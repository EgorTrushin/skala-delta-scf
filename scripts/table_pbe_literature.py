"""Print CP2K versus published PBE excitation energies and export the comparison."""

import argparse
import csv
import math
from pathlib import Path

from table_mae import ROOT, number, read_csv

SOURCES = {'QUEST1': '10.1021/acs.jctc.1c00157', 'QUEST6': '10.1021/acs.jctc.3c01319'}
FIELDS = ('dataset', 'molecule', 'quantity', 'literature_eV', 'cp2k_eV',
          'difference_eV', 'doi', 'table', 'printed_page', 'column')


def compare(calculations, literature):
    cp = {}
    for r in read_csv(calculations, ('dataset','molecule','method','functional','basis','singlet_eV','triplet_eV')):
        if (r['basis'], r['method'], r['functional']) != ('aug-cc-pVTZ', 'UKS', 'PBE'):
            continue
        key = r['dataset'], r['molecule']
        if key in cp:
            raise ValueError(f'Duplicate CP2K result: {key}')
        cp[key] = r
    rows, seen = [], set()
    for r in read_csv(literature, ('dataset','molecule','quantity','literature_eV','doi','table','printed_page','column')):
        d, m, q = r['dataset'], r['molecule'], r['quantity']
        key = d, m, q
        if d not in SOURCES or q not in ('singlet','triplet') or (d == 'QUEST6' and q != 'singlet'):
            raise ValueError(f'Unsupported literature target: {key}')
        if key in seen or r['doi'] != SOURCES[d]:
            raise ValueError(f'Duplicate target or incorrect source: {key}')
        seen.add(key)
        if (d, m) not in cp:
            raise ValueError(f'Missing CP2K result: {key}')
        published = number(r['literature_eV'], str(key))
        calculated = number(cp[d,m][q+'_eV'], str(key))
        rows.append(dict(r, literature_eV=published, cp2k_eV=calculated,
                         difference_eV=calculated-published))
    quest1 = {m for d,m,q in seen if d == 'QUEST1'}
    quest6 = {m for d,m,q in seen if d == 'QUEST6'}
    if len(quest1) != 7 or len(quest6) != 11:
        raise ValueError('Expected seven QUEST1 and eleven QUEST6 matched systems')
    expected = {('QUEST1',m,q) for m in quest1 for q in ('singlet','triplet')}
    expected |= {('QUEST6',m,'singlet') for m in quest6}
    if seen != expected:
        raise ValueError('Incomplete literature singlet/triplet coverage')
    if quest6 != {m for d,m in cp if d=='QUEST6'}:
        raise ValueError('QUEST6 literature and CP2K subsets differ')
    return rows


def summary(rows):
    summaries=[]
    for dataset, quantity in [('QUEST1','singlet'),('QUEST1','triplet'),('QUEST1','all'),('QUEST6','singlet')]:
        values=[abs(r['difference_eV']) for r in rows if r['dataset']==dataset
                and (quantity=='all' or r['quantity']==quantity)]
        summaries.append(dict(dataset=dataset,quantity=quantity,count=len(values),
                              mad_eV=math.fsum(values)/len(values),maxad_eV=max(values)))
    return summaries


def print_table(rows):
    widths=[max(len(str(row[i])) for row in rows) for i in range(len(rows[0]))]
    for row in rows:
        print('  '.join(str(v).ljust(widths[i]) for i,v in enumerate(row)))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calculations',type=Path,default=ROOT/'data/csv/cp2k_excitations.csv')
    parser.add_argument('--literature',type=Path,default=ROOT/'data/references/pbe_literature.csv')
    parser.add_argument('--output',type=Path,default=ROOT/'output/tables/pbe_literature.csv')
    args=parser.parse_args(argv)
    try:
        rows=compare(args.calculations.expanduser(),args.literature.expanduser())
        stats=summary(rows)
        output=args.output.expanduser()
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=FIELDS);writer.writeheader()
            for r in rows:
                writer.writerow({k:f'{r[k]:.12f}' if isinstance(r[k],float) else r[k] for k in FIELDS})
        summary_path=output.with_name(output.stem+'_summary.csv')
        with summary_path.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(stats[0]));writer.writeheader();writer.writerows(stats)
    except (OSError,ValueError) as exc:
        parser.exit(1,f'{exc}\n')
    print('CP2K UKS/PBE/aug-cc-pVTZ versus published orbital-optimized PBE results (eV)')
    print('Singlets are spin purified. Delta = CP2K - literature.\n')
    for dataset in SOURCES:
        print(dataset + (' — Ivanov et al. (2021)' if dataset=='QUEST1' else ' — Selenius et al. (2024)'))
        print_table([['Molecule','State','Literature','CP2K','Delta']]+[
            [r['molecule'],r['quantity'],f"{r['literature_eV']:.2f}",f"{r['cp2k_eV']:.4f}",f"{r['difference_eV']:+.4f}"]
            for r in rows if r['dataset']==dataset])
        print()
    print_table([['Dataset','States','Count','MAD (eV)','MaxAD (eV)']]+[
        [r['dataset'],r['quantity'],r['count'],f"{r['mad_eV']:.6f}",f"{r['maxad_eV']:.6f}"] for r in stats])
    print('\nPublished energies have two decimal places; numerical protocols differ.')
    print('Agreement supports consistency of the targets, not independent state assignment.')
    print(f'\nWrote {output}\nWrote {summary_path}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
