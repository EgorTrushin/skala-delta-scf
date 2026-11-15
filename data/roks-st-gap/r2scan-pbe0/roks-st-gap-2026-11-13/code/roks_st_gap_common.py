"""Portable frozen-ROKS primitives (NumPy only).

Adapted from the validated November 6 CP2K procedure. Frozen diagnostics are
not converged UKS states. AO metrics require complete, well-conditioned MOs.
"""
from io import BytesIO
from pathlib import Path
import struct
import re
import json
import hashlib
import gzip
import numpy as np

ORTHO_TOL = 1e-7
METRIC_CONDITION_LIMIT = 1e6
DENSITY_TOL = 1e-8
HARTREE_TO_EV = 27.211396132

def read_wfn(path):
    """Read the real, little-endian, 4-byte-record restart used by our build."""
    return read_wfn_bytes(Path(path).read_bytes())

def read_wfn_bytes(data):
    """Read one immutable restart snapshot, including its serialized AO layout."""
    stream = BytesIO(data)

    def record(dtype, count=None):
        marker = stream.read(4)
        if len(marker) != 4:
            raise ValueError("truncated WFN record marker")
        size = struct.unpack('<i', marker)[0]
        itemsize = np.dtype(dtype).itemsize
        if size <= 0 or size % itemsize or (count is not None and size != count * itemsize):
            raise ValueError("invalid WFN record size or unsupported restart format")
        data = stream.read(size)
        if len(data) != size or stream.read(4) != marker:
            raise ValueError("truncated WFN or mismatched Fortran record markers")
        return np.frombuffer(data, dtype=dtype).copy()

    header = record('<i4', 5)
    natom, nspin, nao, nset, nshell = map(int, header)
    if min(natom, nao, nset, nshell) <= 0 or nspin not in (1, 2):
        raise ValueError("unsupported WFN header")
    basis = [record('<i4') for _ in range(3)]
    spins = []
    for _ in range(nspin):
        nmo, homo, lfomo, nelectron = map(int, record('<i4', 4))
        if not 0 < nmo <= nao or nelectron < 0:
            raise ValueError("invalid WFN orbital/electron count")
        values = record('<f8', 2 * nmo)
        coeff = np.column_stack([record('<f8', nao) for _ in range(nmo)])
        if not np.isfinite(values).all() or not np.isfinite(coeff).all():
            raise ValueError("non-finite WFN coefficients, energies or occupations")
        occ = values[nmo:]
        full = 2 if nspin == 1 else 1
        if not np.all((abs(occ) < 1e-8) | (abs(occ - full) < 1e-8)):
            raise ValueError("fractional occupations are not supported by this screen")
        if abs(float(occ.sum()) - nelectron) > 1e-8:
            raise ValueError("WFN occupations disagree with its electron count")
        spins.append({'C': coeff, 'eps': values[:nmo], 'occ': occ,
                      'homo': homo, 'lfomo': lfomo, 'nelectron': nelectron})
    if stream.read(1):
        raise ValueError("unparsed trailing WFN data; unsupported restart format")
    result = {'header': header, 'basis': basis, 'spins': spins}
    ao_atom_slices(result)  # Validate the serialized basis layout too.
    return result

def wfn_bytes(wfn):
    """Serialize a prepared real WFN in the same format as the restart reader."""
    def record(array, dtype):
        payload = np.asarray(array, dtype=dtype).tobytes()
        marker = struct.pack('<i', len(payload))
        return marker + payload + marker

    records = [record(wfn['header'], '<i4')]
    records.extend(record(array, '<i4') for array in wfn['basis'])
    for spin in wfn['spins']:
        nmo = spin['C'].shape[1]
        records.append(record([nmo, spin['homo'], spin['lfomo'], spin['nelectron']], '<i4'))
        records.append(record(np.r_[spin['eps'], spin['occ']], '<f8'))
        records.extend(record(column, '<f8') for column in spin['C'].T)
    return b''.join(records)

def ao_atom_slices(wfn):
    """AO intervals in CP2K atom order, from the Fortran-order basis records."""
    natom, _, nao, nset, nshell = map(int, wfn['header'])
    sets, shells, functions = wfn['basis']
    if [len(a) for a in wfn['basis']] != [natom, nset * natom, nshell * nset * natom]:
        raise ValueError('invalid WFN basis record dimensions')
    shells = shells.reshape((nset, natom), order='F')
    functions = functions.reshape((nshell, nset, natom), order='F')
    if (np.any(sets < 0) or np.any(sets > nset) or np.any(shells < 0)
            or np.any(shells > nshell) or np.any(functions < 0)):
        raise ValueError('invalid WFN basis shell counts')
    slices, start = [], 0
    for atom, count in enumerate(sets):
        if np.any(shells[count:, atom]) or np.any(functions[:, count:, atom]):
            raise ValueError('nonzero unused WFN basis sets')
        size = 0
        for block in range(count):
            active = shells[block, atom]
            if (np.any(functions[:active, block, atom] <= 0)
                    or np.any(functions[active:, block, atom])):
                raise ValueError('invalid WFN basis shell layout')
            size += int(functions[:active, block, atom].sum())
        slices.append(slice(start, start + size))
        start += size
    if start != nao:
        raise ValueError('WFN basis layout disagrees with AO count')
    return slices

def recover_overlap(wfn):
    """Recover S = C^-T C^-1 from a complete, well-conditioned real MO basis.

    A rectangular matrix cannot determine the full metric. No pseudoinverse
    is used. Orthonormality of the source is an assumption, not independent
    validation; other saved states must be cross-checked against this metric.
    """
    failures = []
    for spin, data in enumerate(wfn['spins']):
        c = data['C']
        if c.shape[0] != c.shape[1]:
            continue
        singular = np.linalg.svd(c, compute_uv=False)
        condition = float(singular[0] / singular[-1]) if singular[-1] > 0 else float('inf')
        if not np.isfinite(condition) or condition > METRIC_CONDITION_LIMIT:
            failures.append(f'spin {spin + 1}: singular/ill-conditioned MO matrix')
            continue
        inverse = np.linalg.solve(c, np.eye(len(c)))
        overlap = inverse.T @ inverse
        try:
            error = orthonormality_error(wfn, overlap)
        except ValueError as exc:
            failures.append(str(exc))
            continue
        return overlap, {
            'kind': 'recovered_from_wfn', 'spin_channel': spin + 1,
            'coefficient_condition_number': condition,
            'source_orthonormality_error': error,
            'assumption': 'Source contains a complete AO-metric orthonormal MO basis; source self-check is not independent validation.',
        }
    raise ValueError('; '.join(failures) if failures else 'no complete square MO matrix; supply an AO overlap cache')

def orthonormality_error(wfn, overlap):
    errors = []
    for spin in wfn['spins']:
        c = spin['C']
        if c.shape[0] != overlap.shape[0]:
            raise ValueError("WFN / AO overlap dimensions differ")
        errors.append(float(np.max(abs(c.T @ overlap @ c - np.eye(c.shape[1])))))
    error = max(errors)
    if error > ORTHO_TOL:
        raise ValueError(f"MO orthonormality error {error:.3g} exceeds {ORTHO_TOL:g}")
    return error

def energy_components(text, scf_energy):
    """Full-precision breakdown, checked against the first SCF row."""
    labels = {
        'total': r'Total energy',
        'xc_grid': r'Exchange-correlation energy',
        'xc_atomic': r'GAPW\| Exc from hard and soft atomic rho1',
        'hartree_grid': r'Hartree energy',
        'hartree_atomic': r'GAPW\| local Eh = 1 center integrals',
        'core_hamiltonian': r'Core Hamiltonian energy',
        'core_overlap': r'Overlap energy of the core charge distribution',
        'core_self': r'Self energy of the core charge distribution',
    }
    values = {}
    for name, label in labels.items():
        matches = re.findall(r'^\s*' + label + r':\s*([-+\d.E]+)\s*$', text, re.M)
        if len(matches) != 1:
            raise ValueError(f'Expected one {name} energy component, got {len(matches)}')
        values[name] = float(matches[0])
    hf = re.findall(r'^\s*Hartree-Fock Exchange energy:\s*([-+\d.E]+)\s*$', text, re.M)
    if len(hf) > 1:
        raise ValueError('Repeated exact-exchange energy')
    values['exact_exchange'] = float(hf[0]) if hf else 0.0
    if not all(np.isfinite(v) for v in values.values()):
        raise ValueError('Nonfinite energy component')
    if abs(values['total'] - scf_energy) > 1e-9:
        raise ValueError('Energy breakdown differs from the first SCF row')
    if abs(values['total'] - sum(v for k, v in values.items() if k != 'total')) > 1e-9:
        raise ValueError('Incomplete energy decomposition')
    return values

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def explicit_determinant(source, configuration):
    """Preserve labeled shared orbitals and construct explicit M or T spins."""
    if configuration not in ('M', 'T') or len(source['spins']) != 2:
        raise ValueError('Expected a two-spin ROKS checkpoint and M or T')
    alpha, beta = source['spins']
    nc = beta['nelectron']
    if alpha['nelectron'] != nc + 2:
        raise ValueError('Expected exactly two open-shell orbitals')
    if not np.all(alpha['occ'][:nc + 2] == 1) or not np.all(beta['occ'][:nc] == 1):
        raise ValueError('Expected contiguous occupied ROKS orbitals')
    c = alpha['C']
    if np.max(abs(c[:, :nc] @ c[:, :nc].T -
                  beta['C'][:, :nc] @ beta['C'][:, :nc].T)) > 1e-7:
        raise ValueError('ROKS closed spin densities do not match')
    groups = (list(range(nc)) + [nc], list(range(nc)) + [nc + 1]) if configuration == 'M' else (
        list(range(nc + 2)), list(range(nc)))
    result = dict(header=source['header'].copy(), basis=source['basis'], spins=[])
    for indices in groups:
        n = len(indices)
        result['spins'].append(dict(C=c[:, indices].copy(), eps=alpha['eps'][indices].copy(),
                                   occ=np.ones(n), homo=n, lfomo=n + 1, nelectron=n))
    original = sum((s['C'] * s['occ']) @ s['C'].T for s in source['spins'])
    prepared = sum((s['C'] * s['occ']) @ s['C'].T for s in result['spins'])
    if np.max(abs(original - prepared)) > 1e-7:
        raise ValueError('Prepared determinant changes the total ROKS density')
    return result

def replace_section(text, name, replacement):
    pattern = rf'(?ms)^([ \t]*)&{name}\b.*?^\1&END {name}[ \t]*\n'
    updated, count = re.subn(pattern, replacement, text, count=1)
    if count != 1:
        raise ValueError(f'Expected section {name}')
    return updated

def probe_input(template, configuration):
    text = replace_section(template, 'LOW_SPIN_ROKS', '')
    text = re.sub(r'(?m)^\s*WFN_RESTART_FILE_NAME.*$', '    WFN_RESTART_FILE_NAME seed.wfn', text)
    text = text.replace('ROKS T', 'UKS T')
    text = re.sub(r'(?m)^\s*MULTIPLICITY.*$', '    MULTIPLICITY ' + ('1' if configuration == 'M' else '3'), text)
    text = text.replace('SKALA_RUNTIME SELF', 'SKALA_RUNTIME MPI')
    text = replace_section(text, 'SCF', '''    &SCF
      SCF_GUESS RESTART
      MAX_SCF 1
      EPS_SCF 1.0E-30
      &OT
        MINIMIZER CG
        PRECONDITIONER NONE
        STEPSIZE 0.001
      &END OT
    &END SCF
''')
    text = replace_section(text, 'PRINT', '''    &PRINT
      &AO_MATRICES ON
        FILENAME =density.Log
        NDIGITS 16
        DENSITY T
        &EACH
          QS_SCF 1
        &END EACH
      &END AO_MATRICES
    &END PRINT
''')
    text = re.sub(r'(?m)^\s*PROJECT\s+.*$', '  PROJECT frozen', text)
    return text

def density_matrices(path, nao):
    """Read a single pre-update alpha/beta density dump without filling gaps."""
    matrices, counts = {}, {}
    key, columns = None, None
    for line in density_bytes(path).decode().splitlines():
        match = re.fullmatch(r'\s*DENSITY MATRIX FOR (ALPHA|BETA) SPIN\s*', line)
        if match:
            key = match[1]
            if key in matrices:
                raise ValueError('More than one density evaluation')
            matrices[key] = np.full((nao, nao), np.nan)
            counts[key] = np.zeros((nao, nao), dtype=np.uint8)
            columns = None
        elif not line.strip():
            continue
        elif re.fullmatch(r'\s*\d+(?:\s+\d+)*\s*', line):
            columns = np.array([int(v) - 1 for v in line.split()])
            if key is None or len(columns) not in (1, 2) or not np.all((columns >= 0) & (columns < nao)):
                raise ValueError('Invalid density column header')
        else:
            fields = line.split()
            if key is None or columns is None or len(fields) != 4 + len(columns) or not fields[0].isdigit():
                raise ValueError('Unexpected matrix line: ' + line[:120])
            row = int(fields[0]) - 1
            if not 0 <= row < nao or counts[key][row, columns].any():
                raise ValueError('Invalid or repeated density element')
            matrices[key][row, columns] = [float(v.replace('D', 'E')) for v in fields[4:]]
            counts[key][row, columns] += 1
    if set(matrices) != {'ALPHA', 'BETA'}:
        raise ValueError('Missing spin density')
    for key, matrix in matrices.items():
        if not np.all(counts[key] == 1) or not np.isfinite(matrix).all() or np.max(abs(matrix - matrix.T)) > 1e-12:
            raise ValueError('Incomplete or asymmetric density matrix')
    return matrices

def density_bytes(path):
    path = Path(path)
    data = path.read_bytes()
    return gzip.decompress(data) if path.suffix == '.gz' else data

def validate(directory):
    directory = Path(directory)
    info = json.loads((directory / 'trial.json').read_text())
    if sha(directory / 'seed.wfn') != info['seed_sha256'] or sha(directory / 'probe.inp') != info['input_sha256']:
        raise ValueError('Input/seed changed')
    source = read_wfn(directory / 'seed.wfn')
    density_path = directory / 'density.Log'
    if not density_path.exists():
        density_path = directory / 'density.Log.gz'
    matrices = density_matrices(density_path, int(source['header'][2]))
    differences, relative_errors = {}, {}
    for spin, label in zip(source['spins'], ('ALPHA', 'BETA')):
        expected = (spin['C'] * spin['occ']) @ spin['C'].T
        differences[label] = float(np.max(abs(expected - matrices[label])))
        relative_errors[label] = float(np.linalg.norm(expected - matrices[label]) /
                                       np.linalg.norm(expected))
    if max(*differences.values(), *relative_errors.values()) > DENSITY_TOL:
        raise ValueError('Evaluated density differs from frozen seed: ' + str(differences))
    text = (directory / 'probe.out').read_text()
    matches = re.findall(r'^\s*(\d+)\s+OT\s+\S+\s+\S+\s+\S+\s+\S+\s+([-+\d.E]+)\s+[-+\d.E]+\s*$', text, re.M)
    if len(matches) != 1 or matches[0][0] != '1':
        raise ValueError(f'Expected exactly one pre-update OT energy: {matches}')
    if 'SCF run NOT converged' not in text and 'SCF run converged in' not in text:
        raise ValueError('Missing SCF termination diagnostic')
    energy = float(matches[0][1])
    if not np.isfinite(energy):
        raise ValueError('Nonfinite frozen energy')
    result = dict(**info, energy_Ha=energy, density_max_difference=differences,
                  density_relative_error=relative_errors, density_tolerance=DENSITY_TOL,
                  validated=True, output_sha256=sha(directory / 'probe.out'),
                  density_sha256=hashlib.sha256(density_bytes(density_path)).hexdigest(),
                  qualification='First pre-update energy with verified source density; not an SCF-converged state.')
    (directory / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(id=info['id'], energy_Ha=energy, density_max_difference=differences)))

