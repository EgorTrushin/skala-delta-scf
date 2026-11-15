"""Shared plotting styles and routines for the four manuscript benchmark figures."""

import argparse
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from matplotlib import patheffects, rc_context
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator, MultipleLocator

from table_mae import (BASIS, ROOT, calculate_maes, read_csv)

FUNCTIONALS = ('pbe', 'r2scan', 'pbe0', 'skala')
METHODS = ('UKS', 'ROKS')
QUANTITIES = ('singlet', 'triplet', 'gap')

LABELS = ("PBE", "r²SCAN", "PBE0", "Skala")

COLORS = ("orangered", "dodgerblue", "tab:green", "orange")

TITLES = ("Singlet excitation", "Triplet excitation", "Singlet-triplet gap")

POINT_TEXT_SCALE = 1.44

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "axes.unicode_minus": False,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "hatch.linewidth": 0.6,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
}

SIGNED_ERROR_STYLE = {
    **STYLE,
    "font.family": ["Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "axes.linewidth": 0.55,
    "xtick.labelsize": 7.7,
    "ytick.labelsize": 7.7,
}

def tint_color(color, strength):
    """Mix a color with white without PDF transparency."""
    return tuple(1 - strength * (1 - channel) for channel in to_rgb(color))

def _swarm(ax, values, size, stroke):
    """Separate circular markers in display coordinates without changing errors."""
    diameter = (np.sqrt(size) + stroke + 0.10) * ax.figure.dpi / 72
    origin = ax.transData.transform((0, 0))
    xscale = ax.transData.transform((1, 0))[0] - origin[0]
    y = ax.transData.transform(np.column_stack((np.zeros(len(values)), values)))[:, 1]
    placed = []
    positions = np.zeros(len(values))
    for index in np.argsort(values, kind="stable"):
        neighbours = [(xj, yj) for xj, yj in placed if abs(y[index] - yj) < diameter]
        candidates = [0.0]
        for xj, yj in neighbours:
            separation = diameter * np.sqrt(1 - ((y[index] - yj) / diameter) ** 2)
            candidates.extend((xj - separation, xj + separation))
        for x in sorted(candidates, key=lambda value: (abs(value), value)):
            extent = [x, *(xj for xj, _ in placed)]
            if max(extent) - min(extent) > 0.86 * xscale:
                continue
            if all(
                ((x - xj) / diameter) ** 2 + ((y[index] - yj) / diameter) ** 2
                >= 1 - 1e-9
                for xj, yj in neighbours
            ):
                positions[index] = x / xscale
                placed.append((x, y[index]))
                break
        else:
            raise ValueError("marker footprints exceed the functional's swarm width")
    positions -= (positions.min() + positions.max()) / 2
    return positions

def signed_errors(rows, statistics, *, ylim=None, yticks=(-1.5, -1.0, -0.5, 0)):
    """Pale, outlined circles and rounded mean bars, with matching method styles."""
    fig = Figure(figsize=(7.2, 4.85))
    axes = fig.subplots(2, 3, sharey=True)
    fig.subplots_adjust(
        left=0.09, right=0.985, top=0.95, bottom=0.10, hspace=0.42, wspace=0.08
    )
    stats = {(r["method"], r["quantity"], r["functional"]): r for r in statistics}
    if ylim is None:
        errors = [r["error_ev"] for r in rows]
        ylim = (min(-1.53, min(errors) - 0.1), max(0.4, max(errors) + 0.1))
    marker_size, marker_stroke = 12, 0.5
    for row_index, method in enumerate(METHODS):
        for col_index, (quantity, title) in enumerate(zip(QUANTITIES, TITLES)):
            ax = axes[row_index, col_index]
            ax.set(xlim=(-0.5, 3.5), ylim=ylim, yticks=yticks)
            for position, (xc, color) in enumerate(zip(FUNCTIONALS, COLORS)):
                values = np.array(
                    [
                        r["error_ev"]
                        for r in rows
                        if (r["method"], r["quantity"], r["functional"])
                        == (method, quantity, xc)
                    ]
                )
                ax.scatter(
                    position + _swarm(ax, values, marker_size, marker_stroke),
                    values,
                    s=marker_size,
                    edgecolors=color,
                    facecolors=tint_color(color, 0.22),
                    linewidths=marker_stroke,
                    marker="o",
                    zorder=3,
                )
                mean = stats[method, quantity, xc]["mean_error_ev"]
                ax.plot(
                    [position - 0.23, position + 0.23],
                    [mean, mean],
                    color=color,
                    lw=2.1,
                    solid_capstyle="round",
                    path_effects=[
                        patheffects.Stroke(linewidth=3.1, foreground="white"),
                        patheffects.Normal(),
                    ],
                    zorder=4,
                )
            ax.axhline(0, color="0.55", lw=0.6, ls=(0, (2.5, 2.5)), zorder=1)
            ax.set(xticks=range(4), xticklabels=LABELS)
            ax.grid(False)
            ax.set_axisbelow(True)
            ax.tick_params(
                axis="both",
                width=0.55,
                length=2.8,
                pad=3,
                labelsize=SIGNED_ERROR_STYLE["xtick.labelsize"] * POINT_TEXT_SCALE,
            )
            for spine in ax.spines.values():
                spine.set_color("0.3")
                spine.set_linewidth(0.55)
            panel = 3 * row_index + col_index
            method_label = "ΔSCF-UKS" if method == "UKS" else "ΔSCF-ROKS"
            ax.set_title(
                f"({chr(97 + panel)}) {method_label}\n{title}",
                fontsize=8.3 * POINT_TEXT_SCALE,
                pad=2,
            )
            if col_index == 0:
                ax.set_ylabel(
                    "Calculated - TBE (eV)",
                    fontsize=SIGNED_ERROR_STYLE["axes.labelsize"] * POINT_TEXT_SCALE,
                )
    return fig

GAP_STYLE = {
    **SIGNED_ERROR_STYLE,
    "axes.titlesize": 9,
    "axes.edgecolor": "0.3",
    "xtick.major.width": 0.55,
    "ytick.major.width": 0.55,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
}

METHOD_LABELS = ("ΔSCF-UKS", "ΔSCF-ROKS")

GAP_LABEL = r"$\Delta E_{\mathrm{ST}}$ (eV)"

NEGATIVE_FILL = "#f1f3f5"

def _points(ax, x, y, color):
    return ax.scatter(
        x,
        y,
        s=18,
        marker="o",
        facecolors=tint_color(color, 0.22),
        edgecolors=color,
        linewidths=0.65,
        zorder=3,
    )

def _mean(ax, x, value, color, width=0.23):
    return ax.plot(
        [x - width, x + width],
        [value, value],
        color=color,
        linewidth=2.0,
        solid_capstyle="round",
        zorder=4,
        path_effects=[
            patheffects.Stroke(linewidth=3.0, foreground="white"),
            patheffects.Normal(),
        ],
    )

def _gap_background(ax, limits):
    ax.grid(False)
    ax.axhspan(limits[0], 0, color=NEGATIVE_FILL, zorder=0)
    ax.axhline(0, color="0.4", linewidth=0.7, linestyle=(0, (3, 3)), zorder=1)

def gap_swarm(report):
    """Show signed gaps directly, alongside the eight reference gaps."""
    fig = Figure(figsize=(7.2, 3.45))
    axes = fig.subplots(1, 2, sharey=True)
    fig.subplots_adjust(left=0.095, right=0.98, bottom=0.24, top=0.90, wspace=0.14)
    for m, ax in enumerate(axes):
        ax.set(
            xlim=(-0.55, 4.5),
            ylim=report.gap_limits,
            ylabel=GAP_LABEL if m == 0 else "",
        )
        _gap_background(ax, report.gap_limits)
        values = (*report.gaps[m], report.reference)
        for position, (data, color) in enumerate(
            zip(values, (*COLORS, "0.2"), strict=True)
        ):
            _points(ax, position + _swarm(ax, data, 18, 0.65), data, color)
            _mean(ax, position, data.mean(), color)
        ax.set(xticks=range(5), xticklabels=(*LABELS, "TBE"))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=8))
        ax.tick_params(
            axis="both", labelsize=GAP_STYLE["xtick.labelsize"] * POINT_TEXT_SCALE
        )
        ax.yaxis.label.set_fontsize(GAP_STYLE["axes.labelsize"] * POINT_TEXT_SCALE)
        ax.set_title(
            f"({chr(97 + m)}) {METHOD_LABELS[m]}",
            fontsize=GAP_STYLE["axes.titlesize"] * POINT_TEXT_SCALE,
            pad=5,
        )
    return fig

MARKERS = ("o", "s", "D", "^")

SCATTER_STYLE = {
    **STYLE,
    "font.size": 8.3,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
}

def energy_limits(rows):
    """Use the same padded range for calculated and reference energies."""
    values = [r[k] for r in rows for k in ("value_ev", "reference_ev")]
    return (np.floor((min(values) - 0.05) * 2) / 2,
            np.ceil((max(values) + 0.05) * 2) / 2)

def scatter_plot(rows, statistics):
    """Show exact energies with pale markers, equal axes and eleven-target MAEs."""
    # Size at the intended paper column width so labels need no downscaling.
    fig = Figure(figsize=(3.5, 3.6))
    ax = fig.subplots()
    fig.subplots_adjust(left=0.16, right=0.95, bottom=0.18, top=0.96)
    limits = energy_limits(rows)
    ax.plot(limits, limits, color="0.5", lw=0.7, ls=(0, (4, 4)), zorder=1)
    stats = {r["functional"]: r for r in statistics}
    for xc, label, color, marker in zip(FUNCTIONALS, LABELS, COLORS, MARKERS, strict=True):
        group = [r for r in rows if r["functional"] == xc]
        ax.scatter(
            [r["reference_ev"] for r in group], [r["value_ev"] for r in group],
            label=f"{label}  MAE={stats[xc]['mae_ev']:.3f}",
            s=25, marker=marker, facecolor=tint_color(color, 0.23),
            edgecolor=color, linewidth=0.85,
            zorder=3,
        )
    ax.set(xlabel="Reference excitation energy (eV)",
           ylabel="Calculated excitation energy (eV)",
           xlim=limits, ylim=limits, aspect="equal")
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.grid(False)
    ax.tick_params(direction="out", length=3, width=0.65, pad=3)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelspacing=0.4,
              handlelength=1.0, handletextpad=0.5, columnspacing=1.0, borderaxespad=0.5)
    return fig


def load_data(calculations, references):
    """Use the same full-subset validation and MAEs as the main-text table."""
    maes = calculate_maes(calculations, references)
    reference_rows = read_csv(references, ('dataset', 'molecule', 'singlet_eV', 'triplet_eV'))
    targets = {}
    for r in reference_rows:
        ref = {'singlet': float(r['singlet_eV'])}
        if r['triplet_eV']:
            ref['triplet'] = float(r['triplet_eV'])
            ref['gap'] = ref['singlet'] - ref['triplet']
        targets[r['dataset'], r['molecule']] = ref
    raw = read_csv(calculations, ('dataset', 'basis', 'method', 'functional', 'molecule'))
    lookup = {(r['dataset'], r['method'], r['functional'].lower(), r['molecule']): r
              for r in raw if r['basis'] == BASIS}
    reports = {}
    for dataset in ('QUEST1', 'ST_gap', 'QUEST6'):
        systems = [r['molecule'] for r in reference_rows if r['dataset'] == dataset]
        methods = ('UKS',) if dataset == 'QUEST6' else METHODS
        quantities = ('singlet',) if dataset == 'QUEST6' else QUANTITIES
        rows, stats = [], []
        for method in methods:
            for xc in FUNCTIONALS:
                for quantity in quantities:
                    group = []
                    for system in systems:
                        r = lookup[dataset, method, xc, system]
                        value = float(r[quantity + '_eV'])
                        ref = targets[dataset, system][quantity]
                        group.append(dict(system=system, method=method, functional=xc,
                                          quantity=quantity, value_ev=value, reference_ev=ref,
                                          error_ev=value-ref))
                    rows.extend(group)
                    label = {'pbe':'PBE','r2scan':'r2SCAN','pbe0':'PBE0','skala':'Skala'}[xc]
                    stats.append(dict(method=method, functional=xc, quantity=quantity,
                                      mean_error_ev=float(np.mean([r['error_ev'] for r in group])),
                                      mae_ev=maes[dataset, method, quantity, label]))
        reports[dataset] = SimpleNamespace(rows=rows, statistics=stats, systems=systems)
    report = reports['ST_gap']
    report.gaps = np.array([[[float(lookup['ST_gap', method, xc, system]['gap_eV'])
                             for system in report.systems] for xc in FUNCTIONALS] for method in METHODS])
    report.reference = np.array([targets['ST_gap', system]['gap'] for system in report.systems])
    if not np.all(report.reference < 0):
        raise ValueError('ST-gap distribution requires negative reference gaps')
    values = np.append(report.gaps.ravel(), report.reference)
    report.gap_limits = (np.floor((min(0, values.min()) - 0.025) * 10) / 10,
                         np.ceil((max(0, values.max()) + 0.025) * 10) / 10)
    return reports


def main(kind, argv=None):
    parser = argparse.ArgumentParser(description=f'Reproduce {kind} figures from CP2K / {BASIS} CSV data.')
    parser.add_argument('--calculations', type=Path, default=ROOT / 'data/csv/cp2k_excitations.csv')
    parser.add_argument('--references', type=Path, default=ROOT / 'data/references/excitations.csv')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'output/figures')
    parser.add_argument('--png', action='store_true', help='Also save 300-dpi PNG previews')
    args = parser.parse_args(argv)
    try:
        calculations, references = args.calculations.expanduser(), args.references.expanduser()
        reports = load_data(calculations, references)
        output = args.output_dir.expanduser()
        output.mkdir(parents=True, exist_ok=True)
        metadata = dict(Creator='skala-delta-scf', CreationDate=None, ModDate=None,
                        Subject=f'CP2K / {BASIS}; calculation CSV SHA256 '
                        f'{hashlib.sha256(calculations.read_bytes()).hexdigest()}; reference CSV SHA256 '
                        f'{hashlib.sha256(references.read_bytes()).hexdigest()}')
        # Keep the reviewed dimensions, colors, marker shapes and tight cropping.
        with rc_context(STYLE):
            if kind == 'quest1':
                r = reports['QUEST1']
                with rc_context(SIGNED_ERROR_STYLE):
                    figures = [('quest1_signed_errors', signed_errors(r.rows, r.statistics))]
            elif kind == 'st_gap':
                r = reports['ST_gap']
                bound = np.ceil((max(abs(row['error_ev']) for row in r.rows) + 0.025) * 10) / 10
                tick_count = int(np.floor(bound / 0.2))
                with rc_context(SIGNED_ERROR_STYLE):
                    errors = signed_errors(r.rows, r.statistics, ylim=(-bound, bound),
                                           yticks=0.2*np.arange(-tick_count, tick_count+1))
                with rc_context(GAP_STYLE):
                    gaps = gap_swarm(r)
                figures = [('st_gap_signed_errors', errors), ('st_gap_gap_distributions', gaps)]
            elif kind == 'quest6':
                r = reports['QUEST6']
                with rc_context(SCATTER_STYLE):
                    figures = [('quest6_scatter', scatter_plot(r.rows, r.statistics))]
            else:
                raise ValueError(f'Unknown figure group: {kind}')
            for name, fig in figures:
                path = output / (name + '.pdf')
                fig.savefig(path, metadata={**metadata, 'Title':name}, bbox_inches='tight', pad_inches=0.05)
                print(f'Wrote {path}')
                if args.png:
                    path = output / (name + '.png')
                    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.05)
                    print(f'Wrote {path}')
                fig.clear()
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f'{exc}\n')
    return 0
