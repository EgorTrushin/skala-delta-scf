"""Reproduce the frozen-orbital ROKS manuscript figure from the archived CSV."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from matplotlib import rc_context
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from plot_common import SIGNED_ERROR_STYLE
from table_roks import ROOT, MOLECULES, load_data

OUTPUT_STEM = "st_gap_roks_analysis"
PUBLICATION_STYLE = {
    **SIGNED_ERROR_STYLE,
    "font.size": 8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,
    "mathtext.fontset": "dejavusans",
    "hatch.linewidth": 0.45,
}
STYLES = (("PBE@PBE", "orangered"), ("Skala@PBE", "dodgerblue"),
          ("Skala@Skala", "orange"))
PAIRS = (("pbe", "pbe"), ("skala", "pbe"), ("skala", "skala"))


def load_gaps(path):
    rows = load_data(path)
    gaps = np.array([[rows[m, evaluator, donor]["gap_eV"]
                      for evaluator, donor in PAIRS] for m in MOLECULES])
    if np.any(gaps < -0.085) or np.any(gaps > 0.175):
        raise ValueError("Values exceed the publication axis limits")
    return SimpleNamespace(systems=MOLECULES, gaps=gaps)


def grouped_gap_figure(data, width=3.4):
    """Plot signed gaps in the original benchmark molecule and fixed stage order."""
    fig = Figure(figsize=(width, 3.1))
    ax = fig.subplots()
    fig.subplots_adjust(left=0.58 / width, right=1 - 0.07 / width,
                        bottom=0.28, top=0.89)
    order = np.arange(len(data.systems))
    values = data.gaps[order]
    x = np.arange(len(order))
    bar_width = 0.25
    for stage, (label, color) in enumerate(STYLES):
        ax.bar(x + (stage - 1) * bar_width, values[:, stage], width=bar_width,
               color=color, edgecolor="black", linewidth=0.45,
               label=label, zorder=2)
    ax.axhline(0, color="0.2", linewidth=0.65, zorder=3)
    ax.set(xlim=(-0.55, len(order) - 0.45), ylim=(-0.085, 0.175),
           ylabel=r"$\Delta E_{\mathrm{ST}}$ (eV)", yticks=(-0.05, 0, 0.05, 0.10, 0.15),
           xticks=x, xticklabels=[data.systems[i].replace("Molecule", "Molecule ") for i in order])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.tick_params(axis="x", length=0, pad=1)
    ax.tick_params(axis="y", length=3, width=0.55, pad=3)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=50, ha="right", rotation_mode="anchor")
    ax.spines["bottom"].set_visible(False)
    for side in ("left", "top", "right"):
        ax.spines[side].set_color("0.35")
    ax.grid(False)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3,
              frameon=False, borderaxespad=0, handlelength=1.1,
              handletextpad=0.4, columnspacing=0.7)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/csv/roks_st_gap.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/figures")
    parser.add_argument("--png", action="store_true", help="Also save a PNG preview")
    parser.add_argument("--width", type=float, default=3.4)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    if not np.isfinite(args.width) or args.width < 3.4 or args.dpi < 150:
        parser.error("Use a finite width >= 3.4 inches and dpi >= 150")
    data = load_gaps(args.data)
    digest = hashlib.sha256(args.data.read_bytes()).hexdigest()
    metadata = {
        "Title": "ST-gap ROKS frozen-orbital comparison",
        "Creator": "skala-delta-scf: plot_st_gap_frozen.py",
        "Subject": "CP2K ROKS/aug-cc-pVTZ; fixed geometries; "
                   "PBE@PBE -> Skala@PBE -> Skala@Skala; CSV SHA256 " + digest,
        "CreationDate": None, "ModDate": None,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    with rc_context(PUBLICATION_STYLE):
        fig = grouped_gap_figure(data, args.width)
        for extension in (["pdf", "png"] if args.png else ["pdf"]):
            path = args.output_dir / f"{OUTPUT_STEM}.{extension}"
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight", pad_inches=0.05,
                        metadata=metadata if extension == "pdf" else None)
            outputs.append(path)
        fig.clear()
    report = dict(source_csv_sha256=digest, molecules=data.systems,
                  stages=[s[0] for s in STYLES], gaps_eV=data.gaps.tolist(),
                  mean_gaps_eV=data.gaps.mean(axis=0).tolist(),
                  inversion_counts=(data.gaps < 0).sum(axis=0).tolist(),
                  canvas_width_inches=args.width, canvas_height_inches=3.1,
                  bbox_inches="tight", pad_inches=0.05,
                  outputs_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs})
    (args.output_dir / f"{OUTPUT_STEM}_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print("\n".join(str(p) for p in outputs))
    print("Mean gaps (eV):", report["mean_gaps_eV"])
    print("Inversion counts:", report["inversion_counts"])

if __name__ == "__main__":
    main()
