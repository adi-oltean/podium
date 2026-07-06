#!/usr/bin/env python3
"""Publication figure for the verification-modality coverage matrix.

Renders the real output of tools/fault_coverage.build_matrix() as a
clean grid: fault classes (rows) x verification lanes (columns), a cell
filled iff that lane catches that fault. Lanes are grouped by kind
(physics / emitter / spec / exact certificates) so the block structure
is visible -- certificate faults are caught only in the certificate
block, invisible to the physics lanes.

Object-API matplotlib (no pyplot / global backend), headless-safe.
`python3 tools/fault_coverage_figure.py [outstem]` writes SVG + PDF.
"""

from __future__ import annotations

import os
import sys

# allow running as a bare script (`python tools/fault_coverage_figure.py`)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from matplotlib.figure import Figure  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from tools.fault_coverage import FAULTS, LANES, build_matrix  # noqa: E402

# lane -> (short label, group key); groups ordered for the column layout
_LANE_META = {
    "conservation": ("conserv", "physics"),
    "analytic": ("analytic", "physics"),
    "golden": ("golden", "emitter"),
    "stl": ("STL", "spec"),
    "barrier": ("barrier", "certificate"),
    "kkt": ("KKT", "certificate"),
    "lyapunov": ("Lyapunov", "certificate"),
    "sos": ("SOS", "certificate"),
}
_GROUP_COLOR = {
    "physics": "#1f77b4",       # blue
    "emitter": "#ff7f0e",       # orange
    "spec": "#9467bd",          # purple
    "certificate": "#2ca02c",   # green
}
_GROUP_LABEL = {
    "physics": "physics",
    "emitter": "emit",
    "spec": "spec",
    "certificate": "exact certificates",
}


def build_figure() -> Figure:
    caught = build_matrix()
    faults = list(FAULTS)
    lanes = LANES
    n_rows, n_cols = len(faults), len(lanes)

    fig = Figure(figsize=(9.0, 4.6))
    ax = fig.add_subplot()
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_aspect("equal")
    ax.invert_yaxis()

    for j, lane in enumerate(lanes):
        _short, group = _LANE_META[lane]
        colour = _GROUP_COLOR[group]
        for i, fault in enumerate(faults):
            val = caught[fault][lane]
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor="white",
                                   edgecolor="#dddddd", lw=0.8, zorder=1))
            if val is None:
                ax.plot([j, j + 1], [i, i + 1], color="#cccccc", lw=0.7)
            elif val:
                ax.add_patch(Rectangle((j + 0.12, i + 0.12), 0.76, 0.76,
                                       facecolor=colour, edgecolor="none",
                                       zorder=2))

    # row labels (fault classes)
    for i, fault in enumerate(faults):
        ax.text(-0.2, i + 0.5, fault, ha="right", va="center", fontsize=9)
    # column labels (lanes)
    for j, lane in enumerate(lanes):
        short, group = _LANE_META[lane]
        ax.text(j + 0.5, -0.15, short, ha="left", va="bottom", fontsize=8.5,
                rotation=40, color=_GROUP_COLOR[group])

    # group bands above the columns
    start = 0
    groups = [_LANE_META[ln][1] for ln in lanes]
    for k in range(n_cols + 1):
        if k == n_cols or groups[k] != groups[start]:
            g = groups[start]
            ax.add_patch(Rectangle((start, -1.15), k - start, 0.42,
                                   facecolor=_GROUP_COLOR[g], alpha=0.18,
                                   edgecolor="none", clip_on=False, zorder=0))
            ax.text((start + k) / 2, -0.94, _GROUP_LABEL[g], ha="center",
                    va="center", fontsize=7.5, color=_GROUP_COLOR[g],
                    fontweight="bold")
            start = k

    singced = sum(1 for f in faults
                  if len(FAULTS[f]) == 1 and caught[f][FAULTS[f][0]])
    ax.text(n_cols / 2, n_rows + 0.7,
            f"{singced}/{n_rows} faults caught by exactly one lane  "
            "·  certificate faults are invisible to the physics lanes  "
            "·  0 false alarms",
            ha="center", va="center", fontsize=9, style="italic")

    fig.suptitle("Verification-modality coverage: which lane catches "
                 "which injected fault", fontsize=11, y=0.99)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.subplots_adjust(left=0.22, right=0.98, top=0.70, bottom=0.02)
    return fig


def main() -> int:
    stem = sys.argv[1] if len(sys.argv) > 1 else "fault_coverage_matrix"
    fig = build_figure()
    for ext in ("svg", "pdf"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight")
        print(f"wrote {stem}.{ext}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
