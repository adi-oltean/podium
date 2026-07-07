"""Generate the exact optimality-gap recovery-rate figure for the paper.

Reproduces, from the real exact-rational certificate code in
``podium.verify.bracket``, the quadratic (nonsingular) vs linear (trust-region
hard case) recovery of the certified lower bound as the dual multiplier
approaches its optimum. Every plotted value is an exact rational verified by
``certify_lower_bound``. Output is strict black-and-white (series distinguished
by line style and marker, not color).

Run:  python tools/recovery_figure.py docs/paper/recovery.pdf
"""
import sys
from fractions import Fraction as F

sys.path.insert(0, "src")
from podium.verify import bracket  # noqa: E402


def smooth_series():
    # min ||x||^2 s.t. ||x - (3,0)|| >= 5 ; lambda* = 2/5, J* = 4 (nonsingular).
    p0, q0, r0, p1, q1, r1 = bracket.keepout_qcqp((F(3), F(0)), F(5))
    lam_star, j_star = F(2, 5), F(4)
    pts = []
    for k in range(1, 7):
        d = F(1, 10 ** k)
        lam = lam_star - d
        g = bracket.dual_value(p0, q0, r0, p1, q1, r1, lam)
        assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, lam, g)
        pts.append((lam, d, j_star - g))
    return pts


def hard_series():
    # min -x1^2 + 2 x2^2 + 2 x2 s.t. 1 - ||x||^2 >= 0 ; lambda* = 1, J* = -4/3
    # (trust-region hard case: A(lambda*) singular).
    p0 = [[F(-1), F(0)], [F(0), F(2)]]
    q0, r0 = [F(0), F(2)], F(0)
    p1 = [[F(-1), F(0)], [F(0), F(-1)]]
    q1, r1 = [F(0), F(0)], F(1)
    lam_star, j_star = F(1), F(-4, 3)
    pts = []
    for k in range(1, 7):
        d = F(1, 10 ** k)
        lam = lam_star + d                       # approach from A > 0 side
        g = bracket.dual_value(p0, q0, r0, p1, q1, r1, lam)
        assert bracket.certify_lower_bound(p0, q0, r0, p1, q1, r1, lam, g)
        pts.append((lam, d, j_star - g))
    return pts


def main():
    smooth, hard = smooth_series(), hard_series()

    print("SMOOTH (nonsingular, lambda*=2/5, J*=4):")
    for lam, d, gap in smooth:
        print(f"  lam={str(lam):>12}  |lam-lam*|={str(d):>10}  gap={gap}")
    print("HARD (trust-region, lambda*=1, J*=-4/3):")
    for lam, d, gap in hard:
        print(f"  lam={str(lam):>12}  |lam-lam*|={str(d):>10}  gap={gap}")

    out = sys.argv[1] if len(sys.argv) > 1 else None
    if out is not None:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.7})
    fig, ax = plt.subplots(figsize=(5.4, 3.4))

    xs = [float(d) for _, d, _ in smooth]
    ys_s = [float(gap) for _, _, gap in smooth]
    ys_h = [float(gap) for _, _, gap in hard]

    ax.loglog(xs, ys_s, color="black", linestyle="-", marker="o",
              markerfacecolor="white", markeredgecolor="black", markersize=6,
              linewidth=1.2, label="nonsingular (slope 2)")
    ax.loglog(xs, ys_h, color="black", linestyle="--", marker="s",
              markerfacecolor="black", markeredgecolor="black", markersize=5,
              linewidth=1.2, label="trust-region hard case (slope 1)")

    # slope reference guides (dotted / dash-dot), anchored for visibility
    gx = [1e-6, 1e-1]
    ax.loglog(gx, [40 * x ** 2 for x in gx], color="black", linestyle=":",
              linewidth=0.9)
    ax.loglog(gx, [0.9 * x for x in gx], color="black", linestyle="-.",
              linewidth=0.9)
    ax.text(2e-6, 40 * (2e-6) ** 2 * 1.6, "slope 2", fontsize=8)
    ax.text(2e-6, 0.9 * 2e-6 * 1.6, "slope 1", fontsize=8)

    ax.set_xlabel(r"dual-multiplier error $|\lambda - \lambda^\star|$")
    ax.set_ylabel(r"certified gap $J^\star - g(\lambda)$")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, color="0.7")
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    fig.tight_layout()

    if out is not None:
        fig.savefig(out)
        print(f"\nSaved figure to {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
