"""Schematics of the Lie-factor architecture (print figures for the manuals)."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent
NAVY = "#1f3a5f"
RUST = "#b85c38"
SLATE = "#5a6a7a"
INK = "#1a1a1a"
PAPER = "#f4f1ea"
CREAM = "#fffdf8"
CHECK = "#6e6e6e"
BLUEF = "#e7eef6"
RUSTF = "#f6ebe4"


def box(ax, x, y, w, h, fc, ec=NAVY, lw=1.2, r=0.08, ls="-"):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.015,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, linestyle=ls, zorder=2,
    ))


def arrow(ax, a, b, color=NAVY, lw=1.35):
    ax.add_patch(FancyArrowPatch(
        a, b, arrowstyle="-|>", mutation_scale=11,
        linewidth=lw, color=color, zorder=3,
    ))


def T(ax, x, y, s, size=10, color=INK, weight="normal", ha="center", va="center"):
    ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va,
            fontweight=weight, zorder=4)


def elbow(ax, x0, y0, x1, y1, color=NAVY, lw=1.2):
    """Vertical then horizontal then vertical, ending with an arrowhead."""
    mid = 0.5 * (y0 + y1)
    ax.plot([x0, x0], [y0, mid], color=color, lw=lw, zorder=3, solid_capstyle="round")
    ax.plot([x0, x1], [mid, mid], color=color, lw=lw, zorder=3, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch(
        (x1, mid), (x1, y1), arrowstyle="-|>", mutation_scale=11,
        linewidth=lw, color=color, zorder=3,
    ))


def fig_object():
    fig, ax = plt.subplots(figsize=(7.5, 4.55), dpi=200)
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 7.35)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    T(ax, 6.1, 7.12, "The ring is one Lie factor", size=13, weight="medium")

    hs = [(0.35, r"$H_{\mathrm{drift}}$"), (1.5, r"$H_{\mathrm{quad}}$"),
          (2.65, r"$H_{\mathrm{bend}}$"), (3.8, r"$H_{\mathrm{sext}}$")]
    for x, lab in hs:
        box(ax, x, 5.55, 1.05, 0.85, PAPER)
        T(ax, x + 0.525, 5.97, lab, size=10)
    T(ax, 2.6, 5.25, r"same $H$, expanded to a polynomial", size=8, color=SLATE)

    arrow(ax, (4.95, 5.97), (5.45, 5.97))

    box(ax, 5.45, 5.4, 2.05, 1.15, BLUEF)
    T(ax, 6.475, 6.15, "product", size=10)
    T(ax, 6.475, 5.8, "of factors", size=10)
    T(ax, 6.475, 5.52, r"$R,\;f_3,\;f_4$", size=8, color=SLATE)

    arrow(ax, (7.55, 5.97), (8.05, 5.97))

    box(ax, 8.05, 4.80, 3.8, 2.15, CREAM, lw=1.5)
    T(ax, 9.95, 6.72, r"one-turn factor $\mathcal{F}$", size=10, weight="medium")
    layers = [
        (6.10, r"$R$", "linear", BLUEF, NAVY),
        (5.50, r"$e^{:f_3:}$", "cubic", RUSTF, RUST),
        (4.90, r"$e^{:f_4:}$", "applied first", PAPER, SLATE),
    ]
    for y, title, sub, fc, ec in layers:
        box(ax, 8.25, y, 3.4, 0.50, fc, ec=ec, r=0.06)
        T(ax, 9.35, y + 0.25, title, size=11)
        T(ax, 10.85, y + 0.25, sub, size=8, color=SLATE)

    T(ax, 6.1, 4.42, r"three uses of the same $\mathcal{F}$", size=10, color=SLATE)

    box(ax, 0.35, 1.45, 3.55, 1.55, BLUEF, ec=NAVY)
    T(ax, 2.125, 2.62, "read", size=11, weight="medium")
    T(ax, 2.125, 2.18, r"Twiss, $D$, $\eta$, $\xi$", size=10)
    T(ax, 2.125, 1.78, "observables(lat)", size=8, color=SLATE)

    box(ax, 4.3, 1.45, 3.55, 1.55, RUSTF, ec=RUST)
    T(ax, 6.075, 2.62, "track", size=11, weight="medium")
    T(ax, 6.075, 2.18, "flow of each factor,", size=10)
    T(ax, 6.075, 1.78, "or unexpanded $H$", size=10)

    box(ax, 8.25, 1.45, 3.6, 1.55, "#f0f0f0", ec=CHECK, ls="--", lw=1.05)
    T(ax, 10.05, 2.62, "check", size=11, weight="medium", color=CHECK)
    T(ax, 10.05, 2.18, "RK4 polynomial map", size=9, color=CHECK)
    T(ax, 10.05, 1.78, "not the ring object", size=8, color=CHECK)

    elbow(ax, 9.95, 4.80, 2.125, 3.08)
    elbow(ax, 9.95, 4.80, 6.075, 3.08)
    elbow(ax, 9.95, 4.80, 10.05, 3.08, color=CHECK, lw=1.05)

    T(ax, 6.1, 0.55, "MaryLie: concatenate in the algebra; analyse the product; track a point through it.",
      size=8, color=SLATE)

    fig.tight_layout(pad=0.2)
    fig.savefig(OUT / "architecture_object.png", dpi=200, bbox_inches="tight",
                pad_inches=0.18, facecolor="white")
    plt.close(fig)


def fig_orders():
    fig, ax = plt.subplots(figsize=(7.5, 3.55), dpi=200)
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 5.35)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor("white")

    T(ax, 6.1, 5.0, r"One polynomial, three homogeneous readings", size=13, weight="medium")
    T(ax, 6.1, 4.55, r"$M=\mathrm{factor\_to\_map}(\mathcal{F})$  —  not a stack of independent theories",
      size=9, color=SLATE)

    rows = [
        (3.15, BLUEF, NAVY, "degree 1", r"$R$",
         r"Twiss, $D_1$, $\eta_0$",
         r"$R=\exp(L\,S\,\mathrm{Hess}\,H)$ is $\exp(:f_2:)$ closed"),
        (1.75, RUSTF, RUST, "degree 2", r"$M^{[2]}$",
         r"$D_2$, $\eta_1$",
         r"a $6\times 6$ matrix cannot see this"),
        (0.35, PAPER, SLATE, "degree 3", r"$M^{[3]}$",
         r"$\xi_1=\mathrm{d}Q/\mathrm{d}\delta$",
         r"Jacobian of $M$ on the off-momentum orbit"),
    ]
    for y, fc, ec, deg, left, mid, note in rows:
        box(ax, 0.35, y, 11.5, 1.22, fc, ec=ec, lw=1.25)
        T(ax, 1.55, y + 0.61, deg, size=9, color=ec)
        T(ax, 3.85, y + 0.61, left, size=13)
        T(ax, 7.15, y + 0.78, mid, size=11)
        T(ax, 7.15, y + 0.38, note, size=8, color=SLATE)

    fig.tight_layout(pad=0.2)
    fig.savefig(OUT / "architecture_orders.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    fig_object()
    fig_orders()
    print("wrote", OUT / "architecture_object.png")
    print("wrote", OUT / "architecture_orders.png")
