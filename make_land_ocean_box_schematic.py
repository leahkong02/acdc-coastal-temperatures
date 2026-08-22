"""Generate a compact schematic of the minimal land--ocean box model."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "figures"

TEXT = "#1F2D36"
MUTED = "#586A75"
LAND = "#A65A3A"
LAND_FILL = "#F3E3D7"
OCEAN = "#347DA0"
OCEAN_FILL = "#DDECF4"
ATM_FILL = "#EAF0F3"
FT_FILL = "#F2F2EF"


def add_box(
    ax,
    xy,
    width,
    height,
    title,
    state,
    facecolor,
    edgecolor,
    linestyle="-",
):
    """Add a labeled rounded model compartment."""
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.025,rounding_size=0.06",
        linewidth=1.45,
        edgecolor=edgecolor,
        facecolor=facecolor,
        linestyle=linestyle,
        zorder=2,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + 0.68 * height,
        title,
        ha="center",
        va="center",
        fontsize=10.2,
        fontweight="semibold",
        color=TEXT,
        zorder=3,
    )
    ax.text(
        x + width / 2,
        y + 0.31 * height,
        state,
        ha="center",
        va="center",
        fontsize=11.2,
        color=TEXT,
        zorder=3,
    )


def add_arrow(
    ax,
    start,
    end,
    color,
    label=None,
    label_xy=None,
    arrowstyle="-|>",
    linewidth=1.35,
):
    """Add a process arrow and an optional direct label."""
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle=arrowstyle,
        mutation_scale=11,
        linewidth=linewidth,
        color=color,
        shrinkA=2,
        shrinkB=2,
        connectionstyle="arc3,rad=0",
        zorder=4,
    )
    ax.add_patch(arrow)
    if label is not None and label_xy is not None:
        ax.text(
            label_xy[0],
            label_xy[1],
            label,
            ha="center",
            va="center",
            fontsize=9.2,
            color=color,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8},
            zorder=5,
        )


def main():
    """Render PNG, SVG, and PDF versions of the schematic."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.unicode_minus": False,
            "svg.fonttype": "none",
        }
    )

    fig, ax = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 7.2)
    ax.set_ylim(0, 5.0)
    ax.axis("off")

    ax.text(
        3.6,
        4.86,
        "Minimal land–ocean box model",
        ha="center",
        va="top",
        fontsize=14.0,
        fontweight="semibold",
        color=TEXT,
    )

    add_box(
        ax,
        (1.35, 3.85),
        4.50,
        0.62,
        "Free troposphere",
        r"$\theta_{FT}(t),\ q_{FT}$  (prescribed)",
        FT_FILL,
        MUTED,
        linestyle="--",
    )
    add_box(
        ax,
        (0.62, 2.08),
        2.42,
        1.00,
        "Land atmosphere",
        r"$\theta_L,\ q_L$",
        ATM_FILL,
        LAND,
    )
    add_box(
        ax,
        (4.16, 2.08),
        2.42,
        1.00,
        "Ocean atmosphere",
        r"$\theta_o,\ q_o$",
        ATM_FILL,
        OCEAN,
    )
    add_box(
        ax,
        (0.62, 0.66),
        2.42,
        0.78,
        "Land surface",
        r"$T_s,\ m$",
        LAND_FILL,
        LAND,
    )
    add_box(
        ax,
        (4.16, 0.66),
        2.42,
        0.78,
        "Ocean surface",
        r"$T_o$  (prescribed)",
        OCEAN_FILL,
        OCEAN,
        linestyle="--",
    )

    add_arrow(
        ax,
        (3.08, 2.58),
        (4.12, 2.58),
        MUTED,
        label=r"mixing  $\tau_{\rm mix}$",
        label_xy=(3.60, 2.82),
        arrowstyle="<|-|>",
        linewidth=1.45,
    )

    add_arrow(
        ax,
        (2.55, 3.83),
        (2.08, 3.11),
        MUTED,
        label=r"$W_{eL}$",
        label_xy=(2.18, 3.52),
    )
    add_arrow(
        ax,
        (4.65, 3.83),
        (5.12, 3.11),
        MUTED,
        label=r"$W_{eo}$",
        label_xy=(5.02, 3.52),
    )

    add_arrow(
        ax,
        (1.48, 1.47),
        (1.48, 2.05),
        LAND,
        label=r"$H_L$",
        label_xy=(1.28, 1.76),
    )
    add_arrow(
        ax,
        (2.18, 1.47),
        (2.18, 2.05),
        LAND,
        label=r"$E_L$",
        label_xy=(2.38, 1.76),
    )
    add_arrow(
        ax,
        (5.02, 1.47),
        (5.02, 2.05),
        OCEAN,
        label=r"$H_o$",
        label_xy=(4.82, 1.76),
    )
    add_arrow(
        ax,
        (5.72, 1.47),
        (5.72, 2.05),
        OCEAN,
        label=r"$E_o$",
        label_xy=(5.92, 1.76),
    )

    add_arrow(
        ax,
        (0.08, 1.18),
        (0.59, 1.18),
        LAND,
        label=r"$F(t)$",
        label_xy=(0.31, 1.38),
    )
    add_arrow(
        ax,
        (0.08, 0.87),
        (0.59, 0.87),
        LAND,
        label=r"$P(t)$",
        label_xy=(0.31, 0.67),
    )

    ax.text(
        3.60,
        0.18,
        r"Prognostic: $T_s,m,\theta_L,q_L,\theta_o,q_o$     "
        r"Prescribed: $T_o,\theta_{FT},q_{FT}$",
        ha="center",
        va="center",
        fontsize=8.8,
        color=MUTED,
    )

    outputs = {
        "png": OUTPUT_DIR / "land_ocean_box_model_schematic.png",
        "svg": OUTPUT_DIR / "land_ocean_box_model_schematic.svg",
        "pdf": OUTPUT_DIR / "land_ocean_box_model_schematic.pdf",
    }
    fig.savefig(outputs["png"], dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs["svg"], bbox_inches="tight", facecolor="white")
    fig.savefig(outputs["pdf"], bbox_inches="tight", facecolor="white")
    plt.close(fig)

    svg_text = outputs["svg"].read_text(encoding="utf-8")
    clean_svg = "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n"
    outputs["svg"].write_text(clean_svg, encoding="utf-8")

    for output in outputs.values():
        print(output)


if __name__ == "__main__":
    main()
