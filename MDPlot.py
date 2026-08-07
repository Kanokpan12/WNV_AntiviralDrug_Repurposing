"""
NS5_TEC — MD Simulation Analysis
pip install matplotlib numpy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import os

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────────────────────
WORK_DIR = "/Users/kanokpant.sriwong/Desktop/md_dock"

FILES = {
    "rmsd_protein": "rmsd_protein.xvg",
    "rmsd_ligand":  "rmsd_ligand.xvg",
    "rmsf":         "rmsf.xvg",
    "hbond":        "hbond.xvg",
}

# Publication color palette (colorblind-friendly)
COLORS = {
    "protein": "#2166AC",
    "ligand":  "#D6604D",
    "rmsf":    "#4DAC26",
    "flex":    "#B2182B",
    "hbond":   "#762A83",
    "mean":    "#000000",
}

FONT = "Arial"
plt.rcParams.update({
    "font.family":       FONT,
    "font.size":         9,
    "axes.linewidth":    0.8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "xtick.direction":   "out",
    "ytick.direction":   "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "pdf.fonttype":      42,   # editable text in Illustrator
    "svg.fonttype":      "none",
})

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def parse_xvg(filepath):
    x, y = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            cols = line.split()
            if len(cols) >= 2:
                try:
                    x.append(float(cols[0]))
                    y.append(float(cols[1]))
                except ValueError:
                    continue
    return np.array(x), np.array(y)


def fp(key):
    return os.path.join(WORK_DIR, FILES[key])


def style_ax(ax, xlabel="", ylabel=""):
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(False)


def panel_label(ax, label):
    """Add bold panel label e.g. A, B, C, D"""
    ax.text(-0.12, 1.06, label, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="left")


def missing_panel(ax, key):
    ax.text(0.5, 0.5, f"{FILES[key]}\nnot found",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=8, color="#AAAAAA",
            bbox=dict(boxstyle="round", facecolor="#F5F5F5",
                      edgecolor="#CCCCCC", linewidth=0.8))
    ax.set_xticks([])
    ax.set_yticks([])


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    os.chdir(WORK_DIR)

    fig = plt.figure(figsize=(10, 8))
    fig.patch.set_facecolor("white")

    gs = gridspec.GridSpec(
        2, 2,
        left=0.10, right=0.97,
        top=0.93, bottom=0.09,
        hspace=0.52, wspace=0.38
    )

    # ── Figure title ─────────────────────────────────────────────────────────
    fig.suptitle(
        "MD Simulation of WNV NS5–TEC Complex",
        fontsize=11, fontweight="bold", y=0.99, x=0.53
    )

    # ── A: Protein RMSD ──────────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    panel_label(ax1, "A")
    if os.path.exists(fp("rmsd_protein")):
        t, r = parse_xvg(fp("rmsd_protein"))
        r_a = r * 10
        ax1.plot(t, r_a, color=COLORS["protein"], lw=0.9, alpha=0.85, rasterized=True)
        ax1.fill_between(t, r_a, alpha=0.10, color=COLORS["protein"])
        m = np.mean(r_a)
        ax1.axhline(m, color=COLORS["mean"], ls="--", lw=1.0,
                    label=f"Mean = {m:.2f} Å")
        ax1.legend(fontsize=8, frameon=False)
        ax1.set_ylim(bottom=0)
    else:
        missing_panel(ax1, "rmsd_protein")
    style_ax(ax1, xlabel="Time (ns)", ylabel="RMSD (Å)")
    ax1.set_title("Protein Backbone RMSD", fontsize=9, pad=4)

    # ── B: Ligand RMSD ───────────────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    panel_label(ax2, "B")
    if os.path.exists(fp("rmsd_ligand")):
        t, r = parse_xvg(fp("rmsd_ligand"))
        r_a = r * 10
        ax2.plot(t, r_a, color=COLORS["ligand"], lw=0.9, alpha=0.85, rasterized=True)
        ax2.fill_between(t, r_a, alpha=0.10, color=COLORS["ligand"])
        m = np.mean(r_a)
        ax2.axhline(m, color=COLORS["mean"], ls="--", lw=1.0,
                    label=f"Mean = {m:.2f} Å")
        ax2.legend(fontsize=8, frameon=False)
        ax2.set_ylim(bottom=0)
    else:
        missing_panel(ax2, "rmsd_ligand")
    style_ax(ax2, xlabel="Time (ns)", ylabel="RMSD (Å)")
    ax2.set_title("Ligand (TEC) RMSD", fontsize=9, pad=4)

    # ── C: RMSF ──────────────────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 0])
    panel_label(ax3, "C")
    if os.path.exists(fp("rmsf")):
        res, r = parse_xvg(fp("rmsf"))
        r_a = r * 10
        thresh = np.mean(r_a) + np.std(r_a)
        bar_colors = [COLORS["flex"] if v > thresh else COLORS["rmsf"] for v in r_a]
        ax3.bar(res, r_a, color=bar_colors, width=1,
                edgecolor="none", alpha=0.85)
        m = np.mean(r_a)
        ax3.axhline(m, color=COLORS["mean"], ls="--", lw=1.0,
                    label=f"Mean = {m:.2f} Å")
        ax3.legend(handles=[
            Patch(facecolor=COLORS["rmsf"],  alpha=0.85, label="Normal"),
            Patch(facecolor=COLORS["flex"],  alpha=0.85, label="Flexible"),
            plt.Line2D([0], [0], color=COLORS["mean"], ls="--",
                       lw=1.0, label=f"Mean = {m:.2f} Å"),
        ], fontsize=8, frameon=False)
        ax3.set_ylim(bottom=0)
    else:
        missing_panel(ax3, "rmsf")
    style_ax(ax3, xlabel="Residue Number", ylabel="RMSF (Å)")
    ax3.set_title("Per-Residue Flexibility (RMSF)", fontsize=9, pad=4)

    # ── D: H-bonds ───────────────────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 1])
    panel_label(ax4, "D")
    if os.path.exists(fp("hbond")):
        t, h = parse_xvg(fp("hbond"))
        ax4.plot(t, h, color=COLORS["hbond"], lw=0.9, alpha=0.85, rasterized=True)
        ax4.fill_between(t, h, alpha=0.10, color=COLORS["hbond"])
        m = np.mean(h)
        ax4.axhline(m, color=COLORS["mean"], ls="--", lw=1.0,
                    label=f"Mean = {m:.1f}")
        ax4.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax4.legend(fontsize=8, frameon=False, loc="upper left")
        ax4.set_ylim(bottom=0)
    else:
        missing_panel(ax4, "hbond")
        ax4.text(0.5, 0.18,
                 "Generate with:\necho '1 13' | gmx hbond -s md.tpr\n-f md.xtc -n index.ndx -num hbond.xvg",
                 ha="center", va="center", transform=ax4.transAxes,
                 fontsize=7, color="#888888", family="monospace")
    style_ax(ax4, xlabel="Time (ns)", ylabel="Number of H-bonds")
    ax4.set_title("Protein–Ligand Hydrogen Bonds", fontsize=9, pad=4)

    # ── Save ─────────────────────────────────────────────────────────────────
    for fmt, dpi in [("png", 300)]:
        out = os.path.join(WORK_DIR, f"md_publication.{fmt}")
        kwargs = dict(bbox_inches="tight", facecolor="white")
        if dpi:
            kwargs["dpi"] = dpi
        plt.savefig(out, **kwargs)
        print(f"✅ Saved: {out}")

    plt.show()


if __name__ == "__main__":
    main()


