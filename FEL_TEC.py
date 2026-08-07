"""
NS5 Free Energy Landscape — Apo vs Holo (NS5+TEC)
Full pipeline: Trajectory → PCA → FEL → Figure
pip install MDAnalysis scikit-learn matplotlib numpy scipy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sklearn.decomposition import PCA
import MDAnalysis as mda
from MDAnalysis.analysis import align
import warnings
warnings.filterwarnings("ignore")

# ── PATHS ─────────────────────────────────────────────────────────────────────
APO_TPR  = "/Users/kanokpant.sriwong/Desktop/md_dock/apo_md/npt.gro"
APO_XTC  = "/Users/kanokpant.sriwong/Desktop/md_dock/apo_md/md.xtc"
HOLO_TPR = "/Users/kanokpant.sriwong/Desktop/md_dock/npt.gro"
HOLO_XTC = "/Users/kanokpant.sriwong/Desktop/md_dock/md.xtc"
OUTFILE  = "/Users/kanokpant.sriwong/Desktop/md_dock/FEL_publication.png"
# ─────────────────────────────────────────────────────────────────────────────

kBT = 2.479  # kJ/mol at 300 K

# Rainbow colormap: blue (stable) → red (unstable)
RAINBOW = LinearSegmentedColormap.from_list("fel", [
    "#0000FF",   # blue   — minimum energy
    "#00AAFF",   # cyan
    "#00FF88",   # green
    "#AAFF00",   # yellow-green
    "#FFFF00",   # yellow
    "#FF8800",   # orange
    "#FF0000",   # red    — maximum energy
], N=512)

plt.rcParams.update({
    "font.family":   "Arial",
    "font.size":     10,
    "axes.linewidth": 1.0,
    "pdf.fonttype":  42,
    "svg.fonttype":  "none",
})

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Load trajectory and run PCA
# ─────────────────────────────────────────────────────────────────────────────

def load_and_pca(tpr, xtc, label=""):
    print(f"\n{'─'*50}")
    print(f"  Loading {label} trajectory...")
    u = mda.Universe(tpr, xtc)
    backbone = u.select_atoms("backbone")
    print(f"  Atoms selected (backbone): {len(backbone)}")
    print(f"  Total frames: {len(u.trajectory)}")

    # Align all frames to first frame
    print(f"  Aligning trajectory...")
    ref = mda.Universe(tpr, xtc)
    aligner = align.AlignTraj(u, ref,
                               select="backbone",
                               in_memory=True).run()

    # Collect coordinates for all frames
    print(f"  Collecting coordinates...")
    coords = []
    for ts in u.trajectory:
        coords.append(backbone.positions.flatten())
    coords = np.array(coords)
    print(f"  Coordinate matrix shape: {coords.shape}")

    # PCA
    print(f"  Running PCA...")
    pca = PCA(n_components=2)
    projected = pca.fit_transform(coords)
    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100
    print(f"  PC1 variance explained: {var1:.1f}%")
    print(f"  PC2 variance explained: {var2:.1f}%")

    return projected[:, 0], projected[:, 1], var1, var2


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Compute Gibbs FEL
# ─────────────────────────────────────────────────────────────────────────────

def compute_fel(pc1, pc2, bins=60, smooth=1.8):
    """
    ΔG = -kBT ln P(PC1, PC2)
    Normalized so global minimum = 0
    """
    H, xe, ye = np.histogram2d(pc1, pc2, bins=bins, density=True)
    H = gaussian_filter(H.astype(float), sigma=smooth)
    H[H <= 0] = np.nan
    G = -kBT * np.log(H / np.nanmax(H))
    G[G > 8.5] = np.nan
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    return X, Y, G


def get_minimum(X, Y, G):
    idx = np.unravel_index(np.nanargmin(G), G.shape)
    return X[idx], Y[idx], G[idx]


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Plotting functions
# ─────────────────────────────────────────────────────────────────────────────

def plot_3d(ax, X, Y, G, title, norm, label, var1, var2, fig):
    Gf = np.nan_to_num(G, nan=float(np.nanmax(G[~np.isnan(G)])))

    # Surface
    ax.plot_surface(X, Y, Gf,
                    facecolors=RAINBOW(norm(Gf)),
                    alpha=0.92, linewidth=0,
                    antialiased=True, shade=False)

    # Contour projection at bottom
    z_floor = -0.5
    ax.contourf(X, Y, Gf,
                zdir="z", offset=z_floor,
                levels=40, cmap=RAINBOW, norm=norm, alpha=0.85)

    # Axes
    ax.set_zlim(z_floor, 8.5)
    ax.set_xlabel(f"PC1 ({var1:.1f}%)", fontsize=8, labelpad=8)
    ax.set_ylabel(f"PC2 ({var2:.1f}%)", fontsize=8, labelpad=8)
    ax.set_zlabel("G (kJ/mol)", fontsize=8, labelpad=8)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=14)
    ax.tick_params(labelsize=6)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor("#CCCCCC")
    ax.yaxis.pane.set_edgecolor("#CCCCCC")
    ax.zaxis.pane.set_edgecolor("#CCCCCC")
    ax.grid(True, color="#DDDDDD", linewidth=0.4)
    ax.view_init(elev=30, azim=-55)

    # Panel label
    ax.text2D(-0.06, 1.04, label, transform=ax.transAxes,
              fontsize=14, fontweight="bold", va="top")

    # Individual colorbar for this 3D panel
    sm = plt.cm.ScalarMappable(cmap=RAINBOW, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.1, shrink=0.6)
    cb.set_label("G (kJ/mol)", fontsize=8)
    cb.ax.tick_params(labelsize=7)


def plot_2d(ax, X, Y, G, title, norm, label, var1, var2, fig):
    Gf = np.nan_to_num(G, nan=float(np.nanmax(G[~np.isnan(G)])))

    cf = ax.contourf(X, Y, Gf, levels=60,
                     cmap=RAINBOW, norm=norm, alpha=0.95)
    ax.contour(X, Y, Gf, levels=15,
               colors="white", linewidths=0.3, alpha=0.3)

    # Individual colorbar for this 2D panel
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("G (kJ/mol)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    ax.set_xlabel(f"PC1 ({var1:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({var2:.1f}%)", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel label
    ax.text(-0.10, 1.06, label, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="top")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    # ── PCA ──────────────────────────────────────────────────────────────────
    apc1, apc2, avar1, avar2 = load_and_pca(APO_TPR,  APO_XTC,  "Apo-NS5")
    hpc1, hpc2, hvar1, hvar2 = load_and_pca(HOLO_TPR, HOLO_XTC, "Holo NS5+TEC")

    # ── FEL ──────────────────────────────────────────────────────────────────
    print("\nComputing FEL...")
    Xa, Ya, Ga = compute_fel(apc1, apc2)
    Xh, Yh, Gh = compute_fel(hpc1, hpc2)

    xma, yma, zma = get_minimum(Xa, Ya, Ga)
    xmh, ymh, zmh = get_minimum(Xh, Yh, Gh)
    print(f"  Apo  minimum: PC1={xma:.2f}  PC2={yma:.2f}  ΔG={zma:.3f} kJ/mol")
    print(f"  Holo minimum: PC1={xmh:.2f}  PC2={ymh:.2f}  ΔG={zmh:.3f} kJ/mol")

    norm = Normalize(vmin=0, vmax=8)

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 12))
    fig.patch.set_facecolor("white")

    fig.suptitle(
        "Free Energy Landscape of NS5  —  TEC",
        fontsize=11, fontweight="bold", y=0.995
    )

    gs = gridspec.GridSpec(
        2, 2,
        left=0.07, right=0.96,
        top=0.93, bottom=0.05,
        hspace=0.38, wspace=0.32
    )

    # ── Row 0: 3D panels ─────────────────────────────────────────────────────
    ax3a = fig.add_subplot(gs[0, 0], projection="3d")
    ax3h = fig.add_subplot(gs[0, 1], projection="3d")

    plot_3d(ax3a, Xa, Ya, Ga,
            title="Apo structure - NS5",
            norm=norm, label="A",
            var1=avar1, var2=avar2, fig=fig)

    plot_3d(ax3h, Xh, Yh, Gh,
            title="Holo structure - NS5_TEC complex",
            norm=norm, label="B",
            var1=hvar1, var2=hvar2, fig=fig)

    # ── Row 1: 2D panels ─────────────────────────────────────────────────────
    ax2a = fig.add_subplot(gs[1, 0])
    ax2h = fig.add_subplot(gs[1, 1])

    plot_2d(ax2a, Xa, Ya, Ga,
            title="Apo structure - NS5",
            norm=norm, label="C",
            var1=avar1, var2=avar2, fig=fig)

    plot_2d(ax2h, Xh, Yh, Gh,
            title="Holo structure - NS5_TEC complex",
            norm=norm, label="D",
            var1=hvar1, var2=hvar2, fig=fig)

    # ── Footer ───────────────────────────────────────────────────────────────
    fig.text(
        0.5, 0.005,
        "",
        ha="center", fontsize=7.5, color="#666666"
    )

    # ── Save ─────────────────────────────────────────────────────────────────
    for fmt, dpi in [("png", 300)]:
        out = OUTFILE.replace(".png", f".{fmt}")
        kw = dict(bbox_inches="tight", facecolor="white")
        if dpi:
            kw["dpi"] = dpi
        plt.savefig(out, **kw)
        print(f"✅ Saved: {out}")

    plt.show()


if __name__ == "__main__":
    main()


