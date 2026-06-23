"""
Per-gene-family expression heatmaps — Supplementary Figs. S17–S23.

Layout
------
  Columns : Ant_Vm | Ant_VF | Ant_MF ‖ Palp_Vm | Palp_VF | Palp_MF ‖ Tarsi_Vm | Tarsi_VF | Tarsi_MF
  Rows    : genes grouped by appendage of peak expression (Ant → Palp → Tarsi),
            sorted by mean expression within each group (descending)
  Colour  : log2(counts + 1) centred at threshold = 10 norm counts; global scale across all families
  Page    : hard A4 (210 × 297 mm); content within 180 × 260 mm margins
  Split   : families with >80 genes rendered in two side-by-side columns
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches   # used for legend patches
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE       = pathlib.Path(__file__).parent / "Submission_ready" / "app" / "data"
OUT_APP    = BASE / "figures" / "supplementary"
OUT_SUB    = pathlib.Path(__file__).parent / "Submission_ready" / "Supplementary_Figures"
COUNTS_FILE = BASE / "BSF_normalized_counts_nameX.csv"
# BSF_Olfactory_ID.csv — used to add TRP gene names (and any others missing from Name col)
OLFACTORY_ID_FILE = pathlib.Path(__file__).parent / "BSF2026" / "Supplemetary_Files" / "csv" / "BSF_Olfactory_ID.csv"
for p in [OUT_APP, OUT_SUB]:
    p.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
THRESHOLD  = 10.0          # normalised-count expression threshold (white on heatmap)
THR_LOG    = float(np.log2(THRESHOLD + 1.0))   # = log2(11) ≈ 3.459
SPLIT_AT   = 80            # families with more genes → 2-column layout

FAM_FIGNUMS = {"OR": "S17", "GR": "S18", "IR": "S19",
               "OBP": "S20", "PPK": "S21", "CSP": "S22", "TRP": "S23"}
FAMILIES = list(FAM_FIGNUMS.keys())

# A4 exact size
A4_W, A4_H = 8.268, 11.693          # inches (210 × 297 mm)
# Margins as fractions of figure size  (15 mm L/R, 18 mm top, 19 mm bottom)
MAR_L, MAR_R = 15/210, 15/210
MAR_T, MAR_B = 18/297, 19/297

# Colour map: blue – white – red
CMAP = LinearSegmentedColormap.from_list("bwr", ["#1f4aa8", "#ffffff", "#b2182b"])

# Column spec: (data prefix, display label, label colour)
COLS = [
    ("Ant_Vm",  "Ant\nVm",    "#4D9ACD"),
    ("Ant_VF",  "Ant\nVF",    "#E878A0"),
    ("Ant_MF",  "Ant\nMF",    "#7B1FA2"),
    ("P_Vm",    "Palp\nVm",   "#4D9ACD"),
    ("P_VF",    "Palp\nVF",   "#E878A0"),
    ("P_MF",    "Palp\nMF",   "#7B1FA2"),
    ("Leg_Vm",  "Tarsi\nVm",  "#4D9ACD"),
    ("Leg_VF",  "Tarsi\nVF",  "#E878A0"),
    ("Leg_MF",  "Tarsi\nMF",  "#7B1FA2"),
]
N_COLS = len(COLS)   # 9


# ── Load & prepare data ────────────────────────────────────────────────────────
print("Loading normalised counts …")
raw = pd.read_csv(COUNTS_FILE)

# Enrich Name column using BSF_Olfactory_ID.csv (fills TRP and any others
# whose names were not already in the count matrix Name column)
if OLFACTORY_ID_FILE.exists():
    olf = pd.read_csv(OLFACTORY_ID_FILE)
    olf["XM"] = olf["Transcript"].str.replace(r"^rna-", "", regex=True)
    # Build mapping XM_ → gene name, only for entries not already named
    xm_to_name = dict(zip(olf["XM"], olf["Name"]))
    mask_unnamed = raw["Name"].isna()
    raw.loc[mask_unnamed, "Name"] = raw.loc[mask_unnamed, "Gene"].map(xm_to_name)
    n_filled = raw["Name"].notna().sum() - (~mask_unnamed).sum()
    print(f"  Filled {n_filled} gene names from BSF_Olfactory_ID.csv")
else:
    print(f"  WARNING: {OLFACTORY_ID_FILE} not found — TRP family will be empty")

for prefix, _, _ in COLS:
    rep_cols = [c for c in raw.columns if c.startswith(prefix) and c[-1].isdigit()]
    raw[f"{prefix}_mean"] = raw[rep_cols].apply(pd.to_numeric, errors="coerce").mean(axis=1).fillna(0)

MEAN_COLS = [f"{p}_mean" for p, _, _ in COLS]

def assign_fam(name):
    n = str(name).upper()
    # Check TRP first (before OR/GR/IR) to avoid false positives from substring matches
    if n.startswith("TRP"):
        return "TRP"
    for tag in ["OBP", "PPK", "CSP", "GR", "IR", "OR"]:   # OBP before OR/IR
        if tag in n:
            return tag
    return None

raw["gene_fam"] = raw["Name"].apply(assign_fam)

# ── Compute GLOBAL colour scale (same for all families) ───────────────────────
# Anchored to the original 6-family set (S17-S22) so that adding TRP (S23) does
# not shift the shared colour norm retroactively.  TRP uses this same scale.
ANCHOR_FAMILIES = ["OR", "GR", "IR", "OBP", "PPK", "CSP"]
all_centered = []
for fam in ANCHOR_FAMILIES:
    sub = raw[raw["gene_fam"] == fam]
    if sub.empty:
        continue
    mat = sub[MEAN_COLS].to_numpy(dtype=float)
    centered = np.log2(mat + 1.0) - THR_LOG
    all_centered.append(centered.ravel())

all_vals = np.concatenate(all_centered)
# Use 99th-percentile absolute value so a few extreme outliers don't crush the scale
GLOBAL_LIM = float(max(np.percentile(np.abs(all_vals), 99), 0.5))
NORM = TwoSlopeNorm(vmin=-GLOBAL_LIM, vcenter=0.0, vmax=GLOBAL_LIM)

GLOBAL_LIM = GLOBAL_LIM / 4        # user requested 4× tighter scale
print(f"Global colour limit: ±{GLOBAL_LIM:.4f}  (log₂(counts+1) − log₂(11), scaled to 99th-pct/4)")

# ── Helpers ───────────────────────────────────────────────────────────────────
def order_by_appendage(mat_centered):
    """Group rows by appendage of peak expression (Ant→Palp→Tarsi),
    sorted by mean expression within each group (descending)."""
    ant   = mat_centered[:, 0:3].mean(axis=1)
    palp  = mat_centered[:, 3:6].mean(axis=1)
    tarsi = mat_centered[:, 6:9].mean(axis=1)
    peak  = np.argmax(np.column_stack([ant, palp, tarsi]), axis=1)
    order = []
    for t in range(3):
        idx  = np.where(peak == t)[0]
        vals = np.column_stack([ant, palp, tarsi])[idx, t]
        order.extend(idx[np.argsort(vals)[::-1]])
    assigned = set(order)
    order += [i for i in range(len(mat_centered)) if i not in assigned]
    return np.array(order, dtype=int)


def font_pt(n_rows, axes_h_mm):
    """Font size (pt) that avoids row-label overlap given axes height in mm."""
    mm_per_row = axes_h_mm / max(n_rows, 1)
    pt = mm_per_row * 2.835 * 0.62
    return float(np.clip(pt, 3.5, 9.0))



def draw_panel(ax, mat, names, col_labels, col_colors,
               ytick_side="left", fsize=7.0):
    """Render one heatmap panel; return the AxesImage for the shared colourbar."""
    im = ax.imshow(mat, aspect="auto", cmap=CMAP, norm=NORM, interpolation="nearest")

    ax.set_xticks(np.arange(N_COLS))
    ax.set_xticklabels(col_labels, fontsize=7, ha="center")
    for tick, col in zip(ax.get_xticklabels(), col_colors):
        tick.set_color(col)

    ax.set_yticks(np.arange(len(names)))
    if ytick_side == "right":
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
    ax.set_yticklabels(names, fontsize=fsize, va="center")
    ax.tick_params(axis="y", length=0, pad=2)

    # White separators between tissue groups
    for xpos in [2.5, 5.5]:
        ax.axvline(xpos, color="white", linewidth=1.5, zorder=3)

    return im


# ── Per-family loop ────────────────────────────────────────────────────────────
col_labels = [lbl for _, lbl, _ in COLS]
col_colors = [col for _, _, col in COLS]

for fam in FAMILIES:
    fignum = FAM_FIGNUMS[fam]
    sub = raw[raw["gene_fam"] == fam].copy().reset_index(drop=True)
    if sub.empty:
        print(f"  {fam}: no genes, skipping.")
        continue

    n_genes  = len(sub)
    mat_raw  = sub[MEAN_COLS].to_numpy(dtype=float)
    mat_cent = np.log2(mat_raw + 1.0) - THR_LOG      # centred at threshold

    gene_order    = order_by_appendage(mat_cent)
    mat_ordered   = mat_cent[gene_order]
    names_ordered = sub["Name"].iloc[gene_order].astype(str).tolist()

    two_col = n_genes > SPLIT_AT
    half    = (n_genes + 1) // 2 if two_col else n_genes

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(A4_W, A4_H))

    cx0, cy0 = MAR_L, MAR_B
    cx1, cy1 = 1 - MAR_R, 1 - MAR_T
    cw, ch   = cx1 - cx0, cy1 - cy0

    # Vertical budget (fractions of ch)
    title_frac  = 0.050
    xlabel_frac = 0.065
    legend_frac = 0.040
    hm_frac     = 1 - title_frac - xlabel_frac - legend_frac

    hm_y0   = cy0 + (legend_frac + xlabel_frac) * ch
    hm_y1   = cy1 - title_frac * ch
    hm_h_mm = hm_frac * ch * A4_H * 25.4

    if two_col:
        # Both panels have gene names on the LEFT.
        # [lbl: left labels][left hm][gap wide enough for right labels][right hm][cbar]
        lbl_f  = 0.20;  hm_f = 0.25;  gap_f = 0.22;  cbar_f = 0.03

        lhm_x0 = cx0 + lbl_f * cw
        lhm_x1 = lhm_x0 + hm_f * cw
        rhm_x0 = lhm_x1 + gap_f * cw
        rhm_x1 = rhm_x0 + hm_f * cw
        cbar_x  = rhm_x1 + 0.02 * cw

        ax_l = fig.add_axes([lhm_x0, hm_y0, hm_f * cw, hm_y1 - hm_y0])
        ax_r = fig.add_axes([rhm_x0, hm_y0, hm_f * cw, hm_y1 - hm_y0])

        n_l, n_r = half, n_genes - half
        fl = font_pt(n_l, hm_h_mm)
        fr = font_pt(n_r, hm_h_mm)

        im = draw_panel(ax_l, mat_ordered[:n_l], names_ordered[:n_l],
                        col_labels, col_colors, ytick_side="left", fsize=fl)
        draw_panel(ax_r, mat_ordered[n_l:], names_ordered[n_l:],
                   col_labels, col_colors, ytick_side="left", fsize=fr)

        cbar_ax = fig.add_axes([cbar_x, hm_y0 + 0.10*(hm_y1-hm_y0),
                                 cbar_f * cw, 0.80*(hm_y1-hm_y0)])

    else:
        # Single column: left_labels | heatmap | gap | cbar | right_margin
        lbl_f  = 0.36;  hm_f = 0.52;  cbar_gap = 0.02;  cbar_f = 0.03

        lhm_x0 = cx0 + lbl_f * cw
        lhm_w  = hm_f * cw
        cbar_x  = lhm_x0 + lhm_w + cbar_gap * cw

        ax_l = fig.add_axes([lhm_x0, hm_y0, lhm_w, hm_y1 - hm_y0])
        fl   = font_pt(n_genes, hm_h_mm)

        im = draw_panel(ax_l, mat_ordered, names_ordered,
                        col_labels, col_colors, ytick_side="left", fsize=fl)

        cbar_ax = fig.add_axes([cbar_x, hm_y0 + 0.10*(hm_y1-hm_y0),
                                 cbar_f * cw, 0.80*(hm_y1-hm_y0)])

    # ── Colourbar ─────────────────────────────────────────────────────────────
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("log₂(counts+1) − log₂(11)\n[white = 10 counts threshold]",
                   fontsize=6.5, rotation=90, labelpad=4)
    # Show a few intuitive tick values
    tick_vals = [v for v in [-GLOBAL_LIM, -2, -1, 0, 1, 2, GLOBAL_LIM]
                 if -GLOBAL_LIM <= v <= GLOBAL_LIM]
    cbar.set_ticks(tick_vals)
    cbar.ax.set_yticklabels([f"{v:+.1f}" for v in tick_vals], fontsize=6)
    cbar.ax.axhline(0, color="black", linewidth=0.8)

    # ── Title ─────────────────────────────────────────────────────────────────
    title_y = cy1 - (title_frac / 2) * ch
    fig.text(0.5, title_y,
             f"Figure {fignum}  —  {fam} family expression heatmap   (n = {n_genes} genes)",
             ha="center", va="center", fontsize=9, fontweight="bold",
             transform=fig.transFigure)

    # ── State colour legend ────────────────────────────────────────────────────
    patches = [
        mpatches.Patch(color="#4D9ACD", label="Virgin male"),
        mpatches.Patch(color="#E878A0", label="Virgin female"),
        mpatches.Patch(color="#7B1FA2", label="Mated female"),
    ]
    fig.legend(handles=patches,
               loc="lower center", bbox_to_anchor=(0.5, cy0 - 0.005),
               ncol=3, fontsize=8, frameon=False,
               bbox_transform=fig.transFigure)

    # ── Save (exact A4, no bbox_inches adjustment) ─────────────────────────────
    for out_dir in [OUT_APP, OUT_SUB]:
        out_path = out_dir / f"Figure_{fignum}.pdf"
        fig.savefig(out_path, format="pdf", bbox_inches=None)
        print(f"    saved: {out_path}")

    plt.close(fig)
    layout = f"2-col ({half}+{n_genes-half})" if two_col else "1-col"
    print(f"  {fam:4s}  {n_genes:3d} genes  {layout}  font={fl:.1f}pt  → {fignum}")

print(f"\nDone. Global scale: ±{GLOBAL_LIM:.2f} log₂ units (centred at 10 counts).")
