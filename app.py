# =============================================================================
# app.py  –  BSF Transcriptome Explorer  (merged from app_t1, c1, s1, h1)
# =============================================================================

# === IMPORTS ===
import base64
from urllib.parse import quote
import os
import re
import math
import pathlib
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from PIL import Image
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from typing import List

# === PAGE CONFIG (must be first Streamlit call) ===
st.set_page_config(
    page_title="Appendage identity dominates chemosensory transcriptomic "
               "organization in the black soldier fly | Perets et al. 2026",
    page_icon="🪰",
    layout="wide",
)
st.markdown("# Appendage identity dominates chemosensory transcriptomic "
            "organization in the black soldier fly")
st.caption("Perets *et al.* 2026 · data browser")

st.markdown("""
<style>
/* Target the tab labels */
div[data-baseweb="tab-list"] button:nth-child(1) p::before {
  content: "● ";
  color: #01045A;
  font-weight: 900;
}
div[data-baseweb="tab-list"] button:nth-child(2) p::before {
  content: "● ";
  color: #3F6ADE;
  font-weight: 900;
}
div[data-baseweb="tab-list"] button:nth-child(3) p::before {
  content: "● ";
  color: #9142E0;
  font-weight: 900;
}
div[data-baseweb="tab-list"] button:nth-child(4) p::before {
  content: "● ";
  color: #03B5E2;
  font-weight: 900;
}
</style>
""", unsafe_allow_html=True)

# =============================================================================
# === SHARED CONSTANTS ===
# =============================================================================

BASE_DIR = str(pathlib.Path(__file__).resolve().parent / "data")

# GO domain colors (identical across all files)
GO_COLS = {"MF": "#3FC498", "CC": "#4D50DB", "BP": "#A860E3", "Unknown": "#666666"}
GO_COLS_WITH_CHEMO = {**GO_COLS, "Chemosensory": "#8B0000"}

# Gene counts from per-family IQ-TREE Newick files (authoritative — some seqs pruned vs FASTA)
# H. illucens (Hill) total in trees: 546 (OR=190 not 192, TRP=92 not 111 after deduplication)
FAMILY_SPECIES_COUNTS = {
    "CSP": {"Hill": 10, "Dmel": 4, "Mdom": 5, "GmmC": 5},
    "GR":  {"Hill": 39, "Aaeg": 107, "Dmel": 68},
    "IR":  {"Hill": 101, "Aaeg": 135, "Dmel": 59, "Csty": 21},
    "OBP": {"Hill": 75, "Dmel": 51, "Mdom": 87, "Aaeg": 56},
    "OR":  {"Hill": 190, "Dmel": 62, "Csty": 41, "Mdom": 80, "Aaeg": 117, "Bdor": 1},
    "PPK": {"Hill": 39},
    "TRP": {"Hill": 92, "Aaeg": 10, "Dmel": 16},
}
SPECIES_FULL_NAMES = {
    "Hill": "Hermetia illucens",
    "Dmel": "Drosophila melanogaster",
    "Aaeg": "Aedes aegypti",
    "Mdom": "Musca domestica",
    "Csty": "Calliphora stygia",
    "GmmC": "Glossina morsitans",
    "Bdor": "Bactrocera dorsalis",
}
SPECIES_COLORS = {
    "Hill": "#000000",  # black (extracted from PDF)
    "Dmel": "#7F59A5",  # purple
    "Aaeg": "#EF4538",  # red
    "Mdom": "#F6891F",  # orange
    "Csty": "#94C57E",  # light green
    "GmmC": "#2A86C7",  # blue
    "Bdor": "#BE922D",  # gold
}

# Chemosensory pie colors (identical across t1 and c1)
CHEMO_PIE_COLORS = {
    "Appendage-specific": "#A60000",
    "Appendage-biased": "#1C1B8D",
    "Expressed": "#555555",
    "Not expressed": "#D6D6D6",
}

# Class columns (identical across files)
CLASS_COLS = ["Antenna_Class", "Palp_Class", "Tarsi_Class"]

# Tissue prefixes in normalized counts (identical across t1 and c1)
TISSUE_PREFIXES = {
    "Antenna": ["Ant_"],
    "Maxillary palp": ["P_", "Palp_"],
    "Tarsi": ["Leg_", "Tar_"],
}

# Figure 2C partitions variance three ways - appendage, sex and mating status -
# plus the unexplained remainder. The app used to read the 2-way model from
# Additional file 32, which merged sex and mating into one slice.
VP_PIE_COLORS = {
    "Appendage identity": "#07045F",
    "Sex": "#00C9FB",
    "Mating status": "#7B1FA2",
    "Shared / confounded": "#9E9E9E",
    "Unexplained": "#D6D6D6",
    # retained so the Python fallback keeps its colours
    "Appendage (unique effect)": "#07045F",
    "Reproductive state (unique effect)": "#00C9FB",
    "Unexplained variance (other factors)": "#D6D6D6",
}


# =============================================================================
# === APP-SPECIFIC CONSTANTS ===
# =============================================================================

# --- T1 constants ---
T1_DEFAULT_PATHS = {
    "BASE_DIR": BASE_DIR,
    "FILE_ANT_P": "condition_vs_Ant_vs_P_name.csv",
    "FILE_ANT_LEG": "condition_vs_Ant_vs_Leg_name.csv",
    "FILE_LEG_P": "condition_vs_Leg_vs_P_name.csv",
    "FILE_GO": "BSF_all-rna_GO_ID_annotated.csv",
    "FILE_NORM": "BSF_normalized_counts_nameX.csv",
    "CHEMO_DIR": BASE_DIR,
    "FILE_LENGTHS": "",
    "FILE_FLY_BASE": "fly_base.png",
    "FILE_FLY_ANT": "fly_antenna.png",
    "FILE_FLY_PALP": "fly_palp.png",
    "FILE_FLY_TARSI": "fly_tarsi.png",
}

TISSUE_LABELS = ["Antenna", "Maxillary palp", "Tarsi"]

GO_FULL = {
    "MF": "Molecular function",
    "CC": "Cellular component",
    "BP": "Biological process",
    "Unknown": "Unknown",
}

CHEMO_FAM_ORDER = ["Or", "Gr", "Ir", "Obp", "Csp", "Ppk", "Trp"]

STATE_COLORS = {
    "Vm": "#ADD8E6",
    "VF": "#FFC0CB",
    "MF": "#800080",
}

DEFAULT_EXPR_THR = 10.0

# --- C1 constants ---
C1_DEFAULT_PATHS = {
    "BASE_DIR": BASE_DIR,
    "FILE_ANT_P": "condition_vs_Ant_vs_P_name.csv",
    "FILE_ANT_LEG": "condition_vs_Ant_vs_Leg_name.csv",
    "FILE_LEG_P": "condition_vs_Leg_vs_P_name.csv",
    "FILE_GO": "BSF_all-rna_GO_ID_annotated.csv",
    "FILE_NORM": "BSF_normalized_counts_nameX.csv",
    "CHEMO_DIR": BASE_DIR,
}

PAL_GO = {**GO_COLS, "Labeled": "#A60000", "Not Significant": "#D9D9D9"}

DOMAIN_FULL = {
    "MF": "Molecular function",
    "CC": "Cellular component",
    "BP": "Biological process",
    "Unknown": "Unknown",
}

CHEMO_TAGS = ("OR", "IR", "GR", "OBP", "CSP", "PPK", "ORCO")

C1_TISSUE_LABELS = {
    "Antenna": "Antenna",
    "Palp": "Maxillary palp",
    "Tarsi": "Tarsi",
}

# --- S1 constants ---
S1_DEFAULT_PATHS = {
    "BASE_DIR": BASE_DIR,
    "DE_FILES": {
        "Antenna - MF vs VF": "results_ant_MF_vs_VF.csv",
        "Antenna - VF vs Vm": "results_ant_VF_vs_Vm.csv",
        "Antenna - MF vs Vm": "results_ant_MF_vs_Vm.csv",
        "Palp - MF vs VF":    "results_palp_MF_vs_VF.csv",
        "Palp - VF vs Vm":    "results_palp_VF_vs_Vm.csv",
        "Palp - MF vs Vm":    "results_palp_MF_vs_Vm.csv",
        "Tarsi - MF vs VF":   "results_leg_MF_vs_VF.csv",
        "Tarsi - VF vs Vm":   "results_leg_VF_vs_Vm.csv",
        "Tarsi - MF vs Vm":   "results_leg_MF_vs_Vm.csv",
    },
    "FILE_GO":   "BSF_all-rna_GO_ID_annotated.csv",
    "FILE_NORM": "BSF_normalized_counts_nameX.csv",
}

TISSUE_DISPLAY = {
    "Antenna": "Antenna",
    "Palp": "Maxillary palp",
    "Tarsi": "Tarsi",
}

COND_FULL = {
    "MF": "Mated female",
    "VF": "Virgin female",
    "Vm": "Virgin male",
    "VM": "Virgin male",
}

COL_MF = dict(light="#C77CFF", dark="#7B1FA2")
COL_VF = dict(light="#FF66B3", dark="#99004D")
COL_Vm = dict(light="#99B3FF", dark="#003399")

S1_PAL_GO = {**GO_COLS, "Not Significant": "#D9D9D9"}

VENN_PADJ_THR = 0.001
VENN_LFC_THR = 1.0

# --- H1 constants ---
H1_NORM_COUNTS_FILE = pathlib.Path(BASE_DIR) / "BSF_normalized_counts_nameX.csv"

H1_THRESHOLD = 10.0


# =============================================================================
# === SHARED HELPER FUNCTIONS ===
# =============================================================================

def make_joinkey(df: pd.DataFrame) -> pd.Series:
    """Create a JoinKey from Gene or Name column."""
    if "Gene" in df.columns:
        k = df["Gene"].astype(str)
    elif "Name" in df.columns:
        k = df["Name"].astype(str)
    else:
        raise ValueError("Neither 'Gene' nor 'Name' present in table.")
    return k.str.strip()


def clean_go_domain(x):
    x0 = (
        x.astype(str)
        .str.strip()
        .str.lower()
        .str.replace("[^a-z]", "", regex=True)
    )
    mapd = {
        "mf": "MF",
        "molecularfunction": "MF",
        "cc": "CC",
        "cellularcomponent": "CC",
        "bp": "BP",
        "biologicalprocess": "BP",
    }
    return x0.map(mapd).fillna("Unknown")


@st.cache_data(show_spinner=False)
def load_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


@st.cache_data(show_spinner=False)
def load_go_map(path: str) -> pd.DataFrame:
    go_raw = pd.read_csv(path)
    nm = [c for c in go_raw.columns]
    first_col = nm[0]

    def pick(colnames, opts):
        low = {c.lower(): c for c in colnames}
        for k in opts:
            if k in low:
                return low[k]
        return None

    go_name_col = pick(nm, [
        "go_name", "go term name", "go_term_name", "go term", "go_description",
        "go_desc", "go label", "goname", "go"
    ])
    go_domain_col = pick(nm, [
        "go_domain", "domain", "aspect", "goaspect", "go_domain_name", "go_aspect"
    ])

    out = pd.DataFrame({
        "Gene": go_raw[first_col].astype(str).str.strip(),
        "GO_Name": go_raw[go_name_col].astype(str).str.strip() if go_name_col else "Unknown",
        "GO_Domain": go_raw[go_domain_col] if go_domain_col else "Unknown",
    })
    out["GO_Domain"] = clean_go_domain(out["GO_Domain"])
    out.loc[out["GO_Name"].isna() | (out["GO_Name"] == ""), "GO_Name"] = "Unknown"
    out = out.drop_duplicates(subset=["Gene"], keep="first").reset_index(drop=True)
    return out


def normalize_class_column(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.strip()
    out = s.copy()
    mask_spec = s.str.contains("specific", case=False, na=False)
    out[mask_spec] = "Appendage-specific"
    mask_bias = s.str.contains("biased", case=False, na=False)
    out[mask_bias] = "Appendage-biased"
    mask_not = s.str.contains("not expressed", case=False, na=False)
    out[mask_not] = "Not expressed"
    mask_expr = s.str.match(r"^expressed$", case=False, na=False)
    out[mask_expr] = "Expressed"
    return out


def deduplicate_family_for_pies(
    fam_df: pd.DataFrame,
    gene_col: str = "Gene_ID",
    cluster_col: str = "Cluster_Rep",
) -> pd.DataFrame:
    df = fam_df.copy()
    if df.empty:
        df["Weight"] = pd.Series(dtype="int64")
        return df
    if gene_col not in df.columns:
        df["Weight"] = 1
        return df
    gene_id_clean = df[gene_col].astype(str).str.strip()
    if cluster_col in df.columns:
        rep = df[cluster_col].astype(str).str.strip()
        bad = rep.isna() | (rep == "") | (rep.str.upper() == "NA")
        rep_key = rep.where(~bad, gene_id_clean)
    else:
        rep_key = gene_id_clean
    df["_rep_key"] = rep_key
    sort_cols = ["_rep_key"]
    if "Is_Duplicate" in df.columns:
        sort_cols.append("Is_Duplicate")
    df_sorted = df.sort_values(sort_cols, ascending=True)
    rep_df = (
        df_sorted
        .groupby("_rep_key", as_index=False)
        .head(1)
        .copy()
    )
    rep_df["Weight"] = 1
    return rep_df


def hex_to_rgba(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# =============================================================================
# === APP_T1 FUNCTIONS ===
# =============================================================================

@st.cache_data(show_spinner=False)
def t1_load_norm_counts(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["JoinKey"] = make_joinkey(df)
    return df


def t1_compute_tissue_means(matrix_df: pd.DataFrame) -> pd.DataFrame:
    df = matrix_df.copy()
    if "JoinKey" not in df.columns:
        df["JoinKey"] = make_joinkey(df)
    out = df[["JoinKey"]].copy()
    for tissue, prefixes in TISSUE_PREFIXES.items():
        cols = [c for c in df.columns if any(c.startswith(p) for p in prefixes)]
        if cols:
            out[tissue + "_mean"] = df[cols].apply(
                pd.to_numeric, errors="coerce"
            ).mean(axis=1)
    out = out.drop_duplicates(subset=["JoinKey"])
    return out


@st.cache_data(show_spinner=False)
def t1_load_gene_lengths(path: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    df = pd.read_csv(path)
    nm = list(df.columns)

    def pick(colnames, opts):
        low = {c.lower(): c for c in colnames}
        for k in opts:
            if k in low:
                return low[k]
        return None

    gene_col = pick(nm, ["gene", "name", "gene_id", "id"])
    len_col = pick(nm, ["length_bp", "length", "len", "gene_length_bp", "gene_length"])
    if gene_col is None or len_col is None:
        return pd.DataFrame()
    out = df[[gene_col, len_col]].copy()
    out.columns = ["Gene", "Length_bp"]
    out["Gene"] = out["Gene"].astype(str).str.strip()
    out = out.dropna(subset=["Length_bp"])
    out = out[out["Length_bp"] > 0]
    return out


def t1_derive_gene_lengths_from_counts(norm_counts: pd.DataFrame) -> pd.DataFrame:
    df = norm_counts.copy()
    if "Length" not in df.columns:
        return pd.DataFrame()
    if "JoinKey" not in df.columns:
        df["JoinKey"] = make_joinkey(df)
    out = df[["JoinKey", "Length"]].copy()
    out.columns = ["Gene", "Length_bp"]
    out["Gene"] = out["Gene"].astype(str).str.strip()
    out["Length_bp"] = pd.to_numeric(out["Length_bp"], errors="coerce")
    out = out.dropna(subset=["Length_bp"])
    out = out[out["Length_bp"] > 0]
    out = out.drop_duplicates(subset=["Gene"])
    return out


@st.cache_data(show_spinner=False)
def t1_compute_fpkm_matrix(norm_counts: pd.DataFrame,
                            gene_lengths: pd.DataFrame) -> pd.DataFrame:
    if gene_lengths.empty:
        return pd.DataFrame()
    df = norm_counts.copy()
    df["JoinKey"] = make_joinkey(df)
    gl = gene_lengths.copy()
    gl["Gene"] = gl["Gene"].astype(str).str.strip()
    merged = df.merge(gl, left_on="JoinKey", right_on="Gene", how="inner")
    if merged.empty:
        return pd.DataFrame()
    length_kb = merged["Length_bp"].astype(float) / 1000.0
    sample_cols = [
        c for c in merged.columns
        if c not in ["Gene", "Name", "JoinKey", "Length_bp", "Length"]
           and not c.startswith("Unnamed")
    ]
    fpkm = merged[["JoinKey"]].copy()
    for col in sample_cols:
        counts = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
        lib_size = counts.sum()
        if lib_size <= 0:
            fpkm[col] = 0.0
            continue
        fpkm[col] = (counts / length_kb) / (lib_size / 1e6)
    return fpkm


@st.cache_data(show_spinner=False)
def t1_load_chemo_tables(chemo_dir: str):
    family_tables = {}
    for fam in CHEMO_FAM_ORDER:
        csv_path = os.path.join(chemo_dir, f"chemo_{fam}.csv")
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        family_tables[fam] = df
    return family_tables


@st.cache_data(show_spinner=False)
def t1_load_fly_images(defaults):
    base_dir = defaults["BASE_DIR"]

    def load_one(filename):
        if not filename:
            return None
        path = os.path.join(base_dir, filename)
        if not os.path.exists(path):
            return None
        img = Image.open(path).convert("RGBA")
        return img

    img_base = load_one(defaults.get("FILE_FLY_BASE", ""))
    img_ant = load_one(defaults.get("FILE_FLY_ANT", ""))
    img_palp = load_one(defaults.get("FILE_FLY_PALP", ""))
    img_tarsi = load_one(defaults.get("FILE_FLY_TARSI", ""))
    return {
        "base": img_base,
        "Antenna": img_ant,
        "Maxillary palp": img_palp,
        "Tarsi": img_tarsi,
    }


def t1_make_tissue_fly(tissue: str, fly_images: dict):
    base = fly_images.get("base")
    overlay = fly_images.get(tissue)
    if base is None and overlay is None:
        return None
    if base is None:
        return overlay
    if overlay is None:
        return base
    canvas = base.copy()
    canvas.alpha_composite(overlay)
    return canvas


@st.cache_data(show_spinner=False)
def t1_load_all_data(defaults):
    base_dir = defaults["BASE_DIR"]
    go_path = os.path.join(base_dir, defaults["FILE_GO"])
    norm_path = os.path.join(base_dir, defaults["FILE_NORM"])
    chemo_dir = defaults.get("CHEMO_DIR", base_dir)
    len_path = os.path.join(base_dir, defaults.get("FILE_LENGTHS", "")) if defaults.get("FILE_LENGTHS") else ""

    go_map = load_go_map(go_path)
    norm_counts = t1_load_norm_counts(norm_path)
    chemo_tables = t1_load_chemo_tables(chemo_dir)

    gene_lengths_file = t1_load_gene_lengths(len_path) if len_path else pd.DataFrame()
    if gene_lengths_file.empty:
        gene_lengths = t1_derive_gene_lengths_from_counts(norm_counts)
    else:
        gene_lengths = gene_lengths_file

    ant_p_path = os.path.join(base_dir, defaults["FILE_ANT_P"])
    ant_leg_path = os.path.join(base_dir, defaults["FILE_ANT_LEG"])
    leg_p_path = os.path.join(base_dir, defaults["FILE_LEG_P"])

    AntP = load_csv(ant_p_path) if os.path.exists(ant_p_path) else pd.DataFrame()
    AntLeg = load_csv(ant_leg_path) if os.path.exists(ant_leg_path) else pd.DataFrame()
    LegP = load_csv(leg_p_path) if os.path.exists(leg_p_path) else pd.DataFrame()

    return {
        "GO_MAP": go_map,
        "NORM_COUNTS": norm_counts,
        "GENE_LENGTHS": gene_lengths,
        "CHEMO_TABLES": chemo_tables,
        "DE": {"AntP": AntP, "AntLeg": AntLeg, "LegP": LegP},
    }


def t1_summarise_expression_by_tissue(norm_means, go_map, expr_thr):
    summaries = {}
    for tissue in TISSUE_LABELS:
        col = tissue + "_mean"
        if col not in norm_means.columns:
            continue
        df_t = norm_means[norm_means[col] >= expr_thr].copy()
        n_total = df_t.shape[0]
        merged = df_t.merge(go_map, left_on="JoinKey", right_on="Gene", how="left")
        merged["GO_Domain"] = clean_go_domain(merged["GO_Domain"])
        counts = (
            merged["GO_Domain"]
            .value_counts()
            .reindex(["MF", "CC", "BP", "Unknown"], fill_value=0)
        )
        summaries[tissue] = {
            "total_expressed": int(n_total),
            "go_counts": counts.to_dict(),
        }
    return summaries


def t1_summarise_chemo_expression_by_threshold(chemo_tables, expr_table, expr_thr):
    summary = {"Antenna": {}, "Maxillary palp": {}, "Tarsi": {}}
    if expr_table is None or expr_table.empty:
        return summary
    if "JoinKey" not in expr_table.columns:
        return summary
    expr_table = expr_table.copy()
    expr_table["JoinKey_clean"] = expr_table["JoinKey"].astype(str).str.strip().str.upper()
    if "Name" in expr_table.columns:
        expr_table["Name_clean"] = expr_table["Name"].astype(str).str.strip().str.upper()
    expressed_ids = {}
    for tissue in TISSUE_LABELS:
        col = tissue + "_mean"
        if col not in expr_table.columns:
            continue
        mask = expr_table[col] >= expr_thr
        keys = set(expr_table.loc[mask, "JoinKey_clean"])
        if "Name_clean" in expr_table.columns:
            keys |= set(expr_table.loc[mask, "Name_clean"])
        expressed_ids[tissue] = keys
    for fam in CHEMO_FAM_ORDER:
        fam_df = chemo_tables.get(fam, pd.DataFrame()).copy()
        if fam_df.empty or "Gene_ID" not in fam_df.columns:
            continue
        fam_rep = deduplicate_family_for_pies(fam_df)
        fam_rep["Gene_ID_clean"] = fam_rep["Gene_ID"].astype(str).str.strip().str.upper()
        fam_counts = {}
        for tissue in TISSUE_LABELS:
            ids_t = expressed_ids.get(tissue, set())
            if not ids_t:
                fam_counts[tissue] = 0
                continue
            mask = fam_rep["Gene_ID_clean"].isin(ids_t)
            n_expr = int(fam_rep.loc[mask, "Weight"].sum())
            fam_counts[tissue] = n_expr
        for tissue in TISSUE_LABELS:
            if tissue not in summary:
                summary[tissue] = {}
            summary[tissue][fam.upper()] = fam_counts[tissue]
    return summary


def t1_build_all_summaries(go_map, norm_means, expr_table, norm_counts, chemo_tables, expr_thr):
    expr_summary = t1_summarise_expression_by_tissue(norm_means, go_map, expr_thr)
    chemo_summary = t1_summarise_chemo_expression_by_threshold(chemo_tables, expr_table, expr_thr)
    return expr_summary, chemo_summary


def t1_get_go_gene_table(go_map, norm_means, tissue, expr_thr):
    col = tissue + "_mean"
    if col not in norm_means.columns:
        return pd.DataFrame()
    df_t = norm_means[norm_means[col] >= expr_thr].copy()
    if df_t.empty:
        return pd.DataFrame()
    merged = df_t.merge(go_map, left_on="JoinKey", right_on="Gene", how="left")
    merged["GO_Domain"] = clean_go_domain(merged["GO_Domain"])
    merged = merged.rename(columns={col: "Expression"})
    keep_cols = ["Gene", "GO_Name", "GO_Domain", "Expression"]
    keep_cols = [c for c in keep_cols if c in merged.columns]
    merged = merged[keep_cols].drop_duplicates()
    return merged


def t1_get_chemo_gene_table(chemo_tables, expr_table, tissue, family_key, expr_thr):
    if expr_table is None or expr_table.empty:
        return pd.DataFrame()
    expr_table = expr_table.copy()
    expr_table["JoinKey_clean"] = expr_table["JoinKey"].astype(str).str.strip().str.upper()
    if "Name" in expr_table.columns:
        expr_table["Name_clean"] = expr_table["Name"].astype(str).str.strip().str.upper()
    col = tissue + "_mean"
    if col not in expr_table.columns:
        return pd.DataFrame()
    mask_expr = expr_table[col] >= expr_thr
    expr_sub = expr_table.loc[mask_expr].copy()
    ids = set(expr_sub["JoinKey_clean"])
    if "Name_clean" in expr_sub.columns:
        ids |= set(expr_sub["Name_clean"])
    fam = family_key.capitalize()
    fam_df = chemo_tables.get(fam, pd.DataFrame()).copy()
    if fam_df.empty or "Gene_ID" not in fam_df.columns:
        return pd.DataFrame()
    fam_rep = deduplicate_family_for_pies(fam_df)
    fam_rep["Gene_ID_clean"] = fam_rep["Gene_ID"].astype(str).str.strip().str.upper()
    mask = fam_rep["Gene_ID_clean"].isin(ids)
    fam_sel = fam_rep.loc[mask].copy()
    if fam_sel.empty:
        return pd.DataFrame()
    drop_cols = {"_rep_key", "Weight", "Gene_ID_clean"}
    cols = [c for c in fam_sel.columns if c not in drop_cols]
    fam_sel = fam_sel[cols].drop_duplicates()
    return fam_sel


def t1_parse_sample_meta(sample_name: str):
    tissue = None
    rest = None
    if sample_name.startswith("Ant_"):
        tissue = "Antenna"
        rest = sample_name[4:]
    elif sample_name.startswith("Leg_"):
        tissue = "Tarsi"
        rest = sample_name[4:]
    elif sample_name.startswith("P_"):
        tissue = "Maxillary palp"
        rest = sample_name[2:]
    elif sample_name.startswith("Palp_"):
        tissue = "Maxillary palp"
        rest = sample_name[5:]
    else:
        return None, None, None, None
    if rest.startswith("VF"):
        state = "Virgin female"
        short = "VF"
    elif rest.startswith("Vm"):
        state = "Virgin male"
        short = "Vm"
    elif rest.startswith("MF"):
        state = "Mated female"
        short = "MF"
    else:
        return None, None, None, None
    group = f"{tissue} - {state}"
    return tissue, state, short, group


@st.cache_data(show_spinner=False)
def t1_compute_pca_by_group(norm_counts: pd.DataFrame, n_components: int = 2):
    df = norm_counts.copy()
    candidate_cols = [
        c for c in df.columns
        if c not in ["Gene", "Name", "JoinKey", "Length"]
           and not c.startswith("Unnamed")
    ]
    sample_cols = []
    tissues = []
    states = []
    shorts = []
    groups = []
    for c in candidate_cols:
        t, s, short_s, g = t1_parse_sample_meta(c)
        if t is None:
            continue
        sample_cols.append(c)
        tissues.append(t)
        states.append(s)
        shorts.append(short_s)
        groups.append(g)
    if len(sample_cols) < 2:
        return pd.DataFrame()
    expr = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X = np.log2(expr.values + 1.0)
    X = X.T  # samples × genes
    # Centre then scale each gene (matches R prcomp(scale.=TRUE))
    X = X - X.mean(axis=0, keepdims=True)
    gene_std = X.std(axis=0, ddof=1)
    gene_std[gene_std == 0] = 1.0
    X = X / gene_std
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    PCs = U[:, :n_components] * S[:n_components]
    pc_sample = pd.DataFrame({
        "Sample": sample_cols,
        "PC1": PCs[:, 0],
        "PC2": PCs[:, 1] if n_components > 1 else 0.0,
        "Tissue": tissues,
        "State": states,
        "State_short": shorts,
        "Group": groups,
    })
    pc_group = (
        pc_sample
        .groupby(["Tissue", "State", "State_short", "Group"], as_index=False)[["PC1", "PC2"]]
        .mean()
    )
    return pc_group


@st.cache_data(show_spinner=False)
def t1_compute_pca_loadings(norm_counts: pd.DataFrame, n_components: int = 5) -> pd.DataFrame:
    df = norm_counts.copy()
    if "JoinKey" not in df.columns:
        df["JoinKey"] = make_joinkey(df)
    candidate_cols = [
        c for c in df.columns
        if c not in ["Gene", "Name", "JoinKey", "Length"]
           and not c.startswith("Unnamed")
    ]
    sample_cols = []
    for c in candidate_cols:
        t, s, short_s, g = t1_parse_sample_meta(c)
        if t is None:
            continue
        sample_cols.append(c)
    if len(sample_cols) < 2:
        return pd.DataFrame()
    expr = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X = np.log2(expr.values + 1.0)
    X = X.T  # samples × genes
    X = X - X.mean(axis=0, keepdims=True)
    gene_std = X.std(axis=0, ddof=1)
    gene_std[gene_std == 0] = 1.0
    X = X / gene_std  # scale to match R prcomp(scale.=TRUE)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    k = min(n_components, Vt.shape[0])
    Vt_k = Vt[:k, :]
    load_df = pd.DataFrame({"JoinKey": df["JoinKey"].values})
    if "Name" in df.columns:
        load_df["Name"] = df["Name"].astype(str).values
    for i in range(k):
        load_df[f"PC{i+1}"] = Vt_k[i, :]
    return load_df


@st.cache_data(show_spinner=False)
def t1_compute_group_means_for_corr(norm_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Match the R script exactly:
    - restrict to the 27 library columns
    - log2(count + 1) FIRST
    - mean in log space per pattern (9 groups)
    """
    df = norm_counts.copy()

    # Keep JoinKey for downstream tables if you need it
    if "JoinKey" not in df.columns:
        df["JoinKey"] = make_joinkey(df)

    # Exact libraries from the R script
    libs = [
        "Ant_MF1", "Ant_MF2", "Ant_MF3",
        "Ant_VF1", "Ant_VF2", "Ant_VF3",
        "Ant_Vm1", "Ant_Vm2", "Ant_Vm3",
        "Leg_MF1", "Leg_MF2", "Leg_MF3",
        "Leg_VF1", "Leg_VF2", "Leg_VF3",
        "Leg_Vm1", "Leg_Vm2", "Leg_Vm3",
        "P_MF1",   "P_MF2",   "P_MF3",
        "P_VF1",   "P_VF2",   "P_VF3",
        "P_Vm1",   "P_Vm2",   "P_Vm3",
    ]
    libs_present = [c for c in libs if c in df.columns]
    if not libs_present:
        return pd.DataFrame()

    # Numeric, keep NaN as NaN (R uses pairwise.complete.obs later)
    counts_sel = df[libs_present].apply(pd.to_numeric, errors="coerce")

    # IMPORTANT: log transform before averaging (R behavior)
    log_counts = np.log2(counts_sel + 1.0)

    out = df[["JoinKey"]].copy()

    group_patterns = [
        "Ant_Vm", "Ant_VF", "Ant_MF",
        "Leg_Vm", "Leg_VF", "Leg_MF",
        "P_Vm",   "P_VF",   "P_MF",
    ]

    for pattern in group_patterns:
        cols = [c for c in log_counts.columns if c.startswith(pattern)]
        if not cols:
            continue
        out[pattern + "_mean"] = log_counts[cols].mean(axis=1, skipna=True)

    return out


def t1_compute_correlation_matrix(group_means: pd.DataFrame) -> pd.DataFrame:
    """
    Match R:
    cor(group_means_mat, method="pearson", use="pairwise.complete.obs")
    Pandas .corr uses pairwise deletion by default.
    """
    if group_means is None or group_means.empty:
        return pd.DataFrame()

    expr = group_means.drop(columns=["JoinKey"], errors="ignore")
    expr = expr.dropna(axis=1, how="all")

    if expr.shape[1] < 2:
        return pd.DataFrame()

    # Pairwise complete obs is the default behavior for pandas corr.
    return expr.corr(method="pearson", min_periods=1)

def t1_parse_corr_label(label: str):
    parts = label.split("_")
    if len(parts) < 2:
        return label, label, "#000000"
    prefix = parts[0]
    state_code = parts[1]
    if prefix == "Ant":
        tissue = "Antenna"
    elif prefix == "Leg":
        tissue = "Tarsi"
    elif prefix in ("P", "Palp"):
        tissue = "Maxillary palp"
    else:
        tissue = prefix
    if state_code == "Vm":
        state_full = "Virgin male"
    elif state_code == "VF":
        state_full = "Virgin female"
    elif state_code == "MF":
        state_full = "Mated female"
    else:
        state_full = state_code
    color = STATE_COLORS.get(state_code, "#000000")
    return tissue, state_full, color


def t1_build_correlation_figure(corr: pd.DataFrame):
    if corr.empty:
        return None

    def short_tissue_label(tissue: str) -> str:
        if tissue == "Maxillary palp":
            return "M. Palp"
        return tissue

    corr_plot = corr.T.copy()
    labels = list(corr_plot.columns)
    rows = []
    for i, row_lab in enumerate(labels):
        for j, col_lab in enumerate(labels):
            if j < i:
                continue
            r = float(corr_plot.loc[row_lab, col_lab])
            tissue_row, state_full_row, color_row = t1_parse_corr_label(row_lab)
            tissue_col, state_full_col, color_col = t1_parse_corr_label(col_lab)
            parts_row = row_lab.split("_")
            state_code_row = parts_row[1] if len(parts_row) > 1 else ""
            parts_col = col_lab.split("_")
            state_code_col = parts_col[1] if len(parts_col) > 1 else ""
            rows.append({
                "Row": row_lab, "Col": col_lab, "corr": r, "abs_corr": abs(r),
                "Row_tissue_short": short_tissue_label(tissue_row),
                "Row_state_code": state_code_row,
                "Row_state_full": state_full_row,
                "Row_color": color_row,
                "Col_tissue_short": short_tissue_label(tissue_col),
                "Col_state_code": state_code_col,
                "Col_state_full": state_full_col,
                "Col_color": color_col,
            })
    df_long = pd.DataFrame(rows)
    x_order = labels[::-1]
    y_order = labels[::-1]
    fig = px.scatter(
        df_long, x="Col", y="Row", size="abs_corr", color="corr",
        color_continuous_scale="Blues", range_color=(0.0, 1.0), size_max=11,
        custom_data=["Row_tissue_short", "Row_state_full", "Row_color",
                     "Col_tissue_short", "Col_state_full", "Col_color", "corr"],
    )
    fig.update_traces(
        hovertemplate=(
            "<b><span style='color:black'>%{customdata[0]}</span> "
            "<span style='color:%{customdata[2]}'>%{customdata[1]}</span></b>"
            " vs "
            "<b><span style='color:black'>%{customdata[3]}</span> "
            "<span style='color:%{customdata[5]}'>%{customdata[4]}</span></b>"
            "<br>Pearson r = %{customdata[6]:.2f}<extra></extra>"
        ),
        marker=dict(line=dict(width=0)),
    )
    x_ticktext = []
    y_ticktext = []
    for lab in x_order:
        tissue, state_full, color = t1_parse_corr_label(lab)
        tissue_short = short_tissue_label(tissue)
        parts = lab.split("_")
        state_code = parts[1] if len(parts) > 1 else ""
        txt = (f"<span style='color:black'>{tissue_short}</span> "
               f"<span style='color:{color}'>{state_code}</span>")
        x_ticktext.append(txt)
    for lab in y_order:
        tissue, state_full, color = t1_parse_corr_label(lab)
        tissue_short = short_tissue_label(tissue)
        parts = lab.split("_")
        state_code = parts[1] if len(parts) > 1 else ""
        txt = (f"<span style='color:black'>{tissue_short}</span> "
               f"<span style='color:{color}'>{state_code}</span>")
        y_ticktext.append(txt)
    fig.update_xaxes(title="", categoryorder="array", categoryarray=x_order,
                     showgrid=False, ticks="outside", tickmode="array",
                     tickvals=x_order, ticktext=x_ticktext)
    fig.update_yaxes(title="", categoryorder="array", categoryarray=y_order,
                     showgrid=False, ticks="outside", tickmode="array",
                     tickvals=y_order, ticktext=y_ticktext)
    fig.update_layout(
        height=600, width=600, yaxis_scaleanchor="x", yaxis_scaleratio=1,
        coloraxis_colorbar=dict(title="Pearson\ncorrelation"), showlegend=False,
    )
    return fig


def t1_fig_chemo_pie(summary, fam_key: str, tissue_label: str, show_legend: bool):
    fam_dict = summary.get(fam_key, {})
    counts = fam_dict.get(tissue_label, {})
    if not counts:
        return px.scatter()
    df = pd.DataFrame({"Class": list(counts.keys()), "Count": list(counts.values())})
    df = df[df["Count"] > 0]
    if df.empty:
        return px.scatter()
    fig = px.pie(df, names="Class", values="Count", color="Class",
                 color_discrete_map=CHEMO_PIE_COLORS, hole=0)
    fig.update_traces(textinfo="value",
                      hovertemplate="%{label}: %{value} genes (%{percent:.1%})<extra></extra>",
                      showlegend=False)
    base_layout = dict(margin=dict(l=5, r=5, t=10, b=60), height=320)
    if show_legend:
        legend_items = [("Appendage-specific", "Specific"), ("Appendage-biased", "Biased"),
                        ("Expressed", "Expressed"), ("Not expressed", "Not expressed")]
        for key, label in legend_items:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                     marker=dict(symbol="square", size=10,
                                                 color=CHEMO_PIE_COLORS[key]),
                                     name=label, showlegend=True, hoverinfo="skip"))
        base_layout["showlegend"] = True
        base_layout["legend"] = dict(orientation="h", yanchor="top", y=-0.05,
                                     xanchor="center", x=0.5, font=dict(size=10))
    else:
        base_layout["showlegend"] = False
    fig.update_layout(**base_layout)
    return fig


def t1_compute_chemo_pie_summary_from_expr(expr_table, chemo_tables, expr_thr):
    if expr_table is None or expr_table.empty or not chemo_tables:
        return {}
    et = expr_table.copy()
    for tissue in TISSUE_LABELS:
        col = tissue + "_mean"
        if col not in et.columns:
            et[col] = 0.0
        et[col] = pd.to_numeric(et[col], errors="coerce").fillna(0.0)
    et["JoinKey_clean"] = et["JoinKey"].astype(str).str.strip().str.upper()
    if "Name" in et.columns:
        et["Name_clean"] = et["Name"].astype(str).str.strip().str.upper()
    expr_map = {}
    for _, row in et.iterrows():
        vals = {tissue: float(row.get(tissue + "_mean", 0.0)) for tissue in TISSUE_LABELS}
        jk = row["JoinKey_clean"]
        if jk:
            expr_map[jk] = vals
        if "Name_clean" in et.columns:
            nm = row["Name_clean"]
            if isinstance(nm, str) and nm:
                expr_map[nm] = vals

    def classify_gene(vals):
        out = {}
        expr_flags = {tissue: (vals.get(tissue, 0.0) >= expr_thr) for tissue in TISSUE_LABELS}
        means = {t: vals.get(t, 0.0) for t in TISSUE_LABELS}
        max_mean = max(means.values()) if means else 0.0
        for tissue in TISSUE_LABELS:
            v = means.get(tissue, 0.0)
            is_expr = expr_flags[tissue]
            others = [t for t in TISSUE_LABELS if t != tissue]
            others_expr = any(expr_flags[o] for o in others)
            if is_expr and not any(expr_flags[o] for o in others):
                cls = "Appendage-specific"
            elif is_expr and others_expr and (v > max(means[o] for o in others)):
                cls = "Appendage-biased"
            elif is_expr:
                cls = "Expressed"
            else:
                cls = "Not expressed"
            out[tissue] = cls
        return out

    chemo_summary = {}
    for fam in CHEMO_FAM_ORDER:
        fam_df = chemo_tables.get(fam, pd.DataFrame()).copy()
        if fam_df.empty or "Gene_ID" not in fam_df.columns:
            continue
        fam_rep = deduplicate_family_for_pies(fam_df)
        fam_rep["Gene_ID_clean"] = fam_rep["Gene_ID"].astype(str).str.strip().str.upper()
        fam_counts = {
            tissue: {"Appendage-specific": 0, "Appendage-biased": 0, "Expressed": 0, "Not expressed": 0}
            for tissue in TISSUE_LABELS
        }
        for _, row in fam_rep.iterrows():
            gid = row["Gene_ID_clean"]
            vals = expr_map.get(gid, {t: 0.0 for t in TISSUE_LABELS})
            classes = classify_gene(vals)
            for tissue in TISSUE_LABELS:
                cls = classes[tissue]
                fam_counts[tissue][cls] += 1
        chemo_summary[fam] = fam_counts
    return chemo_summary
def t1_get_gene_family_from_counts(norm_counts: pd.DataFrame) -> pd.Series:
    """
    Returns a 'gene_fam' Series aligned to norm_counts rows.
    If norm_counts already has gene_fam, use it.
    Else infer from Name/name column.
    """
    if norm_counts is None or norm_counts.empty:
        return pd.Series(dtype=str)

    if "gene_fam" in norm_counts.columns:
        return norm_counts["gene_fam"].astype(str)

    name_col = None
    for c in ["Name", "name", "Gene", "gene"]:
        if c in norm_counts.columns:
            name_col = c
            break
    if name_col is None:
        return pd.Series(["nan"] * len(norm_counts), index=norm_counts.index, dtype=str)

    s = (norm_counts[name_col].astype("object")
         .where(norm_counts[name_col].notna(), "").astype(str).str.upper())

    def infer(x: str) -> str:
        # conservative prefix inference
        if x.startswith("OR"): return "OR"
        if x.startswith("GR"): return "GR"
        if x.startswith("IR"): return "IR"
        if x.startswith("PPK"): return "PPK"
        if x.startswith("OBP"): return "OBP"
        if x.startswith("CSP"): return "CSP"
        if x.startswith("TRP"): return "TRP"
        return "nan"

    return s.apply(infer)


def t1_filter_chemoreceptors(norm_counts: pd.DataFrame) -> pd.DataFrame:
    """
    Chemoreceptors only (match your R intention): OR, GR, IR, PPK, TRP.
    If you want to include OBP/CSP, add them below.
    """
    if norm_counts is None or norm_counts.empty:
        return pd.DataFrame()

    fam = t1_get_gene_family_from_counts(norm_counts)
    keep_fams = {"OR", "GR", "IR", "PPK", "TRP", "OBP", "CSP"}
    keep = fam.isin(list(keep_fams))
    return norm_counts.loc[keep].copy()


def _t1_extract_library_cols(norm_counts: pd.DataFrame) -> list:
    """
    Keep only the 27 RNA libs using the same parsing logic you already use.
    """
    if norm_counts is None or norm_counts.empty:
        return []

    candidate_cols = [c for c in norm_counts.columns if not c.startswith("Unnamed")]
    lib_cols = []
    for c in candidate_cols:
        t, s, short_s, g = t1_parse_sample_meta(c)
        if t is None:
            continue
        lib_cols.append(c)
    return lib_cols


def t1_variance_partition_pie(norm_counts: pd.DataFrame) -> pd.DataFrame:
    # Exact 27 libraries (same as your R script)
    libs = [
        "Ant_MF1", "Ant_MF2", "Ant_MF3",
        "Ant_VF1", "Ant_VF2", "Ant_VF3",
        "Ant_Vm1", "Ant_Vm2", "Ant_Vm3",
        "Leg_MF1", "Leg_MF2", "Leg_MF3",
        "Leg_VF1", "Leg_VF2", "Leg_VF3",
        "Leg_Vm1", "Leg_Vm2", "Leg_Vm3",
        "P_MF1",   "P_MF2",   "P_MF3",
        "P_VF1",   "P_VF2",   "P_VF3",
        "P_Vm1",   "P_Vm2",   "P_Vm3",
    ]
    libs = [c for c in libs if c in norm_counts.columns]
    if len(libs) < 2:
        return pd.DataFrame()

    # counts -> numeric, keep NaN (do not fill with 0)
    counts = norm_counts[libs].apply(pd.to_numeric, errors="coerce")

    # log2 transform first (same as R)
    X = np.log2(counts + 1.0).to_numpy(dtype=float).T  # samples x genes

    # remove zero-variance genes (same as R log_counts_pca)
    gene_sd = np.nanstd(X, axis=0)
    keep = np.isfinite(gene_sd) & (gene_sd > 0)
    X = X[:, keep]
    if X.shape[1] < 2:
        return pd.DataFrame()

    # build metadata in the same way as your app parser
    meta_rows = []
    for c in libs:
        tissue, state, state_short, group = t1_parse_sample_meta(c)
        meta_rows.append({"Library": c, "Tissue": tissue, "State_short": state_short})
    meta = pd.DataFrame(meta_rows)

    # one-hot encodings (drop_first like standard regression)
    Z_tissue = pd.get_dummies(meta["Tissue"], drop_first=True).to_numpy(dtype=float)
    Z_state = pd.get_dummies(meta["State_short"], drop_first=True).to_numpy(dtype=float)

    def adj_r2_multivariate(Y: np.ndarray, Z: np.ndarray) -> float:
        n = Y.shape[0]
        if Z.size == 0:
            return 0.0

        # intercept + predictors
        Z1 = np.column_stack([np.ones((n, 1)), Z])
        # use rank, not just number of columns (closer to vegan behavior)
        rank = np.linalg.matrix_rank(Z1) - 1  # exclude intercept
        if rank < 1:
            return 0.0

        B, _, _, _ = np.linalg.lstsq(Z1, Y, rcond=None)
        Yhat = Z1 @ B

        ss_tot = np.nansum((Y - np.nanmean(Y, axis=0, keepdims=True)) ** 2)
        ss_res = np.nansum((Y - Yhat) ** 2)
        if ss_tot <= 0:
            return 0.0

        r2 = 1.0 - (ss_res / ss_tot)

        denom = max(n - rank - 1, 1)
        adj = 1.0 - (1.0 - r2) * (n - 1) / denom
        return float(max(0.0, min(1.0, adj)))

    # R logic:
    # unique_tissue = adjR2(full) - adjR2(state_only)
    # unique_state  = adjR2(full) - adjR2(tissue_only)
    Z_full = np.column_stack([Z_tissue, Z_state])
    adj_full = adj_r2_multivariate(X, Z_full)
    adj_state_only = adj_r2_multivariate(X, Z_state)
    adj_tissue_only = adj_r2_multivariate(X, Z_tissue)

    tissue_unique = max(0.0, adj_full - adj_state_only)
    state_unique = max(0.0, adj_full - adj_tissue_only)
    unexplained = max(0.0, 1.0 - (tissue_unique + state_unique))

    return pd.DataFrame({
        "Component": [
            "Appendage (unique effect)",
            "Reproductive state (unique effect)",
            "Unexplained variance (other factors)",
        ],
        "Fraction": [tissue_unique, state_unique, unexplained],
    })


def t1_style_pca_fig(fig: go.Figure, height: int = 600, width: int = 600) -> go.Figure:
    fig.update_layout(
        height=height,
        width=width,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=40, r=20, t=40, b=40),
        showlegend=False,
    )
    fig.update_xaxes(showgrid=False, zeroline=False, showline=True, linecolor="black")
    fig.update_yaxes(showgrid=False, zeroline=False, showline=True, linecolor="black")
    return fig
def render_t1_tab():
    # ---------- Load data ----------
    data = t1_load_all_data(T1_DEFAULT_PATHS)
    go_map = data["GO_MAP"]
    norm_counts = data["NORM_COUNTS"]
    gene_lengths = data["GENE_LENGTHS"]
    chemo_tables = data["CHEMO_TABLES"]
    fly_images = t1_load_fly_images(T1_DEFAULT_PATHS)

    # ---------- Sidebar options ----------
    with st.sidebar:
        st.header("Transcriptome Overview")
        unit_choice = st.radio(
            "Expression unit for threshold",
            ["Normalized counts", "FPKM"],
            index=0,
            key="t1_unit_choice",
        )
        expr_thr = st.number_input(
            "Expression threshold",
            min_value=0.0, max_value=1e6, value=DEFAULT_EXPR_THR, step=1.0,
            help="Genes with mean expression above this are counted as expressed.",
            key="t1_expr_thr",
        )
        st.markdown("---")
        st.subheader("FPKM conversion")
        st.markdown(
            "FPKM is computed for each sample as:\n\n"
            "`FPKM = (counts / length_kb) / (library_size / 1e6)`\n\n"
            "where `length_kb` is gene length in kilobases and `library_size` "
            "is the sum of counts for that sample."
        )

    use_fpkm = (unit_choice == "FPKM")
    if use_fpkm:
        fpkm_matrix = t1_compute_fpkm_matrix(norm_counts, gene_lengths)
        if fpkm_matrix.empty:
            st.warning(
                "FPKM selected but could not be computed (no valid gene lengths). "
                "Using normalized counts instead."
            )
            matrix_for_means = norm_counts
            unit_label = "normalized counts"
        else:
            matrix_for_means = fpkm_matrix
            unit_label = "FPKM"
    else:
        matrix_for_means = norm_counts
        unit_label = "normalized counts"

    norm_means = t1_compute_tissue_means(matrix_for_means)

    if "Name" in norm_counts.columns:
        expr_table = norm_means.merge(
            norm_counts[["JoinKey", "Name"]].drop_duplicates(),
            on="JoinKey", how="left",
        )
    else:
        expr_table = norm_means.copy()

    expr_summary, chemo_summary = t1_build_all_summaries(
        go_map, norm_means, expr_table, norm_counts, chemo_tables, expr_thr
    )

    tissue_sel = st.radio("Choose tissue", TISSUE_LABELS, horizontal=True, key="t1_tissue_sel")

    expr = expr_summary.get(tissue_sel, {})
    chemo = chemo_summary.get(tissue_sel, {})

    total = expr.get("total_expressed", 0)
    go_counts = expr.get("go_counts", {})
    mf = go_counts.get("MF", 0)
    cc = go_counts.get("CC", 0)
    bp = go_counts.get("BP", 0)
    unk = go_counts.get("Unknown", 0)

    st.markdown(
        f"**{tissue_sel}:** {total} genes expressed "
        f"(mean {unit_label} >= {expr_thr:.2f})."
    )

    col_fly, col_right = st.columns([1, 2])

    with col_fly:
        img_fly = t1_make_tissue_fly(tissue_sel, fly_images)
        if img_fly is not None:
            st.markdown("**Location on fly**")
            st.image(img_fly)
        else:
            st.info("Fly silhouette images not found or not configured.")

    with col_right:
        st.markdown("**GO domain composition (expressed genes)**")
        go_df = pd.DataFrame({
            "GO_Code": ["MF", "CC", "BP", "Unknown"],
            "GO_Domain": [GO_FULL[c] for c in ["MF", "CC", "BP", "Unknown"]],
            "Count": [mf, cc, bp, unk],
        })
        max_go = max(go_df["Count"].max(), 1)
        fig_go = px.bar(go_df, x="GO_Domain", y="Count", color="GO_Code",
                        color_discrete_map=GO_COLS)
        fig_go.update_layout(showlegend=False, xaxis_title="GO domain",
                             yaxis_title="Number of genes", height=400)
        fig_go.update_yaxes(range=[0, max_go * 1.25])
        st.plotly_chart(fig_go, use_container_width=True)

        with st.expander("Show GO gene table and download CSV"):
            go_table = t1_get_go_gene_table(go_map, norm_means, tissue_sel, expr_thr)
            if go_table.empty:
                st.info("No expressed genes for this tissue at current threshold.")
            else:
                domain_choice_full = st.selectbox(
                    "GO domain",
                    [GO_FULL[c] for c in ["MF", "CC", "BP", "Unknown"]],
                    index=0, key="t1_go_domain_sel",
                )
                rev_map = {v: k for k, v in GO_FULL.items()}
                domain_code = rev_map[domain_choice_full]
                sub = go_table[go_table["GO_Domain"] == domain_code].copy()
                st.write(f"{sub.shape[0]} genes in {domain_choice_full}")
                st.dataframe(sub, use_container_width=True)
                csv = sub.to_csv(index=False).encode("utf-8")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"{tissue_sel}_{domain_code}_GO_genes.csv",
                                   mime="text/csv", key="t1_go_dl")

        st.markdown("---")
        st.markdown("**Chemosensory families (deduplicated, expressed)**")
        chemo_counts = []
        for fam in CHEMO_FAM_ORDER:
            fam_key = fam.upper()
            chemo_counts.append({"Family": fam, "Count": int(chemo.get(fam_key, 0))})
        chemo_df = pd.DataFrame(chemo_counts)
        max_chemo = max(chemo_df["Count"].max(), 1)
        fig_chemo_bar = px.bar(chemo_df, x="Family", y="Count")
        fig_chemo_bar.update_traces(marker_color="black")
        fig_chemo_bar.update_layout(showlegend=False, xaxis_title="Family",
                                    yaxis_title="Number of genes", height=400)
        fig_chemo_bar.update_yaxes(range=[0, max_chemo * 1.25])
        fig_chemo_bar.update_xaxes(tickmode="array", tickvals=CHEMO_FAM_ORDER,
                                   ticktext=[f"<i>{f}</i>" for f in CHEMO_FAM_ORDER])
        st.plotly_chart(fig_chemo_bar, use_container_width=True)

        with st.expander("Show chemosensory genes and download CSV"):
            fam_choice_label = st.selectbox("Chemosensory family", CHEMO_FAM_ORDER,
                                            index=0, key="t1_chemo_fam_sel")
            fam_choice_key = fam_choice_label.upper()
            chemo_table = t1_get_chemo_gene_table(
                chemo_tables, expr_table, tissue_sel, fam_choice_key, expr_thr
            )
            if chemo_table.empty:
                st.info("No chemosensory genes above threshold for this tissue and family.")
            else:
                st.write(f"{chemo_table.shape[0]} genes in {fam_choice_label} for {tissue_sel}")
                st.dataframe(chemo_table, use_container_width=True)
                csv = chemo_table.to_csv(index=False).encode("utf-8")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"{tissue_sel}_{fam_choice_label}_chemo_genes.csv",
                                   mime="text/csv", key="t1_chemo_dl")

    st.markdown("---")
    st.subheader("Chemosensory pies")
    st.caption("Manuscript: **Figure 3D**")

    chemo_pie_summary = t1_compute_chemo_pie_summary_from_expr(expr_table, chemo_tables, expr_thr)

    if not chemo_pie_summary:
        st.info("No chemosensory tables available to build pies.")
    else:
        st.markdown(
            f"<h3 style='text-align:left;'>Tissue: {'M. palp' if tissue_sel == 'Maxillary palp' else tissue_sel}</h3>",
            unsafe_allow_html=True,
        )
        col_p0, col_p1, col_p2 = st.columns(3)
        col_map = {
            "Or": col_p0, "Gr": col_p0,
            "Ir": col_p1, "Obp": col_p1,
            "Csp": col_p2, "Ppk": col_p2, "Trp": col_p2,
        }
        for fam in CHEMO_FAM_ORDER:
            col_here = col_map[fam]
            with col_here:
                fig_p = t1_fig_chemo_pie(chemo_pie_summary, fam_key=fam,
                                         tissue_label=tissue_sel, show_legend=(fam == "Or"))
                if fig_p.data:
                    st.plotly_chart(fig_p, use_container_width=True,
                                    key=f"t1_chemo_pie_{tissue_sel}_{fam}")
                    st.markdown(f"<p style='text-align:center; font-size:0.95rem;'><i>{fam}</i></p>",
                                unsafe_allow_html=True)
                else:
                    st.caption(f"No data for {fam} in {tissue_sel}")

    st.markdown("---")
    st.subheader("PCA")
    st.caption("Manuscript: **Figure 2A**")

    pc_group = t1_compute_pca_by_group(norm_counts)
    group_means = t1_compute_group_means_for_corr(norm_counts)
    pca_loadings = t1_compute_pca_loadings(norm_counts, n_components=5)

    # Load R pre-computed Pearson correlation (matches manuscript Fig 2B exactly)
    _corr_csv = str(pathlib.Path(__file__).resolve().parent / "data" / "stats_reports" / "05_pearson_correlation_long.csv")
    corr_matrix = t1_compute_correlation_matrix(group_means)  # Python fallback
    if os.path.exists(_corr_csv):
        try:
            _corr_raw = pd.read_csv(_corr_csv)
            _needed = [c for c in ["Row", "Col", "corr", "dataset"] if c in _corr_raw.columns]
            if "dataset" in _corr_raw.columns:
                _corr_raw = _corr_raw[_corr_raw["dataset"] == "ALL GENES"]
            if {"Row", "Col", "corr"}.issubset(_corr_raw.columns):
                _rows = _corr_raw["Row"].str.replace("_mean", "", regex=False)
                _cols = _corr_raw["Col"].str.replace("_mean", "", regex=False)
                _labels = sorted(set(_rows) | set(_cols))
                _mat = pd.DataFrame(index=_labels, columns=_labels, dtype=float)
                for _, row in _corr_raw.iterrows():
                    r_lbl = str(row["Row"]).replace("_mean", "")
                    c_lbl = str(row["Col"]).replace("_mean", "")
                    _mat.loc[r_lbl, c_lbl] = float(row["corr"])
                    _mat.loc[c_lbl, r_lbl] = float(row["corr"])
                np.fill_diagonal(_mat.values, 1.0)
                if not _mat.isna().all().all():
                    corr_matrix = _mat.dropna(how="all").dropna(axis=1, how="all")
        except Exception:
            pass  # keep Python-computed fallback

    # Load precomputed % variance from R (matches manuscript Fig 2A)
    _pca_ov_csv = str(pathlib.Path(__file__).resolve().parent / "data" / "stats_reports" / "01_pca_overview.csv")
    _pca_pc1_pct, _pca_pc2_pct = 39.6, 16.6  # R defaults (all genes, global PCA)
    if os.path.exists(_pca_ov_csv):
        try:
            _pca_ov = pd.read_csv(_pca_ov_csv)
            _row = _pca_ov[
                (_pca_ov["dataset"] == "ALL GENES") &
                (_pca_ov["analysis_scope"] == "global_appendage_pca")
            ]
            if not _row.empty:
                _pca_pc1_pct = round(float(_row["pc1_variance"].iloc[0]), 1)
                _pca_pc2_pct = round(float(_row["pc2_variance"].iloc[0]), 1)
        except Exception:
            pass

    col_pca, col_corr = st.columns(2)

    with col_pca:
        if pc_group.empty:
            st.info(
                "Could not find enough samples matching patterns like Ant_VF1, Leg_Vm2, P_MF3.\n"
                "Check your column names if you expect 9 groups."
            )
        else:
            fig_pca = px.scatter(
                pc_group, x="PC1", y="PC2", color="State_short",
                color_discrete_map=STATE_COLORS,
                hover_data={"Tissue": True, "State": True, "State_short": False},
            )
            fig_pca.update_traces(mode="markers", marker=dict(size=12), showlegend=False)
            x_min = pc_group["PC1"].min()
            x_max = pc_group["PC1"].max()
            y_min = pc_group["PC2"].min()
            y_max = pc_group["PC2"].max()
            pad_x = max((x_max - x_min) * 0.2, 0.5)
            pad_y = max((y_max - y_min) * 0.2, 0.5)
            x_range = [x_min - pad_x, x_max + pad_x]
            y_range = [y_min - pad_y, y_max + pad_y]
            fig_pca.update_xaxes(range=x_range, showgrid=False, zeroline=False,
                                 title=f"PC1 ({_pca_pc1_pct}%)")
            fig_pca.update_yaxes(range=y_range, showgrid=False, zeroline=False,
                                 scaleanchor="x", scaleratio=1,
                                 title=f"PC2 ({_pca_pc2_pct}%)")
            for tissue in pc_group["Tissue"].unique():
                sub = pc_group[pc_group["Tissue"] == tissue]
                if sub.empty:
                    continue
                cx = sub["PC1"].mean()
                cy = sub["PC2"].mean()
                dx = (sub["PC1"] - cx).abs().max()
                dy = (sub["PC2"] - cy).abs().max()
                r = max(dx, dy)
                if r == 0:
                    r = max(x_max - x_min, y_max - y_min) / 10.0 or 0.5
                r *= 1.3
                fig_pca.add_shape(type="circle", xref="x", yref="y",
                                  x0=cx - r, x1=cx + r, y0=cy - r, y1=cy + r,
                                  line=dict(width=0), fillcolor="rgba(128,128,128,0.15)",
                                  layer="below")
                fig_pca.add_annotation(x=cx, y=cy - r * 1.05, text=tissue, showarrow=False,
                                       font=dict(size=24, color="black"), align="center",
                                       yanchor="top")
            fig_pca.add_shape(type="line", x0=0, x1=0, y0=y_range[0], y1=y_range[1],
                              line=dict(color="black", width=1, dash="dot"))
            fig_pca.add_shape(type="line", x0=x_range[0], x1=x_range[1], y0=0, y1=0,
                              line=dict(color="black", width=1, dash="dot"))
            fig_pca = t1_style_pca_fig(fig_pca, height=600, width=600)
            st.plotly_chart(fig_pca, use_container_width=False)

        if not pca_loadings.empty:
            with st.expander("Show PCA loadings table"):
                pc_options = [c for c in pca_loadings.columns if c.startswith("PC")]
                pc_choice = st.selectbox("Choose principal component", pc_options, index=0,
                                         key="t1_pc_choice")
                tbl = pca_loadings.copy()
                tbl["abs_loading"] = tbl[pc_choice].abs()
                tbl = tbl.sort_values("abs_loading", ascending=False)
                cols = ["JoinKey"]
                if "Name" in tbl.columns:
                    cols.append("Name")
                cols += [pc_choice, "abs_loading"]
                st.dataframe(tbl[cols], use_container_width=True)

    with col_corr:
        st.subheader("Correlation")
        st.caption("Manuscript: **Figure 2B**")
        if corr_matrix.empty:
            st.info(
                "Could not compute correlation matrix. "
                "Expected columns like Ant_VF1, Leg_MF1, P_Vm1 in the counts table."
            )
        else:
            fig_corr = t1_build_correlation_figure(corr_matrix)
            st.plotly_chart(fig_corr, use_container_width=False)

    st.markdown("---")
    st.markdown("---")
    st.subheader("Chemoreceptors-only PCA and Pearson correlation")
    st.caption("Manuscript: **Figure 2E** (PCA) · **Figure 2F** (correlation)")

    norm_counts_chemo = t1_filter_chemoreceptors(norm_counts)

    if norm_counts_chemo.empty:
        st.info("No chemoreceptor genes detected in the counts table (OR/GR/IR/PPK/TRP).")
    else:
        col_pca2, col_corr2 = st.columns(2)

        with col_pca2:
            pc_group_chemo = t1_compute_pca_by_group(norm_counts_chemo)
            if pc_group_chemo.empty:
                st.info("Could not compute PCA for chemoreceptors (check library columns).")
            else:
                fig_pca_chemo = px.scatter(
                    pc_group_chemo, x="PC1", y="PC2", color="State_short",
                    color_discrete_map=STATE_COLORS,
                    hover_data={"Tissue": True, "State": True, "State_short": False},
                )
                fig_pca_chemo.update_traces(mode="markers", marker=dict(size=12), showlegend=False)

                # Axis padding to keep circles/labels inside the frame
                x_min = pc_group_chemo["PC1"].min()
                x_max = pc_group_chemo["PC1"].max()
                y_min = pc_group_chemo["PC2"].min()
                y_max = pc_group_chemo["PC2"].max()
                pad_x = max((x_max - x_min) * 0.2, 0.5)
                pad_y = max((y_max - y_min) * 0.2, 0.5)
                x_range = [x_min - pad_x, x_max + pad_x]
                y_range = [y_min - pad_y, y_max + pad_y]

                # Load chemo-specific % variance from R
                _chemo_pc1_pct, _chemo_pc2_pct = 36.1, 30.3
                if os.path.exists(_pca_ov_csv):
                    try:
                        _pca_ov2 = pd.read_csv(_pca_ov_csv)
                        _row2 = _pca_ov2[
                            (_pca_ov2["dataset"].str.contains("CHEMOSENSORY", na=False)) &
                            (_pca_ov2["analysis_scope"] == "global_appendage_pca")
                        ]
                        if not _row2.empty:
                            _chemo_pc1_pct = round(float(_row2["pc1_variance"].iloc[0]), 1)
                            _chemo_pc2_pct = round(float(_row2["pc2_variance"].iloc[0]), 1)
                    except Exception:
                        pass

                fig_pca_chemo.update_xaxes(range=x_range, showgrid=False, zeroline=False,
                                           title=f"PC1 ({_chemo_pc1_pct}%)")
                fig_pca_chemo.update_yaxes(
                    range=y_range, showgrid=False, zeroline=False,
                    scaleanchor="x", scaleratio=1, title=f"PC2 ({_chemo_pc2_pct}%)"
                )

                # Grey "clouds" per tissue (same style as the global PCA)
                for tissue in pc_group_chemo["Tissue"].unique():
                    sub = pc_group_chemo[pc_group_chemo["Tissue"] == tissue]
                    if sub.empty:
                        continue
                    cx = sub["PC1"].mean()
                    cy = sub["PC2"].mean()
                    dx = (sub["PC1"] - cx).abs().max()
                    dy = (sub["PC2"] - cy).abs().max()
                    r = max(dx, dy)
                    if r == 0:
                        r = max(x_max - x_min, y_max - y_min) / 10.0 or 0.5
                    r *= 1.3
                    fig_pca_chemo.add_shape(
                        type="circle", xref="x", yref="y",
                        x0=cx - r, x1=cx + r, y0=cy - r, y1=cy + r,
                        line=dict(width=0),
                        fillcolor="rgba(128,128,128,0.15)",
                        layer="below",
                    )
                    fig_pca_chemo.add_annotation(
                        x=cx, y=cy - r * 1.05, text=tissue, showarrow=False,
                        font=dict(size=24, color="black"), align="center",
                    )

                # Reference lines at 0
                fig_pca_chemo.add_shape(
                    type="line", x0=0, x1=0, y0=y_range[0], y1=y_range[1],
                    line=dict(color="black", width=1, dash="dot")
                )
                fig_pca_chemo.add_shape(
                    type="line", x0=x_range[0], x1=x_range[1], y0=0, y1=0,
                    line=dict(color="black", width=1, dash="dot")
                )

                fig_pca_chemo = t1_style_pca_fig(fig_pca_chemo, height=600, width=600)
                st.plotly_chart(fig_pca_chemo, use_container_width=False)

                pca_loadings_chemo = t1_compute_pca_loadings(norm_counts_chemo, n_components=5)
                if not pca_loadings_chemo.empty:
                    with st.expander("Show PCA loadings table (chemoreceptors)"):
                        pc_options = [c for c in pca_loadings_chemo.columns if c.startswith("PC")]
                        pc_choice = st.selectbox(
                            "Choose principal component",
                            pc_options,
                            index=0,
                            key="t1_pca_loading_pc_chemo",
                        )
                        n_top = st.slider(
                            "Top loadings to show",
                            min_value=10,
                            max_value=200,
                            value=50,
                            step=10,
                            key="t1_pca_loading_n_chemo",
                        )
                        tmp = pca_loadings_chemo.sort_values(
                            pc_choice,
                            key=lambda s: s.abs(),
                            ascending=False,
                        ).head(n_top)
                        st.dataframe(tmp, use_container_width=True, hide_index=True)

        with col_corr2:
            group_means_chemo = t1_compute_group_means_for_corr(norm_counts_chemo)
            corr_matrix_chemo = t1_compute_correlation_matrix(group_means_chemo)
            if corr_matrix_chemo.empty:
                st.info("Could not compute Pearson correlation (chemoreceptors).")
            else:
                fig_corr_chemo = t1_build_correlation_figure(corr_matrix_chemo)
                st.plotly_chart(fig_corr_chemo, use_container_width=False)

    st.markdown("---")
    st.subheader("Multivariate variance partition (pie)")
    st.caption("Manuscript: **Figure 2C** (all genes) · **Figure 2G** (chemosensory genes)")

    # Load precomputed R values (vegan varpart / partial RDA — matches manuscript)
    _vp_csv = str(pathlib.Path(__file__).resolve().parent / "data" / "stats_reports" / "07_variation_partitioning.csv")

    def _load_vp_from_csv(csv_path, dataset_filter):
        if not os.path.exists(csv_path):
            return pd.DataFrame()
        df = pd.read_csv(csv_path)
        mask = df["dataset"].str.strip() == dataset_filter
        sub = df[mask & df["model"].str.contains("3-way", na=False)].copy()
        if sub.empty:
            return pd.DataFrame()
        comp_map = {
            "Appendage unique": "Appendage identity",
            "Sex unique": "Sex",
            "Mating unique": "Mating status",
            "Shared/confounded": "Shared / confounded",
            "Unexplained": "Unexplained",
        }
        order = ["Appendage identity", "Sex", "Mating status",
                 "Shared / confounded", "Unexplained"]
        rows = []
        for _, r in sub.iterrows():
            comp = comp_map.get(r["component"].strip(), r["component"].strip())
            frac = max(0.0, float(r["adj_r2_fraction"]))
            if comp == "Shared / confounded" and frac <= 0:
                continue          # it is exactly zero here; an empty wedge helps nobody
            rows.append({"Component": comp, "Fraction": frac})
        out = pd.DataFrame(rows)
        out["_o"] = out["Component"].map({c: i for i, c in enumerate(order)})
        return out.sort_values("_o").drop(columns="_o").reset_index(drop=True)

    pie_all   = _load_vp_from_csv(_vp_csv, "ALL GENES")
    pie_chemo = _load_vp_from_csv(_vp_csv, "CHEMOSENSORY GENES (Name non-empty)")

    # Fall back to Python computation if CSV not found
    if pie_all.empty:
        pie_all = t1_variance_partition_pie(norm_counts)
    if pie_chemo.empty:
        pie_chemo = t1_variance_partition_pie(norm_counts_chemo) if not norm_counts_chemo.empty else pd.DataFrame()

    col_vp1, col_vp2 = st.columns(2)

    def _render_vp_pie(pie_df, title_caption, key_suffix):
        st.caption(title_caption)
        if pie_df.empty:
            st.info("Could not compute variance partition.")
            return
        fig = px.pie(pie_df, names="Component", values="Fraction",
                     color="Component", color_discrete_map=VP_PIE_COLORS)
        fig.update_traces(textinfo="percent+label",
                          hovertemplate="%{label}: %{percent:.1%}<extra></extra>")
        fig.update_layout(height=420, width=520, showlegend=False)
        st.plotly_chart(fig, use_container_width=False, key=f"vp_pie_{key_suffix}")

    with col_vp1:
        _render_vp_pie(pie_all, "All genes (n = 24,608; Ezekiel adj. R²; RDA, vegan). "
                                "Three-way model: appendage, sex and mating status.", "all")
    with col_vp2:
        _render_vp_pie(pie_chemo, "Chemosensory genes only. "
                                  "Three-way model: appendage, sex and mating status.", "chemo")

# =============================================================================
# === APP_C1 FUNCTIONS ===
# =============================================================================

@st.cache_data(show_spinner=False)
def c1_load_all(ant_p, ant_leg, leg_p, go_path):
    return load_csv(ant_p), load_csv(ant_leg), load_csv(leg_p), load_csv(go_path)


def c1_extract_chemo_tag(s: pd.Series) -> pd.Series:
    s = s.astype(str)
    hit = pd.Series(index=s.index, dtype="object")
    for tag in CHEMO_TAGS:
        mask = s.str.contains(fr"\b{tag}[0-9A-Za-z._-]*\b", case=False, regex=True)
        hit.loc[mask] = tag.upper()
    return hit


def c1_annotate_with_go(df: pd.DataFrame, go_map: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["JoinKey"] = make_joinkey(out)
    out = out.merge(go_map, how="left", left_on="JoinKey", right_on="Gene")
    out["GO_Domain"] = clean_go_domain(out["GO_Domain"])
    out.loc[out["GO_Name"].isna() | (out["GO_Name"] == ""), "GO_Name"] = "Unknown"
    return out


def c1_prep_volcano_df(df, go_map, left_name, right_name, padj_thr, lfc_thr, strong_lfc=2.0):
    v = c1_annotate_with_go(df, go_map).copy()
    v["padj"] = pd.to_numeric(v.get("padj", 1.0), errors="coerce").fillna(1.0)
    v["log2FoldChange"] = pd.to_numeric(v.get("log2FoldChange", 0.0), errors="coerce").fillna(0.0)
    v["padj_safe"] = np.maximum(v["padj"].astype(float).values, np.finfo(float).tiny)
    v["-log10"] = -np.log10(v["padj_safe"])
    abs_lfc = np.abs(v["log2FoldChange"].astype(float))
    v["is_sig"] = (v["padj"].astype(float) < padj_thr) & (abs_lfc >= lfc_thr)
    v["Direction"] = np.where(
        (v["is_sig"]) & (v["log2FoldChange"].astype(float) > 0), f"{left_name} up",
        np.where(
            (v["is_sig"]) & (v["log2FoldChange"].astype(float) < 0),
            f"{right_name} up", "Not sig"
        )
    )
    name_col = v["Name"] if "Name" in v.columns else v["JoinKey"]
    v["ChemoName"] = c1_extract_chemo_tag(name_col).fillna(c1_extract_chemo_tag(v["JoinKey"]))
    v["is_chemo"] = v["ChemoName"].notna()

    def colorkey(row):
        if not row["is_sig"]:
            return "Not Significant"
        return row["GO_Domain"] if pd.notna(row["GO_Domain"]) else "Unknown"

    v["ColorKey"] = v.apply(colorkey, axis=1)
    v["alpha_pt"] = np.where(
        (v["is_sig"]) & (abs_lfc >= strong_lfc), 0.95,
        np.where(v["is_sig"], 0.65, 0.25)
    )
    v.attrs["labels"] = {"left": left_name, "right": right_name}
    return v


def c1_make_percent_table(vdf: pd.DataFrame) -> pd.DataFrame:
    labs = vdf.attrs.get("labels", {"left": "Left", "right": "Right"})
    d = vdf[vdf["is_sig"]].copy()

    d["Dir"] = np.where(
        d["log2FoldChange"] > 0,
        f"{labs['left']} up",
        f"{labs['right']} up"
    )

    d["Tissue"] = d["Dir"].str.replace(" up", " biased", regex=False)
    d["GO_Domain"] = pd.Categorical(
        d["GO_Domain"].fillna("Unknown"),
        ["MF", "CC", "BP", "Unknown"],
        ordered=True
    )

    g = d.groupby(["Dir", "Tissue", "GO_Domain"], dropna=False).size().reset_index(name="N")
    g["Percent"] = g.groupby("Dir")["N"].transform(lambda x: 100 * x / x.sum())
    g["GO_Domain_full"] = g["GO_Domain"].map(DOMAIN_FULL).fillna("Unknown")
    g["TextN"] = "N=" + g["N"].astype(int).astype(str)
    return g


def c1_make_percent_table_chemo(vdf: pd.DataFrame, show_chemo: bool) -> pd.DataFrame:
    """Like c1_make_percent_table but with optional Chemosensory category for DEGs with is_chemo==True."""
    if not show_chemo:
        return c1_make_percent_table(vdf)
    labs = vdf.attrs.get("labels", {"left": "Left", "right": "Right"})
    d = vdf[vdf["is_sig"]].copy()
    d["Dir"] = np.where(
        d["log2FoldChange"] > 0,
        f"{labs['left']} up",
        f"{labs['right']} up"
    )
    d["Tissue"] = d["Dir"].str.replace(" up", " biased", regex=False)
    is_chemo_mask = d["is_chemo"].fillna(False).astype(bool) if "is_chemo" in d.columns else pd.Series(False, index=d.index)
    d["GO_Domain"] = np.where(
        is_chemo_mask,
        "Chemosensory",
        d["GO_Domain"].fillna("Unknown").astype(str)
    )
    cat_order = ["Chemosensory", "MF", "CC", "BP", "Unknown"]
    d["GO_Domain"] = pd.Categorical(d["GO_Domain"], cat_order, ordered=True)
    g = d.groupby(["Dir", "Tissue", "GO_Domain"], dropna=False, observed=True).size().reset_index(name="N")
    g["Percent"] = g.groupby("Dir")["N"].transform(lambda x: 100 * x / x.sum())
    domain_full_chemo = {**DOMAIN_FULL, "Chemosensory": "Chemosensory"}
    g["GO_Domain_full"] = g["GO_Domain"].astype(str).map(domain_full_chemo).fillna("Unknown")
    g["TextN"] = "N=" + g["N"].astype(int).astype(str)
    return g


def c1_volcano_title(vdf: pd.DataFrame) -> str:
    labs = vdf.attrs.get("labels", {})
    left = labs.get("left", "Left")
    right = labs.get("right", "Right")
    n_left_side = (vdf["Direction"] == f"{right} up").sum()
    n_right_side = (vdf["Direction"] == f"{left} up").sum()
    return f"{right} : {n_left_side}    {left} : {n_right_side}"


@st.cache_data(show_spinner=False)
def c1_load_norm_means(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["JoinKey"] = make_joinkey(df)
    out = df[["JoinKey"]].copy()
    for tissue, prefixes in TISSUE_PREFIXES.items():
        cols = [c for c in df.columns if any(c.startswith(p) for p in prefixes)]
        if cols:
            out[tissue + "_mean"] = df[cols].mean(axis=1)
    out = out.drop_duplicates(subset=["JoinKey"])
    return out


@st.cache_data(show_spinner=False)
def c1_build_name_map(antp, antleg, legp):
    parts = []
    for df in (antp, antleg, legp):
        if "Name" in df.columns:
            tmp = df.copy()
            tmp["JoinKey"] = make_joinkey(tmp)
            parts.append(tmp[["JoinKey", "Name"]])
    if not parts:
        return pd.DataFrame(columns=["JoinKey", "Name"])
    nm = pd.concat(parts, ignore_index=True)
    nm = nm.dropna(subset=["Name"])
    nm = nm.drop_duplicates(subset=["JoinKey"])
    return nm


def c1_build_go_group_detail(vdf, dir_label, go_domain, norm_means):
    sub = vdf[
        (vdf["is_sig"]) & (vdf["Direction"] == dir_label) & (vdf["GO_Domain"] == go_domain)
    ].copy()
    if sub.empty:
        return sub
    cols_keep = [c for c in ["Gene", "Name", "JoinKey", "log2FoldChange", "padj",
                              "GO_Name", "GO_Domain", "Direction"] if c in sub.columns]
    sub = sub[cols_keep]
    if norm_means is not None:
        sub = sub.merge(norm_means, on="JoinKey", how="left")
    if "JoinKey" in sub.columns:
        if "Gene" in sub.columns:
            sub = sub.drop(columns=["Gene"])
        sub = sub.rename(columns={"JoinKey": "Gene"})
    base_order = ["Gene", "Name", "log2FoldChange", "padj", "GO_Name", "GO_Domain", "Direction"]
    other_cols = [c for c in sub.columns if c not in base_order]
    sub = sub[[c for c in base_order if c in sub.columns] + other_cols]
    if "padj" in sub.columns:
        sub["padj"] = pd.to_numeric(sub["padj"], errors="coerce")
        sub["padj"] = sub["padj"].map(
            lambda x: f"{x:.2e}" if x is not None and not pd.isna(x) else ""
        )
    return sub


def c1_find_hits(vdf: pd.DataFrame, query: str) -> pd.DataFrame:
    if not query or str(query).strip() == "":
        return vdf.iloc[0:0]
    q = str(query).strip()
    cols = [c for c in ["Gene", "Name", "JoinKey"] if c in vdf.columns]
    if not cols:
        return vdf.iloc[0:0]
    mask_exact = False
    for c in cols:
        m = vdf[c].astype(str).str.strip() == q
        mask_exact = (mask_exact | m) if isinstance(mask_exact, pd.Series) else m
    hits = vdf[mask_exact]
    if not hits.empty:
        return hits.iloc[[0]]
    pattern = rf"\b{re.escape(q)}\b"
    mask_word = False
    for c in cols:
        m = vdf[c].astype(str).str.contains(pattern, case=False, regex=True, na=False)
        mask_word = (mask_word | m) if isinstance(mask_word, pd.Series) else m
    hits = vdf[mask_word]
    if not hits.empty:
        return hits.iloc[[0]]
    qlow = q.lower()
    mask_cont = False
    for c in cols:
        m = vdf[c].astype(str).str.contains(qlow, case=False, na=False)
        mask_cont = (mask_cont | m) if isinstance(mask_cont, pd.Series) else m
    hits = vdf[mask_cont]
    if hits.empty:
        return hits
    return hits.iloc[[0]]


# The published Figure 3A colours chemosensory genes red, every other
# significant gene light blue, and the non-significant cloud pale grey, with
# the chemosensory points drawn last so they sit on top. The GO-domain
# colouring the app used to default to is kept behind a tick.
ALT_VOLCANO_COLS = {
    "Chemosensory": "#A60000",
    "All other genes": "#9ECAE1",
    "Not significant": "#E4E4E4",
}

def c1_fig_volcano(vdf, padj_thr, lfc_thr, overlay_chemo, highlight_query,
                   color_by_go=False):
    if not color_by_go:
        chemo = chemo_transcript_set()
        key = (vdf["Gene"] if "Gene" in vdf.columns else vdf["JoinKey"]).astype(str).str.strip()
        is_ch = key.isin(chemo) if chemo else vdf["is_chemo"].fillna(False).astype(bool)
        vdf = vdf.copy()
        vdf["ColorKey"] = np.where(vdf["is_sig"] & is_ch, "Chemosensory",
                           np.where(vdf["is_sig"], "All other genes", "Not significant"))
        # chemosensory last so the genes this paper is about are never buried
        order = {"Not significant": 0, "All other genes": 1, "Chemosensory": 2}
        vdf = vdf.assign(_z=vdf["ColorKey"].map(order)).sort_values("_z")
        palette = ALT_VOLCANO_COLS
    else:
        palette = PAL_GO
    fig = px.scatter(
        vdf, x="log2FoldChange", y="-log10", color="ColorKey",
        color_discrete_map=palette,
        hover_data={
            "JoinKey": True,
            "Name": True if "Name" in vdf.columns else False,
            "padj": ":.2e", "log2FoldChange": ":.3f",
            "GO_Name": True, "GO_Domain": True,
        },
        render_mode="webgl", height=600,
    )
    fig.update_traces(marker=dict(size=5, opacity=0.6), selector=dict(mode="markers"))
    fig.add_hline(y=-np.log10(padj_thr), line=dict(dash="dash", width=0.5))
    fig.add_vline(x=-lfc_thr, line=dict(dash="dash", width=0.5))
    fig.add_vline(x=lfc_thr, line=dict(dash="dash", width=0.5))
    fig.update_layout(
        legend_title_text="",
        xaxis_title="Log2(FoldChange)",
        yaxis_title="-log10(<i>padj</i>)",
        uirevision="volcano"
    )

    if overlay_chemo:
        chemo = vdf[(vdf["is_sig"]) & (vdf["is_chemo"])].copy()
        if not chemo.empty:
            fig.add_trace(go.Scattergl(
                x=chemo["log2FoldChange"], y=chemo["-log10"],
                mode="markers", name="Labeled (chemo)",
                marker=dict(color="#8B0000", size=7, opacity=0.9),
                hovertext=chemo["JoinKey"], hoverinfo="text"
            ))

    hits = c1_find_hits(vdf, highlight_query)
    if not hits.empty:
        fig.add_trace(go.Scattergl(
            x=hits["log2FoldChange"], y=hits["-log10"],
            mode="markers", name=f"Search: {highlight_query}",
            marker=dict(color="red", size=10, line=dict(color="white", width=1.5), opacity=1.0),
            hovertext=(hits["Gene"] if "Gene" in hits.columns else hits["JoinKey"]),
            hoverinfo="text"
        ))
    return fig

@st.cache_data(show_spinner=False)
def c1_load_chemo_tables(chemo_dir: str):
    fam_order = ["Or", "Gr", "Ir", "Obp", "Csp", "Ppk", "Trp"]
    family_tables = {}
    for fam in fam_order:
        csv_path = os.path.join(chemo_dir, f"chemo_{fam}.csv")
        if not os.path.exists(csv_path):
            continue
        df = pd.read_csv(csv_path)
        family_tables[fam] = df
    if not family_tables:
        raise FileNotFoundError(
            f"No chemo_*.csv files found in {chemo_dir}. "
            f"Expected files like chemo_Or.csv, chemo_Gr.csv, etc."
        )
    return family_tables


def c1_fig_chemo_pie(summary, fam_key: str, tissue_key: str, show_legend: bool) -> go.Figure:
    counts = summary.get(fam_key, {}).get(tissue_key, {})
    if not counts:
        return go.Figure()
    df = pd.DataFrame({"Class": list(counts.keys()), "Count": list(counts.values())})
    df = df[df["Count"] > 0]
    if df.empty:
        return go.Figure()
    fig = px.pie(df, names="Class", values="Count", color="Class",
                 color_discrete_map=CHEMO_PIE_COLORS, hole=0)
    fig.update_traces(textinfo="value",
                      hovertemplate="%{label}: %{value} genes (%{percent:.1%})<extra></extra>",
                      showlegend=False)
    base_layout = dict(margin=dict(l=5, r=5, t=10, b=60), height=330)
    if show_legend:
        legend_items = [("Appendage-specific", "Specific"), ("Appendage-biased", "Biased"),
                        ("Expressed", "Expressed"), ("Not expressed", "Not expressed")]
        for key, label in legend_items:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                     marker=dict(symbol="square", size=10,
                                                 color=CHEMO_PIE_COLORS[key]),
                                     name=label, showlegend=True, hoverinfo="skip"))
        base_layout["showlegend"] = True
        base_layout["legend"] = dict(orientation="h", yanchor="top", y=-0.05,
                                     xanchor="center", x=0.5, font=dict(size=10))
    else:
        base_layout["showlegend"] = False
    fig.update_layout(**base_layout)
    return fig


def c1_up_keys(df, sign, padj_thr, lfc_thr):
    df = df.copy()
    df["padj"] = pd.to_numeric(df.get("padj", 1.0), errors="coerce").fillna(1.0)
    df["log2FoldChange"] = pd.to_numeric(df.get("log2FoldChange", 0.0), errors="coerce").fillna(0.0)
    if sign == +1:
        sub = df[(df["padj"] < padj_thr) & (df["log2FoldChange"] > lfc_thr)]
    else:
        sub = df[(df["padj"] < padj_thr) & (df["log2FoldChange"] < -lfc_thr)]
    return set(make_joinkey(sub).dropna().unique())


def c1_go_name_table_overlap(keys, go_map, domain, drop_unknown=False):
    gm = go_map.copy()
    gm["GO_Domain"] = clean_go_domain(gm["GO_Domain"])
    annot = gm[gm["Gene"].astype(str).isin(keys) & (gm["GO_Domain"] == domain)].copy()
    if drop_unknown:
        annot = annot[annot["GO_Name"].astype(str) != "Unknown"]
    return (
        annot.groupby("GO_Name").size().reset_index(name="N")
        .sort_values(["N", "GO_Name"], ascending=[False, True])
    )


def c1_fig_go_names(tbl, domain, x_max=None):
    if tbl.empty:
        return go.Figure()
    tbl = tbl.copy()
    if x_max is None:
        x_max = int(tbl["N"].max())
    label = DOMAIN_FULL.get(domain, domain).lower()
    fig = px.bar(tbl, x="N", y="GO_Name", orientation="h",
                 title=f"Overlap GO Names in {label}")
    fig.update_traces(marker_color=GO_COLS.get(domain, "#CCCCCC"))
    fig.update_layout(yaxis=dict(tickfont=dict(size=10)),
                      xaxis_title="Gene count", yaxis_title="GO Name")
    fig.update_xaxes(range=[0, max(1, int(x_max * 1.1))])
    return fig


def c1_domain_keys_for_tissue(tissue, domain, antp_v, antleg_v, legp_v):
    if tissue == "Antenna":
        ma = (antp_v["is_sig"]) & (antp_v["Direction"].astype(str).str.contains("Antenna", na=False)) & (antp_v["GO_Domain"] == domain)
        mb = (antleg_v["is_sig"]) & (antleg_v["Direction"].astype(str).str.contains("Antenna", na=False)) & (antleg_v["GO_Domain"] == domain)
        set1 = set(antp_v.loc[ma, "JoinKey"].astype(str))
        set2 = set(antleg_v.loc[mb, "JoinKey"].astype(str))
        label_left, label_right = "Maxillary palp", "Tarsi"
    elif tissue == "Maxillary palp":
        ma = (antp_v["is_sig"]) & (antp_v["Direction"].astype(str).str.contains("Maxillary palp", na=False)) & (antp_v["GO_Domain"] == domain)
        mb = (legp_v["is_sig"]) & (legp_v["Direction"].astype(str).str.contains("Maxillary palp", na=False)) & (legp_v["GO_Domain"] == domain)
        set1 = set(antp_v.loc[ma, "JoinKey"].astype(str))
        set2 = set(legp_v.loc[mb, "JoinKey"].astype(str))
        label_left, label_right = "Antenna", "Tarsi"
    else:
        ma = (antleg_v["is_sig"]) & (antleg_v["Direction"].astype(str).str.contains("Tarsi", na=False)) & (antleg_v["GO_Domain"] == domain)
        mb = (legp_v["is_sig"]) & (legp_v["Direction"].astype(str).str.contains("Tarsi", na=False)) & (legp_v["GO_Domain"] == domain)
        set1 = set(antleg_v.loc[ma, "JoinKey"].astype(str))
        set2 = set(legp_v.loc[mb, "JoinKey"].astype(str))
        label_left, label_right = "Antenna", "Maxillary palp"
    overlap = set1 & set2
    return label_left, label_right, set1, set2, overlap


def c1_fig_venn_two_sets(label_left, label_right, n_left, n_right, n_overlap, title, domain_color):
    fig = go.Figure()
    r = 1.0
    c1_coord = (-0.5, 0)
    c2_coord = (0.5, 0)
    fill_rgba = hex_to_rgba(domain_color, 0.35)
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c1_coord[0]-r, x1=c1_coord[0]+r, y0=c1_coord[1]-r, y1=c1_coord[1]+r,
                  line=dict(color="rgba(0,0,0,0)"), fillcolor=fill_rgba)
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c2_coord[0]-r, x1=c2_coord[0]+r, y0=c2_coord[1]-r, y1=c2_coord[1]+r,
                  line=dict(color="rgba(0,0,0,0)"), fillcolor=fill_rgba)
    fig.add_annotation(x=c1_coord[0]-0.35, y=0, text=str(n_left), showarrow=False, font=dict(size=12))
    fig.add_annotation(x=c2_coord[0]+0.35, y=0, text=str(n_right), showarrow=False, font=dict(size=12))
    fig.add_annotation(x=0, y=0.15, text=str(n_overlap), showarrow=False, font=dict(size=14))
    fig.add_annotation(x=c1_coord[0], y=-0.9, text=label_left, showarrow=False, font=dict(size=10))
    fig.add_annotation(x=c2_coord[0], y=-0.9, text=label_right, showarrow=False, font=dict(size=10))
    fig.update_xaxes(visible=False, range=[-2, 2])
    fig.update_yaxes(visible=False, range=[-1.5, 1.5])
    fig.update_layout(title=title, showlegend=False, margin=dict(l=10, r=10, t=40, b=10), height=220)
    return fig


def c1_compute_tissue_classes(norm_means, AntP, AntLeg, LegP,
                               expr_thr=10.0, lfc_thr=1.0, padj_thr=0.001):
    if norm_means is None:
        return pd.DataFrame(columns=["JoinKey", "Antenna_Class", "Palp_Class", "Tarsi_Class"])
    base = norm_means.copy()
    antp = AntP.copy(); antp["JoinKey"] = make_joinkey(antp)
    antleg = AntLeg.copy(); antleg["JoinKey"] = make_joinkey(antleg)
    legp = LegP.copy(); legp["JoinKey"] = make_joinkey(legp)
    for df, lfc_col, padj_col in [
        (antp, "LFC_AntP", "padj_AntP"),
        (antleg, "LFC_AntLeg", "padj_AntLeg"),
        (legp, "LFC_LegP", "padj_LegP"),
    ]:
        df["log2FoldChange"] = pd.to_numeric(df.get("log2FoldChange", 0.0), errors="coerce").fillna(0.0)
        df["padj"] = pd.to_numeric(df.get("padj", 1.0), errors="coerce").fillna(1.0)
        base = base.merge(
            df[["JoinKey", "log2FoldChange", "padj"]].rename(
                columns={"log2FoldChange": lfc_col, "padj": padj_col}
            ), on="JoinKey", how="left",
        )
    for c in ["Antenna_mean", "Maxillary palp_mean", "Tarsi_mean"]:
        if c in base.columns:
            base[c] = base[c].fillna(0.0)
    exprA = base.get("Antenna_mean", pd.Series(0.0, index=base.index))
    exprP = base.get("Maxillary palp_mean", pd.Series(0.0, index=base.index))
    exprT = base.get("Tarsi_mean", pd.Series(0.0, index=base.index))
    LFC_AntP = base["LFC_AntP"]; padj_AntP = base["padj_AntP"]
    LFC_AntLeg = base["LFC_AntLeg"]; padj_AntLeg = base["padj_AntLeg"]
    LFC_LegP = base["LFC_LegP"]; padj_LegP = base["padj_LegP"]
    cond_exprA = exprA >= expr_thr; cond_exprP = exprP >= expr_thr; cond_exprT = exprT >= expr_thr
    cond_DE_A_vs_P = (LFC_AntP > lfc_thr) & (padj_AntP < padj_thr)
    cond_DE_P_vs_A = (LFC_AntP < -lfc_thr) & (padj_AntP < padj_thr)
    cond_DE_A_vs_T = (LFC_AntLeg > lfc_thr) & (padj_AntLeg < padj_thr)
    cond_DE_T_vs_A = (LFC_AntLeg < -lfc_thr) & (padj_AntLeg < padj_thr)
    cond_DE_T_vs_P = (LFC_LegP > lfc_thr) & (padj_LegP < padj_thr)
    cond_DE_P_vs_T = (LFC_LegP < -lfc_thr) & (padj_LegP < padj_thr)
    antenna_specific = cond_exprA & (exprP < expr_thr) & (exprT < expr_thr)
    antenna_biased = cond_exprA & ((exprP >= expr_thr) | (exprT >= expr_thr)) & cond_DE_A_vs_P & cond_DE_A_vs_T
    antenna_expressed = cond_exprA & ~antenna_specific & ~antenna_biased
    antenna_not = ~cond_exprA
    base["Antenna_Class"] = np.select(
        [antenna_specific, antenna_biased, antenna_expressed, antenna_not],
        ["Appendage-specific", "Appendage-biased", "Expressed", "Not expressed"], default="Not expressed")
    palp_specific = cond_exprP & (exprA < expr_thr) & (exprT < expr_thr)
    palp_biased = cond_exprP & ((exprA >= expr_thr) | (exprT >= expr_thr)) & cond_DE_P_vs_A & cond_DE_P_vs_T
    palp_expressed = cond_exprP & ~palp_specific & ~palp_biased
    palp_not = ~cond_exprP
    base["Palp_Class"] = np.select(
        [palp_specific, palp_biased, palp_expressed, palp_not],
        ["Appendage-specific", "Appendage-biased", "Expressed", "Not expressed"], default="Not expressed")
    tarsi_specific = cond_exprT & (exprA < expr_thr) & (exprP < expr_thr)
    tarsi_biased = cond_exprT & ((exprA >= expr_thr) | (exprP >= expr_thr)) & cond_DE_T_vs_A & cond_DE_T_vs_P
    tarsi_expressed = cond_exprT & ~tarsi_specific & ~tarsi_biased
    tarsi_not = ~cond_exprT
    base["Tarsi_Class"] = np.select(
        [tarsi_specific, tarsi_biased, tarsi_expressed, tarsi_not],
        ["Appendage-specific", "Appendage-biased", "Expressed", "Not expressed"], default="Not expressed")
    return base[["JoinKey", "Antenna_Class", "Palp_Class", "Tarsi_Class"]]


def c1_load_all_data(default_paths, padj_thr, lfc_thr, strong_lfc):
    BASE_DIR_C1 = default_paths["BASE_DIR"]
    ant_p_path = os.path.join(BASE_DIR_C1, default_paths["FILE_ANT_P"])
    ant_leg_path = os.path.join(BASE_DIR_C1, default_paths["FILE_ANT_LEG"])
    leg_p_path = os.path.join(BASE_DIR_C1, default_paths["FILE_LEG_P"])
    go_path = os.path.join(BASE_DIR_C1, default_paths["FILE_GO"])
    norm_path = os.path.join(BASE_DIR_C1, default_paths["FILE_NORM"])
    try:
        AntP, AntLeg, LegP, _GO_raw = c1_load_all(ant_p_path, ant_leg_path, leg_p_path, go_path)
    except Exception as e:
        st.error(f"Error loading DE or GO files: {e}")
        st.stop()
    GO_MAP = load_go_map(go_path)
    try:
        NORM_MEANS = c1_load_norm_means(norm_path)
    except Exception as e:
        st.warning(f"Could not load normalized counts from {norm_path}: {e}")
        NORM_MEANS = None
    NAME_MAP = c1_build_name_map(AntP, AntLeg, LegP)
    chemo_dir = default_paths.get("CHEMO_DIR", BASE_DIR_C1)
    try:
        CHEMO_TABLES = c1_load_chemo_tables(chemo_dir)
        CHEMO_AVAILABLE = True
        chemo_error = ""
    except Exception as e:
        CHEMO_TABLES = {}
        CHEMO_AVAILABLE = False
        chemo_error = str(e)
    def _neg_lfc(df):
        d = df.copy()
        if "log2FoldChange" in d.columns:
            d["log2FoldChange"] = -pd.to_numeric(d["log2FoldChange"], errors="coerce")
        return d
    antp_v = c1_prep_volcano_df(_neg_lfc(AntP), GO_MAP, "Maxillary palp", "Antenna", padj_thr, lfc_thr, strong_lfc)
    antleg_v = c1_prep_volcano_df(_neg_lfc(AntLeg), GO_MAP, "Tarsi", "Antenna", padj_thr, lfc_thr, strong_lfc)
    legp_v = c1_prep_volcano_df(LegP, GO_MAP, "Tarsi", "Maxillary palp", padj_thr, lfc_thr, strong_lfc)
    return {
        "paths": {"BASE_DIR": BASE_DIR_C1, "chemo_dir": chemo_dir},
        "de": {"AntP": AntP, "AntLeg": AntLeg, "LegP": LegP},
        "GO_MAP": GO_MAP, "NORM_MEANS": NORM_MEANS, "NAME_MAP": NAME_MAP,
        "CHEMO_TABLES": CHEMO_TABLES, "CHEMO_AVAILABLE": CHEMO_AVAILABLE,
        "CHEMO_ERROR": chemo_error,
        "VOLCANO": {"antp_v": antp_v, "antleg_v": antleg_v, "legp_v": legp_v},
    }


def c1_render_volcano_tab(volcano_dfs, padj_thr, lfc_thr, show_chemo, highlight_query):
    antp_v = volcano_dfs["antp_v"]
    antleg_v = volcano_dfs["antleg_v"]
    legp_v = volcano_dfs["legp_v"]
    st.subheader("Volcano plots")
    st.caption("Manuscript: **Figure 3A**")
    go_3a = st.checkbox(
        "Colour by GO domain instead", value=False, key="c1_3a_go",
        help=("The published panel colours chemosensory genes red and every other "
              "significant gene light blue."))
    cols = st.columns(3)
    with cols[0]:
        st.markdown(f"**{c1_volcano_title(antp_v)}**")
        fig = c1_fig_volcano(antp_v, padj_thr, lfc_thr, show_chemo, highlight_query, go_3a)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=antp_v.to_csv(index=False).encode(),
                           file_name="antenna_vs_maxillary_palp_volcano_table.csv",
                           mime="text/csv", key="c1_dl_antp")
    with cols[1]:
        st.markdown(f"**{c1_volcano_title(antleg_v)}**")
        fig = c1_fig_volcano(antleg_v, padj_thr, lfc_thr, show_chemo, highlight_query, go_3a)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=antleg_v.to_csv(index=False).encode(),
                           file_name="antenna_vs_tarsi_volcano_table.csv",
                           mime="text/csv", key="c1_dl_antleg")
    with cols[2]:
        st.markdown(f"**{c1_volcano_title(legp_v)}**")
        fig = c1_fig_volcano(legp_v, padj_thr, lfc_thr, show_chemo, highlight_query, go_3a)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=legp_v.to_csv(index=False).encode(),
                           file_name="tarsi_vs_maxillary_palp_volcano_table.csv",
                           mime="text/csv", key="c1_dl_legp")


# =============================================================================
# === FIGURE 3B: CHEMOSENSORY VERSUS THE REST ===
# The published panel used to split each DEG set by GO domain. Additional file
# 11 carries exactly one GO term per transcript, so that split described the
# annotation pipeline rather than the biology, and it was replaced. The panel
# now shows chemosensory against everything else, with an enrichment test
# against the genes DESeq2 actually tested in that contrast. The GO split is
# still reachable behind a tick, for anyone who wants it.
# =============================================================================
@st.cache_data(show_spinner=False)
def chemo_transcript_set():
    """The 567 chemosensory transcript models, from Additional file 37.

    Membership is taken from the published map rather than inferred from the
    Name string: name matching misses genes whose annotation carries no family
    tag, which undercounts every set (215 -> 202 for the antenna, for example)
    and shifts the odds ratios away from the published ones.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                        "supplementary",
                        "Additional_file_37_transcript_model_to_gene_locus.csv")
    if not os.path.exists(path):
        return set()
    d = pd.read_csv(path)
    return set(d["Transcript_ID"].astype(str).str.strip())


CHEMO_SPLIT_COLS = {"Chemosensory": "#A60000", "All other genes": "#BFBFBF"}


def c1_odds_ratio(k_set, n_set, k_bg, n_bg):
    """Odds ratio with a Woolf 95% interval and a normal-approximation p.

    The manuscript reports a one-sided Fisher exact test (Additional file 38);
    SciPy is not available in the deployed environment, so this is the standard
    log-odds approximation. It agrees with the published values to two decimals
    for every set in Figure 3.
    """
    import math
    a = k_set                      # chemosensory, in the set
    b = n_set - k_set              # other, in the set
    c = k_bg - k_set               # chemosensory, outside
    d = (n_bg - n_set) - c         # other, outside
    if min(a, b, c, d) <= 0:       # Haldane-Anscombe correction
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orv = (a * d) / (b * c)
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    lo, hi = math.exp(math.log(orv) - 1.96*se), math.exp(math.log(orv) + 1.96*se)
    z = abs(math.log(orv)) / se
    pv = math.erfc(z / math.sqrt(2))          # two-sided
    return orv, lo, hi, pv


def c1_stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def c1_chemo_composition(vdf):
    """One row per DEG set (direction) of one pairwise comparison."""
    labs = vdf.attrs.get("labels", {"left": "Left", "right": "Right"})
    chemo = chemo_transcript_set()
    key = vdf["Gene"].astype(str).str.strip() if "Gene" in vdf.columns \
        else vdf["JoinKey"].astype(str).str.strip()
    is_ch = key.isin(chemo) if chemo else vdf["is_chemo"].fillna(False).astype(bool)
    vdf = vdf.assign(_is_ch=is_ch.values)
    n_bg = len(vdf)
    k_bg = int(vdf["_is_ch"].sum())
    out = []
    sig = vdf[vdf["is_sig"]]
    for side, lab in ((+1, labs["left"]), (-1, labs["right"])):
        d = sig[sig["log2FoldChange"] > 0] if side > 0 else sig[sig["log2FoldChange"] < 0]
        n_set = len(d)
        if n_set == 0:
            continue
        k_set = int(d["_is_ch"].sum())
        orv, lo, hi, pv = c1_odds_ratio(k_set, n_set, k_bg, n_bg)
        out.append({"Set": f"{lab}-biased", "DEGs": n_set, "Chemosensory": k_set,
                    "All other genes": n_set - k_set,
                    "% chemosensory": round(100*k_set/n_set, 2),
                    "Odds ratio": round(orv, 2), "CI95 low": round(lo, 2),
                    "CI95 high": round(hi, 2), "p": pv, "sig": c1_stars(pv)})
    return pd.DataFrame(out)


def c1_render_composition_tab(volcano_dfs, show_go, norm_means):
    st.subheader("Chemosensory composition of each DEG set")
    st.caption("Manuscript: **Figure 3B**")
    show_go_here = st.checkbox(
        "Break down by GO domain instead", value=show_go, key="c1_3b_go",
        help=("The manuscript shows chemosensory against all other genes. Additional "
              "file 11 holds one GO term per transcript, so a GO split reflects the "
              "annotation pipeline more than the biology."))
    if show_go_here:
        c1_render_go_tab(volcano_dfs, norm_means, show_chemo_go=True)
        return
    st.markdown(
        "Bar height is the number of differentially expressed genes; the value inside "
        "each segment is its share of that set. Genes in the seven chemosensory "
        "families count as chemosensory whatever their GO domain. Above each bar is "
        "the odds ratio for chemosensory enrichment against the genes tested in that "
        "contrast, with significance stars."
    )
    cols = st.columns(3)
    tables = []
    for col, key, name in [
        (cols[0], "antp_v", "Antenna vs Maxillary palp"),
        (cols[1], "antleg_v", "Antenna vs Tarsi"),
        (cols[2], "legp_v", "Tarsi vs Maxillary palp"),
    ]:
        vdf = volcano_dfs[key]
        tbl = c1_chemo_composition(vdf)
        if tbl.empty:
            with col:
                st.info(f"No DEGs for {name}.")
            continue
        tbl.insert(0, "Comparison", name)
        tables.append(tbl)
        long = tbl.melt(id_vars=["Set", "DEGs", "Odds ratio", "sig"],
                        value_vars=["Chemosensory", "All other genes"],
                        var_name="Class", value_name="N")
        long["pct"] = long["N"] / long["DEGs"] * 100
        fig = px.bar(long, x="Set", y="N", color="Class",
                     color_discrete_map=CHEMO_SPLIT_COLS,
                     category_orders={"Class": ["Chemosensory", "All other genes"]},
                     text=long["pct"].map(lambda v: f"{v:.1f}%"))
        fig.update_traces(textposition="inside", insidetextanchor="middle",
                          textfont=dict(size=10, color="white"))
        for _, r in tbl.iterrows():
            fig.add_annotation(x=r["Set"], y=r["DEGs"], yshift=12, showarrow=False,
                               text=f"OR {r['Odds ratio']:.2f} {r['sig']}",
                               font=dict(size=10, color="#333333"))
        fig.update_layout(height=430, plot_bgcolor="white", bargap=0.45,
                          margin=dict(l=10, r=10, t=46, b=10),
                          yaxis_title="DEGs", xaxis_title=None,
                          legend=dict(orientation="h", yanchor="bottom", y=1.04,
                                      x=0, font=dict(size=9)),
                          title=dict(text=name, font=dict(size=12)))
        with col:
            st.plotly_chart(fig, use_container_width=True, key=f"c1_3b_{key}")
    if tables:
        allt = pd.concat(tables, ignore_index=True)
        allt["p"] = allt["p"].map(lambda v: f"{v:.3g}")
        st.dataframe(allt, use_container_width=True, hide_index=True)
        st.download_button("Download this table (CSV)",
                           allt.to_csv(index=False).encode("utf-8"),
                           file_name="Figure_3B_chemosensory_composition.csv",
                           mime="text/csv", key="c1_3b_dl")
        st.caption(
            "The manuscript's published values come from a one-sided Fisher exact "
            "test (Additional file 38); these are the log-odds approximation, which "
            "agrees to two decimals."
        )


# =============================================================================
# === FIGURE 3C: CONSENSUS APPENDAGE-BIASED SETS ===
# A gene is consensus appendage-biased when it is significantly higher in that
# appendage than in BOTH of the others, i.e. it sits in the intersection of the
# two pairwise contrasts that involve it. Targets: 494 / 539 / 4161.
#
# Note the sign convention differs from the sex and mating tables: in the
# appendage files a positive log2 fold change already means the first-named
# appendage is higher, so these are read unnegated.
# =============================================================================
CONSENSUS_SPEC = {
    "Antenna":        [("AntLeg", +1), ("AntP", +1)],
    "Maxillary palp": [("AntP", -1), ("LegP", -1)],
    "Tarsi":          [("LegP", +1), ("AntLeg", -1)],
}
CONSENSUS_PAIR_LABEL = {
    ("AntLeg", +1): "vs tarsi", ("AntP", +1): "vs maxillary palp",
    ("AntP", -1): "vs antenna", ("LegP", -1): "vs tarsi",
    ("LegP", +1): "vs maxillary palp", ("AntLeg", -1): "vs antenna",
}
CONSENSUS_COLS = {"Antenna": "#1C1B8D", "Maxillary palp": "#A60000", "Tarsi": "#2E7D32"}


def c1_consensus_sets(de_data, padj_thr, lfc_thr):
    out = {}
    for app, spec in CONSENSUS_SPEC.items():
        parts = []
        for key, sign in spec:
            d = de_data[key]
            padj = pd.to_numeric(d["padj"], errors="coerce")
            lfc = pd.to_numeric(d["log2FoldChange"], errors="coerce")
            m = (padj < padj_thr) & ((lfc >= lfc_thr) if sign > 0 else (lfc <= -lfc_thr))
            gcol = "Gene" if "Gene" in d.columns else d.columns[0]
            parts.append(set(d.loc[m, gcol].astype(str).str.strip()))
        out[app] = {"sets": parts, "labels": [CONSENSUS_PAIR_LABEL[t] for t in spec],
                    "consensus": parts[0] & parts[1]}
    return out


def c1_render_consensus_tab(de_data, padj_thr, lfc_thr, go_fallback):
    st.subheader("Consensus appendage-biased gene sets")
    st.caption("Manuscript: **Figure 3C**")
    show_go = st.checkbox("Show the GO-name bars and overlaps instead", value=False,
                          key="c1_3c_go")
    if show_go:
        go_fallback()
        return
    st.markdown(
        "A gene is consensus appendage-biased when it is significantly higher in "
        "that appendage than in **both** of the others, so it falls in the overlap "
        "of the two pairwise contrasts that involve it. These are the sets the GO "
        "over-representation in Figures S9–S11 is run on."
    )
    cons = c1_consensus_sets(de_data, padj_thr, lfc_thr)
    chemo = chemo_transcript_set()
    cols = st.columns(3)
    rows = []
    for i, (app, d) in enumerate(cons.items()):
        A, B = d["sets"]
        inter = d["consensus"]
        colour = CONSENSUS_COLS[app]
        with cols[i]:
            fig = go.Figure()
            for cx, lab in ((-0.55, d["labels"][0]), (0.55, d["labels"][1])):
                fig.add_shape(type="circle", x0=cx-1.15, x1=cx+1.15, y0=-1.15, y1=1.15,
                              fillcolor=colour, opacity=0.18, line=dict(color=colour, width=1))
            for x, txt in ((-1.15, len(A-B)), (1.15, len(B-A)), (0, len(inter))):
                fig.add_annotation(x=x, y=0, text=f"<b>{txt:,}</b>", showarrow=False,
                                   font=dict(size=15 if x == 0 else 12,
                                             color="#111111" if x == 0 else "#555555"))
            fig.add_annotation(x=-0.55, y=1.42, text=d["labels"][0], showarrow=False,
                               font=dict(size=9, color=colour))
            fig.add_annotation(x=0.55, y=-1.42, text=d["labels"][1], showarrow=False,
                               font=dict(size=9, color=colour))
            fig.update_layout(
                height=320, margin=dict(l=4, r=4, t=34, b=4),
                xaxis=dict(range=[-2.1, 2.1], visible=False),
                yaxis=dict(range=[-1.8, 1.8], visible=False, scaleanchor="x"),
                plot_bgcolor="white", showlegend=False,
                title=dict(text=f"{app}-biased", font=dict(size=12, color=colour), x=0.5))
            st.plotly_chart(fig, use_container_width=True, key=f"c1_3c_{app}")
        n_ch = len({g for g in inter if g in chemo})
        rows.append({"Appendage": app, "Consensus genes": len(inter),
                     "Chemosensory": n_ch,
                     "% chemosensory": round(100*n_ch/max(len(inter), 1), 2)})
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    pick = st.selectbox("List the consensus set for", list(cons.keys()), key="c1_3c_pick")
    genes = sorted(cons[pick]["consensus"])
    gdf = pd.DataFrame({"Gene": genes})
    gdf["Chemosensory"] = gdf["Gene"].isin(chemo)
    ex_table(gdf, f"cons_{pick}", f"Figure_3C_consensus_{pick.replace(' ','_')}.csv",
             height=320)


def c1_render_go_tab(volcano_dfs, norm_means, show_chemo_go=False):
    antp_v = volcano_dfs["antp_v"]
    antleg_v = volcano_dfs["antleg_v"]
    legp_v = volcano_dfs["legp_v"]
    go_col_map = GO_COLS_WITH_CHEMO if show_chemo_go else GO_COLS
    go_cat_order = ["Chemosensory", "MF", "CC", "BP", "Unknown"] if show_chemo_go else ["MF", "CC", "BP", "Unknown"]
    cols = st.columns(3)
    for col, vdf, title_stub, fname in [
        (cols[0], antp_v, c1_volcano_title(antp_v), "go_percent_antenna_vs_maxillary_palp.csv"),
        (cols[1], antleg_v, c1_volcano_title(antleg_v), "go_percent_antenna_vs_tarsi.csv"),
        (cols[2], legp_v, c1_volcano_title(legp_v), "go_percent_tarsi_vs_maxillary_palp.csv"),
    ]:
        with col:
            tbl = c1_make_percent_table_chemo(vdf, show_chemo_go)
            figp = px.bar(
                tbl, x="Dir", y="Percent", color="GO_Domain",
                color_discrete_map=go_col_map,
                category_orders={"GO_Domain": go_cat_order},
                text="TextN", barmode="stack", height=420,
                title=f"GO composition - {title_stub}",
                custom_data=["GO_Domain_full", "Tissue", "N", "Percent"],
            )
            figp.update_traces(
                textposition="inside",
                hovertemplate=(
                    "GO domain=%{customdata[0]}"
                    "<br>Tissue=%{customdata[1]}"
                    "<br>N=%{customdata[2]}"
                    "<br>Percent=%{customdata[3]:.1f}%<extra></extra>"
                ),
            )
            figp.update_yaxes(range=[0, 100])
            figp.update_layout(legend_title_text="", xaxis_title="",
                               yaxis_title="% of significant DEGs")
            st.plotly_chart(figp, use_container_width=True)
            st.download_button("Download %", data=tbl.to_csv(index=False).encode(),
                               file_name=fname, mime="text/csv",
                               key=f"c1_go_dl_{fname}")
    st.markdown("---")
    st.markdown("#### Inspect genes by contrast, direction and GO domain")
    contrast_label = st.selectbox(
        "Contrast",
        ["Antenna vs Maxillary palp", "Antenna vs Tarsi", "Tarsi vs Maxillary palp"],
        key="c1_go_filter_contrast"
    )
    if contrast_label == "Antenna vs Maxillary palp":
        vdf_sel = antp_v; contrast_code = "antp"
    elif contrast_label == "Antenna vs Tarsi":
        vdf_sel = antleg_v; contrast_code = "antleg"
    else:
        vdf_sel = legp_v; contrast_code = "legp"
    dir_opts = sorted([d for d in vdf_sel["Direction"].unique() if d != "Not sig"])
    if not dir_opts:
        st.info("No significant directions for this contrast with current thresholds.")
        return
    dir_label = st.selectbox("Direction (tissue up)", dir_opts, key="c1_go_filter_direction")
    present_domains = (
        vdf_sel[(vdf_sel["is_sig"]) & (vdf_sel["Direction"] == dir_label)]
        ["GO_Domain"].dropna().unique().tolist()
    )
    domain_order = ["MF", "CC", "BP", "Unknown"]
    present_domains = [d for d in domain_order if d in present_domains]
    if not present_domains:
        st.info("No GO domains for this direction with current thresholds.")
        return
    go_domain = st.selectbox("GO domain", present_domains, key="c1_go_filter_domain")
    detail_df = c1_build_go_group_detail(vdf_sel, dir_label, go_domain, norm_means)
    if detail_df.empty:
        st.write("No genes found for this combination.")
    else:
        st.dataframe(detail_df)
        st.download_button("Download gene table (CSV)",
                           data=detail_df.to_csv(index=False).encode(),
                           file_name=f"go_detail_{contrast_code}_{dir_label.replace(' ', '_')}_{go_domain}.csv",
                           mime="text/csv", key="c1_go_detail_dl")


def c1_render_overlap_tab(de_data, volcano_dfs, go_map, name_map, norm_means, padj_thr, lfc_thr):
    AntP = de_data["AntP"]; AntLeg = de_data["AntLeg"]; LegP = de_data["LegP"]
    antp_v = volcano_dfs["antp_v"]; antleg_v = volcano_dfs["antleg_v"]; legp_v = volcano_dfs["legp_v"]
    st.subheader("GO Name bars and overlaps per tissue")
    st.caption("Manuscript: **Figure 3B** | Overlap = genes upregulated in both relevant contrasts for that tissue.")
    drop_unknown = st.checkbox("Drop GO_Name='Unknown' in bars and Venns", value=True, key="c1_drop_unk")
    top_n = st.number_input("Top N GO Names per domain (bars)", min_value=1, value=20, key="c1_top_n")
    ant_keys = c1_up_keys(AntP, +1, padj_thr, lfc_thr) & c1_up_keys(AntLeg, +1, padj_thr, lfc_thr)
    palp_keys = c1_up_keys(AntP, -1, padj_thr, lfc_thr) & c1_up_keys(LegP, -1, padj_thr, lfc_thr)
    tarsi_keys = c1_up_keys(AntLeg, -1, padj_thr, lfc_thr) & c1_up_keys(LegP, +1, padj_thr, lfc_thr)
    domains_for_venn = ["MF", "CC", "BP", "Unknown"]
    if drop_unknown:
        domains_for_venn = [d for d in domains_for_venn if d != "Unknown"]
    for tissue, keys in [("Antenna", ant_keys), ("Maxillary palp", palp_keys), ("Tarsi", tarsi_keys)]:
        st.markdown(f"### {tissue}")
        venn_cols = st.columns(len(domains_for_venn))
        for i, dom in enumerate(domains_for_venn):
            with venn_cols[i]:
                lbl_left, lbl_right, s1, s2, ov = c1_domain_keys_for_tissue(
                    tissue, dom, antp_v, antleg_v, legp_v)
                n1, n2, n_overlap = len(s1), len(s2), len(ov)
                if n1 == 0 and n2 == 0:
                    st.caption(f"No {DOMAIN_FULL.get(dom, dom).lower()} genes with current thresholds.")
                else:
                    col_d = GO_COLS.get(dom, "#CCCCCC")
                    title = DOMAIN_FULL.get(dom, dom)
                    fvenn = c1_fig_venn_two_sets(lbl_left, lbl_right, n1, n2, n_overlap,
                                                 title=title, domain_color=col_d)
                    st.plotly_chart(fvenn, use_container_width=True,
                                    key=f"c1_venn_{tissue}_{dom}")
        bar_cols = st.columns(3)
        t_mf = c1_go_name_table_overlap(keys, go_map, "MF", drop_unknown).head(int(top_n))
        t_cc = c1_go_name_table_overlap(keys, go_map, "CC", drop_unknown).head(int(top_n))
        t_bp = c1_go_name_table_overlap(keys, go_map, "BP", drop_unknown).head(int(top_n))
        x_max = max([t["N"].max() if not t.empty else 0 for t in (t_mf, t_cc, t_bp)])
        with bar_cols[0]:
            st.plotly_chart(c1_fig_go_names(t_mf, "MF", x_max), use_container_width=True,
                            key=f"c1_gobar_mf_{tissue}")
        with bar_cols[1]:
            st.plotly_chart(c1_fig_go_names(t_cc, "CC", x_max), use_container_width=True,
                            key=f"c1_gobar_cc_{tissue}")
        with bar_cols[2]:
            st.plotly_chart(c1_fig_go_names(t_bp, "BP", x_max), use_container_width=True,
                            key=f"c1_gobar_bp_{tissue}")
        with st.expander(f"Download GO Name tables - {tissue}"):
            cx1, cx2, cx3 = st.columns(3)
            cx1.download_button("MF table (CSV)", data=t_mf.to_csv(index=False).encode(),
                                file_name=f"{tissue}_MF_overlap_top{int(top_n)}.csv",
                                mime="text/csv", key=f"c1_mf_dl_{tissue}")
            cx2.download_button("CC table (CSV)", data=t_cc.to_csv(index=False).encode(),
                                file_name=f"{tissue}_CC_overlap_top{int(top_n)}.csv",
                                mime="text/csv", key=f"c1_cc_dl_{tissue}")
            cx3.download_button("BP table (CSV)", data=t_bp.to_csv(index=False).encode(),
                                file_name=f"{tissue}_BP_overlap_top{int(top_n)}.csv",
                                mime="text/csv", key=f"c1_bp_dl_{tissue}")
    st.markdown("---")
    st.markdown("#### Inspect overlapping genes by tissue, GO domain and GO name")
    tissue_sel = st.selectbox("Focal tissue", ["Antenna", "Maxillary palp", "Tarsi"],
                               key="c1_overlap_tissue")
    key_set = {"Antenna": ant_keys, "Maxillary palp": palp_keys, "Tarsi": tarsi_keys}[tissue_sel]
    domain_sel = st.selectbox("GO domain for overlap table", ["MF", "CC", "BP", "Unknown"],
                               key="c1_overlap_domain")
    base = go_map.copy()
    base["GO_Domain"] = clean_go_domain(base["GO_Domain"])
    base = base[base["Gene"].astype(str).isin(key_set) & (base["GO_Domain"] == domain_sel)].copy()
    if base.empty:
        st.info("No overlapping genes for this combination.")
        return
    base = base.merge(name_map.rename(columns={"JoinKey": "Gene"}), on="Gene", how="left")
    go_name_opts = sorted(base["GO_Name"].dropna().unique().tolist())
    go_name_sel = st.multiselect("Filter by GO_Name (leave empty for all)", go_name_opts,
                                  default=[], key="c1_overlap_go_name")
    if go_name_sel:
        base = base[base["GO_Name"].isin(go_name_sel)]
    gene_filter = st.text_input("Optional gene filter (substring on Gene)", value="",
                                 key="c1_overlap_gene_filter")
    if gene_filter.strip() != "":
        gf = gene_filter.strip().lower()
        base = base[base["Gene"].astype(str).str.lower().str.contains(gf)]
    if norm_means is not None:
        base = base.merge(norm_means.rename(columns={"JoinKey": "Gene"}), on="Gene", how="left")
    col_order = ["Gene", "Name", "GO_Name", "GO_Domain"]
    other_cols = [c for c in base.columns if c not in col_order]
    base = base[col_order + other_cols]
    st.dataframe(base)
    st.download_button("Download overlap table (CSV)",
                       data=base.to_csv(index=False).encode(),
                       file_name=f"overlap_{tissue_sel}_{domain_sel}.csv",
                       mime="text/csv", key="c1_overlap_dl")


def c1_render_chemo_tab(chemo_tables, chemo_available, chemo_error, chemo_dir,
                         de_data, norm_means, padj_thr, lfc_thr, name_map):
    expr_thr = 10.0
    st.subheader("Chemosensory gene family classification")
    st.caption("Manuscript: **Figure 3D**")
    if not chemo_available:
        st.warning(f"Could not load chemoperception gene tables from {chemo_dir}: {chemo_error}")
        return
    has_preclassified = any(
        any(col in df.columns for col in CLASS_COLS) for df in chemo_tables.values()
    )
    fam_order = ["Or", "Gr", "Ir", "Obp", "Csp", "Ppk", "Trp"]
    tissue_keys = ["Antenna", "Palp", "Tarsi"]
    class_cols_map = {"Antenna": "Antenna_Class", "Palp": "Palp_Class", "Tarsi": "Tarsi_Class"}

    if has_preclassified:
        st.caption(
            "Using tissue classes already present in chemo_*.csv (e.g. Antenna-specific, Antenna-biased). "
            "Labels are normalized to Appendage-specific / Appendage-biased / Expressed / Not expressed."
        )
        chemo_summary = {}
        for fam in fam_order:
            fam_df = chemo_tables.get(fam, pd.DataFrame()).copy()
            if fam_df.empty:
                continue
            for col in CLASS_COLS:
                if col in fam_df.columns:
                    fam_df[col] = normalize_class_column(fam_df[col])
            fam_rep = deduplicate_family_for_pies(fam_df)
            fam_counts = {}
            for tkey in tissue_keys:
                col = class_cols_map[tkey]
                if col not in fam_rep.columns:
                    fam_counts[tkey] = {"Appendage-specific": 0, "Appendage-biased": 0,
                                        "Expressed": 0, "Not expressed": 0}
                    continue
                if "Weight" in fam_rep.columns:
                    vc = fam_rep.groupby(col)["Weight"].sum()
                else:
                    vc = fam_rep[col].value_counts(dropna=True)
                fam_counts[tkey] = {
                    "Appendage-specific": int(vc.get("Appendage-specific", 0)),
                    "Appendage-biased": int(vc.get("Appendage-biased", 0)),
                    "Expressed": int(vc.get("Expressed", 0)),
                    "Not expressed": int(vc.get("Not expressed", 0)),
                }
            chemo_summary[fam] = fam_counts
        cols = st.columns(3)
        for ti, tkey in enumerate(["Antenna", "Palp", "Tarsi"]):
            with cols[ti]:
                st.markdown(f"<h3 style='text-align:center;'>{C1_TISSUE_LABELS[tkey]}</h3>",
                            unsafe_allow_html=True)
                for fi, fam in enumerate(fam_order):
                    if tkey not in chemo_summary.get(fam, {}):
                        continue
                    show_leg = (fi == 0)
                    fig_p = c1_fig_chemo_pie(chemo_summary, fam, tkey, show_legend=show_leg)
                    if fig_p.data:
                        st.plotly_chart(fig_p, use_container_width=True,
                                        key=f"c1_pie_pre_{tkey}_{fam}")
                    else:
                        st.caption(f"No data available for {fam}")
                    st.markdown(f"<p style='text-align:center; font-size:1.0rem;'>{fam}</p>",
                                unsafe_allow_html=True)
        st.markdown("---")
        st.markdown("#### Inspect chemoperception genes")
        fam_sel = st.selectbox("Gene family", fam_order, key="c1_chemo_family_pre")
        tissue_label_sel = st.selectbox("Tissue", ["Antenna", "Maxillary palp", "Tarsi"],
                                         key="c1_chemo_tissue_pre")
        tissue_key_sel = "Palp" if tissue_label_sel == "Maxillary palp" else tissue_label_sel
        class_opts = ["Appendage-specific", "Appendage-biased", "Expressed", "Not expressed"]
        class_sel = st.multiselect("Class (leave empty for all)", class_opts, default=[],
                                    key="c1_chemo_class_pre")
        df_base = chemo_tables.get(fam_sel, pd.DataFrame()).copy()
        if df_base.empty:
            st.info("No data for this family.")
            return
        for col in CLASS_COLS:
            if col in df_base.columns:
                df_base[col] = normalize_class_column(df_base[col])
        col_class = class_cols_map[tissue_key_sel]
        if col_class not in df_base.columns:
            st.info(f"No class information for tissue {tissue_label_sel} in this family.")
            return
        s = df_base[col_class].astype(str)
        mask = pd.Series(True, index=df_base.index)
        if class_sel:
            mask = s.isin(class_sel)
        df_view = df_base.loc[mask].copy()
        front_cols = ["Gene_ID", "Family", "Antenna_Class", "Palp_Class", "Tarsi_Class"]
        front_cols = [c for c in front_cols if c in df_view.columns]
        other_cols = [c for c in df_view.columns if c not in front_cols]
        df_view = df_view[front_cols + other_cols]
        st.dataframe(df_view)
        st.download_button("Download table (CSV)", data=df_view.to_csv(index=False).encode(),
                           file_name=f"chemo_{fam_sel}_{tissue_key_sel}.csv",
                           mime="text/csv", key="c1_chemo_pre_dl")
        return

    # Branch B: compute classes from DE + expression
    if norm_means is None:
        st.warning("Normalized counts are not available - cannot compute tissue-specific/biased classes.")
        return
    st.caption(
        "Each pie shows numbers of genes that are specific, biased, "
        "expressed or not expressed per tissue and gene family. "
        f"Classification uses normalized counts threshold = {expr_thr}, "
        f"and DE thresholds padj < {padj_thr}, |log2FC| >= {lfc_thr} (strict for biased)."
    )
    AntP = de_data["AntP"]; AntLeg = de_data["AntLeg"]; LegP = de_data["LegP"]
    class_table = c1_compute_tissue_classes(norm_means, AntP, AntLeg, LegP,
                                             expr_thr=expr_thr, lfc_thr=lfc_thr, padj_thr=padj_thr)
    class_table = class_table.merge(name_map, on="JoinKey", how="left")
    with st.expander("Debug - overall class counts (all genes)"):
        for tissue, col in [("Antenna", "Antenna_Class"), ("Palp", "Palp_Class"), ("Tarsi", "Tarsi_Class")]:
            if col in class_table.columns:
                st.write(tissue, class_table[col].value_counts())
    chemo_summary = {}
    for fam in fam_order:
        fam_df = chemo_tables.get(fam, pd.DataFrame()).copy()
        if fam_df.empty or "Gene_ID" not in fam_df.columns:
            continue
        fam_df_clean = fam_df.copy()
        fam_df_clean["Gene_ID_clean"] = fam_df_clean["Gene_ID"].astype(str).str.strip().str.upper()
        class_table_clean = class_table.copy()
        if "Name" in class_table_clean.columns:
            class_table_clean["Name_clean"] = class_table_clean["Name"].astype(str).str.strip().str.upper()
        class_table_clean["JoinKey_clean"] = class_table_clean["JoinKey"].astype(str).str.strip().str.upper()
        merged = None; best_match_count = -1
        if "Name_clean" in class_table_clean.columns:
            merged_name = fam_df_clean.merge(class_table_clean, left_on="Gene_ID_clean",
                                              right_on="Name_clean", how="left", suffixes=("", "_name"))
            match_count = sum(merged_name[col].notna().sum()
                               for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                               if col in merged_name.columns)
            if match_count > best_match_count:
                merged = merged_name; best_match_count = match_count
        merged_join = fam_df_clean.merge(class_table_clean, left_on="Gene_ID_clean",
                                          right_on="JoinKey_clean", how="left", suffixes=("", "_join"))
        match_count = sum(merged_join[col].notna().sum()
                           for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                           if col in merged_join.columns)
        if match_count > best_match_count:
            merged = merged_join; best_match_count = match_count
        if merged is None or best_match_count == 0:
            merged_raw = fam_df.merge(class_table, left_on="Gene_ID", right_on="Name", how="left")
            match_count_raw = sum(merged_raw[col].notna().sum()
                                   for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                                   if col in merged_raw.columns)
            if match_count_raw > best_match_count:
                merged = merged_raw; best_match_count = match_count_raw
            if "JoinKey" in class_table.columns:
                merged_raw2 = fam_df.merge(class_table, left_on="Gene_ID", right_on="JoinKey", how="left")
                match_count_raw2 = sum(merged_raw2[col].notna().sum()
                                        for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                                        if col in merged_raw2.columns)
                if match_count_raw2 > best_match_count:
                    merged = merged_raw2; best_match_count = match_count_raw2
        merged_rep = deduplicate_family_for_pies(merged)
        fam_counts = {}
        for tkey in tissue_keys:
            col = class_cols_map[tkey]
            if col not in merged_rep.columns:
                fam_counts[tkey] = {"Appendage-specific": 0, "Appendage-biased": 0,
                                    "Expressed": 0, "Not expressed": 0}
                continue
            if "Weight" in merged_rep.columns:
                vc = merged_rep.groupby(col)["Weight"].sum()
            else:
                vc = merged_rep[col].value_counts(dropna=True)
            fam_counts[tkey] = {
                "Appendage-specific": int(vc.get("Appendage-specific", 0)),
                "Appendage-biased": int(vc.get("Appendage-biased", 0)),
                "Expressed": int(vc.get("Expressed", 0)),
                "Not expressed": int(vc.get("Not expressed", 0)),
            }
        chemo_summary[fam] = fam_counts
    cols = st.columns(3)
    for ti, tkey in enumerate(["Antenna", "Palp", "Tarsi"]):
        with cols[ti]:
            st.markdown(f"<h3 style='text-align:center;'>{C1_TISSUE_LABELS[tkey]}</h3>",
                        unsafe_allow_html=True)
            for fi, fam in enumerate(fam_order):
                if tkey not in chemo_summary.get(fam, {}):
                    continue
                show_leg = (fi == 0)
                fig_p = c1_fig_chemo_pie(chemo_summary, fam, tkey, show_legend=show_leg)
                if fig_p.data:
                    st.plotly_chart(fig_p, use_container_width=True,
                                    key=f"c1_pie_fb_{tkey}_{fam}")
                else:
                    st.caption(f"No data available for {fam}")
                st.markdown(f"<p style='text-align:center; font-size:1.0rem;'>{fam}</p>",
                            unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("#### Inspect chemoperception genes")
    fam_sel = st.selectbox("Gene family", fam_order, key="c1_chemo_family_fb")
    tissue_label_sel = st.selectbox("Tissue", ["Antenna", "Maxillary palp", "Tarsi"],
                                     key="c1_chemo_tissue_fb")
    tissue_key_sel = "Palp" if tissue_label_sel == "Maxillary palp" else tissue_label_sel
    class_opts = ["Appendage-specific", "Appendage-biased", "Expressed", "Not expressed"]
    class_sel = st.multiselect("Class (leave empty for all)", class_opts, default=[],
                                key="c1_chemo_class_fb")
    df_base = chemo_tables.get(fam_sel, pd.DataFrame()).copy()
    if df_base.empty or "Gene_ID" not in df_base.columns:
        st.info("No data for this family.")
        return
    fam_df_clean = df_base.copy()
    fam_df_clean["Gene_ID_clean"] = fam_df_clean["Gene_ID"].astype(str).str.strip().str.upper()
    class_table_clean = class_table.copy()
    if "Name" in class_table_clean.columns:
        class_table_clean["Name_clean"] = class_table_clean["Name"].astype(str).str.strip().str.upper()
    class_table_clean["JoinKey_clean"] = class_table_clean["JoinKey"].astype(str).str.strip().str.upper()
    merged = None; best_match_count = -1
    if "Name_clean" in class_table_clean.columns:
        merged_name = fam_df_clean.merge(class_table_clean, left_on="Gene_ID_clean",
                                          right_on="Name_clean", how="left", suffixes=("", "_name"))
        match_count = sum(merged_name[col].notna().sum()
                           for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                           if col in merged_name.columns)
        if match_count > best_match_count:
            merged = merged_name; best_match_count = match_count
    merged_join = fam_df_clean.merge(class_table_clean, left_on="Gene_ID_clean",
                                      right_on="JoinKey_clean", how="left", suffixes=("", "_join"))
    match_count = sum(merged_join[col].notna().sum()
                       for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                       if col in merged_join.columns)
    if match_count > best_match_count:
        merged = merged_join; best_match_count = match_count
    if merged is None or best_match_count == 0:
        merged_raw = df_base.merge(class_table, left_on="Gene_ID", right_on="Name", how="left")
        match_count_raw = sum(merged_raw[col].notna().sum()
                               for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                               if col in merged_raw.columns)
        if match_count_raw > best_match_count:
            merged = merged_raw; best_match_count = match_count_raw
        if "JoinKey" in class_table.columns:
            merged_raw2 = df_base.merge(class_table, left_on="Gene_ID", right_on="JoinKey", how="left")
            match_count_raw2 = sum(merged_raw2[col].notna().sum()
                                    for col in ["Antenna_Class", "Palp_Class", "Tarsi_Class"]
                                    if col in merged_raw2.columns)
            if match_count_raw2 > best_match_count:
                merged = merged_raw2
    col_class = class_cols_map[tissue_key_sel]
    if col_class not in merged.columns:
        st.info(f"No class information for tissue {tissue_label_sel} in this family.")
        return
    s = merged[col_class].astype(str)
    mask = pd.Series(True, index=merged.index)
    if class_sel:
        mask = s.isin(class_sel)
    df_view = merged.loc[mask].copy()
    front_cols = ["Gene_ID", "Family", "Antenna_Class", "Palp_Class", "Tarsi_Class"]
    front_cols = [c for c in front_cols if c in df_view.columns]
    other_cols = [c for c in df_view.columns if c not in front_cols]
    df_view = df_view[front_cols + other_cols]
    st.dataframe(df_view)
    st.download_button("Download table (CSV)", data=df_view.to_csv(index=False).encode(),
                       file_name=f"chemo_{fam_sel}_{tissue_key_sel}.csv",
                       mime="text/csv", key="c1_chemo_fb_dl")


def render_c1_tab():
    with st.sidebar:
        st.header("Appendage Comparison")
        padj_thr = st.number_input("padj threshold", value=0.001, min_value=1e-12,
                                    max_value=1.0, step=0.0005, format="%.6f",
                                    key="c1_padj_thr")
        lfc_thr = st.number_input("|log2FC| threshold", value=1.0, min_value=0.0,
                                   max_value=10.0, step=0.1, key="c1_lfc_thr")
        strong_lfc = st.number_input("Strong |log2FC| (alpha)", value=2.0, min_value=0.0,
                                      max_value=10.0, step=0.5, key="c1_strong_lfc")
        st.markdown("---")
        show_chemo = st.checkbox("Overlay chemosensory points (dark red)", value=True,
                                  key="c1_show_chemo")
        show_chemo_go = st.checkbox("Show chemosensory % separately in GO bars", value=False,
                                     key="c1_show_chemo_go")
        highlight_query = st.text_input("Search gene (Gene / Name / JoinKey)", value="",
                                         placeholder="e.g. XM_038046469.1 or OR120",
                                         key="c1_highlight_query")

    data = c1_load_all_data(C1_DEFAULT_PATHS, padj_thr, lfc_thr, strong_lfc)
    volcano_dfs = data["VOLCANO"]
    de_data = data["de"]
    go_map = data["GO_MAP"]
    norm_means = data["NORM_MEANS"]
    name_map = data["NAME_MAP"]
    chemo_tables = data["CHEMO_TABLES"]
    chemo_available = data["CHEMO_AVAILABLE"]
    chemo_error = data["CHEMO_ERROR"]
    chemo_dir = data["paths"]["chemo_dir"]

    vol_tab, go_tab, overlap_tab, chemo_tab = st.tabs(
        ["Volcano (Fig. 3A)", "Composition (Fig. 3B)", "Consensus sets (Fig. 3C)",
         "Chemosensory classes (Fig. 3D)"]
    )
    with vol_tab:
        c1_render_volcano_tab(volcano_dfs, padj_thr, lfc_thr, show_chemo, highlight_query)
    with go_tab:
        c1_render_composition_tab(volcano_dfs, show_chemo_go, norm_means)
    with overlap_tab:
        c1_render_consensus_tab(
            de_data, padj_thr, lfc_thr,
            lambda: c1_render_overlap_tab(de_data, volcano_dfs, go_map, name_map,
                                          norm_means, padj_thr, lfc_thr))
    with chemo_tab:
        c1_render_chemo_tab(chemo_tables, chemo_available, chemo_error, chemo_dir,
                             de_data, norm_means, padj_thr, lfc_thr, name_map)


# =============================================================================
# === APP_S1 FUNCTIONS ===
# =============================================================================

def s1_get_cond_colors(cond: str):
    cond = str(cond)
    if "MF" in cond:
        return COL_MF
    if "VF" in cond and "Vm" not in cond:
        return COL_VF
    if "Vm" in cond:
        return COL_Vm
    return dict(light="#CCCCCC", dark="#666666")


def s1_get_full_label_color(full_name: str) -> str:
    if full_name == COND_FULL.get("Vm"):
        return COL_Vm["dark"]
    if full_name == COND_FULL.get("VF"):
        return COL_VF["dark"]
    if full_name == COND_FULL.get("MF"):
        return COL_MF["dark"]
    return "#000000"


def s1_extract_chemo_tag(text: str):
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return np.nan
    m = re.search(r"(?i)\b(OR|IR|GR|OBP|CSP|PPK)[0-9A-Za-z._-]*\b", str(text))
    return m.group(0).upper() if m else np.nan


def s1_build_join_key(row, name_col="Name", gene_col="Gene"):
    gene_val = str(row[gene_col]).strip() if gene_col in row and pd.notna(row[gene_col]) else ""
    name_val = str(row[name_col]).strip() if name_col in row and pd.notna(row[name_col]) else ""
    if gene_val:
        return gene_val
    if name_val:
        return name_val
    return np.nan


def s1_first_non_nan_string(*vals):
    for v in vals:
        if isinstance(v, str) and v.strip() != "":
            return v
    return np.nan


@st.cache_data(show_spinner=False)
def s1_load_de_table(path):
    df = pd.read_csv(path)
    if "Name" not in df.columns:
        df["Name"] = np.nan
    if "Gene" not in df.columns:
        df["Gene"] = np.nan
    if "padj" not in df.columns:
        st.error(f"File {os.path.basename(path)} is missing column 'padj'")
    if "log2FoldChange" not in df.columns:
        st.error(f"File {os.path.basename(path)} is missing column 'log2FoldChange'")
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = -pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["JoinKey"] = df.apply(s1_build_join_key, axis=1)
    df["ChemoName"] = df.apply(
        lambda r: s1_first_non_nan_string(
            s1_extract_chemo_tag(r.get("Name")),
            s1_extract_chemo_tag(r.get("Gene")),
            s1_extract_chemo_tag(r.get("JoinKey")),
        ), axis=1,
    )
    df["is_chemo"] = df["ChemoName"].notna()
    padj_safe = df["padj"].clip(lower=np.finfo(float).tiny)
    df["neglog10_padj"] = -np.log10(padj_safe)
    return df


def s1_classify_de(df, cond1, cond2, padj_thr, lfc_thr, strong_lfc):
    df = df.copy()
    df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    df["is_sig"] = (df["padj"] < padj_thr) & (df["log2FoldChange"].abs() >= lfc_thr)
    df["Side"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df.loc[df["is_sig"] & (df["log2FoldChange"] >= lfc_thr), "Side"] = cond1
    df.loc[df["is_sig"] & (df["log2FoldChange"] <= -lfc_thr), "Side"] = cond2
    df["Strength"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df.loc[df["is_sig"] & (df["log2FoldChange"].abs() >= strong_lfc), "Strength"] = "Strong"
    df.loc[df["is_sig"] & df["Strength"].isna(), "Strength"] = "Moderate"
    return df


@st.cache_data(show_spinner=False)
def s1_load_go_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        st.warning(f"GO map file not found at {path}")
        return pd.DataFrame(columns=["Gene", "GO_Name", "GO_Domain"])
    go_raw = pd.read_csv(path)
    nm = list(go_raw.columns)
    if not nm:
        return pd.DataFrame(columns=["Gene", "GO_Name", "GO_Domain"])
    first_col = nm[0]

    def pick(colnames, opts):
        low = {c.lower(): c for c in colnames}
        for k in opts:
            if k in low:
                return low[k]
        return None

    go_name_col = pick(nm, ["go_name", "go term name", "go_term_name", "go term",
                             "go_description", "go_desc", "go label", "goname", "go"])
    go_domain_col = pick(nm, ["go_domain", "domain", "aspect", "goaspect",
                               "go_domain_name", "go_aspect"])
    out = pd.DataFrame({
        "Gene": go_raw[first_col].astype(str).str.strip(),
        "GO_Name": go_raw[go_name_col].astype(str).str.strip() if go_name_col else "Unknown",
        "GO_Domain": go_raw[go_domain_col] if go_domain_col else "Unknown",
    })
    out["GO_Domain"] = clean_go_domain(out["GO_Domain"])
    out.loc[out["GO_Name"].isna() | (out["GO_Name"] == ""), "GO_Name"] = "Unknown"
    out = out.drop_duplicates(subset=["Gene"], keep="first").reset_index(drop=True)
    return out


def s1_annotate_with_go(df, go_map):
    if go_map is None or go_map.empty:
        out = df.copy()
        out["GO_Name"] = "Unknown"
        out["GO_Domain"] = "Unknown"
        return out
    out = df.copy()
    if "JoinKey" not in out.columns:
        out["JoinKey"] = out.apply(s1_build_join_key, axis=1)
    out = out.merge(go_map, how="left", left_on="JoinKey", right_on="Gene", suffixes=("", "_go"))
    out["GO_Domain"] = clean_go_domain(out["GO_Domain"])
    out.loc[out["GO_Name"].isna() | (out["GO_Name"] == ""), "GO_Name"] = "Unknown"
    if "Gene_go" in out.columns:
        out["Gene"] = out.apply(
            lambda r: s1_first_non_nan_string(r.get("Gene"), r.get("Gene_go")), axis=1)
        out = out.drop(columns=["Gene_go"])
    return out


@st.cache_data(show_spinner=False)
def s1_load_norm_table(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        st.warning(f"Normalized counts file not found at {path}")
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "Name" not in df.columns:
        df["Name"] = np.nan
    if "Gene" not in df.columns:
        df["Gene"] = np.nan
    df["JoinKey"] = df.apply(s1_build_join_key, axis=1)
    return df


def s1_compute_tissue_state_means(norm_df: pd.DataFrame) -> pd.DataFrame:
    if norm_df.empty:
        return pd.DataFrame(columns=["JoinKey"])
    norm_df = norm_df.copy()
    num_cols = norm_df.select_dtypes(include=[np.number]).columns
    tissue_tags = {"Antenna": ["ant", "a_"], "Palp": ["palp", "pal", "p_"],
                   "Tarsi": ["tars", "tar", "leg", "l_"]}
    state_tags = {"MF": ["mf"], "VF": ["vf"], "Vm": ["vm"]}
    for tissue, ttags in tissue_tags.items():
        t_tags = [t.lower() for t in ttags]
        for state, stags in state_tags.items():
            s_tags = [s.lower() for s in stags]
            cols = [c for c in num_cols
                    if any(t in c.lower() for t in t_tags) and any(s in c.lower() for s in s_tags)]
            if cols:
                new_col = f"mean_{tissue}_{state}"
                norm_df[new_col] = norm_df[cols].mean(axis=1)
    mean_cols = [c for c in norm_df.columns if c.startswith("mean_")]
    return norm_df[["JoinKey"] + mean_cols] if mean_cols else norm_df[["JoinKey"]]


def s1_de_counts_caption(df, cond1, cond2):
    n1 = int(((df["is_sig"]) & (df["Side"] == cond1)).sum())
    n2 = int(((df["is_sig"]) & (df["Side"] == cond2)).sum())
    label1 = COND_FULL.get(cond1, cond1)
    label2 = COND_FULL.get(cond2, cond2)
    return f"{label2}: {n2} | {label1}: {n1}"


def s1_find_hits(vdf: pd.DataFrame, query: str) -> pd.DataFrame:
    if not query or str(query).strip() == "":
        return vdf.iloc[0:0]
    q = str(query).strip().lower()
    cols = [c for c in ["Gene", "Name", "JoinKey", "ChemoName"] if c in vdf.columns]
    if not cols:
        return vdf.iloc[0:0]
    mask_total = None
    for c in cols:
        s = vdf[c].astype(str).str.strip().str.lower()
        m = s == q
        mask_total = m if mask_total is None else (mask_total | m)
    if mask_total is None:
        return vdf.iloc[0:0]
    hits = vdf[mask_total]
    if hits.empty:
        return vdf.iloc[0:0]
    return hits.iloc[[0]]


def s1_build_volcano_figure(df, cond1, cond2, padj_thr, lfc_thr,
                             search_term="", overlay_chemo=True, color_by_go=True):
    disp_left = COND_FULL.get(cond2, cond2)
    disp_right = COND_FULL.get(cond1, cond1)
    vdf = df.copy()
    if color_by_go and "ColorKey_GO" in vdf.columns:
        vdf["ColorKey_plot"] = vdf["ColorKey_GO"]
        palette = S1_PAL_GO
    else:
        def classify_row(r):
            if not r.get("is_sig", False) or pd.isna(r.get("Side")):
                return "Not significant"
            side = r["Side"]
            strength = r.get("Strength", "Moderate")
            if strength not in ["Moderate", "Strong"]:
                strength = "Moderate"
            return f"{side} {strength}"
        vdf["ColorKey_plot"] = vdf.apply(classify_row, axis=1)
        col1 = s1_get_cond_colors(cond1)
        col2 = s1_get_cond_colors(cond2)
        palette = {
            f"{cond1} Moderate": col1["light"], f"{cond1} Strong": col1["dark"],
            f"{cond2} Moderate": col2["light"], f"{cond2} Strong": col2["dark"],
            "Not significant": "#D9D9D9",
        }
    col_left = s1_get_full_label_color(disp_left)
    col_right = s1_get_full_label_color(disp_right)
    x_label_html = (
        f"Log2 fold change (<span style='color:{col_left}'>{disp_left}</span> vs "
        f"<span style='color:{col_right}'>{disp_right}</span>)"
    )
    hover_data = {
        "JoinKey": True, "Name": "Name" in vdf.columns,
        "ChemoName": "ChemoName" in vdf.columns,
        "GO_Name": "GO_Name" in vdf.columns, "GO_Domain": "GO_Domain" in vdf.columns,
        "padj": True, "log2FoldChange": True,
    }
    fig = px.scatter(
        vdf, x="log2FoldChange", y="neglog10_padj",
        color="ColorKey_plot", color_discrete_map=palette, hover_data=hover_data,
        labels={"log2FoldChange": x_label_html, "neglog10_padj": "-log10(<i>padj</i>)"},
        render_mode="webgl",
    )
    fig.update_traces(marker=dict(size=5, line=dict(width=0)))
    fig.add_hline(y=-np.log10(padj_thr), line_dash="dash", line_width=0.8)
    fig.add_vline(x=lfc_thr, line_dash="dash", line_width=0.8)
    fig.add_vline(x=-lfc_thr, line_dash="dash", line_width=0.8)
    if overlay_chemo:
        chemo = vdf[(vdf["is_sig"]) & (vdf["is_chemo"]) & vdf["ChemoName"].notna()]
        if not chemo.empty:
            fig.add_trace(go.Scattergl(
                x=chemo["log2FoldChange"], y=chemo["neglog10_padj"],
                mode="markers", marker=dict(color="#8B0000", size=7),
                name="Chemosensory", hovertext=chemo["ChemoName"], hoverinfo="text",
            ))
    if search_term:
        hits = s1_find_hits(vdf, search_term)
        if not hits.empty:
            hit = hits.iloc[[0]]
            fig.add_trace(go.Scattergl(
                x=hit["log2FoldChange"], y=hit["neglog10_padj"], mode="markers",
                marker=dict(size=7, symbol="circle", color="red",
                            line=dict(width=1, color="white")),
                name="Search hit", hovertext=hit["JoinKey"], hoverinfo="text",
            ))
    fig.update_layout(template="simple_white", legend_title=None)
    fig.update_yaxes(range=[0, max(5, vdf["neglog10_padj"].max() * 1.05)])
    fig.update_xaxes(range=[-10, 10])
    return fig


def s1_build_bar_figure(df, cond1, cond2):
    total_genes = len(df)
    df_counts = df[df["is_sig"] & df["Side"].notna() & df["Strength"].notna()].copy()
    if df_counts.empty:
        return None
    counts = (df_counts.groupby(["Side", "Strength"], as_index=False)
               .size().rename(columns={"size": "N"}))
    all_idx = pd.MultiIndex.from_product([[cond1, cond2], ["Moderate", "Strong"]],
                                          names=["Side", "Strength"])
    counts = counts.set_index(["Side", "Strength"]).reindex(all_idx, fill_value=0).reset_index()
    if total_genes > 0:
        counts["Percent"] = counts["N"].apply(lambda n: 100 * n / total_genes)
    else:
        counts["Percent"] = 0.0
    counts["FillKey"] = counts["Side"] + " " + counts["Strength"]
    counts["Label"] = counts.apply(lambda r: f"{int(r['N'])} ({r['Percent']:.1f}%)", axis=1)
    COL1 = s1_get_cond_colors(cond1)
    COL2 = s1_get_cond_colors(cond2)
    fill_cols = {
        f"{cond1} Moderate": COL1["light"], f"{cond1} Strong": COL1["dark"],
        f"{cond2} Moderate": COL2["light"], f"{cond2} Strong": COL2["dark"],
    }
    x_order = [cond2, cond1]
    fig = px.bar(counts, x="Side", y="N", color="FillKey", color_discrete_map=fill_cols,
                 text="Label", category_orders={"Side": x_order})
    fig.update_layout(barmode="stack", template="simple_white", legend_title=None,
                      yaxis_title="Number of significant genes", xaxis_title="")
    fig.update_traces(textposition="inside")
    y_max = max(1, counts["N"].max())
    fig.update_yaxes(range=[0, y_max * 1.2])
    return fig


@st.cache_data(show_spinner=False)
def s1_load_all_contrasts(base_dir, padj_thr, lfc_thr, strong_lfc):
    go_path = os.path.join(base_dir, S1_DEFAULT_PATHS["FILE_GO"])
    norm_path = os.path.join(base_dir, S1_DEFAULT_PATHS["FILE_NORM"])
    go_map = s1_load_go_table(go_path)
    norm_raw = s1_load_norm_table(norm_path)
    norm_means = s1_compute_tissue_state_means(norm_raw) if not norm_raw.empty else pd.DataFrame(columns=["JoinKey"])
    by_tissue = {}
    for label, fname in S1_DEFAULT_PATHS["DE_FILES"].items():
        tissue_label, contrast_str = label.split(" - ", 1)
        tissue = tissue_label
        path = os.path.join(base_dir, fname)
        if not os.path.exists(path):
            st.warning(f"DE file not found: {path}")
            continue
        df_raw = s1_load_de_table(path)
        parts = contrast_str.split("vs")
        if len(parts) != 2:
            continue
        cond1 = parts[0].strip()
        cond2 = parts[1].strip()
        df_cls = s1_classify_de(df_raw, cond1, cond2, padj_thr, lfc_thr, strong_lfc)
        df_go = s1_annotate_with_go(df_cls, go_map)
        df_go["ColorKey_GO"] = np.where(df_go["is_sig"], df_go["GO_Domain"], "Not Significant")
        if not norm_means.empty:
            col1_raw = f"mean_{tissue}_{cond1}"
            col2_raw = f"mean_{tissue}_{cond2}"
            cols_to_merge = [c for c in [col1_raw, col2_raw] if c in norm_means.columns]
            if cols_to_merge:
                tmp = norm_means[["JoinKey"] + cols_to_merge]
                df_go = df_go.merge(tmp, on="JoinKey", how="left")
                rename_map = {}
                if col1_raw in cols_to_merge:
                    rename_map[col1_raw] = f"Mean {COND_FULL.get(cond1, cond1)} ({TISSUE_DISPLAY[tissue]})"
                if col2_raw in cols_to_merge:
                    rename_map[col2_raw] = f"Mean {COND_FULL.get(cond2, cond2)} ({TISSUE_DISPLAY[tissue]})"
                if rename_map:
                    df_go = df_go.rename(columns=rename_map)
        by_tissue.setdefault(tissue, []).append({"cond1": cond1, "cond2": cond2, "df": df_go})
    return by_tissue, norm_means


@st.cache_data(show_spinner=False)
def s1_build_master_annotation(by_tissue, norm_means):
    frames = []
    for tissue, contrast_list in by_tissue.items():
        for item in contrast_list:
            df = item["df"]
            cols = [c for c in ["JoinKey", "Gene", "Name", "GO_Domain", "GO_Name"] if c in df.columns]
            frames.append(df[cols])
    if not frames:
        base = pd.DataFrame(columns=["JoinKey", "Gene", "Name", "GO_Domain", "GO_Name"])
    else:
        base = pd.concat(frames, ignore_index=True)
        base = base.drop_duplicates(subset=["JoinKey"], keep="first")
    if not norm_means.empty:
        base = base.merge(norm_means, on="JoinKey", how="left")
    return base


def s1_compute_venn_regions(setA, setB, setC):
    A, B, C = setA, setB, setC
    return {
        "A only": A - B - C, "B only": B - A - C, "C only": C - A - B,
        "A & B": (A & B) - C, "A & C": (A & C) - B, "B & C": (B & C) - A,
        "A & B & C": A & B & C,
    }


def s1_get_up_sets_venn(by_tissue, state_code, partner_code=None):
    tissue_sets = {}
    for tissue, contrast_list in by_tissue.items():
        genes = set()
        for item in contrast_list:
            cond1 = item["cond1"]; cond2 = item["cond2"]; df = item["df"].copy()
            if partner_code is None:
                if state_code not in {cond1, cond2}:
                    continue
            else:
                if {state_code, partner_code} != {cond1, cond2}:
                    continue
            df["padj"] = pd.to_numeric(df["padj"], errors="coerce")
            df["log2FoldChange"] = pd.to_numeric(df["log2FoldChange"], errors="coerce")
            base_mask = df["padj"] < VENN_PADJ_THR
            if cond1 == state_code:
                mask = base_mask & (df["log2FoldChange"] >= VENN_LFC_THR)
            elif cond2 == state_code:
                mask = base_mask & (df["log2FoldChange"] <= -VENN_LFC_THR)
            else:
                continue
            m2 = mask & df["JoinKey"].notna()
            g = df.loc[m2, "JoinKey"].astype(str)
            genes |= set(g)
        tissue_sets[tissue] = genes
    return tissue_sets


def s1_build_venn_figure(state_full_label, state_color, labels, regions_counts):
    def hex_to_rgba_local(hex_color, alpha=0.20):
        if not isinstance(hex_color, str):
            return hex_color
        hc = hex_color.lstrip("#")
        if len(hc) != 6:
            return hex_color
        r = int(hc[0:2], 16)
        g = int(hc[2:4], 16)
        b = int(hc[4:6], 16)
        return f"rgba({r},{g},{b},{alpha})"

    fig = go.Figure()
    r = 1.5
    circles = [
        {"x0": 0.0, "y0": 1.8, "name": "A"},
        {"x0": -r, "y0": 0.0, "name": "B"},
        {"x0": r, "y0": 0.0, "name": "C"},
    ]
    fill_rgba = hex_to_rgba_local(state_color, alpha=0.20)
    for c in circles:
        fig.add_shape(type="circle", x0=c["x0"]-r, y0=c["y0"]-r, x1=c["x0"]+r, y1=c["y0"]+r,
                      line=dict(width=0), fillcolor=fill_rgba, layer="below")
    text_positions = {
        "A only": (0.0, 3.0), "B only": (-2.0, -0.2), "C only": (2.0, -0.2),
        "A & B": (-1.0, 1.2), "A & C": (1.0, 1.2), "B & C": (0.0, -0.4), "A & B & C": (0.0, 0.9),
    }
    for key, pos in text_positions.items():
        count = regions_counts.get(key, 0)
        fig.add_trace(go.Scatter(x=[pos[0]], y=[pos[1]], text=[str(count)], mode="text",
                                  textfont=dict(size=16), showlegend=False, hoverinfo="skip"))
    labA, labB, labC = labels
    fig.add_annotation(x=0.0, y=3.6, text=labA, showarrow=False,
                       font=dict(size=14, color=state_color))
    fig.add_annotation(x=-2.5, y=-1.8, text=labB, showarrow=False,
                       font=dict(size=14, color=state_color))
    fig.add_annotation(x=2.5, y=-1.8, text=labC, showarrow=False,
                       font=dict(size=14, color=state_color))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(template="simple_white", showlegend=False,
                      margin=dict(l=10, r=10, t=10, b=10), width=380, height=320)
    return fig


def s1_make_region_table(state_code, gene_set, master_annot):
    if not gene_set:
        return pd.DataFrame(columns=["Gene", "Name", "GO_Domain", "GO_Name"])
    df = master_annot.copy()
    df["JoinKey_str"] = df["JoinKey"].astype(str)
    df = df[df["JoinKey_str"].isin(set(gene_set))].drop(columns=["JoinKey_str"])
    state_full = COND_FULL.get(state_code, state_code)
    tissue_order = ["Antenna", "Palp", "Tarsi"]
    for tissue in tissue_order:
        raw_col = f"mean_{tissue}_{state_code}"
        friendly = f"Mean {state_full} ({TISSUE_DISPLAY[tissue]})"
        if raw_col in df.columns:
            df[friendly] = df[raw_col]
    cols = []
    for c in ["Gene", "Name"]:
        if c in df.columns:
            cols.append(c)
    for tissue in tissue_order:
        friendly = f"Mean {state_full} ({TISSUE_DISPLAY[tissue]})"
        if friendly in df.columns:
            cols.append(friendly)
    for c in ["GO_Domain", "GO_Name"]:
        if c in df.columns:
            cols.append(c)
    if cols:
        df = df[cols]
    return df


def s1_sort_key(item):
    c1, c2 = item["cond1"], item["cond2"]
    if c1 == "VF" and c2 == "Vm":
        return 0
    if c1 == "MF" and c2 == "VF":
        return 1
    return 2


def s1_fig_venn_sex_mating(n_sex_only, n_mat_only, n_overlap, tissue_label):
    """2-circle Venn: sex-biased (left, pink) vs mating-regulated (right, purple) for one tissue."""
    fill_sex = "rgba(180, 180, 180, 0.30)"
    fill_mat = "rgba(180, 180, 180, 0.30)"
    r = 1.0
    c1x, c2x = -0.5, 0.5
    fig = go.Figure()
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c1x-r, x1=c1x+r, y0=-r, y1=r,
                  line=dict(color="#888888", width=1.5), fillcolor=fill_sex)
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c2x-r, x1=c2x+r, y0=-r, y1=r,
                  line=dict(color="#888888", width=1.5), fillcolor=fill_mat)
    fig.add_annotation(x=c1x-0.4, y=0.05, text=str(n_sex_only), showarrow=False,
                       font=dict(size=15, color="black"))
    fig.add_annotation(x=c2x+0.4, y=0.05, text=str(n_mat_only), showarrow=False,
                       font=dict(size=15, color="black"))
    fig.add_annotation(x=0, y=0.05, text=str(n_overlap), showarrow=False,
                       font=dict(size=15, color="white"))
    fig.add_annotation(x=c1x, y=-1.3, text=f"Sex-biased<br>(total {n_sex_only + n_overlap})",
                       showarrow=False, font=dict(size=11, color="#444444"))
    fig.add_annotation(x=c2x, y=-1.3, text=f"Mating-regulated<br>(total {n_mat_only + n_overlap})",
                       showarrow=False, font=dict(size=11, color="#444444"))
    fig.update_xaxes(visible=False, range=[-2.1, 2.1])
    fig.update_yaxes(visible=False, range=[-1.8, 1.4])
    fig.update_layout(
        title=dict(text=f"<b>{tissue_label}</b>", x=0.5, font=dict(size=13)),
        showlegend=False, margin=dict(l=5, r=5, t=35, b=5),
        height=280, width=310, plot_bgcolor="white", paper_bgcolor="white",
    )
    return fig


def s1_assign_domain_chemo_priority(df):
    """Return per-row GO domain with chemosensory genes overriding to 'Chemosensory'."""
    if "ChemoName" in df.columns:
        return np.where(df["ChemoName"].notna(), "Chemosensory", df["GO_Domain"].fillna("Unknown"))
    return df["GO_Domain"].fillna("Unknown").values


def s1_compute_domain_venns(tissue_cross, tissue):
    """Per-GO-domain Venn counts (sex-only | both | mating-only) with chemosensory priority."""
    if tissue not in tissue_cross:
        return {}
    td = tissue_cross[tissue]
    sex_df = td["sex_df"]
    mat_df = td["mat_df"]
    all_sex = td["sex_vf"] | td["sex_vm"]
    all_mat = td["mat_mf"] | td["mat_vf"]
    ov = td["overlap"]
    sex_only = all_sex - ov
    mat_only = all_mat - ov

    def domain_map(df, gene_set):
        gs = {str(g) for g in gene_set}
        sub = df[df["JoinKey"].astype(str).isin(gs)].drop_duplicates("JoinKey").copy()
        if sub.empty:
            return {}
        sub["_dom"] = s1_assign_domain_chemo_priority(sub)
        return dict(zip(sub["JoinKey"].astype(str), sub["_dom"]))

    sex_dom = domain_map(sex_df, sex_only | ov)
    mat_dom = domain_map(mat_df, mat_only | ov)

    results = {}
    for dom in ["Chemosensory", "MF", "CC", "BP", "Unknown"]:
        n_sex = sum(1 for g in sex_only if sex_dom.get(str(g)) == dom)
        n_mat = sum(1 for g in mat_only if mat_dom.get(str(g)) == dom)
        n_ov  = sum(1 for g in ov if (sex_dom.get(str(g)) or mat_dom.get(str(g))) == dom)
        results[dom] = (n_sex, n_mat, n_ov)
    return results


def s1_fig_domain_venn(n_sex_only, n_mat_only, n_overlap, domain, dom_color):
    """2-circle Venn for a single GO domain, using domain color."""
    fill = hex_to_rgba(dom_color, 0.30)
    r = 1.0; c1x, c2x = -0.5, 0.5
    fig = go.Figure()
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c1x-r, x1=c1x+r, y0=-r, y1=r,
                  line=dict(color="rgba(0,0,0,0)"), fillcolor=fill)
    fig.add_shape(type="circle", xref="x", yref="y",
                  x0=c2x-r, x1=c2x+r, y0=-r, y1=r,
                  line=dict(color="rgba(0,0,0,0)"), fillcolor=fill)
    for x, y, txt in [(c1x-0.35, 0, str(n_sex_only)), (0, 0, str(n_overlap)), (c2x+0.35, 0, str(n_mat_only))]:
        fig.add_annotation(x=x, y=y, text=txt, showarrow=False,
                           font=dict(size=13, color="white" if x == 0 else "black"))
    fig.update_xaxes(visible=False, range=[-2, 2])
    fig.update_yaxes(visible=False, range=[-1.5, 1.5])
    dom_full = {**DOMAIN_FULL, "Chemosensory": "Chemosensory"}
    fig.update_layout(
        title=dict(text=dom_full.get(domain, domain), x=0.5,
                   font=dict(size=11, color=dom_color)),
        showlegend=False, margin=dict(l=2, r=2, t=30, b=2),
        height=210, width=220, plot_bgcolor="white", paper_bgcolor="white",
    )
    return fig


# =============================================================================
# === FIGURE 4D: REORGANISATION OF SEX BIAS BY MATING ===
# Sex bias measured in virgins (x, VF vs Vm) against sex bias measured in mated
# females (y, MF vs Vm).  Both axes are female-minus-male against the SAME
# reference group, virgin males, so a change of sign between them is a genuine
# reversal rather than an artefact of a shared denominator.
#
# Two layers are deliberately kept apart, as in the manuscript:
#   POSITION is direction  - which sex is favoured, in each group
#   COLOUR   is evidence   - in which group the bias reached significance
# "Reversed" is an evidence call, so it gets its own mark (an open ring) and
# never appears in the quadrant wording.
# =============================================================================
CLS_CHEMO = "Chemosensory"
CLS_BOTH  = "Non-chemosensory, sex-biased in both states"
CLS_ONE   = "Non-chemosensory, sex-biased in one state"
RETENTION_COLORS = {CLS_ONE: "#D8D8D8", CLS_BOTH: "#4A76C4", CLS_CHEMO: "#B2182B"}
RETENTION_LIM = 10.0

QUAD_TEXT = {
    "x- y+": "Male-biased in virgins,<br>female-biased in mated females",
    "x+ y+": "Female-biased in virgins<br>and in mated females",
    "x- y-": "Male-biased in virgins<br>and in mated females",
    "x+ y-": "Female-biased in virgins,<br>male-biased in mated females",
}


def s1_build_retention(tissue_cross, tissue):
    """Figure 4D data: every gene significant in EITHER sex contrast."""
    if tissue not in tissue_cross:
        return None
    td = tissue_cross[tissue]
    msex_df = td.get("msex_df")
    if msex_df is None:
        return None
    keep = [c for c in ("JoinKey", "log2FoldChange", "padj", "is_sig", "ChemoName",
                        "Name", "GO_Name") if c in td["sex_df"].columns]
    sx = td["sex_df"][keep].copy()
    sx = sx.drop_duplicates("JoinKey").rename(
        columns={"log2FoldChange": "x", "padj": "padj_virgin", "is_sig": "sig_v"})
    mx = msex_df[["JoinKey", "log2FoldChange", "padj", "is_sig"]].copy()
    mx = mx.drop_duplicates("JoinKey").rename(
        columns={"log2FoldChange": "y", "padj": "padj_mated", "is_sig": "sig_m"})
    d = sx.merge(mx, on="JoinKey", how="inner")
    d["x"] = pd.to_numeric(d["x"], errors="coerce")
    d["y"] = pd.to_numeric(d["y"], errors="coerce")
    d = d.dropna(subset=["x", "y"])
    d = d[d["sig_v"] | d["sig_m"]].copy()
    if d.empty:
        return None
    d["reversed"] = d["sig_v"] & d["sig_m"] & (np.sign(d["x"]) != np.sign(d["y"]))
    is_chemo = d["ChemoName"].notna()
    d["Class"] = np.where(is_chemo, CLS_CHEMO,
                 np.where(d["sig_v"] & d["sig_m"], CLS_BOTH, CLS_ONE))
    # clamp the view; nothing is dropped, out-of-range points sit on the boundary
    d["oor"] = (d["x"].abs() > RETENTION_LIM) | (d["y"].abs() > RETENTION_LIM)
    d["px"] = d["x"].clip(-RETENTION_LIM, RETENTION_LIM)
    d["py"] = d["y"].clip(-RETENTION_LIM, RETENTION_LIM)
    d["Quad"] = np.select(
        [(d.x < 0) & (d.y > 0), (d.x > 0) & (d.y > 0),
         (d.x < 0) & (d.y < 0), (d.x > 0) & (d.y < 0)],
        ["x- y+", "x+ y+", "x- y-", "x+ y-"], default="axis")
    # family membership from Additional file 37, which is authoritative
    fam, sym = ex_family_map()
    gk = d["JoinKey"].astype(str)
    d["Family"] = gk.map(fam).fillna("")
    d["Symbol"] = gk.map(sym).fillna("")
    d["is_chemo"] = d["Family"].astype(bool)
    if "Name" not in d.columns:
        d["Name"] = ""
    d["Label"] = np.where(d["Symbol"].astype(bool), d["Symbol"],
                  np.where(d["ChemoName"].notna(), d["ChemoName"].astype(str), gk))
    d["hover"] = (
        "<b>" + d["Label"].astype(str) + "</b><br>" + gk
        + np.where(d["Family"].astype(bool), "  ·  " + d["Family"].astype(str), "")
        + "<br>virgins "   + d["x"].round(2).astype(str)
        + "  (q " + d["padj_virgin"].map(lambda v: f"{v:.1e}" if pd.notna(v) else "NA") + ")"
        + "<br>mated "     + d["y"].round(2).astype(str)
        + "  (q " + d["padj_mated"].map(lambda v: f"{v:.1e}" if pd.notna(v) else "NA") + ")"
        + "<br>shift "     + (d["y"] - d["x"]).round(2).astype(str)
        + np.where(d["reversed"], "<br><b>sex bias reverses after mating</b>", "")
        + np.where(d["Name"].astype(str).str.len() > 0,
                   "<br>" + d["Name"].astype(str).str.slice(0, 60), "")
    )
    return d


def s1_scatter_overlays(fig, d, mark_chemo, search):
    """Optional marks on the 4D / S26 scatter: chemosensory genes, and a search.

    Both are drawn as open rings rather than new colours, so they stack on top
    of whatever the points already encode instead of competing with it.
    """
    if mark_chemo and "is_chemo" in d.columns:
        ch = d[d["is_chemo"]]
        if not ch.empty:
            fig.add_trace(go.Scattergl(
                x=ch["px"], y=ch["py"], mode="markers", name="Chemosensory",
                marker=dict(size=9, color="rgba(0,0,0,0)",
                            line=dict(width=1.2, color="#A60000")),
                text=ch["hover"], hovertemplate="%{text}<extra></extra>"))
    q = (search or "").strip().lower()
    if q:
        hay = (d["JoinKey"].astype(str) + " " + d["Label"].astype(str) + " "
               + d.get("Name", pd.Series("", index=d.index)).astype(str) + " "
               + d.get("Family", pd.Series("", index=d.index)).astype(str)).str.lower()
        hit = d[hay.str.contains(q, regex=False, na=False)]
        if not hit.empty:
            fig.add_trace(go.Scattergl(
                x=hit["px"], y=hit["py"], mode="markers+text",
                name=f"matches “{search.strip()}”",
                marker=dict(size=13, color="rgba(0,0,0,0)",
                            line=dict(width=2.0, color="#00A878"), symbol="circle"),
                text=hit["Label"].astype(str), textposition="top center",
                textfont=dict(size=9, color="#00634A"),
                hovertext=hit["hover"], hovertemplate="%{hovertext}<extra></extra>"))


def s1_fig_retention(d, tissue_label, show_quad_text=False,
                     mark_chemo=False, search=""):
    L = RETENTION_LIM
    fig = go.Figure()
    # shade the two quadrants where the direction of bias differs between groups
    for x0, x1, y0, y1 in [(-L, 0, 0, L), (0, L, -L, 0)]:
        fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1,
                      fillcolor="#F0F0F0", opacity=0.6, line_width=0, layer="below")
    fig.add_shape(type="line", x0=-L, x1=L, y0=-L, y1=L,
                  line=dict(color="#BBBBBB", width=1, dash="dot"), layer="below")
    for ln in (dict(x0=-L, x1=L, y0=0, y1=0), dict(x0=0, x1=0, y0=-L, y1=L)):
        fig.add_shape(type="line", line=dict(color="#555555", width=1), layer="below", **ln)
    # grey backdrop first, chemosensory last, so the genes the paper is about sit on top
    for cls in (CLS_ONE, CLS_BOTH, CLS_CHEMO):
        sub = d[d["Class"] == cls]
        if sub.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=sub["px"], y=sub["py"], mode="markers", name=cls,
            marker=dict(size=4.5, color=RETENTION_COLORS[cls], line=dict(width=0)),
            text=sub["hover"], hovertemplate="%{text}<extra></extra>",
        ))
    rev = d[d["reversed"]]
    if not rev.empty:
        fig.add_trace(go.Scattergl(
            x=rev["px"], y=rev["py"], mode="markers", name="Direction reversed",
            marker=dict(size=10, color="rgba(0,0,0,0)",
                        line=dict(width=1.3, color="#111111")),
            text=rev["hover"], hovertemplate="%{text}<extra></extra>",
        ))
    s1_scatter_overlays(fig, d, mark_chemo, search)
    # quadrant counts; the wording is spelled out only in the first panel
    for key, (qx, qy, ax_, ay) in {
        "x- y+": (-L * 0.96,  L * 0.96, "left", "top"),
        "x+ y+": ( L * 0.96,  L * 0.96, "right", "top"),
        "x- y-": (-L * 0.96, -L * 0.96, "left", "bottom"),
        "x+ y-": ( L * 0.96, -L * 0.96, "right", "bottom"),
    }.items():
        n = int((d["Quad"] == key).sum())
        txt = f"{QUAD_TEXT[key]}<br>n = {n}" if show_quad_text else f"n = {n}"
        fig.add_annotation(x=qx, y=qy, text=txt, showarrow=False,
                           xanchor=ax_, yanchor=ay,
                           font=dict(size=9, color="#666666"), align=ax_)
    fig.update_layout(
        title=dict(text=tissue_label, font=dict(size=13)),
        xaxis=dict(title="Sex bias in virgins (log₂FC, VF vs Vm)",
                   range=[-L, L], zeroline=False),
        yaxis=dict(title="Sex bias in mated females (log₂FC, MF vs Vm)",
                   range=[-L, L], zeroline=False, scaleanchor="x", scaleratio=1),
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0, font=dict(size=9)),
        margin=dict(l=10, r=10, t=35, b=10), height=520,
        plot_bgcolor="white",
    )
    return fig


MFS_COLS = {"Induced in mated females": "#E08214",
            "Reduced in mated females": "#2166AC",
            "Not mating-responsive": "#D4D4D4"}


def s1_fig_mfs(d, tissue_label, up, dn, mark_chemo=False, search=""):
    """Figure S26: the same axes as 4D, coloured by the mating response.

    The diagonal is the null for the mating contrast: on it a gene has the same
    sex bias in both groups, so y = x means no mating response.
    """
    L = RETENTION_LIM
    d = d.copy()
    key = d["JoinKey"].astype(str)
    d["MFS"] = np.where(key.isin(up), "Induced in mated females",
               np.where(key.isin(dn), "Reduced in mated females",
                        "Not mating-responsive"))
    fig = go.Figure()
    fig.add_shape(type="line", x0=-L, x1=L, y0=-L, y1=L,
                  line=dict(color="#555555", width=1.4), layer="below")
    for cls in ("Not mating-responsive", "Reduced in mated females",
                "Induced in mated females"):
        sub = d[d["MFS"] == cls]
        if sub.empty:
            continue
        fig.add_trace(go.Scattergl(
            x=sub["px"], y=sub["py"], mode="markers", name=cls,
            marker=dict(size=4.5, color=MFS_COLS[cls]),
            text=sub["hover"], hovertemplate="%{text}<extra></extra>"))
    s1_scatter_overlays(fig, d, mark_chemo, search)
    fig.add_annotation(x=-L*0.55, y=-L*0.55, text="y = x, no mating response",
                       showarrow=False, textangle=-45,
                       font=dict(size=9, color="#555555"), yshift=10)
    fig.update_layout(
        title=dict(text=tissue_label, font=dict(size=13)),
        xaxis=dict(title="Sex bias in virgins (log₂FC)", range=[-L, L]),
        yaxis=dict(title="Sex bias in mated females (log₂FC)", range=[-L, L],
                   scaleanchor="x", scaleratio=1),
        legend=dict(orientation="h", yanchor="top", y=-0.18, x=0, font=dict(size=9)),
        margin=dict(l=10, r=10, t=35, b=10), height=520, plot_bgcolor="white")
    return fig, d


def s1_build_scatter_sex_mating(tissue_cross, tissue):
    """Scatter of LFC_sex vs LFC_mating for genes significant in BOTH contrasts."""
    if tissue not in tissue_cross:
        return None, None
    td = tissue_cross[tissue]
    ov = td["overlap"]
    if not ov:
        return None, None
    ov_str = {str(g) for g in ov}
    sex_df = td["sex_df"].copy()
    mat_df = td["mat_df"].copy()
    sex_sub = sex_df[sex_df["JoinKey"].astype(str).isin(ov_str)].drop_duplicates("JoinKey")
    mat_sub = mat_df[mat_df["JoinKey"].astype(str).isin(ov_str)].drop_duplicates("JoinKey")[
        ["JoinKey", "log2FoldChange", "padj"]
    ].rename(columns={"log2FoldChange": "LFC_mating", "padj": "padj_mating"})
    merged = sex_sub.merge(mat_sub, on="JoinKey", how="inner")
    merged["LFC_sex"] = pd.to_numeric(merged["log2FoldChange"], errors="coerce")
    merged["LFC_mating"] = pd.to_numeric(merged["LFC_mating"], errors="coerce")
    merged = merged.dropna(subset=["LFC_sex", "LFC_mating"])
    if len(merged) < 2:
        return None, None
    merged["PlotDomain"] = s1_assign_domain_chemo_priority(merged)
    r_val = float(merged["LFC_sex"].corr(merged["LFC_mating"]))
    return merged, r_val


def s1_fig_scatter(merged, r_val, tissue_label):
    """Build the LFC_sex vs LFC_mating scatter (Figure 4E)."""
    ax = max(merged["LFC_sex"].abs().max(), merged["LFC_mating"].abs().max()) * 1.15
    ax = max(float(ax), 2.5)

    n_conc = int(((merged["LFC_sex"] > 0) == (merged["LFC_mating"] > 0)).sum())
    n_disc = len(merged) - n_conc

    fig = go.Figure()
    # Quadrant fills
    for x0, x1, y0, y1, col in [
        (0, ax, 0, ax, "#FFE6F0"),    # Q1 female+mating-induced
        (-ax, 0, -ax, 0, "#E6EEFF"),  # Q3 male+mating-suppressed
        (0, ax, -ax, 0, "#FFF3E6"),   # Q4 female+mating-suppressed
        (-ax, 0, 0, ax, "#FFF3E6"),   # Q2 male+mating-induced
    ]:
        fig.add_shape(type="rect", x0=x0, x1=x1, y0=y0, y1=y1,
                      fillcolor=col, opacity=0.35, line_width=0, layer="below")
    # Quadrant labels
    qlabs = [
        (ax*0.97, ax*0.97, "right", "top", "Female-biased;<br>mating-induced"),
        (-ax*0.97, ax*0.97, "left", "top", "Male-biased;<br>mating-induced"),
        (-ax*0.97, -ax*0.97, "left", "bottom", "Male-biased;<br>mating-suppressed"),
        (ax*0.97, -ax*0.97, "right", "bottom", "Female-biased;<br>mating-suppressed"),
    ]
    for x, y, xa, ya, txt in qlabs:
        fig.add_annotation(x=x, y=y, text=f"<i style='color:grey;font-size:9px'>{txt}</i>",
                           showarrow=False, xanchor=xa, yanchor=ya)
    # Axis lines
    fig.add_hline(y=0, line=dict(color="black", width=0.7))
    fig.add_vline(x=0, line=dict(color="black", width=0.7))
    # Diagonal reference
    fig.add_shape(type="line", x0=-ax, x1=ax, y0=-ax, y1=ax,
                  line=dict(color="red", width=0.8, dash="dash"))
    # Regression line
    coef = np.polyfit(merged["LFC_sex"], merged["LFC_mating"], 1)
    xs = np.linspace(-ax, ax, 120)
    fig.add_trace(go.Scatter(x=xs, y=np.polyval(coef, xs), mode="lines",
                              line=dict(color="black", width=1.2), showlegend=False, hoverinfo="skip"))
    # Points by domain
    hover_col = "Name" if "Name" in merged.columns else "JoinKey"
    for dom in ["Chemosensory", "MF", "CC", "BP", "Unknown"]:
        mask = merged["PlotDomain"] == dom
        if not mask.any():
            continue
        sub = merged[mask]
        col = GO_COLS_WITH_CHEMO.get(dom, "#888888")
        sz = 9 if dom == "Chemosensory" else 6
        op = 0.9 if dom == "Chemosensory" else 0.65
        fig.add_trace(go.Scattergl(
            x=sub["LFC_sex"], y=sub["LFC_mating"], mode="markers",
            marker=dict(color=col, size=sz, opacity=op),
            name={**DOMAIN_FULL, "Chemosensory": "Chemosensory"}.get(dom, dom),
            hovertext=sub[hover_col].fillna(sub["JoinKey"]).astype(str),
            hovertemplate="<b>%{hovertext}</b><br>LFC sex=%{x:.2f}<br>LFC mating=%{y:.2f}<extra></extra>",
        ))
    # Stats box
    fig.add_annotation(
        x=ax*0.97, y=-ax*0.50, xanchor="right",
        text=f"<i>r</i> = {r_val:.3f}<br>n = {len(merged)}<br>reinforced: {n_conc}<br>reversed: {n_disc}",
        showarrow=False, font=dict(size=11), bgcolor="rgba(255,255,255,0.85)",
        bordercolor="grey", borderwidth=1,
    )
    fig.update_xaxes(range=[-ax, ax], showgrid=False, zeroline=False,
                     title="log₂FC sex (VF/Vm) — positive: female-biased")
    fig.update_yaxes(range=[-ax, ax], showgrid=False, zeroline=False,
                     scaleanchor="x", scaleratio=1,
                     title="log₂FC mating (MF/VF) — positive: mating-induced")
    fig.update_layout(
        title=dict(text=f"<b>{tissue_label}</b> (n={len(merged)}, <i>r</i>={r_val:.3f})",
                   x=0.5, font=dict(size=12)),
        height=520, legend_title_text="GO domain",
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="h", y=-0.15, itemsizing="constant"),
    )
    return fig


def s1_go_pct_for_set(gene_set, source_df, label, show_chemo=False):
    """GO domain % breakdown for a gene set, using annotations already in source_df."""
    if not gene_set:
        return pd.DataFrame()
    gene_set_str = {str(g) for g in gene_set}
    sub = source_df[source_df["JoinKey"].astype(str).isin(gene_set_str)].drop_duplicates(subset=["JoinKey"]).copy()
    if sub.empty:
        return pd.DataFrame()
    total = len(sub)
    has_chemo = show_chemo and "ChemoName" in sub.columns
    if has_chemo:
        sub["GO_Cat"] = np.where(sub["ChemoName"].notna(), "Chemosensory", sub["GO_Domain"].fillna("Unknown"))
        cat_order = ["Chemosensory", "MF", "CC", "BP", "Unknown"]
    else:
        sub["GO_Cat"] = sub["GO_Domain"].fillna("Unknown")
        cat_order = ["MF", "CC", "BP", "Unknown"]
    counts = sub["GO_Cat"].value_counts().reindex(cat_order, fill_value=0)
    domain_full_chemo = {**DOMAIN_FULL, "Chemosensory": "Chemosensory"}
    result = pd.DataFrame({
        "Category": [label] * len(cat_order),
        "GO_Domain": cat_order,
        "N": counts.values,
        "Percent": counts.values / max(total, 1) * 100,
    })
    result["GO_Domain_full"] = result["GO_Domain"].map(domain_full_chemo).fillna(result["GO_Domain"])
    result["TextN"] = "N=" + result["N"].astype(str)
    return result


def s1_compute_sex_mating_overlap(by_tissue):
    result = {}
    for tissue in ["Antenna", "Palp", "Tarsi"]:
        if tissue not in by_tissue:
            continue
        contrast_list = by_tissue[tissue]
        sex_item  = next((x for x in contrast_list if x["cond1"] == "VF" and x["cond2"] == "Vm"), None)
        mat_item  = next((x for x in contrast_list if x["cond1"] == "MF" and x["cond2"] == "VF"), None)
        msex_item = next((x for x in contrast_list if x["cond1"] == "MF" and x["cond2"] == "Vm"), None)
        if sex_item is None or mat_item is None:
            continue
        sex_df = sex_item["df"]
        mat_df = mat_item["df"]
        sex_vf = set(sex_df.loc[sex_df["is_sig"] & (sex_df["Side"] == "VF"), "JoinKey"].astype(str))
        sex_vm = set(sex_df.loc[sex_df["is_sig"] & (sex_df["Side"] == "Vm"), "JoinKey"].astype(str))
        mat_mf = set(mat_df.loc[mat_df["is_sig"] & (mat_df["Side"] == "MF"), "JoinKey"].astype(str))
        mat_vf = set(mat_df.loc[mat_df["is_sig"] & (mat_df["Side"] == "VF"), "JoinKey"].astype(str))
        # The manuscript calls a gene mating-responsive only when it moves in the
        # same direction against BOTH virgin groups, so the mated-female versus
        # virgin-male contrast gates the mated-female versus virgin-female one.
        # Without this gate the app reports larger sets than Figure 4B.
        if msex_item is not None:
            msex_df = msex_item["df"]
            msex_up = set(msex_df.loc[msex_df["is_sig"] & (msex_df["Side"] == "MF"), "JoinKey"].astype(str))
            msex_dn = set(msex_df.loc[msex_df["is_sig"] & (msex_df["Side"] == "Vm"), "JoinKey"].astype(str))
            mat_mf &= msex_up
            mat_vf &= msex_dn
        overlap = (sex_vf | sex_vm) & (mat_mf | mat_vf)
        result[tissue] = {
            "sex_vf": sex_vf, "sex_vm": sex_vm,
            "mat_mf": mat_mf, "mat_vf": mat_vf,
            "overlap": overlap,
            "ov_female": sex_vf & overlap,
            "ov_male": sex_vm & overlap,
            "sex_df": sex_df, "mat_df": mat_df,
            "msex_df": msex_item["df"] if msex_item is not None else None,
        }
    return result


def s1_build_sex_mating_bar(data):
    """Grouped bar: Sex / Mating / Overlap (concordant + reversed split)."""
    # Overlap split: concordant = same LFC sign in both contrasts
    sex_df = data["sex_df"]
    mat_df = data["mat_df"]
    ov = data["overlap"]
    # Determine concordant/reversed using LFC sign
    ov_str = {str(g) for g in ov}
    sex_lfc = dict(zip(
        sex_df["JoinKey"].astype(str),
        pd.to_numeric(sex_df["log2FoldChange"], errors="coerce")
    ))
    mat_lfc = dict(zip(
        mat_df["JoinKey"].astype(str),
        pd.to_numeric(mat_df["log2FoldChange"], errors="coerce")
    ))
    n_conc = sum(1 for g in ov_str
                 if pd.notna(sex_lfc.get(g)) and pd.notna(mat_lfc.get(g))
                 and np.sign(sex_lfc.get(g, 0)) == np.sign(mat_lfc.get(g, 0)))
    n_rev  = len(ov) - n_conc

    rows = [
        {"group": "Sex-biased", "category": "Female-biased (VF>Vm)", "n": len(data["sex_vf"])},
        {"group": "Sex-biased", "category": "Male-biased (Vm>VF)",   "n": len(data["sex_vm"])},
        {"group": "Mating-regulated", "category": "Mating-induced (MF>VF)",    "n": len(data["mat_mf"])},
        {"group": "Mating-regulated", "category": "Mating-suppressed (VF>MF)", "n": len(data["mat_vf"])},
        {"group": "Both contrasts", "category": "Reinforced (same direction)", "n": n_conc},
        {"group": "Both contrasts", "category": "Reversed (opposite direction)","n": n_rev},
    ]
    df_bar = pd.DataFrame(rows)
    color_map = {
        "Female-biased (VF>Vm)":       COL_VF["dark"],
        "Male-biased (Vm>VF)":         COL_Vm["dark"],
        "Mating-induced (MF>VF)":      COL_MF["dark"],
        "Mating-suppressed (VF>MF)":   COL_VF["light"],
        "Reinforced (same direction)":  "#7B1FA2",
        "Reversed (opposite direction)":"#E65100",
    }
    cat_order = list(color_map.keys())
    fig = px.bar(
        df_bar, x="group", y="n", color="category", barmode="group", text="n",
        color_discrete_map=color_map,
        category_orders={
            "group": ["Sex-biased", "Mating-regulated", "Both contrasts"],
            "category": cat_order,
        },
        labels={"group": "", "n": "Number of significant DEGs", "category": ""},
        height=440,
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(template="simple_white", legend_title_text="",
                      legend=dict(orientation="h", y=-0.2))
    y_max = max(1, df_bar["n"].max())
    fig.update_yaxes(range=[0, y_max * 1.35])
    return fig


def s1_make_overlap_download_table(gene_set, sex_df, mat_df):
    if not gene_set:
        return pd.DataFrame()
    gene_set_str = set(str(g) for g in gene_set)
    sex_sub = sex_df[sex_df["JoinKey"].astype(str).isin(gene_set_str)].copy()
    keep_sex = ["JoinKey"]
    if "Name" in sex_sub.columns:
        keep_sex.append("Name")
    keep_sex += [c for c in ["log2FoldChange", "padj"] if c in sex_sub.columns]
    expr_sex_cols = [c for c in sex_sub.columns if c.startswith("Mean ")]
    go_cols = [c for c in sex_sub.columns if c in ("GO_Domain", "GO_Name")]
    keep_sex += expr_sex_cols + go_cols
    sex_out = sex_sub[[c for c in keep_sex if c in sex_sub.columns]].rename(
        columns={"log2FoldChange": "log2FC_sex", "padj": "padj_sex"}
    )
    mat_sub = mat_df[mat_df["JoinKey"].astype(str).isin(gene_set_str)].copy()
    keep_mat = ["JoinKey", "log2FoldChange", "padj"]
    expr_mat_cols = [c for c in mat_sub.columns if c.startswith("Mean ") and c not in expr_sex_cols]
    keep_mat += expr_mat_cols
    mat_out = mat_sub[[c for c in keep_mat if c in mat_sub.columns]].rename(
        columns={"log2FoldChange": "log2FC_mating", "padj": "padj_mating"}
    )
    merged = sex_out.merge(mat_out, on="JoinKey", how="outer")
    return merged.sort_values("padj_sex", na_position="last").reset_index(drop=True)


def render_s1_tab():
    base_dir = S1_DEFAULT_PATHS["BASE_DIR"]

    with st.sidebar:
        st.header("Sex & State Analysis")
        padj_thr = st.sidebar.number_input(
            "padj threshold", value=0.001, min_value=1e-10, max_value=0.1,
            format="%.4g", key="s1_padj_thr"
        )
        lfc_thr = st.sidebar.number_input(
            "|log2FC| threshold ", value=1.0, min_value=0.0, max_value=10.0,
            step=0.25, key="s1_lfc_thr"
        )
        strong_lfc = st.sidebar.number_input(
            "Strong |log2FC|", value=2.5, min_value=0.0, max_value=10.0,
            step=0.25, key="s1_strong_lfc"
        )
        search_term = st.sidebar.text_input(
            "Search gene ", value="",
            placeholder="e.g. OR120 or XM_012345678",
            key="s1_search_term",
        )
        color_by_go = st.sidebar.checkbox(
            "Color by GO domain ", value=False, key="s1_color_by_go"
        )
        overlay_chemo = st.sidebar.checkbox(
            "Overlay chemosensory genes", value=True, key="s1_overlay_chemo"
        )

    by_tissue, norm_means = s1_load_all_contrasts(base_dir, padj_thr, lfc_thr, strong_lfc)
    master_annot = s1_build_master_annotation(by_tissue, norm_means)

    if search_term.strip():
        q = search_term.strip().lower()
        found_any = False
        for tissue, contrast_list in by_tissue.items():
            for item in contrast_list:
                df_all = item["df"]
                cols = [c for c in ["Gene", "Name", "JoinKey", "ChemoName"] if c in df_all.columns]
                for c in cols:
                    s = df_all[c].astype(str).str.strip().str.lower()
                    if (s == q).any():
                        found_any = True
                        break
                if found_any:
                    break
            if found_any:
                break
        if not found_any:
            st.sidebar.warning(f"'{search_term}' not found in any loaded contrast.")

    # Compute cross-contrast overlap once; shared across Sex×Mating, Fig 4D/E, and GO Names tabs
    tissue_cross = s1_compute_sex_mating_overlap(by_tissue)
    _go_path_s1 = os.path.join(base_dir, S1_DEFAULT_PATHS["FILE_GO"])
    go_map_s1 = s1_load_go_table(_go_path_s1)

    tab_volc, tab_bar, tab_venn, tab_4de, tab_gonames = st.tabs(
        ["Volcano plots (Fig. S25)", "DE summary bars + DEG tables",
         "Venn diagrams (Fig. 4C)", "Sex-bias reorganisation (Figs. 4D, S26)",
         "GO Term Browser (Figs. S22–S24)"]
    )

    tissue_labels_s1 = {"Antenna": "Antenna", "Palp": "Maxillary palp", "Tarsi": "Tarsi"}

    # -------------------- VOLCANO TAB --------------------
    with tab_volc:
        st.subheader("Volcano plots (three per appendage)")
        st.caption("Manuscript: **Figure S25**")
        st.markdown(
            "Figure 4A and 4B of the manuscript are lollipop plots; the volcanoes were "
            "moved to Figure S25 because they show the full distribution behind the "
            "counts, which the lollipops cannot. Per-gene responses are in "
            "**Data Explorer → Response (lollipops)**. The third panel, mated female "
            "versus virgin male, is the second of the two contrasts a gene must pass "
            "to count as mating-responsive."
        )
        for tissue in ["Antenna", "Palp", "Tarsi"]:
            if tissue not in by_tissue:
                continue
            st.markdown(f"### {TISSUE_DISPLAY[tissue]}")
            contrasts = sorted(by_tissue[tissue], key=s1_sort_key)
            cols = st.columns(len(contrasts))
            for i, item in enumerate(contrasts):
                cond1, cond2, df_cls = item["cond1"], item["cond2"], item["df"]
                with cols[i]:
                    st.caption(s1_de_counts_caption(df_cls, cond1, cond2))
                    fig = s1_build_volcano_figure(
                        df_cls, cond1, cond2, padj_thr, lfc_thr,
                        search_term=search_term, overlay_chemo=overlay_chemo,
                        color_by_go=color_by_go,
                    )
                    st.plotly_chart(fig, use_container_width=True,
                                    key=f"s1_volc_{tissue}_{cond1}_{cond2}")

    # -------------------- BARS TAB AND DEG TABLES --------------------
    with tab_bar:
        st.subheader("Upregulated DEGs summary and DEG tables")
        for tissue in ["Antenna", "Palp", "Tarsi"]:
            if tissue not in by_tissue:
                continue
            st.markdown(f"### {TISSUE_DISPLAY[tissue]}")
            contrasts = sorted(by_tissue[tissue], key=s1_sort_key)
            cols = st.columns(len(contrasts))
            for i, item in enumerate(contrasts):
                cond1, cond2, df_cls = item["cond1"], item["cond2"], item["df"]
                with cols[i]:
                    label1 = COND_FULL.get(cond1, cond1)
                    label2 = COND_FULL.get(cond2, cond2)
                    col1_c = s1_get_full_label_color(label1)
                    col2_c = s1_get_full_label_color(label2)
                    title_html = (
                        f"<span style='color:{col2_c}; font-weight:bold'>{label2}</span> "
                        f"vs "
                        f"<span style='color:{col1_c}; font-weight:bold'>{label1}</span>"
                    )
                    st.markdown(title_html, unsafe_allow_html=True)
                    fig_bar = s1_build_bar_figure(df_cls, cond1, cond2)
                    if fig_bar is None:
                        st.info("No significant genes for current thresholds.")
                    else:
                        st.plotly_chart(fig_bar, use_container_width=True,
                                        key=f"s1_bar_{tissue}_{cond1}_{cond2}")
                    with st.expander("Show DEG table for this contrast"):
                        deg = df_cls[df_cls["is_sig"]].copy()
                        if deg.empty:
                            st.info("No DEGs for this contrast.")
                        else:
                            deg["Up_state"] = deg["Side"].map(COND_FULL)
                            state_options = ["All states"]
                            for lab in [label1, label2]:
                                if lab not in state_options:
                                    state_options.append(lab)
                            state_choice = st.selectbox(
                                "Filter by sex / physiological state (upregulated side)",
                                state_options, key=f"s1_state_{tissue}_{cond1}_{cond2}",
                            )
                            if state_choice != "All states":
                                deg = deg[deg["Up_state"] == state_choice]
                            strength_choice = st.selectbox(
                                "Filter by DEG strength",
                                ["All strengths", "Moderate only", "Strong only"],
                                key=f"s1_strength_{tissue}_{cond1}_{cond2}",
                            )
                            if strength_choice == "Moderate only":
                                deg = deg[deg["Strength"] == "Moderate"]
                            elif strength_choice == "Strong only":
                                deg = deg[deg["Strength"] == "Strong"]
                            if deg.empty:
                                st.info("No DEGs for the selected filters.")
                            else:
                                col_expr1 = f"Mean {COND_FULL.get(cond1, cond1)} ({TISSUE_DISPLAY[tissue]})"
                                col_expr2 = f"Mean {COND_FULL.get(cond2, cond2)} ({TISSUE_DISPLAY[tissue]})"
                                columns_to_show = []
                                for c in ["JoinKey", "Gene", "Name"]:
                                    if c in deg.columns:
                                        columns_to_show.append(c)
                                if col_expr1 in deg.columns:
                                    columns_to_show.append(col_expr1)
                                if col_expr2 in deg.columns:
                                    columns_to_show.append(col_expr2)
                                for c in ["Up_state", "Strength", "log2FoldChange", "padj",
                                          "GO_Domain", "GO_Name"]:
                                    if c in deg.columns:
                                        columns_to_show.append(c)
                                if not columns_to_show:
                                    if "padj" in deg.columns:
                                        st.dataframe(deg.style.format({"padj": "{:.1e}"}))
                                    else:
                                        st.dataframe(deg)
                                else:
                                    df_to_show = deg[columns_to_show]
                                    if "padj" in df_to_show.columns:
                                        st.dataframe(df_to_show.style.format({"padj": "{:.1e}"}))
                                    else:
                                        st.dataframe(df_to_show)

    # -------------------- VENN TAB --------------------
    with tab_venn:
        st.subheader("Genes shared among appendages")
        st.caption("Manuscript: **Figure 4C**")
        st.markdown(
            "Four transcriptional programmes, each as a three-way Venn across the "
            "appendages. Sex bias is measured in virgins. A gene counts as "
            "mating-responsive only when it moves in the same direction against "
            "**both** virgin groups, which is the manuscript's definition."
        )
        # colours as in the manuscript: virgin male, virgin female, mated female
        PROGRAMMES = [
            ("male_biased",  "Male-biased in virgins",
             "significant in virgin male vs virgin female", "#003399", "sex_vm"),
            ("female_biased", "Female-biased in virgins",
             "significant in virgin female vs virgin male", "#99004D", "sex_vf"),
            ("mating_down",  "Reduced in mated females",
             "lower in mated females than both virgin groups", "#2471A3", "mat_vf"),
            ("mating_up",    "Induced in mated females",
             "higher in mated females than both virgin groups", "#7B1FA2", "mat_mf"),
        ]
        region_options = [
            "Antenna only", "Maxillary palp only", "Tarsi only",
            "Antenna & Maxillary palp", "Antenna & Tarsi",
            "Maxillary palp & Tarsi", "All three appendages",
        ]
        region_key_map = {
            "Antenna only": "A only", "Maxillary palp only": "B only",
            "Tarsi only": "C only", "Antenna & Maxillary palp": "A & B",
            "Antenna & Tarsi": "A & C", "Maxillary palp & Tarsi": "B & C",
            "All three appendages": "A & B & C",
        }
        summary = []
        grid = [st.columns(2), st.columns(2)]
        for i, (pid, title, sub, colour, setkey) in enumerate(PROGRAMMES):
            cell = grid[i // 2][i % 2]
            with cell:
                st.markdown(f"**{title}**")
                st.caption(sub)
                sets = {}
                for t in ["Antenna", "Palp", "Tarsi"]:
                    td = tissue_cross.get(t)
                    sets[t] = {str(g) for g in td[setkey]} if td else set()
                A, B, C = sets["Antenna"], sets["Palp"], sets["Tarsi"]
                regions = s1_compute_venn_regions(A, B, C)
                counts = {k: len(v) for k, v in regions.items()}
                fig = s1_build_venn_figure(
                    title, colour, ["Antenna", "Maxillary palp", "Tarsi"], counts)
                st.plotly_chart(fig, use_container_width=False, key=f"s1_venn_{pid}")
                summary.append({
                    "Programme": title, "Antenna": len(A), "Maxillary palp": len(B),
                    "Tarsi": len(C), **{k: counts[v] for k, v in region_key_map.items()}})
                choice = st.selectbox("Region to list", region_options,
                                      key=f"s1_venn_region_{pid}")
                genes = sorted(regions[region_key_map[choice]])
                st.caption(f"{len(genes):,} genes in {choice.lower()}")
                if genes:
                    gdf = pd.DataFrame({"Gene": genes})
                    nm = tissue_cross.get("Antenna", {}).get("sex_df")
                    if nm is not None and "ChemoName" in nm.columns:
                        lut = (nm.drop_duplicates("JoinKey")
                                 .set_index(nm.drop_duplicates("JoinKey")["JoinKey"].astype(str))
                                 ["ChemoName"])
                        gdf["Chemosensory name"] = gdf["Gene"].map(lut)
                    st.download_button(
                        "Download these genes (CSV)",
                        gdf.to_csv(index=False).encode(),
                        file_name=f"Figure_4C_{pid}_{region_key_map[choice].replace(' ','')}.csv",
                        mime="text/csv", key=f"s1_venn_dl_{pid}")
        if summary:
            st.markdown("---")
            st.dataframe(pd.DataFrame(summary), use_container_width=True, hide_index=True)
            st.download_button(
                "Download all region counts (CSV)",
                pd.DataFrame(summary).to_csv(index=False).encode(),
                file_name="Figure_4C_regions.csv", mime="text/csv", key="s1_venn_dl_all")

    with tab_4de:
        if not tissue_cross:
            st.warning("No overlap data — ensure both contrasts are loaded for each tissue.")
        else:
            # ── Figure 4D: reorganisation of sex bias by mating ─────────────
            st.subheader("Reorganisation of sex bias by mating")
            st.caption("Manuscript: **Figure 4D**")
            st.markdown(
                "Sex bias measured in virgins (x, virgin female vs virgin male) against sex bias "
                "measured in mated females (y, mated female vs virgin male), for every gene "
                "significant in **either** contrast. Both axes are female-minus-male against the "
                "same reference group, virgin males, so a change of sign between them is a genuine "
                "reversal of sex bias rather than an artefact of a shared denominator. "
                "Points on the dotted diagonal are unchanged by mating; vertical distance from it "
                "is the mating effect. Shaded quadrants are those where the direction of bias "
                "differs between the two groups. Axes are clipped at ±10 log₂FC and any gene "
                "beyond is drawn on the boundary."
            )
            view = st.radio(
                "Colour points by",
                ["Sex-bias class (Figure 4D)", "Mating response (Figure S26)"],
                horizontal=True, key="s1_4d_view")
            as_s26 = view.startswith("Mating")
            o1, o2 = st.columns([1, 2])
            with o1:
                mark_chemo = st.checkbox("Mark chemosensory genes", value=False,
                                         key="s1_4d_chemo")
            with o2:
                scatter_q = st.text_input(
                    "Find genes in the scatter", "", key="s1_4d_search",
                    placeholder="ORco, OBP34, XM_038066019.1, or a family such as OR")
            st.caption(
                "Hover any point for its accession, family, both fold changes with "
                "their adjusted p values, and the shift between them. Matches to the "
                "search are ringed in green and labelled."
            )
            ret_rows = {}
            ret_cols = st.columns(3)
            for ci, tissue in enumerate(["Antenna", "Palp", "Tarsi"]):
                d_ret = s1_build_retention(tissue_cross, tissue)
                ret_rows[tissue] = d_ret
                with ret_cols[ci]:
                    if d_ret is None:
                        st.info(f"Mated-female contrast not loaded for {TISSUE_DISPLAY[tissue]}.")
                    elif as_s26:
                        td = tissue_cross[tissue]
                        fig_s26, _ = s1_fig_mfs(d_ret, TISSUE_DISPLAY[tissue],
                                                td["mat_mf"], td["mat_vf"],
                                                mark_chemo, scatter_q)
                        st.plotly_chart(fig_s26, use_container_width=True,
                                        key=f"s1_mfs_{tissue}")
                    else:
                        st.plotly_chart(
                            s1_fig_retention(d_ret, TISSUE_DISPLAY[tissue],
                                             show_quad_text=(ci == 0),
                                             mark_chemo=mark_chemo, search=scatter_q),
                            use_container_width=True, key=f"s1_retention_{tissue}")
            if scatter_q.strip():
                hits = []
                for tissue, d_ret in ret_rows.items():
                    if d_ret is None:
                        continue
                    hay = (d_ret["JoinKey"].astype(str) + " " + d_ret["Label"].astype(str)
                           + " " + d_ret["Family"].astype(str)).str.lower()
                    h = d_ret[hay.str.contains(scatter_q.strip().lower(), regex=False, na=False)]
                    if not h.empty:
                        hits.append(h.assign(Appendage=TISSUE_DISPLAY[tissue])[
                            ["Appendage", "JoinKey", "Label", "Family", "x", "y",
                             "padj_virgin", "padj_mated", "reversed"]])
                if hits:
                    hv = pd.concat(hits, ignore_index=True).rename(columns={
                        "JoinKey": "Gene", "Label": "Name", "x": "log₂FC virgins",
                        "y": "log₂FC mated females", "padj_virgin": "q virgins",
                        "padj_mated": "q mated", "reversed": "Reversed"})
                    st.dataframe(hv.round(3), use_container_width=True, hide_index=True)
                    st.download_button(
                        "Download these matches (CSV)", hv.to_csv(index=False).encode(),
                        file_name="Figure_4D_search_matches.csv", mime="text/csv",
                        key="s1_4d_search_dl")
                else:
                    st.info(f"Nothing in the scatter matches “{scatter_q.strip()}”.")
            summary = []
            for tissue, d_ret in ret_rows.items():
                if d_ret is None:
                    continue
                summary.append({
                    "Appendage": TISSUE_DISPLAY[tissue],
                    "Sex-biased in virgins": int(d_ret["sig_v"].sum()),
                    "Sex-biased in mated females": int(d_ret["sig_m"].sum()),
                    "Significant in both": int((d_ret["sig_v"] & d_ret["sig_m"]).sum()),
                    "Direction reversed": int(d_ret["reversed"].sum()),
                    "Chemosensory": int((d_ret["Class"] == CLS_CHEMO).sum()),
                })
            if summary:
                st.dataframe(pd.DataFrame(summary), use_container_width=True, hide_index=True)
                rev_all = pd.concat(
                    [d.assign(Appendage=TISSUE_DISPLAY[t]) for t, d in ret_rows.items()
                     if d is not None and d["reversed"].any()],
                    ignore_index=True) if any(
                    d is not None and d["reversed"].any() for d in ret_rows.values()) else None
                if rev_all is not None:
                    with st.expander("Genes whose sex bias reverses after mating"):
                        rv = rev_all[rev_all["reversed"]][
                            ["Appendage", "JoinKey", "ChemoName", "x", "y"]].copy()
                        rv.columns = ["Appendage", "Gene", "Chemosensory name",
                                      "log₂FC virgins", "log₂FC mated females"]
                        st.dataframe(rv.round(2), use_container_width=True, hide_index=True)

            st.markdown("---")
            # ── Figure 4D: Venns ────────────────────────────────────────────
            st.subheader("Overlap between sex-biased and mating-responsive gene sets")
            st.caption(
                "Within each appendage. This was a panel of Figure 4 in an earlier "
                "version and is kept here because it is the clearest view of the "
                "overlap; the current Figure 4D is the scatter above."
            )
            st.markdown("**Total** (all genes) and **per GO domain** (with chemosensory priority).")
            for tissue in ["Antenna", "Palp", "Tarsi"]:
                if tissue not in tissue_cross:
                    continue
                td = tissue_cross[tissue]
                all_sex = td["sex_vf"] | td["sex_vm"]
                all_mat = td["mat_mf"] | td["mat_vf"]
                ov = td["overlap"]
                sex_only = all_sex - ov
                mat_only = all_mat - ov

                st.markdown(f"### {TISSUE_DISPLAY[tissue]}")
                # Total Venn
                vcol_total, vcol_domains = st.columns([1, 3])
                with vcol_total:
                    st.caption("**All genes**")
                    fig_vt = s1_fig_venn_sex_mating(
                        len(sex_only), len(mat_only), len(ov), TISSUE_DISPLAY[tissue]
                    )
                    st.plotly_chart(fig_vt, use_container_width=False,
                                    key=f"s1_4d_total_{tissue}")
                with vcol_domains:
                    st.caption("**Per GO domain** (chemosensory priority)")
                    dom_venns = s1_compute_domain_venns(tissue_cross, tissue)
                    dom_cols = st.columns(5)
                    for di, dom in enumerate(["Chemosensory", "MF", "CC", "BP", "Unknown"]):
                        ns, nm, no = dom_venns.get(dom, (0, 0, 0))
                        col_d = GO_COLS_WITH_CHEMO.get(dom, GO_COLS.get(dom, "#888888"))
                        with dom_cols[di]:
                            fig_dv = s1_fig_domain_venn(ns, nm, no, dom, col_d)
                            st.plotly_chart(fig_dv, use_container_width=False,
                                            key=f"s1_4d_dom_{tissue}_{dom}")

                with st.expander(f"Download gene sets — {TISSUE_DISPLAY[tissue]}"):
                    region_opts = {
                        "Sex-only (sex DEG, not mating)": sex_only,
                        "Mating-only (mating DEG, not sex)": mat_only,
                        "Overlap (DEG in both contrasts)": ov,
                    }
                    rc = st.selectbox("Select region", list(region_opts.keys()),
                                      key=f"s1_4d_region_{tissue}")
                    gene_set_4d = region_opts[rc]
                    if not gene_set_4d:
                        st.info("No genes in this region for current thresholds.")
                    else:
                        dl4d = s1_make_overlap_download_table(gene_set_4d, td["sex_df"], td["mat_df"])
                        st.write(f"{len(dl4d)} genes")
                        fmt4d = {c: "{:.2e}" for c in ("padj_sex", "padj_mating") if c in dl4d.columns}
                        st.dataframe(dl4d.style.format(fmt4d) if fmt4d else dl4d, use_container_width=True)
                        csv4d = dl4d.to_csv(index=False).encode("utf-8")
                        safe_4dt = tissue.replace(" ", "_")
                        safe_4dr = rc.replace(" ", "_").replace("(", "").replace(")", "")[:30]
                        st.download_button("Download CSV", data=csv4d,
                                           file_name=f"BSF_{safe_4dt}_{safe_4dr}.csv",
                                           mime="text/csv", key=f"s1_4d_dl_{tissue}_{safe_4dr}")

    # -------------------- GO NAMES TAB (S22–S24 equivalent) --------------------
    with tab_gonames:
        st.subheader("GO term enrichment browser")
        st.markdown(
            "Browse the top enriched GO terms for each gene set from the Venn diagrams. "
            "Supplementary Figures S22 (sex-biased only), S23 (mating-responsive only), S24 (overlap)."
        )
        if not tissue_cross:
            st.warning("No data — ensure both contrasts are loaded for each tissue.")
        else:
            drop_unk_gn = True  # always drop genes with no GO annotation
            top_n_gn = int(st.number_input("Top N GO names per domain", min_value=1,
                                            max_value=50, value=20, key="s1_gn_top_n"))

            for tissue in ["Antenna", "Palp", "Tarsi"]:
                if tissue not in tissue_cross:
                    continue
                td = tissue_cross[tissue]
                all_sex = td["sex_vf"] | td["sex_vm"]
                all_mat = td["mat_mf"] | td["mat_vf"]
                ov = td["overlap"]
                sex_only_gn = all_sex - ov
                mat_only_gn = all_mat - ov

                st.markdown(f"### {TISSUE_DISPLAY[tissue]}")
                reg_tabs = st.tabs(["Sex-biased only (S22)", "Mating-responsive only (S23)", "Sex & mating overlap (S24)"])

                for ri, (region_label, region_keys, rtab) in enumerate([
                    ("Sex-only", sex_only_gn, reg_tabs[0]),
                    ("Mating-only", mat_only_gn, reg_tabs[1]),
                    ("Overlap", ov, reg_tabs[2]),
                ]):
                    with rtab:
                        if not region_keys:
                            st.info(f"No {region_label.lower()} genes for current thresholds.")
                            continue
                        bar_cs = st.columns(3)
                        gn_tbls = {}
                        x_max_gn = 1
                        for dom in ["MF", "CC", "BP"]:
                            t_gn = c1_go_name_table_overlap(region_keys, go_map_s1, dom, drop_unk_gn).head(top_n_gn)
                            gn_tbls[dom] = t_gn
                            if not t_gn.empty:
                                x_max_gn = max(x_max_gn, int(t_gn["N"].max()))
                        for j, dom in enumerate(["MF", "CC", "BP"]):
                            with bar_cs[j]:
                                t_gn = gn_tbls[dom]
                                if t_gn.empty:
                                    st.caption(f"No {DOMAIN_FULL.get(dom, dom)} terms")
                                else:
                                    fig_gn = c1_fig_go_names(t_gn, dom, x_max=x_max_gn)
                                    st.plotly_chart(fig_gn, use_container_width=True,
                                                    key=f"s1_gn_{tissue}_{ri}_{dom}")
                        dl_gn_cs = st.columns(3)
                        for j, dom in enumerate(["MF", "CC", "BP"]):
                            t_gn = gn_tbls[dom]
                            if not t_gn.empty:
                                with dl_gn_cs[j]:
                                    st.download_button(
                                        f"Download {dom}",
                                        data=t_gn.to_csv(index=False).encode("utf-8"),
                                        file_name=f"BSF_{tissue.replace(' ','_')}_{region_label.replace('-','_')}_{dom}_GO_names.csv",
                                        mime="text/csv",
                                        key=f"s1_gn_dl_{tissue}_{ri}_{dom}",
                                    )
                st.markdown("---")


# =============================================================================
# === APP_H1 FUNCTIONS ===
# =============================================================================

def h1_detect_species(label: str) -> str:
    if label is None:
        return "Hill"
    for prefix in ["Aaeg", "Dmel", "Mdom", "Csty", "Bdor", "Gmor", "Hill"]:
        if label.startswith(prefix):
            return prefix
    return "Hill"


def h1_add_hill_prefix(label: str) -> str:
    if label is None:
        return label
    species = h1_detect_species(label)
    if species == "Hill":
        if label.startswith("Hill"):
            return label
        return f"Hill{label}"
    return label


@st.cache_data
def h1_load_matrix(unit_choice: str) -> pd.DataFrame:
    if not H1_NORM_COUNTS_FILE.exists():
        return pd.DataFrame()

    norm_counts = pd.read_csv(H1_NORM_COUNTS_FILE, header=0)
    norm_counts["JoinKey"] = make_joinkey(norm_counts)

    if unit_choice == "FPKM":
        gene_lengths = t1_derive_gene_lengths_from_counts(norm_counts)
        fpkm_matrix = t1_compute_fpkm_matrix(norm_counts, gene_lengths)
        if not fpkm_matrix.empty:
            if "Name" in norm_counts.columns:
                fpkm_matrix = fpkm_matrix.merge(
                    norm_counts[["JoinKey", "Name"]].drop_duplicates(),
                    on="JoinKey", how="left"
                )
            return fpkm_matrix

    return norm_counts


@st.cache_data
def h1_load_norm_counts(unit_choice: str = "Normalized counts") -> pd.DataFrame:
    df = h1_load_matrix(unit_choice)
    if df.empty:
        return df

    def row_mean_safe(row, cols):
        vals = [pd.to_numeric(row[c], errors="coerce") for c in cols if c in row.index]
        vals = [v for v in vals if pd.notna(v)]
        if not vals:
            return 0.0
        return float(sum(vals)) / float(len(vals))

    for base in ["Ant_MF", "Ant_VF", "Ant_Vm", "Leg_MF", "Leg_VF", "Leg_Vm",
                 "P_MF", "P_VF", "P_Vm"]:
        cols = [f"{base}{i}" for i in [1, 2, 3]]
        df[f"{base}_mean"] = df.apply(lambda r, cc=cols: row_mean_safe(r, cc), axis=1)

    df["gene_fam"] = "nan"
    name_col = "Name" if "Name" in df.columns else "name"

    def assign_family(name: str) -> str:
        if "OR" in name: return "OR"
        if "GR" in name: return "GR"
        if "IR" in name: return "IR"
        if "OBP" in name: return "OBP"
        if "CSP" in name: return "CSP"
        if "PPK" in name: return "PPK"
        return "nan"

    names = df[name_col].astype("object").where(df[name_col].notna(), "").astype(str)
    df["gene_fam"] = names.apply(assign_family)
    df["Name_prefixed"] = names.apply(h1_add_hill_prefix)
    return df

def h1_make_chemo_heatmap_matplotlib(
    df: pd.DataFrame,
    gene_fam: str,
    tissues_keep: List[str],
    states_keep: List[str],
    expr_thr: float,
    search_term: str = "",
    unit_choice: str = "Normalized counts",
):
    if df.empty:
        return None

    fam = gene_fam.upper().strip()
    sub = df[df["gene_fam"].astype(str).str.upper() == fam].copy()
    if sub.empty:
        return None

    q = (search_term or "").strip().lower()
    if q:
        name_series = sub["Name_prefixed"].astype(str)
        mask = (
            name_series.str.lower().str.contains(q)
            | name_series.str.lower().str.replace("hill", "", regex=False).str.contains(q)
        )
        sub = sub[mask]
        if sub.empty:
            return None

    tissue_to_prefix = {"Antenna": "Ant", "Maxillary palp": "P", "Tarsi": "Leg"}
    state_order = ["MF", "VF", "Vm"]
    state_full = {
        "MF": "Mated Female",
        "VF": "Virgin Female",
        "Vm": "Virgin Male",
    }
    state_colors = {
        "MF": "#800080",
        "VF": "#FFC0CB",
        "Vm": "#ADD8E6",
    }

    cols = []
    col_labels = []
    col_states = []

    for tissue_label in ["Antenna", "Maxillary palp", "Tarsi"]:
        if tissue_label not in tissues_keep:
            continue
        prefix = tissue_to_prefix[tissue_label]
        for st_code in state_order:
            if st_code not in states_keep:
                continue
            c = f"{prefix}_{st_code}_mean"
            if c in sub.columns:
                cols.append(c)
                col_labels.append(f"{tissue_label}\n{state_full[st_code]}")
                col_states.append(st_code)

    if not cols:
        return None

    mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    mat_log = np.log2(mat + 1.0)

    thr_log = float(np.log2(expr_thr + 1.0))
    mat_centered = mat_log - thr_log

    row_order = np.argsort(np.nanmax(mat_centered, axis=1))[::-1]
    sub = sub.iloc[row_order].copy()
    mat_centered = mat_centered[row_order, :]

    cmap = LinearSegmentedColormap.from_list("bwr_custom", ["#1f4aa8", "#ffffff", "#b2182b"])
    vmax = float(np.nanmax(mat_centered)) if np.isfinite(np.nanmax(mat_centered)) else 1.0
    vmin = float(np.nanmin(mat_centered)) if np.isfinite(np.nanmin(mat_centered)) else -1.0
    lim = max(abs(vmin), abs(vmax), 1e-6)
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)

    n_rows = mat_centered.shape[0]
    n_cols = mat_centered.shape[1]
    fig_w = 6.5 + (n_cols * 0.7)
    fig_h = 4.0 + min(16.0, n_rows * 0.20)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=140)

    im = ax.imshow(mat_centered, aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(col_labels, fontsize=9)
    for tick, st_code in zip(ax.get_xticklabels(), col_states):
        tick.set_color(state_colors.get(st_code, "black"))

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(sub["Name_prefixed"].astype(str).tolist(), fontsize=8)

    ax.set_title(
        f"{fam} expression ({unit_choice}; log2(x+1) centered at threshold {expr_thr})",
        fontsize=12,
        pad=10,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(f"log2({unit_choice}+1) - log2(threshold+1)", rotation=90)

    fig.tight_layout()
    return fig


def render_h1_tab():
    st.subheader("Chemosensory heatmap")
    st.caption("Manuscript: **Figs. S13–S19** (per gene family, A4 PDF — OR, GR, IR, OBP, PPK, CSP, TRP)")

    unit_choice = st.radio(
        "Expression unit",
        ["Normalized counts", "FPKM"],
        index=0,
        key="h1_unit_choice",
    )

    df_norm = h1_load_norm_counts(unit_choice)
    if df_norm.empty:
        st.info(
            "No normalized counts file found.\n\n"
            f"Expected at: `{H1_NORM_COUNTS_FILE}`.\n"
            "Use the same file you used in R: BSF_normalized_counts_nameX.csv."
        )
        return

    fam_choice = st.selectbox(
        "Gene family",
        ["OR", "GR", "IR", "OBP", "CSP", "PPK", "TRP"],
        index=0,
        key="h1_fam_choice",
    )

    tissue_choices = st.multiselect(
        "Tissues (columns)",
        ["Antenna", "Maxillary palp", "Tarsi"],
        default=["Antenna", "Maxillary palp", "Tarsi"],
        key="h1_tissues_keep",
    )

    state_choices = st.multiselect(
        "Sex/state (columns)",
        ["MF", "VF", "Vm"],
        default=["MF", "VF", "Vm"],
        key="h1_states_keep",
    )

    expr_thr = st.number_input(
        f"Expression threshold in {unit_choice} (used for blue-white-red centering)",
        min_value=0.0,
        max_value=1e6,
        value=float(H1_THRESHOLD),
        step=1.0,
        key="h1_expr_thr",
    )

    search_term = st.text_input(
        "Optional search (filters rows, e.g. OR2, GR1, IR, OBP)",
        value="",
        key="h1_search_term2",
    )

    fig = h1_make_chemo_heatmap_matplotlib(
        df=df_norm,
        gene_fam=fam_choice,
        tissues_keep=tissue_choices,
        states_keep=state_choices,
        expr_thr=expr_thr,
        search_term=search_term,
        unit_choice=unit_choice,
    )

    if fig is None:
        st.info("No genes/columns match your current filters.")
        return

    st.pyplot(fig, clear_figure=True)


# =============================================================================
# === PDF DISPLAY HELPER ===
# =============================================================================

@st.cache_data(show_spinner=False)
def _pdf_to_png_bytes(file_path: str, dpi: int = 150) -> bytes | None:
    """Render first page of a PDF to PNG bytes using PyMuPDF (cached per path)."""
    try:
        import fitz  # pymupdf
        doc = fitz.open(file_path)
        mat = fitz.Matrix(dpi / 72, dpi / 72)
        pix = doc[0].get_pixmap(matrix=mat, alpha=False)
        return pix.tobytes("png")
    except Exception:
        return None


def show_pdf(file_path: str, height: int = 850):
    if not os.path.exists(file_path):
        st.warning(f"Figure file not found: {os.path.basename(file_path)}")
        return
    img_bytes = _pdf_to_png_bytes(file_path)
    if img_bytes:
        st.image(img_bytes, use_column_width=True)
    else:
        st.info(f"Could not render preview for {os.path.basename(file_path)}. "
                "Use the download button above.")


# =============================================================================
# === APP_TREES FUNCTIONS ===
# =============================================================================

TREE_FAMILY_MAP = {
    "All families": ("Figure_S1.pdf", None, None, None),
    "CSP (10 genes)":  ("Figure_S2.pdf", "Additional_file_26_CSP_orthology.csv", "CSP_tree_clean.txt", "CSP_sequences_clean.fasta"),
    "GR (39 genes)":   ("Figure_S3.pdf", "Additional_file_21_GR_orthology.csv",  "GR_tree_clean.txt",  "GR_sequences_clean.fasta"),
    "IR (101 genes)":  ("Figure_S4.pdf", "Additional_file_22_IR_orthology.csv",  "IR_tree_clean.txt",  "IR_sequences_clean.fasta"),
    "OBP (75 genes)":  ("Figure_S5.pdf", "Additional_file_23_OBP_orthology.csv", "OBP_tree_clean.txt", "OBP_sequences_clean.fasta"),
    "OR (192 genes)":  ("Figure_S6.pdf", "Additional_file_20_OR_orthology.csv",  "OR_tree_clean.txt",  "OR_sequences_clean.fasta"),
    "PPK (39 genes)":  ("Figure_S7.pdf", "Additional_file_24_PPK_orthology.csv", "PPK_tree_clean.txt", "PPK_sequences_clean.fasta"),
    "TRP (111 genes)": ("Figure_S8.pdf", "Additional_file_25_TRP_orthology.csv", "TRP_tree_clean.txt", "TRP_sequences_clean.fasta"),
}

ORTHOLOGY_TYPE_COLORS = {
    "1:1_ortholog": "#2E7D32",
    "co-ortholog": "#1565C0",
    "ortholog_group": "#E65100",
    "H.illucens_specific": "#7B1FA2",
}

TREE_KEY_FINDINGS = {
    "CSP (10 genes)": "CSP10 (→ *DmelCSP4*) is a cysteine string protein by naming convergence, not a canonical chemosensory protein.",
    "GR (39 genes)": "GR1–GR5 form a CO₂ receptor clade (*Gr21a/Gr63a* orthologs, BS = 87–98). GR40 → *Gr43a* fructose receptor (BS = 99).",
    "IR (101 genes)": "IR7 = *Ir8a* ortholog (BS = 86); IR3 = *Ir25a* ortholog (BS = 95). 51 antennal tuning receptor paralogues.",
    "OBP (75 genes)": "Diverse family; includes Obp31 (mating-induced → mated-female-biased reversal) and Obp27/30/48 (antenna-biased).",
    "OR (192 genes)": "Largest family. ORco AND Or_putative_Orco-like both in Orco clade (BS = 98); provisional annotation retained.",
    "PPK (39 genes)": "PPK3 = *ppk23* ortholog (BS = 100). PPK27/28 in *ppk15* neighbourhood (NOT *ppk25* as initially annotated). PPK18 = Nanchung co-ortholog (BS = 100).",
    "TRP (111 genes)": "TRPN expansion: 44 genes (NompC clade) — larger than TRPm (15 genes). Mechanosensory roles likely.",
}


def render_plotly_tree(newick_str: str, highlight_gene: str = "", height: int = 700):
    """Interactive rectangular phylogram via Plotly. No CDN — pure Python/Plotly/BioPython."""
    import re as _re
    import sys as _sys
    from Bio import Phylo as _Phylo
    from io import StringIO as _StringIO
    import plotly.graph_objects as go

    _sys.setrecursionlimit(5000)

    # IQ-TREE format: ):length[bootstrap] → )bootstrap:length
    nwk = _re.sub(r'\):([0-9]+(?:\.[0-9]+)?(?:[eE][+\-]?[0-9]+)?)\[([0-9]+)\]',
                  r')\2:\1', newick_str)
    try:
        tree = _Phylo.read(_StringIO(nwk), "newick")
    except Exception as exc:
        st.error(f"Tree parsing failed: {exc}")
        return

    # ── Layout: x = cumulative branch length; y = leaf rank, internals = mean ──
    def assign_x(clade, x=0.0):
        clade._x = x
        for child in clade.clades:
            assign_x(child, x + (child.branch_length or 0.0))

    _ctr = [0]
    def assign_y(clade):
        if clade.is_terminal():
            clade._y = _ctr[0]
            _ctr[0] += 1
        else:
            for child in clade.clades:
                assign_y(child)
            clade._y = sum(c._y for c in clade.clades) / len(clade.clades)

    assign_x(tree.root)
    assign_y(tree.root)
    n_leaves = _ctr[0]

    # ── Branch lines ──
    xs, ys = [], []
    def add_lines(clade, parent=None):
        if parent is not None:
            xs.extend([parent._x, clade._x, None])
            ys.extend([clade._y, clade._y, None])
        if not clade.is_terminal():
            cy = [c._y for c in clade.clades]
            xs.extend([clade._x, clade._x, None])
            ys.extend([min(cy), max(cy), None])
            for child in clade.clades:
                add_lines(child, clade)

    add_lines(tree.root)

    # ── Leaf nodes ──
    def get_sp(name):
        for sp in SPECIES_COLORS:
            if name and name.startswith(sp):
                return sp
        return None

    hl = (highlight_gene or "").strip().lower()
    lx, ly, lcolor, lsize, lsymbol, lhover = [], [], [], [], [], []
    hl_clades = []

    for clade in tree.get_terminals():
        nm = clade.name or ""
        sp = get_sp(nm)
        is_hl = bool(hl) and hl in nm.lower()
        lx.append(clade._x)
        ly.append(clade._y)
        lhover.append(nm)
        if is_hl:
            lcolor.append("#FFD700")
            lsize.append(14)
            lsymbol.append("star")
            hl_clades.append(clade)
        else:
            lcolor.append(SPECIES_COLORS.get(sp, "#888888"))
            lsize.append(9 if sp == "Hill" else 7)
            lsymbol.append("circle")

    show_labels = n_leaves <= 60

    # ── Bootstrap labels (≥70 only) ──
    bx, by, bt = [], [], []
    for clade in tree.get_nonterminals():
        if clade.confidence is not None and int(clade.confidence) >= 70:
            bx.append(clade._x)
            by.append(clade._y)
            bt.append(str(int(clade.confidence)))

    # ── Build figure ──
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs, y=ys,
        mode="lines",
        line=dict(color="#444444", width=0.8),
        hoverinfo="none",
        showlegend=False,
    ))

    leaf_labels = [c.name or "" for c in tree.get_terminals()] if show_labels else [""] * n_leaves
    fig.add_trace(go.Scatter(
        x=lx, y=ly,
        mode=("markers+text" if show_labels else "markers"),
        text=leaf_labels,
        textposition="middle right",
        textfont=dict(size=8, family="monospace"),
        marker=dict(color=lcolor, size=lsize, symbol=lsymbol),
        customdata=lhover,
        hovertemplate="%{customdata}<extra></extra>",
        showlegend=False,
    ))

    if hl_clades:
        fig.add_trace(go.Scatter(
            x=[c._x for c in hl_clades],
            y=[c._y for c in hl_clades],
            mode="markers+text",
            text=[c.name for c in hl_clades],
            textposition="middle right",
            textfont=dict(size=10, color="#B8860B", family="monospace"),
            marker=dict(color="#FFD700", size=14, symbol="star",
                        line=dict(color="#B8860B", width=1.5)),
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
        ))

    if bt:
        fig.add_trace(go.Scatter(
            x=bx, y=by,
            mode="text",
            text=bt,
            textfont=dict(size=7, color="#999999"),
            hoverinfo="none",
            showlegend=False,
        ))

    max_x = max((c._x for c in tree.get_terminals()), default=1.0)
    fig.update_layout(
        height=height,
        margin=dict(l=10, r=(200 if show_labels else 80), t=15, b=40),
        plot_bgcolor="white",
        paper_bgcolor="white",
        xaxis=dict(
            title="Branch length",
            range=[-max_x * 0.02, max_x * (1.35 if show_labels else 1.08)],
            showgrid=True,
            gridcolor="#F0F0F0",
            zeroline=False,
        ),
        yaxis=dict(visible=False, range=[-1, n_leaves]),
        hoverlabel=dict(bgcolor="white", font_size=11, font_family="monospace"),
    )

    st.plotly_chart(fig, use_container_width=True)

    if hl_clades:
        names = ", ".join(c.name for c in hl_clades[:5])
        extra = f" (+{len(hl_clades) - 5} more)" if len(hl_clades) > 5 else ""
        st.success(f"Highlighted: **{names}**{extra}")
    elif hl:
        st.warning(f"No gene matching '{highlight_gene}' found in this tree.")


# render_phylocanvas_tree removed — replaced by render_plotly_tree above


def render_trees_tab():
    SUPP_FIG_DIR  = str(pathlib.Path(__file__).resolve().parent / "data" / "figures" / "supplementary")
    ORTH_DIR      = str(pathlib.Path(__file__).resolve().parent / "data" / "orthology")
    TREE_DIR      = str(pathlib.Path(__file__).resolve().parent / "data" / "trees")
    SEQ_TABLE_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "sequence_tables")
    FASTA_DIR     = str(pathlib.Path(__file__).resolve().parent / "data" / "sequences")

    st.header("Phylogenetic Trees")
    st.caption("Manuscript: **Figure 1** (main text, circular all families) | **Fig. S1** (full circular PDF) · **Figs. S2–S8** (per-family trees)")
    st.markdown(
        "Maximum-likelihood trees of *H. illucens* chemosensory gene families inferred with IQ-TREE "
        "(ultrafast bootstrap, 1000 replicates). Outgroups: *Drosophila melanogaster*, *Musca domestica*, "
        "*Calliphora stygia*, *Aedes aegypti* (GR), *Bombyx mori* (OBP)."
    )

    fam_choice = st.selectbox("Select gene family", list(TREE_FAMILY_MAP.keys()), key="trees_fam_sel")
    pdf_name, orth_name, tree_name, fasta_name = TREE_FAMILY_MAP[fam_choice]
    pdf_path  = os.path.join(SUPP_FIG_DIR, pdf_name)
    tree_path = os.path.join(TREE_DIR, tree_name) if tree_name else None
    fasta_path = os.path.join(FASTA_DIR, fasta_name) if fasta_name else None
    fam_key = fam_choice.split(" ")[0]  # e.g., "OR"

    if fam_choice == "All families":
        hill_total = sum(v.get("Hill", 0) for v in FAMILY_SPECIES_COUNTS.values())
        all_total  = sum(n for v in FAMILY_SPECIES_COUNTS.values() for n in v.values())
        st.markdown(
            f"Supplementary Fig. S1 — circular ML tree of all seven *H. illucens* chemosensory gene families. "
            f"**{hill_total} *H. illucens* genes** (546 retained in trees after deduplication; 567 annotated) "
            f"+ {all_total - hill_total} outgroup sequences = **{all_total} sequences total**."
        )

        # PDF + download
        tab_pdf_all, tab_interactive_all = st.tabs(["📄 Figure S1 (PDF)", "🌿 Browse by family"])
        with tab_pdf_all:
            if os.path.exists(pdf_path):
                with open(pdf_path, "rb") as f:
                    st.download_button("Download Fig. S1 PDF", data=f.read(),
                                       file_name="Figure_S1_all_families.pdf", mime="application/pdf",
                                       key="trees_dl_all_pdf")
            # Figure S1 is 4.4 MB — render via PyMuPDF to avoid iframe issues
            show_pdf(pdf_path, height=900)

        with tab_interactive_all:
            st.markdown(
                "The combined Newick for all families has corrupt gene labels (IQ-TREE intermediate format). "
                "Browse each family's interactive tree individually below — same trees as in Fig. S1."
            )
            st.subheader("Gene counts per species across all families")
            rows = []
            for fam, sp_counts in FAMILY_SPECIES_COUNTS.items():
                for sp, n in sp_counts.items():
                    rows.append({"Family": fam, "Species_code": sp,
                                 "Species": SPECIES_FULL_NAMES.get(sp, sp), "N_genes": n})
            df_all = pd.DataFrame(rows)
            pivot = df_all.pivot_table(index="Species", columns="Family", values="N_genes",
                                       aggfunc="sum", fill_value=0)
            pivot["Total"] = pivot.sum(axis=1)
            pivot = pivot.sort_values("Total", ascending=False)
            st.dataframe(pivot, use_container_width=True)

            st.markdown("---")
            st.subheader("Interactive trees — select family")
            all_fam_keys = [k for k in TREE_FAMILY_MAP if k != "All families"]
            sel_fam = st.selectbox("Family to view", all_fam_keys, key="trees_all_fam_sel")
            _, _, sel_tree_name, sel_fasta_name = TREE_FAMILY_MAP[sel_fam]
            sel_fam_key = sel_fam.split(" ")[0]
            sel_tree_path  = os.path.join(TREE_DIR, sel_tree_name) if sel_tree_name else None
            sel_fasta_path = os.path.join(FASTA_DIR, sel_fasta_name) if sel_fasta_name else None
            sel_seq_path   = os.path.join(SEQ_TABLE_DIR, f"{sel_fam_key}_sequence_table.csv")
            sel_seq_df     = pd.read_csv(sel_seq_path) if os.path.exists(sel_seq_path) else None

            search_all = st.text_input("Search gene", placeholder="e.g. OR5, HillGR14, Q9W5G6",
                                       key="trees_all_search")
            if sel_fasta_path and os.path.exists(sel_fasta_path):
                with open(sel_fasta_path, "rb") as f:
                    fasta_b = f.read()
                st.download_button(f"Download {sel_fam_key} sequences (FASTA)", data=fasta_b,
                                   file_name=sel_fasta_name, mime="text/plain",
                                   key="trees_all_fasta_dl")

            sp_c = FAMILY_SPECIES_COUNTS.get(sel_fam_key, {})
            if sp_c:
                sp_cols = st.columns(len(sp_c))
                for ci, (sp, n) in enumerate(sp_c.items()):
                    with sp_cols[ci]:
                        col = SPECIES_COLORS.get(sp, "#888")
                        st.markdown(
                            f"<div style='text-align:center'>"
                            f"<span style='color:{col};font-size:1.3em;font-weight:bold'>{n}</span><br>"
                            f"<i style='font-size:0.8em'>{SPECIES_FULL_NAMES.get(sp, sp)}</i></div>",
                            unsafe_allow_html=True,
                        )

            if sel_tree_path and os.path.exists(sel_tree_path):
                sel_nwk = pathlib.Path(sel_tree_path).read_text(encoding="utf-8").strip()
                st.caption("🖱 Drag to pan · Scroll to zoom · Hover leaf for gene name · Search to highlight")
                render_plotly_tree(sel_nwk, highlight_gene=search_all, height=650)
                leg_cols = st.columns(len(SPECIES_COLORS))
                for ci, (sp, col) in enumerate(SPECIES_COLORS.items()):
                    with leg_cols[ci]:
                        st.markdown(f"<span style='color:{col}'>●</span> <b>{sp}</b>",
                                    unsafe_allow_html=True)
        return

    # --- Per-family view ---
    finding = TREE_KEY_FINDINGS.get(fam_choice, "")
    if finding:
        st.info(f"**Key finding:** {finding}")

    # Species counts for this family
    sp_counts = FAMILY_SPECIES_COUNTS.get(fam_key, {})
    if sp_counts:
        st.markdown("**Genes in tree per species:**")
        sp_cols = st.columns(len(sp_counts))
        for ci, (sp, n) in enumerate(sp_counts.items()):
            with sp_cols[ci]:
                col = SPECIES_COLORS.get(sp, "#888")
                st.markdown(
                    f"<div style='text-align:center'>"
                    f"<span style='color:{col};font-size:1.4em;font-weight:bold'>{n}</span><br>"
                    f"<span style='font-size:0.8em'><i>{SPECIES_FULL_NAMES.get(sp, sp)}</i></span></div>",
                    unsafe_allow_html=True,
                )

    st.markdown("---")

    fig_num_label = pdf_name.replace("Figure_S", "S").replace(".pdf", "")  # e.g. "S2"
    tab_pdf, tab_interactive = st.tabs([f"📄 Fig. {fig_num_label} (PDF)", "🌿 Interactive Tree"])

    # ── PDF tab ──────────────────────────────────────────────────────────────
    with tab_pdf:
        col_tree, col_orth = st.columns([3, 2])
        with col_tree:
            if os.path.exists(pdf_path):
                with open(pdf_path, "rb") as f:
                    st.download_button("Download PDF", data=f.read(),
                                       file_name=pdf_name, mime="application/pdf",
                                       key=f"trees_dl_{fam_choice}_pdf")
            show_pdf(pdf_path, height=850)

        with col_orth:
            st.subheader("Orthology table")
            if orth_name:
                orth_path = os.path.join(ORTH_DIR, orth_name)
                if os.path.exists(orth_path):
                    orth_df = load_csv(orth_path)
                    if "Orthology_type" in orth_df.columns:
                        type_counts = orth_df["Orthology_type"].value_counts()
                        st.caption("Orthology type breakdown:")
                        for ot, n in type_counts.items():
                            color = ORTHOLOGY_TYPE_COLORS.get(ot, "#444")
                            st.markdown(
                                f"<span style='color:{color}'>●</span> **{ot}**: {n}",
                                unsafe_allow_html=True,
                            )
                        st.markdown("---")
                    search_q = st.text_input("Filter genes", key=f"trees_search_{fam_choice}",
                                             placeholder="e.g. OR1, Gr43a, ppk23, Drosophila")
                    if search_q.strip():
                        sq = search_q.strip().lower()
                        mask = orth_df.apply(
                            lambda col: col.astype(str).str.lower().str.contains(sq, na=False), axis=1
                        ).any(axis=1)
                        orth_show = orth_df[mask]
                    else:
                        orth_show = orth_df
                    display_cols = ["Hill_gene_symbol", "Orthology_type", "Putative_function_from_Dmel",
                                    "Bootstrap_support", "Nearest_outgroup_species"]
                    display_cols = [c for c in display_cols if c in orth_show.columns]
                    st.caption(f"{len(orth_show)} of {len(orth_df)} genes shown")
                    st.dataframe(orth_show[display_cols] if display_cols else orth_show,
                                 use_container_width=True, height=580)
                    with open(orth_path, "rb") as f:
                        st.download_button("Download orthology CSV", data=f.read(),
                                           file_name=orth_name, mime="text/csv",
                                           key=f"trees_orth_dl_{fam_choice}")
                else:
                    st.info(f"Orthology file not found: {orth_name}")

    # ── Interactive tree tab ──────────────────────────────────────────────────
    with tab_interactive:
        # Load sequence table for color mapping
        seq_tbl_path = os.path.join(SEQ_TABLE_DIR, f"{fam_key}_sequence_table.csv")
        seq_table_df = None
        if os.path.exists(seq_tbl_path):
            seq_table_df = pd.read_csv(seq_tbl_path)

        # Search / highlight
        search_gene = st.text_input(
            "Search gene (highlights in tree and filters table)",
            placeholder="e.g. OR5, Gr43a, HillGR14, Q9W5G6",
            key=f"trees_interactive_search_{fam_choice}",
        )

        # FASTA download
        if fasta_path and os.path.exists(fasta_path):
            with open(fasta_path, "rb") as f:
                fasta_bytes = f.read()
            n_seqs = fasta_bytes.count(b">")
            st.download_button(
                f"Download sequences ({n_seqs} genes, FASTA)",
                data=fasta_bytes,
                file_name=fasta_name,
                mime="text/plain",
                key=f"trees_fasta_dl_{fam_choice}",
            )

        # Interactive tree
        if tree_path and os.path.exists(tree_path):
            newick_str = pathlib.Path(tree_path).read_text(encoding="utf-8").strip()
            st.caption("🖱 Drag to pan · Scroll to zoom · Hover leaf for gene name · Search to highlight")
            render_plotly_tree(newick_str, highlight_gene=search_gene.strip(), height=680)
            # Species color legend
            st.markdown("**Species color key:**")
            leg_cols = st.columns(len(SPECIES_COLORS))
            for ci, (sp, col) in enumerate(SPECIES_COLORS.items()):
                with leg_cols[ci]:
                    full = SPECIES_FULL_NAMES.get(sp, sp)
                    st.markdown(
                        f"<span style='color:{col};font-size:1.1em'>●</span> "
                        f"<b>{sp}</b> <i style='font-size:0.8em'>{full}</i>",
                        unsafe_allow_html=True,
                    )
        else:
            st.info("Interactive tree file not found.")

        # Sequence table
        if seq_table_df is not None and not seq_table_df.empty:
            st.markdown("---")
            st.subheader("Gene metadata table")
            st.caption("Includes protein accession, NCBI description, and species for every gene in the tree.")
            show_cols = [c for c in ["Tree_ID", "Gene_symbol", "Species", "Protein_accession",
                                     "NCBI_description", "Database_source"] if c in seq_table_df.columns]
            if search_gene.strip():
                sq = search_gene.strip().lower()
                mask = seq_table_df.apply(
                    lambda col: col.astype(str).str.lower().str.contains(sq, na=False), axis=1
                ).any(axis=1)
                tbl_show = seq_table_df[mask]
            else:
                tbl_show = seq_table_df
            st.caption(f"{len(tbl_show)} of {len(seq_table_df)} genes shown")
            st.dataframe(tbl_show[show_cols] if show_cols else tbl_show,
                         use_container_width=True, height=400)


# =============================================================================
# === APP_FIGURES FUNCTIONS ===
# =============================================================================

MAIN_FIGURE_CAPTIONS = {
    "Figure 1 — Phylogeny of chemosensory families": (
        "**Figure 1.** Phylogenetic analysis of H. illucens chemosensory gene families. A) Adult H. "
        "illucens highlighting the three chemosensory appendages analyzed in this study: antennae, "
        "maxillary palps, and tarsi (image adapted from Beta Bugs Ltd., UK). B) Simplified phylogeny of "
        "representative dipteran families included in the comparative analyses, illustrating the basal "
        "phylogenetic position of H. illucens (Stratiomyidae) within the Brachycera. The time scale (Myr) "
        "indicates the approximate divergence of major dipteran lineages based on Wiegmann et al. (2011) "
        "[69]. C) Unrooted maximum-likelihood phylogeny of 1,620 chemosensory proteins, comprising 567 H. "
        "illucens sequences and 1,053 homologs from Aedes aegypti, Drosophila melanogaster, Musca "
        "domestica, Calliphora stygia, Bactrocera dorsalis, and Glossina morsitans. Proteins are grouped "
        "into the seven conserved chemosensory gene families: odorant receptors (OR), gustatory receptors "
        "(GR), ionotropic receptors (IR), odorant-binding proteins (OBP), chemosensory proteins (CSP), "
        "pickpocket channels (PPK), and transient receptor potential channels (TRP). Branches are "
        "color-coded by gene family. Expanded phylogenies with protein names, orthology assignments, and "
        "bootstrap support values are provided in Figures S1-S8. "
    ),
    "Figure 2 — Transcriptome overview": (
        "**Figure 2.** Appendage identity dominates transcriptomic organization in H. illucens. A) "
        "Principal component analysis (PCA) of global transcriptomic profiles from antennae, maxillary "
        "palps, and tarsi of virgin males (VM), virgin females (VF), and mated females (MF) (n = 3 "
        "biological replicates per group). Samples cluster primarily by appendage. Grey polygons "
        "delineate appendage clusters, numbers within polygons indicate cluster area in PC1–PC2 space, "
        "and grey lines connect appendage centroids. Numbers along the connecting lines indicate "
        "Euclidean distances between appendage centroids in PCA space. B) Pearson correlation matrix of "
        "mean global gene expression profiles across appendages and reproductive groups. Circle size and "
        "color intensity are proportional to the Pearson correlation coefficient. C) Variation "
        "partitioning of the global transcriptome showing the proportion of multivariate variance "
        "explained by appendage identity, sex, mating status, and unexplained variance. D) PCA of the "
        "global transcriptome within each appendage, illustrating variation among virgin males, virgin "
        "females, and mated females. E) PCA of the 567 annotated chemosensory transcript models showing "
        "stronger appendage-specific separation than observed for the global transcriptome. Grey "
        "polygons, centroid connections, cluster areas, and centroid distances are displayed as in panel "
        "A. F) Pearson correlation matrix of chemosensory gene expression profiles. G) Variation "
        "partitioning of chemosensory gene expression. H) Within-appendage PCA of chemosensory gene "
        "expression across reproductive groups. Together, these analyses show that appendage identity is "
        "the primary determinant of both global and chemosensory transcriptomic variation, whereas sex "
        "and mating status contribute comparatively modest effects. "
    ),
    "Figure 3 — Appendage-specific divergence": (
        "**Figure 3.** Appendage-specific transcriptional divergence and chemosensory specialization in "
        "H. illucens. A) Pairwise differential expression analyses comparing antennae, maxillary palps, "
        "and tarsi (|log₂FC| ≥ 1; adjusted p < 0.001). Volcano plots show significantly up- and "
        "downregulated genes for each comparison. Chemosensory genes are highlighted in red. Numbers "
        "above each plot indicate the number of chemosensory differentially expressed genes (DEGs) "
        "relative to the total number of DEGs, and the chemosensory genes with the largest fold changes "
        "are labelled. B) Chemosensory and non-chemosensory composition of the DEGs identified in each "
        "pairwise comparison. Bar heights represent DEG counts and the values within the bars the "
        "percentage of the DEG set; genes belonging to the seven chemosensory families are counted as "
        "chemosensory irrespective of their GO domain. C) Identification of consensus appendage-biased "
        "gene sets by intersecting pairwise DEG comparisons. Venn diagrams show genes consistently "
        "enriched in the antenna, maxillary palp, or tarsi relative to the other two appendages, "
        "providing appendage-specific transcriptional signatures. D) Classification of the 510 annotated "
        "chemosensory genes into appendage-specific, appendage-biased, broadly expressed, and "
        "non-expressed categories across the seven chemosensory gene families (Or, Gr, Ir, Obp, Csp, ppk, "
        "and Trp). Classification criteria are described in the Materials and Methods. Together, these "
        "analyses demonstrate that each appendage possesses a distinct transcriptional identity supported "
        "by characteristic enrichment of chemosensory gene families, with Ors predominating in the "
        "antennae, Obps in the maxillary palps, and ppk, Ir, Gr, and Trp genes in the tarsi. "
    ),
    "Figure 4 — Sex & mating remodelling": (
        "**Figure 4.** Sex- and mating-dependent transcriptional remodelling across adult chemosensory "
        "appendages. Genes are called differentially expressed at |log₂FC| ≥ 1 and adjusted p < 0.001. "
        "Mating-responsive denotes significance in the same direction against both virgin females and "
        "virgin males. A) Chemosensory loci responding to sex (left) or mating (right), one row per gene "
        "locus, grouped by family. Stems run from no change to the log₂ fold change; dot area is mean "
        "normalised expression in that appendage. Colour denotes appendage: antenna, dark blue; maxillary "
        "palp, mid blue; tarsi, grey. Panels A and B share one dot-size scale. B) Non-chemosensory genes, "
        "the top 20 within each contrast ranked by |log₂FC| × −log₁₀(adjusted p), sex above and mating "
        "below. Colours as in A. Names are abbreviated; an asterisk marks a name inferred from the best "
        "Drosophila melanogaster BLASTp match rather than a H. illucens annotation. C) Venn diagrams of "
        "genes shared among appendages for four programmes: male-biased and female-biased in virgins, and "
        "reduced or induced in mated females. Representative genes are named beside the regions they "
        "occupy. D) Reorganisation of sex bias by mating. Sex bias measured in virgins (virgin female "
        "versus virgin male, x) against sex bias measured in mated females (mated female versus virgin "
        "male, y), for every gene significant in either contrast. Because both axes share virgin males as "
        "the reference, a change of sign between them is a genuine reversal of sex bias. Genes are "
        "classified by their behaviour across the two reproductive states: sex-biased in both states "
        "(blue), sex-biased in one state only (grey), or sex-biased in both with the direction reversed "
        "(ringed). Red marks chemosensory genes. Quadrant labels give the direction of bias in each "
        "group, and the counts give the genes in each quadrant. Axes are clipped at ±10 log2FC, with any "
        "gene beyond drawn on the boundary. "
    ),
    "Figure 5 — Hierarchical model": (
        "**Figure 5.** A hierarchical model for the evolution and modulation of peripheral chemosensory "
        "systems in Diptera. A conceptual model summarizing the evolutionary and physiological "
        "organization of the Hermetia illucens chemosensory system inferred from this study. The model "
        "proposes three hierarchical levels operating across distinct evolutionary timescales. (1) "
        "Evolutionary conservation (deep evolutionary timescale; proposed, not tested here): the antenna, "
        "maxillary palp and tarsi carry distinct anatomical and functional identities, which we "
        "hypothesise became stabilised early in dipteran evolution. The present data establish that these "
        "identities are invariant across sex and reproductive state in H. illucens; dating their origin "
        "requires appendage-resolved data from further lineages. (2) Molecular diversification "
        "(intermediate evolutionary timescale): within these conserved appendages, chemosensory receptor "
        "gene families (OR, GR, IR, OBP, CSP, TRP, and PPK) diversified through lineage-specific "
        "expansion, contraction, and differential deployment, enabling ecological specialization without "
        "fundamentally altering appendage identity. (3) Physiological modulation (adult timescale): sex, "
        "reproductive state, and circadian context modulate the activity of conserved sensory circuits "
        "through regulation of shared signaling pathways, including sex determination (transformer), "
        "circadian regulation (daywake), neuropeptide signalling (the SIFamide receptor) and remodelling "
        "of the cuticle and extracellular matrix, rather than through wholesale remodeling of peripheral "
        "receptor repertoires. In the tarsi this modulation abolishes constitutive sex-biased expression "
        "and installs a female-specific programme of induced genes, with no accompanying change in "
        "receptor-gene expression. Together, these findings support a hierarchical model in which stable "
        "appendage-specific sensory architectures provide the anatomical setting in which receptor-family "
        "diversification and reversible physiological regulation generate behavioral flexibility. "
    ),
}

SUPP_FIGURE_CAPTIONS = {
    "Figure S1 — Tree: all families": (
        "**Figure S1.** phylogenetic reconstruction of the H. illucens chemosensory repertoire. Circular "
        "maximum-likelihood phylogeny of all 1,620 chemosensory proteins analyzed in this study, "
        "comprising 567 annotated H. illucens sequences and 1,053 homologous proteins from Aedes aegypti, "
        "Drosophila melanogaster, Musca domestica, Calliphora stygia, Bactrocera dorsalis, and Glossina "
        "morsitans. Branches are color-coded according to the seven conserved chemosensory protein "
        "families: odorant receptors (OR), gustatory receptors (GR), ionotropic receptors (IR), transient "
        "receptor potential channels (TRP), odorant-binding proteins (OBP), chemosensory proteins (CSP), "
        "and pickpocket channels (PPK). Terminal labels identify individual proteins, with species "
        "indicated by label color. Grey circles at internal nodes indicate bootstrap support (70–100%), "
        "with node size proportional to support; values below 70 are omitted. This overview summarizes "
        "the phylogenetic relationships among all annotated chemosensory proteins and provides the basis "
        "for the detailed family-specific phylogenies presented in Figures S2–S8. "
    ),
    "Figure S2 — Tree: CSP": (
        "**Figure S2.** Maximum-likelihood phylogeny of the H. illucens chemosensory protein (CSP) "
        "family. Circular maximum-likelihood phylogeny of the 10 annotated H. illucens CSP proteins "
        "reconstructed together with homologous CSP sequences from Drosophila melanogaster, Musca "
        "domestica, Bactrocera dorsalis, and Glossina morsitans. Terminal labels identify individual "
        "proteins, with species indicated by label color. Grey circles at internal nodes indicate "
        "bootstrap support (50–100%), with node size proportional to support; values below 50 are "
        "omitted. The phylogeny identifies orthologous relationships and lineage-specific expansions "
        "within the CSP family, providing the basis for the annotation and nomenclature of H. illucens "
        "CSP genes. "
    ),
    "Figure S3 — Tree: GR": (
        "**Figure S3.** Maximum-likelihood phylogeny of the H. illucens gustatory receptor (GR) family. "
        "Circular maximum-likelihood phylogeny of the 39 annotated H. illucens gustatory receptor (GR) "
        "proteins reconstructed together with homologous GR sequences from Aedes aegypti and Drosophila "
        "melanogaster. Terminal labels identify individual receptors, with species indicated by label "
        "color. Grey circles at internal nodes indicate bootstrap support (50–100%), with node size "
        "proportional to support; values below 50 are omitted. The phylogeny identifies orthologous "
        "relationships between H. illucens and well-characterized dipteran GRs, including the conserved "
        "carbon dioxide receptor clade (HillGR1–HillGR5), the fructose receptor ortholog HillGR40 "
        "(Gr43a), and the bitter receptor orthologs HillGR24 (Gr66a) and HillGR39 (Gr33a), providing the "
        "basis for functional annotation of the H. illucens GR repertoire. "
    ),
    "Figure S4 — Tree: IR": (
        "**Figure S4.** Maximum-likelihood phylogeny of the H. illucens ionotropic receptor (IR) family. "
        "Circular maximum-likelihood phylogeny of the 101 annotated H. illucens ionotropic receptor (IR) "
        "proteins reconstructed together with homologous IR sequences from Aedes aegypti, Drosophila "
        "melanogaster, and Calliphora stygia. Terminal labels identify individual receptors, with species "
        "indicated by label color. Grey circles at internal nodes indicate bootstrap support (53–100%), "
        "with node size proportional to support; values below 53 are omitted. The phylogeny identifies "
        "conserved co-receptors (Ir8a, Ir25a and Ir76b), antennal tuning receptors, and lineage-specific "
        "expansions within the H. illucens IR repertoire, providing the basis for functional annotation "
        "and gene nomenclature. "
    ),
    "Figure S5 — Tree: OBP": (
        "**Figure S5.** Maximum-likelihood phylogeny of the H. illucens odorant-binding protein (OBP) "
        "family. Circular maximum-likelihood phylogeny of the 75 annotated H. illucens odorant-binding "
        "protein (OBP) sequences reconstructed together with homologous OBPs from Aedes aegypti, "
        "Drosophila melanogaster, and Musca domestica. Terminal labels identify individual proteins, with "
        "species indicated by label color. Grey circles at internal nodes indicate bootstrap support "
        "(50–100%), with node size proportional to support; values below 50 are omitted. The phylogeny "
        "identifies orthologous relationships and lineage-specific expansions within the H. illucens OBP "
        "family and provides the basis for gene annotation and nomenclature. "
    ),
    "Figure S6 — Tree: OR": (
        "**Figure S6.** Maximum-likelihood phylogeny of the H. illucens odorant receptor (OR) family. "
        "Circular maximum-likelihood phylogeny of the 192 annotated H. illucens odorant receptor (OR) "
        "proteins reconstructed together with homologous OR sequences from Aedes aegypti, Drosophila "
        "melanogaster, Musca domestica, and Calliphora stygia. Terminal labels identify individual "
        "receptors, with species indicated by label color. Grey circles at internal nodes indicate "
        "bootstrap support (50–100%), with node size proportional to support; values below 50 are "
        "omitted. The phylogeny identifies the conserved odorant receptor co-receptors (HillORco and "
        "HillORco2) together with extensive lineage-specific expansion of the OR repertoire, providing "
        "the basis for functional annotation and evolutionary comparisons. "
    ),
    "Figure S7 — Tree: PPK": (
        "**Figure S7.** Maximum-likelihood phylogeny of the H. illucens pickpocket (PPK) family. Circular "
        "maximum-likelihood phylogeny of the 39 annotated H. illucens pickpocket (PPK) proteins "
        "reconstructed together with homologous PPK sequences from Aedes aegypti and Drosophila "
        "melanogaster. Terminal labels identify individual proteins, with species indicated by label "
        "color. Grey circles at internal nodes indicate bootstrap support (50–100%), with node size "
        "proportional to support; values below 50 are omitted. The phylogeny identifies orthologous "
        "relationships within the PPK family, including the ppk23 pheromone-sensing lineage and the Nach "
        "lineage, providing the basis for functional annotation of the H. illucens PPK repertoire. "
    ),
    "Figure S8 — Tree: TRP": (
        "**Figure S8.** Maximum-likelihood phylogeny of the H. illucens transient receptor potential "
        "(TRP) family. Circular maximum-likelihood phylogeny of the 111 annotated H. illucens transient "
        "receptor potential (TRP) proteins reconstructed together with homologous TRP sequences from "
        "Aedes aegypti and Drosophila melanogaster. Terminal labels identify individual proteins, with "
        "species indicated by label color. Grey circles at internal nodes indicate bootstrap support "
        "(50–100%), with node size proportional to support; values below 50 are omitted. The phylogeny "
        "resolves the major TRP subfamilies, including TRPA, TRPM, TRPN (NompC), TRPC, TRPV, TRPML, PKD, "
        "and Brivido channels, highlighting lineage-specific expansions within the H. illucens TRP "
        "repertoire and providing the basis for functional annotation. "
    ),
    "Figure S9 — GO: antenna-biased": (
        "**Figure S9.** Gene Ontology enrichment of antenna-biased genes. with bar length giving the "
        "significance of over-representation, −log₁₀ of the Benjamini–Hochberg adjusted p, and gene "
        "counts printed beside each bar. Terms are ordered by significance; bars that do not reach q < "
        "0.05 are drawn pale. GO terms shown in bold and enclosed in a black box denote the biologically "
        "salient category discussed in the Results: olfactory receptor activity, the dominant molecular "
        "function of the antenna. "
    ),
    "Figure S10 — GO: palp-biased": (
        "**Figure S10.** Gene Ontology enrichment of maxillary palp-biased genes. Unlike the antenna, the "
        "maxillary palp exhibits a heterogeneous functional profile: the most strongly over-represented "
        "terms are solute:inorganic anion antiporter activity and odorant binding, with structural "
        "constituents of the ribosome and cuticle and monooxygenase activity also enriched.with bar "
        "length giving the significance of over-representation, −log₁₀ of the Benjamini–Hochberg adjusted "
        "p, and gene counts printed beside each bar. Terms are ordered by significance; bars that do not "
        "reach q < 0.05 are drawn pale. solute:inorganic anion antiporter activity and odorant binding, "
        "the two most strongly over-represented terms. "
    ),
    "Figure S11 — GO: tarsi-biased": (
        "**Figure S11.** Gene Ontology enrichment of tarsus-biased genes. The tarsus has the largest "
        "appendage-biased repertoire, but over-representation is confined to signalling and regulatory "
        "functions, principally G protein-coupled receptor activity, DNA-binding transcription factor "
        "activity, regulation of alternative mRNA splicing and neuropeptide hormone activity.with bar "
        "length giving the significance of over-representation, −log₁₀ of the Benjamini–Hochberg adjusted "
        "p, and gene counts printed beside each bar. Terms are ordered by significance; bars that do not "
        "reach q < 0.05 are drawn pale. G protein-coupled receptor activity. Protein binding and "
        "nucleotide binding are the numerically largest categories in this set but are not "
        "over-represented against the annotated background (q = 0.77 and q = 1.00), and are therefore not "
        "boxed. "
    ),
    "Figure S12 — GO domain composition": (
        "**Figure S12.** Gene Ontology domain composition of appendage-biased genes. For each pairwise "
        "appendage comparison, the differentially expressed genes biased to either appendage (adjusted p "
        "< 0.001, |log₂FC| ≥ 1) are broken down by Gene Ontology domain, with chemosensory genes shown as "
        "a separate category. This breakdown is descriptive: Additional file 11 assigns exactly one GO "
        "term to each transcript, so the domain split reflects how the genome was annotated rather than a "
        "tested property of the gene sets. Statistically tested over-representation is reported in "
        "Figures S9–S11, and the chemosensory enrichment tests, with odds ratios and 95% confidence "
        "intervals, in Additional file 38. "
    ),
    "Figure S13 — Heatmap: OR": (
        "**Figure S13.** Expression heatmap of the H. illucens odorant receptor (OR) family. Heatmap "
        "showing normalized expression of the 192 annotated odorant receptor (OR) genes across the "
        "antennae, maxillary palps, and tarsi of virgin males (VM), virgin females (VF), and mated "
        "females (MF). Expression values are displayed as log₂(counts + 1) and centered on a threshold of "
        "10 normalized counts (white = 10 counts, blue = lower expression, red = higher expression). "
        "Genes are ordered according to the appendage showing the highest mean expression (antenna → "
        "maxillary palp → tarsi), and subsequently by expression level within each appendage. The heatmap "
        "illustrates the pronounced enrichment of OR expression in the antennae, with a smaller subset of "
        "appendage-specific receptors expressed in the maxillary palps and relatively few enriched in the "
        "tarsi. "
    ),
    "Figure S14 — Heatmap: GR": (
        "**Figure S14.** Expression heatmap of the H. illucens gustatory receptor (GR) family. Heatmap "
        "showing normalized expression of the 39 annotated gustatory receptor (GR) genes across the "
        "antennae, maxillary palps, and tarsi of virgin males (VM), virgin females (VF), and mated "
        "females (MF). Expression values are displayed as log₂(counts + 1) and centered on a threshold of "
        "10 normalized counts (white = 10 counts). Genes are ordered by appendage of highest mean "
        "expression and then by expression level. The heatmap highlights distinct appendage-specific GR "
        "subsets, including receptors preferentially expressed in the maxillary palps and tarsi, "
        "consistent with their roles in CO₂ detection and contact chemosensation. "
    ),
    "Figure S15 — Heatmap: IR": (
        "**Figure S15.** Expression heatmap of the H. illucens ionotropic receptor (IR) family. Heatmap "
        "showing normalized expression of the 101 annotated ionotropic receptor (IR) genes across the "
        "antennae, maxillary palps, and tarsi of virgin males (VM), virgin females (VF), and mated "
        "females (MF). Expression values are displayed as log₂(counts + 1) and centered on a threshold of "
        "10 normalized counts (white = 10 counts). Genes are ordered by appendage of highest mean "
        "expression and then by expression level. The heatmap illustrates broad expression of the "
        "conserved IR co-receptors together with appendage-specific tuning receptor subsets enriched in "
        "the antennae and tarsi. "
    ),
    "Figure S16 — Heatmap: OBP": (
        "**Figure S16.** Expression heatmap of the H. illucens odorant-binding protein (OBP) family. "
        "Heatmap showing normalized expression of the 75 annotated odorant-binding protein (OBP) genes "
        "across the antennae, maxillary palps, and tarsi of virgin males (VM), virgin females (VF), and "
        "mated females (MF). Expression values are displayed as log₂(counts + 1) and centered on a "
        "threshold of 10 normalized counts (white = 10 counts). Genes are ordered by appendage of highest "
        "mean expression and then by expression level. Distinct antennal-, maxillary palp-, and "
        "tarsus-enriched OBP clusters are apparent, illustrating the extensive appendage-specific "
        "diversification of odorant-binding proteins. "
    ),
    "Figure S17 — Heatmap: PPK": (
        "**Figure S17.** Expression heatmap of the H. illucens pickpocket (PPK) family. Heatmap showing "
        "normalized expression of the 39 annotated pickpocket (PPK) genes across the antennae, maxillary "
        "palps, and tarsi of virgin males (VM), virgin females (VF), and mated females (MF). Expression "
        "values are displayed as log₂(counts + 1) and centered on a threshold of 10 normalized counts "
        "(white = 10 counts). Genes are ordered by appendage of highest mean expression and then by "
        "expression level. Most PPK genes exhibit preferential expression in the tarsi, consistent with "
        "their established roles in contact chemosensation and mechanosensory function. "
    ),
    "Figure S18 — Heatmap: CSP": (
        "**Figure S18.** Expression heatmap of the H. illucens chemosensory protein (CSP) family. Heatmap "
        "showing normalized expression of the 10 annotated chemosensory protein (CSP) genes across the "
        "antennae, maxillary palps, and tarsi of virgin males (VM), virgin females (VF), and mated "
        "females (MF). Expression values are displayed as log₂(counts + 1) and centered on a threshold of "
        "10 normalized counts (white = 10 counts). Genes are ordered by appendage of highest mean "
        "expression and then by expression level. In contrast to other chemosensory families, most CSP "
        "genes exhibit broad expression across all appendages, suggesting more generalized functions in "
        "the peripheral sensory system. "
    ),
    "Figure S19 — Heatmap: TRP": (
        "**Figure S19.** Expression heatmap of the H. illucens transient receptor potential (TRP) family. "
        "Heatmap showing normalized expression of the 111 annotated transient receptor potential (TRP) "
        "channel genes across the antennae, maxillary palps, and tarsi of virgin males (VM), virgin "
        "females (VF), and mated females (MF). Expression values are displayed as log₂(counts + 1) and "
        "centered on a threshold of 10 normalized counts (white = 10 counts). Genes are ordered by "
        "appendage of highest mean expression and then by expression level. The heatmap reveals "
        "widespread expression of TRP channels together with appendage-specific subsets, particularly "
        "within the tarsi, consistent with roles in thermo-, mechano-, and polymodal sensory signaling. "
    ),
    "Figure S20 — GO: sex-biased": (
        "**Figure S20.** Gene Ontology enrichment of sex-biased genes across chemosensory appendages. Top "
        "20 enriched Gene Ontology (GO) terms associated with sex-biased differentially expressed genes "
        "(virgin female vs. virgin male; adjusted Bar length gives the significance of "
        "over-representation, −log₁₀ of the Benjamini–Hochberg adjusted p; terms are ordered by "
        "significance and bars that do not reach q < 0.05 are drawn pale.p < 0.001, |log₂FC| ≥ 1) in the "
        "antennae (top row), maxillary palps (middle row), and tarsi (bottom row). GO terms are presented "
        "separately for Molecular Function (MF), Cellular Component (CC), and Biological Process (BP). "
        "Bar lengths indicate the number of genes assigned to each GO term. Panel titles indicate the "
        "appendage (rows: antennae, maxillary palps, tarsi) and Gene Ontology domain (columns) shown. GO "
        "terms shown in bold and enclosed in a black box denote the categories discussed in the main "
        "text: chemosensory (olfactory receptor activity, odorant binding), pheromone and detoxification "
        "metabolism (monooxygenase activity), cuticular structure (structural constituent of cuticle), "
        "and reproduction-, signalling- and immunity-related processes (lipid metabolic process, MAPK "
        "cascade, defense response, developmental process involved in reproduction). "
    ),
    "Figure S21 — GO: mating-responsive": (
        "**Figure S21.** Gene Ontology enrichment of mating-responsive genes across chemosensory "
        "appendages. Bar length gives the significance of over-representation, −log₁₀ of the "
        "Benjamini–Hochberg adjusted p; terms are ordered by significance and bars that do not reach q < "
        "0.05 are drawn pale.(mated female vs. virgin female and mated female vs. virgin male, "
        "significant in the same direction in both; adjusted p < 0.001, |log₂FC| ≥ 1) in the antennae "
        "(top row), maxillary palps (middle row), and tarsi (bottom row). GO terms are presented "
        "separately for Molecular Function (MF), Cellular Component (CC), and Biological Process (BP). "
        "Bar lengths indicate the number of genes assigned to each GO term. Panel titles indicate the "
        "appendage (rows: antennae, maxillary palps, tarsi) and Gene Ontology domain (columns) shown. GO "
        "terms shown in bold and enclosed in a black box denote the mating-responsive categories "
        "discussed in the main text: monooxygenase activity, odorant binding and G protein-coupled "
        "receptor activity, together with the biological processes MAPK cascade, defense response and "
        "lipid metabolic process. "
    ),
    "Figure S22 — GO: sex-biased only": (
        "**Figure S22.** Gene Ontology enrichment of constitutively sex-biased genes. Top enriched Gene "
        "Ontology (GO) terms associated with genes showing constitutive sex-biased expression (virgin "
        "female vs. virgin male) but no significant mating response within each chemosensory appendage "
        "(adjusted Bar length gives the significance of over-representation, −log₁₀ of the "
        "Benjamini–Hochberg adjusted p; terms are ordered by significance and bars that do not reach q < "
        "0.05 are drawn pale (p< 0.001, |log₂FC| ≥ 1). Results are shown separately for the antennae (top "
        "row), maxillary palps (middle row), and tarsi (bottom row). For each appendage, the top 15 "
        "enriched GO terms are presented for Molecular Function (MF), Cellular Component (CC), and "
        "Biological Process (BP), with bar lengths representing the number of genes assigned to each GO "
        "term. These analyses identify appendage-specific biological functions associated with "
        "constitutive sexual dimorphism). GO terms shown in bold and enclosed in a black box denote the "
        "categories underlying constitutive sexual dimorphism discussed in the main text: olfactory "
        "receptor activity and monooxygenase activity in the antenna, odorant binding in the palp and "
        "tarsi, and developmental process involved in reproduction in the tarsi. "
    ),
    "Figure S23 — GO: mating-responsive only": (
        "**Figure S23.** Gene Ontology enrichment of mating-responsive genes. Bar length gives the "
        "significance of over-representation, −log₁₀ of the Benjamini–Hochberg adjusted p; terms are "
        "ordered by significance and bars that do not reach q < 0.05 are drawn pale.genes responding to "
        "mating (significant in the same direction against both virgin females and virgin males) but not "
        "significantly sex GO terms shown in bold and enclosed in a black box denote the post-mating "
        "categories discussed in the main text: monooxygenase activity and odorant binding, and the "
        "biological processes immune system process and lipid metabolic process, most pronounced in the "
        "tarsi. "
    ),
    "Figure S24 — GO: sex & mating": (
        "**Figure S24.** Gene Ontology enrichment of genes regulated by both sex and mating status. Top "
        "enriched Gene Ontology (GO) terms associated with genes that were significantly differentially "
        "expressed in both the sex comparison (virgin female vs. virgin male) and the mating comparison "
        "(mated female vs. virgin female) within the same appendage (adjusted Bar length gives the "
        "significance of over-representation, −log₁₀ of the Benjamini–Hochberg adjusted p; terms are "
        "ordered by significance and bars that do not reach q < 0.05 are drawn pale.p < 0.001, |log₂FC| ≥ "
        "1). Results are shown separately for the antennae (top row), maxillary palps (middle row), and "
        "tarsi (bottom row). For each appendage, the top 15 enriched GO terms are presented for Molecular "
        "Function (MF), Cellular Component (CC), and Biological Process (BP), with bar lengths "
        "representing the number of genes assigned to each GO term. Absence of enriched terms is "
        "indicated where no GO category met the enrichment criteria. These analyses identify biological "
        "pathways jointly influenced by constitutive sexual dimorphism and reproductive state. GO terms "
        "shown in bold and enclosed in a black box denote the categories jointly regulated by sex and "
        "mating discussed in the main text: the tarsal terms monooxygenase activity, G protein-coupled "
        "receptor activity and odorant binding, and the biological processes MAPK cascade and defense "
        "response. "
    ),
    "Figure S25 — DE across appendages": (
        "**Figure S25.** Differential expression across appendages, by sex and by mating. Volcano plots "
        "for every gene tested in each appendage. Left column, sex bias in virgins (virgin female versus "
        "virgin male); right column, the mating response, which requires significance in the same "
        "direction against both virgin females and virgin males. Points are coloured grey where not "
        "significant, blue where significant, and red where significant and belonging to one of the seven "
        "chemosensory families. Dashed lines mark the significance thresholds (adjusted p < 0.001, "
        "|log₂FC| ≥ 1). Axes are clipped at |log₂FC| 12 and −log₁₀(adjusted p) 300; 38 of 148,968 points "
        "lie beyond and are drawn on the boundary. These are the distributions underlying the counts in "
        "Figure 4. "
    ),
    "Figure S26 — Mated-female-specific": (
        "**Figure S26.** Mated-female-specific expression against constitutive sex bias. Axes as in "
        "Figure 4D: sex bias measured in virgins (virgin female versus virgin male) against sex bias "
        "measured in mated females (mated female versus virgin male). Because both axes share virgin "
        "males as the reference, vertical displacement from the diagonal equals the log2 fold change "
        "between mated and virgin females; genes on the line show sex bias unchanged by mating. Coloured "
        "points are the mating-responsive set defined in the Methods, significant in the same direction "
        "against both virgin groups: orange where induced in mated females, blue where reduced, grey "
        "where not mating-responsive. Red rings mark chemosensory genes. Counts in each panel give the "
        "mating-responsive genes, and the chemosensory genes within each category, as a fraction of those "
        "plotted. Axes are clipped at ±10 log2FC. Because mated females were older than virgin males at "
        "collection, this contrast carries the same age confound as the mating comparison (Methods). "
    ),
}

_SUPP_FILES_ORDERED = [
    ("Additional_file_1_Appendage_DE_Antenna_vs_Tarsi.csv.gz",
     "Additional file 1 — Appendage DE: antenna vs tarsi (DESeq2) (0.8 MB) [gzip]", "application/gzip"),
    ("Additional_file_2_Appendage_DE_Antenna_vs_MaxillaryPalp.csv.gz",
     "Additional file 2 — Appendage DE: antenna vs maxillary palp (DESeq2) (0.8 MB) [gzip]", "application/gzip"),
    ("Additional_file_3_Appendage_DE_Tarsi_vs_MaxillaryPalp.csv.gz",
     "Additional file 3 — Appendage DE: tarsi vs maxillary palp (DESeq2) (0.8 MB) [gzip]", "application/gzip"),
    ("Additional_file_4_Antenna_sex_VirginFemale_vs_VirginMale.csv.gz",
     "Additional file 4 — Antenna, sex: virgin female vs virgin male (DESeq2) (1.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_5_Antenna_mating_MatedFemale_vs_VirginFemale.csv.gz",
     "Additional file 5 — Antenna, mating: mated female vs virgin female (DESeq2) (1.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_6_MaxillaryPalp_sex_VirginFemale_vs_VirginMale.csv.gz",
     "Additional file 6 — Maxillary palp, sex: virgin female vs virgin male (DESeq2) (1.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_7_MaxillaryPalp_mating_MatedFemale_vs_VirginFemale.csv.gz",
     "Additional file 7 — Maxillary palp, mating: mated female vs virgin female (DESeq2) (1.3 MB) [gzip]", "application/gzip"),
    ("Additional_file_8_Tarsi_sex_VirginFemale_vs_VirginMale.csv.gz",
     "Additional file 8 — Tarsi, sex: virgin female vs virgin male (DESeq2) (1.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_9_Tarsi_mating_MatedFemale_vs_VirginFemale.csv.gz",
     "Additional file 9 — Tarsi, mating: mated female vs virgin female (DESeq2) (1.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_10_normalized_counts_all_samples.csv.gz",
     "Additional file 10 — Normalised count matrix, all 27 libraries (2.8 MB) [gzip]", "application/gzip"),
    ("Additional_file_11_GO_annotations_all_genes.csv.gz",
     "Additional file 11 — GO annotations for all 24,828 genes (0.3 MB) [gzip]", "application/gzip"),
    ("Additional_file_12_AA_identity_matrix_CSP.csv",
     "Additional file 12 — Amino-acid identity matrix — CSP (1 kB)", "text/csv"),
    ("Additional_file_13_AA_identity_matrix_GR.csv",
     "Additional file 13 — Amino-acid identity matrix — GR (11 kB)", "text/csv"),
    ("Additional_file_14_AA_identity_matrix_IR.csv",
     "Additional file 14 — Amino-acid identity matrix — IR (71 kB)", "text/csv"),
    ("Additional_file_15_AA_identity_matrix_OBP.csv",
     "Additional file 15 — Amino-acid identity matrix — OBP (39 kB)", "text/csv"),
    ("Additional_file_16_AA_identity_matrix_OR.csv",
     "Additional file 16 — Amino-acid identity matrix — OR (0.2 MB)", "text/csv"),
    ("Additional_file_17_AA_identity_matrix_PPK.csv",
     "Additional file 17 — Amino-acid identity matrix — PPK (11 kB)", "text/csv"),
    ("Additional_file_18_AA_identity_matrix_TRP.csv",
     "Additional file 18 — Amino-acid identity matrix — TRP (82 kB)", "text/csv"),
    ("Additional_file_19_Diptera_chemosensory_comparison.csv",
     "Additional file 19 — Chemosensory family sizes across Diptera (0 kB)", "text/csv"),
    ("Additional_file_20_OR_orthology.csv",
     "Additional file 20 — Orthology assignments — OR (82 kB)", "text/csv"),
    ("Additional_file_21_GR_orthology.csv",
     "Additional file 21 — Orthology assignments — GR (15 kB)", "text/csv"),
    ("Additional_file_22_IR_orthology.csv",
     "Additional file 22 — Orthology assignments — IR (0.2 MB)", "text/csv"),
    ("Additional_file_23_OBP_orthology.csv",
     "Additional file 23 — Orthology assignments — OBP (34 kB)", "text/csv"),
    ("Additional_file_24_PPK_orthology.csv",
     "Additional file 24 — Orthology assignments — PPK (7 kB)", "text/csv"),
    ("Additional_file_25_TRP_orthology.csv",
     "Additional file 25 — Orthology assignments — TRP (18 kB)", "text/csv"),
    ("Additional_file_26_CSP_orthology.csv",
     "Additional file 26 — Orthology assignments — CSP (4 kB)", "text/csv"),
    ("Additional_file_27_Unknown_BLAST.csv",
     "Additional file 27 — BLASTp hits for unknown-domain genes (0.2 MB)", "text/csv"),
    ("Additional_file_28_AA_identity_1to1_orthologs.csv",
     "Additional file 28 — Amino-acid identity of 1:1 orthologs (4 kB)", "text/csv"),
    ("Additional_file_29_Antenna_matedsex_MatedFemale_vs_VirginMale.csv.gz",
     "Additional file 29 — Antenna, mated sex: mated female vs virgin male (DESeq2) (1.6 MB) [gzip]", "application/gzip"),
    ("Additional_file_30_MaxillaryPalp_matedsex_MatedFemale_vs_VirginMale.csv.gz",
     "Additional file 30 — Maxillary palp, mated sex: mated female vs virgin male (DESeq2) (1.7 MB) [gzip]", "application/gzip"),
    ("Additional_file_31_Tarsi_matedsex_MatedFemale_vs_VirginMale.csv.gz",
     "Additional file 31 — Tarsi, mated sex: mated female vs virgin male (DESeq2) (1.7 MB) [gzip]", "application/gzip"),
    ("Additional_file_32_multivariate_statistics.zip",
     "Additional file 32 — Multivariate statistics: PCA, PERMANOVA, variation partitioning (25 kB) [zip]", "application/zip"),
    ("Additional_file_33_GO_overrepresentation_significant.csv",
     "Additional file 33 — GO over-representation — significant terms (24 kB)", "text/csv"),
    ("Additional_file_34_GO_overrepresentation_all_tests.csv.gz",
     "Additional file 34 — GO over-representation — all tests (0.2 MB) [gzip]", "application/gzip"),
    ("Additional_file_35_identical_protein_clusters.csv",
     "Additional file 35 — Identical-protein clusters (17 kB)", "text/csv"),
    ("Additional_file_36_newick_trees.zip",
     "Additional file 36 — Newick trees for all seven families and the combined tree (48 kB) [zip]", "application/zip"),
    ("Additional_file_37_transcript_model_to_gene_locus.csv",
     "Additional file 37 — Transcript model to gene locus map (567 models, 393 loci) (21 kB)", "text/csv"),
    ("Additional_file_38_Figure3_statistical_tests.csv",
     "Additional file 38 — Statistical tests underlying Figure 3 (1 kB)", "text/csv"),
    ("Additional_file_39_FPKM_all_samples.csv.gz",
     "Additional file 39 — FPKM matrix for all 27 libraries (derived from Additional file 10) (1.7 MB) [gzip]", "application/gzip"),
]


def render_figures_tab():
    MAIN_FIG_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "figures" / "main")
    SUPP_FIG_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "figures" / "supplementary")
    SUPP_FILE_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "supplementary")

    st.header("Figures & Supplementary Data")

    fig_tab1, fig_tab2, fig_tab3, fig_tab4 = st.tabs(
        ["Main Figures", "Supplementary Figures",
         "Download Additional Files (AF1–AF39)", "Analysis code"]
    )

    with fig_tab1:
        fig_choice = st.radio(
            "Select figure",
            list(MAIN_FIGURE_CAPTIONS.keys()),
            horizontal=True,
            key="fig_main_choice",
        )
        label = fig_choice.split("—")[0].strip()
        fig_num = label.replace("Figure ", "")
        pdf_path = os.path.join(MAIN_FIG_DIR, f"Figure_{fig_num}.pdf")
        if os.path.exists(pdf_path):
            with open(pdf_path, "rb") as f:
                st.download_button(
                    f"Download {label} PDF",
                    data=f.read(),
                    file_name=f"Figure_{fig_num}.pdf",
                    mime="application/pdf",
                    key=f"fig_main_dl_{fig_num}",
                )
        show_pdf(pdf_path, height=900)
        st.markdown(MAIN_FIGURE_CAPTIONS[fig_choice])

    with fig_tab2:
        supp_choice = st.radio(
            "Select supplementary figure",
            list(SUPP_FIGURE_CAPTIONS.keys()),
            horizontal=True,
            key="fig_supp_choice",
        )
        s_label = supp_choice.split("—")[0].strip()
        s_num = s_label.replace("Figure S", "")
        pdf_path_s = os.path.join(SUPP_FIG_DIR, f"Figure_S{s_num}.pdf")
        if os.path.exists(pdf_path_s):
            with open(pdf_path_s, "rb") as f:
                st.download_button(
                    f"Download {s_label} PDF",
                    data=f.read(),
                    file_name=f"Figure_S{s_num}.pdf",
                    mime="application/pdf",
                    key=f"fig_supp_dl_{s_num}",
                )
        show_pdf(pdf_path_s, height=850)
        st.markdown(SUPP_FIGURE_CAPTIONS[supp_choice])

    with fig_tab3:
        st.subheader("Additional Files (AF1–AF39)")
        st.caption(
            "All additional files from Perets *et al.* 2026. "
            "AF1–AF9: DESeq2 results per contrast. AF10: count matrix. AF11: GO annotations. "
            "AF12–AF18: amino-acid identity matrices. AF19: Diptera comparison. "
            "AF20–AF26: per-family orthology tables. AF27: unknown-domain BLASTp. "
            "AF28: 1:1 ortholog identity. AF29–AF31: mated-female versus virgin-male contrasts. "
            "AF32: multivariate statistics. AF33–AF34: GO over-representation. "
            "AF35: identical-protein clusters. AF36: Newick trees. "
            "AF37: transcript-to-locus map. AF38: Figure 3 statistical tests. "
            "AF39: FPKM matrix. "
            "Large tables are served gzipped."
        )
        for fname, desc, mime in _SUPP_FILES_ORDERED:
            fpath = os.path.join(SUPP_FILE_DIR, fname)
            col_desc, col_btn = st.columns([5, 1])
            with col_desc:
                st.markdown(f"**{desc}**")
            with col_btn:
                if os.path.exists(fpath):
                    ext = fname.rsplit(".", 1)[-1].upper()
                    with open(fpath, "rb") as f:
                        st.download_button(
                            ext,
                            data=f.read(),
                            file_name=fname,
                            mime=mime,
                            key=f"supp_dl_{fname}",
                        )
                else:
                    st.caption("—")

    with fig_tab4:
        st.subheader("Analysis code")
        st.caption(
            "The R scripts that produce the figures from the Additional files. "
            "Under the current numbering these are no longer Additional files "
            "themselves, so they are published here and in the repository."
        )
        code_dir = str(pathlib.Path(__file__).resolve().parent / "data" / "code")
        notes = {
            "Figure_2.R": "Figure 2 and the multivariate statistics",
            "Figure_3.R": "Figure 3, Figures S9\u2013S19",
            "Figure_4.R": "Figure 4, Figures S20\u2013S26",
            "Figure_4A_lollipop.R": "Figure 4A, chemosensory lollipops",
            "Figure_4B_lollipop.R": "Figure 4B, non-chemosensory lollipops",
            "Figure_4_names.R": "shared gene-name shortening for panels B and D",
        }
        if not os.path.isdir(code_dir):
            st.info("Analysis code not bundled with this deployment.")
        else:
            for fn in sorted(os.listdir(code_dir)):
                fp = os.path.join(code_dir, fn)
                if not os.path.isfile(fp):
                    continue
                kb = os.path.getsize(fp) / 1e3
                col_d, col_b = st.columns([5, 1])
                with col_d:
                    st.markdown(f"**{fn}** — {notes.get(fn, 'analysis script')} ({kb:.0f} kB)")
                with col_b:
                    with open(fp, "rb") as fh:
                        st.download_button("R", fh.read(), file_name=fn,
                                           mime="text/plain", key=f"code_dl_{fn}")


# =============================================================================
# === MAIN APP ===
# =============================================================================


# =============================================================================
# === DATA EXPLORER ===
# The app exists so a reader of the manuscript can interrogate the data behind
# the figures, not just look at the figures again.  This section exposes the
# layers that have no figure of their own, or whose figure shows only the top
# of a much longer table: the GO over-representation tests, the per-locus
# chemosensory response, the BLAST evidence behind unnamed genes, and the
# comparative counts.  Everything here reads an Additional file directly.
# =============================================================================
SUPP_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "supplementary")

EXPLORER_FILES = {
    "go_sig":   "Additional_file_33_GO_overrepresentation_significant.csv",
    "go_all":   "Additional_file_34_GO_overrepresentation_all_tests.csv.gz",
    "diptera":  "Additional_file_19_Diptera_chemosensory_comparison.csv",
    "unknown":  "Additional_file_27_Unknown_BLAST.csv",
    "clusters": "Additional_file_35_identical_protein_clusters.csv",
    "locus":    "Additional_file_37_transcript_model_to_gene_locus.csv",
    "fig3stat": "Additional_file_38_Figure3_statistical_tests.csv",
}


@st.cache_data(show_spinner=False)
def ex_load(key):
    path = os.path.join(SUPP_DIR, EXPLORER_FILES[key])
    if not os.path.exists(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def ex_download(df, label, fname, key=None):
    st.download_button(label, df.to_csv(index=False).encode("utf-8"),
                       file_name=fname, mime="text/csv",
                       key=f"dl_{key or fname}")


def ex_table(df, key, fname, height=360, caption=None):
    """Show a dataframe with a free-text filter across every column, and a
    download of exactly what the filter leaves."""
    if df is None or df.empty:
        st.info("Nothing to show.")
        return df
    q = st.text_input("Filter rows (searches every column)", "",
                      key=f"flt_{key}", placeholder="gene, term, family, accession…")
    sub = df
    if q.strip():
        terms = [t for t in q.lower().split() if t]
        hay = df.astype(str).apply(lambda c: c.str.lower())
        mask = pd.Series(True, index=df.index)
        for t in terms:                      # every term must appear somewhere
            mask &= hay.apply(lambda c: c.str.contains(t, regex=False, na=False)).any(axis=1)
        sub = df[mask]
    st.caption(caption or f"{len(sub):,} of {len(df):,} rows")
    st.dataframe(sub, use_container_width=True, hide_index=True, height=height)
    ex_download(sub, "Download these rows (CSV)", fname, key=f"tbl_{key}")
    return sub


# -------------------------------------------------------------- GO over-rep
def ex_tab_go():
    st.subheader("Gene Ontology over-representation")
    st.markdown(
        "Every GO term tested against its gene set, with the odds ratio and the "
        "Benjamini–Hochberg adjusted p value. Figures S9–S12 and S20–S24 plot the "
        "strongest terms; this is the whole table. Tests are one-sided Fisher "
        "exact tests, corrected within each gene set × GO domain."
    )
    sig = ex_load("go_sig")
    if sig.empty:
        st.info("Additional file 33 not found.")
        return
    show_all = st.checkbox(
        "Include tests that did not reach significance (Additional file 34)",
        value=False, key="ex_go_all")
    df = ex_load("go_all") if show_all else sig
    if df.empty:
        st.info("Additional file 34 not found.")
        return
    c1, c2, c3 = st.columns(3)
    with c1:
        sets = sorted(df["Gene_set"].dropna().unique())
        pick = st.multiselect("Gene set", sets, default=sets[:1], key="ex_go_set")
    with c2:
        doms = sorted(df["GO_domain"].dropna().unique())
        dpick = st.multiselect("GO domain", doms, default=doms, key="ex_go_dom")
    with c3:
        qmax = st.number_input("Maximum q (BH)", value=0.05, min_value=1e-12,
                               max_value=1.0, format="%.4g", key="ex_go_q")
    term_q = st.text_input("Filter GO term contains", "", key="ex_go_term")
    sub = df[df["Gene_set"].isin(pick) & df["GO_domain"].isin(dpick)]
    sub = sub[pd.to_numeric(sub["q_BH"], errors="coerce") <= qmax]
    if term_q.strip():
        sub = sub[sub["GO_term"].str.contains(term_q.strip(), case=False, na=False)]
    if sub.empty:
        st.info("No terms match these filters.")
        return
    sub = sub.sort_values("q_BH")
    st.caption(f"{len(sub):,} terms")
    top = sub.head(25).iloc[::-1]
    fig = go.Figure(go.Bar(
        x=top["odds_ratio"], y=top["GO_term"], orientation="h",
        marker=dict(color=top["odds_ratio"], colorscale="Blues", showscale=False),
        customdata=np.stack([top["q_BH"], top["k_in_set"], top["n_set"]], axis=-1),
        hovertemplate=("<b>%{y}</b><br>odds ratio %{x:.2f}<br>"
                       "q = %{customdata[0]:.3g}<br>"
                       "%{customdata[1]} of %{customdata[2]} genes<extra></extra>"),
    ))
    fig.update_layout(height=max(320, 22 * len(top) + 90), plot_bgcolor="white",
                      margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title="Odds ratio", yaxis_title=None,
                      title=dict(text="Strongest 25 terms", font=dict(size=12)))
    st.plotly_chart(fig, use_container_width=True, key="ex_go_fig")
    ex_table(sub, "go", "GO_overrepresentation_filtered.csv", height=360)


# ------------------------------------------------- chemosensory per-locus view
# Every contrast in the paper, restricted to the 567 chemosensory transcript
# models.  Family membership comes from Additional file 37 rather than from
# string-matching a name, so the counts here are the published ones.
DE_CONTRASTS = {
    "Antenna vs Tarsi":              ("condition_vs_Ant_vs_Leg_name.csv", "appendage"),
    "Antenna vs Maxillary palp":     ("condition_vs_Ant_vs_P_name.csv",   "appendage"),
    "Tarsi vs Maxillary palp":       ("condition_vs_Leg_vs_P_name.csv",   "appendage"),
    "Antenna — sex (VF vs Vm)":          ("results_ant_VF_vs_Vm.csv",  "sex"),
    "Maxillary palp — sex (VF vs Vm)":   ("results_palp_VF_vs_Vm.csv", "sex"),
    "Tarsi — sex (VF vs Vm)":            ("results_leg_VF_vs_Vm.csv",  "sex"),
    "Antenna — mating (MF vs VF)":        ("results_ant_MF_vs_VF.csv",  "mating"),
    "Maxillary palp — mating (MF vs VF)": ("results_palp_MF_vs_VF.csv", "mating"),
    "Tarsi — mating (MF vs VF)":          ("results_leg_MF_vs_VF.csv",  "mating"),
    "Antenna — mated female vs virgin male":        ("results_ant_MF_vs_Vm.csv",  "matedsex"),
    "Maxillary palp — mated female vs virgin male": ("results_palp_MF_vs_Vm.csv", "matedsex"),
    "Tarsi — mated female vs virgin male":          ("results_leg_MF_vs_Vm.csv",  "matedsex"),
}


@st.cache_data(show_spinner=False)
def ex_family_map():
    loci = ex_load("locus")
    if loci.empty:
        return {}, {}
    fam = dict(zip(loci["Transcript_ID"].astype(str), loci["Family"].astype(str)))
    sym = dict(zip(loci["Transcript_ID"].astype(str), loci["Hill_gene_symbol"].astype(str)))
    return fam, sym


@st.cache_data(show_spinner=False)
def ex_load_contrast(fname, padj_thr, lfc_thr, chemo_only):
    path = os.path.join(BASE_DIR, fname)
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, usecols=lambda c: c in
                    ("Gene", "Name", "log2FoldChange", "padj", "baseMean"))
    d["Gene"] = d["Gene"].astype(str)
    d["lfc"] = -pd.to_numeric(d["log2FoldChange"], errors="coerce")   # see README
    d["padj"] = pd.to_numeric(d["padj"], errors="coerce")
    d["sig"] = (d["padj"] < padj_thr) & (d["lfc"].abs() >= lfc_thr)
    fam, sym = ex_family_map()
    d["Family"] = d["Gene"].map(fam)
    d["Symbol"] = d["Gene"].map(sym)
    if chemo_only:
        d = d[d["Family"].notna()]
    return d.drop(columns=["log2FoldChange"], errors="ignore")


# Figure 4A/4B colour every stem by the appendage it belongs to:
# antenna dark blue, maxillary palp mid blue, tarsi grey.
APPENDAGE_COLS = {"Antenna": "#08306B", "Maxillary palp": "#4292C6", "Tarsi": "#737373"}

def ex_lollipop(d, title, n_max, colour_appendage=None):
    """Stems from no change to the log2 fold change, as in Figure 4A.

    Dot colour is the appendage, matching the published panels. For an
    appendage-versus-appendage contrast there is no single appendage, so the
    dot takes the colour of whichever appendage the gene is higher in.
    """
    d = d.sort_values("lfc")
    if len(d) > n_max:
        keep = pd.concat([d.head(n_max // 2), d.tail(n_max - n_max // 2)])
        d = keep.sort_values("lfc")
    lab = d["Symbol"].fillna(d["Gene"]).astype(str)
    fig = go.Figure()
    if colour_appendage is None:
        cols = pd.Series("#737373", index=d.index)
    elif isinstance(colour_appendage, tuple):          # (positive side, negative side)
        hi, lo = colour_appendage
        cols = pd.Series(np.where(d["lfc"] > 0, APPENDAGE_COLS.get(hi, "#737373"),
                                  APPENDAGE_COLS.get(lo, "#737373")), index=d.index)
    else:
        cols = pd.Series(APPENDAGE_COLS.get(colour_appendage, "#737373"), index=d.index)
    for x, y, c in zip(d["lfc"], lab, cols):
        fig.add_shape(type="line", x0=0, x1=x, y0=y, y1=y,
                      line=dict(color="#C9C9C9", width=1))
    size = 8.0
    if "baseMean" in d.columns and d["baseMean"].notna().any():
        bm = pd.to_numeric(d["baseMean"], errors="coerce").fillna(0).clip(lower=1)
        size = 5 + 9 * (np.log10(bm) - np.log10(bm).min()) / max(
            np.log10(bm).max() - np.log10(bm).min(), 1e-9)
    fig.add_trace(go.Scatter(
        x=d["lfc"], y=lab, mode="markers",
        marker=dict(size=size, color=list(cols),
                    line=dict(width=0.5, color="#3A3A3A")),
        customdata=np.stack([d["Gene"], d["padj"].fillna(1),
                             d.get("Family", pd.Series("", index=d.index))], axis=-1),
        hovertemplate=("<b>%{y}</b><br>log₂FC %{x:.2f}<br>"
                       "%{customdata[0]} · %{customdata[2]}<br>"
                       "adjusted p %{customdata[1]:.3g}<extra></extra>")))
    fig.add_vline(x=0, line_width=1, line_color="#555555")
    fig.update_layout(height=max(340, 15 * len(d) + 110), plot_bgcolor="white",
                      margin=dict(l=10, r=10, t=40, b=10), showlegend=False,
                      xaxis_title="log₂ fold change",
                      title=dict(text=title, font=dict(size=12)))
    return fig


def ex_tab_chemo(by_tissue, padj_thr, lfc_thr):
    st.subheader("Response locus by locus")
    st.markdown(
        "The data behind the **Figure 4A** lollipops, extended to every contrast "
        "in the paper — including the three appendage comparisons. Each stem runs "
        "from no change to the gene's log₂ fold change; dot area is mean "
        "normalised expression. Pick a contrast, a family, and whether to look "
        "only at chemosensory genes."
    )
    c1, c2 = st.columns([3, 2])
    with c1:
        picks = st.multiselect("Contrast", list(DE_CONTRASTS.keys()),
                               default=["Antenna — sex (VF vs Vm)"], key="ex_ll_con")
    with c2:
        chemo_only = st.checkbox("Chemosensory genes only", value=True, key="ex_ll_chemo")
    if not picks:
        st.info("Choose at least one contrast.")
        return
    frames = []
    for name in picks:
        fname, layer = DE_CONTRASTS[name]
        d = ex_load_contrast(fname, padj_thr, lfc_thr, chemo_only)
        if d.empty:
            continue
        d = d.copy()
        d["Appendage"] = name.split(" — ")[0]
        d["Contrast"] = name
        frames.append(d)
    if not frames:
        st.info("Those contrasts are not available.")
        return
    D = pd.concat(frames, ignore_index=True)
    c3, c4, c5 = st.columns(3)
    with c3:
        fams = sorted(D["Family"].dropna().unique())
        fpick = st.multiselect("Family", fams, default=fams, key="ex_ll_fam") if fams else []
    with c4:
        sig_only = st.checkbox("Significant only", value=True, key="ex_ll_sig")
    with c5:
        n_max = int(st.number_input("Maximum loci to plot", min_value=10, max_value=200,
                                    value=50, step=10, key="ex_ll_n"))
    sub = D.copy()
    if fpick:
        sub = sub[sub["Family"].isin(fpick) | sub["Family"].isna()]
    if chemo_only and fpick:
        sub = sub[sub["Family"].isin(fpick)]
    if sig_only:
        sub = sub[sub["sig"]]
    if sub.empty:
        st.info("Nothing matches these filters. Try clearing 'Significant only'.")
        return
    st.caption(f"{len(sub):,} loci across {sub['Contrast'].nunique()} contrast(s)")
    pair = {"Antenna vs Tarsi": ("Antenna", "Tarsi"),
            "Antenna vs Maxillary palp": ("Antenna", "Maxillary palp"),
            "Tarsi vs Maxillary palp": ("Tarsi", "Maxillary palp")}
    for name in picks:
        part = sub[sub["Contrast"] == name]
        if part.empty:
            continue
        colour = pair.get(name) or name.split(" \u2014 ")[0]
        st.plotly_chart(
            ex_lollipop(part.drop(columns=["Appendage"]), name, n_max, colour),
            use_container_width=True, key=f"ex_ll_{name}")
    st.caption(
        "Dot colour is the appendage, as in Figure 4A and 4B: antenna dark blue, "
        "maxillary palp mid blue, tarsi grey. For an appendage-versus-appendage "
        "contrast the dot takes the colour of the appendage the gene is higher in. "
        "Dot area is mean normalised expression."
    )
    out = sub[["Contrast", "Gene", "Symbol", "Family", "lfc", "padj", "sig", "baseMean"]]
    out = out.rename(columns={"lfc": "log2FC", "padj": "adjusted_p",
                              "sig": "significant", "baseMean": "mean_normalised_count"})
    ex_table(out, "lolli", "locus_response.csv", height=340)


# ----------------------------------------- S26: mated-female-specific expression
def ex_tab_mfs(tissue_cross):
    st.subheader("Mated-female-specific expression")
    st.markdown(
        "The interactive form of **Figure S26**. Sex bias in virgins (x) against "
        "sex bias in mated females (y). The diagonal is the null for the mating "
        "contrast: a gene sitting on it has the same sex bias in both groups, so "
        "**y = x means no mating response**. Vertical displacement from the "
        "diagonal is the mating effect, and colour says whether that displacement "
        "passed the strict test against both virgin groups."
    )
    if not tissue_cross:
        st.info("Contrasts not loaded.")
        return
    L = RETENTION_LIM
    cols = st.columns(3)
    tables, genes = [], []
    for ci, tissue in enumerate(["Antenna", "Palp", "Tarsi"]):
        td = tissue_cross.get(tissue)
        d = s1_build_retention(tissue_cross, tissue)
        if d is None or td is None:
            with cols[ci]:
                st.info(f"Not available for {TISSUE_DISPLAY[tissue]}.")
            continue
        up, dn = td["mat_mf"], td["mat_vf"]
        d = d.copy()
        key = d["JoinKey"].astype(str)
        d["MFS"] = np.where(key.isin(up), "Induced in mated females",
                   np.where(key.isin(dn), "Reduced in mated females",
                            "Not mating-responsive"))
        colmap = {"Induced in mated females": "#C0392B",
                  "Reduced in mated females": "#2471A3",
                  "Not mating-responsive": "#D5D5D5"}
        fig = go.Figure()
        fig.add_shape(type="line", x0=-L, x1=L, y0=-L, y1=L,
                      line=dict(color="#555555", width=1.4), layer="below")
        for cls in ("Not mating-responsive", "Reduced in mated females",
                    "Induced in mated females"):
            part = d[d["MFS"] == cls]
            if part.empty:
                continue
            fig.add_trace(go.Scattergl(
                x=part["px"], y=part["py"], mode="markers", name=cls,
                marker=dict(size=4.5, color=colmap[cls]),
                customdata=np.stack([part["JoinKey"], part["x"], part["y"]], axis=-1),
                hovertemplate=("<b>%{customdata[0]}</b><br>virgins %{customdata[1]:.2f}"
                               "<br>mated %{customdata[2]:.2f}<extra></extra>")))
        fig.update_layout(
            title=dict(text=TISSUE_DISPLAY[tissue], font=dict(size=13)),
            xaxis=dict(title="Sex bias in virgins (log₂FC)", range=[-L, L]),
            yaxis=dict(title="Sex bias in mated females (log₂FC)", range=[-L, L],
                       scaleanchor="x", scaleratio=1),
            legend=dict(orientation="h", yanchor="top", y=-0.18, x=0, font=dict(size=9)),
            margin=dict(l=10, r=10, t=35, b=10), height=520, plot_bgcolor="white")
        with cols[ci]:
            st.plotly_chart(fig, use_container_width=True, key=f"ex_mfs_{tissue}")
        t = d.groupby("MFS").agg(
            genes=("JoinKey", "size"),
            chemosensory=("Class", lambda c: int((c == CLS_CHEMO).sum()))).reset_index()
        t.insert(0, "Appendage", TISSUE_DISPLAY[tissue])
        tables.append(t)
        genes.append(d.assign(Appendage=TISSUE_DISPLAY[tissue])[
            ["Appendage", "JoinKey", "ChemoName", "x", "y", "MFS"]].rename(
            columns={"JoinKey": "Gene", "ChemoName": "Chemosensory name",
                     "x": "log2FC_virgins", "y": "log2FC_mated"}))
    if tables:
        summ = pd.concat(tables, ignore_index=True)
        st.dataframe(summ, use_container_width=True, hide_index=True)
        ex_download(summ, "Download these counts (CSV)",
                    "S26_matedfemale_specific_counts.csv")
    if genes:
        allg = pd.concat(genes, ignore_index=True)
        with st.expander(f"All {len(allg):,} plotted genes"):
            st.dataframe(allg.round(3), use_container_width=True, hide_index=True,
                         height=360)
            ex_download(allg, "Download every plotted gene (CSV)",
                        "S26_matedfemale_specific_genes.csv")


# ------------------------------------------------------------- gene lookup
def ex_tab_gene(by_tissue):
    st.subheader("Look up a gene")
    st.markdown(
        "Every contrast for one gene, side by side, with its BLAST evidence. "
        "Search by accession, chemosensory name or description."
    )
    unk = ex_load("unknown")
    loci = ex_load("locus")
    q = st.text_input("Gene accession, name or description", "", key="ex_gene_q")
    if not q.strip():
        st.caption("Try ORco, OBP34, or an XM_ accession.")
        return
    qq = q.strip().lower()
    hits = set()
    for _, items in by_tissue.items():
        for it in items:
            d = it["df"]
            m = (d["JoinKey"].astype(str).str.lower().str.contains(qq, na=False)
                 | d["ChemoName"].astype(str).str.lower().str.contains(qq, na=False)
                 | d.get("Name", pd.Series("", index=d.index)).astype(str)
                   .str.lower().str.contains(qq, na=False))
            hits.update(d.loc[m, "JoinKey"].astype(str).tolist())
    if not hits:
        st.warning("No gene matches that search.")
        return
    hits = sorted(hits)
    if len(hits) > 1:
        gene = st.selectbox(f"{len(hits)} matches", hits, key="ex_gene_pick")
    else:
        gene = hits[0]
    rows = []
    for tissue, items in by_tissue.items():
        for it in items:
            d = it["df"]
            r = d[d["JoinKey"].astype(str) == gene]
            if r.empty:
                continue
            r = r.iloc[0]
            rows.append({
                "Appendage": TISSUE_DISPLAY[tissue],
                "Contrast": f'{it["cond1"]} vs {it["cond2"]}',
                "log₂FC": round(float(r["log2FoldChange"]), 3),
                "adjusted p": float(r["padj"]) if pd.notna(r["padj"]) else None,
                "Significant": bool(r["is_sig"]),
            })
    if rows:
        rdf = pd.DataFrame(rows)
        st.dataframe(rdf, use_container_width=True, hide_index=True)
        ex_download(rdf, "Download this gene's contrasts (CSV)", f"{gene}_contrasts.csv")
    if not loci.empty:
        lr = loci[loci["Transcript_ID"].astype(str) == gene]
        if not lr.empty:
            r = lr.iloc[0]
            st.caption(f"Family **{r['Family']}** · symbol **{r['Hill_gene_symbol']}** "
                       f"· locus **{r['gene_id']}**")
    if not unk.empty:
        ur = unk[unk["Gene"].astype(str) == gene]
        if not ur.empty:
            with st.expander("BLAST evidence and category"):
                st.dataframe(ur.T.rename(columns=lambda c: "value"),
                             use_container_width=True)


# ------------------------------------------------------------- comparative
def ex_tab_comparative():
    st.subheader("Comparative and structural data")
    dip = ex_load("diptera")
    if not dip.empty:
        st.markdown("**Chemosensory family sizes across Diptera** (Additional file 19)")
        long = dip.melt(id_vars="Family", var_name="Species", value_name="Genes")
        long["Species"] = long["Species"].str.replace("_", ". ", regex=False)
        fig = px.bar(long, x="Family", y="Genes", color="Species", barmode="group",
                     color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(height=380, plot_bgcolor="white",
                          margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(fig, use_container_width=True, key="ex_dip_fig")
        ex_table(dip, "dip", "AF19_Diptera_comparison.csv", height=260)
    cl = ex_load("clusters")
    if not cl.empty:
        st.markdown("**Identical protein clusters** (Additional file 35)")
        st.caption(
            f"{len(cl)} clusters of transcript models encoding identical proteins. "
            "These are why the paper reports 567 transcript models but 393 gene loci."
        )
        ex_table(cl, "clusters", "AF35_identical_protein_clusters.csv", height=260)
    f3 = ex_load("fig3stat")
    if not f3.empty:
        st.markdown("**Enrichment tests behind Figure 3** (Additional file 38)")
        ex_table(f3, "fig3", "AF38_Figure3_statistical_tests.csv", height=300)



# --------------------------------------------------------------- statistics
STATS_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "stats_reports")

STATS_TITLES = {
    "00_SUMMARY_permanova_silhouette_betadisper.csv": "Summary: PERMANOVA, silhouette, dispersion",
    "01_pca_overview.csv": "PCA: variance explained per component",
    "02_pca_group_geometry_and_silhouette.csv": "PCA: group geometry and silhouette",
    "03_permanova_summary.csv": "PERMANOVA: main effects",
    "04a_pairwise_permanova_tests.csv": "PERMANOVA: pairwise tests",
    "04b_pairwise_mean_distances.csv": "Pairwise mean distances",
    "05_pearson_correlation_long.csv": "Pearson correlations, long form",
    "06_pearson_correlation_summary.csv": "Pearson correlations, summary",
    "07_variation_partitioning.csv": "Variation partitioning",
}


def ex_tab_stats():
    st.subheader("Multivariate statistics")
    st.markdown(
        "The statistical reports behind Figure 2 (Additional file 32): how much "
        "variance appendage identity, sex and mating status each explain, whether "
        "the groups separate, and by how much. Every table downloads as CSV."
    )
    if not os.path.isdir(STATS_DIR):
        st.info("Additional file 32 not found.")
        return
    avail = [f for f in sorted(os.listdir(STATS_DIR)) if f.endswith(".csv")]
    if not avail:
        st.info("No statistical reports found.")
        return
    labels = {STATS_TITLES.get(f, f): f for f in avail}
    pick = st.selectbox("Report", list(labels.keys()), key="ex_stats_pick")
    fname = labels[pick]
    df = pd.read_csv(os.path.join(STATS_DIR, fname))
    ex_table(df, f"stats_{fname}", fname, height=420)
    with st.expander("All reports in one download"):
        for f in avail:
            d = pd.read_csv(os.path.join(STATS_DIR, f))
            st.caption(f"{STATS_TITLES.get(f, f)} — {len(d):,} rows")
            ex_download(d, f"Download {f}", f, key=f"all_{f}")


# ------------------------------------------------- amino-acid identity matrices
IDENTITY_FILES = {
    "CSP": "Additional_file_12_AA_identity_matrix_CSP.csv",
    "GR":  "Additional_file_13_AA_identity_matrix_GR.csv",
    "IR":  "Additional_file_14_AA_identity_matrix_IR.csv",
    "OBP": "Additional_file_15_AA_identity_matrix_OBP.csv",
    "OR":  "Additional_file_16_AA_identity_matrix_OR.csv",
    "PPK": "Additional_file_17_AA_identity_matrix_PPK.csv",
    "TRP": "Additional_file_18_AA_identity_matrix_TRP.csv",
    "1:1 orthologs": "Additional_file_28_AA_identity_1to1_orthologs.csv",
}


@st.cache_data(show_spinner=False)
def ex_load_identity(fname):
    path = os.path.join(SUPP_DIR, fname)
    if not os.path.exists(path):
        return pd.DataFrame()
    return pd.read_csv(path)


def ex_tab_identity():
    st.subheader("Amino-acid identity")
    st.markdown(
        "Pairwise amino-acid identity within each chemosensory family "
        "(Additional files 12\u201318) and between the 1:1 orthologs "
        "(Additional file 28). Use it to ask how similar two paralogs really are: "
        "a bright block is a recent expansion, a dark field is an old one."
    )
    c1, c2 = st.columns([2, 3])
    with c1:
        fam = st.selectbox("Family", list(IDENTITY_FILES.keys()), key="ex_id_fam")
    df = ex_load_identity(IDENTITY_FILES[fam])
    if df.empty:
        st.info("That matrix is not available.")
        return
    first = df.columns[0]
    names = df[first].astype(str)
    mat = df.drop(columns=[first]).apply(pd.to_numeric, errors="coerce")
    with c2:
        pick = st.text_input(
            "Restrict to proteins whose name contains", "", key="ex_id_sel",
            placeholder="e.g. OR1  — leave blank for the whole family")
    if pick.strip():
        keep = names.str.contains(pick.strip(), case=False, na=False)
        if keep.sum() < 2:
            st.warning("Fewer than two proteins match; showing the whole family.")
        else:
            cols = [c for c in mat.columns if str(c) in set(names[keep])]
            mat = mat.loc[keep.values, cols] if cols else mat.loc[keep.values]
            names = names[keep]
    n = len(names)
    if n > 160:
        st.caption(f"{n} proteins — the heat map is drawn without labels at this size.")
    fig = px.imshow(
        mat.values,
        x=list(names) if n <= 160 else None,
        y=list(names) if n <= 160 else None,
        color_continuous_scale="Viridis", origin="upper", aspect="auto",
        labels=dict(color="% identity"),
    )
    fig.update_traces(hovertemplate="%{y} vs %{x}<br>%{z:.1f}% identity<extra></extra>")
    fig.update_layout(height=max(420, min(820, 9 * n + 140)),
                      margin=dict(l=10, r=10, t=36, b=10),
                      title=dict(text=f"{fam} — pairwise amino-acid identity",
                                 font=dict(size=12)),
                      xaxis=dict(tickfont=dict(size=7), showticklabels=n <= 160),
                      yaxis=dict(tickfont=dict(size=7), showticklabels=n <= 160))
    st.plotly_chart(fig, use_container_width=True, key="ex_id_fig")
    vals = mat.where(~np.eye(len(mat), dtype=bool)[:len(mat), :mat.shape[1]]) \
        if mat.shape[0] == mat.shape[1] else mat
    flat = pd.to_numeric(pd.Series(vals.values.ravel()), errors="coerce").dropna()
    if len(flat):
        st.caption(
            f"{n} proteins \u00b7 off-diagonal identity: median {flat.median():.1f}%, "
            f"range {flat.min():.1f}\u2013{flat.max():.1f}%"
        )
    ex_table(df, f"ident_{fam}", IDENTITY_FILES[fam], height=320)


def render_explorer_tab():
    st.header("Data Explorer")
    st.caption(
        "Every table behind the manuscript — searchable, filterable and "
        "downloadable as CSV. Significance thresholds are yours to set."
    )
    padj_thr = st.sidebar.number_input(
        "padj threshold", value=0.001, min_value=1e-10, max_value=0.1,
        format="%.4g", key="ex_padj_thr")
    lfc_thr = st.sidebar.number_input(
        "|log2FC| threshold", value=1.0, min_value=0.0, max_value=10.0,
        step=0.25, key="ex_lfc_thr")
    st.sidebar.caption("Manuscript defaults: adjusted p < 0.001, |log₂FC| ≥ 1.")
    by_tissue, _ = s1_load_all_contrasts(BASE_DIR, padj_thr, lfc_thr, 2.5)
    tissue_cross = s1_compute_sex_mating_overlap(by_tissue)
    t1, t2, t3, t4, t5, t6, t7 = st.tabs(
        ["Response (lollipops)", "Mated-female-specific", "GO over-representation",
         "Gene lookup", "Multivariate statistics", "Amino-acid identity",
         "Comparative data"])
    with t1:
        ex_tab_chemo(by_tissue, padj_thr, lfc_thr)
    with t2:
        ex_tab_mfs(tissue_cross)
    with t3:
        ex_tab_go()
    with t4:
        ex_tab_gene(by_tissue)
    with t5:
        ex_tab_stats()
    with t6:
        ex_tab_identity()
    with t7:
        ex_tab_comparative()


# =============================================================================
# === PUBLIC ADDRESS ===
# The deployed address, also cited in the manuscript.  Override it with the
# APP_URL environment variable when running a private or local instance.
# =============================================================================
APP_URL = os.environ.get(
    "APP_URL",
    "https://hillucensolfactomeapp-rk3chuhczmmvbpfq2ctvan.streamlit.app/",
)


def render_open_on_desktop():
    """Small sidebar helper: the browser is cramped on a phone, so offer an
    easy way to move the link to a computer."""
    with st.sidebar.expander("Open on a computer", expanded=True):
        st.caption(
            "The heatmaps, trees and volcano plots need a wide screen. "
            "Mail yourself the link and open it on a desktop."
        )
        subject = "H. illucens olfactome browser"
        body = (
            "Open this on a computer for the full interactive view:\n\n"
            f"{APP_URL}\n\n"
            "Data browser for Perets et al. 2026, chemosensory transcriptomes of "
            "the antenna, maxillary palp and tarsi of Hermetia illucens."
        )
        href = f"mailto:?subject={quote(subject)}&body={quote(body)}"
        st.markdown(
            f'<a href="{href}" target="_blank" rel="noopener" '
            'style="display:inline-block;padding:0.45rem 0.9rem;border-radius:0.4rem;'
            'background:#01045A;color:#fff;text-decoration:none;font-weight:600;">'
            'Email me this link</a>',
            unsafe_allow_html=True,
        )
        st.caption("Or copy it:")
        st.code(APP_URL, language=None)


render_open_on_desktop()

section = st.sidebar.radio(
    "Section",
    [
        "Phylogenetic Trees",
        "Transcriptome Overview",
        "Appendage Comparison",
        "Sex & Mating Analysis",
        "Chemosensory Heatmap",
        "Data Explorer",
        "Figures & Data",
    ],
    key="main_section",
)

if section == "Phylogenetic Trees":
    render_trees_tab()
elif section == "Transcriptome Overview":
    render_t1_tab()
elif section == "Appendage Comparison":
    render_c1_tab()
elif section == "Sex & Mating Analysis":
    render_s1_tab()
elif section == "Chemosensory Heatmap":
    render_h1_tab()
elif section == "Data Explorer":
    render_explorer_tab()
else:
    render_figures_tab()

    
st.markdown("---")
st.markdown(
    "<div style='font-size:18px; color:#444;'>"
    "<b>All rights of this data are for Jonathan Bohbot Lab.</b> "
    "Use of this data for research or commercial use must contact Jonathan Bohbot Lab."
    "</div>",
    unsafe_allow_html=True,
)
