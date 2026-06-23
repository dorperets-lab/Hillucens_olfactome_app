# =============================================================================
# app.py  –  BSF Transcriptome Explorer  (merged from app_t1, c1, s1, h1)
# =============================================================================

# === IMPORTS ===
import base64
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
    page_title="H. illucens Olfactome | Perets et al. 2026",
    page_icon="🪰",
    layout="wide",
)
st.markdown("# Appendage-resolved transcriptomics reveals constitutive sex-biased expression and reproductive-state plasticity in adult *Hermetia illucens* chemosensory appendages")
st.caption("Perets et al. 2026 · *Genome Biology* · Interactive data companion")

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

VP_PIE_COLORS = {
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
        "Palp - MF vs VF":    "results_palp_MF_vs_VF.csv",
        "Palp - VF vs Vm":    "results_palp_VF_vs_Vm.csv",
        "Tarsi - MF vs VF":   "results_leg_MF_vs_VF.csv",
        "Tarsi - VF vs Vm":   "results_leg_VF_vs_Vm.csv",
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
        color_continuous_scale="Blues", range_color=(0.0, 1.0), size_max=20,
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

    s = norm_counts[name_col].astype(str).str.upper()

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
            "Csp": col_p2, "Ppk": col_p2,
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
        sub = df[mask & df["model"].str.contains("2-way", na=False)].copy()
        if sub.empty:
            return pd.DataFrame()
        comp_map = {
            "Appendage unique": "Appendage (unique effect)",
            "Sex / mating-status unique": "Reproductive state (unique effect)",
            "Unexplained": "Unexplained variance (other factors)",
        }
        rows = []
        for _, r in sub.iterrows():
            comp = comp_map.get(r["component"].strip(), r["component"].strip())
            rows.append({"Component": comp, "Fraction": max(0.0, float(r["adj_r2_fraction"]))})
        return pd.DataFrame(rows)

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
        _render_vp_pie(pie_all, "All genes (n = 24,608; Ezekiel adj. R²; RDA, vegan)", "all")
    with col_vp2:
        _render_vp_pie(pie_chemo, "Chemosensory genes only (n ≈ 553)", "chemo")

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


def c1_fig_volcano(vdf, padj_thr, lfc_thr, overlay_chemo, highlight_query):
    fig = px.scatter(
        vdf, x="log2FoldChange", y="-log10", color="ColorKey",
        color_discrete_map=PAL_GO,
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
    cols = st.columns(3)
    with cols[0]:
        st.markdown(f"**{c1_volcano_title(antp_v)}**")
        fig = c1_fig_volcano(antp_v, padj_thr, lfc_thr, show_chemo, highlight_query)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=antp_v.to_csv(index=False).encode(),
                           file_name="antenna_vs_maxillary_palp_volcano_table.csv",
                           mime="text/csv", key="c1_dl_antp")
    with cols[1]:
        st.markdown(f"**{c1_volcano_title(antleg_v)}**")
        fig = c1_fig_volcano(antleg_v, padj_thr, lfc_thr, show_chemo, highlight_query)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=antleg_v.to_csv(index=False).encode(),
                           file_name="antenna_vs_tarsi_volcano_table.csv",
                           mime="text/csv", key="c1_dl_antleg")
    with cols[2]:
        st.markdown(f"**{c1_volcano_title(legp_v)}**")
        fig = c1_fig_volcano(legp_v, padj_thr, lfc_thr, show_chemo, highlight_query)
        st.plotly_chart(fig, use_container_width=True)
        st.download_button("Download table (filtered)",
                           data=legp_v.to_csv(index=False).encode(),
                           file_name="tarsi_vs_maxillary_palp_volcano_table.csv",
                           mime="text/csv", key="c1_dl_legp")


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
        ["Volcano (Fig. 3A)", "GO domain % (Fig. 3B)", "GO Names (Fig. 3C)", "Chemosensory pies (Fig. 3D)"]
    )
    with vol_tab:
        c1_render_volcano_tab(volcano_dfs, padj_thr, lfc_thr, show_chemo, highlight_query)
    with go_tab:
        c1_render_go_tab(volcano_dfs, norm_means, show_chemo_go)
    with overlap_tab:
        c1_render_overlap_tab(de_data, volcano_dfs, go_map, name_map, norm_means, padj_thr, lfc_thr)
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
        sex_item = next((x for x in contrast_list if x["cond1"] == "VF" and x["cond2"] == "Vm"), None)
        mat_item = next((x for x in contrast_list if x["cond1"] == "MF" and x["cond2"] == "VF"), None)
        if sex_item is None or mat_item is None:
            continue
        sex_df = sex_item["df"]
        mat_df = mat_item["df"]
        sex_vf = set(sex_df.loc[sex_df["is_sig"] & (sex_df["Side"] == "VF"), "JoinKey"].astype(str))
        sex_vm = set(sex_df.loc[sex_df["is_sig"] & (sex_df["Side"] == "Vm"), "JoinKey"].astype(str))
        mat_mf = set(mat_df.loc[mat_df["is_sig"] & (mat_df["Side"] == "MF"), "JoinKey"].astype(str))
        mat_vf = set(mat_df.loc[mat_df["is_sig"] & (mat_df["Side"] == "VF"), "JoinKey"].astype(str))
        overlap = (sex_vf | sex_vm) & (mat_mf | mat_vf)
        result[tissue] = {
            "sex_vf": sex_vf, "sex_vm": sex_vm,
            "mat_mf": mat_mf, "mat_vf": mat_vf,
            "overlap": overlap,
            "ov_female": sex_vf & overlap,
            "ov_male": sex_vm & overlap,
            "sex_df": sex_df, "mat_df": mat_df,
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
        ["Volcano (Figs. 4A–B)", "DE summary bars + DEG tables", "Venn diagrams (Fig. 4C)",
         "Scatter & Overlap (Figs. 4D–E)", "GO Term Browser (Figs. S14–S16)"]
    )

    tissue_labels_s1 = {"Antenna": "Antenna", "Palp": "Maxillary palp", "Tarsi": "Tarsi"}

    # -------------------- VOLCANO TAB --------------------
    with tab_volc:
        st.subheader("Volcano plots (two per tissue)")
        st.caption("Manuscript: **Figure 4A** (VF vs Vm, sex-biased) · **Figure 4B** (MF vs VF, mating-regulated)")
        for tissue in ["Antenna", "Palp", "Tarsi"]:
            if tissue not in by_tissue:
                continue
            st.markdown(f"### {TISSUE_DISPLAY[tissue]}")
            cols = st.columns(2)
            contrasts = sorted(by_tissue[tissue], key=s1_sort_key)
            for i, item in enumerate(contrasts[:2]):
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
            cols = st.columns(2)
            contrasts = sorted(by_tissue[tissue], key=s1_sort_key)
            for i, item in enumerate(contrasts[:2]):
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
        st.subheader("Venn diagrams - tissue overlap of upregulated genes")
        st.caption("Manuscript: **Figure 4C**")
        st.markdown(
            "_Venns use fixed thresholds: padj < 0.001 and |log2FC| >= 1.0 for upregulation._"
        )
        region_options = [
            "Antenna only", "Palp only", "Tarsi only",
            "Antenna & Palp", "Antenna & Tarsi", "Palp & Tarsi", "All 3 tissues",
        ]
        region_key_map = {
            "Antenna only": "A only", "Palp only": "B only", "Tarsi only": "C only",
            "Antenna & Palp": "A & B", "Antenna & Tarsi": "A & C",
            "Palp & Tarsi": "B & C", "All 3 tissues": "A & B & C",
        }

        row1 = st.columns(2)
        with row1[0]:
            state_code = "Vm"
            state_full = COND_FULL.get(state_code, state_code)
            state_col = s1_get_full_label_color(state_full)
            st.markdown(f"### {state_full}")
            sets_vm = s1_get_up_sets_venn(by_tissue, state_code="Vm", partner_code=None)
            setA = sets_vm.get("Antenna", set()); setB = sets_vm.get("Palp", set()); setC = sets_vm.get("Tarsi", set())
            regions = s1_compute_venn_regions(setA, setB, setC)
            counts = {k: len(v) for k, v in regions.items()}
            st.caption(f"Antenna: {len(setA)} | Maxillary palp: {len(setB)} | Tarsi: {len(setC)} | All 3: {len(regions['A & B & C'])}")
            fig_venn = s1_build_venn_figure(state_full, state_col,
                                             (tissue_labels_s1["Antenna"], tissue_labels_s1["Palp"], tissue_labels_s1["Tarsi"]),
                                             counts)
            st.plotly_chart(fig_venn, use_container_width=False, key="s1_venn_Vm")
            st.markdown("**Genes in selected tissue combination:**")
            region_choice = st.selectbox("Tissue overlap region", region_options, key="s1_venn_region_Vm")
            region_key = region_key_map[region_choice]
            gene_set = regions.get(region_key, set())
            df_region = s1_make_region_table(state_code, gene_set, master_annot)
            if df_region.empty:
                st.info("No genes in this region for these thresholds.")
            else:
                st.dataframe(df_region)
                csv = df_region.to_csv(index=False).encode("utf-8")
                safe_state = state_full.replace(" ", "_")
                safe_region = region_choice.replace(" ", "_").replace("&", "and")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"BSF_{safe_state}_venn_{safe_region}.csv",
                                   mime="text/csv", key=f"s1_dl_Vm_{region_key}")

        with row1[1]:
            state_code = "VF"; partner = "Vm"
            state_full = COND_FULL.get(state_code, state_code)
            state_col = s1_get_full_label_color(state_full)
            st.markdown(f"### {state_full} (VF vs Vm)")
            sets_vf_vm = s1_get_up_sets_venn(by_tissue, state_code="VF", partner_code="Vm")
            setA = sets_vf_vm.get("Antenna", set()); setB = sets_vf_vm.get("Palp", set()); setC = sets_vf_vm.get("Tarsi", set())
            regions = s1_compute_venn_regions(setA, setB, setC)
            counts = {k: len(v) for k, v in regions.items()}
            st.caption(f"Antenna: {len(setA)} | Maxillary palp: {len(setB)} | Tarsi: {len(setC)} | All 3: {len(regions['A & B & C'])}")
            fig_venn = s1_build_venn_figure(f"{state_full} vs {COND_FULL.get(partner, partner)}", state_col,
                                             (tissue_labels_s1["Antenna"], tissue_labels_s1["Palp"], tissue_labels_s1["Tarsi"]),
                                             counts)
            st.plotly_chart(fig_venn, use_container_width=False, key="s1_venn_VF_Vm")
            st.markdown("**Genes in selected tissue combination:**")
            region_choice = st.selectbox("Tissue overlap region", region_options, key="s1_venn_region_VF_Vm")
            region_key = region_key_map[region_choice]
            gene_set = regions.get(region_key, set())
            df_region = s1_make_region_table(state_code, gene_set, master_annot)
            if df_region.empty:
                st.info("No genes in this region for these thresholds.")
            else:
                st.dataframe(df_region)
                csv = df_region.to_csv(index=False).encode("utf-8")
                safe_state = f"{state_full}_vs_{COND_FULL.get(partner, partner)}".replace(" ", "_")
                safe_region = region_choice.replace(" ", "_").replace("&", "and")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"BSF_{safe_state}_venn_{safe_region}.csv",
                                   mime="text/csv", key=f"s1_dl_VF_Vm_{region_key}")

        st.markdown("---")
        row2 = st.columns(2)
        with row2[0]:
            state_code = "VF"; partner = "MF"
            state_full = COND_FULL.get(state_code, state_code)
            state_col = s1_get_full_label_color(state_full)
            st.markdown(f"### {state_full} (MF vs VF)")
            sets_vf_mf = s1_get_up_sets_venn(by_tissue, state_code="VF", partner_code="MF")
            setA = sets_vf_mf.get("Antenna", set()); setB = sets_vf_mf.get("Palp", set()); setC = sets_vf_mf.get("Tarsi", set())
            regions = s1_compute_venn_regions(setA, setB, setC)
            counts = {k: len(v) for k, v in regions.items()}
            st.caption(f"Antenna: {len(setA)} | Maxillary palp: {len(setB)} | Tarsi: {len(setC)} | All 3: {len(regions['A & B & C'])}")
            fig_venn = s1_build_venn_figure(f"{state_full} vs {COND_FULL.get(partner, partner)}", state_col,
                                             (tissue_labels_s1["Antenna"], tissue_labels_s1["Palp"], tissue_labels_s1["Tarsi"]),
                                             counts)
            st.plotly_chart(fig_venn, use_container_width=False, key="s1_venn_VF_MF")
            st.markdown("**Genes in selected tissue combination:**")
            region_choice = st.selectbox("Tissue overlap region", region_options, key="s1_venn_region_VF_MF")
            region_key = region_key_map[region_choice]
            gene_set = regions.get(region_key, set())
            df_region = s1_make_region_table(state_code, gene_set, master_annot)
            if df_region.empty:
                st.info("No genes in this region for these thresholds.")
            else:
                st.dataframe(df_region)
                csv = df_region.to_csv(index=False).encode("utf-8")
                safe_state = f"{state_full}_vs_{COND_FULL.get(partner, partner)}".replace(" ", "_")
                safe_region = region_choice.replace(" ", "_").replace("&", "and")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"BSF_{safe_state}_venn_{safe_region}.csv",
                                   mime="text/csv", key=f"s1_dl_VF_MF_{region_key}")

        with row2[1]:
            state_code = "MF"
            state_full = COND_FULL.get(state_code, state_code)
            state_col = s1_get_full_label_color(state_full)
            st.markdown(f"### {state_full}")
            sets_mf = s1_get_up_sets_venn(by_tissue, state_code="MF", partner_code=None)
            setA = sets_mf.get("Antenna", set()); setB = sets_mf.get("Palp", set()); setC = sets_mf.get("Tarsi", set())
            regions = s1_compute_venn_regions(setA, setB, setC)
            counts = {k: len(v) for k, v in regions.items()}
            st.caption(f"Antenna: {len(setA)} | Maxillary palp: {len(setB)} | Tarsi: {len(setC)} | All 3: {len(regions['A & B & C'])}")
            fig_venn = s1_build_venn_figure(state_full, state_col,
                                             (tissue_labels_s1["Antenna"], tissue_labels_s1["Palp"], tissue_labels_s1["Tarsi"]),
                                             counts)
            st.plotly_chart(fig_venn, use_container_width=False, key="s1_venn_MF")
            st.markdown("**Genes in selected tissue combination:**")
            region_choice = st.selectbox("Tissue overlap region", region_options, key="s1_venn_region_MF")
            region_key = region_key_map[region_choice]
            gene_set = regions.get(region_key, set())
            df_region = s1_make_region_table(state_code, gene_set, master_annot)
            if df_region.empty:
                st.info("No genes in this region for these thresholds.")
            else:
                st.dataframe(df_region)
                csv = df_region.to_csv(index=False).encode("utf-8")
                safe_state = state_full.replace(" ", "_")
                safe_region = region_choice.replace(" ", "_").replace("&", "and")
                st.download_button("Download CSV", data=csv,
                                   file_name=f"BSF_{safe_state}_venn_{safe_region}.csv",
                                   mime="text/csv", key=f"s1_dl_MF_{region_key}")

    # -------------------- FIG 4 A-E TAB --------------------
    with tab_4de:
        if not tissue_cross:
            st.warning("No overlap data — ensure both contrasts are loaded for each tissue.")
        else:
            # ── Figure 4 A-C: Scatter plots ──────────────────────────────────
            st.subheader("Sex × mating LFC correlation — genes significant in both contrasts")
            st.caption("Manuscript: **Figure 4E**")
            st.markdown(
                "Each point is a gene significant in **both** the sex (VF vs Vm) "
                "and mating (MF vs VF) contrasts. Colored by GO domain (chemosensory genes dark red). "
                "Quadrants: Q1 (pink) = female-biased + mating-induced; Q3 (blue) = male-biased + mating-suppressed; "
                "Q4/Q2 (orange) = reversed."
            )
            scatter_cols = st.columns(3)
            for ci, tissue in enumerate(["Antenna", "Palp", "Tarsi"]):
                merged_sc, r_sc = s1_build_scatter_sex_mating(tissue_cross, tissue)
                with scatter_cols[ci]:
                    if merged_sc is None:
                        st.info(f"No overlap genes for {TISSUE_DISPLAY[tissue]}.")
                    else:
                        fig_sc = s1_fig_scatter(merged_sc, r_sc, TISSUE_DISPLAY[tissue])
                        st.plotly_chart(fig_sc, use_container_width=True,
                                        key=f"s1_scatter_{tissue}")

            st.markdown("---")
            # ── Figure 4D: Venns ────────────────────────────────────────────
            st.subheader("Overlap between sex-biased and mating-regulated gene sets")
            st.caption("Manuscript: **Figure 4D** (total Venn) · per-GO-domain Venns")
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

    # -------------------- GO NAMES TAB (S14–S16 equivalent) --------------------
    with tab_gonames:
        st.subheader("GO term enrichment browser")
        st.markdown(
            "Browse the top enriched GO terms for each gene set from the Venn diagrams. "
            "Supplementary Figures S14 (sex-biased only), S15 (mating-regulated only), S16 (overlap)."
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
                reg_tabs = st.tabs(["Sex-biased only (S14)", "Mating-regulated only (S15)", "Sex & mating overlap (S16)"])

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

    df["gene_fam"] = df[name_col].astype(str).apply(assign_family)
    df["Name_prefixed"] = df[name_col].astype(str).apply(h1_add_hill_prefix)
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
    st.caption("Manuscript: **Figs. S17–S23** (per gene family, A4 PDF — OR, GR, IR, OBP, PPK, CSP, TRP)")

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
        st.image(img_bytes, use_container_width=True)
    else:
        st.info(f"Could not render preview for {os.path.basename(file_path)}. "
                "Use the download button above.")


# =============================================================================
# === APP_TREES FUNCTIONS ===
# =============================================================================

TREE_FAMILY_MAP = {
    "All families": ("Figure_S1.pdf", None, None, None),
    "CSP (10 genes)":  ("Figure_S2.pdf", "S25_CSP_orthology.csv", "CSP_tree_clean.txt", "CSP_sequences_clean.fasta"),
    "GR (39 genes)":   ("Figure_S3.pdf", "S20_GR_orthology.csv",  "GR_tree_clean.txt",  "GR_sequences_clean.fasta"),
    "IR (101 genes)":  ("Figure_S4.pdf", "S21_IR_orthology.csv",  "IR_tree_clean.txt",  "IR_sequences_clean.fasta"),
    "OBP (75 genes)":  ("Figure_S5.pdf", "S22_OBP_orthology.csv", "OBP_tree_clean.txt", "OBP_sequences_clean.fasta"),
    "OR (192 genes)":  ("Figure_S6.pdf", "S19_OR_orthology.csv",  "OR_tree_clean.txt",  "OR_sequences_clean.fasta"),
    "PPK (39 genes)":  ("Figure_S7.pdf", "S23_PPK_orthology.csv", "PPK_tree_clean.txt", "PPK_sequences_clean.fasta"),
    "TRP (111 genes)": ("Figure_S8.pdf", "S24_TRP_orthology.csv", "TRP_tree_clean.txt", "TRP_sequences_clean.fasta"),
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
    "Figure 1 — Phylogenetic trees": (
        "**Figure 1.** Phylogenetic characterisation of *H. illucens* chemosensory gene families. "
        "ML trees of all seven families (OR 192 genes, IR 101, TRP 111, OBP 75, GR 39, PPK 39, CSP 10) "
        "inferred with IQ-TREE. Bootstrap support ≥ 70 shown on branches. "
        "Scale bar = amino acid substitutions per site."
    ),
    "Figure 2 — Transcriptome overview": (
        "**Figure 2.** Transcriptome-wide overview of chemosensory appendage gene expression. "
        "(A) PCA of all 27 RNA-seq libraries (log₂ normalised counts). "
        "(B) Pearson correlation matrix of 9 group means (convex hull areas: Antenna 2764.9, "
        "Palp 394.7, Tarsi 2043.9 PC-space units²). "
        "(C) Variance partitioning pie — all genes (n = 24,608; Ezekiel adj. R²). "
        "(E) Chemoreceptors-only PCA. "
        "(F) Chemoreceptors-only Pearson correlation. "
        "(G) Variance partitioning pie — chemosensory genes only."
    ),
    "Figure 3 — Appendage comparison": (
        "**Figure 3.** Appendage-specific transcriptional profiles. "
        "(A) Volcano plots of pairwise appendage DESeq2 contrasts (padj < 0.001, |log₂FC| ≥ 1). "
        "(B) Stacked GO domain % bars for appendage-biased DEGs. "
        "(C) GO names overlap per tissue. "
        "(D) Chemosensory gene classification pies per appendage. "
        "See also Figs. S9–S11 (GO enrichment bars) and Figs. S17–S23 (per-family heatmaps)."
    ),
    "Figure 4 — Sex & mating analysis": (
        "**Figure 4.** Sex- and mating-induced transcriptional regulation across appendages. "
        "(A) Volcano plots — sex contrast (VF vs Vm) per appendage. "
        "(B) Volcano plots — mating contrast (MF vs VF) per appendage. "
        "(C) Three-tissue Venn diagrams of upregulated genes per sex/mating state. "
        "(D) Two-circle Venn diagrams showing overlap between sex-biased and mating-regulated DEGs "
        "per appendage (padj < 0.001, |log₂FC| ≥ 1), with per-GO-domain breakdown. "
        "(E) Scatter plots of log₂FC(VF/Vm) vs log₂FC(MF/VF) for genes significant in both contrasts. "
        "Chemosensory genes highlighted (dark red). "
        "Note: tarsi *r* = −0.903 reflects a shared-VF contrast design artifact (see Limitations). "
        "See also Figs. S14–S16 (GO term names for sex/mating gene sets)."
    ),
}

SUPP_FIGURE_CAPTIONS = {
    "Figure S1 — All families (circular)": (
        "**Figure S1.** Circular ML tree of all seven *H. illucens* chemosensory gene families "
        "combined, with gene name labels. Bootstrap ≥ 70 shown. Families colour-coded."
    ),
    "Figure S2 — CSP": (
        "**Figure S2.** ML tree of the *H. illucens* Chemosensory Protein (CSP) family (10 genes). "
        "Note: CSP10 (→ *DmelCSP4*) is a cysteine string protein by naming convergence."
    ),
    "Figure S3 — GR": (
        "**Figure S3.** ML tree of the *H. illucens* Gustatory Receptor (GR) family (39 genes). "
        "GR1–5 form a CO₂ receptor clade; GR40 → *Gr43a* fructose receptor (BS = 99)."
    ),
    "Figure S4 — IR": (
        "**Figure S4.** ML tree of the *H. illucens* Ionotropic Receptor (IR) family (101 genes). "
        "IR7 = *Ir8a* ortholog (BS = 86); IR3 = *Ir25a* ortholog (BS = 95)."
    ),
    "Figure S5 — OBP": (
        "**Figure S5.** ML tree of the *H. illucens* Odorant Binding Protein (OBP) family (75 genes)."
    ),
    "Figure S6 — OR": (
        "**Figure S6.** ML tree of the *H. illucens* Olfactory Receptor (OR) family (192 genes). "
        "ORco and Or_putative_Orco-like both cluster in the Orco clade (BS = 98)."
    ),
    "Figure S7 — PPK": (
        "**Figure S7.** ML tree of the *H. illucens* Pickpocket (PPK) family (39 genes). "
        "PPK3 = *ppk23* ortholog (BS = 100); PPK18 = Nanchung co-ortholog (BS = 100)."
    ),
    "Figure S8 — TRP": (
        "**Figure S8.** ML tree of the *H. illucens* Transient Receptor Potential (TRP) family (111 genes). "
        "TRPN expansion: 44 genes (NompC clade), larger than TRPm (15 genes)."
    ),
    "Figure S9 — GO bars: Antenna": (
        "**Figure S9.** GO enrichment stacked bars for antenna-biased genes "
        "(genes significantly up in Antenna vs both Tarsi and Maxillary palp)."
    ),
    "Figure S10 — GO bars: Maxillary palp": (
        "**Figure S10.** GO enrichment stacked bars for maxillary palp-biased genes."
    ),
    "Figure S11 — GO bars: Tarsi": (
        "**Figure S11.** GO enrichment stacked bars for tarsi-biased genes."
    ),
    "Figure S12 — GO bars: Sex-biased genes": (
        "**Figure S12.** GO term enrichment of sex-biased genes (virgin female vs. virgin male) "
        "per chemosensory appendage. Top 20 GO terms per domain (MF, CC, BP) among significant "
        "DEGs (adjusted *p* < 0.001, |log₂FC| ≥ 1)."
    ),
    "Figure S13 — GO bars: Mating-regulated genes": (
        "**Figure S13.** GO term enrichment of mating-regulated genes (mated female vs. virgin "
        "female) per chemosensory appendage. Layout as in Supplementary Fig. S12."
    ),
    "Figure S14 — GO names: Sex-biased genes (only)": (
        "**Figure S14.** Top GO term names (per domain: MF, CC, BP) for genes with sex-biased expression "
        "(VF vs Vm) that are not also mating-regulated — per appendage (padj < 0.001, |log₂FC| ≥ 1). "
        "Interactive version: Sex & Mating Analysis → GO Term Browser → Sex-biased only."
    ),
    "Figure S15 — GO names: Mating-only genes": (
        "**Figure S15.** Top GO term names for genes differentially expressed in the mating contrast "
        "(MF vs VF) only — not in the sex contrast — per appendage. Layout as in Supplementary Fig. S14."
    ),
    "Figure S16 — GO names: Sex/mating overlap genes": (
        "**Figure S16.** Top GO term names for genes differentially expressed in **both** the sex "
        "(VF vs Vm) and mating (MF vs VF) contrasts within the same appendage. "
        "These represent genes jointly regulated by both sex and reproductive state."
    ),
    "Figure S17 — Heatmap: OR family": (
        "**Figure S17.** Expression heatmap of the Olfactory Receptor (OR) family (192 genes). "
        "Rows grouped by appendage of peak expression (Antenna → Maxillary palp → Tarsi), "
        "sorted by expression level within each group. Columns: Virgin male → Virgin female → Mated female "
        "per appendage. Colour scale: log₂(normalised counts + 1) centred at threshold of 10 counts; "
        "same scale across Figs. S17–S23."
    ),
    "Figure S18 — Heatmap: GR family": (
        "**Figure S18.** Expression heatmap of the Gustatory Receptor (GR) family (39 genes). "
        "Layout and scale as in Fig. S17."
    ),
    "Figure S19 — Heatmap: IR family": (
        "**Figure S19.** Expression heatmap of the Ionotropic Receptor (IR) family (101 genes). "
        "Layout and scale as in Fig. S17."
    ),
    "Figure S20 — Heatmap: OBP family": (
        "**Figure S20.** Expression heatmap of the Odorant Binding Protein (OBP) family (75 genes). "
        "Layout and scale as in Fig. S17."
    ),
    "Figure S21 — Heatmap: PPK family": (
        "**Figure S21.** Expression heatmap of the Pickpocket (PPK) family (39 genes). "
        "Layout and scale as in Fig. S17."
    ),
    "Figure S22 — Heatmap: CSP family": (
        "**Figure S22.** Expression heatmap of the Chemosensory Protein (CSP) family (10 genes). "
        "Layout and scale as in Fig. S17."
    ),
    "Figure S23 — Heatmap: TRP family": (
        "**Figure S23.** Expression heatmap of the Transient Receptor Potential channel (Trp) family (111 genes). "
        "Layout and colour scale as in Supplementary Fig. S17. Columns represent the nine group-mean normalised counts "
        "(Antenna: Vm, VF, MF; Palp: Vm, VF, MF; Tarsi: Vm, VF, MF). "
        "Rows sorted by maximum expression across all conditions (descending)."
    ),
}

_SUPP_FILES_ORDERED = [
    # AF 1–9: DESeq2 results per contrast
    ("S1_Appendage_DE_Antenna_vs_Tarsi.csv",
     "Additional file 1 — Appendage DE: Antenna vs Tarsi (DESeq2, 24,828 genes, 8 cols)", "text/csv"),
    ("S2_Appendage_DE_Antenna_vs_MaxillaryPalp.csv",
     "Additional file 2 — Appendage DE: Antenna vs Maxillary Palp", "text/csv"),
    ("S3_Appendage_DE_Tarsi_vs_MaxillaryPalp.csv",
     "Additional file 3 — Appendage DE: Tarsi vs Maxillary Palp", "text/csv"),
    ("S4_Antenna_sex_VirginFemale_vs_VirginMale.csv",
     "Additional file 4 — Antenna: sex effect VF vs Vm (14 cols, per-replicate counts)", "text/csv"),
    ("S5_Antenna_mating_MatedFemale_vs_VirginFemale.csv",
     "Additional file 5 — Antenna: mating effect MF vs VF", "text/csv"),
    ("S6_MaxillaryPalp_sex_VirginFemale_vs_VirginMale.csv",
     "Additional file 6 — Maxillary palp: sex effect VF vs Vm", "text/csv"),
    ("S7_MaxillaryPalp_mating_MatedFemale_vs_VirginFemale.csv",
     "Additional file 7 — Maxillary palp: mating effect MF vs VF", "text/csv"),
    ("S8_Tarsi_sex_VirginFemale_vs_VirginMale.csv",
     "Additional file 8 — Tarsi: sex effect VF vs Vm", "text/csv"),
    ("S9_Tarsi_mating_MatedFemale_vs_VirginFemale.csv",
     "Additional file 9 — Tarsi: mating effect MF vs VF", "text/csv"),
    # AF 10–11: Expression and annotation matrices
    ("S10_normalized_counts_all_samples.csv",
     "Additional file 10 — Normalized count matrix (24,828 genes × 27 samples + genomic coords)", "text/csv"),
    ("S11_GO_annotations_all_genes.csv",
     "Additional file 11 — GO annotation database (24,828 genes, 6 cols)", "text/csv"),
    # AF 12–18: Amino acid identity matrices
    ("S12_AA_identity_matrix_CSP.csv",
     "Additional file 12 — Pairwise AA identity matrix: Csp (10 × 10)", "text/csv"),
    ("S13_AA_identity_matrix_GR.csv",
     "Additional file 13 — Pairwise AA identity matrix: Gr (39 × 39)", "text/csv"),
    ("S14_AA_identity_matrix_IR.csv",
     "Additional file 14 — Pairwise AA identity matrix: Ir (101 × 101)", "text/csv"),
    ("S15_AA_identity_matrix_OBP.csv",
     "Additional file 15 — Pairwise AA identity matrix: Obp (75 × 75)", "text/csv"),
    ("S16_AA_identity_matrix_OR.csv",
     "Additional file 16 — Pairwise AA identity matrix: Or (192 × 192)", "text/csv"),
    ("S17_AA_identity_matrix_PPK.csv",
     "Additional file 17 — Pairwise AA identity matrix: Ppk (39 × 39)", "text/csv"),
    ("S18_AA_identity_matrix_TRP.csv",
     "Additional file 18 — Pairwise AA identity matrix: Trp (111 × 111)", "text/csv"),
    # AF 19: Diptera chemosensory comparison
    ("AF19_Diptera_chemosensory_comparison.csv",
     "Additional file 19 — Diptera chemosensory gene family comparison (6 species)", "text/csv"),
    # AF 20–26: Per-family phylogenetic orthology tables
    ("S19_OR_orthology.csv",
     "Additional file 20 — Or orthology classification (ML tree-based, 192 genes)", "text/csv"),
    ("S20_GR_orthology.csv",
     "Additional file 21 — Gr orthology classification (39 genes)", "text/csv"),
    ("S21_IR_orthology.csv",
     "Additional file 22 — Ir orthology classification (101 genes)", "text/csv"),
    ("S22_OBP_orthology.csv",
     "Additional file 23 — Obp orthology classification (75 genes)", "text/csv"),
    ("S23_PPK_orthology.csv",
     "Additional file 24 — Ppk orthology classification (39 genes)", "text/csv"),
    ("S24_TRP_orthology.csv",
     "Additional file 25 — Trp orthology classification (111 genes)", "text/csv"),
    ("S25_CSP_orthology.csv",
     "Additional file 26 — Csp orthology classification (10 genes)", "text/csv"),
    # AF 27–28: BLASTp and 1:1 ortholog identity
    ("AF27_Unknown_BLAST.csv",
     "Additional file 27 — BLASTp results for Unknown-domain DEGs vs D. melanogaster proteome (538 genes)", "text/csv"),
    ("AF28_AA_identity_1to1_orthologs.csv",
     "Additional file 28 — Pairwise AA identity between H. illucens 1:1 orthologs and nearest Dmel homolog", "text/csv"),
    # AF 29–33: Analysis scripts
    ("AF29_Figure_1.R",
     "Additional file 29 — R script: multivariate statistics (PCA, PERMANOVA, variance partitioning, correlation)", "text/plain"),
    ("AF30_Figure_2.R",
     "Additional file 30 — R script: multivariate transcriptome overview (identical to AF29)", "text/plain"),
    ("AF31_Figure_3.R",
     "Additional file 31 — R script: appendage DE analysis, GO enrichment, chemosensory classification (Fig. 3)", "text/plain"),
    ("AF32_Figure_4.R",
     "Additional file 32 — R script: sex- and mating-state DE, volcano plots, Venn diagrams, scatter (Fig. 4)", "text/plain"),
    ("AF33_generate_heatmaps.py",
     "Additional file 33 — Python script: per-family expression heatmaps (Figs. S17–S23)", "text/plain"),
]


def render_figures_tab():
    MAIN_FIG_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "figures" / "main")
    SUPP_FIG_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "figures" / "supplementary")
    SUPP_FILE_DIR = str(pathlib.Path(__file__).resolve().parent / "data" / "supplementary")

    st.header("Figures & Supplementary Data")

    fig_tab1, fig_tab2, fig_tab3 = st.tabs(
        ["Main Figures", "Supplementary Figures", "Download Additional Files (AF1–AF33)"]
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
        st.subheader("Additional Files (AF1–AF33)")
        st.caption(
            "All additional files from Perets *et al.* 2026 (v5.2). "
            "AF1–AF9: DESeq2 DE results per contrast. AF10: count matrix. AF11: GO database. "
            "AF12–AF18: amino acid identity matrices. AF19: Diptera comparison. "
            "AF20–AF26: per-family orthology tables. AF27: Unknown-domain BLASTp. "
            "AF28: 1:1 ortholog AA identity. AF29–AF33: R/Python analysis scripts."
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


# =============================================================================
# === MAIN APP ===
# =============================================================================

section = st.sidebar.radio(
    "Section",
    [
        "Phylogenetic Trees",
        "Transcriptome Overview",
        "Appendage Comparison",
        "Sex & Mating Analysis",
        "Chemosensory Heatmap",
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
