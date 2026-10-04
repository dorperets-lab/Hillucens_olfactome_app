################################################################################
# Figure 3 — appendage-specific transcriptional divergence
#
# WHAT THIS SCRIPT COVERS — main figure 3 and every supplementary figure that
# belongs to it. Nothing else is needed to reproduce them.
#
#   3A   pairwise volcanoes, now with the chemosensory genes labelled
#   3B   chemosensory versus rest, with a Fisher test  (REPLACES the GO bars)
#   3C   consensus appendage-biased Venn                (UNCHANGED)
#   3D   classification of the 510 non-redundant genes  (UNCHANGED)
#   S9   GO over-representation, antenna-biased genes
#   S10  GO over-representation, maxillary palp-biased genes
#   S11  GO over-representation, tarsus-biased genes
#   S12  GO domain composition of appendage-biased genes (descriptive)
#   S13-S19  per-family expression heatmaps (OR, GR, IR, OBP, PPK, CSP, TRP)
#
# HOW THIS FILE IS PUT TOGETHER — READ THIS BEFORE EDITING
#   The long middle section is the PUBLISHED analysis, with only its input and
#   output paths changed. That is deliberate: panels 3C and 3D are protected
#   numbers, and reimplementing them would risk moving values that are already
#   verified and quoted in the Results. Everything the revision changes is added
#   in clearly marked sections at the END of the file, which redraw 3A and 3B
#   from the same data objects and then build the supplementary figures.
#
#   So: to change a revision panel, edit the sections at the bottom. Do not
#   "tidy" the middle.
#
# WHAT CHANGED FROM THE PUBLISHED FILE
#   1. Inputs now come from 04_Additional_Files/ instead of the deleted
#      C:/Users/dorpe/Downloads/BSF2026/... path. The old input CSVs were
#      verified field by field against the Additional files before being
#      deleted; they differed only in NA versus empty cell and 0 versus 0.0.
#   2. The one input that is not an Additional file, BSF_Olfactory_ID.csv, now
#      lives in 05_R_Scripts/inputs/.
#   3. The count-ranked GO bar PDFs the old script wrote as S9-S11 are no longer
#      written. They are rebuilt properly at the end of this file.
#   4. Panels 3A and 3B are redrawn at the end, in their revised form.
#
# SIGN CONVENTION — VERIFIED, DO NOT "TIDY"
#   Additional files 1-3 store positive log2FC = higher in the FIRST-named
#   appendage and are NOT negated. Verified by correlating the stored values
#   against log2(mean ratios) recomputed from Additional file 10:
#   r = +0.996, +0.984, +0.995; ORco = +11.46 with antenna 113,064 vs tarsi 40.
#   Reversing this silently reports the antenna as chemosensory-DEPLETED.
#   (Additional files 4-9 and 29-31, used by Figure_4.R, go the other way.)
#
# OUTPUTS
#   02_Main_Figures/panels_new/            Figure_3A_*.pdf, Figure_3B_*.pdf,
#                                          Figure_3B_enrichment_stats.csv,
#                                          Figure_3_report.txt, classification CSV
#   03_Supplementary_Figures/panels_new/   Figure_S9/S10/S11.pdf,
#                                          Figure_S17 ... S23.pdf,
#                                          ORA_appendage_consensus.csv
#   Panels 3C and 3D are drawn to the plot device, as in the published script.
#
# PACKAGES  ggplot2, dplyr, tidyr, ggrepel, ggforce, gridExtra, cowplot, scales
#
# HOW TO RUN
#   setwd(".../Submission_v5_2026-09/05_R_Scripts")
#   source("Figure_3.R")
#
# CHECK THESE — printed at the end of the run
#   3A  Antenna 215/2372, Maxillary palp 35/1350 (antenna vs palp)
#   3B  antenna-biased chemosensory enrichment OR 6.26 vs palp, 5.33 vs tarsi;
#       palp and tarsi NOT enriched
#   3C  494 / 539 / 4161
#   3D  antenna 95 specific + 54 biased; palp 13 + 15; tarsi 23 + 75
#   S9  antenna MF olfactory receptor activity k = 127/278, OR 166.6, q = 1.17e-171
#   If a run disagrees, stop and say so rather than adjusting the figure.
################################################################################

################################################################################
# Figure 2 — single-file analysis script
#
# This file is the ONE script to source for everything behind Figure 2
# of Perets et al. 2026 (appendage DE, GO domain composition, Venn
# overlaps, chemosensory pie charts).
#
# Usage
# -----
#   1. Open R / RStudio.
#   2. Edit BASE_DIR below to point at the folder containing the
#      condition_vs_*, chemo_*.csv, and GO-annotated CSV files
#      (e.g. Supplemetary_Files/csv/ for this repository).
#   3. Source this file.  All plots are printed in order and a full
#      numerical report is written to Figure2_report.txt.
#
# Faithful reproduction of all plots from the Streamlit app
# (render_c1_tab / "Appendage Comparison" section).
#
# 5 sub-tabs reproduced:
#   1. Volcano (3 volcano plots)
#   2. GO domain % (3 stacked bar charts, with Chemosensory category)
#   3. GO Names overlap (Venn diagrams + horizontal bar charts per appendage)
#   4. Chemosensory pies (7 families x 3 appendages = 21 pie charts)
#   5. Chemosensory gene identity report (specific/biased gene lists)
#
# Convention: Antenna is always on the NEGATIVE x-axis side of volcanos
################################################################################

# ── Libraries ────────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(scales)
library(ggforce)   # for geom_circle in Venn
library(cowplot)   # for get_legend

# ── Paths — everything comes from this package ───────────────────────────────
PKG      <- "C:/Users/dorpe/OneDrive/Desktop/PhD_Projects/BSF/Updated_GDrive_version/Submission_v5_2026-09"
AF_DIR   <- file.path(PKG, "04_Additional_Files")
MAIN_OUT <- file.path(PKG, "02_Main_Figures", "panels_new")
SUPP_OUT <- file.path(PKG, "03_Supplementary_Figures", "panels_new")
for (d in c(MAIN_OUT, SUPP_OUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

BASE_DIR <- AF_DIR          # kept so downstream references still resolve
FIG_DIR  <- MAIN_OUT

ant_p_file    <- file.path(AF_DIR, "Additional_file_2_Appendage_DE_Antenna_vs_MaxillaryPalp.csv")
ant_leg_file  <- file.path(AF_DIR, "Additional_file_1_Appendage_DE_Antenna_vs_Tarsi.csv")
leg_p_file    <- file.path(AF_DIR, "Additional_file_3_Appendage_DE_Tarsi_vs_MaxillaryPalp.csv")
go_file       <- file.path(AF_DIR, "Additional_file_11_GO_annotations_all_genes.csv")
norm_file     <- file.path(AF_DIR, "Additional_file_10_normalized_counts_all_samples.csv")
chemo_id_file <- file.path(PKG, "05_R_Scripts", "inputs", "BSF_Olfactory_ID.csv")

report_file <- file.path(MAIN_OUT, "Figure_3_report.txt")

# ── Thresholds (matching app defaults) ───────────────────────────────────────────
PADJ_THR   <- 0.001
LFC_THR    <- 1.0
STRONG_LFC <- 2.0
EXPR_THR   <- 10.0
TOP_N_GO   <- 20

# ── Color palettes (exact hex from app) ──────────────────────────────────────────
GO_COLS <- c(Chemosensory = "#A60000", MF = "#3FC498", CC = "#4D50DB",
             BP = "#A860E3", Unknown = "#666666")
PAL_GO  <- c(GO_COLS, `Not Significant` = "#D9D9D9")

CHEMO_PIE_COLORS <- c(
  `Appendage-specific` = "#A60000",
  `Appendage-biased`   = "#1C1B8D",
  Expressed          = "#555555",
  `Not expressed`    = "#D6D6D6"
)

DOMAIN_FULL <- c(Chemosensory = "Chemosensory", MF = "Molecular function",
                 CC = "Cellular component", BP = "Biological process",
                 Unknown = "Unknown")

# ── Appendage prefixes (column name prefixes in normalized counts) ────────────────
TISSUE_PREFIXES <- list(
  Antenna          = c("Ant_"),
  `Maxillary palp` = c("P_", "Palp_"),
  Tarsi            = c("Leg_", "Tar_")
)

# =============================================================================
# REPORT — open file for writing
# =============================================================================
rpt <- file(report_file, open = "wt")
write_rpt <- function(...) {
  msg <- paste0(...)
  writeLines(msg, rpt)
}

write_rpt("================================================================")
write_rpt("BSF Appendage Comparison — Detailed Report")
write_rpt("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
write_rpt("================================================================")
write_rpt("")
write_rpt("Thresholds:")
write_rpt("  padj < ", PADJ_THR)
write_rpt("  |log2FoldChange| >= ", LFC_THR)
write_rpt("  Expression threshold >= ", EXPR_THR, " (normalized counts)")
write_rpt("  Top N GO names shown: ", TOP_N_GO)
write_rpt("")

# =============================================================================
# READ & PREPARE DATA
# =============================================================================

AntP   <- read.csv(ant_p_file,   stringsAsFactors = FALSE)
AntLeg <- read.csv(ant_leg_file, stringsAsFactors = FALSE)
LegP   <- read.csv(leg_p_file,   stringsAsFactors = FALSE)
go_raw <- read.csv(go_file,      stringsAsFactors = FALSE)
norm   <- read.csv(norm_file,    stringsAsFactors = FALSE)

write_rpt("Input files:")
write_rpt("  Ant vs Palp DE:   ", ant_p_file,   " (", nrow(AntP), " genes)")
write_rpt("  Ant vs Leg DE:    ", ant_leg_file,  " (", nrow(AntLeg), " genes)")
write_rpt("  Leg vs Palp DE:   ", leg_p_file,    " (", nrow(LegP), " genes)")
write_rpt("  GO annotations:   ", go_file,       " (", nrow(go_raw), " entries)")
write_rpt("  Norm counts:      ", norm_file,     " (", nrow(norm), " genes)")
write_rpt("")

# ── GO map (deduplicated, first per gene) ────────────────────────────────────────
clean_go_domain <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z]", "", x)
  mapd <- c(mf = "MF", molecularfunction = "MF",
            cc = "CC", cellularcomponent = "CC",
            bp = "BP", biologicalprocess = "BP")
  out <- mapd[x]
  out[is.na(out)] <- "Unknown"
  return(unname(out))
}

go_map <- data.frame(
  Gene      = trimws(as.character(go_raw[[1]])),
  GO_Name   = trimws(as.character(go_raw$GO_Name)),
  GO_Domain = clean_go_domain(go_raw$GO_Domain),
  stringsAsFactors = FALSE
)
go_map$GO_Name[is.na(go_map$GO_Name) | go_map$GO_Name == ""] <- "Unknown"
go_map <- go_map[!duplicated(go_map$Gene), ]

write_rpt("GO map: ", nrow(go_map), " unique genes with annotations")
write_rpt("  GO domain distribution:")
for (d in c("MF", "CC", "BP", "Unknown")) {
  write_rpt("    ", d, ": ", sum(go_map$GO_Domain == d))
}
write_rpt("")

# ── Helper: make JoinKey ─────────────────────────────────────────────────────────
make_joinkey <- function(df) trimws(as.character(df$Gene))

# ── Annotate DE tables with GO ───────────────────────────────────────────────────
annotate_with_go <- function(df) {
  df$JoinKey <- make_joinkey(df)
  merged <- merge(df, go_map, by.x = "JoinKey", by.y = "Gene", all.x = TRUE)
  merged$GO_Domain[is.na(merged$GO_Domain)] <- "Unknown"
  merged$GO_Name[is.na(merged$GO_Name) | merged$GO_Name == ""] <- "Unknown"
  return(merged)
}

# ── Chemo tag extraction ─────────────────────────────────────────────────────────
CHEMO_TAGS <- c("OR", "IR", "GR", "OBP", "CSP", "PPK", "ORCO", "TRP")

extract_chemo_tag <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  for (tag in CHEMO_TAGS) {
    pattern <- paste0("\\b", tag, "[0-9A-Za-z._-]*\\b")
    hits <- grepl(pattern, x, ignore.case = TRUE)
    out[hits & is.na(out)] <- tag
  }
  return(out)
}

# ── Prep volcano data frame ─────────────────────────────────────────────────────
# log2FoldChange is kept UNCHANGED (matching the app exactly), so that every
# downstream filter on LFC sign (Venn, GO %, up_keys, Direction) is
# semantically identical to c1_prep_volcano_df in app2.py.
#
# `plot_lfc` is a display-only mirror of log2FoldChange used by the volcano
# plot: when left_name == "Antenna" we negate it so Antenna ends up on the
# NEGATIVE x side, as you requested earlier. Nothing else reads plot_lfc.
#
# left_name / right_name attrs preserve app semantics:
#   LFC > 0  →  "<left_name> up"
#   LFC < 0  →  "<right_name> up"
# disp_left / disp_right attrs describe what sits where on the plot:
#   disp_left  = tissue shown on positive plot x
#   disp_right = tissue shown on negative plot x
prep_volcano <- function(df, left_name, right_name) {
  v <- annotate_with_go(df)
  v$padj <- as.numeric(v$padj); v$padj[is.na(v$padj)] <- 1.0
  v$log2FoldChange <- as.numeric(v$log2FoldChange); v$log2FoldChange[is.na(v$log2FoldChange)] <- 0.0

  # Display-only flip so Antenna ends up on the negative x axis.
  if (left_name == "Antenna") {
    v$plot_lfc <- -v$log2FoldChange
    disp_left  <- right_name    # positive x on plot
    disp_right <- left_name     # negative x on plot  (Antenna)
  } else {
    v$plot_lfc <- v$log2FoldChange
    disp_left  <- left_name
    disp_right <- right_name
  }

  v$padj_safe <- pmax(v$padj, .Machine$double.xmin)
  v$neg_log10 <- -log10(v$padj_safe)
  abs_lfc <- abs(v$log2FoldChange)
  v$is_sig <- (v$padj < PADJ_THR) & (abs_lfc >= LFC_THR)

  # Direction uses ORIGINAL log2FoldChange — matches app convention.
  v$Direction <- ifelse(v$is_sig & v$log2FoldChange > 0, paste0(left_name,  " up"),
                 ifelse(v$is_sig & v$log2FoldChange < 0, paste0(right_name, " up"),
                        "Not sig"))

  # Chemo tag
  name_col <- if ("Name" %in% names(v)) v$Name else v$JoinKey
  v$ChemoName <- extract_chemo_tag(name_col)
  v$is_chemo <- !is.na(v$ChemoName)

  # Color key: Chemosensory takes priority, then GO domain, then Not Significant
  v$ColorKey <- ifelse(!v$is_sig, "Not Significant",
                ifelse(v$is_chemo, "Chemosensory",
                       v$GO_Domain))

  attr(v, "left_name")  <- left_name     # semantic
  attr(v, "right_name") <- right_name    # semantic
  attr(v, "disp_left")  <- disp_left     # visual: positive x
  attr(v, "disp_right") <- disp_right    # visual: negative x
  return(v)
}

antp_v   <- prep_volcano(AntP,   "Antenna", "Maxillary palp")
antleg_v <- prep_volcano(AntLeg, "Antenna", "Tarsi")
legp_v   <- prep_volcano(LegP,   "Tarsi",   "Maxillary palp")

# ── Norm means per tissue ────────────────────────────────────────────────────────
norm$JoinKey <- make_joinkey(norm)

compute_tissue_means <- function(df) {
  out <- data.frame(JoinKey = df$JoinKey, stringsAsFactors = FALSE)
  for (tissue in names(TISSUE_PREFIXES)) {
    prefixes <- TISSUE_PREFIXES[[tissue]]
    cols <- names(df)[sapply(names(df), function(cn) any(startsWith(cn, prefixes)))]
    if (length(cols) > 0) {
      out[[paste0(tissue, "_mean")]] <- rowMeans(df[, cols, drop = FALSE], na.rm = TRUE)
    }
  }
  out <- out[!duplicated(out$JoinKey), ]
  return(out)
}

norm_means <- compute_tissue_means(norm)

# =============================================================================
# SUB-TAB 1: VOLCANO PLOTS (3 plots, one at a time)
# =============================================================================

write_rpt("================================================================")
write_rpt("1. VOLCANO PLOTS")
write_rpt("================================================================")
write_rpt("Convention: Antenna is always on the negative x-axis side.")
write_rpt("")

volcano_title <- function(vdf) {
  disp_left  <- attr(vdf, "disp_left")   # appendage on positive plot x
  disp_right <- attr(vdf, "disp_right")  # appendage on negative plot x
  # Count on plot coordinates, so the number next to each label is exactly
  # what the viewer sees on that side of the volcano.
  n_neg_side     <- sum(vdf$is_sig & vdf$plot_lfc < 0, na.rm = TRUE)
  n_pos_side     <- sum(vdf$is_sig & vdf$plot_lfc > 0, na.rm = TRUE)
  chemo_neg_side <- sum(vdf$is_sig & vdf$plot_lfc < 0 & vdf$is_chemo, na.rm = TRUE)
  chemo_pos_side <- sum(vdf$is_sig & vdf$plot_lfc > 0 & vdf$is_chemo, na.rm = TRUE)
  paste0(disp_right, " : ", chemo_neg_side, "/", n_neg_side,
         "    ",
         disp_left,  " : ", chemo_pos_side, "/", n_pos_side)
}

make_volcano_plot <- function(vdf) {
  vdf$ColorKey <- factor(vdf$ColorKey,
    levels = c("Not Significant", "MF", "CC", "BP", "Unknown", "Chemosensory"))
  vdf <- vdf[order(vdf$ColorKey), ]

  title_text <- volcano_title(vdf)

  p <- ggplot(vdf, aes(x = plot_lfc, y = neg_log10, color = ColorKey)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(
      values = PAL_GO,
      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown", "Not Significant"),
      labels = c("Chemosensory", "Molecular function", "Cellular component",
                 "Biological process", "Unknown", "Not Significant"),
      na.value = "#D9D9D9"
    ) +
    labs(
      title = title_text,
      x = expression(Log[2](FoldChange)),
      y = expression(-log[10](italic(padj))),
      color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  return(p)
}

# Report helper for volcano
report_volcano <- function(vdf, label) {
  disp_left  <- attr(vdf, "disp_left")
  disp_right <- attr(vdf, "disp_right")
  n_sig <- sum(vdf$is_sig)
  # Count by direction in the ORIGINAL LFC space, then relabel to visual sides.
  n_pos_side <- sum(vdf$is_sig & vdf$plot_lfc > 0, na.rm = TRUE)
  n_neg_side <- sum(vdf$is_sig & vdf$plot_lfc < 0, na.rm = TRUE)
  n_chemo_sig <- sum(vdf$is_sig & vdf$is_chemo)

  write_rpt("  --- ", label, " ---")
  write_rpt("  Comparison: ", disp_left, " (positive x) vs ", disp_right, " (negative x)")
  write_rpt("  Total significant DEGs (padj<", PADJ_THR, ", |LFC|>=", LFC_THR, "): ", n_sig)
  write_rpt("    ", disp_left,  " up (positive x): ", n_pos_side)
  write_rpt("    ", disp_right, " up (negative x): ", n_neg_side)
  write_rpt("    Chemosensory among significant: ", n_chemo_sig)
  write_rpt("  GO domain breakdown of significant DEGs:")
  if (n_sig > 0) {
    dom_tbl <- table(vdf$GO_Domain[vdf$is_sig])
    for (d in names(dom_tbl)) {
      write_rpt("    ", d, ": ", dom_tbl[d])
    }
  }
  write_rpt("")
}

# Plot and report each volcano individually
print(make_volcano_plot(antp_v))
report_volcano(antp_v, "Antenna vs Maxillary palp")

print(make_volcano_plot(antleg_v))
report_volcano(antleg_v, "Antenna vs Tarsi")

print(make_volcano_plot(legp_v))
report_volcano(legp_v, "Tarsi vs Maxillary palp")

cat("Done: Volcano plots\n")


# =============================================================================
# SUB-TAB 2: GO DOMAIN % (3 stacked bar charts, with Chemosensory category)
# =============================================================================

write_rpt("================================================================")
write_rpt("2. GO DOMAIN % STACKED BAR CHARTS")
write_rpt("================================================================")
write_rpt("Chemosensory genes are assigned their own category, taking priority")
write_rpt("over any GO domain assignment.")
write_rpt("")

make_go_percent_table <- function(vdf) {
  left  <- attr(vdf, "left_name")
  right <- attr(vdf, "right_name")

  d <- vdf[vdf$is_sig, , drop = FALSE]
  d$Dir <- ifelse(d$log2FoldChange > 0, paste0(left, " up"), paste0(right, " up"))

  # Chemosensory overrides GO domain
  d$PlotDomain <- ifelse(d$is_chemo, "Chemosensory", as.character(d$GO_Domain))
  d$PlotDomain <- factor(d$PlotDomain, levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))

  g <- d %>%
    group_by(Dir, PlotDomain, .drop = FALSE) %>%
    summarise(N = n(), .groups = "drop") %>%
    group_by(Dir) %>%
    mutate(Percent = 100 * N / sum(N)) %>%
    ungroup() %>%
    mutate(TextLabel = ifelse(Percent >= 2,
                              paste0("N=", N, "\n", sprintf("%.1f%%", Percent)),
                              ""))

  return(g)
}

make_go_bar <- function(vdf) {
  tbl <- make_go_percent_table(vdf)
  title_text <- paste0("GO composition - ", volcano_title(vdf))

  p <- ggplot(tbl, aes(x = Dir, y = Percent, fill = PlotDomain)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    geom_text(aes(label = TextLabel), position = position_stack(vjust = 0.5),
              size = 2.5, color = "white", fontface = "bold", lineheight = 0.85) +
    scale_fill_manual(
      values = GO_COLS,
      labels = DOMAIN_FULL,
      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown")
    ) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(title = title_text, x = NULL, y = "% of significant DEGs", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      legend.position = "bottom",
      axis.text.x = element_text(size = 8, angle = 15, hjust = 1)
    )
  return(p)
}

# Report helper for GO bars
report_go_bar <- function(vdf, label) {
  tbl <- make_go_percent_table(vdf)
  write_rpt("  --- ", label, " ---")
  for (dir_val in unique(tbl$Dir)) {
    sub <- tbl[tbl$Dir == dir_val, ]
    total_n <- sum(sub$N)
    write_rpt("  Direction: ", dir_val, " (total N=", total_n, ")")
    for (i in seq_len(nrow(sub))) {
      write_rpt("    ", as.character(sub$PlotDomain[i]), ": N=", sub$N[i],
                " (", sprintf("%.1f", sub$Percent[i]), "%)")
    }
  }
  write_rpt("")
}

# Plot and report each GO bar individually
print(make_go_bar(antp_v))
report_go_bar(antp_v, "Antenna vs Maxillary palp")

print(make_go_bar(antleg_v))
report_go_bar(antleg_v, "Antenna vs Tarsi")

print(make_go_bar(legp_v))
report_go_bar(legp_v, "Tarsi vs Maxillary palp")

cat("Done: GO domain % bar charts\n")


# =============================================================================
# SUB-TAB 3: GO NAMES OVERLAP -- Venn diagrams + horizontal bar charts
# =============================================================================

write_rpt("================================================================")
write_rpt("3. GO NAMES OVERLAP (Venn diagrams + horizontal bar charts)")
write_rpt("================================================================")
write_rpt("Overlap = genes significantly up in both relevant contrasts for an appendage.")
write_rpt("  Antenna:  up in AntP AND up in AntLeg")
write_rpt("  Palp:     down in AntP AND down in LegP")
write_rpt("  Tarsi:    down in AntLeg AND up in LegP")
write_rpt("")

# ── Helper: get significant gene keys ────────────────────────────────────────────
up_keys <- function(df, sign) {
  df$padj <- as.numeric(df$padj); df$padj[is.na(df$padj)] <- 1.0
  df$log2FoldChange <- as.numeric(df$log2FoldChange); df$log2FoldChange[is.na(df$log2FoldChange)] <- 0.0
  df$JoinKey <- make_joinkey(df)
  if (sign == 1) {
    sub <- df[df$padj < PADJ_THR & df$log2FoldChange > LFC_THR, ]
  } else {
    sub <- df[df$padj < PADJ_THR & df$log2FoldChange < -LFC_THR, ]
  }
  return(unique(trimws(sub$JoinKey)))
}

ant_keys   <- intersect(up_keys(AntP, +1),  up_keys(AntLeg, +1))
palp_keys  <- intersect(up_keys(AntP, -1),  up_keys(LegP, -1))
tarsi_keys <- intersect(up_keys(AntLeg, -1), up_keys(LegP, +1))

write_rpt("Overlap gene counts:")
write_rpt("  Antenna overlap:  ", length(ant_keys), " genes")
write_rpt("  Palp overlap:     ", length(palp_keys), " genes")
write_rpt("  Tarsi overlap:    ", length(tarsi_keys), " genes")
write_rpt("")

# ── Domain-specific keys for Venn ────────────────────────────────────────────────
domain_keys_for_tissue <- function(tissue, domain, antp_v, antleg_v, legp_v) {
  if (tissue == "Antenna") {
    s1 <- antp_v$JoinKey[antp_v$is_sig & antp_v$log2FoldChange > 0 & antp_v$GO_Domain == domain]
    s2 <- antleg_v$JoinKey[antleg_v$is_sig & antleg_v$log2FoldChange > 0 & antleg_v$GO_Domain == domain]
    lbl_left <- "Maxillary palp"; lbl_right <- "Tarsi"
  } else if (tissue == "Maxillary palp") {
    s1 <- antp_v$JoinKey[antp_v$is_sig & antp_v$log2FoldChange < 0 & antp_v$GO_Domain == domain]
    s2 <- legp_v$JoinKey[legp_v$is_sig & legp_v$log2FoldChange < 0 & legp_v$GO_Domain == domain]
    lbl_left <- "Antenna"; lbl_right <- "Tarsi"
  } else {
    s1 <- antleg_v$JoinKey[antleg_v$is_sig & antleg_v$log2FoldChange < 0 & antleg_v$GO_Domain == domain]
    s2 <- legp_v$JoinKey[legp_v$is_sig & legp_v$log2FoldChange > 0 & legp_v$GO_Domain == domain]
    lbl_left <- "Maxillary palp"; lbl_right <- "Antenna"
  }
  s1 <- unique(as.character(s1)); s2 <- unique(as.character(s2))
  overlap <- intersect(s1, s2)
  list(lbl_left = lbl_left, lbl_right = lbl_right,
       s1 = s1, s2 = s2, overlap = overlap,
       n1 = length(s1), n2 = length(s2), n_overlap = length(overlap))
}

# ── Draw Venn ────────────────────────────────────────────────────────────────────
draw_venn_gg <- function(lbl_left, lbl_right, n_left, n_right, n_overlap,
                         title, domain_color) {
  circles <- data.frame(x0 = c(-0.5, 0.5), y0 = c(0, 0), r = c(1, 1))
  labels <- data.frame(
    x   = c(-0.85, 0.85, 0, -0.5, 0.5),
    y   = c(0, 0, 0.15, -0.9, -0.9),
    txt = c(as.character(n_left), as.character(n_right),
            as.character(n_overlap), lbl_left, lbl_right),
    sz  = c(4, 4, 5, 3.5, 3.5)
  )
  fill_col <- adjustcolor(domain_color, alpha.f = 0.35)

  p <- ggplot() +
    ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r),
                         fill = fill_col, color = NA, linewidth = 0) +
    geom_text(data = labels, aes(x = x, y = y, label = txt), size = labels$sz) +
    coord_fixed(xlim = c(-2, 2), ylim = c(-1.5, 1.5)) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
  return(p)
}

# ── GO name horizontal bar chart ─────────────────────────────────────────────────
go_name_overlap_table <- function(keys, domain, drop_unknown = FALSE) {
  gm <- go_map[go_map$Gene %in% keys & go_map$GO_Domain == domain, ]
  if (drop_unknown) gm <- gm[gm$GO_Name != "Unknown", ]
  if (nrow(gm) == 0) return(data.frame(GO_Name = character(0), N = integer(0)))
  tbl <- gm %>%
    group_by(GO_Name) %>%
    summarise(N = n(), .groups = "drop") %>%
    arrange(desc(N), GO_Name) %>%
    head(TOP_N_GO)
  return(tbl)
}

make_go_name_bar <- function(tbl, domain, x_max = NULL) {
  if (nrow(tbl) == 0) return(ggplot() + theme_void() + labs(title = paste("No", domain, "genes")))
  if (is.null(x_max)) x_max <- max(tbl$N)
  tbl$GO_Name <- factor(tbl$GO_Name, levels = rev(tbl$GO_Name))
  label <- tolower(DOMAIN_FULL[domain])

  p <- ggplot(tbl, aes(x = N, y = GO_Name)) +
    geom_bar(stat = "identity", fill = GO_COLS[domain], width = 0.7) +
    scale_x_continuous(limits = c(0, max(1, x_max * 1.1))) +
    labs(title = paste0("Overlap GO Names in ", label), x = "Gene count", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          axis.text.y = element_text(size = 7))
  return(p)
}

# ── Generate all Venn + bar plots per tissue (one at a time) ─────────────────────
tissues_overlap <- list(
  list(name = "Antenna",          keys = ant_keys),
  list(name = "Maxillary palp",   keys = palp_keys),
  list(name = "Tarsi",            keys = tarsi_keys)
)

domains_for_venn <- c("MF", "CC", "BP")

for (tinfo in tissues_overlap) {
  tissue <- tinfo$name
  keys   <- tinfo$keys

  write_rpt("  --- ", tissue, " ---")

  # Venn diagrams — one per domain
  for (dom in domains_for_venn) {
    dk <- domain_keys_for_tissue(tissue, dom, antp_v, antleg_v, legp_v)
    vp <- draw_venn_gg(dk$lbl_left, dk$lbl_right,
                       dk$n1, dk$n2, dk$n_overlap,
                       title = paste0(tissue, " - ", DOMAIN_FULL[dom]),
                       domain_color = GO_COLS[dom])
    print(vp)

    write_rpt("  Venn ", dom, ": ", dk$lbl_left, "=", dk$n1,
              ", ", dk$lbl_right, "=", dk$n2, ", overlap=", dk$n_overlap)
  }

  # GO name bars — one per domain
  t_mf <- go_name_overlap_table(keys, "MF")
  t_cc <- go_name_overlap_table(keys, "CC")
  t_bp <- go_name_overlap_table(keys, "BP")
  x_max <- max(c(if (nrow(t_mf) > 0) max(t_mf$N) else 0,
                  if (nrow(t_cc) > 0) max(t_cc$N) else 0,
                  if (nrow(t_bp) > 0) max(t_bp$N) else 0))

  print(make_go_name_bar(t_mf, "MF", x_max))
  print(make_go_name_bar(t_cc, "CC", x_max))
  print(make_go_name_bar(t_bp, "BP", x_max))

  # Report top GO names
  for (dom_info in list(list("MF", t_mf), list("CC", t_cc), list("BP", t_bp))) {
    dom <- dom_info[[1]]; tbl <- dom_info[[2]]
    write_rpt("  Top GO Names (", dom, ", N=", nrow(tbl), "):")
    if (nrow(tbl) > 0) {
      for (i in seq_len(min(nrow(tbl), TOP_N_GO))) {
        write_rpt("    ", tbl$GO_Name[i], ": ", tbl$N[i])
      }
    }
  }
  write_rpt("")
}
cat("Done: GO Names overlap (Venns + bars)\n")

# ── Save supplementary GO bar PDFs (S1=Antenna, S2=Palp, S3=Tarsi) ──────────
sup_fig_labels <- c(Antenna = "S1", `Maxillary palp` = "S2", Tarsi = "S3")

for (tinfo in tissues_overlap) {
  tissue <- tinfo$name
  keys   <- tinfo$keys
  fig_label <- sup_fig_labels[tissue]

  t_mf <- go_name_overlap_table(keys, "MF")
  t_cc <- go_name_overlap_table(keys, "CC")
  t_bp <- go_name_overlap_table(keys, "BP")
  x_max <- max(c(if (nrow(t_mf) > 0) max(t_mf$N) else 0,
                  if (nrow(t_cc) > 0) max(t_cc$N) else 0,
                  if (nrow(t_bp) > 0) max(t_bp$N) else 0))

  p_mf <- make_go_name_bar(t_mf, "MF", x_max) +
    labs(title = paste0(tissue, " — Molecular function"))
  p_cc <- make_go_name_bar(t_cc, "CC", x_max) +
    labs(title = paste0(tissue, " — Cellular component"))
  p_bp <- make_go_name_bar(t_bp, "BP", x_max) +
    labs(title = paste0(tissue, " — Biological process"))

  combined <- gridExtra::arrangeGrob(
    p_mf, p_cc, p_bp,
    ncol = 1,
    top = grid::textGrob(
      paste0("Figure ", fig_label, ". GO term enrichment of ",
             tissue, "-biased genes"),
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    )
  )

  # SUPERSEDED. These count-ranked GO bars were the submitted S9-S11, and they
  # are exactly what Phi objected to (C20/C24/C61): a long bar meant a large GO
  # term, not an enriched one. They are rebuilt as a proper over-representation
  # analysis in the S9-S11 section at the end of this script, so nothing is
  # written here any more. `combined` is still built so the on-screen report and
  # the overlap tables above are unchanged.
  cat("  (count-ranked GO bars for ", tissue,
      " not written - superseded by the S9-S11 ORA section)\n", sep = "")
}


# =============================================================================
# SUB-TAB 4: CHEMOSENSORY PIES (7 families x 3 tissues = up to 21 pies)
# =============================================================================

write_rpt("================================================================")
write_rpt("4. CHEMOSENSORY PIE CHARTS")
write_rpt("================================================================")
write_rpt("Classification per appendage per gene family:")
write_rpt("  Appendage-specific: expressed only in that appendage (>= ", EXPR_THR, ")")
write_rpt("  Appendage-biased:   highest expression + significantly DE vs others")
write_rpt("  Expressed:       expressed but not specific or biased")
write_rpt("  Not expressed:   below threshold in that appendage")
write_rpt("")

# ── Read chemosensory gene ID table ─────────────────────────────────────────────
chemo_id <- read.csv(chemo_id_file, stringsAsFactors = FALSE)
# Clean BOM if present
names(chemo_id) <- gsub("^\uFEFF", "", names(chemo_id))
chemo_id$Name <- trimws(as.character(chemo_id$Name))
chemo_id$Transcript <- trimws(as.character(chemo_id$Transcript))

# Extract family from gene name (e.g. "OR123" -> "Or", "ORCO" -> "Or")
extract_family <- function(name) {
  name_upper <- toupper(name)
  ifelse(grepl("^ORCO", name_upper), "Or",
  ifelse(grepl("^OR",   name_upper), "Or",
  ifelse(grepl("^GR",   name_upper), "Gr",
  ifelse(grepl("^IR",   name_upper), "Ir",
  ifelse(grepl("^OBP",  name_upper), "Obp",
  ifelse(grepl("^CSP",  name_upper), "Csp",
  ifelse(grepl("^PPK",  name_upper), "Ppk",
  ifelse(grepl("^TRP",  name_upper), "Trp",
         NA_character_))))))))
}
chemo_id$Family <- extract_family(chemo_id$Name)
chemo_id <- chemo_id[!is.na(chemo_id$Family), ]

write_rpt("  Chemosensory ID file: ", nrow(chemo_id), " genes")
write_rpt("  Family distribution:")
for (f in c("Or","Gr","Ir","Obp","Csp","Ppk","Trp")) {
  write_rpt("    ", f, ": ", sum(chemo_id$Family == f))
}
write_rpt("")

# ── Map chemosensory names to transcript IDs in norm counts ─────────────────────
# norm$Gene contains transcript IDs (e.g. "XM_038060289.1")
# chemo_id$Transcript has "rna-" prefix (e.g. "rna-XM_038060289.1") — strip it
chemo_id$NormKey <- sub("^rna-", "", chemo_id$Transcript)

# ── Deduplicate chemosensory genes using AA similarity matrices ─────────────────
# Genes with 100% AA identity are paralogs/duplications and should be counted as
# one representative. We read each family's similarity matrix, find transitive
# clusters of 100%-identical genes, and keep one representative per cluster.

sim_matrix_files <- list(
  Or  = file.path(AF_DIR, "Additional_file_16_AA_identity_matrix_OR.csv"),
  Gr  = file.path(AF_DIR, "Additional_file_13_AA_identity_matrix_GR.csv"),
  Ir  = file.path(AF_DIR, "Additional_file_14_AA_identity_matrix_IR.csv"),
  Obp = file.path(AF_DIR, "Additional_file_15_AA_identity_matrix_OBP.csv"),
  Csp = file.path(AF_DIR, "Additional_file_12_AA_identity_matrix_CSP.csv"),
  Ppk = file.path(AF_DIR, "Additional_file_17_AA_identity_matrix_PPK.csv"),
  Trp = file.path(AF_DIR, "Additional_file_18_AA_identity_matrix_TRP.csv")
)

find_duplicate_clusters <- function(sim_file) {
  if (!file.exists(sim_file)) return(list())
  mat <- read.csv(sim_file, row.names = 1, check.names = FALSE)
  genes <- rownames(mat)
  # Find all pairs with >= 100% identity
  pairs <- list()
  for (i in seq_along(genes)) {
    for (j in seq_along(genes)) {
      if (i >= j) next
      val <- suppressWarnings(as.numeric(mat[i, j]))
      if (!is.na(val) && val >= 100.0) {
        pairs[[length(pairs) + 1]] <- c(genes[i], genes[j])
      }
    }
  }
  if (length(pairs) == 0) return(list())
  # Build adjacency and find connected components (transitive clusters)
  adj <- list()
  for (p in pairs) {
    adj[[p[1]]] <- unique(c(adj[[p[1]]], p[2]))
    adj[[p[2]]] <- unique(c(adj[[p[2]]], p[1]))
  }
  visited <- character(0)
  clusters <- list()
  for (node in names(adj)) {
    if (node %in% visited) next
    cluster <- character(0)
    stack <- node
    while (length(stack) > 0) {
      n <- stack[1]; stack <- stack[-1]
      if (n %in% visited) next
      visited <- c(visited, n)
      cluster <- c(cluster, n)
      neighbors <- adj[[n]]
      if (!is.null(neighbors)) stack <- c(stack, setdiff(neighbors, visited))
    }
    clusters[[length(clusters) + 1]] <- sort(cluster)
  }
  return(clusters)
}

# Build a lookup: gene name -> representative (first alphabetically in cluster)
dup_rep_map <- list()  # gene_name -> representative_name
dup_cluster_members <- list()  # gene_name -> all OTHER members in its 100% cluster
dup_clusters_all <- list()
for (fam in names(sim_matrix_files)) {
  clusters <- find_duplicate_clusters(sim_matrix_files[[fam]])
  dup_clusters_all[[fam]] <- clusters
  for (cl in clusters) {
    rep_name <- cl[1]  # first alphabetically (already sorted)
    for (gene in cl) {
      dup_rep_map[[toupper(gene)]] <- toupper(rep_name)
      dup_cluster_members[[toupper(gene)]] <- setdiff(cl, gene)
    }
  }
}

# Mark representatives in chemo_id — keep only one per cluster
chemo_id$RepName <- toupper(chemo_id$Name)
for (i in seq_len(nrow(chemo_id))) {
  key <- toupper(chemo_id$Name[i])
  if (!is.null(dup_rep_map[[key]])) {
    chemo_id$RepName[i] <- dup_rep_map[[key]]
  }
}
chemo_id$is_rep <- (toupper(chemo_id$Name) == chemo_id$RepName)

n_before <- nrow(chemo_id)
chemo_id_dedup <- chemo_id[chemo_id$is_rep, ]
n_after <- nrow(chemo_id_dedup)

write_rpt("  AA similarity deduplication (100% identity clusters):")
write_rpt("    Before: ", n_before, " genes")
write_rpt("    After:  ", n_after, " genes (removed ", n_before - n_after, " duplicates)")
for (fam in names(dup_clusters_all)) {
  cls <- dup_clusters_all[[fam]]
  if (length(cls) > 0) {
    write_rpt("    ", toupper(fam), ": ", length(cls), " cluster(s)")
    for (cl in cls) {
      write_rpt("      ", paste(cl, collapse = ", "), " -> keep ", cl[1])
    }
  }
}
write_rpt("")

# Use deduplicated set for classification
chemo_id <- chemo_id_dedup

# ── Compute appendage means for chemosensory genes ──────────────────────────────
ant_cols  <- grep("^Ant_",  names(norm), value = TRUE)
palp_cols <- grep("^P_",    names(norm), value = TRUE)
tarsi_cols <- grep("^Leg_", names(norm), value = TRUE)

norm$Ant_mean  <- rowMeans(norm[, ant_cols,   drop = FALSE], na.rm = TRUE)
norm$Palp_mean <- rowMeans(norm[, palp_cols,  drop = FALSE], na.rm = TRUE)
norm$Tarsi_mean <- rowMeans(norm[, tarsi_cols, drop = FALSE], na.rm = TRUE)

# ── Classify each chemosensory gene per appendage ───────────────────────────────
# Appendage-specific: expressed (>= EXPR_THR) ONLY in that appendage
# Appendage-biased: expressed in that appendage AND significantly DE (higher)
#                   vs BOTH other appendages (padj < PADJ_THR, |LFC| >= LFC_THR)
# Expressed: expressed (>= EXPR_THR) in that appendage but not specific/biased
# Not expressed: mean < EXPR_THR in that appendage

classify_chemo_genes <- function(chemo_id, norm, AntP, AntLeg, LegP) {
  # Build lookup: Gene -> row in norm
  norm_lookup <- setNames(seq_len(nrow(norm)), trimws(as.character(norm$Gene)))

  # Build DE lookup tables (Gene -> padj, LFC)
  make_de_lookup <- function(de_df) {
    de_df$Gene <- trimws(as.character(de_df$Gene))
    list(
      padj = setNames(as.numeric(de_df$padj), de_df$Gene),
      lfc  = setNames(as.numeric(de_df$log2FoldChange), de_df$Gene)
    )
  }
  de_ant_p   <- make_de_lookup(AntP)    # LFC > 0 = Ant up; LFC < 0 = Palp up
  de_ant_leg <- make_de_lookup(AntLeg)  # LFC > 0 = Ant up; LFC < 0 = Leg up
  de_leg_p   <- make_de_lookup(LegP)    # LFC > 0 = Leg up; LFC < 0 = Palp up

  results <- data.frame(
    Name = chemo_id$Name,
    Family = chemo_id$Family,
    Transcript = chemo_id$Transcript,
    Ant_mean = NA_real_, Palp_mean = NA_real_, Tarsi_mean = NA_real_,
    Antenna_Class = NA_character_, Palp_Class = NA_character_, Tarsi_Class = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(chemo_id))) {
    # Try matching by NormKey (transcript without rna- prefix) first, then by Name
    tkey <- chemo_id$NormKey[i]
    nkey <- chemo_id$Name[i]
    idx <- norm_lookup[tkey]
    if (is.na(idx)) idx <- norm_lookup[nkey]
    # Also try matching via norm$Name column
    if (is.na(idx) && "Name" %in% names(norm)) {
      nm_match <- which(trimws(as.character(norm$Name)) == nkey)
      if (length(nm_match) > 0) idx <- nm_match[1]
    }

    if (is.na(idx)) {
      results$Antenna_Class[i] <- "Not expressed"
      results$Palp_Class[i]    <- "Not expressed"
      results$Tarsi_Class[i]   <- "Not expressed"
      next
    }

    ant_m  <- norm$Ant_mean[idx]
    palp_m <- norm$Palp_mean[idx]
    tar_m  <- norm$Tarsi_mean[idx]
    results$Ant_mean[i]  <- ant_m
    results$Palp_mean[i] <- palp_m
    results$Tarsi_mean[i] <- tar_m

    expressed_ant  <- ant_m  >= EXPR_THR
    expressed_palp <- palp_m >= EXPR_THR
    expressed_tar  <- tar_m  >= EXPR_THR

    # Get DE stats using the Gene (transcript) key
    gene_key <- trimws(as.character(norm$Gene[idx]))

    # Antenna classification
    if (!expressed_ant) {
      results$Antenna_Class[i] <- "Not expressed"
    } else if (expressed_ant & !expressed_palp & !expressed_tar) {
      results$Antenna_Class[i] <- "Appendage-specific"
    } else {
      # Check if biased: Ant significantly higher than BOTH others
      p_vs_palp <- de_ant_p$padj[gene_key]; lfc_vs_palp <- de_ant_p$lfc[gene_key]
      p_vs_tar  <- de_ant_leg$padj[gene_key]; lfc_vs_tar <- de_ant_leg$lfc[gene_key]
      biased <- FALSE
      if (!is.na(p_vs_palp) && !is.na(lfc_vs_palp) &&
          !is.na(p_vs_tar)  && !is.na(lfc_vs_tar)) {
        biased <- (p_vs_palp < PADJ_THR & lfc_vs_palp >= LFC_THR &
                   p_vs_tar  < PADJ_THR & lfc_vs_tar  >= LFC_THR)
      }
      results$Antenna_Class[i] <- ifelse(biased, "Appendage-biased", "Expressed")
    }

    # Maxillary palp classification
    if (!expressed_palp) {
      results$Palp_Class[i] <- "Not expressed"
    } else if (expressed_palp & !expressed_ant & !expressed_tar) {
      results$Palp_Class[i] <- "Appendage-specific"
    } else {
      # Palp higher than Ant: AntP LFC < 0 (Palp up)
      # Palp higher than Tar: LegP LFC < 0 (Palp up)
      p_vs_ant <- de_ant_p$padj[gene_key]; lfc_vs_ant <- de_ant_p$lfc[gene_key]
      p_vs_tar <- de_leg_p$padj[gene_key]; lfc_vs_tar <- de_leg_p$lfc[gene_key]
      biased <- FALSE
      if (!is.na(p_vs_ant) && !is.na(lfc_vs_ant) &&
          !is.na(p_vs_tar) && !is.na(lfc_vs_tar)) {
        biased <- (p_vs_ant < PADJ_THR & lfc_vs_ant <= -LFC_THR &
                   p_vs_tar < PADJ_THR & lfc_vs_tar <= -LFC_THR)
      }
      results$Palp_Class[i] <- ifelse(biased, "Appendage-biased", "Expressed")
    }

    # Tarsi classification
    if (!expressed_tar) {
      results$Tarsi_Class[i] <- "Not expressed"
    } else if (expressed_tar & !expressed_ant & !expressed_palp) {
      results$Tarsi_Class[i] <- "Appendage-specific"
    } else {
      # Tar higher than Ant: AntLeg LFC < 0 (Leg up)
      # Tar higher than Palp: LegP LFC > 0 (Leg up)
      p_vs_ant  <- de_ant_leg$padj[gene_key]; lfc_vs_ant  <- de_ant_leg$lfc[gene_key]
      p_vs_palp <- de_leg_p$padj[gene_key];   lfc_vs_palp <- de_leg_p$lfc[gene_key]
      biased <- FALSE
      if (!is.na(p_vs_ant)  && !is.na(lfc_vs_ant) &&
          !is.na(p_vs_palp) && !is.na(lfc_vs_palp)) {
        biased <- (p_vs_ant  < PADJ_THR & lfc_vs_ant  <= -LFC_THR &
                   p_vs_palp < PADJ_THR & lfc_vs_palp >= LFC_THR)
      }
      results$Tarsi_Class[i] <- ifelse(biased, "Appendage-biased", "Expressed")
    }
  }
  return(results)
}

chemo_classified <- classify_chemo_genes(chemo_id, norm, AntP, AntLeg, LegP)
write_rpt("  Classification complete: ", nrow(chemo_classified), " chemosensory genes classified")
write_rpt("")

fam_order <- c("Or", "Gr", "Ir", "Obp", "Csp", "Ppk", "Trp")
tissue_keys <- c("Antenna", "Palp", "Tarsi")
tissue_labels_map <- c(Antenna = "Antenna", Palp = "Maxillary palp", Tarsi = "Tarsi")
class_cols_map <- c(Antenna = "Antenna_Class", Palp = "Palp_Class", Tarsi = "Tarsi_Class")

# ── Build chemo summary from classified data ─────────────────────────────────────
chemo_summary <- list()
for (fam in fam_order) {
  fam_df <- chemo_classified[chemo_classified$Family == fam, ]
  if (nrow(fam_df) == 0) next
  fam_counts <- list()
  for (tkey in tissue_keys) {
    col <- class_cols_map[tkey]
    vc <- table(fam_df[[col]])
    fam_counts[[tkey]] <- c(
      `Appendage-specific` = as.integer(ifelse("Appendage-specific" %in% names(vc), vc["Appendage-specific"], 0)),
      `Appendage-biased`   = as.integer(ifelse("Appendage-biased" %in% names(vc), vc["Appendage-biased"], 0)),
      Expressed            = as.integer(ifelse("Expressed" %in% names(vc), vc["Expressed"], 0)),
      `Not expressed`      = as.integer(ifelse("Not expressed" %in% names(vc), vc["Not expressed"], 0))
    )
  }
  chemo_summary[[fam]] <- fam_counts
}

# ── Report chemo summary table ───────────────────────────────────────────────────
write_rpt("  Chemosensory classification summary:")
write_rpt("  ", sprintf("%-6s", "Family"),
          sprintf("%-10s", "Appendage"),
          sprintf("%10s", "Specific"),
          sprintf("%10s", "Biased"),
          sprintf("%10s", "Expressed"),
          sprintf("%12s", "Not expr"))
write_rpt("  ", paste(rep("-", 58), collapse = ""))
for (fam in fam_order) {
  if (is.null(chemo_summary[[fam]])) next
  for (tkey in tissue_keys) {
    cv <- chemo_summary[[fam]][[tkey]]
    write_rpt("  ", sprintf("%-6s", toupper(fam)),
              sprintf("%-10s", tkey),
              sprintf("%10d", cv["Appendage-specific"]),
              sprintf("%10d", cv["Appendage-biased"]),
              sprintf("%10d", cv["Expressed"]),
              sprintf("%12d", cv["Not expressed"]))
  }
}
write_rpt("")

# ── Draw pie chart ───────────────────────────────────────────────────────────────
make_chemo_pie <- function(counts_vec, fam_label, tissue_label) {
  df <- data.frame(Class = names(counts_vec), Count = as.integer(counts_vec),
                   stringsAsFactors = FALSE)
  df <- df[df$Count > 0, ]
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(fam_label, "-", tissue_label, "(no data)")))
  }
  class_order <- c("Appendage-specific", "Appendage-biased", "Expressed", "Not expressed")
  df$Class <- factor(df$Class, levels = class_order)
  total <- sum(df$Count)

  p <- ggplot(df, aes(x = "", y = Count, fill = Class)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CHEMO_PIE_COLORS, drop = FALSE) +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = paste0(fam_label, " - ", tissue_label, " (n=", total, ")")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      legend.position = "bottom",
      legend.text = element_text(size = 8)
    )
  return(p)
}

# ── Plot each pie individually ───────────────────────────────────────────────────
for (fam in fam_order) {
  if (is.null(chemo_summary[[fam]])) next
  for (tkey in tissue_keys) {
    tissue_label <- tissue_labels_map[tkey]
    cv <- chemo_summary[[fam]][[tkey]]
    print(make_chemo_pie(cv, toupper(fam), tissue_label))
  }
}

cat("Done: Chemosensory pie charts\n")


# =============================================================================
# 5. CHEMOSENSORY GENE IDENTITY REPORT
# =============================================================================

write_rpt("================================================================")
write_rpt("5. CHEMOSENSORY GENE IDENTITY REPORT")
write_rpt("================================================================")
write_rpt("Listing all genes classified as Appendage-specific or Appendage-biased")
write_rpt("per family per appendage.")
write_rpt("")

for (fam in fam_order) {
  fam_df <- chemo_classified[chemo_classified$Family == fam, ]
  if (nrow(fam_df) == 0) next

  write_rpt("  ── ", toupper(fam), " (", nrow(fam_df), " genes) ──")

  for (tkey in tissue_keys) {
    col <- class_cols_map[tkey]
    tissue_label <- tissue_labels_map[tkey]

    specific <- fam_df$Name[fam_df[[col]] == "Appendage-specific"]
    biased   <- fam_df$Name[fam_df[[col]] == "Appendage-biased"]
    expressed <- fam_df$Name[fam_df[[col]] == "Expressed"]

    write_rpt("  ", tissue_label, ":")
    if (length(specific) > 0) {
      write_rpt("    Appendage-specific (", length(specific), "): ",
                paste(sort(specific), collapse = ", "))
    } else {
      write_rpt("    Appendage-specific: none")
    }
    if (length(biased) > 0) {
      write_rpt("    Appendage-biased (", length(biased), "): ",
                paste(sort(biased), collapse = ", "))
    } else {
      write_rpt("    Appendage-biased: none")
    }
    write_rpt("    Expressed (not biased): ", length(expressed))
    write_rpt("    Not expressed: ",
              sum(fam_df[[col]] == "Not expressed"))
  }
  write_rpt("")
}

# ── Summary across all families ──────────────────────────────────────────────────
write_rpt("  ── SUMMARY ACROSS ALL FAMILIES ──")
for (tkey in tissue_keys) {
  col <- class_cols_map[tkey]
  tissue_label <- tissue_labels_map[tkey]
  all_specific <- chemo_classified$Name[chemo_classified[[col]] == "Appendage-specific"]
  all_biased   <- chemo_classified$Name[chemo_classified[[col]] == "Appendage-biased"]
  all_expressed <- chemo_classified$Name[chemo_classified[[col]] == "Expressed"]
  all_not      <- chemo_classified$Name[chemo_classified[[col]] == "Not expressed"]

  write_rpt("")
  write_rpt("  ", tissue_label, " TOTAL:")
  write_rpt("    Appendage-specific: ", length(all_specific))
  if (length(all_specific) > 0) {
    write_rpt("      ", paste(sort(all_specific), collapse = ", "))
  }
  write_rpt("    Appendage-biased: ", length(all_biased))
  if (length(all_biased) > 0) {
    write_rpt("      ", paste(sort(all_biased), collapse = ", "))
  }
  write_rpt("    Expressed: ", length(all_expressed))
  write_rpt("    Not expressed: ", length(all_not))
}
write_rpt("")

# ── Add AA identity columns to classification table ─────────────────────────────
chemo_classified$Has_Identical <- FALSE
chemo_classified$Identical_To  <- NA_character_
for (i in seq_len(nrow(chemo_classified))) {
  key <- toupper(chemo_classified$Name[i])
  members <- dup_cluster_members[[key]]
  if (!is.null(members) && length(members) > 0) {
    chemo_classified$Has_Identical[i] <- TRUE
    chemo_classified$Identical_To[i]  <- paste(sort(members), collapse = ", ")
  }
}

# ── Save classification table as CSV ─────────────────────────────────────────────
class_csv_file <- file.path(FIG_DIR, "Figure2_chemosensory_classification.csv")
write.csv(chemo_classified, class_csv_file, row.names = FALSE)
write_rpt("Classification table saved to: ", class_csv_file)
n_ident <- sum(chemo_classified$Has_Identical)
write_rpt("  Genes with 100% AA identical partner(s): ", n_ident,
          " out of ", nrow(chemo_classified))
write_rpt("")

# =============================================================================
# 6. SUPPLEMENTARY GENE TABLES PER APPENDAGE (with GO, expression, DE stats)
# =============================================================================

write_rpt("================================================================")
write_rpt("6. SUPPLEMENTARY GENE TABLES PER APPENDAGE")
write_rpt("================================================================")

# Helper: build DE lookup keyed on Gene (transcript ID)
make_de_df <- function(de_raw) {
  data.frame(
    Gene            = trimws(as.character(de_raw$Gene)),
    log2FoldChange  = as.numeric(de_raw$log2FoldChange),
    padj            = as.numeric(de_raw$padj),
    stringsAsFactors = FALSE
  )
}
de_antp   <- make_de_df(AntP)
de_antleg <- make_de_df(AntLeg)
de_legp   <- make_de_df(LegP)

# Name lookup from norm (Gene -> Name column if present)
norm_name_map <- if ("Name" %in% names(norm)) {
  setNames(trimws(as.character(norm$Name)), trimws(as.character(norm$Gene)))
} else { NULL }

# Chemo family lookup (transcript key -> family)
chemo_fam_lookup <- setNames(chemo_id$Family, chemo_id$NormKey)

# Per-appendage config: which two DE tables to join + column naming
appendage_de_config <- list(
  Antenna = list(
    de1 = de_antp,   lbl1 = "Ant_vs_Palp",
    de2 = de_antleg, lbl2 = "Ant_vs_Tarsi"
  ),
  `Maxillary palp` = list(
    de1 = de_antp,  lbl1 = "Ant_vs_Palp",
    de2 = de_legp,  lbl2 = "Tarsi_vs_Palp"
  ),
  Tarsi = list(
    de1 = de_antleg, lbl1 = "Ant_vs_Tarsi",
    de2 = de_legp,   lbl2 = "Tarsi_vs_Palp"
  )
)

for (tinfo in tissues_overlap) {
  tissue <- tinfo$name
  keys   <- tinfo$keys
  cfg    <- appendage_de_config[[tissue]]

  base <- data.frame(Gene = keys, stringsAsFactors = FALSE)

  # Gene name
  if (!is.null(norm_name_map)) {
    base$Name <- norm_name_map[base$Gene]
  }

  # GO annotation
  go_sub <- go_map[, c("Gene", "GO_Name", "GO_Domain")]
  base <- merge(base, go_sub, by = "Gene", all.x = TRUE)
  base$GO_Name[is.na(base$GO_Name)]     <- "Unknown"
  base$GO_Domain[is.na(base$GO_Domain)] <- "Unknown"

  # Mean expression per appendage
  base <- merge(base, norm_means, by.x = "Gene", by.y = "JoinKey", all.x = TRUE)

  # DE stats — comparison 1 (copy to avoid mutating shared df)
  de1 <- data.frame(cfg$de1, check.names = FALSE)
  names(de1)[names(de1) == "log2FoldChange"] <- paste0("LFC_", cfg$lbl1)
  names(de1)[names(de1) == "padj"]           <- paste0("padj_", cfg$lbl1)
  base <- merge(base, de1, by = "Gene", all.x = TRUE)

  # DE stats — comparison 2
  de2 <- data.frame(cfg$de2, check.names = FALSE)
  names(de2)[names(de2) == "log2FoldChange"] <- paste0("LFC_", cfg$lbl2)
  names(de2)[names(de2) == "padj"]           <- paste0("padj_", cfg$lbl2)
  base <- merge(base, de2, by = "Gene", all.x = TRUE)

  # Chemosensory family
  base$Chemosensory_family <- chemo_fam_lookup[base$Gene]

  # Sort by GO domain then GO name
  base <- base %>% arrange(GO_Domain, GO_Name, Gene)

  tissue_clean <- gsub(" ", "_", tissue)
  go_csv_file <- file.path(FIG_DIR,
                           paste0("Table_S_GO_names_", tissue_clean, ".csv"))
  write.csv(base, go_csv_file, row.names = FALSE)

  write_rpt("  ", tissue, ": ", nrow(base), " genes -> ", go_csv_file)
  for (dom in c("BP", "CC", "MF", "Unknown")) {
    n_dom   <- sum(base$GO_Domain == dom)
    n_terms <- length(unique(base$GO_Name[base$GO_Domain == dom &
                                           base$GO_Name != "Unknown"]))
    write_rpt("    ", dom, ": ", n_dom, " genes, ", n_terms, " unique GO terms")
  }
}
write_rpt("")

cat("Saved supplementary gene tables per appendage.\n")

# =============================================================================
# CLOSE REPORT
# =============================================================================

write_rpt("")
write_rpt("================================================================")
write_rpt("END OF REPORT")
write_rpt("================================================================")
close(rpt)

cat("\nAll plots displayed.")
cat("\nReport written to:", report_file, "\n")


# ── Session info (run to capture package versions for Methods reporting) ─────
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())

################################################################################
################################################################################
##
##  REVISION SECTIONS — everything below was added for the v5 revision.
##  Panels 3C and 3D above are untouched and reproduce the submitted numbers.
##
################################################################################
################################################################################

suppressPackageStartupMessages({ library(ggrepel); library(tidyr) })

AF_FILE <- function(x) file.path(AF_DIR, x)

save_pdf <- function(p, dir, filename, w, h) {
  dev_fun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(dir, filename), plot = p, width = w, height = h,
         units = "in", device = dev_fun)
  cat("written: ", file.path(dir, filename), "\n", sep = "")
  invisible(p)
}

# =============================================================================
# PANEL 3A (revised) — same volcanoes, with the chemosensory genes labelled
#
# Phi's C18: the panel showed which genes were chemosensory by colour but never
# said which ones. Geometry, counts, colours and the Antenna-on-negative-x
# convention are unchanged; only repel labels are added. The caption now also
# says |log2FC| >= 1 rather than "> 1", which is what the code always did.
# =============================================================================
N_LABEL_PER_SIDE <- 6

# ── Exact volcano geometry, to match the Illustrator layout ──────────────────
# The PLOT PANEL (the bordered box, excluding ticks, axis text and legend) is
# forced to exactly 53.9988 x 52.9798 mm. ggsave sizes the whole device, not the
# panel, so the panel is pinned in the gtable and the device is then sized to
# whatever the finished table needs. Every volcano therefore drops into the
# Illustrator layout at identical size with no rescaling.
#
# Axes are fixed too: x -20..20, y 0..300, shared by every volcano in Figures 3
# and 4 so the appendages can be compared by eye.
#
# A handful of points fall outside those ranges: genes with padj stored as 0
# (DESeq2 underflow) and 13 points beyond |log2FC| 20 - none of the latter is
# significant. They are CAPPED onto the boundary rather than dropped, so nothing
# silently disappears. Say so in the caption.
PANEL_W_MM <- 53.9988
PANEL_H_MM <- 52.9798
X_MAX <- 20
Y_MAX <- 300

cap_y <- function(v) pmin(v, Y_MAX)
cap_x <- function(v) pmax(pmin(v, X_MAX), -X_MAX)

volcano_x_scale <- function()
  scale_x_continuous(limits = c(-X_MAX, X_MAX), breaks = seq(-X_MAX, X_MAX, 10),
                     expand = expansion(mult = 0.02))
volcano_y_scale <- function()
  scale_y_continuous(limits = c(0, Y_MAX), breaks = seq(0, Y_MAX, 100),
                     expand = expansion(mult = c(0.01, 0.04)))

# theme_bw, grid removed
volcano_theme <- function(base = 9)
  theme_bw(base_size = base) +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 9),
        plot.subtitle    = element_text(hjust = 0.5, size = 7, colour = "grey30"),
        axis.text        = element_text(size = 7),
        axis.title       = element_text(size = 8),
        legend.position  = "bottom",
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(3.2, "mm"),
        legend.margin    = margin(t = -2))

# Pin the panel to the exact millimetre size, then size the device to suit.

# ── Volcano gene labels ──────────────────────────────────────────────────────
# Two kinds of label are drawn, in ONE repel layer so they avoid each other:
#   chemosensory  - the Name column (Or182 -> Or182), red, italic
#   everything else - the description from Additional file 11, dark grey, upright
#                     ("alpha-amylase A-like", "farnesol dehydrogenase-like", ...)
# Non-chemosensory genes have no Name - that is how is_chemo is defined - so the
# description is the only usable source. Genes whose description is missing or
# just "uncharacterized LOCxxxxxxx" are skipped and the next-ranked gene is taken
# instead, so a label slot is never wasted on an accession number.
#
# TUNE THESE if the panel is too crowded or too bare. Genes are picked on BOTH
# axes and the two sets are combined, so a gene that is only moderately shifted
# but hugely significant gets a label, and so does a strongly shifted one that is
# less significant. "_X" = ranked by |log2FoldChange|, "_Y" = ranked by -log10(padj).
# A gene must FIRST clear one of these two floors to be eligible for a label at
# all - it has to be genuinely extreme on one axis or the other. Then the caps
# below take the top few of whatever qualifies. Raise the floors to label less.
LAB_MIN_LFC <- 10    # |log2FoldChange| ...
LAB_MIN_Y   <- 100   # ... OR -log10(padj)
# Labels are ranked by a combined score (see lab_score below) and capped.
LAB_TOP_CHEMO <- 30   # per side, best-scoring chemosensory loci
LAB_TOP_OTHER <- 30   # per side, best-scoring other genes past the threshold
LAB_MAXCHAR <- 20
LAB_SIZE    <- 1     # text size of the labels
# TRUE = label one dot per locus; FALSE = label every dot, duplicates
# sharing the same representative name.
LAB_ONE_PER_LOCUS <- TRUE
LAB_OTHER_COL <- "#141414"   # near-black

.desc_raw <- read.csv(file.path(AF_DIR, "Additional_file_11_GO_annotations_all_genes.csv"),
                      stringsAsFactors = FALSE)
DESC_MAP <- setNames(trimws(sub("\\s*\\[Source:.*$", "", .desc_raw$description)),
                     toupper(trimws(.desc_raw$Gene)))

# Second name source: Additional file 27 carries NCBI_gene_name for exactly the
# genes that have no GO term, which is the main reason a description is missing.
.af27 <- tryCatch(read.csv(file.path(AF_DIR, "Additional_file_27_Unknown_BLAST.csv"),
                           stringsAsFactors = FALSE), error = function(e) NULL)
NCBI_MAP <- if (is.null(.af27)) character(0) else {
  k <- toupper(trimws(.af27$Gene)); v <- trimws(.af27$NCBI_gene_name)
  keep <- nzchar(k) & nzchar(v) & !is.na(v)
  setNames(v[keep], k[keep])[!duplicated(k[keep])]
}

other_label <- function(keys) {
  k <- toupper(trimws(keys))
  d <- unname(DESC_MAP[k]); d[is.na(d)] <- ""
  d[grepl("^uncharacterized", d, ignore.case = TRUE)] <- ""
  # fall back to the NCBI name where the description is missing or useless
  f <- unname(NCBI_MAP[k]); f[is.na(f)] <- ""
  f[grepl("^uncharacterized", f, ignore.case = TRUE)] <- ""
  d <- ifelse(nzchar(d), d, f)
  # If a gene has no name in Additional file 11 and none in Additional file 27,
  # it gets NO label - an accession tells the reader nothing. Its dot is still
  # plotted and still counted, it simply goes unnamed. 53 such genes in Antenna
  # vs Tarsi, 69 in Tarsi vs Palp.
  ifelse(nchar(d) > LAB_MAXCHAR, paste0(substr(d, 1, LAB_MAXCHAR - 1), "\u2026"), d)
}

# ── One name per gene locus ──────────────────────────────────────────────────
# Several chemosensory "genes" are redundant RefSeq transcript models of ONE
# locus - Obp8 and Obp9 are both LOC119659806, for instance. Labelling both
# would show the same gene twice under two names. Additional file 37 maps all
# 567 models onto their 393 loci; each locus is labelled once, using the
# lowest-numbered symbol as the representative (Obp8, not Obp9).
.loc <- read.csv(file.path(AF_DIR, "Additional_file_37_transcript_model_to_gene_locus.csv"),
                 stringsAsFactors = FALSE)
.loc$sym <- trimws(.loc$Hill_gene_symbol)
LOCUS_ID <- setNames(trimws(.loc$gene_id), toupper(.loc$sym))
.n   <- suppressWarnings(as.numeric(gsub("^[^0-9]*", "", .loc$sym)))
.ord <- order(.loc$gene_id, ifelse(is.na(.n), Inf, .n), .loc$sym)
.f   <- .loc[.ord, ]
LOCUS_REP <- setNames(.f$sym[!duplicated(.f$gene_id)], .f$gene_id[!duplicated(.f$gene_id)])

# symbol -> the representative symbol for its locus (unchanged if not mapped)
locus_name <- function(sym) {
  k <- toupper(trimws(as.character(sym)))
  id <- unname(LOCUS_ID[k])
  rep <- unname(LOCUS_REP[id])
  ifelse(is.na(rep), as.character(sym), rep)
}
# symbol -> locus id, for collapsing duplicates (falls back to the symbol)
locus_key <- function(sym) {
  k <- toupper(trimws(as.character(sym)))
  id <- unname(LOCUS_ID[k])
  ifelse(is.na(id), k, id)
}

# Ranking score: significance AND effect size together. pi-score (Xiao et al.
# 2014, Bioinformatics 30:801) = |log2FoldChange| x -log10(padj). Ranking on padj
# alone favours high-expression genes with trivial shifts; ranking on fold change
# alone favours noisy low-count genes. The product demands both.
lab_score <- function(x, y) abs(x) * y

# Build the label table. For each side of the volcano and each kind of gene, take
# the top N by |x| AND the top N by y, then combine. Non-chemosensory genes with
# no usable description are dropped BEFORE the top-N cut, so a slot is never
# spent on an accession number.
build_labels <- function(df, x, key, sig = df$is_sig) {
  # EVERY significant chemosensory locus is labelled. The y/|x| threshold exists
  # to thin out the thousands of non-chemosensory genes, and applying it to the
  # chemosensory set as well was discarding almost all of them - 156 of 187 loci
  # in Antenna vs Palp, and every one of them in the Figure 4 panels, because
  # most chemosensory genes sit below y = 100 with fold changes under 10.
  ok_all <- sig & is.finite(x) & is.finite(df$neg_log10)
  ok_thr <- ok_all & (abs(x) >= LAB_MIN_LFC | df$neg_log10 >= LAB_MIN_Y)
  out <- list()
  for (sgn in c(-1, 1)) {
    for (kind in c("chemo", "other")) {
      sel <- (if (kind == "chemo") ok_all else ok_thr) &
             (sign(x) == sgn) & (if (kind == "chemo") df$is_chemo else !df$is_chemo)
      if (!any(sel)) next
      idx <- which(sel)
      sc  <- lab_score(x, df$neg_log10)
      idx <- idx[order(-sc[idx])]          # best-scoring first, both kinds
      if (kind == "chemo") {
        # One label per locus. Rank the members by prominence FIRST, so the name
        # lands on the locus's most significant dot. Keeping whichever row came
        # first in the file instead put "Or3" on the OR7 dot at y = 161 and left
        # the real OR3 dot at y = 208 bare.
        txt <- pretty_gene(locus_name(df$Name[idx]))
        # EVERY chemosensory dot is labelled. Duplicate transcript models of one
        # locus all carry the SAME representative name (Obp8 twice, never Obp8
        # and Obp9), which is the "one name per gene" rule. Dropping the
        # duplicates instead left 63 of 250 red dots bare in Antenna vs Palp.
        if (LAB_ONE_PER_LOCUS) {
          dup <- duplicated(locus_key(df$Name[idx]))
          idx <- idx[!dup]; txt <- txt[!dup]
        }
        ntop <- LAB_TOP_CHEMO
      } else {
        txt <- other_label(df[[key]][idx])
        keep <- nzchar(txt); idx <- idx[keep]; txt <- txt[keep]
        ntop <- LAB_TOP_OTHER
      }
      if (!length(idx)) next
      take <- idx[seq_len(min(ntop, length(idx)))]
      out[[length(out) + 1]] <- data.frame(
        xx = x[take], yy = df$neg_log10[take],
        lab = txt[match(take, idx)], kind = kind, stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  lab <- do.call(rbind, out)
  lab <- lab[nzchar(lab$lab), , drop = FALSE]
  lab[!duplicated(paste(lab$lab, round(lab$xx, 4), round(lab$yy, 4))), , drop = FALSE]
}

volcano_labels <- function(lab) {
  if (is.null(lab) || !nrow(lab)) return(NULL)
  ggrepel::geom_text_repel(
    data = lab, aes(x = xx, y = yy, label = lab), inherit.aes = FALSE,
    colour   = ifelse(lab$kind == "chemo", GO_COLS[["Chemosensory"]], LAB_OTHER_COL),
    fontface = ifelse(lab$kind == "chemo", "italic", "plain"),
    size = LAB_SIZE, segment.size = 0.12, segment.colour = "grey60",
    min.segment.length = 0, box.padding = 0.28, point.padding = 0.12,
    max.overlaps = Inf, seed = 1)
}


# A CSV left open in Excel is locked on Windows, and write.csv then aborts the
# whole script - so panels written earlier survive while everything after the
# failure is silently missing. This writes to the intended file when it can, and
# to "<name>_NEW.csv" when it cannot, warning loudly instead of halting.
safe_write_csv <- function(x, path, ...) {
  ok <- tryCatch({ write.csv(x, path, ...); TRUE }, error = function(e) FALSE,
                 warning = function(w) FALSE)
  if (!ok) {
    alt <- sub("\\.csv$", "_NEW.csv", path)
    tryCatch({
      write.csv(x, alt, ...)
      cat("  !! ", basename(path), " is LOCKED (open in Excel?) - wrote ",
          basename(alt), " instead\n", sep = "")
    }, error = function(e)
      cat("  !! could not write ", basename(path), " or its fallback\n", sep = ""))
  }
  invisible(NULL)
}

save_volcano <- function(p, dir, filename) {
  g   <- ggplot2::ggplotGrob(p)
  pan <- g$layout[g$layout$name == "panel", , drop = FALSE]
  g$widths[unique(pan$l)]  <- grid::unit(PANEL_W_MM, "mm")
  g$heights[unique(pan$t)] <- grid::unit(PANEL_H_MM, "mm")
  grDevices::pdf(NULL)                       # device needed to resolve text units
  # The panel is pinned, but the title and the legend can be WIDER than the panel
  # plus its axes. If the page is sized to the table alone they get clipped at the
  # edge. Measure them and pad the outer margin columns so nothing is cut, leaving
  # the panel size untouched.
  tbl_w <- grid::convertWidth(sum(g$widths), "mm", valueOnly = TRUE)
  extras <- c(0)
  for (nm in c("title", "subtitle", "guide-box-bottom", "guide-box")) {
    k <- which(g$layout$name == nm)
    if (length(k) == 1 && !inherits(g$grobs[[k]], "zeroGrob"))
      extras <- c(extras, grid::convertWidth(grid::grobWidth(g$grobs[[k]]),
                                             "mm", valueOnly = TRUE))
  }
  need <- max(extras) + 2                  # 1 mm breathing room each side
  if (need > tbl_w) {
    pad <- grid::unit((need - tbl_w) / 2, "mm")
    g$widths[1]                <- g$widths[1] + pad
    g$widths[length(g$widths)] <- g$widths[length(g$widths)] + pad
  }
  w <- grid::convertWidth(sum(g$widths),   "mm", valueOnly = TRUE)
  h <- grid::convertHeight(sum(g$heights), "mm", valueOnly = TRUE)
  grDevices::dev.off()
  dev_fun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(dir, filename), plot = g, width = w, height = h,
         units = "mm", device = dev_fun)
  cat(sprintf("written: %-34s  panel %.4f x %.4f mm, page %.2f x %.2f mm
",
              filename, PANEL_W_MM, PANEL_H_MM, w, h))
  invisible(g)
}

# Gene labels follow the manuscript convention: Or182, not OR182. The Name column
# of the Additional files stores them upper-case, so the leading family prefix is
# converted to sentence case for display only. Nothing downstream reads this.
pretty_gene <- function(x) {
  x   <- trimws(as.character(x))
  pre <- sub("^([A-Za-z]+).*$", "\\1", x)
  rst <- sub("^[A-Za-z]+", "", x)
  ifelse(pre == "" | is.na(x), x,
         paste0(toupper(substring(pre, 1, 1)), tolower(substring(pre, 2)), rst))
}

volcano_labelled <- function(vdf, key) {
  vdf$neg_log10 <- cap_y(vdf$neg_log10)
  vdf$plot_lfc  <- cap_x(vdf$plot_lfc)
  lab <- build_labels(vdf, vdf$plot_lfc, "JoinKey")
  p <- make_volcano_plot(vdf) +
    volcano_labels(lab) +
    volcano_x_scale() + volcano_y_scale() + volcano_theme() +
    labs(title = NULL)          # set titles in Illustrator; counts are in the CSV
  cat("   3A ", key, "  ", volcano_title(vdf), "\n", sep = "")
  print(p)
  save_volcano(p, MAIN_OUT, paste0("Figure_3A_", key, ".pdf"))
}

cat("\n--- Panel 3A (revised, labelled) ---\n")
volcano_labelled(antp_v,   "Ant_vs_Palp")
volcano_labelled(antleg_v, "Ant_vs_Tarsi")
volcano_labelled(legp_v,   "Tarsi_vs_Palp")

# =============================================================================
# PANEL 3B (replaced) — chemosensory versus the rest, with a test
#
# The published panel split each DEG set into GO domains. Additional file 11
# carries exactly ONE GO term per transcript (0 of 24,828 have more than one),
# so that split describes the annotation pipeline, not the biology — which is
# what C18 objected to, and C20/C22 asked for statistics the percentage stack
# could not carry. The panel now shows chemosensory versus everything else with
# a one-sided Fisher's exact test against the full annotated background.
#
# Target: chemosensory enrichment among antenna-biased genes OR 6.26 vs palp and
# OR 5.33 vs tarsi; palp and tarsi are NOT enriched.
# =============================================================================
CHEMO_SPLIT_COLS <- c(Chemosensory = "#A60000", `All other genes` = "#BFBFBF")

fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 2.2e-16) "p < 2.2e-16" else
  paste0("p = ", format(p, digits = 2, scientific = TRUE))

# Odds ratio above each bar, significance as stars. The p-values themselves are
# either astronomically small or plainly null, so printing them adds nothing the
# stars do not; the confidence intervals live in the caption and the stats CSV.
sig_stars <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else
  if (p < 0.01) "**" else if (p < 0.05) "*" else " n.s."

chemo_fisher <- function(gene_keys, universe_df) {
  inset <- universe_df$JoinKey %in% gene_keys
  ch    <- universe_df$is_chemo
  ft <- fisher.test(matrix(c(sum(inset & ch), sum(inset & !ch),
                             sum(!inset & ch), sum(!inset & !ch)), 2, byrow = TRUE))
  list(n_set = sum(inset), n_chemo = sum(inset & ch),
       odds = unname(ft$estimate), p = ft$p.value)
}

cat("\n--- Panel 3B (replaced: chemosensory vs rest) ---\n")
b_rows <- list()
make_chemo_bar <- function(vdf, key) {
  disp_pos <- attr(vdf, "disp_left"); disp_neg <- attr(vdf, "disp_right")
  sig <- vdf[vdf$is_sig, ]
  sig$Dir <- ifelse(sig$plot_lfc > 0, paste0(disp_pos, "-biased"),
                                      paste0(disp_neg, "-biased"))
  dir_order <- c(paste0(disp_neg, "-biased"), paste0(disp_pos, "-biased"))
  sig$Split <- factor(ifelse(sig$is_chemo, "Chemosensory", "All other genes"),
                      levels = c("Chemosensory", "All other genes"))

  tbl <- sig %>% group_by(Dir, Split, .drop = FALSE) %>%
    summarise(N = dplyr::n(), .groups = "drop") %>%
    group_by(Dir) %>% mutate(Percent = 100 * N / sum(N)) %>% ungroup() %>%
    mutate(TextLabel = ifelse(N > 0, paste0("N = ", N, "\n", sprintf("%.1f%%", Percent)), ""))
  chemo_lab <- tbl %>% filter(Split == "Chemosensory")
  tbl$TextLabel[tbl$Split == "Chemosensory"] <- ""   # drawn above the bar instead

  ann <- character(0)
  for (d in dir_order) {
    ft <- chemo_fisher(sig$JoinKey[sig$Dir == d], vdf)
    b_rows[[length(b_rows) + 1]] <<- data.frame(Contrast = key, Direction = d,
      n_DEG = ft$n_set, n_chemosensory = ft$n_chemo,
      pct_chemosensory = round(100 * ft$n_chemo / max(ft$n_set, 1), 2),
      odds_ratio = round(ft$odds, 3), p_value = ft$p, stringsAsFactors = FALSE)
    ann <- c(ann, sprintf("%.2f%s", ft$odds, sig_stars(ft$p)))
    cat(sprintf("  %-28s %-24s OR = %6.2f  %s\n", key, d, ft$odds, fmt_p(ft$p)))
  }
  ann_df <- data.frame(Dir = factor(dir_order, levels = dir_order), lab = ann)
  n_tot  <- sapply(dir_order, function(d) sum(sig$Dir == d))

  p <- ggplot(tbl, aes(x = factor(Dir, levels = dir_order), y = Percent, fill = Split)) +
    geom_bar(stat = "identity", position = "stack", width = 0.62) +
    geom_text(aes(label = TextLabel), position = position_stack(vjust = 0.5),
              size = 2.8, colour = "white", fontface = "bold", lineheight = 0.85) +
    geom_text(data = chemo_lab, aes(x = Dir, y = 101.5, label = TextLabel),
              inherit.aes = FALSE, size = 2.8, fontface = "bold", vjust = 0,
              lineheight = 0.85, colour = CHEMO_SPLIT_COLS[["Chemosensory"]]) +
    geom_text(data = ann_df, aes(x = Dir, y = 114, label = lab),
              inherit.aes = FALSE, size = 2.5, colour = "grey25") +
    scale_fill_manual(values = CHEMO_SPLIT_COLS, name = NULL) +
    scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = paste0(disp_neg, " vs ", disp_pos), x = NULL, y = "% of significant DEGs",
         subtitle = paste0(dir_order[1], ": ", n_tot[1], " DEGs\n",
                           dir_order[2], ": ", n_tot[2], " DEGs")) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(hjust = .5, face = "bold", size = 10),
          plot.subtitle = element_text(hjust = .5, size = 7.5, colour = "grey30", lineheight = 1.1),
          legend.position = "bottom", legend.text = element_text(size = 8),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  print(p)
  save_pdf(p, MAIN_OUT, paste0("Figure_3B_", key, ".pdf"), 4.0, 4.8)
}
make_chemo_bar(antp_v,   "Ant_vs_Palp")
make_chemo_bar(antleg_v, "Ant_vs_Tarsi")
make_chemo_bar(legp_v,   "Tarsi_vs_Palp")
safe_write_csv(bind_rows(b_rows), file.path(MAIN_OUT, "Figure_3B_enrichment_stats.csv"),
          row.names = FALSE)

# =============================================================================
# SUPPLEMENTARY FIGURES S9-S11 — GO over-representation, appendage-biased sets
#
# The submitted panels ranked GO terms by raw gene count with no test at all,
# which is what C20/C24/C61 objected to. This is a real over-representation
# analysis: one-sided Fisher's exact per term, Benjamini-Hochberg corrected
# within each gene set x domain, non-significant bars drawn pale.
#
# S11 is the important one. Its submitted Illustrator boxes sit on protein
# binding (OR 1.10, q = 0.77) and nucleotide binding (OR 0.69, q = 1.00 — that
# is DEPLETED). Redo the boxes on q < 0.05 terms only.
#
# Validate against 07_Revision_Record/GO_enrichment_tables/:
#   Antenna MF olfactory receptor activity  k = 127/278, K = 186/12026,
#   OR 166.6, q = 1.17e-171
# =============================================================================
# ── Supplementary page geometry (Perets et al., 2026-10) ─────────────────────
# Every supplementary figure except the family trees (S1-S8) is an A4 page with
# the content occupying 180 x 260 mm, centred: side margins (210-180)/2 = 15 mm,
# top and bottom (297-260)/2 = 18.5 mm. Minimum text size is 6 pt throughout;
# geom_text takes millimetres, so 6 pt = 6/2.845 = 2.11.
SUPP_W_IN  <- 8.268            # 210 mm
SUPP_H_IN  <- 11.693           # 297 mm
SUPP_MARGIN <- ggplot2::margin(18.5, 15, 18.5, 15, unit = "mm")
PT_MIN     <- 6                # smallest theme text, in points
MM_MIN     <- 2.11             # smallest geom_text, in millimetres (= 6 pt)

Q_THR <- 0.05

read_raw_de <- function(f) {
  d <- read.csv(AF_FILE(f), stringsAsFactors = FALSE)
  data.frame(JoinKey = trimws(as.character(d$Gene)),
             lfc  = suppressWarnings(as.numeric(d$log2FoldChange)),
             padj = suppressWarnings(as.numeric(d$padj)), stringsAsFactors = FALSE)
}
.AT <- read_raw_de("Additional_file_1_Appendage_DE_Antenna_vs_Tarsi.csv")          # + = Antenna
.AP <- read_raw_de("Additional_file_2_Appendage_DE_Antenna_vs_MaxillaryPalp.csv")  # + = Antenna
.TP <- read_raw_de("Additional_file_3_Appendage_DE_Tarsi_vs_MaxillaryPalp.csv")    # + = Tarsi
.up <- function(d) d$JoinKey[!is.na(d$padj) & d$padj < PADJ_THR & !is.na(d$lfc) & d$lfc >=  LFC_THR]
.dn <- function(d) d$JoinKey[!is.na(d$padj) & d$padj < PADJ_THR & !is.na(d$lfc) & d$lfc <= -LFC_THR]

CONSENSUS <- list(
  Antenna          = intersect(.up(.AT), .up(.AP)),
  `Maxillary palp` = intersect(.dn(.AP), .dn(.TP)),
  Tarsi            = intersect(.up(.TP), .dn(.AT)))
ORA_UNIVERSE <- Reduce(intersect, list(
  .AT$JoinKey[!is.na(.AT$padj)], .AP$JoinKey[!is.na(.AP$padj)], .TP$JoinKey[!is.na(.TP$padj)]))

cat("\n--- Consensus appendage-biased sets (Figure 3C targets 494 / 539 / 4161) ---\n")
for (a in names(CONSENSUS)) cat(sprintf("  %-15s %d\n", a, length(CONSENSUS[[a]])))

# One real GO annotation per gene
GO_ANN <- go_map[go_map$GO_Domain %in% c("MF", "CC", "BP") &
                 !go_map$GO_Name %in% c("", "Unknown"), c("Gene", "GO_Domain", "GO_Name")]
names(GO_ANN)[1] <- "JoinKey"          # go_map keys on the column named Gene
GO_ANN <- GO_ANN[!duplicated(GO_ANN$JoinKey), ]

ora_panel <- function(set_keys, universe_keys, dom, top_n = 20) {
  bg  <- GO_ANN[GO_ANN$GO_Domain == dom & GO_ANN$JoinKey %in% universe_keys, ]
  ins <- bg[bg$JoinKey %in% set_keys, ]
  n <- nrow(ins); N <- nrow(bg)
  if (n == 0 || N == 0) return(NULL)
  res <- bind_rows(lapply(unique(ins$GO_Name), function(t) {
    k <- sum(ins$GO_Name == t); K <- sum(bg$GO_Name == t)
    ft <- fisher.test(matrix(c(k, n - k, K - k, N - n - (K - k)), 2, byrow = TRUE),
                      alternative = "greater")
    data.frame(GO_Name = t, k = k, n = n, K = K, N = N,
               odds_ratio = unname(ft$estimate), p = ft$p.value, stringsAsFactors = FALSE)
  }))
  res$q <- p.adjust(res$p, method = "BH")
  head(res[order(res$q, -res$odds_ratio), ], top_n)
}
wrap_lab <- function(x, w) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"),
                                  character(1), USE.NAMES = FALSE)

# Domain colours as in the originally submitted figures - molecular function
# green, cellular component blue, biological process purple - with a pale tint
# of the same hue where the term does not reach q < 0.05.
DOM_LAB  <- c("Molecular function", "Cellular component", "Biological process")
DOM_FILL <- setNames(
  c("#3FC498", "#4D50DB", "#9B5DE5", "#CCEDE1", "#CDCEF4", "#E4D2F7"),
  c(DOM_LAB, paste0(DOM_LAB, " (n.s.)")))

cat("\n--- Supplementary figures S9-S11 ---\n")
S9_11_FIG <- c(Antenna = "S9", `Maxillary palp` = "S10", Tarsi = "S11")
ora_all <- list()
for (a in names(CONSENSUS)) {
  rows <- list()
  for (dom in c("MF", "CC", "BP")) {
    r <- ora_panel(CONSENSUS[[a]], ORA_UNIVERSE, dom)
    if (is.null(r)) next
    r$Appendage <- a; r$Domain <- dom; rows[[dom]] <- r
  }
  tbl <- bind_rows(rows)
  if (!nrow(tbl)) next
  ora_all[[a]] <- tbl
  lead <- tbl[tbl$Domain == "MF", ][1, ]
  cat(sprintf("  %-15s -> %-4s  leading MF term: %s (k=%d/%d, OR=%.1f, q=%.2e); %d terms at q<0.05\n",
              a, S9_11_FIG[[a]], lead$GO_Name, lead$k, lead$n, lead$odds_ratio, lead$q,
              sum(tbl$q < Q_THR)))

  pl <- tbl
  pl$Domain <- factor(pl$Domain, levels = c("MF", "CC", "BP"),
                      labels = c("Molecular function", "Cellular component",
                                 "Biological process"))
  pl$score <- -log10(pmax(pl$q, .Machine$double.xmin))
  pl$sig   <- pl$q < Q_THR
  # A SHARED y scale cannot label facets that hold different terms at the same
  # integer position - every term name piled on top of every other, which is why
  # the previous build was unreadable. Make y a factor of facet + term, keep the
  # facets free, and strip the sort prefix in the labels.
  pl <- pl %>% group_by(Domain) %>% arrange(q, .by_group = TRUE) %>%
    mutate(rank = rev(seq_len(dplyr::n()))) %>% ungroup()
  pl$ykey <- sprintf("%d@@%02d@@%s", as.integer(pl$Domain), pl$rank, pl$GO_Name)
  pl$ykey <- factor(pl$ykey, levels = pl$ykey[order(as.integer(pl$Domain), pl$rank)])
  pl$fillkey <- paste0(as.character(pl$Domain), ifelse(pl$sig, "", " (n.s.)"))
  pl$fillkey <- factor(pl$fillkey, levels = as.vector(rbind(DOM_LAB, paste0(DOM_LAB, " (n.s.)"))))

  p <- ggplot(pl, aes(x = score, y = ykey, fill = fillkey)) +
    geom_col(width = 0.7, position = "identity") +
    geom_vline(xintercept = -log10(Q_THR), linetype = "dashed",
               linewidth = 0.4, colour = "#C00000") +
    geom_text(aes(label = k), hjust = -0.25, size = MM_MIN, colour = "grey25") +
    facet_wrap(~ Domain, ncol = 1, scales = "free") +
    scale_y_discrete(labels = function(x) wrap_lab(sub("^\\d+@@\\d+@@", "", x), 44),
                     expand = expansion(add = 0.6)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    scale_fill_manual(values = DOM_FILL, guide = "none") +
    labs(title = paste0("Figure ", S9_11_FIG[[a]],
                        ". Gene Ontology over-representation among ", tolower(a),
                        "-biased genes"),
         subtitle = paste0("One-sided Fisher's exact test against all annotated genes tested in ",
                           "all three pairwise contrasts, Benjamini-Hochberg corrected within each ",
                           "domain. Bar = -log10(q); numbers are gene counts; dashed line marks ",
                           "q = 0.05 and pale bars did not reach it."),
         x = expression(-log[10]*"(BH-adjusted "*italic(p)*")"), y = NULL) +
    theme_bw(base_size = 8) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 6.4, colour = "grey30"),
          axis.text.y   = element_text(size = PT_MIN),
          strip.text    = element_text(size = 8, face = "bold"),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  print(p)
  save_pdf(p + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
           paste0("Figure_", S9_11_FIG[[a]], ".pdf"), SUPP_W_IN, SUPP_H_IN)
}
safe_write_csv(bind_rows(ora_all), file.path(SUPP_OUT, "ORA_appendage_consensus.csv"),
          row.names = FALSE)

# =============================================================================
# SUPPLEMENTARY FIGURES S13-S19 — per-family expression heatmaps
#
# Ported from generate_heatmaps.py, which could no longer run: it pointed at
# Submission_ready/app/data and BSF2026/Supplemetary_Files/, both deleted. The
# plots are unchanged in design; only the captions gain "transcript models".
#
#   colour   log2(mean normalised count + 1) - log2(11), so white = the
#            expression threshold of 10 used throughout
#   scale    shared across families, set from the 99th percentile of |value|
#            over the six original families then tightened 4-fold, so adding TRP
#            never shifts the earlier panels
#   rows     grouped by appendage of peak expression, descending within group
# =============================================================================
THRESHOLD <- 10.0
THR_LOG   <- log2(THRESHOLD + 1)
FAM_FIG   <- c(Or = "S13", Gr = "S14", Ir = "S15", Obp = "S16",
               Ppk = "S17", Csp = "S18", Trp = "S19")   # renumbered 2026-09-27
FAM_LONG  <- c(Or = "odorant receptor (OR)", Gr = "gustatory receptor (GR)",
               Ir = "ionotropic receptor (IR)", Obp = "odorant-binding protein (OBP)",
               Ppk = "pickpocket (PPK)", Csp = "chemosensory protein (CSP)",
               Trp = "transient receptor potential (TRP)")
GROUPS <- c("Ant_Vm", "Ant_VF", "Ant_MF", "P_Vm", "P_VF", "P_MF",
            "Leg_Vm", "Leg_VF", "Leg_MF")
GROUP_LAB <- c("Ant\nVm", "Ant\nVF", "Ant\nMF", "Palp\nVm", "Palp\nVF",
               "Palp\nMF", "Tarsi\nVm", "Tarsi\nVF", "Tarsi\nMF")

hm <- read.csv(norm_file, stringsAsFactors = FALSE, check.names = FALSE)
hm$Name <- trimws(ifelse(is.na(hm$Name), "", as.character(hm$Name)))
hm <- hm[hm$Name != "", ]
hm$Family <- extract_family(hm$Name)

grp_mean <- sapply(GROUPS, function(g) {
  cols <- grep(paste0("^", g, "[0-9]+$"), names(hm), value = TRUE)
  if (!length(cols)) return(rep(NA_real_, nrow(hm)))
  rowMeans(sapply(hm[cols], function(x) suppressWarnings(as.numeric(x))), na.rm = TRUE)
})
colnames(grp_mean) <- GROUPS
CENT <- log2(grp_mean + 1) - THR_LOG

anchor <- hm$Family %in% c("Or", "Gr", "Ir", "Obp", "Ppk", "Csp")
GLOBAL_LIM <- max(quantile(abs(CENT[anchor, ]), 0.99, na.rm = TRUE), 0.5) / 4
cat(sprintf("\n--- Supplementary figures S13-S19 ---\n  colour limit +/-%.4f\n", GLOBAL_LIM))

for (fam in names(FAM_FIG)) {
  idx <- which(hm$Family == fam)
  if (!length(idx)) { cat("  ", fam, ": no genes\n"); next }
  m <- CENT[idx, , drop = FALSE]; nm <- hm$Name[idx]

  ant <- rowMeans(m[, 1:3, drop = FALSE]); palp <- rowMeans(m[, 4:6, drop = FALSE])
  tar <- rowMeans(m[, 7:9, drop = FALSE])
  peak <- max.col(cbind(ant, palp, tar), ties.method = "first")
  peakval <- cbind(ant, palp, tar)[cbind(seq_along(peak), peak)]
  ord <- unlist(lapply(1:3, function(t) { i <- which(peak == t); i[order(-peakval[i])] }))
  m <- m[ord, , drop = FALSE]; nm <- nm[ord]

  n_genes <- length(nm)
  two_col <- n_genes > 80
  half    <- if (two_col) ceiling(n_genes / 2) else n_genes
  block   <- ifelse(seq_len(n_genes) <= half, "A", "B")
  rowpos  <- ifelse(block == "A", seq_len(n_genes), seq_len(n_genes) - half)

  df <- data.frame(gene = rep(nm, times = length(GROUPS)),
                   block = rep(block, times = length(GROUPS)),
                   rowpos = rep(rowpos, times = length(GROUPS)),
                   grp = rep(GROUPS, each = n_genes),
                   val = as.vector(m), stringsAsFactors = FALSE)
  df$grp <- factor(df$grp, levels = GROUPS, labels = GROUP_LAB)
  lab_map <- setNames(nm, paste(block, rowpos))
  fs <- max(3.2, min(7, 260 / max(half, 1)))

  p <- ggplot(df, aes(x = grp, y = rowpos, fill = pmax(pmin(val, GLOBAL_LIM), -GLOBAL_LIM))) +
    geom_tile(colour = "grey92", linewidth = 0.05) +
    facet_wrap(~ block, nrow = 1, scales = "free_y") +
    scale_y_reverse(breaks = df$rowpos,
                    labels = function(b) lab_map[paste(df$block[match(b, df$rowpos)], b)],
                    expand = expansion(add = 0.5)) +
    scale_fill_gradient2(low = "#1f4aa8", mid = "#ffffff", high = "#b2182b",
                         midpoint = 0, limits = c(-GLOBAL_LIM, GLOBAL_LIM),
                         name = expression(log[2]*"(counts+1) - "*log[2]*"(11)")) +
    labs(title = paste0("Figure ", FAM_FIG[[fam]], ". Expression heatmap of the H. illucens ",
                        FAM_LONG[[fam]], " family"),
         subtitle = paste0(n_genes, " annotated transcript models across the three appendages ",
                           "and three reproductive states. White marks the expression threshold ",
                           "of 10 normalised counts; the colour scale is shared by S13-S19."),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 6.4, colour = "grey30"),
          axis.text.y   = element_text(size = fs),
          axis.text.x   = element_text(size = 6, lineheight = 0.85),
          strip.text    = element_blank(),
          legend.position = "bottom", legend.key.width = unit(1.4, "cm"),
          panel.grid = element_blank())
  print(p)
  save_pdf(p + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
             paste0("Figure_", FAM_FIG[[fam]], ".pdf"), SUPP_W_IN, SUPP_H_IN)
  cat(sprintf("  %-4s %3d models -> %s%s\n", fam, n_genes, FAM_FIG[[fam]],
              ifelse(two_col, "  (two columns)", "")))
}


# =============================================================================
# PANEL 3A — ALTERNATIVE VERSION, no GO colouring
#
# Requested 2026-09-17. Same geometry, same counts, same labels as the panels
# above; only the colour scheme changes. Chemosensory genes stay red, every
# other significant gene is light blue, and non-significant genes stay pale
# grey so the volcano still reads as a volcano.
#
# These are written ALONGSIDE the standard panels as Figure_3A_alt_*.pdf.
# Nothing is overwritten — pick whichever set you prefer in Illustrator.
#
# If you want literally every non-chemosensory point light blue, including the
# non-significant cloud, set ALT_BLUE_INCLUDES_NONSIG <- TRUE below.
# =============================================================================
ALT_CHEMO_COL   <- "#A60000"
ALT_OTHER_COL   <- "#9ECAE1"   # light blue
ALT_NONSIG_COL  <- "#E4E4E4"
ALT_BLUE_INCLUDES_NONSIG <- FALSE

volcano_alt <- function(vdf, key, labels = TRUE) {
  vdf$neg_log10 <- cap_y(vdf$neg_log10)
  vdf$plot_lfc  <- cap_x(vdf$plot_lfc)
  vdf$AltKey <- ifelse(vdf$is_chemo & vdf$is_sig, "Chemosensory",
                ifelse(vdf$is_sig | ALT_BLUE_INCLUDES_NONSIG,
                       "All other genes", "Not significant"))
  vdf$AltKey <- factor(vdf$AltKey,
                       levels = c("Not significant", "All other genes", "Chemosensory"))
  vdf <- vdf[order(vdf$AltKey), ]          # chemosensory drawn last, on top

  lab <- build_labels(vdf, vdf$plot_lfc, "JoinKey")

  p <- ggplot(vdf, aes(x = plot_lfc, y = neg_log10, colour = AltKey)) +
    geom_point(size = 0.8, alpha = 0.65) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", linewidth = 0.3) +
    (if (labels) volcano_labels(lab)) +
    scale_colour_manual(
      values = c(`Not significant` = ALT_NONSIG_COL,
                 `All other genes` = ALT_OTHER_COL,
                 Chemosensory      = ALT_CHEMO_COL),
      breaks = c("Chemosensory", "All other genes", "Not significant"),
      drop = FALSE, name = NULL) +
    volcano_y_scale() +
    labs(title = NULL,
         x = expression(Log[2](FoldChange)),
         y = expression(-log[10](italic(padj)))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "bottom", legend.text = element_text(size = 8),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    volcano_x_scale() + volcano_theme()
  tag <- if (labels) "alt_" else "nolab_"
  cat("   3A ", tag, key, "  ", volcano_title(vdf), "\n", sep = "")
  print(p)
  save_volcano(p, MAIN_OUT, paste0("Figure_3A_", tag, key, ".pdf"))
}

cat("\n--- Panel 3A (alternative: chemosensory vs light blue, no GO) ---\n")
volcano_alt(antp_v,   "Ant_vs_Palp")
volcano_alt(antleg_v, "Ant_vs_Tarsi")
volcano_alt(legp_v,   "Tarsi_vs_Palp")

# Same panels with NO gene labels at all, written alongside the labelled ones.
cat("\n--- Panel 3A (no labels) ---\n")
volcano_alt(antp_v,   "Ant_vs_Palp",   labels = FALSE)
volcano_alt(antleg_v, "Ant_vs_Tarsi",  labels = FALSE)
volcano_alt(legp_v,   "Tarsi_vs_Palp", labels = FALSE)


# =============================================================================
# SUPPLEMENTARY FIGURE S12 + ADDITIONAL FILE 38
#
# S12  GO-domain composition of each appendage-biased DEG set - the breakdown
#      that used to be panel 3B. DESCRIPTIVE ONLY: Additional file 11 assigns
#      exactly one GO term per transcript, so the split describes the annotation,
#      not the biology (Phi, C18). Tested enrichment is in Figures S9-S11 and
#      Additional file 38.
#
# AF38 Statistical report behind panel 3B: chemosensory enrichment per
#      appendage-biased set with odds ratios, 95% confidence intervals and exact
#      p-values, plus the between-appendage comparison.
#
# Uses this script's own prepped frames (antp_v / antleg_v / legp_v); positive
# log2FoldChange = higher in the left-named appendage, per prep_volcano().
# =============================================================================
cat("\n--- Supplementary figure S12 and Additional file 38 ---\n")

AF_OUT <- file.path(PKG, "04_Additional_Files")
S12_SETS <- list(list(v = antp_v,   key = "Ant_vs_Palp"),
                 list(v = antleg_v, key = "Ant_vs_Tarsi"),
                 list(v = legp_v,   key = "Tarsi_vs_Palp"))

go_rows <- list(); enr_rows <- list()
for (it in S12_SETS) {
  v  <- it$v
  Lname <- attr(v, "left_name"); Rname <- attr(v, "right_name")
  sig <- v[v$is_sig, ]
  for (s in c(1, -1)) {
    nm  <- if (s > 0) Lname else Rname
    set <- sig[sign(sig$log2FoldChange) == s, ]
    dom <- ifelse(set$is_chemo, "Chemosensory", set$GO_Domain)
    tb  <- as.data.frame(table(factor(dom,
             levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))))
    names(tb) <- c("Category", "N")
    tb$Contrast <- it$key; tb$Appendage <- nm
    tb$Percent  <- 100 * tb$N / max(sum(tb$N), 1)
    go_rows[[length(go_rows) + 1]] <- tb

    inset <- v$JoinKey %in% set$JoinKey
    tab <- matrix(c(sum(inset & v$is_chemo),  sum(inset & !v$is_chemo),
                    sum(!inset & v$is_chemo), sum(!inset & !v$is_chemo)), 2, byrow = TRUE)
    ft <- fisher.test(tab)
    enr_rows[[length(enr_rows) + 1]] <- data.frame(
      Comparison = paste0(Lname, " vs ", Rname), Gene_set = paste0(nm, "-biased"),
      Test = "chemosensory vs all expressed genes",
      n_DEG = sum(inset), n_chemosensory = tab[1, 1],
      pct_chemosensory = round(100 * tab[1, 1] / max(sum(inset), 1), 2),
      odds_ratio = round(unname(ft$estimate), 3),
      CI95_low = round(ft$conf.int[1], 3), CI95_high = round(ft$conf.int[2], 3),
      p_value = ft$p.value, stringsAsFactors = FALSE)
  }
  a <- sig[sign(sig$log2FoldChange) > 0, ]; b <- sig[sign(sig$log2FoldChange) < 0, ]
  tab <- matrix(c(sum(a$is_chemo), nrow(a) - sum(a$is_chemo),
                  sum(b$is_chemo), nrow(b) - sum(b$is_chemo)), 2, byrow = TRUE)
  ft <- fisher.test(tab)
  enr_rows[[length(enr_rows) + 1]] <- data.frame(
    Comparison = paste0(Lname, " vs ", Rname),
    Gene_set = paste0(Lname, "-biased vs ", Rname, "-biased"),
    Test = "between appendages", n_DEG = nrow(a) + nrow(b),
    n_chemosensory = sum(a$is_chemo) + sum(b$is_chemo), pct_chemosensory = NA,
    odds_ratio = round(unname(ft$estimate), 3),
    CI95_low = round(ft$conf.int[1], 3), CI95_high = round(ft$conf.int[2], 3),
    p_value = ft$p.value, stringsAsFactors = FALSE)
}

ENR <- bind_rows(enr_rows)
safe_write_csv(ENR, file.path(AF_OUT, "Additional_file_38_Figure3_statistical_tests.csv"),
               row.names = FALSE)
cat("written: Additional_file_38_Figure3_statistical_tests.csv\n")
print(ENR[, c("Gene_set", "Test", "odds_ratio", "CI95_low", "CI95_high", "p_value")])

GO <- bind_rows(go_rows)
GO$Category  <- factor(GO$Category, levels = c("Chemosensory","MF","CC","BP","Unknown"))
GO$Contrast  <- factor(GO$Contrast, levels = c("Ant_vs_Palp","Ant_vs_Tarsi","Tarsi_vs_Palp"),
                       labels = c("Antenna vs maxillary palp","Antenna vs tarsi",
                                  "Tarsi vs maxillary palp"))
p24 <- ggplot(GO, aes(x = Appendage, y = Percent, fill = Category)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(aes(label = ifelse(Percent >= 4, sprintf("%.1f", Percent), "")),
            position = position_stack(vjust = 0.5), size = 2.4, colour = "white") +
  facet_wrap(~ Contrast, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(Chemosensory = "#A60000", MF = "#3FC498",
                               CC = "#4D50DB", BP = "#A860E3", Unknown = "#666666"),
                    labels = c("Chemosensory","Molecular function","Cellular component",
                               "Biological process","No GO term"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Figure S12. Gene Ontology domain composition of appendage-biased genes",
       subtitle = paste("Descriptive only. Additional file 11 assigns exactly one GO term per",
                        "transcript, so this\nsplit reflects the annotation, not a tested",
                        "property of the gene sets. Tested enrichment is\nin Figures S9-S11;",
                        "chemosensory enrichment statistics are in Additional file 38."),
       x = NULL, y = "% of differentially expressed genes") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 6.6, colour = "grey30"),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", colour = NA))
print(p24)
save_pdf(p24 + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
         "Figure_S12.pdf", SUPP_W_IN, SUPP_H_IN)

cat("\n================ CHECK AGAINST THE TARGETS ================\n")
cat("3A titles       : Antenna 215/2372, Maxillary palp 35/1350 (Ant vs Palp)\n")
cat("3B enrichment   : antenna-biased OR 6.26 vs palp, OR 5.33 vs tarsi; palp and tarsi not enriched\n")
cat("3C consensus    : ", paste(sapply(CONSENSUS, length), collapse = " / "),
    "   (target 494 / 539 / 4161)\n", sep = "")
cat("3D classification: antenna 95 specific + 54 biased; palp 13 + 15; tarsi 23 + 75\n")
cat("S9  antenna MF  : olfactory receptor activity k = 127/278, OR 166.6, q = 1.17e-171\n")
cat("S10 palp MF     : solute:inorganic anion antiporter activity, 12 genes, OR 64.9\n")
cat("S11 tarsi MF    : GPCR activity 129 genes, OR 6.9 — protein binding and nucleotide\n")
cat("                  binding are NOT significant; redo the Illustrator boxes accordingly\n")
cat("\nIf any line disagrees with its target, stop and report it.\n")
cat("Done: Figure 3 and Supplementary Figures S9-S12 and S13-S19.\n")
