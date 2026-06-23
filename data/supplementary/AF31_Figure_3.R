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

# ── Paths ────────────────────────────────────────────────────────────────────────
BASE_DIR   <- "C:/Users/dorpe/Downloads/BSF2026/BSF2026/Supplemetary_Files/csv"
FIG_DIR    <- "C:/Users/dorpe/Downloads/BSF2026/BSF2026/Supplemetary_Files/Figure_2"

ant_p_file   <- file.path(BASE_DIR, "Appendage_DE", "condition_vs_Ant_vs_P_name1.csv")
ant_leg_file <- file.path(BASE_DIR, "Appendage_DE", "condition_vs_Ant_vs_Leg_name1.csv")
leg_p_file   <- file.path(BASE_DIR, "Appendage_DE", "condition_vs_Leg_vs_P_name1.csv")
go_file      <- file.path(BASE_DIR, "BSF_all-rna_GO_ID_annotated.csv")
norm_file    <- file.path(BASE_DIR, "BSF_normalized_counts_master.csv")
chemo_id_file <- file.path(BASE_DIR, "BSF_Olfactory_ID.csv")

report_file <- file.path(FIG_DIR, "Figure2_report.txt")

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

  sup_pdf <- file.path(FIG_DIR, paste0("Figure_", fig_label, "_GO_bars_",
                                        gsub(" ", "_", tissue), ".pdf"))
  ggsave(sup_pdf, combined, width = 8.27, height = 11.69, units = "in")  # A4
  cat("Saved supplementary figure:", sup_pdf, "\n")
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

sim_matrix_dir <- file.path(BASE_DIR, "AA_chemosensory_similarity_matrix")
sim_matrix_files <- list(
  Or  = file.path(sim_matrix_dir, "BSF_OR_matrix.csv"),
  Gr  = file.path(sim_matrix_dir, "BSF_GR_matrix.csv"),
  Ir  = file.path(sim_matrix_dir, "BSF_IR_matrix.csv"),
  Obp = file.path(sim_matrix_dir, "BSF_OBP_matrix.csv"),
  Csp = file.path(sim_matrix_dir, "BSF_CSP_matrix.csv"),
  Ppk = file.path(sim_matrix_dir, "BSF_PPK_matrix.csv"),
  Trp = file.path(sim_matrix_dir, "BSF_TRP_matrix.csv")
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
