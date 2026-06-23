################################################################################
# Figure 3 — Sex & Reproductive State DE analysis
#
# Perets et al. 2026 — adult Hermetia illucens chemosensory appendages
#
# Usage:
#   1. Open R / RStudio
#   2. Source this file — all plots print to screen, report to Figure3_report.txt
#
# 6 sections:
#   1. Volcano plots (6: two per appendage, GO-domain colored, unified axes)
#   2. DE summary bars (a: moderate/strong condition-colored, b: GO-domain-colored)
#   3. Venn diagrams (4: tissue overlap of state-upregulated genes)
#   4. Horizontal violin plots (LFC distributions of significant DEGs)
#   5. Within-appendage overlap (sex-dimorphic vs mating-responsive DEGs)
#   6. Supplementary CSVs (annotated DEG tables per contrast)
################################################################################

# ── Libraries ────────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(scales)
library(ggforce)
library(cowplot)
library(ggrepel)
library(jsonlite)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR <- "C:/Users/dorpe/Downloads/BSF2026/BSF2026/Supplemetary_Files/csv"
DE_DIR   <- file.path(BASE_DIR, "Sex_Reproductive_State _DE")
FIG_DIR  <- "C:/Users/dorpe/Downloads/BSF2026/BSF2026/Supplemetary_Files/Figure_3"

go_file       <- file.path(BASE_DIR, "BSF_all-rna_GO_ID_annotated.csv")
norm_file     <- file.path(BASE_DIR, "BSF_normalized_counts_master.csv")
chemo_id_file <- file.path(BASE_DIR, "BSF_Olfactory_ID.csv")

de_files <- list(
  list(tissue = "Antenna",        cond1 = "VF", cond2 = "Vm",
       file = file.path(DE_DIR, "results_ant_VF_vs_Vm.csv")),
  list(tissue = "Antenna",        cond1 = "MF", cond2 = "VF",
       file = file.path(DE_DIR, "results_ant_MF_vs_VF.csv")),
  list(tissue = "Maxillary palp", cond1 = "VF", cond2 = "Vm",
       file = file.path(DE_DIR, "results_palp_VF_vs_Vm.csv")),
  list(tissue = "Maxillary palp", cond1 = "MF", cond2 = "VF",
       file = file.path(DE_DIR, "results_palp_MF_vs_VF.csv")),
  list(tissue = "Tarsi",          cond1 = "VF", cond2 = "Vm",
       file = file.path(DE_DIR, "results_leg_VF_vs_Vm.csv")),
  list(tissue = "Tarsi",          cond1 = "MF", cond2 = "VF",
       file = file.path(DE_DIR, "results_leg_MF_vs_VF.csv"))
)

report_file <- file.path(FIG_DIR, "Figure3_report.txt")

# ── Thresholds ───────────────────────────────────────────────────────────────
PADJ_THR   <- 0.001
LFC_THR    <- 1.0
STRONG_LFC <- 2.5

# ── Full condition labels ────────────────────────────────────────────────────
COND_FULL <- c(MF = "Mated female", VF = "Virgin female", Vm = "Virgin male")

# ── Color palettes ───────────────────────────────────────────────────────────
COL_MF <- c(light = "#C77CFF", dark = "#7B1FA2")
COL_VF <- c(light = "#FF66B3", dark = "#99004D")
COL_VM <- c(light = "#99B3FF", dark = "#003399")

COND_COLORS <- list(MF = COL_MF, VF = COL_VF, Vm = COL_VM)

GO_COLS <- c(Chemosensory = "#A60000", MF = "#3FC498", CC = "#4D50DB",
             BP = "#A860E3", Unknown = "#666666")
PAL_GO  <- c(GO_COLS, `Not Significant` = "#D9D9D9")

DOMAIN_FULL <- c(Chemosensory = "Chemosensory", MF = "Molecular function",
                 CC = "Cellular component", BP = "Biological process",
                 Unknown = "Unknown")

# ── Chemosensory tags ────────────────────────────────────────────────────────
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

# ── Tissue / contrast ordering ───────────────────────────────────────────────
TISSUES_ORDER  <- c("Antenna", "Maxillary palp", "Tarsi")
CONTRAST_ORDER <- c("VF_vs_Vm", "MF_vs_VF")

# =============================================================================
# REPORT — open
# =============================================================================
rpt <- file(report_file, open = "wt")
write_rpt <- function(...) writeLines(paste0(...), rpt)

write_rpt("================================================================")
write_rpt("BSF Sex & Reproductive State DE — Figure 3 Report")
write_rpt("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
write_rpt("================================================================")
write_rpt("")
write_rpt("Thresholds:")
write_rpt("  padj < ", PADJ_THR)
write_rpt("  |log2FoldChange| >= ", LFC_THR)
write_rpt("  Strong |LFC| >= ", STRONG_LFC)
write_rpt("")

# =============================================================================
# READ & PREPARE DATA
# =============================================================================

# ── GO map ───────────────────────────────────────────────────────────────────
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

go_raw <- read.csv(go_file, stringsAsFactors = FALSE)
go_map <- data.frame(
  Gene      = trimws(as.character(go_raw[[1]])),
  GO_Name   = trimws(as.character(go_raw$GO_Name)),
  GO_Domain = clean_go_domain(go_raw$GO_Domain),
  stringsAsFactors = FALSE
)
go_map$GO_Name[is.na(go_map$GO_Name) | go_map$GO_Name == ""] <- "Unknown"
go_map <- go_map[!duplicated(go_map$Gene), ]

# ── Normalized counts ────────────────────────────────────────────────────────
norm <- read.csv(norm_file, stringsAsFactors = FALSE)
norm$JoinKey <- trimws(as.character(norm$Gene))

TISSUE_PREFIXES <- list(Antenna = "Ant_", `Maxillary palp` = "P_", Tarsi = "Leg_")
STATE_CODES     <- c("MF", "VF", "Vm")

for (tissue in names(TISSUE_PREFIXES)) {
  prefix <- TISSUE_PREFIXES[[tissue]]
  for (state in STATE_CODES) {
    pattern <- paste0("^", prefix, state, "[0-9]")
    cols <- grep(pattern, names(norm), value = TRUE)
    if (length(cols) > 0) {
      col_name <- paste0("mean_", gsub(" ", "_", tissue), "_", state)
      norm[[col_name]] <- rowMeans(norm[, cols, drop = FALSE], na.rm = TRUE)
    }
  }
}

mean_cols_all <- grep("^mean_", names(norm), value = TRUE)
norm_means    <- norm[, c("JoinKey", mean_cols_all), drop = FALSE]
norm_means    <- norm_means[!duplicated(norm_means$JoinKey), ]

norm_name_map <- if ("Name" %in% names(norm)) {
  setNames(trimws(as.character(norm$Name)), norm$JoinKey)
} else { NULL }

write_rpt("GO map: ", nrow(go_map), " unique genes")
write_rpt("Norm counts: ", nrow(norm), " genes, ", length(mean_cols_all), " mean columns")
write_rpt("")

# ── Chemosensory gene ID table ──────────────────────────────────────────────
chemo_id <- read.csv(chemo_id_file, stringsAsFactors = FALSE)
names(chemo_id) <- gsub("^﻿", "", names(chemo_id))
chemo_id$Name      <- trimws(as.character(chemo_id$Name))
chemo_id$Transcript <- trimws(as.character(chemo_id$Transcript))
chemo_id$NormKey   <- sub("^rna-", "", chemo_id$Transcript)

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

chemo_name_lookup <- setNames(chemo_id$Name, chemo_id$NormKey)
chemo_fam_lookup  <- setNames(chemo_id$Family, chemo_id$NormKey)

# ── Dmel BLAST lookup (from pre-computed S1 supplementary) ───────────────────
s1_path <- file.path(FIG_DIR, "S1_all_DEGs.csv")
dmel_lookup <- NULL
if (file.exists(s1_path)) {
  s1_raw <- read.csv(s1_path, stringsAsFactors = FALSE)
  blast_cols_avail <- intersect(
    c("Gene","Dmel_hit","Dmel_gene_symbol","Dmel_protein_description",
      "Dmel_evalue","Dmel_pct_aa_identity","Dmel_query_coverage_pct","Dmel_bitscore"),
    names(s1_raw))
  dmel_lookup <- s1_raw[!duplicated(s1_raw$Gene), blast_cols_avail]
  # Strip accession prefix (e.g. "NP_xxx.x| adult cuticle protein 1" -> "adult cuticle protein 1")
  dmel_lookup$Dmel_clean <- ifelse(
    grepl("\\|", dmel_lookup$Dmel_gene_symbol),
    trimws(sub("^[^|]+\\|\\s*", "", dmel_lookup$Dmel_gene_symbol)),
    dmel_lookup$Dmel_gene_symbol
  )
  dmel_lookup$Dmel_clean[dmel_lookup$Dmel_hit != "TRUE" |
                           dmel_lookup$Dmel_clean == ""] <- "No Dmel hit"
  cat("Dmel BLAST lookup: ", nrow(dmel_lookup), " genes\n", sep = "")
}

# ── Read & annotate all DE tables ────────────────────────────────────────────
# LFC is NEGATED on read (matching app convention):
#   positive = cond1 up, negative = cond2 up

read_de <- function(filepath) {
  df <- read.csv(filepath, stringsAsFactors = FALSE)
  df$Gene <- trimws(as.character(df$Gene))
  df$padj <- as.numeric(df$padj); df$padj[is.na(df$padj)] <- 1.0
  df$log2FoldChange <- -as.numeric(df$log2FoldChange)
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0.0
  if (!"Name" %in% names(df)) df$Name <- NA_character_
  df$Name <- trimws(as.character(df$Name))
  return(df)
}

annotate_de <- function(df, cond1, cond2) {
  df <- merge(df, go_map, by = "Gene", all.x = TRUE)
  df$GO_Domain[is.na(df$GO_Domain)] <- "Unknown"
  df$GO_Name[is.na(df$GO_Name) | df$GO_Name == ""] <- "Unknown"

  abs_lfc   <- abs(df$log2FoldChange)
  df$is_sig <- (df$padj < PADJ_THR) & (abs_lfc >= LFC_THR)

  df$Side <- NA_character_
  df$Side[df$is_sig & df$log2FoldChange >=  LFC_THR] <- cond1
  df$Side[df$is_sig & df$log2FoldChange <= -LFC_THR] <- cond2

  df$Strength <- NA_character_
  df$Strength[df$is_sig & abs_lfc >= STRONG_LFC] <- "Strong"
  df$Strength[df$is_sig & is.na(df$Strength)]    <- "Moderate"

  name_col <- ifelse(!is.na(df$Name) & df$Name != "", df$Name, df$Gene)
  df$ChemoName <- extract_chemo_tag(name_col)
  df$is_chemo  <- !is.na(df$ChemoName)

  df$padj_safe <- pmax(df$padj, .Machine$double.xmin)
  df$neg_log10 <- -log10(df$padj_safe)

  df$ColorKey <- ifelse(!df$is_sig, "Not Significant",
                 ifelse(df$is_chemo, "Chemosensory",
                        df$GO_Domain))
  return(df)
}

# Container for all 6 contrasts
all_contrasts <- list()
for (item in de_files) {
  df_raw <- read_de(item$file)
  df_ann <- annotate_de(df_raw, item$cond1, item$cond2)

  short_key <- paste0(item$tissue, " ", item$cond1, "_vs_", item$cond2)
  label     <- paste0(item$tissue, ": ", COND_FULL[item$cond1],
                      " vs ", COND_FULL[item$cond2])

  all_contrasts[[short_key]] <- list(
    tissue = item$tissue, cond1 = item$cond1, cond2 = item$cond2,
    label = label, short = short_key, df = df_ann
  )
}

write_rpt("Contrasts loaded:")
for (key in names(all_contrasts)) {
  it <- all_contrasts[[key]]
  n_sig <- sum(it$df$is_sig)
  n_c1  <- sum(it$df$Side == it$cond1, na.rm = TRUE)
  n_c2  <- sum(it$df$Side == it$cond2, na.rm = TRUE)
  n_ch  <- sum(it$df$is_sig & it$df$is_chemo)
  write_rpt("  ", it$label, ": ", nrow(it$df), " genes, ", n_sig, " DEGs (",
            COND_FULL[it$cond1], "=", n_c1, ", ", COND_FULL[it$cond2], "=", n_c2,
            "), ", n_ch, " chemosensory")
}
write_rpt("")

# =============================================================================
# GRAPHICS DEVICE — redirect print() output in non-interactive (Rscript) runs
# =============================================================================
PDF_DIR <- file.path(FIG_DIR, "Figure_3_PDFs")
dir.create(PDF_DIR, showWarnings = FALSE, recursive = TRUE)
if (!interactive()) {
  pdf(file.path(tempdir(), "fig3_preview_sink.pdf"), onefile = TRUE, width = 8, height = 6)
}

# =============================================================================
# 1. VOLCANO PLOTS (6 x 2: condition-colored + GO-colored, unified axes)
# =============================================================================

write_rpt("================================================================")
write_rpt("1. VOLCANO PLOTS")
write_rpt("================================================================")

VOLCANO_Y_MAX <- max(sapply(all_contrasts, function(it) {
  max(it$df$neg_log10, na.rm = TRUE)
})) * 1.05

volcano_base <- function(item, y_max) {
  df <- item$df; cond1 <- item$cond1; cond2 <- item$cond2
  n_pos  <- sum(df$is_sig & df$log2FoldChange > 0, na.rm = TRUE)
  n_neg  <- sum(df$is_sig & df$log2FoldChange < 0, na.rm = TRUE)
  ch_pos <- sum(df$is_sig & df$log2FoldChange > 0 & df$is_chemo, na.rm = TRUE)
  ch_neg <- sum(df$is_sig & df$log2FoldChange < 0 & df$is_chemo, na.rm = TRUE)
  title_text <- paste0(item$tissue, ":  ",
                       COND_FULL[cond2], " ", ch_neg, "/", n_neg,
                       "    ",
                       COND_FULL[cond1], " ", ch_pos, "/", n_pos)
  list(df = df, cond1 = cond1, cond2 = cond2, title = title_text)
}

make_volcano_condition <- function(item, y_max) {
  v <- volcano_base(item, y_max); df <- v$df
  col1 <- COND_COLORS[[v$cond1]]; col2 <- COND_COLORS[[v$cond2]]
  df$CondColor <- ifelse(!df$is_sig, "NS",
                  ifelse(df$log2FoldChange > 0,
                    ifelse(abs(df$log2FoldChange) >= STRONG_LFC,
                           paste0(v$cond1, "_Strong"), paste0(v$cond1, "_Moderate")),
                    ifelse(abs(df$log2FoldChange) >= STRONG_LFC,
                           paste0(v$cond2, "_Strong"), paste0(v$cond2, "_Moderate"))))
  df$CondColor <- factor(df$CondColor,
    levels = c("NS", paste0(v$cond1, "_Moderate"), paste0(v$cond1, "_Strong"),
               paste0(v$cond2, "_Moderate"), paste0(v$cond2, "_Strong")))
  df <- df[order(df$CondColor), ]
  fill_vals <- c(NS = "#D9D9D9",
    setNames(col1["light"], paste0(v$cond1, "_Moderate")),
    setNames(col1["dark"],  paste0(v$cond1, "_Strong")),
    setNames(col2["light"], paste0(v$cond2, "_Moderate")),
    setNames(col2["dark"],  paste0(v$cond2, "_Strong")))
  fill_labels <- c("Not Significant",
    paste(COND_FULL[v$cond1], "Moderate"), paste(COND_FULL[v$cond1], "Strong"),
    paste(COND_FULL[v$cond2], "Moderate"), paste(COND_FULL[v$cond2], "Strong"))
  names(fill_labels) <- names(fill_vals)

  ggplot(df, aes(x = log2FoldChange, y = neg_log10, color = CondColor)) +
    geom_point(size = 5, alpha = 0.5) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = fill_vals, labels = fill_labels, name = NULL) +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = v$title,
         x = expression(Log[2](FoldChange)),
         y = expression(-log[10](italic(padj)))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = "bottom", legend.text = element_text(size = 8),
          panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

make_volcano_go <- function(item, y_max) {
  v <- volcano_base(item, y_max); df <- v$df
  df$ColorKey <- factor(df$ColorKey,
    levels = c("Not Significant", "MF", "CC", "BP", "Unknown", "Chemosensory"))
  df <- df[order(df$ColorKey), ]

  ggplot(df, aes(x = log2FoldChange, y = neg_log10, color = ColorKey)) +
    geom_point(size = 5, alpha = 0.5) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(
      values = PAL_GO,
      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown", "Not Significant"),
      labels = c("Chemosensory", "Molecular function", "Cellular component",
                 "Biological process", "Unknown", "Not Significant"),
      na.value = "#D9D9D9") +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = v$title,
         x = expression(Log[2](FoldChange)),
         y = expression(-log[10](italic(padj)))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          legend.position = "bottom", legend.text = element_text(size = 8),
          panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
}

# ── 1a. Condition-colored volcanos ──────────────────────────────────────────
write_rpt("  1a. Condition-colored volcanos (Moderate / Strong LFC):")
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      print(make_volcano_condition(it, VOLCANO_Y_MAX))
      n_sig <- sum(it$df$is_sig); n_ch <- sum(it$df$is_sig & it$df$is_chemo)
      write_rpt("  ", it$label, ": ", n_sig, " DEGs (", n_ch, " chemosensory)")
    }
  }
}
write_rpt("")

# ── 1b. GO-domain-colored volcanos ─────────────────────────────────────────
write_rpt("  1b. GO-domain-colored volcanos:")
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      print(make_volcano_go(it, VOLCANO_Y_MAX))
    }
  }
}
write_rpt("")
cat("Done: Volcano plots\n")


# =============================================================================
# 2. DE SUMMARY BARS (moderate / strong)
# =============================================================================

write_rpt("================================================================")
write_rpt("2. DE SUMMARY BARS")
write_rpt("================================================================")

# Bar x-axis order:
#   VF vs Vm → Vm left, VF right (male first)
#   MF vs VF → VF left, MF right (virgin first)
bar_x_order <- function(cond1, cond2) {
  if (setequal(c(cond1, cond2), c("VF", "Vm")))
    return(c(COND_FULL["Vm"], COND_FULL["VF"]))
  if (cond2 == "VF") return(c(COND_FULL["VF"], COND_FULL[cond1]))
  if (cond1 == "VF") return(c(COND_FULL["VF"], COND_FULL[cond2]))
  return(c(COND_FULL[cond1], COND_FULL[cond2]))
}

# Pre-compute global y-axis max for condition-colored bars
BAR_Y_MAX_COND <- 0
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      sig <- all_contrasts[[key]]$df[all_contrasts[[key]]$df$is_sig, ]
      for (side in unique(sig$Side)) {
        n_side <- sum(sig$Side == side, na.rm = TRUE)
        BAR_Y_MAX_COND <- max(BAR_Y_MAX_COND, n_side)
      }
    }
  }
}
BAR_Y_MAX_COND <- BAR_Y_MAX_COND * 1.08

make_de_bar <- function(item, y_max) {
  df <- item$df; cond1 <- item$cond1; cond2 <- item$cond2

  sig <- df[df$is_sig & !is.na(df$Side) & !is.na(df$Strength), ]
  if (nrow(sig) == 0) return(NULL)

  all_combos <- expand.grid(Side = c(cond1, cond2),
                             Strength = c("Moderate", "Strong"),
                             stringsAsFactors = FALSE)
  counts <- sig %>%
    group_by(Side, Strength) %>%
    summarise(N = n(), .groups = "drop")
  counts <- merge(all_combos, counts, all.x = TRUE)
  counts$N[is.na(counts$N)] <- 0

  total <- nrow(df)
  counts$Pct   <- 100 * counts$N / total
  counts$Label <- ifelse(counts$N > 0,
                         paste0(counts$N, " (", sprintf("%.1f%%", counts$Pct), ")"), "")
  counts$Side_full <- COND_FULL[counts$Side]
  counts$FillKey   <- paste0(counts$Side, " ", counts$Strength)

  col1 <- COND_COLORS[[cond1]]; col2 <- COND_COLORS[[cond2]]
  fill_map <- setNames(
    c(col1["light"], col1["dark"], col2["light"], col2["dark"]),
    c(paste(cond1, "Moderate"), paste(cond1, "Strong"),
      paste(cond2, "Moderate"), paste(cond2, "Strong"))
  )
  label_map <- setNames(
    c(paste(COND_FULL[cond1], "Moderate"), paste(COND_FULL[cond1], "Strong"),
      paste(COND_FULL[cond2], "Moderate"), paste(COND_FULL[cond2], "Strong")),
    c(paste(cond1, "Moderate"), paste(cond1, "Strong"),
      paste(cond2, "Moderate"), paste(cond2, "Strong"))
  )

  x_order <- bar_x_order(cond1, cond2)

  ggplot(counts, aes(x = factor(Side_full, levels = x_order), y = N, fill = FillKey)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65) +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = fill_map, labels = label_map) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0))) +
    labs(title = paste0(item$tissue, ": ", COND_FULL[cond1], " vs ", COND_FULL[cond2]),
         x = NULL, y = "Number of significant DEGs", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "bottom", legend.text = element_text(size = 8))
}

# ── 2a. Condition-colored bars (moderate / strong) ──────────────────────────
write_rpt("  2a. Condition-colored bars (moderate / strong):")
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      p <- make_de_bar(it, BAR_Y_MAX_COND)
      if (!is.null(p)) print(p)
      sig <- it$df[it$df$is_sig, ]
      n_mod <- sum(sig$Strength == "Moderate", na.rm = TRUE)
      n_str <- sum(sig$Strength == "Strong",   na.rm = TRUE)
      write_rpt("  ", it$label, ": Moderate=", n_mod, ", Strong=", n_str)
    }
  }
}
write_rpt("")

# ── 2b. GO-domain-colored bars ──────────────────────────────────────────────
write_rpt("  2b. GO-domain-colored bars:")

make_go_de_bar <- function(item) {
  df <- item$df; cond1 <- item$cond1; cond2 <- item$cond2
  sig <- df[df$is_sig, ]
  if (nrow(sig) == 0) return(NULL)

  sig$Dir <- ifelse(sig$log2FoldChange > 0,
                    paste0(COND_FULL[cond1], " up"),
                    paste0(COND_FULL[cond2], " up"))
  sig$PlotDomain <- ifelse(sig$is_chemo, "Chemosensory", sig$GO_Domain)
  sig$PlotDomain <- factor(sig$PlotDomain,
                           levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))

  tbl <- sig %>%
    group_by(Dir, PlotDomain, .drop = FALSE) %>%
    summarise(N = n(), .groups = "drop") %>%
    group_by(Dir) %>%
    mutate(Percent = 100 * N / sum(N)) %>%
    ungroup() %>%
    mutate(TextLabel = ifelse(Percent >= 2,
                              paste0("N=", N, "\n", sprintf("%.1f%%", Percent)),
                              ""))

  # VF vs Vm → Vm up left, VF up right; MF vs VF → VF up left, MF up right
  if (setequal(c(cond1, cond2), c("VF", "Vm"))) {
    dir_order <- c(paste0(COND_FULL["Vm"], " up"),
                   paste0(COND_FULL["VF"], " up"))
  } else if (cond2 == "VF") {
    dir_order <- c(paste0(COND_FULL["VF"], " up"),
                   paste0(COND_FULL[cond1], " up"))
  } else {
    dir_order <- c(paste0(COND_FULL["VF"], " up"),
                   paste0(COND_FULL[cond2], " up"))
  }

  ggplot(tbl, aes(x = factor(Dir, levels = dir_order), y = Percent, fill = PlotDomain)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65) +
    geom_text(aes(label = TextLabel), position = position_stack(vjust = 0.5),
              size = 2.5, color = "white", fontface = "bold", lineheight = 0.85) +
    scale_fill_manual(values = GO_COLS, labels = DOMAIN_FULL,
                      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown")) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
    labs(title = paste0(item$tissue, ": ", COND_FULL[cond1], " vs ", COND_FULL[cond2]),
         x = NULL, y = "% of significant DEGs", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "bottom", legend.text = element_text(size = 8),
          axis.text.x = element_text(size = 9))
}

for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      p <- make_go_de_bar(it)
      if (!is.null(p)) print(p)

      sig <- it$df[it$df$is_sig, ]
      sig$PlotDomain <- ifelse(sig$is_chemo, "Chemosensory", sig$GO_Domain)
      for (dom in c("Chemosensory", "MF", "CC", "BP", "Unknown")) {
        write_rpt("    ", it$short, " ", dom, ": ", sum(sig$PlotDomain == dom))
      }
    }
  }
}
write_rpt("")
cat("Done: DE summary bars\n")


# =============================================================================
# 3. VENN DIAGRAMS (4: tissue overlap of state-upregulated genes)
# =============================================================================

write_rpt("================================================================")
write_rpt("3. VENN DIAGRAMS")
write_rpt("================================================================")
write_rpt("Overlap of upregulated genes across appendages per state.")
write_rpt("")

get_up_genes <- function(contrast_suffix, side) {
  sets <- list()
  for (tissue in TISSUES_ORDER) {
    key <- paste0(tissue, " ", contrast_suffix)
    if (is.null(all_contrasts[[key]])) next
    df <- all_contrasts[[key]]$df
    if (side == "cond1") {
      genes <- df$Gene[df$padj < PADJ_THR & df$log2FoldChange >= LFC_THR]
    } else {
      genes <- df$Gene[df$padj < PADJ_THR & df$log2FoldChange <= -LFC_THR]
    }
    sets[[tissue]] <- unique(genes[!is.na(genes)])
  }
  return(sets)
}

compute_venn_regions <- function(sA, sB, sC) {
  list(
    `A only`    = setdiff(setdiff(sA, sB), sC),
    `B only`    = setdiff(setdiff(sB, sA), sC),
    `C only`    = setdiff(setdiff(sC, sA), sB),
    `A & B`     = setdiff(intersect(sA, sB), sC),
    `A & C`     = setdiff(intersect(sA, sC), sB),
    `B & C`     = setdiff(intersect(sB, sC), sA),
    `A & B & C` = Reduce(intersect, list(sA, sB, sC))
  )
}

draw_venn_3way <- function(title, color, labels, regions) {
  counts <- sapply(regions, length)
  circles <- data.frame(x0 = c(0, -1.5, 1.5), y0 = c(1.8, 0, 0),
                         r = c(1.5, 1.5, 1.5))
  fill_col <- adjustcolor(color, alpha.f = 0.2)

  num_df <- data.frame(
    x   = c(0, -2.0, 2.0, -1.0, 1.0, 0, 0),
    y   = c(3.0, -0.2, -0.2, 1.2, 1.2, -0.4, 0.9),
    txt = as.character(counts[c("A only", "B only", "C only",
                                 "A & B", "A & C", "B & C", "A & B & C")]),
    stringsAsFactors = FALSE
  )
  lbl_df <- data.frame(
    x = c(0, -2.5, 2.5), y = c(3.8, -1.8, -1.8), txt = labels,
    stringsAsFactors = FALSE
  )

  ggplot() +
    ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r),
                          fill = fill_col, color = NA, linewidth = 0) +
    geom_text(data = num_df, aes(x = x, y = y, label = txt), size = 5) +
    geom_text(data = lbl_df, aes(x = x, y = y, label = txt),
              size = 4, color = color, fontface = "bold") +
    coord_fixed(xlim = c(-4, 4), ylim = c(-2.5, 4.5)) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                                     color = color))
}

venn_configs <- list(
  list(title = "Virgin male upregulated",
       contrast = "VF_vs_Vm", side = "cond2", color = COL_VM["dark"]),
  list(title = "Virgin female upregulated (vs Vm)",
       contrast = "VF_vs_Vm", side = "cond1", color = COL_VF["dark"]),
  list(title = "Virgin female upregulated (vs MF)",
       contrast = "MF_vs_VF", side = "cond2", color = COL_VF["dark"]),
  list(title = "Mated female upregulated",
       contrast = "MF_vs_VF", side = "cond1", color = COL_MF["dark"])
)

tissue_labels_venn <- c("Antenna", "Maxillary palp", "Tarsi")

for (vc in venn_configs) {
  sets <- get_up_genes(vc$contrast, vc$side)
  sA <- if (!is.null(sets[["Antenna"]]))        sets[["Antenna"]]        else character(0)
  sB <- if (!is.null(sets[["Maxillary palp"]])) sets[["Maxillary palp"]] else character(0)
  sC <- if (!is.null(sets[["Tarsi"]]))          sets[["Tarsi"]]          else character(0)

  regions <- compute_venn_regions(sA, sB, sC)
  print(draw_venn_3way(vc$title, vc$color, tissue_labels_venn, regions))

  write_rpt("  ", vc$title, ":")
  write_rpt("    Antenna=", length(sA), ", Maxillary palp=", length(sB),
            ", Tarsi=", length(sC))
  write_rpt("    A only=", length(regions[["A only"]]),
            ", B only=", length(regions[["B only"]]),
            ", C only=", length(regions[["C only"]]))
  write_rpt("    A&B=", length(regions[["A & B"]]),
            ", A&C=", length(regions[["A & C"]]),
            ", B&C=", length(regions[["B & C"]]),
            ", A&B&C=", length(regions[["A & B & C"]]))
  write_rpt("")
}
cat("Done: Venn diagrams\n")


# =============================================================================
# 4. HORIZONTAL VIOLIN PLOTS (LFC distributions)
# =============================================================================

write_rpt("================================================================")
write_rpt("4. VIOLIN PLOTS — LFC DISTRIBUTIONS")
write_rpt("================================================================")

# Build background (all genes) and foreground (sig only) data
violin_bg <- do.call(rbind, lapply(names(all_contrasts), function(key) {
  it <- all_contrasts[[key]]
  data.frame(LFC = it$df$log2FoldChange,
             Tissue = it$tissue,
             Contrast = paste0(COND_FULL[it$cond1], " vs\n", COND_FULL[it$cond2]),
             stringsAsFactors = FALSE)
}))

violin_fg <- do.call(rbind, lapply(names(all_contrasts), function(key) {
  it <- all_contrasts[[key]]
  sig <- it$df[it$df$is_sig, ]
  if (nrow(sig) == 0) return(NULL)
  data.frame(LFC = sig$log2FoldChange,
             Tissue = it$tissue,
             Contrast = paste0(COND_FULL[it$cond1], " vs\n", COND_FULL[it$cond2]),
             Up_in = COND_FULL[sig$Side],
             stringsAsFactors = FALSE)
}))

cond_fill <- setNames(
  c(COL_MF["dark"], COL_VF["dark"], COL_VM["dark"]),
  c("Mated female", "Virgin female", "Virgin male")
)

for (tissue in TISSUES_ORDER) {
  bg_sub <- violin_bg[violin_bg$Tissue == tissue, ]
  fg_sub <- violin_fg[violin_fg$Tissue == tissue, ]
  if (nrow(fg_sub) == 0) next

  # Count DEGs per contrast x direction for annotation
  jitter_counts <- fg_sub %>%
    group_by(Contrast, Up_in) %>%
    summarise(N = n(), .groups = "drop")
  jitter_counts$x_pos <- ifelse(
    jitter_counts$Up_in %in% c("Virgin female", "Mated female"), 8, -8)
  jitter_counts$Label <- paste0("n=", jitter_counts$N)

  p <- ggplot() +
    geom_violin(data = bg_sub, aes(y = Contrast, x = LFC),
                fill = "grey90", color = "grey70", alpha = 0.5,
                scale = "width", width = 0.8) +
    geom_jitter(data = fg_sub, aes(y = Contrast, x = LFC, color = Up_in),
                height = 0.15, size = 1, alpha = 0.6) +
    geom_text(data = jitter_counts,
              aes(y = Contrast, x = x_pos, label = Label, color = Up_in),
              size = 3, fontface = "bold", show.legend = FALSE) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
               linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = cond_fill, name = "Upregulated in") +
    scale_x_continuous(limits = c(-10, 10)) +
    labs(title = tissue,
         x = expression(Log[2](FoldChange)), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.position = "bottom", panel.grid.minor = element_blank())

  print(p)

  for (contr in unique(fg_sub$Contrast)) {
    sub <- fg_sub[fg_sub$Contrast == contr, ]
    write_rpt("  ", tissue, " — ", gsub("\n", " ", contr), ": ",
              nrow(sub), " DEGs, median LFC=", sprintf("%.2f", median(sub$LFC)),
              " [", sprintf("%.1f", min(sub$LFC)), ", ",
              sprintf("%.1f", max(sub$LFC)), "]")
  }
}
write_rpt("")
cat("Done: Violin plots\n")


# =============================================================================
# 5. SEX-DIMORPHIC & MATING-RESPONSIVE OVERLAP — detailed analysis + CSVs
# =============================================================================
#
# Sex-dimorphic genes: significantly DE between virgin females (VF) and virgin
# males (Vm). These reflect baseline transcriptional differences between sexes
# independent of mating — e.g. genes that underlie sex-specific olfactory
# tuning, pheromone detection, or oviposition-related chemosensation.
#
# Mating-responsive genes: significantly DE between mated females (MF) and
# virgin females (VF). These change after mating — potentially shifting the
# female chemosensory profile toward host-seeking for oviposition.
#
# Sex-and-mating-regulated genes (overlap): DE in BOTH contrasts.
#   - Concordant (same LFC sign): mating reinforces the sex difference.
#     E.g. a female-biased gene further upregulated after mating → the mated
#     female's expression diverges even more from males.
#   - Discordant (opposite LFC sign): mating reverses or attenuates the sex
#     difference. E.g. a gene higher in VF than Vm that drops after mating →
#     the mated female's profile shifts toward a more male-like state for
#     that gene, potentially reflecting a post-reproductive sensory switch.
#
# Ranking: genes are ranked by |LFC_sex| (= |LFC in VF vs Vm|), capturing
# how strongly a gene differs between the sexes. Genes with the largest
# absolute sex effect are the most sexually dimorphic in that appendage.
#
# Connection to violin plots (Section 4): each significant DEG shown in the
# VF vs Vm violin is a sex-dimorphic gene. Sex-and-mating-regulated genes also appear
# in the MF vs VF violin for the same appendage.

write_rpt("================================================================")
write_rpt("5. SEX-DIMORPHIC & MATING-RESPONSIVE OVERLAP")
write_rpt("================================================================")
write_rpt("")
write_rpt("Sex-dimorphic = DE in VF vs Vm (baseline sex difference).")
write_rpt("Mating-responsive = DE in MF vs VF (expression change after mating).")
write_rpt("Sex-and-mating-regulated = DE in both contrasts.")
write_rpt("  Direction reinforced = same LFC direction (mating reinforces sex bias).")
write_rpt("  Direction reversed = opposite LFC direction (mating reverses or attenuates sex bias).")
write_rpt("Ranked by |LFC_sex| — largest absolute sex effect first.")
write_rpt("")

# Shorten a GO name to first 3 words (for scatter labels)
shorten_go <- function(x) {
  sapply(as.character(x), function(s) {
    if (is.na(s) || s == "" || s == "Unknown") return("")
    words <- strsplit(s, "\\s+")[[1]]
    if (length(words) <= 3) s else paste(words[1:3], collapse = " ")
  }, USE.NAMES = FALSE)
}

# Pre-compute global symmetric axis limit for scatter plots across all tissues
SCATTER_AXIS_LIM <- 0
for (tissue in TISSUES_ORDER) {
  key_sex  <- paste0(tissue, " VF_vs_Vm")
  key_mate <- paste0(tissue, " MF_vs_VF")
  if (is.null(all_contrasts[[key_sex]]) || is.null(all_contrasts[[key_mate]])) next
  df_sex   <- all_contrasts[[key_sex]]$df
  df_mate  <- all_contrasts[[key_mate]]$df
  sig_sex  <- df_sex[df_sex$is_sig, ]
  sig_mate <- df_mate[df_mate$is_sig, ]
  both_g   <- intersect(sig_sex$Gene, sig_mate$Gene)
  if (length(both_g) > 0) {
    lfc_sx <- sig_sex$log2FoldChange[match(both_g, sig_sex$Gene)]
    lfc_mt <- df_mate$log2FoldChange[match(both_g, df_mate$Gene)]
    SCATTER_AXIS_LIM <- max(SCATTER_AXIS_LIM, abs(lfc_sx), abs(lfc_mt), na.rm = TRUE)
  }
}
SCATTER_AXIS_LIM <- ceiling(SCATTER_AXIS_LIM * 1.15)

s2_list <- list()  # accumulate per-tissue sex/mating classification rows

for (tissue in TISSUES_ORDER) {
  key_sex  <- paste0(tissue, " VF_vs_Vm")
  key_mate <- paste0(tissue, " MF_vs_VF")
  if (is.null(all_contrasts[[key_sex]]) || is.null(all_contrasts[[key_mate]])) next

  df_sex   <- all_contrasts[[key_sex]]$df
  df_mate  <- all_contrasts[[key_mate]]$df
  sig_sex  <- df_sex[df_sex$is_sig, ]
  sig_mate <- df_mate[df_mate$is_sig, ]

  genes_sex  <- unique(sig_sex$Gene)
  genes_mate <- unique(sig_mate$Gene)
  both       <- intersect(genes_sex, genes_mate)
  only_sex   <- setdiff(genes_sex, genes_mate)
  only_mate  <- setdiff(genes_mate, genes_sex)

  # ── Venn overview (simple counts) ──────────────────────────────────────────
  circles <- data.frame(x0 = c(-0.5, 0.5), y0 = c(0, 0), r = c(1, 1))
  lbl_df <- data.frame(
    x   = c(-0.85, 0.85, 0, -0.5, 0.5),
    y   = c(0, 0, 0.15, -1.0, -1.0),
    txt = c(length(only_sex), length(only_mate), length(both),
            "Sex-dimorphic\nonly", "Mating-regulated\nonly"),
    sz  = c(5, 5, 5, 3.5, 3.5), stringsAsFactors = FALSE
  )
  venn_annot <- data.frame(
    x   = c(-1.6, 0, 1.6),
    y   = c(1.3, -0.55, 1.3),
    txt = c(
      "DE between VF & Vm\nbut NOT after mating\n(constitutive sex\ndifference)",
      "DE in BOTH contrasts\n(mating modulates a\nsex-dimorphic gene)",
      "DE between MF & VF\nbut NOT between sexes\n(mating-specific\nchange)"
    ), stringsAsFactors = FALSE
  )

  p <- ggplot() +
    ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r),
                          fill = adjustcolor("#555555", alpha.f = 0.15), color = NA) +
    geom_text(data = lbl_df, aes(x = x, y = y, label = txt), size = lbl_df$sz) +
    geom_text(data = venn_annot, aes(x = x, y = y, label = txt),
              size = 2.3, color = "grey30", lineheight = 0.9) +
    coord_fixed(xlim = c(-2.5, 2.5), ylim = c(-1.8, 2.0)) +
    labs(title = paste0(tissue, " — sex-dimorphic vs mating-responsive DEGs"),
         subtitle = paste0(
           "Left circle: DE in VF vs Vm (sex-dimorphic, baseline sex difference)    ",
           "Right circle: DE in MF vs VF (mating-responsive, post-mating change)")) +
    theme_void() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 7, color = "grey40"))
  print(p)

  # ── GO-domain breakdown: assign domain with chemosensory priority ──────────
  assign_go_domain <- function(gene_vec) {
    dom <- go_map$GO_Domain[match(gene_vec, go_map$Gene)]
    dom[is.na(dom)] <- "Unknown"
    dom[!is.na(chemo_fam_lookup[gene_vec])] <- "Chemosensory"
    return(dom)
  }

  venn_go_data <- data.frame(
    Gene = c(only_sex, both, only_mate),
    Region = c(rep("Sex-dimorphic only", length(only_sex)),
               rep("Sex-and-mating-regulated", length(both)),
               rep("Mating-regulated only", length(only_mate))),
    stringsAsFactors = FALSE
  )
  venn_go_data$GO_Domain <- assign_go_domain(venn_go_data$Gene)
  venn_go_data$GO_Domain <- factor(venn_go_data$GO_Domain,
    levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))
  venn_go_data$Region <- factor(venn_go_data$Region,
    levels = c("Sex-dimorphic only", "Sex-and-mating-regulated", "Mating-regulated only"))

  if (!is.null(norm_name_map)) {
    venn_go_data$Label <- norm_name_map[venn_go_data$Gene]
    venn_go_data$Label[is.na(venn_go_data$Label) | venn_go_data$Label == ""] <-
      venn_go_data$Gene[is.na(venn_go_data$Label) | venn_go_data$Label == ""]
  } else {
    venn_go_data$Label <- venn_go_data$Gene
  }

  # ── Per-domain Venns (sex-dimorphic vs mating-responsive, split by domain) ──
  for (dom in c("Chemosensory", "MF", "CC", "BP", "Unknown")) {
    dom_sex_genes  <- venn_go_data$Gene[
      venn_go_data$Region != "Mating-regulated only" &
      venn_go_data$GO_Domain == dom]
    dom_mate_genes <- venn_go_data$Gene[
      venn_go_data$Region != "Sex-dimorphic only" &
      venn_go_data$GO_Domain == dom]
    n_only_sex  <- length(setdiff(dom_sex_genes, dom_mate_genes))
    n_only_mate <- length(setdiff(dom_mate_genes, dom_sex_genes))
    n_both      <- length(intersect(dom_sex_genes, dom_mate_genes))

    circles_d <- data.frame(x0 = c(-0.5, 0.5), y0 = c(0, 0), r = c(1, 1))
    lbl_d <- data.frame(
      x   = c(-0.85, 0.85, 0, -0.5, 0.5),
      y   = c(0, 0, 0.15, -0.9, -0.9),
      txt = c(as.character(n_only_sex), as.character(n_only_mate),
              as.character(n_both), "Sex-dimorphic", "Mating-regulated"),
      sz  = c(4, 4, 5, 3.5, 3.5), stringsAsFactors = FALSE
    )
    dom_label <- if (dom == "Chemosensory") "Chemosensory" else DOMAIN_FULL[dom]
    p_dom <- ggplot() +
      ggforce::geom_circle(data = circles_d, aes(x0 = x0, y0 = y0, r = r),
        fill = adjustcolor(GO_COLS[dom], alpha.f = 0.25), color = NA) +
      geom_text(data = lbl_d, aes(x = x, y = y, label = txt), size = lbl_d$sz) +
      coord_fixed(xlim = c(-2, 2), ylim = c(-1.5, 1.5)) +
      labs(title = paste0(tissue, " — ", dom_label)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10,
                                       color = GO_COLS[dom]))
    print(p_dom)
    write_rpt("  Venn ", dom, ": sex-only=", n_only_sex,
              ", mating-only=", n_only_mate, ", both=", n_both)
  }

  # ── BLAST characterization of Unknown-domain genes ────────────────────────
  if (!is.null(dmel_lookup)) {
    unk_venn <- venn_go_data[venn_go_data$GO_Domain == "Unknown", ]
    if (nrow(unk_venn) > 0) {
      unk_venn <- merge(unk_venn,
                        dmel_lookup[, c("Gene", "Dmel_hit", "Dmel_clean")],
                        by = "Gene", all.x = TRUE)
      unk_venn$Dmel_clean[is.na(unk_venn$Dmel_clean)] <- "No Dmel hit"

      n_unk    <- nrow(unk_venn)
      n_hit    <- sum(unk_venn$Dmel_hit == "TRUE", na.rm = TRUE)
      n_no_hit <- n_unk - n_hit
      write_rpt("  Unknown-domain genes (", tissue, "): N=", n_unk,
                " | Dmel hit: ", n_hit, " | No Dmel hit: ", n_no_hit)
      for (reg in levels(unk_venn$Region)) {
        reg_sub <- unk_venn[unk_venn$Region == reg, ]
        if (nrow(reg_sub) == 0) next
        write_rpt("    ", reg, " (N=", nrow(reg_sub), ", with Dmel hit: ",
                  sum(reg_sub$Dmel_hit == "TRUE", na.rm = TRUE), "):")
        top_h <- sort(
          table(reg_sub$Dmel_clean[reg_sub$Dmel_clean != "No Dmel hit"]),
          decreasing = TRUE)
        for (i in seq_len(min(10, length(top_h)))) {
          write_rpt("      ", names(top_h)[i], ": N=", top_h[i])
        }
      }

      # Horizontal bar: top Dmel hits across Unknown genes, coloured by region
      blast_names <- unk_venn$Dmel_clean[unk_venn$Dmel_clean != "No Dmel hit"]
      if (length(blast_names) > 0) {
        top_names <- names(sort(table(blast_names), decreasing = TRUE))[
          seq_len(min(15, length(unique(blast_names))))]
        unk_bar <- unk_venn[unk_venn$Dmel_clean %in% top_names, ] %>%
          group_by(Region, Dmel_clean) %>%
          summarise(N = n(), .groups = "drop")
        unk_bar$Dmel_clean <- factor(
          substr(unk_bar$Dmel_clean, 1, 55),
          levels = rev(substr(top_names, 1, 55)))
        region_cols <- c(
          `Sex-dimorphic only`     = "#003399",
          `Sex-and-mating-regulated` = "#7B1FA2",
          `Mating-regulated only`   = "#99004D")
        p_unk <- ggplot(unk_bar,
                        aes(x = N, y = Dmel_clean, fill = Region)) +
          geom_bar(stat = "identity", position = "dodge", width = 0.7) +
          scale_fill_manual(values = region_cols) +
          labs(
            title    = paste0(tissue, " — Dmel BLAST hits: Unknown-domain DEGs"),
            subtitle = paste0("N=", n_unk, " Unknown | ",
                              n_hit, " with Dmel hit | ",
                              n_no_hit, " no hit"),
            x = "Number of genes", y = NULL, fill = NULL) +
          theme_minimal(base_size = 10) +
          theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 10),
            plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"),
            legend.position = "bottom",
            axis.text.y  = element_text(size = 8)
          )
        print(p_unk)
      }
    }
  }

  # ── GO-domain stacked bar summary ──────────────────────────────────────────
  venn_go_counts <- venn_go_data %>%
    group_by(Region, GO_Domain, .drop = FALSE) %>%
    summarise(N = n(),
              Names = paste(sort(Label), collapse = ", "),
              .groups = "drop") %>%
    group_by(Region) %>%
    mutate(Total = sum(N)) %>%
    ungroup()
  venn_go_counts$NLabel <- ifelse(venn_go_counts$N > 0,
    paste0(venn_go_counts$N), "")

  p_go_venn <- ggplot(venn_go_counts,
      aes(x = Region, y = N, fill = GO_Domain)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(aes(label = NLabel),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = GO_COLS, labels = DOMAIN_FULL,
                      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown")) +
    labs(title = paste0(tissue, " — GO domain breakdown"),
         subtitle = paste0("Sex-dimorphic: ", length(only_sex),
                           "  |  Dual: ", length(both),
                           "  |  Mating only: ", length(only_mate)),
         x = NULL, y = "Number of DEGs", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
          legend.position = "bottom",
          axis.text.x = element_text(size = 9))
  print(p_go_venn)

  # Report chemosensory gene names per region
  for (reg in levels(venn_go_data$Region)) {
    chemo_in_reg <- venn_go_data$Label[
      venn_go_data$Region == reg & venn_go_data$GO_Domain == "Chemosensory"]
    if (length(chemo_in_reg) > 0) {
      write_rpt("  Chemosensory in ", reg, ": ",
                paste(sort(chemo_in_reg), collapse = ", "))
    }
  }

  # ── Save Venn supplementary CSV ──────────────────────────────────────────
  # Enrich venn_go_data with GO_Name and LFC from both contrasts
  venn_csv_data <- venn_go_data
  venn_csv_data$GO_Name <- go_map$GO_Name[match(venn_csv_data$Gene, go_map$Gene)]
  venn_csv_data$GO_Name[is.na(venn_csv_data$GO_Name)] <- "Unknown"
  venn_csv_data$LFC_sex  <- sig_sex$log2FoldChange[match(venn_csv_data$Gene, sig_sex$Gene)]
  venn_csv_data$padj_sex <- sig_sex$padj[match(venn_csv_data$Gene, sig_sex$Gene)]
  venn_csv_data$LFC_mating  <- df_mate$log2FoldChange[match(venn_csv_data$Gene, df_mate$Gene)]
  venn_csv_data$padj_mating <- df_mate$padj[match(venn_csv_data$Gene, df_mate$Gene)]
  venn_csv_data$Chemosensory_family <- chemo_fam_lookup[venn_csv_data$Gene]
  venn_csv_data$Name <- venn_csv_data$Label

  # ── Build full annotated table of ALL sex-dimorphic genes ──────────────────
  # Start from every gene that is DE in VF vs Vm (sex-dimorphic)
  sex_tbl <- data.frame(Gene = genes_sex, stringsAsFactors = FALSE)

  # Name
  if (!is.null(norm_name_map)) sex_tbl$Name <- norm_name_map[sex_tbl$Gene]

  # GO annotation
  sex_tbl <- merge(sex_tbl, go_map, by = "Gene", all.x = TRUE)
  sex_tbl$GO_Name[is.na(sex_tbl$GO_Name)]     <- "Unknown"
  sex_tbl$GO_Domain[is.na(sex_tbl$GO_Domain)] <- "Unknown"

  # Sex effect (VF vs Vm)
  sex_lu_lfc  <- setNames(sig_sex$log2FoldChange, sig_sex$Gene)
  sex_lu_padj <- setNames(sig_sex$padj, sig_sex$Gene)
  sex_lu_side <- setNames(sig_sex$Side, sig_sex$Gene)
  sex_tbl$LFC_sex  <- sex_lu_lfc[sex_tbl$Gene]
  sex_tbl$padj_sex <- sex_lu_padj[sex_tbl$Gene]
  sex_tbl$Sex_bias <- ifelse(sex_tbl$LFC_sex > 0, "Female-biased", "Male-biased")

  # Mating effect (MF vs VF) — for ALL genes, not just significant
  mate_lu_lfc  <- setNames(df_mate$log2FoldChange, df_mate$Gene)
  mate_lu_padj <- setNames(df_mate$padj, df_mate$Gene)
  mate_lu_side <- setNames(df_mate$Side, df_mate$Gene)
  sex_tbl$LFC_mating  <- mate_lu_lfc[sex_tbl$Gene]
  sex_tbl$padj_mating <- mate_lu_padj[sex_tbl$Gene]
  sex_tbl$Mating_DE   <- sex_tbl$Gene %in% genes_mate

  # Classification
  sex_tbl$Category <- ifelse(
    !sex_tbl$Mating_DE, "Sex-dimorphic only",
    ifelse(sign(sex_tbl$LFC_sex) == sign(sex_tbl$LFC_mating),
           "Sex-and-mating-regulated (reinforced)",
           "Sex-and-mating-regulated (reversed)")
  )

  # Expression means
  sex_tbl <- merge(sex_tbl, norm_means, by.x = "Gene", by.y = "JoinKey", all.x = TRUE)

  # Chemosensory family
  sex_tbl$Chemosensory_family <- chemo_fam_lookup[sex_tbl$Gene]

  # Rank by |LFC_sex| descending
  sex_tbl$abs_LFC_sex <- abs(sex_tbl$LFC_sex)
  sex_tbl <- sex_tbl %>% arrange(desc(abs_LFC_sex))
  sex_tbl$Rank <- seq_len(nrow(sex_tbl))

  # Remove helper column
  sex_tbl$abs_LFC_sex <- NULL

  tissue_clean <- gsub(" ", "_", tissue)

  # ── Save per-appendage supplementary CSV (S3/S4/S5) ─────────────────────
  s_num_map <- c(Antenna = "S3", `Maxillary palp` = "S4", Tarsi = "S5")
  s_num     <- s_num_map[tissue]
  supp_fname <- paste0(s_num, "_", tissue_clean, "_sex_dimorphic_genes.csv")
  write.csv(sex_tbl, file.path(FIG_DIR, supp_fname), row.names = FALSE)
  write_rpt("  Saved: ", supp_fname, " (", nrow(sex_tbl), " sex-dimorphic genes, ranked by |LFC_sex|)")

  # ── Accumulate S2 rows for this appendage ────────────────────────────────
  s2_sex_part <- sex_tbl[, c("Gene", "Name", "Category", "GO_Domain", "GO_Name",
                               "LFC_sex", "padj_sex", "Sex_bias",
                               "LFC_mating", "padj_mating", "Chemosensory_family")]
  if (length(only_mate) > 0) {
    mate_only <- data.frame(Gene = only_mate, stringsAsFactors = FALSE)
    if (!is.null(norm_name_map)) mate_only$Name <- norm_name_map[mate_only$Gene]
    mate_only$Name[is.na(mate_only$Name) | mate_only$Name == ""] <-
      mate_only$Gene[is.na(mate_only$Name) | mate_only$Name == ""]
    mate_only$GO_Domain   <- go_map$GO_Domain[match(mate_only$Gene, go_map$Gene)]
    mate_only$GO_Domain[is.na(mate_only$GO_Domain)] <- "Unknown"
    mate_only$GO_Name     <- go_map$GO_Name[match(mate_only$Gene, go_map$Gene)]
    mate_only$GO_Name[is.na(mate_only$GO_Name)] <- "Unknown"
    mate_only$LFC_sex     <- df_sex$log2FoldChange[match(mate_only$Gene, df_sex$Gene)]
    mate_only$padj_sex    <- df_sex$padj[match(mate_only$Gene, df_sex$Gene)]
    mate_only$Sex_bias    <- NA_character_
    mate_only$LFC_mating  <- sig_mate$log2FoldChange[match(mate_only$Gene, sig_mate$Gene)]
    mate_only$padj_mating <- sig_mate$padj[match(mate_only$Gene, sig_mate$Gene)]
    mate_only$Category    <- "Mating-regulated only"
    mate_only$Chemosensory_family <- chemo_fam_lookup[mate_only$Gene]
    s2_tissue <- rbind(s2_sex_part, mate_only[, names(s2_sex_part)])
  } else {
    s2_tissue <- s2_sex_part
  }
  s2_tissue$Appendage <- tissue
  s2_list[[tissue]] <- s2_tissue[, c("Appendage", "Gene", "Name", "Category",
                                      "LFC_sex", "padj_sex", "Sex_bias",
                                      "LFC_mating", "padj_mating",
                                      "GO_Domain", "GO_Name", "Chemosensory_family")]

  # ── LFC correlation scatter for sex-and-mating-regulated genes ─────────────────────
  dual <- sex_tbl[sex_tbl$Mating_DE, ]
  n_conc <- sum(dual$Category == "Sex-and-mating-regulated (reinforced)")
  n_disc <- sum(dual$Category == "Sex-and-mating-regulated (reversed)")

  if (nrow(dual) > 0) {
    r_val <- cor(dual$LFC_sex, dual$LFC_mating, use = "complete.obs")

    # GO domain with chemosensory priority
    dual$PlotDomain <- dual$GO_Domain
    dual$PlotDomain[!is.na(dual$Chemosensory_family)] <- "Chemosensory"
    dual$PlotDomain <- factor(dual$PlotDomain,
      levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))

    # Quadrant annotations placed at corners of the fixed axis range
    q_pos <- SCATTER_AXIS_LIM * 0.95
    quad_df <- data.frame(
      x     = c( q_pos, -q_pos, -q_pos,  q_pos),
      y     = c( q_pos,  q_pos, -q_pos, -q_pos),
      hjust = c(1, 0, 0, 1),
      vjust = c(1, 1, 0, 0),
      txt = c(
        "Female-biased;\nupregulated in\nmated females",
        "Male-biased;\nupregulated in\nmated females",
        "Male-biased;\ndownregulated in\nmated females",
        "Female-biased;\ndownregulated in\nmated females"
      ), stringsAsFactors = FALSE
    )

    p2 <- ggplot(dual, aes(x = LFC_sex, y = LFC_mating, color = PlotDomain)) +
      annotate("rect",
               xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
               fill = "#FFE6F0", alpha = 0.3) +
      annotate("rect",
               xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
               fill = "#E6EEFF", alpha = 0.3) +
      annotate("rect",
               xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
               fill = "#FFF3E6", alpha = 0.3) +
      annotate("rect",
               xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
               fill = "#FFF3E6", alpha = 0.3) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      geom_vline(xintercept = 0, linewidth = 0.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "red", linewidth = 0.3) +
      geom_point(size = 3.5, alpha = 0.7) +
      geom_text(data = quad_df, aes(x = x, y = y, label = txt,
                hjust = hjust, vjust = vjust),
                inherit.aes = FALSE, size = 2.5, color = "grey35",
                lineheight = 0.85, fontface = "italic") +
      scale_color_manual(values = GO_COLS, labels = DOMAIN_FULL,
                         breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown"),
                         name = NULL) +
      annotate("text", x = SCATTER_AXIS_LIM * 0.95, y = SCATTER_AXIS_LIM * 0.55,
               hjust = 1, size = 3, color = "grey40",
               label = paste0("reinforced: ", n_conc,
                              "\nreversed: ", n_disc)) +
      scale_x_continuous(limits = c(-SCATTER_AXIS_LIM, SCATTER_AXIS_LIM)) +
      scale_y_continuous(limits = c(-SCATTER_AXIS_LIM, SCATTER_AXIS_LIM)) +
      coord_fixed() +
      labs(title = paste0(tissue, " — sex-and-mating-regulated DEGs (n=", nrow(dual),
                          ", r=", sprintf("%.2f", r_val), ")"),
           x = expression(Log[2]*FC~sex~effect~(VF~vs~Vm):
                           ~~positive==female~biased~","~~negative==male~biased),
           y = expression(Log[2]*FC~mating~effect~(MF~vs~VF):
                           ~~positive==up~after~mating~","~~negative==down)) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
            legend.position = "bottom",
            axis.title.x = element_text(size = 7.5),
            axis.title.y = element_text(size = 7.5),
            panel.grid = element_blank())
    print(p2)
  }

  # ── Report ─────────────────────────────────────────────────────────────────
  n_fb <- sum(sex_tbl$Sex_bias == "Female-biased")
  n_mb <- sum(sex_tbl$Sex_bias == "Male-biased")

  write_rpt("  ── ", tissue, " ──")
  write_rpt("  Total sex-dimorphic: ", nrow(sex_tbl),
            " (female-biased: ", n_fb, ", male-biased: ", n_mb, ")")
  write_rpt("  Mating-regulated only: ", length(only_mate))
  write_rpt("  Sex-and-mating-regulated: ", nrow(dual),
            " (reinforced: ", n_conc, ", reversed: ", n_disc, ")")
  if (nrow(dual) > 0) {
    write_rpt("  Pearson r (sex LFC vs mating LFC): ", sprintf("%.3f", r_val))
  }
  write_rpt("")

  # Top 10 most sex-dimorphic genes
  write_rpt("  Top 10 most sex-dimorphic genes (ranked by |LFC_sex|):")
  top10 <- head(sex_tbl, 10)
  for (i in seq_len(nrow(top10))) {
    g <- top10[i, ]
    gene_label <- if (!is.na(g$Name) && g$Name != "") g$Name else g$Gene
    write_rpt("    ", g$Rank, ". ", gene_label,
              " (LFC_sex=", sprintf("%.2f", g$LFC_sex),
              ", ", g$Sex_bias,
              ", GO=", g$GO_Name,
              ", ", g$Category,
              ifelse(g$Mating_DE,
                     paste0(", LFC_mating=", sprintf("%.2f", g$LFC_mating)),
                     ""),
              ifelse(!is.na(g$Chemosensory_family),
                     paste0(", ", g$Chemosensory_family), ""),
              ")")
  }
  write_rpt("")

  # GO name distribution of sex-dimorphic genes
  go_counts <- sex_tbl %>%
    filter(GO_Name != "Unknown") %>%
    group_by(GO_Name) %>%
    summarise(N = n(), mean_abs_LFC = mean(abs(LFC_sex)),
              n_dual = sum(Mating_DE), .groups = "drop") %>%
    arrange(desc(N)) %>%
    head(15)

  if (nrow(go_counts) > 0) {
    write_rpt("  Top GO terms among sex-dimorphic genes:")
    for (i in seq_len(nrow(go_counts))) {
      gc <- go_counts[i, ]
      write_rpt("    ", gc$GO_Name, ": N=", gc$N,
                ", mean|LFC|=", sprintf("%.2f", gc$mean_abs_LFC),
                ", ", gc$n_dual, " also mating-responsive")
    }
    write_rpt("")
  }

  # Chemosensory genes among sex-dimorphic
  chemo_sex <- sex_tbl[!is.na(sex_tbl$Chemosensory_family), ]
  if (nrow(chemo_sex) > 0) {
    write_rpt("  Chemosensory genes among sex-dimorphic (", nrow(chemo_sex), "):")
    for (i in seq_len(nrow(chemo_sex))) {
      g <- chemo_sex[i, ]
      write_rpt("    ", g$Name,
                " (", g$Chemosensory_family, ", LFC_sex=", sprintf("%.2f", g$LFC_sex),
                ", ", g$Sex_bias, ", ", g$Category, ")")
    }
    write_rpt("")
  }
}
cat("Done: Sex-dimorphic overlap analysis\n")

# ── Write S2 ─────────────────────────────────────────────────────────────────
s2_all <- do.call(rbind, s2_list)
rownames(s2_all) <- NULL
# Clean up Category labels for readability
s2_all$Category <- gsub("^Sex-and-mating-regulated \\(reinforced\\)$",
  "Sex-and-mating-regulated: direction reinforced by mating", s2_all$Category)
s2_all$Category <- gsub("^Sex-and-mating-regulated \\(reversed\\)$",
  "Sex-and-mating-regulated: direction reversed by mating", s2_all$Category)
# Placeholder columns for BLAST (blast_drosophila.py will fill these)
s2_all$Dmel_hit                 <- ""
s2_all$Dmel_gene_symbol         <- ""
s2_all$Dmel_protein_description <- ""
s2_all$Dmel_evalue              <- ""
s2_all$Dmel_pct_aa_identity     <- ""
s2_all$Dmel_query_coverage_pct  <- ""
s2_all$Dmel_bitscore            <- ""
write.csv(s2_all, file.path(FIG_DIR, "S2_sex_mating_classification.csv"), row.names = FALSE)
write_rpt("  S2_sex_mating_classification.csv: ", nrow(s2_all), " rows (",
          nrow(s2_all[s2_all$Appendage == "Antenna", ]), " antenna / ",
          nrow(s2_all[s2_all$Appendage == "Maxillary palp", ]), " palp / ",
          nrow(s2_all[s2_all$Appendage == "Tarsi", ]), " tarsi)")


# =============================================================================
# 6. SUPPLEMENTARY CSVs — ALL DEGs COMBINED (S1)
# =============================================================================

write_rpt("================================================================")
write_rpt("6. SUPPLEMENTARY CSVs")
write_rpt("================================================================")

s1_rows <- list()
for (key in names(all_contrasts)) {
  it <- all_contrasts[[key]]
  sig <- it$df[it$df$is_sig, ]
  if (nrow(sig) == 0) next

  out <- data.frame(
    Appendage = it$tissue,
    Contrast  = paste0(it$cond1, "_vs_", it$cond2),
    Gene      = sig$Gene,
    stringsAsFactors = FALSE
  )
  if (!is.null(norm_name_map)) out$Name <- norm_name_map[out$Gene]
  out$Name[is.na(out$Name) | out$Name == ""] <- out$Gene[is.na(out$Name) | out$Name == ""]
  idx <- match(out$Gene, sig$Gene)
  out$log2FoldChange      <- sig$log2FoldChange[idx]
  out$padj                <- sig$padj[idx]
  out$Bias_toward         <- COND_FULL[sig$Side[idx]]
  out$Effect_strength     <- sig$Strength[idx]
  out$GO_Domain           <- sig$GO_Domain[idx]
  out$GO_Name             <- sig$GO_Name[idx]
  out$Chemosensory_family <- chemo_fam_lookup[out$Gene]
  out <- out %>% arrange(padj)
  s1_rows[[key]] <- out
}

s1_all <- do.call(rbind, s1_rows)
rownames(s1_all) <- NULL

# ── Merge BLAST results from cache ──────────────────────────────────────────
blast_cache_path <- file.path(FIG_DIR, "blast_cache.json")
if (file.exists(blast_cache_path) && requireNamespace("jsonlite", quietly = TRUE)) {
  blast_cache <- jsonlite::read_json(blast_cache_path)
  dmel_fields <- c("Dmel_hit", "Dmel_gene_symbol", "Dmel_protein_description",
                   "Dmel_evalue", "Dmel_pct_aa_identity", "Dmel_query_coverage_pct",
                   "Dmel_bitscore")
  for (fld in dmel_fields) {
    s1_all[[fld]] <- sapply(s1_all$Gene, function(g) {
      v <- blast_cache[[g]][[fld]]
      if (is.null(v)) "" else as.character(v)
    })
  }
  cat("BLAST data merged into S1 from cache (", length(blast_cache), "entries)\n")
} else {
  for (fld in c("Dmel_hit", "Dmel_gene_symbol", "Dmel_protein_description",
                "Dmel_evalue", "Dmel_pct_aa_identity", "Dmel_query_coverage_pct",
                "Dmel_bitscore")) {
    s1_all[[fld]] <- ""
  }
  cat("blast_cache.json not found — Dmel columns left empty\n")
}

write.csv(s1_all, file.path(FIG_DIR, "S1_all_DEGs.csv"), row.names = FALSE)
write_rpt("  S1_all_DEGs.csv: ", nrow(s1_all), " rows (6 contrasts × 3 appendages)")

write_rpt("")

# ── 00_csv_index_readme.csv ─────────────────────────────────────────────────
index_data <- data.frame(
  file = c(
    "S1_all_DEGs.csv",
    "S2_sex_mating_classification.csv",
    "S3_Antenna_sex_dimorphic_genes.csv",
    "S4_Maxillary_palp_sex_dimorphic_genes.csv",
    "S5_Tarsi_sex_dimorphic_genes.csv"
  ),
  contains = c(
    paste0("All significant DEGs (padj<", PADJ_THR, ", |LFC|>=", LFC_THR, ") from all 6 contrasts.",
           " Columns: Appendage, Contrast, Gene, Name, log2FoldChange, padj,",
           " Bias_toward, Effect_strength, GO_Domain, GO_Name, Chemosensory_family, Dmel_* (BLAST orthologs)."),
    paste0("Sex/mating classification for all 3 appendages combined.",
           " Category: Sex-dimorphic only / Sex-and-mating-regulated (reinforced or reversed) / Mating-regulated only.",
           " Includes LFC and padj from both VF_vs_Vm and MF_vs_VF contrasts, GO annotation,",
           " Chemosensory_family, and Dmel_* BLAST ortholog columns."),
    paste0("Antenna: all sex-dimorphic genes (VF vs Vm) ranked by |LFC_sex|.",
           " Includes GO annotation, Sex_bias (Female/Male-biased), mating response LFC & padj,",
           " per-group normalized expression means, Category, and Chemosensory_family."),
    paste0("Maxillary palp: all sex-dimorphic genes (VF vs Vm) ranked by |LFC_sex|.",
           " Same columns as S3."),
    paste0("Tarsi: all sex-dimorphic genes (VF vs Vm) ranked by |LFC_sex|.",
           " Same columns as S3.")
  ),
  stringsAsFactors = FALSE
)
write.csv(index_data, file.path(FIG_DIR, "00_csv_index_readme.csv"), row.names = FALSE)
write_rpt("  00_csv_index_readme.csv")

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

# =============================================================================
# SAVE ALL PLOTS AS PDFs
# =============================================================================
# Close the combined-preview PDF device opened at the top (non-interactive runs)
if (!interactive() && dev.cur() > 1) dev.off()
cat("\nSaving individual PDFs to:", PDF_DIR, "\n")

save_pdf <- function(p, filename, w = 8, h = 6) {
  ggsave(file.path(PDF_DIR, filename), plot = p, width = w, height = h,
         units = "in", device = "pdf")
}

# ── 1b. GO-colored volcanos ────────────────────────────────────────────────
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      p <- make_volcano_go(it, VOLCANO_Y_MAX)
      fname <- paste0("1b_volcano_GO_", gsub(" ", "_", tissue), "_",
                       it$cond1, "_vs_", it$cond2, ".pdf")
      save_pdf(p, fname, 8, 6)
    }
  }
}

# ── 2b. GO-colored bars ───────────────────────────────────────────────────
for (tissue in TISSUES_ORDER) {
  for (contrast in CONTRAST_ORDER) {
    key <- paste0(tissue, " ", contrast)
    if (!is.null(all_contrasts[[key]])) {
      it <- all_contrasts[[key]]
      p <- make_go_de_bar(it)
      if (!is.null(p)) {
        fname <- paste0("2b_bar_GO_", gsub(" ", "_", tissue), "_",
                         it$cond1, "_vs_", it$cond2, ".pdf")
        save_pdf(p, fname, 6, 5)
      }
    }
  }
}

# ── 3. Venn diagrams ──────────────────────────────────────────────────────
for (vc in venn_configs) {
  sets <- get_up_genes(vc$contrast, vc$side)
  sA <- if (!is.null(sets[["Antenna"]]))        sets[["Antenna"]]        else character(0)
  sB <- if (!is.null(sets[["Maxillary palp"]])) sets[["Maxillary palp"]] else character(0)
  sC <- if (!is.null(sets[["Tarsi"]]))          sets[["Tarsi"]]          else character(0)
  regions <- compute_venn_regions(sA, sB, sC)
  p <- draw_venn_3way(vc$title, vc$color, tissue_labels_venn, regions)
  fname <- paste0("3_venn_", gsub(" ", "_", gsub("[^A-Za-z0-9 ]", "", vc$title)), ".pdf")
  save_pdf(p, fname, 7, 6)
}

# ── 4. Violin plots ──────────────────────────────────────────────────────
for (tissue in TISSUES_ORDER) {
  bg_sub <- violin_bg[violin_bg$Tissue == tissue, ]
  fg_sub <- violin_fg[violin_fg$Tissue == tissue, ]
  if (nrow(fg_sub) == 0) next

  jitter_counts <- fg_sub %>%
    group_by(Contrast, Up_in) %>%
    summarise(N = n(), .groups = "drop")
  jitter_counts$x_pos <- ifelse(
    jitter_counts$Up_in %in% c("Virgin female", "Mated female"), 8, -8)
  jitter_counts$Label <- paste0("n=", jitter_counts$N)

  p <- ggplot() +
    geom_violin(data = bg_sub, aes(y = Contrast, x = LFC),
                fill = "grey90", color = "grey70", alpha = 0.5,
                scale = "width", width = 0.8) +
    geom_jitter(data = fg_sub, aes(y = Contrast, x = LFC, color = Up_in),
                height = 0.15, size = 1, alpha = 0.6) +
    geom_text(data = jitter_counts,
              aes(y = Contrast, x = x_pos, label = Label, color = Up_in),
              size = 3, fontface = "bold", show.legend = FALSE) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
               linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_color_manual(values = cond_fill, name = "Upregulated in") +
    scale_x_continuous(limits = c(-10, 10)) +
    labs(title = tissue,
         x = expression(Log[2](FoldChange)), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          legend.position = "bottom", panel.grid.minor = element_blank())

}

# ── 5. Sex-dimorphic overlap: Venns + GO-domain Venns + stacked bars + scatters
for (tissue in TISSUES_ORDER) {
  key_sex  <- paste0(tissue, " VF_vs_Vm")
  key_mate <- paste0(tissue, " MF_vs_VF")
  if (is.null(all_contrasts[[key_sex]]) || is.null(all_contrasts[[key_mate]])) next

  df_sex_p  <- all_contrasts[[key_sex]]$df
  df_mate_p <- all_contrasts[[key_mate]]$df
  sig_sex_p  <- df_sex_p[df_sex_p$is_sig, ]
  sig_mate_p <- df_mate_p[df_mate_p$is_sig, ]

  genes_sex_p  <- unique(sig_sex_p$Gene)
  genes_mate_p <- unique(sig_mate_p$Gene)
  both_p       <- intersect(genes_sex_p, genes_mate_p)
  only_sex_p   <- setdiff(genes_sex_p, genes_mate_p)
  only_mate_p  <- setdiff(genes_mate_p, genes_sex_p)

  tc <- gsub(" ", "_", tissue)

  # Overall Venn
  circles <- data.frame(x0 = c(-0.5, 0.5), y0 = c(0, 0), r = c(1, 1))
  lbl_df <- data.frame(
    x = c(-0.85, 0.85, 0, -0.5, 0.5),
    y = c(0, 0, 0.15, -1.0, -1.0),
    txt = c(length(only_sex_p), length(only_mate_p), length(both_p),
            "Sex-dimorphic\nonly", "Mating-regulated\nonly"),
    sz = c(5, 5, 5, 3.5, 3.5), stringsAsFactors = FALSE)
  venn_annot <- data.frame(
    x = c(-1.6, 0, 1.6), y = c(1.3, -0.55, 1.3),
    txt = c("DE between VF & Vm\nbut NOT after mating\n(constitutive sex\ndifference)",
            "DE in BOTH contrasts\n(mating modulates a\nsex-dimorphic gene)",
            "DE between MF & VF\nbut NOT between sexes\n(mating-specific\nchange)"),
    stringsAsFactors = FALSE)
  p_ov <- ggplot() +
    ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r),
      fill = adjustcolor("#555555", alpha.f = 0.15), color = NA) +
    geom_text(data = lbl_df, aes(x = x, y = y, label = txt), size = lbl_df$sz) +
    geom_text(data = venn_annot, aes(x = x, y = y, label = txt),
              size = 2.3, color = "grey30", lineheight = 0.9) +
    coord_fixed(xlim = c(-2.5, 2.5), ylim = c(-1.8, 2.0)) +
    labs(title = paste0(tissue, " — sex-dimorphic vs mating-responsive DEGs"),
         subtitle = paste0(
           "Left: DE in VF vs Vm (sex-dimorphic)    ",
           "Right: DE in MF vs VF (mating-responsive)")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 7, color = "grey40"))
  save_pdf(p_ov, paste0("5a_venn_overview_", tc, ".pdf"), 8, 6)

  # Per-domain Venns (including Unknown)
  assign_go_domain_pdf <- function(gene_vec) {
    dom <- go_map$GO_Domain[match(gene_vec, go_map$Gene)]
    dom[is.na(dom)] <- "Unknown"
    dom[!is.na(chemo_fam_lookup[gene_vec])] <- "Chemosensory"
    return(dom)
  }

  venn_all <- data.frame(
    Gene = c(only_sex_p, both_p, only_mate_p),
    Region = c(rep("Sex-dimorphic only", length(only_sex_p)),
               rep("Sex-and-mating-regulated", length(both_p)),
               rep("Mating-regulated only", length(only_mate_p))),
    stringsAsFactors = FALSE)
  venn_all$GO_Domain <- assign_go_domain_pdf(venn_all$Gene)

  for (dom in c("Chemosensory", "MF", "CC", "BP", "Unknown")) {
    dom_sex  <- venn_all$Gene[venn_all$Region != "Mating-regulated only" &
                               venn_all$GO_Domain == dom]
    dom_mate <- venn_all$Gene[venn_all$Region != "Sex-dimorphic only" &
                               venn_all$GO_Domain == dom]
    n_os <- length(setdiff(dom_sex, dom_mate))
    n_om <- length(setdiff(dom_mate, dom_sex))
    n_b  <- length(intersect(dom_sex, dom_mate))
    circles_d <- data.frame(x0 = c(-0.5, 0.5), y0 = c(0, 0), r = c(1, 1))
    lbl_d <- data.frame(
      x = c(-0.85, 0.85, 0, -0.5, 0.5),
      y = c(0, 0, 0.15, -0.9, -0.9),
      txt = c(as.character(n_os), as.character(n_om), as.character(n_b),
              "Sex-dimorphic", "Mating-regulated"),
      sz = c(4, 4, 5, 3.5, 3.5), stringsAsFactors = FALSE)
    dom_label <- if (dom == "Chemosensory") "Chemosensory" else DOMAIN_FULL[dom]
    p_d <- ggplot() +
      ggforce::geom_circle(data = circles_d, aes(x0 = x0, y0 = y0, r = r),
        fill = adjustcolor(GO_COLS[dom], alpha.f = 0.25), color = NA) +
      geom_text(data = lbl_d, aes(x = x, y = y, label = txt), size = lbl_d$sz) +
      coord_fixed(xlim = c(-2, 2), ylim = c(-1.5, 1.5)) +
      labs(title = paste0(tissue, " — ", dom_label)) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10,
                                       color = GO_COLS[dom]))
    save_pdf(p_d, paste0("5b_venn_", dom, "_", tc, ".pdf"), 6, 5)
  }

  # ── BLAST characterization of Unknown-domain genes (PDF section) ──────────
  if (!is.null(dmel_lookup)) {
    unk_pdf <- venn_all[venn_all$GO_Domain == "Unknown", ]
    unk_pdf <- merge(unk_pdf, dmel_lookup[, c("Gene","Dmel_hit","Dmel_clean",
                                               "Dmel_protein_description")],
                    by = "Gene", all.x = TRUE)
    n_unk  <- nrow(unk_pdf)
    n_hit  <- sum(unk_pdf$Dmel_hit == "TRUE", na.rm = TRUE)
    n_miss <- n_unk - n_hit

    unk_lbl_dir <- file.path(FIG_DIR, "Unknown_BLAST")
    dir.create(unk_lbl_dir, showWarnings = FALSE)

    txt_u <- c(
      strrep("=", 72),
      paste0("  Unknown-domain DEGs — BLAST characterization: ", tissue),
      strrep("=", 72),
      paste0("  Total Unknown DEGs: ", n_unk,
             "  |  Dmel BLAST hit: ", n_hit,
             "  |  No hit: ", n_miss),
      ""
    )
    for (reg in c("Sex-dimorphic only", "Sex-and-mating-regulated", "Mating-regulated only")) {
      sub_reg <- unk_pdf[unk_pdf$Region == reg, ]
      if (nrow(sub_reg) == 0) next
      txt_u <- c(txt_u, paste0("  ", reg, " (n=", nrow(sub_reg), "):"))
      sub_hit <- sub_reg[sub_reg$Dmel_hit == "TRUE" & !is.na(sub_reg$Dmel_hit), ]
      sub_noh <- sub_reg[sub_reg$Dmel_hit != "TRUE" | is.na(sub_reg$Dmel_hit), ]
      if (nrow(sub_hit) > 0) {
        for (i in seq_len(min(nrow(sub_hit), 15))) {
          rr <- sub_hit[i, ]
          nm <- if (!is.null(norm_name_map) && !is.na(norm_name_map[rr$Gene])) {
            norm_name_map[rr$Gene]
          } else rr$Gene
          txt_u <- c(txt_u, sprintf("    • %-18s  Dmel: %-30s  | %s",
                                    nm, rr$Dmel_clean,
                                    substr(rr$Dmel_protein_description, 1, 50)))
        }
        if (nrow(sub_hit) > 15)
          txt_u <- c(txt_u, paste0("    ... and ", nrow(sub_hit)-15, " more"))
      }
      if (nrow(sub_noh) > 0)
        txt_u <- c(txt_u, paste0("    No Dmel hit: ", nrow(sub_noh), " genes"))
      txt_u <- c(txt_u, "")
    }
    writeLines(txt_u, file.path(unk_lbl_dir,
                                paste0("Unknown_BLAST_", tc, ".txt")))
    cat("Unknown BLAST report written:", tc, "(", n_unk, "genes )\n")
  }

  # GO-domain stacked bar
  if (!is.null(norm_name_map)) {
    venn_all$Label <- norm_name_map[venn_all$Gene]
    venn_all$Label[is.na(venn_all$Label) | venn_all$Label == ""] <-
      venn_all$Gene[is.na(venn_all$Label) | venn_all$Label == ""]
  } else { venn_all$Label <- venn_all$Gene }
  venn_all$GO_Domain <- factor(venn_all$GO_Domain,
    levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))
  venn_all$Region <- factor(venn_all$Region,
    levels = c("Sex-dimorphic only", "Sex-and-mating-regulated", "Mating-regulated only"))
  vgc <- venn_all %>%
    group_by(Region, GO_Domain, .drop = FALSE) %>%
    summarise(N = n(), .groups = "drop")
  vgc$NLabel <- ifelse(vgc$N > 0, as.character(vgc$N), "")
  p_sb <- ggplot(vgc, aes(x = Region, y = N, fill = GO_Domain)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    geom_text(aes(label = NLabel), position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = GO_COLS, labels = DOMAIN_FULL,
                      breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown")) +
    labs(title = paste0(tissue, " — GO domain breakdown"),
         x = NULL, y = "Number of DEGs", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "bottom", axis.text.x = element_text(size = 9))
  save_pdf(p_sb, paste0("5c_GO_bar_", tc, ".pdf"), 7, 5)

  # LFC scatter
  sex_tbl_p <- data.frame(Gene = genes_sex_p, stringsAsFactors = FALSE)
  sex_tbl_p$LFC_sex  <- sig_sex_p$log2FoldChange[match(sex_tbl_p$Gene, sig_sex_p$Gene)]
  sex_tbl_p$Sex_bias <- ifelse(sex_tbl_p$LFC_sex > 0, "Female-biased", "Male-biased")
  sex_tbl_p$LFC_mating <- df_mate_p$log2FoldChange[match(sex_tbl_p$Gene, df_mate_p$Gene)]
  sex_tbl_p$padj_mating <- df_mate_p$padj[match(sex_tbl_p$Gene, df_mate_p$Gene)]
  sex_tbl_p$Mating_DE <- sex_tbl_p$Gene %in% genes_mate_p
  sex_tbl_p$Category <- ifelse(!sex_tbl_p$Mating_DE, "Sex-dimorphic only",
    ifelse(sign(sex_tbl_p$LFC_sex) == sign(sex_tbl_p$LFC_mating),
           "Sex-and-mating-regulated (reinforced)", "Sex-and-mating-regulated (reversed)"))
  sex_tbl_p$Chemosensory_family <- chemo_fam_lookup[sex_tbl_p$Gene]
  if (!is.null(norm_name_map)) sex_tbl_p$Name <- norm_name_map[sex_tbl_p$Gene]
  sex_tbl_p$GO_Domain <- go_map$GO_Domain[match(sex_tbl_p$Gene, go_map$Gene)]
  sex_tbl_p$GO_Domain[is.na(sex_tbl_p$GO_Domain)] <- "Unknown"
  sex_tbl_p$GO_Name <- go_map$GO_Name[match(sex_tbl_p$Gene, go_map$Gene)]
  sex_tbl_p$GO_Name[is.na(sex_tbl_p$GO_Name) | sex_tbl_p$GO_Name == ""] <- "Unknown"

  dual_p <- sex_tbl_p[sex_tbl_p$Mating_DE, ]
  n_conc_p <- sum(dual_p$Category == "Sex-and-mating-regulated (reinforced)")
  n_disc_p <- sum(dual_p$Category == "Sex-and-mating-regulated (reversed)")

  if (nrow(dual_p) > 0) {
    r_val_p <- cor(dual_p$LFC_sex, dual_p$LFC_mating, use = "complete.obs")
    dual_p$PlotDomain <- dual_p$GO_Domain
    dual_p$PlotDomain[!is.na(dual_p$Chemosensory_family)] <- "Chemosensory"
    dual_p$PlotDomain <- factor(dual_p$PlotDomain,
      levels = c("Chemosensory", "MF", "CC", "BP", "Unknown"))
    q_pos <- SCATTER_AXIS_LIM * 0.95
    quad_df_p <- data.frame(
      x = c(q_pos, -q_pos, -q_pos, q_pos),
      y = c(q_pos, q_pos, -q_pos, -q_pos),
      hjust = c(1, 0, 0, 1), vjust = c(1, 1, 0, 0),
      txt = c("Female-biased;\nupregulated in\nmated females",
              "Male-biased;\nupregulated in\nmated females",
              "Male-biased;\ndownregulated in\nmated females",
              "Female-biased;\ndownregulated in\nmated females"),
      stringsAsFactors = FALSE)
    p_sc <- ggplot(dual_p, aes(x = LFC_sex, y = LFC_mating, color = PlotDomain)) +
      annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
               fill = "#FFE6F0", alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
               fill = "#E6EEFF", alpha = 0.3) +
      annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
               fill = "#FFF3E6", alpha = 0.3) +
      annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
               fill = "#FFF3E6", alpha = 0.3) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      geom_vline(xintercept = 0, linewidth = 0.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "red", linewidth = 0.3) +
      geom_point(size = 3.5, alpha = 0.7) +
      geom_text(data = quad_df_p, aes(x = x, y = y, label = txt,
                hjust = hjust, vjust = vjust),
                inherit.aes = FALSE, size = 2.5, color = "grey35",
                lineheight = 0.85, fontface = "italic") +
      scale_color_manual(values = GO_COLS, labels = DOMAIN_FULL,
                         breaks = c("Chemosensory", "MF", "CC", "BP", "Unknown"),
                         name = NULL) +
      annotate("text", x = q_pos, y = SCATTER_AXIS_LIM * 0.55,
               hjust = 1, size = 3, color = "grey40",
               label = paste0("reinforced: ", n_conc_p,
                              "\nreversed: ", n_disc_p)) +
      scale_x_continuous(limits = c(-SCATTER_AXIS_LIM, SCATTER_AXIS_LIM)) +
      scale_y_continuous(limits = c(-SCATTER_AXIS_LIM, SCATTER_AXIS_LIM)) +
      coord_fixed() +
      labs(title = paste0(tissue, " — sex-and-mating-regulated DEGs (n=", nrow(dual_p),
                          ", r=", sprintf("%.2f", r_val_p), ")"),
           x = expression(Log[2]*FC~sex~(VF~vs~Vm):~~positive==female~biased~","~~negative==male~biased),
           y = expression(Log[2]*FC~mating~(MF~vs~VF):~~positive==up~in~mated~","~~negative==down~in~mated)) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
            legend.position = "bottom",
            axis.title.x = element_text(size = 7.5),
            axis.title.y = element_text(size = 7.5),
            panel.grid = element_blank())
    save_pdf(p_sc, paste0("5d_scatter_", tc, ".pdf"), 8, 8)

    # ── Assign quadrant ────────────────────────────────────────────────────
    dual_p$Quadrant <- ifelse(
      dual_p$LFC_sex >  0 & dual_p$LFC_mating >  0, "Female-biased_mating-reinforced",
      ifelse(dual_p$LFC_sex <  0 & dual_p$LFC_mating >  0, "Male-biased_mating-reversed",
      ifelse(dual_p$LFC_sex <  0 & dual_p$LFC_mating <  0, "Male-biased_mating-reinforced",
                                                             "Female-biased_mating-reversed")))

    quad_meta <- list(
      `Female-biased_mating-reinforced` = "Female-biased; upregulated in mated females   (top-right)",
      `Male-biased_mating-reversed`     = "Male-biased; upregulated in mated females     (top-left)",
      `Male-biased_mating-reinforced`   = "Male-biased; downregulated in mated females   (bottom-left)",
      `Female-biased_mating-reversed`   = "Female-biased; downregulated in mated females (bottom-right)"
    )
    dom_order_sc <- c("Chemosensory", "MF", "CC", "BP", "Unknown")

    scatter_lbl_dir <- file.path(FIG_DIR, "Scatter_labels")
    dir.create(scatter_lbl_dir, showWarnings = FALSE)

    # ── Quadrant gene-list txt ─────────────────────────────────────────────
    txt_lines <- c(
      strrep("=", 75),
      paste0("  Scatter gene list: ", tissue),
      paste0("  Sex-and-mating-regulated genes: ", nrow(dual_p),
             "  (direction reinforced: ", n_conc_p, " | direction reversed: ", n_disc_p, ")"),
      paste0("  Pearson r = ", sprintf("%.3f", r_val_p)),
      paste0("  Format per line:  Gene_ID (Name) — GO_name  [LFC_sex / LFC_mating]"),
      strrep("=", 75)
    )

    for (qk in names(quad_meta)) {
      qg <- dual_p[dual_p$Quadrant == qk, ]
      if (nrow(qg) == 0) next
      txt_lines <- c(txt_lines, "",
                     paste0("── ", quad_meta[[qk]], "  (n=", nrow(qg), ")"),
                     strrep("─", 75))
      for (dom in dom_order_sc) {
        dg <- qg[as.character(qg$PlotDomain) == dom, ]
        if (nrow(dg) == 0) next
        txt_lines <- c(txt_lines, paste0("  [", dom, "]"))
        dg <- dg[order(dg$GO_Name, dg$Gene), ]
        for (i in seq_len(nrow(dg))) {
          rr  <- dg[i, ]
          nm  <- if (!is.null(norm_name_map)) norm_name_map[[rr$Gene]] else NA
          nm  <- if (!is.na(nm) && nchar(nm) > 0) paste0(" (", nm, ")") else ""
          go_lbl <- if (!is.na(rr$Chemosensory_family) && nchar(rr$Chemosensory_family) > 0) {
            paste0(rr$Chemosensory_family, " receptor")
          } else { rr$GO_Name }
          txt_lines <- c(txt_lines,
            sprintf("    %s%s — %s  [%+.2f / %+.2f]",
                    rr$Gene, nm, go_lbl, rr$LFC_sex, rr$LFC_mating))
        }
      }
    }
    txt_lines <- c(txt_lines, "", strrep("=", 75))
    writeLines(txt_lines,
               file.path(scatter_lbl_dir, paste0("5d_scatter_genes_", tc, ".txt")))
    cat("Scatter gene list written:", tc, "(", nrow(dual_p), "genes)\n")

    # ── BLAST CSV ──────────────────────────────────────────────────────────
    blast_sc <- dual_p[, c("Gene", "Quadrant", "PlotDomain", "GO_Name",
                           "Chemosensory_family", "LFC_sex", "LFC_mating", "Category")]
    names(blast_sc)[names(blast_sc) == "PlotDomain"] <- "Domain"
    nm_vec <- if (!is.null(norm_name_map)) norm_name_map[blast_sc$Gene] else rep("", nrow(blast_sc))
    nm_vec[is.na(nm_vec)] <- ""
    blast_sc$Name <- nm_vec
    if (!is.null(dmel_lookup)) {
      blast_sc <- merge(blast_sc, dmel_lookup, by = "Gene", all.x = TRUE)
    }
    blast_sc <- blast_sc[order(blast_sc$Quadrant, as.character(blast_sc$Domain)), ]
    write.csv(blast_sc,
              file.path(scatter_lbl_dir, paste0("5d_scatter_BLAST_", tc, ".csv")),
              row.names = FALSE)
    cat("Scatter BLAST CSV written:", tc, "\n")
  }
}

cat("All PDFs saved to:", PDF_DIR, "\n")


# ── Session info (run to capture package versions for Methods reporting) ─────
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
