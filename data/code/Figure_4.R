################################################################################
# Figure 4 — sex- and mating-dependent transcriptional remodelling
#
# WHAT THIS SCRIPT COVERS — main figure 4 and every supplementary figure that
# belongs to it. Nothing else is needed to reproduce them.
#
#   4A   sex volcanoes, virgin female vs virgin male, one per appendage
#   4B   mating volcanoes on the STRICT definition, one per appendage
#   4C   cross-appendage Venn diagrams, four transcriptional programmes
#   4D   within-appendage Venn, sex-biased vs mating-responsive
#   4E   sex bias in virgins versus sex bias in mated females
#   S20  GO over-representation, sex-biased genes
#   S21  GO over-representation, mating-responsive genes
#   S22  GO over-representation, sex-biased without a mating response
#   S23  GO over-representation, mating-responsive without a sex bias
#   S24  GO over-representation, genes regulated by both
#
# THIS IS THE FIGURE THE REVISION CHANGED MOST
#   "Mating-responsive" used to mean significant in mated female vs virgin
#   female alone. Because that contrast uses virgin females as its only
#   reference, a gene qualified whenever virgin females were the odd group out,
#   whether or not mated females differed from anything. Phi was right (C26,
#   C34). A gene now has to be significant in the SAME DIRECTION against BOTH
#   virgin groups, which takes the counts from 196/386/230 down to 146/318/33
#   and removes the "shared post-mating silencing programme" entirely.
#
#   Panel 4E was replaced rather than patched (C27-C36). Its old axes were sex
#   bias (VF vs Vm) against mating response (MF vs VF); both used virgin females
#   as reference with opposite sign, which forces a negative correlation by
#   construction. Under a permutation null the published r = -0.903 sat at
#   p = 0.29 — indistinguishable from the artefact. Both axes are now genuine
#   sex contrasts sharing virgin males as reference, so a sign flip really is a
#   reversal of sex bias.
#
# SIGN CONVENTION — VERIFIED, DO NOT "TIDY"
#   Additional files 4-9 and 29-31 all store positive = higher in the SECOND
#   named group. Every loader below negates on read, so downstream:
#       sex contrasts     positive = female-biased
#       mating contrasts  positive = mating-induced (higher in mated females)
#   Verified by correlating the stored log2FoldChange against log2(mean ratios)
#   recomputed from Additional file 10: r = -0.986, -0.999, -1.000.
#   (Additional files 1-3, used by Figure_3.R, go the other way and are NOT
#   negated. Getting that one backwards reports the antenna as
#   chemosensory-depleted.)
#
# INPUTS — all from 04_Additional_Files/, nothing else
#   AF4/6/8    sex contrasts       (virgin female vs virgin male)
#   AF5/7/9    mating contrasts    (mated female vs virgin female)
#   AF29/30/31 mated female vs virgin male — the second contrast that makes the
#              strict definition possible
#   AF11       GO annotation, exactly one term per transcript
#
# OUTPUTS
#   02_Main_Figures/panels_new/        Figure_4A_*.pdf, 4B, 4C, 4D, 4E + 3 CSVs
#   03_Supplementary_Figures/panels_new/  Figure_S20.pdf ... Figure_S24.pdf
#
# PACKAGES  ggplot2, dplyr, ggrepel, ggforce
#
# HOW TO RUN
#   setwd(".../Submission_v5_2026-09/05_R_Scripts")
#   source("Figure_4.R")
#
# CHECK THESE BEFORE ASSEMBLING — printed at the end of the run
#   4A  535 / 555 / 263 DEGs; male-biased 403 of 535 and 470 of 555;
#       tarsi female-biased 168 of 263
#   4B  146 / 318 / 33   (antenna 6 up + 140 down, palp 38 + 280, tarsi 27 + 6)
#   4C  all three appendages = 0 (was 16); antenna and palp share 95, all
#       suppressed; antenna only 51, palp only 219, tarsi only 29
#   4D  sex and mating overlap 41 / 168 / 3   (was 61 / 181 / 49)
#   4E  retention 398/535, 495/555, 15/263; reversals 0 / 8 / 3
#   If a run disagrees, stop and say so rather than adjusting the figure.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggrepel); library(ggforce)
})

# ── Edit this one line if the package folder ever moves ──────────────────────
PKG <- "C:/Users/dorpe/OneDrive/Desktop/PhD_Projects/BSF/Updated_GDrive_version/Submission_v5_2026-09"

AF_DIR   <- file.path(PKG, "04_Additional_Files")
MAIN_OUT <- file.path(PKG, "02_Main_Figures", "panels_new")
SUPP_OUT <- file.path(PKG, "03_Supplementary_Figures", "panels_new")
for (d in c(MAIN_OUT, SUPP_OUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

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

PADJ_THR <- 0.001
LFC_THR  <- 1.0     # >= throughout, matching the code and the revised captions
Q_THR    <- 0.05
APPENDAGES <- c("Antenna", "Maxillary palp", "Tarsi")

GO_COLS <- c(Chemosensory = "#A60000", MF = "#3FC498", CC = "#4D50DB",
             BP = "#A860E3", Unknown = "#666666")
PAL_GO  <- c(GO_COLS, `Not Significant` = "#D9D9D9")
DOMAIN_FULL <- c(MF = "Molecular function", CC = "Cellular component",
                 BP = "Biological process")
COL_MF <- c(light = "#C77CFF", dark = "#7B1FA2")
COL_VF <- c(light = "#FF66B3", dark = "#99004D")
COL_VM <- c(light = "#99B3FF", dark = "#003399")

# ── Which outputs to write ───────────────────────────────────────────────────
# The script computes every panel of Figure 4 and S20-S24 in one pass, which is
# right for a clean rebuild but wrong while iterating: re-running it to change
# one panel rewrites every file, including panels already placed in Illustrator.
# Pass one or more filename fragments on the command line to write only those:
#     Rscript Figure_4.R 4D            # just the scatter
#     Rscript Figure_4.R 4C 4D         # Venns and scatter
#     Rscript Figure_4.R               # everything (default)
# The analysis always runs in full; only the writing is filtered, so the printed
# counts remain a check on the whole script.
WRITE_ONLY <- commandArgs(trailingOnly = TRUE)
if (!length(WRITE_ONLY)) WRITE_ONLY <- "all"
want_write <- function(filename) {
  if ("all" %in% WRITE_ONLY) return(TRUE)
  any(vapply(WRITE_ONLY, function(k) grepl(k, filename, fixed = TRUE), logical(1)))
}

save_pdf <- function(p, dir, filename, w, h) {
  if (!want_write(filename)) { cat("skipped: ", filename, "\n", sep = ""); return(invisible(p)) }
  dev_fun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file.path(dir, filename), plot = p, width = w, height = h,
         units = "in", device = dev_fun)
  cat("written: ", file.path(dir, filename), "\n", sep = "")
  invisible(p)
}

# =============================================================================
# GO annotation — Additional file 11 carries exactly ONE term per transcript
# (24,828 rows, no gene with more than one), so a gene belongs to one term and
# one domain. That is also why a "GO domain composition" bar chart says nothing
# about biology, which is what retired the old panel 3B.
# =============================================================================
clean_domain <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[grepl("^molecular",  x)] <- "MF"
  out[grepl("^cellular",   x)] <- "CC"
  out[grepl("^biological", x)] <- "BP"
  out
}
go_raw <- read.csv(file.path(AF_DIR, "Additional_file_11_GO_annotations_all_genes.csv"),
                   stringsAsFactors = FALSE)
GO_MAP <- data.frame(
  JoinKey   = toupper(trimws(go_raw$Gene)),
  GO_Domain = clean_domain(go_raw$GO_Domain),
  GO_Name   = trimws(ifelse(is.na(go_raw$GO_Name), "", go_raw$GO_Name)),
  stringsAsFactors = FALSE)
GO_MAP <- GO_MAP[!duplicated(GO_MAP$JoinKey), ]
GO_ANN <- GO_MAP[!is.na(GO_MAP$GO_Domain) &
                 !GO_MAP$GO_Name %in% c("", "Unknown"), ]   # real annotations only

# =============================================================================
# Chemosensory identity
# The Name column of Additional files 4-9 is filled for exactly the 567
# annotated chemosensory transcript models and empty otherwise, so membership is
# "Name is non-empty" rather than a regex over symbols. Family is the leading
# alphabetic prefix, ORco folded into Or and every TRP variant into Trp.
# Totals reproduce Table 1: Or 192, Ir 101, Trp 111, Obp 75, Gr 39, Ppk 39, Csp 10.
# =============================================================================
chemo_family <- function(name) {
  n <- toupper(trimws(as.character(name)))
  out <- rep(NA_character_, length(n))
  out[grepl("^TRP",  n)] <- "Trp"      # before OR/IR so TRP* never leaks
  out[grepl("^ORCO", n)] <- "Or"
  out[is.na(out) & grepl("^OR",  n)] <- "Or"
  out[is.na(out) & grepl("^OBP", n)] <- "Obp"
  out[is.na(out) & grepl("^IR",  n)] <- "Ir"
  out[is.na(out) & grepl("^GR",  n)] <- "Gr"
  out[is.na(out) & grepl("^PPK", n)] <- "Ppk"
  out[is.na(out) & grepl("^CSP", n)] <- "Csp"
  out[n == ""] <- NA_character_
  out
}

# =============================================================================
# DE loader. `padj_raw` keeps the true NA so the ORA universe can be defined as
# "genes DESeq2 actually tested", which is the correct background.
# =============================================================================
read_de <- function(fname) {
  d <- read.csv(file.path(AF_DIR, fname), stringsAsFactors = FALSE)
  d$Gene    <- trimws(as.character(d$Gene))
  d$JoinKey <- toupper(d$Gene)
  d$padj_raw <- suppressWarnings(as.numeric(d$padj))
  d$padj <- ifelse(is.na(d$padj_raw), 1.0, d$padj_raw)
  d$lfc  <- -suppressWarnings(as.numeric(d$log2FoldChange))   # NEGATED, see header
  d$lfc[is.na(d$lfc)] <- 0.0
  if (!"Name" %in% names(d)) d$Name <- ""
  d$Name <- trimws(ifelse(is.na(d$Name), "", d$Name))
  d <- merge(d, GO_MAP, by = "JoinKey", all.x = TRUE)
  d$GO_Domain[is.na(d$GO_Domain)] <- "Unknown"
  d$Family    <- chemo_family(d$Name)
  d$is_chemo  <- !is.na(d$Family)
  d$is_sig    <- d$padj < PADJ_THR & abs(d$lfc) >= LFC_THR
  d$neg_log10 <- -log10(pmax(d$padj, .Machine$double.xmin))
  d$ColorKey  <- ifelse(!d$is_sig, "Not Significant",
                 ifelse(d$is_chemo, "Chemosensory", d$GO_Domain))
  d
}

SEX_FILES <- c(
  Antenna          = "Additional_file_4_Antenna_sex_VirginFemale_vs_VirginMale.csv",
  `Maxillary palp` = "Additional_file_6_MaxillaryPalp_sex_VirginFemale_vs_VirginMale.csv",
  Tarsi            = "Additional_file_8_Tarsi_sex_VirginFemale_vs_VirginMale.csv")
MATING_FILES <- c(
  Antenna          = "Additional_file_5_Antenna_mating_MatedFemale_vs_VirginFemale.csv",
  `Maxillary palp` = "Additional_file_7_MaxillaryPalp_mating_MatedFemale_vs_VirginFemale.csv",
  Tarsi            = "Additional_file_9_Tarsi_mating_MatedFemale_vs_VirginFemale.csv")
MATEDSEX_FILES <- c(
  Antenna          = "Additional_file_29_Antenna_matedsex_MatedFemale_vs_VirginMale.csv",
  `Maxillary palp` = "Additional_file_30_MaxillaryPalp_matedsex_MatedFemale_vs_VirginMale.csv",
  Tarsi            = "Additional_file_31_Tarsi_matedsex_MatedFemale_vs_VirginMale.csv")

cat("Reading Additional files 4-9 and 29-31 ...\n")
DAT <- lapply(APPENDAGES, function(a) list(
  sex      = read_de(SEX_FILES[[a]]),
  mating   = read_de(MATING_FILES[[a]]),
  matedsex = read_de(MATEDSEX_FILES[[a]])))
names(DAT) <- APPENDAGES

# ── The strict mating-responsive set ─────────────────────────────────────────
# Significant in the same direction against BOTH virgin groups, so mated females
# must differ from virgin females AND from virgin males.
strict_mating <- function(mating, matedsex) {
  m <- mating[mating$is_sig,   c("Gene", "lfc")]
  n <- matedsex[matedsex$is_sig, c("Gene", "lfc")]
  j <- merge(m, n, by = "Gene", suffixes = c("_m", "_n"))
  j <- j[sign(j$lfc_m) == sign(j$lfc_n), ]
  list(all = j$Gene, up = j$Gene[j$lfc_m > 0], down = j$Gene[j$lfc_m < 0])
}
SETS <- lapply(APPENDAGES, function(a) {
  d  <- DAT[[a]]
  sm <- strict_mating(d$mating, d$matedsex)
  list(sex_any       = unique(d$sex$Gene[d$sex$is_sig]),
       female_biased = unique(d$sex$Gene[d$sex$is_sig & d$sex$lfc > 0]),
       male_biased   = unique(d$sex$Gene[d$sex$is_sig & d$sex$lfc < 0]),
       mating_any    = sm$all, mating_up = sm$up, mating_down = sm$down)
})
names(SETS) <- APPENDAGES

# =============================================================================
# PANELS 4A and 4B — volcanoes
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
# Panel 4A/4B plot box, at FINAL size in the Illustrator layout.
# 2026-09-27: widened from 40.4481 x 24.2208 after the GO stacked bars were
# deleted from panels 4A/4B, which freed 28.8 mm per row (7.6 mm gap + 21.2 mm
# of bars). 60 x 25 uses most of that and leaves ~8 mm of air between the A and
# B blocks.
PANEL_W_MM <- 60
PANEL_H_MM <- 25
# Everything is drawn SCALE times larger - box, type, points, rules - so the PDF
# rasterises cleanly in previews and exports. Scale the placed art to 1/3 in
# Illustrator and the box lands on exactly 40.4481 x 24.2208 mm with the type at
# its intended size. Set SCALE <- 1 to emit final-size art directly.
SCALE <- 3
# x axis for Figure 4. No significant gene in any 4A/4B panel exceeds |log2FC| 12
# (checked: the only points beyond it are 8 + 6 + 2 non-significant ones, which
# are capped onto the edge). Figure 3 keeps its own wider +/-20 range.
X_MAX    <- 12                       # axis limit
X_BREAKS <- seq(-10, 10, 5)          # ticks at -10, -5, 0, 5, 10
Y_MAX    <- 300

cap_y <- function(v) pmin(v, Y_MAX)
cap_x <- function(v) pmax(pmin(v, X_MAX), -X_MAX)

volcano_x_scale <- function()
  scale_x_continuous(limits = c(-X_MAX, X_MAX), breaks = X_BREAKS,
                     expand = expansion(mult = 0.02))
volcano_y_scale <- function()
  scale_y_continuous(limits = c(0, Y_MAX), breaks = seq(0, Y_MAX, 100),
                     expand = expansion(mult = c(0.01, 0.04)))

# theme_bw, grid removed
volcano_theme <- function(base = 9 * SCALE)
  theme_bw(base_size = base) +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA,
                                        linewidth = 0.4 * SCALE),
        plot.title       = element_text(hjust = 0.5, face = "bold", size = 9 * SCALE),
        plot.subtitle    = element_text(hjust = 0.5, size = 7 * SCALE, colour = "grey30"),
        axis.text        = element_text(size = 7 * SCALE),
        axis.title       = element_text(size = 8 * SCALE),
        axis.ticks       = element_line(linewidth = 0.3 * SCALE),
        legend.position  = "bottom",
        legend.text      = element_text(size = 7 * SCALE),
        legend.key.size  = unit(3.2 * SCALE, "mm"),
        legend.margin    = margin(t = -2 * SCALE))

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

# Shared name shortening, also used by panel B (Figure_4B_lollipop.R), so the
# two panels of Figure 4 cannot disagree about what a gene is called.
source(file.path(PKG, "05_R_Scripts", "Figure_4_names.R"))

# raw (un-truncated) description for a gene, for shorten() to work on
raw_desc <- function(keys) {
  k <- toupper(trimws(keys))
  d <- unname(DESC_MAP[k]); d[is.na(d)] <- ""
  d[grepl("^uncharacterized", d, ignore.case = TRUE)] <- ""
  f <- unname(NCBI_MAP[k]);  f[is.na(f)] <- ""
  f[grepl("^uncharacterized", f, ignore.case = TRUE)] <- ""
  ifelse(nzchar(d), d, f)
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
    size = LAB_SIZE * SCALE, segment.size = 0.12 * SCALE, segment.colour = "grey60",
    min.segment.length = 0, box.padding = 0.28, point.padding = 0.12,
    max.overlaps = Inf, seed = 1)
}


# A CSV left open in Excel is locked on Windows, and write.csv then aborts the
# whole script - so panels written earlier survive while everything after the
# failure is silently missing. This writes to the intended file when it can, and
# to "<name>_NEW.csv" when it cannot, warning loudly instead of halting.
safe_write_csv <- function(x, path, ...) {
  if (!want_write(basename(path))) {
    cat("skipped: ", basename(path), "\n", sep = ""); return(invisible(NULL))
  }
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
  if (!want_write(filename)) { cat("skipped: ", filename, "\n", sep = ""); return(invisible(p)) }
  g   <- ggplot2::ggplotGrob(p)
  pan <- g$layout[g$layout$name == "panel", , drop = FALSE]
  g$widths[unique(pan$l)]  <- grid::unit(PANEL_W_MM * SCALE, "mm")
  g$heights[unique(pan$t)] <- grid::unit(PANEL_H_MM * SCALE, "mm")
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
  cat(sprintf("written: %-30s panel %.4f x %.4f mm  (= %.4f x %.4f at 1/%g)  page %.1f x %.1f mm
",
              filename, PANEL_W_MM * SCALE, PANEL_H_MM * SCALE,
              PANEL_W_MM, PANEL_H_MM, SCALE, w, h))
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

draw_volcano <- function(df, keep, title_text, lab_neg, lab_pos, show_excluded = FALSE) {
  df$neg_log10 <- cap_y(df$neg_log10)
  df$lfc       <- cap_x(df$lfc)
  df$PanelKey <- ifelse(!df$is_sig, "Not Significant",
                 ifelse(!keep,      "Excluded by strict definition",
                 ifelse(df$is_chemo, "Chemosensory", df$GO_Domain)))
  lev <- c("Not Significant", "Excluded by strict definition",
           "MF", "CC", "BP", "Unknown", "Chemosensory")
  brk <- c("Chemosensory", "MF", "CC", "BP", "Unknown",
           if (show_excluded) "Excluded by strict definition", "Not Significant")
  lbs <- c("Chemosensory", "Molecular function", "Cellular component",
           "Biological process", "Unknown",
           if (show_excluded) "Significant vs virgin females only", "Not Significant")
  pal <- c(PAL_GO, `Excluded by strict definition` = "#BFBFBF")

  df$.labsig  <- df$is_sig & keep        # store BEFORE reordering, see note
  df$PanelKey <- factor(df$PanelKey, levels = lev)
  df <- df[order(df$PanelKey), ]          # chemosensory drawn last, on top

  # NOTE: `keep` is an argument in the ORIGINAL row order. Reordering df above
  # would misalign it, so the mask is carried as a column and reorders with df.
  lab <- build_labels(df, df$lfc, "Gene", sig = df$.labsig)

  ggplot(df, aes(x = lfc, y = neg_log10, colour = PanelKey)) +
    geom_point(size = 0.8 * SCALE, alpha = 0.6) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed",
               linewidth = 0.3 * SCALE) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
               linewidth = 0.3 * SCALE) +
    volcano_labels(lab) +
    scale_colour_manual(values = pal, breaks = brk, labels = lbs,
                        na.value = "#D9D9D9", drop = FALSE) +
    volcano_y_scale() +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6,
             label = lab_neg, size = 2.6 * SCALE, colour = "grey30") +
    annotate("text", x =  Inf, y = Inf, hjust =  1.05, vjust = 1.6,
             label = lab_pos, size = 2.6 * SCALE, colour = "grey30") +
    labs(title = NULL, colour = NULL,   # titles set in Illustrator; counts in the CSV
         x = expression(Log[2](FoldChange)),
         y = expression(-log[10](italic(padj)))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "bottom", legend.text = element_text(size = 7.5),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 3)) +
    volcano_x_scale() + volcano_theme()
}

cat("\n--- Panels 4A and 4B ---\n")
ab_rows <- list()
for (a in APPENDAGES) {
  d <- DAT[[a]]; s <- d$sex
  n_f <- sum(s$is_sig & s$lfc > 0); n_m <- sum(s$is_sig & s$lfc < 0)
  ch_s <- sum(s$is_sig & s$is_chemo)
  p <- draw_volcano(s, rep(TRUE, nrow(s)),
                    paste0(a, " \u2014 sex bias in virgins : ", ch_s, "/", n_f + n_m),
                    paste0("male-biased\n", n_m), paste0("female-biased\n", n_f))
  print(p); save_volcano(p, MAIN_OUT, paste0("Figure_4A_", gsub(" ", "_", a), ".pdf"))

  m <- d$mating; keep <- m$Gene %in% SETS[[a]]$mating_any
  n_up <- sum(keep & m$is_sig & m$lfc > 0); n_dn <- sum(keep & m$is_sig & m$lfc < 0)
  ch_m <- sum(keep & m$is_sig & m$is_chemo); n_drop <- sum(m$is_sig & !keep)
  p <- draw_volcano(m, keep,
                    paste0(a, " \u2014 mating-responsive : ", ch_m, "/", n_up + n_dn),
                    paste0("mating-suppressed\n", n_dn), paste0("mating-induced\n", n_up),
                    show_excluded = TRUE)
  print(p); save_volcano(p, MAIN_OUT, paste0("Figure_4B_", gsub(" ", "_", a), ".pdf"))

  ab_rows[[a]] <- data.frame(Appendage = a,
    sex_total = n_f + n_m, sex_female = n_f, sex_male = n_m, sex_chemo = ch_s,
    mating_strict = n_up + n_dn, mating_up = n_up, mating_down = n_dn,
    mating_chemo = ch_m, dropped_by_strict = n_drop, stringsAsFactors = FALSE)
  cat(sprintf("  %-15s sex %d (%d male / %d female) | mating strict %d (%d up / %d down), dropped %d\n",
              a, n_f + n_m, n_m, n_f, n_up + n_dn, n_up, n_dn, n_drop))
}
AB <- bind_rows(ab_rows)
safe_write_csv(AB, file.path(MAIN_OUT, "Figure_4AB_counts.csv"), row.names = FALSE)

# =============================================================================
# PANEL 4C — three-way Venn across appendages
# =============================================================================
venn3 <- function(A, B, C) list(
  a = setdiff(setdiff(A, B), C), b = setdiff(setdiff(B, A), C),
  c = setdiff(setdiff(C, A), B), ab = setdiff(intersect(A, B), C),
  ac = setdiff(intersect(A, C), B), bc = setdiff(intersect(B, C), A),
  abc = Reduce(intersect, list(A, B, C)))

draw_venn3 <- function(title, subtitle, colour, regions) {
  cnt <- sapply(regions, length)
  circles <- data.frame(x0 = c(0, -1.5, 1.5), y0 = c(1.8, 0, 0), r = rep(1.5, 3))
  nums <- data.frame(x = c(0, -2.0, 2.0, -1.0, 1.0, 0, 0),
                     y = c(3.0, -0.2, -0.2, 1.2, 1.2, -0.4, 0.9),
                     txt = as.character(cnt[c("a","b","c","ab","ac","bc","abc")]))
  labs_df <- data.frame(x = c(0, -2.5, 2.5), y = c(3.8, -1.8, -1.8), txt = APPENDAGES)
  ggplot() +
    geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r),
                fill = adjustcolor(colour, alpha.f = 0.2), colour = NA) +
    geom_text(data = nums, aes(x, y, label = txt), size = 5) +
    geom_text(data = labs_df, aes(x, y, label = txt), size = 4,
              colour = colour, fontface = "bold") +
    coord_fixed(xlim = c(-4, 4), ylim = c(-2.5, 4.5)) +
    labs(title = title, subtitle = subtitle) + theme_void() +
    theme(plot.title    = element_text(hjust = .5, face = "bold", size = 11, colour = colour),
          plot.subtitle = element_text(hjust = .5, size = 7.5, colour = "grey30"))
}

cat("\n--- Panel 4C ---\n")
VENN_CFG <- list(
  list(id = "male_biased",   colour = unname(COL_VM["dark"]), title = "Male-biased (virgins)",
       sub = "significant in virgin male vs virgin female"),
  list(id = "female_biased", colour = unname(COL_VF["dark"]), title = "Female-biased (virgins)",
       sub = "significant in virgin female vs virgin male"),
  list(id = "mating_down",   colour = unname(COL_VF["dark"]), title = "Mated-depleted",
       sub = "lower in mated females than BOTH virgin groups"),
  list(id = "mating_up",     colour = unname(COL_MF["dark"]), title = "Female-mating-biased",
       sub = "higher in mated females than BOTH virgin groups"))
c_rows <- list()
for (v in VENN_CFG) {
  A <- SETS[["Antenna"]][[v$id]]; B <- SETS[["Maxillary palp"]][[v$id]]; C <- SETS[["Tarsi"]][[v$id]]
  reg <- venn3(A, B, C)
  p <- draw_venn3(v$title, v$sub, v$colour, reg)
  print(p); save_pdf(p, MAIN_OUT, paste0("Figure_4C_", v$id, ".pdf"), 3.5, 3.6)
  c_rows[[v$id]] <- data.frame(Set = v$id, Antenna = length(A), Maxillary_palp = length(B),
    Tarsi = length(C), Ant_only = length(reg$a), Palp_only = length(reg$b),
    Tarsi_only = length(reg$c), Ant_Palp = length(reg$ab), Ant_Tarsi = length(reg$ac),
    Palp_Tarsi = length(reg$bc), All_three = length(reg$abc), stringsAsFactors = FALSE)
  cat(sprintf("  %-22s all three = %d, antenna&palp = %d\n", v$id, length(reg$abc), length(reg$ab)))
}
safe_write_csv(bind_rows(c_rows), file.path(MAIN_OUT, "Figure_4C_regions.csv"), row.names = FALSE)

# =============================================================================
# PANEL 4D — within-appendage Venn, sex-biased vs mating-responsive
# =============================================================================
cat("\n--- Panel 4D ---\n")
d_rows <- list()
for (a in APPENDAGES) {
  s <- SETS[[a]]; both <- length(intersect(s$sex_any, s$mating_any))
  nl <- length(s$sex_any); nr <- length(s$mating_any)
  circles <- data.frame(x0 = c(-0.7, 0.7), y0 = c(0, 0), r = c(1.35, 1.35),
                        fill = c(unname(COL_VF["dark"]), unname(COL_MF["dark"])))
  nums <- data.frame(x = c(-1.35, 1.35, 0), y = c(0, 0, 0),
                     txt = as.character(c(nl - both, nr - both, both)))
  labs_df <- data.frame(x = c(-1.5, 1.5), y = c(1.85, 1.85),
                        txt = c(paste0("Sex-biased\nin virgins\n(", nl, ")"),
                                paste0("Mating-\nresponsive\n(", nr, ")")),
                        col = c(unname(COL_VF["dark"]), unname(COL_MF["dark"])))
  p <- ggplot() +
    geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = I(fill)),
                colour = "grey45", alpha = 0.25, linewidth = 0.3) +
    geom_text(data = nums, aes(x, y, label = txt), size = 5) +
    geom_text(data = labs_df, aes(x, y, label = txt, colour = I(col)),
              size = 3.4, fontface = "bold", lineheight = 0.9) +
    coord_fixed(xlim = c(-3, 3), ylim = c(-1.8, 2.7)) +
    labs(title = a, subtitle = "mating-responsive = significant vs both virgin groups") +
    theme_void() +
    theme(plot.title    = element_text(hjust = .5, face = "bold", size = 11),
          plot.subtitle = element_text(hjust = .5, size = 7.5, colour = "grey30"))
  print(p); save_pdf(p, MAIN_OUT, paste0("Figure_4D_", gsub(" ", "_", a), ".pdf"), 3.4, 3.4)
  d_rows[[a]] <- data.frame(Appendage = a, sex_biased = nl, mating_responsive = nr,
    sex_only = nl - both, both = both, mating_only = nr - both, stringsAsFactors = FALSE)
  cat(sprintf("  %-15s sex %d, mating %d, overlap %d\n", a, nl, nr, both))
}
safe_write_csv(bind_rows(d_rows), file.path(MAIN_OUT, "Figure_4D_regions.csv"), row.names = FALSE)

# =============================================================================
# PANEL 4E — sex bias in virgins versus sex bias in mated females
# Both axes share virgin males as reference, so the comparison is interpretable
# and a sign flip is a genuine reversal. No correlation coefficient is reported.
# =============================================================================
# Class names say exactly what the category IS: in which group the gene is
# sex-biased, and whether the direction is the same. "Mated only" alone did not
# say mated-only WHAT. The legend title carries the shared "Sex bias" stem so
# the four entries read as one sentence each.
# Colour marks chemosensory genes red, as it does in every other panel of this
# paper. Reversal therefore needs its own channel and gets an open ring: there
# are only 11 reversals, against ~98 chemosensory points, so the ring is the
# cheaper mark to spend on the rarer event.
# Vocabulary discipline, because "in both" was previously doing two different
# jobs - direction in the quadrant labels, significance in the legend:
#   "biased"       -> ALWAYS direction (which sex is higher)
#   "significant"  -> ALWAYS evidence  (did the test pass)
#   "virgins"/"mated" -> ALWAYS the two groups, never "contrast" or "state"
# "Other" is not a category. The complement of "chemosensory" is
# "non-chemosensory", and the qualifier uses the Methods' own vocabulary: P71
# classifies each gene by its behaviour across the two reproductive states, and
# "reproductive state" is the term used throughout the manuscript.
CLS_CHEMO <- "Chemosensory"
CLS_BOTH  <- "Non-chemosensory, sex-biased in both states"
CLS_ONE   <- "Non-chemosensory, sex-biased in one state"
CLS_LEVELS <- c(CLS_CHEMO, CLS_BOTH, CLS_ONE)

cat("\n--- Panel 4E ---\n")
e_dat <- list(); e_stats <- list()
for (a in APPENDAGES) {
  s <- DAT[[a]]$sex[, c("Gene", "Name", "lfc", "padj")]
  n <- DAT[[a]]$matedsex[, c("Gene", "lfc", "padj")]
  d <- merge(s, n, by = "Gene", suffixes = c("_virgin", "_mated"))
  d$x <- d$lfc_virgin      # already negated: + = female-biased in virgins
  d$y <- d$lfc_mated       # already negated: + = female-biased in mated
  d$sig_v <- d$padj_virgin < PADJ_THR & abs(d$x) >= LFC_THR
  d$sig_m <- d$padj_mated  < PADJ_THR & abs(d$y) >= LFC_THR
  d <- d[d$sig_v | d$sig_m, ]
  d$reversed <- d$sig_v & d$sig_m & sign(d$x) != sign(d$y)
  d$Class <- ifelse(nzchar(d$Name), CLS_CHEMO,
             ifelse(d$sig_v & d$sig_m, CLS_BOTH, CLS_ONE))
  d$Class <- factor(d$Class, levels = CLS_LEVELS)
  d$Appendage <- factor(a, levels = APPENDAGES)
  e_dat[[a]] <- d
  e_stats[[a]] <- data.frame(Appendage = a,
    n_virgin = sum(d$sig_v), n_mated = sum(d$sig_m),
    retained = sum(d$sig_v & d$sig_m), reversed = sum(d$reversed),
    median_abs_virgin = median(abs(d$x[d$sig_v])),
    median_abs_mated  = median(abs(d$y[d$sig_v])), stringsAsFactors = FALSE)
  cat(sprintf("  %-15s retained %d of %d, reversals %d\n", a,
              sum(d$sig_v & d$sig_m), sum(d$sig_v), sum(d$reversed)))
}
E <- bind_rows(e_dat); ES <- bind_rows(e_stats)
safe_write_csv(ES, file.path(MAIN_OUT, "Figure_4E_stats.csv"), row.names = FALSE)

# Zoom. 99.97% of the 2,977 points lie within +/-10; a single antennal gene at
# y = -21.8 was stretching both axes to +/-24 and squeezing every other point
# into the middle sixth of each panel. Clamp the view to +/-10. One point is
# affected, so it is drawn as an ordinary dot sitting ON the boundary rather
# than spending a legend key on a mark used once; the CAPTION names it and
# gives its true value. Nothing is dropped.
lim <- 10
E$.oor <- abs(E$x) > lim | abs(E$y) > lim
E$.px  <- pmax(pmin(E$x, lim), -lim)
E$.py  <- pmax(pmin(E$y, lim), -lim)
cat("  clamped to +/-", lim, ": ", sum(E$.oor), " of ", nrow(E), " points\n", sep = "")

# ── Labels ───────────────────────────────────────────────────────────────────
# The panel is about how far a gene moved OFF the identity line, so that is what
# selects the labels: every reversal (the rarest and the ones the Results name),
# plus the largest departures |y - x| in each appendage. One label per locus, and
# genes with no usable name in Additional file 11 or 27 stay unlabelled rather
# than carrying an accession.
LAB_N_CHEMO     <- 3   # chemosensory labels per appendage, by |y - x|
LAB_MAX_SCATTER <- 22  # characters; a facet is only ~50 mm wide
E$dev    <- abs(E$y - E$x)
E$.chemo <- nzchar(E$Name)          # Name is filled only for the seven families
E$.nm    <- ifelse(E$.chemo, pretty_gene(locus_name(E$Name)),
                             shorten(raw_desc(E$Gene)))
# drop names too long to sit in a 50 mm facet rather than truncating them to an
# ellipsis that names nothing
E$.nm[nchar(E$.nm) > LAB_MAX_SCATTER] <- ""
# EVERY reversal must be named - these are the genes the Results single out. Four
# of the eight palp reversals have no description in Additional file 11, no NCBI
# name in 27 and no D. melanogaster hit at any coverage, so they are labelled by
# accession: unlovely, but a reader can look an accession up and cannot look up
# a blank.
.rev <- E$reversed
E$.nm[.rev & !nzchar(E$.nm)] <- E$Gene[.rev & !nzchar(E$.nm)]
E$lab <- ""
for (.a in APPENDAGES) {
  i <- which(E$Appendage == .a & nzchar(E$.nm))
  if (!length(i)) next
  i <- i[order(-E$dev[i])]
  # one label per NAME: three adult cuticle protein paralogs all resolve to
  # "Acp1", and printing it three times in one facet names nothing
  i <- i[!duplicated(E$.nm[i])]
  rv <- i[E$reversed[i]]
  ch <- setdiff(i[E$.chemo[i]], rv)
  sel <- c(rv, head(ch, LAB_N_CHEMO))
  E$lab[sel] <- E$.nm[sel]
}
cat("  reversals: ", sum(.rev), " (", sum(.rev & !E$.chemo), " non-chemosensory)\n", sep = "")
cat("  labelled: ", sum(nzchar(E$lab)), " points (",
    paste(sprintf("%s %d", APPENDAGES,
          tapply(nzchar(E$lab), E$Appendage, sum)[APPENDAGES]), collapse = ", "),
    ")\n", sep = "")

# What each quadrant MEANS. Two layers are at work here and they must not be
# encoded twice: POSITION is direction (both axes are female-minus-male, so the
# sign says which sex is favoured in that group) and COLOUR is evidence (in
# which group the bias is significant). "Reversed" belongs to the colour layer
# only - it is called when a gene is significant in BOTH contrasts with
# opposite signs - so the quadrant text must NOT use the word. Written
# geometrically, each label is true of every point in its quadrant whatever its
# significance class. Drawn once, in the leftmost facet, as a key for all three.
.qtext <- c("Male-biased in virgins,\nfemale-biased in mated females",
            "Female-biased in virgins\nand in mated females",
            "Male-biased in virgins\nand in mated females",
            "Female-biased in virgins,\nmale-biased in mated females")
.qkey  <- c("x- y+", "x+ y+", "x- y-", "x+ y-")
E$Quad <- ifelse(E$x < 0 & E$y > 0, "x- y+",
          ifelse(E$x > 0 & E$y > 0, "x+ y+",
          ifelse(E$x < 0 & E$y < 0, "x- y-",
          ifelse(E$x > 0 & E$y < 0, "x+ y-", "axis"))))
# how many genes sit in each quadrant of each facet. The descriptive wording is
# printed once, in the leftmost facet; the count is printed in every facet.
QUAD <- do.call(rbind, lapply(seq_along(APPENDAGES), function(ai) {
  a <- APPENDAGES[ai]
  n <- vapply(.qkey, function(k) sum(E$Appendage == a & E$Quad == k), integer(1))
  data.frame(Appendage = factor(a, levels = APPENDAGES),
             qx  = c(-lim,  lim, -lim,  lim) * 0.96,
             qy  = c( lim,  lim, -lim, -lim) * 0.96,
             hj  = c(0, 1, 0, 1), vj = c(1, 1, 0, 0),
             lab = if (ai == 1) paste0(.qtext, "\nn = ", n) else paste0("n = ", n),
             stringsAsFactors = FALSE)
}))
# the same counts as a table, for the caption and the record
.qc <- as.data.frame(table(Appendage = E$Appendage, Quadrant = E$Quad, Class = E$Class))
.qc <- .qc[.qc$Freq > 0, ]
safe_write_csv(.qc[order(.qc$Appendage, .qc$Quadrant), ],
               file.path(MAIN_OUT, "Figure_4D_quadrant_counts.csv"), row.names = FALSE)
cat("  plotted: ", nrow(E), " gene-by-appendage points (",
    paste(sprintf("%s %d", APPENDAGES, as.integer(table(E$Appendage)[APPENDAGES])),
          collapse = ", "), ")\n", sep = "")

# draw order: grey backdrop, then blue, then chemosensory red on top, so the
# genes this paper is about are never buried under the other 2,900 points
E <- E[order(match(as.character(E$Class), c(CLS_ONE, CLS_BOTH, CLS_CHEMO))), ]

p4e <- ggplot(E, aes(.px, .py, colour = Class)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey92", alpha = .5) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey92", alpha = .5) +
  geom_hline(yintercept = 0, linewidth = .3) + geom_vline(xintercept = 0, linewidth = .3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = .3, colour = "grey45") +
  geom_point(size = 0.7, alpha = .8) +
  # chemosensory genes ringed and drawn on top: colour already carries the
  # evidence layer, so membership of the gene families this paper is about
  # needs its own channel rather than a fifth colour
  geom_point(data = subset(E, reversed), aes(shape = "Sex-bias reversal"),
             size = 1.5, stroke = .32, colour = "grey10", fill = NA) +
  ggrepel::geom_text_repel(
    data = subset(E, nzchar(lab)), aes(label = lab),
    size = MM_MIN, colour = "#141414", segment.colour = "grey55",
    segment.size = 0.16, min.segment.length = 0, box.padding = 0.34,
    point.padding = 0.16, force = 3, max.iter = 20000, max.overlaps = Inf,
    seed = 1, show.legend = FALSE) +
  geom_text(data = QUAD, aes(qx, qy, label = lab, hjust = hj, vjust = vj),
            size = MM_MIN, colour = "grey45", lineheight = 0.95,
            inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~ Appendage, nrow = 1) +
  scale_colour_manual(values = setNames(c(GO_COLS[["Chemosensory"]], "#4A76C4", "#C4C4C4"),
                                        CLS_LEVELS),
                      breaks = CLS_LEVELS, name = NULL) +
  # Both extra marks share one shape scale so each earns a legend key. The
  # triangle previously had none, which left a reader no way to learn that one
  # point sits off the axis rather than on it.
  scale_shape_manual(values = c("Sex-bias reversal" = 21), name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.3)),
         shape  = guide_legend(order = 2,
                   override.aes = list(colour = "grey15", fill = NA,
                                       size = 1.5, stroke = .32))) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  # Short axis titles on purpose: the rotated y title cannot be longer than the
  # panel is tall or it is clipped. The contrasts are spelt out in the caption.
  labs(x = expression("Sex bias, VF vs Vm ("*log[2]*"FC)"),
       y = expression("Sex bias, MF vs Vm ("*log[2]*"FC)")) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.key.size = unit(3, "mm"), legend.margin = margin(t = -2, unit = "mm"),
        legend.text = element_text(size = 6.2), legend.spacing.x = unit(1, "mm"),
        axis.title = element_text(size = 7.5))
# One output only. This panel used to be written twice - once as
# Figure_4E_revised.pdf at 9.5 x 4.2 in under its old panel letter, once at the
# assembled size - and after the rebuild both came from the same ggplot object,
# so the folder held two names for one figure. The scatter is panel D now.
# 165 mm wide, sized to drop straight into the layout.
# This build is PANEL D of Figure 4: it answers "does sex bias survive mating?"
# (retained / lost / reversed). It is the only analysis in the figure that is
# not another view of the mating-responsive set, which panels A, B and C all
# already show - which is why it, and not the mated-female-specific build
# below, earns the main-figure slot. The latter is Figure S26.
save_pdf(p4e, MAIN_OUT, "Figure_4D_sexbias_retention.pdf", 6.5, 3.10)

# =============================================================================
# PANEL 4D, ALTERNATIVE BUILD — the same axes read as the design intends
#
# Both axes share virgin males as the reference, so
#     y - x = (MF - Vm) - (VF - Vm) = MF - VF
# and the VERTICAL DISTANCE FROM THE DIAGONAL IS THE MATING CONTRAST. Verified
# gene by gene against Additional files 5/7/9: median |difference| 0.005-0.038
# log2 units, Pearson r 0.9985-1.0000; the two differ only by DESeq2 fitting
# each contrast separately.
#
# That makes this the natural place to show the mated-female-specific set - the
# genes that differ from BOTH virgin groups - which is what the three-group
# design was built to find. Position along the diagonal is sex bias common to
# both states; displacement off it is the mating response.
#
# Output: Figure_4D_ALT_matedspecific.pdf, for comparison with the build above.
# =============================================================================
cat("\n--- Panel 4D, alternative build ---\n")
alt_rows <- list()
for (a in APPENDAGES) {
  sx <- DAT[[a]]$sex[, c("Gene", "Name", "lfc", "padj")]
  nx <- DAT[[a]]$matedsex[, c("Gene", "lfc", "padj")]
  d  <- merge(sx, nx, by = "Gene", suffixes = c("_virgin", "_mated"))
  d$x <- d$lfc_virgin; d$y <- d$lfc_mated
  d$sig_v <- d$padj_virgin < PADJ_THR & abs(d$x) >= LFC_THR
  d$sig_m <- d$padj_mated  < PADJ_THR & abs(d$y) >= LFC_THR
  d$up   <- d$Gene %in% SETS[[a]]$mating_up
  d$down <- d$Gene %in% SETS[[a]]$mating_down
  # plot anything sex-biased in either state OR mated-female-specific
  d <- d[d$sig_v | d$sig_m | d$up | d$down, ]
  d$Appendage <- factor(a, levels = APPENDAGES)
  alt_rows[[a]] <- d
}
A <- bind_rows(alt_rows)
# The criterion belongs in the legend TITLE, so each key can state a direction
# and nothing else. "Sex-biased only" said only-what; "Not mating-responsive"
# states the actual test result for those genes.
MFS_UP   <- "Induced in mated females"
MFS_DOWN <- "Reduced in mated females"
MFS_NO   <- "Not mating-responsive"
A$MFS <- ifelse(A$up, MFS_UP, ifelse(A$down, MFS_DOWN, MFS_NO))
A$MFS <- factor(A$MFS, levels = c(MFS_UP, MFS_DOWN, MFS_NO))
A$.chemo <- nzchar(A$Name)
A$.px <- pmax(pmin(A$x, lim), -lim); A$.py <- pmax(pmin(A$y, lim), -lim)
for (.a in APPENDAGES) {
  k <- A$Appendage == .a
  ki <- k & A$up; kr <- k & A$down; kn <- k & !A$up & !A$down
  cat(sprintf("  %-15s induced %4d (%2d chemo) | reduced %4d (%2d chemo) | not responsive %5d (%2d chemo) | all %5d (%2d chemo)\n",
      .a, sum(ki), sum(ki & A$.chemo), sum(kr), sum(kr & A$.chemo),
      sum(kn), sum(kn & A$.chemo), sum(k), sum(k & A$.chemo)))
}
cat("  plotted: ", nrow(A), " points; mated-female-specific ",
    sum(A$up | A$down), " (", paste(sprintf("%s %d", APPENDAGES,
      vapply(APPENDAGES, function(a) sum(A$Appendage == a & (A$up | A$down)), integer(1))),
      collapse = ", "), ")\n", sep = "")

# label the chemosensory members of the mated-female-specific set, and the
# largest displacements from the diagonal in each appendage
A$dev  <- abs(A$y - A$x)
A$.nm  <- ifelse(A$.chemo, pretty_gene(locus_name(A$Name)), shorten(raw_desc(A$Gene)))
A$.nm[nchar(A$.nm) > LAB_MAX_SCATTER] <- ""
A$lab <- ""
for (.a in APPENDAGES) {
  i <- which(A$Appendage == .a & nzchar(A$.nm) & (A$up | A$down))
  if (!length(i)) next
  i <- i[order(-A$dev[i])]; i <- i[!duplicated(A$.nm[i])]
  ch <- i[A$.chemo[i]]
  sel <- c(ch, head(setdiff(i, ch), 4))
  A$lab[sel] <- A$.nm[sel]
}
A <- A[order(match(as.character(A$MFS), c("Sex-biased only",
        "Mated-female-specific, reduced", "Mated-female-specific, induced"))), ]

# The diagonal is the null for the mating contrast: on it, a gene's sex bias is
# identical in the two states, so MF = VF. Say that, and say what displacement
# from it measures - otherwise the reader has no way to know the panel contains
# the mating comparison at all.
DIAG <- data.frame(
  Appendage = factor(APPENDAGES[1], levels = APPENDAGES),
  # Only the on-line label is kept: the "above/below the line" wording moved to
  # the caption once the per-facet count table took the corners.
  qx = -lim * 0.60, qy = -lim * 0.60,
  hj = 0.5, vj = -0.55, ang = 45,
  lab = "y = x, no mating response",
  stringsAsFactors = FALSE)

# Per-facet counts. Each panel should state its own result rather than sending
# the reader to the text: how many genes are mated-female-specific out of those
# plotted, and how the set splits by direction. In the leftmost facet these sit
# lower, to clear the "above the line" annotation.
CNT <- do.call(rbind, lapply(seq_along(APPENDAGES), function(ai) {
  a  <- APPENDAGES[ai]
  k  <- A$Appendage == a
  nu <- sum(k & A$up); nd <- sum(k & A$down); nt <- sum(k)
  nns <- nt - nu - nd                       # grey: no significant mating response
  # chemosensory counted WITHIN each category, not just as a total: the question
  # a reader asks of this panel is how many of the induced genes are receptors
  ki  <- k & A$up; kr <- k & A$down; kn <- k & !A$up & !A$down
  f   <- function(x) formatC(x, big.mark = ",", format = "d")
  data.frame(Appendage = factor(a, levels = APPENDAGES),
             qx = lim * 0.97,          # lower-right: empty in all three facets
             qy = -lim * 0.97,
             lab = paste(
               sprintf("%-15s %6s %7s", "",               "genes", "chemos."),
               sprintf("%-15s %6s %7d", "induced",        f(sum(ki)), sum(ki & A$.chemo)),
               sprintf("%-15s %6s %7d", "reduced",        f(sum(kr)), sum(kr & A$.chemo)),
               sprintf("%-15s %6s %7d", "not responsive", f(sum(kn)), sum(kn & A$.chemo)),
               sprintf("%-15s %6s %7d", "all plotted",    f(nt),      sum(k & A$.chemo)),
               sep = "\n"),
             stringsAsFactors = FALSE)
}))

p4d_alt <- ggplot(A, aes(.px, .py, colour = MFS)) +
  geom_abline(slope = 1, intercept = 0, linewidth = .4, colour = "grey35") +
  geom_hline(yintercept = 0, linewidth = .25, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "grey75") +
  geom_point(size = 0.7, alpha = .85) +
  geom_point(data = subset(A, .chemo), aes(shape = "Chemosensory gene"),
             size = MM_MIN, stroke = .32, colour = GO_COLS[["Chemosensory"]], fill = NA) +
  geom_text(data = DIAG, aes(qx, qy, label = lab, hjust = hj, vjust = vj,
                             angle = ang),
            size = MM_MIN, colour = "grey40", lineheight = .95,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = CNT, aes(qx, qy, label = lab), hjust = 1, vjust = 0,
            size = MM_MIN, colour = "grey20", lineheight = 1.05, family = "mono",
            inherit.aes = FALSE, show.legend = FALSE) +
  ggrepel::geom_text_repel(data = subset(A, nzchar(lab)), aes(label = lab),
    size = MM_MIN, colour = "#141414", segment.colour = "grey55",
    segment.size = .16, min.segment.length = 0, box.padding = .34,
    point.padding = .16, force = 3, max.iter = 20000, max.overlaps = Inf,
    seed = 1, show.legend = FALSE) +
  facet_wrap(~ Appendage, nrow = 1) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  scale_colour_manual(values = setNames(c("#E08214", "#2166AC", "#D4D4D4"),
                                        c(MFS_UP, MFS_DOWN, MFS_NO)),
                      breaks = c(MFS_UP, MFS_DOWN, MFS_NO),
                      name = "Mating response (vs both virgin groups):") +
  scale_shape_manual(values = c("Chemosensory gene" = 21), name = NULL) +
  guides(colour = guide_legend(order = 1, title.position = "top",
                               override.aes = list(size = 1.4)),
         shape  = guide_legend(order = 2,
                   override.aes = list(colour = GO_COLS[["Chemosensory"]],
                                       fill = NA, size = 1.5, stroke = .32))) +
  labs(x = expression("Sex bias, VF vs Vm ("*log[2]*"FC)"),
       y = expression("Sex bias, MF vs Vm ("*log[2]*"FC)"),
       # Over the plot rather than inside it: the corners are taken by the count
       # tables, and this applies to all three facets equally.
       subtitle = paste("Above the diagonal: higher in mated than in virgin females;",
                        "below: lower. Distance from the line =",
                        "log2FC, mated vs virgin females.")) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.key.size = unit(3, "mm"), legend.margin = margin(t = -2, unit = "mm"),
        legend.text = element_text(size = 6.2), axis.title = element_text(size = 7.5),
        legend.title = element_text(size = 6.4), legend.box = "horizontal",
        legend.box.just = "bottom",
        plot.subtitle = element_text(size = 6, colour = "grey25",
                                     margin = margin(b = 1.5, unit = "mm")))
save_pdf(p4d_alt + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
         "Figure_S26_matedfemale_specific.pdf", SUPP_W_IN, SUPP_H_IN)

# =============================================================================
# SUPPLEMENTARY FIGURES S20-S24 — GO over-representation
#
# Replaces ranking by raw gene count (C20, C22, C61). One-sided Fisher's exact
# test per term against all annotated genes DESeq2 tested in the same contrast,
# Benjamini-Hochberg corrected within each panel. Bars that do not reach
# q < 0.05 are drawn pale so a long bar can no longer be mistaken for a result.
# =============================================================================
ora <- function(gene_set, universe, dom, top_n) {
  ann <- GO_ANN[GO_ANN$GO_Domain == dom, ]
  key_set  <- toupper(gene_set)
  in_set   <- ann[ann$JoinKey %in% key_set, ]
  in_univ  <- ann[ann$JoinKey %in% toupper(universe), ]
  n <- nrow(in_set); N <- nrow(in_univ)
  if (n == 0 || N == 0) return(NULL)
  terms <- unique(in_set$GO_Name)
  res <- lapply(terms, function(t) {
    k <- sum(in_set$GO_Name  == t)
    K <- sum(in_univ$GO_Name == t)
    if (K == 0) return(NULL)
    ft <- fisher.test(matrix(c(k, n - k, K - k, N - n - (K - k)), 2, byrow = TRUE),
                      alternative = "greater")
    data.frame(GO_Name = t, k = k, n = n, K = K, N = N,
               odds_ratio = unname(ft$estimate), p = ft$p.value, stringsAsFactors = FALSE)
  })
  res <- bind_rows(res)
  if (!nrow(res)) return(NULL)
  res$q <- p.adjust(res$p, method = "BH")     # within this panel
  res <- res[order(res$q, -res$odds_ratio), ]
  head(res, top_n)
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

make_supp_go <- function(sets_by_app, univ_by_app, fignum, subtitle, top_n) {
  rows <- list()
  for (a in APPENDAGES) for (dom in c("MF", "CC", "BP")) {
    r <- ora(sets_by_app[[a]], univ_by_app[[a]], dom, top_n)
    if (is.null(r)) next
    r$Appendage <- a; r$Domain <- dom
    rows[[paste(a, dom)]] <- r
  }
  tbl <- bind_rows(rows)
  if (!nrow(tbl)) { cat("  ", fignum, ": nothing to plot\n"); return(invisible(NULL)) }

  tbl$Appendage <- factor(tbl$Appendage, levels = APPENDAGES)
  tbl$Domain    <- factor(tbl$Domain, levels = c("MF", "CC", "BP"),
                          labels = unname(DOMAIN_FULL[c("MF", "CC", "BP")]))
  tbl$score <- -log10(pmax(tbl$q, .Machine$double.xmin))
  tbl$sig   <- tbl$q < Q_THR
  # unique y per facet cell, ordered by q within the cell
  # Nine facets cannot share one y scale: each holds a different term at the same
  # integer position, so a single label set is drawn over all of them. Make y a
  # factor of facet + term and strip the sort prefix in the labels.
  tbl <- tbl %>% group_by(Appendage, Domain) %>% arrange(q, .by_group = TRUE) %>%
    mutate(rank = rev(seq_len(dplyr::n()))) %>% ungroup()
  tbl$ykey <- sprintf("%d@@%d@@%02d@@%s", as.integer(tbl$Appendage),
                      as.integer(tbl$Domain), tbl$rank, tbl$GO_Name)
  tbl$ykey <- factor(tbl$ykey,
    levels = tbl$ykey[order(as.integer(tbl$Appendage), as.integer(tbl$Domain), tbl$rank)])
  tbl$fillkey <- paste0(as.character(tbl$Domain), ifelse(tbl$sig, "", " (n.s.)"))
  tbl$fillkey <- factor(tbl$fillkey,
    levels = as.vector(rbind(DOM_LAB, paste0(DOM_LAB, " (n.s.)"))))

  p <- ggplot(tbl, aes(x = score, y = ykey, fill = fillkey)) +
    geom_col(width = 0.7, position = "identity") +
    geom_vline(xintercept = -log10(Q_THR), linetype = "dashed",
               linewidth = 0.4, colour = "#C00000") +
    geom_text(aes(label = k), hjust = -0.25, size = MM_MIN, colour = "grey25") +
    # facet_grid frees the y scale per ROW, not per cell, so every panel in a row
    # carried all three domains' labels. facet_wrap gives nine independent panels.
    facet_wrap(~ Appendage + Domain, ncol = 3, scales = "free",
               labeller = function(d) { d$lab <- paste(d$Appendage, "\u2014", d$Domain); d["lab"] }) +
    scale_y_discrete(labels = function(x) wrap_lab(sub("^\\d+@@\\d+@@\\d+@@", "", x), 40),
                     expand = expansion(add = 0.6)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    scale_fill_manual(values = DOM_FILL, guide = "none") +
    labs(title = paste0("Figure ", fignum, ". ", wrap_lab(subtitle, 95)),
         subtitle = paste0("One-sided Fisher's exact test against all annotated genes tested in the ",
                           "same contrast, Benjamini-Hochberg corrected within each panel. Bar = ",
                           "-log10(q); numbers are gene counts; dashed line marks q = ", Q_THR,
                           " and pale bars did not reach it. DEG thresholds: adjusted p < ",
                           PADJ_THR, ", |log2FC| >= 1."),
         x = expression(-log[10]*"(BH-adjusted "*italic(p)*")"), y = NULL) +
    theme_bw(base_size = 8) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 6.4, colour = "grey30"),
          axis.text.y   = element_text(size = PT_MIN),
          strip.text    = element_text(size = 6.6, face = "bold"),
          panel.spacing = unit(2.2, "mm"),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  print(p)
  save_pdf(p + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
           paste0("Figure_", fignum, ".pdf"), SUPP_W_IN, SUPP_H_IN)
  n_sig <- sum(tbl$q < Q_THR)
  cat("  ", fignum, ": ", n_sig, " terms at q < ", Q_THR, "\n", sep = "")
  invisible(tbl)
}

cat("\n--- Supplementary figures S20-S24 ---\n")
# Universe per appendage = every gene DESeq2 actually tested in either contrast
UNIV <- lapply(APPENDAGES, function(a)
  union(DAT[[a]]$sex$Gene[!is.na(DAT[[a]]$sex$padj_raw)],
        DAT[[a]]$mating$Gene[!is.na(DAT[[a]]$mating$padj_raw)]))
names(UNIV) <- APPENDAGES

sex_all <- lapply(SETS, `[[`, "sex_any");    names(sex_all) <- APPENDAGES
mat_all <- lapply(SETS, `[[`, "mating_any"); names(mat_all) <- APPENDAGES
overlap  <- lapply(APPENDAGES, function(a) intersect(sex_all[[a]], mat_all[[a]]))
sex_only <- lapply(APPENDAGES, function(a) setdiff(sex_all[[a]], overlap[[which(APPENDAGES == a)]]))
mat_only <- lapply(APPENDAGES, function(a) setdiff(mat_all[[a]], overlap[[which(APPENDAGES == a)]]))
names(overlap) <- names(sex_only) <- names(mat_only) <- APPENDAGES

supp_tables <- list()
supp_tables$S12 <- make_supp_go(sex_all,  UNIV, "S20",
  "Gene Ontology over-representation among sex-biased genes across chemosensory appendages", 20)
supp_tables$S13 <- make_supp_go(mat_all,  UNIV, "S21",
  "Gene Ontology over-representation among mating-responsive genes across chemosensory appendages", 20)
supp_tables$S14 <- make_supp_go(sex_only, UNIV, "S22",
  "Gene Ontology over-representation among sex-biased genes without a mating response", 15)
supp_tables$S15 <- make_supp_go(mat_only, UNIV, "S23",
  "Gene Ontology over-representation among mating-responsive genes without a sex bias", 15)
supp_tables$S16 <- make_supp_go(overlap,  UNIV, "S24",
  "Gene Ontology over-representation among genes regulated by both sex and mating status", 15)

supp_all <- bind_rows(lapply(names(supp_tables), function(f) {
  t <- supp_tables[[f]]; if (is.null(t)) return(NULL); t$Figure <- f; t }))
safe_write_csv(supp_all, file.path(SUPP_OUT, "ORA_S20_S24_sex_and_mating.csv"), row.names = FALSE)

# =============================================================================
# SUPPLEMENTARY FIGURE S25 — the volcano plots, moved out of the main figure
#
# Figure 4 used to carry six volcanoes of ~24,800 points each to show 18, 12, 6,
# 0, 5 and 4 significant chemosensory genes, and a reader could not name one
# gene from them. They are kept here because they show the full distribution
# behind the counts in Figure 4, which the lollipops necessarily do not.
# Left column: sex bias in virgins. Right column: the mating response, with
# significance requiring both virgin contrasts (Methods).
# =============================================================================
cat("\n--- Supplementary figure S25 (volcanoes) ---\n")
S25_YCAP <- 300
v_rows <- list()
for (a in APPENDAGES) {
  sx <- DAT[[a]]$sex
  sx <- data.frame(Gene = sx$Gene, Name = sx$Name, lfc = sx$lfc, padj = sx$padj,
                   sig = sx$is_sig, Contrast = "Sex bias in virgins (VF vs Vm)",
                   stringsAsFactors = FALSE)
  m <- DAT[[a]]$mating; n <- DAT[[a]]$matedsex
  j <- merge(m[, c("Gene","Name","lfc","padj")], n[, c("Gene","lfc","padj")],
             by = "Gene", suffixes = c("", "_n"))
  j$sig <- j$padj < PADJ_THR & abs(j$lfc) >= LFC_THR &
           j$padj_n < PADJ_THR & abs(j$lfc_n) >= LFC_THR & sign(j$lfc) == sign(j$lfc_n)
  mt <- data.frame(Gene = j$Gene, Name = j$Name, lfc = j$lfc, padj = j$padj,
                   sig = j$sig, Contrast = "Mating response (vs both virgin groups)",
                   stringsAsFactors = FALSE)
  d <- rbind(sx, mt); d$Appendage <- a
  v_rows[[a]] <- d
}
V <- bind_rows(v_rows)
V$Appendage <- factor(V$Appendage, levels = APPENDAGES)
V$Contrast  <- factor(V$Contrast,
                levels = c("Sex bias in virgins (VF vs Vm)",
                           "Mating response (vs both virgin groups)"))
V$chemo <- nzchar(V$Name)
V$y <- pmin(-log10(pmax(V$padj, .Machine$double.xmin)), S25_YCAP)
V$Class <- ifelse(!V$sig, "Not significant",
           ifelse(V$chemo, "Significant, chemosensory", "Significant, other"))
V$Class <- factor(V$Class, levels = c("Not significant","Significant, other",
                                      "Significant, chemosensory"))
V <- V[order(as.integer(V$Class)), ]          # chemosensory drawn on top
VX <- 12
V$.px <- pmax(pmin(V$lfc, VX), -VX)
n_out <- sum(abs(V$lfc) > VX | -log10(pmax(V$padj, .Machine$double.xmin)) > S25_YCAP)
cat("  points:", nrow(V), "| clamped to the plot box:", n_out, "\n")
for (a in APPENDAGES) for (cc in levels(V$Contrast)) {
  k <- V$Appendage == a & V$Contrast == cc
  cat(sprintf("  %-15s %-40s sig %4d (chemosensory %2d)\n", a, cc,
      sum(k & V$sig), sum(k & V$sig & V$chemo)))
}
p25 <- ggplot(V, aes(.px, y, colour = Class)) +
  geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed",
             linewidth = .25, colour = "grey55") +
  geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
             linewidth = .25, colour = "grey55") +
  geom_point(size = 0.55, alpha = .75) +
  facet_grid(Appendage ~ Contrast) +
  scale_colour_manual(values = c("Not significant" = "#D8D8D8",
                                 "Significant, other" = "#4A76C4",
                                 "Significant, chemosensory" = GO_COLS[["Chemosensory"]]),
                      name = NULL) +
  scale_x_continuous(limits = c(-VX, VX), breaks = seq(-10, 10, 5)) +
  coord_cartesian(ylim = c(0, S25_YCAP)) +
  guides(colour = guide_legend(override.aes = list(size = 1.6))) +
  labs(x = expression(log[2]*" fold change"),
       y = expression(-log[10]*"(adjusted "*italic(p)*")"),
       title = "Figure S25. Differential expression across appendages, by sex and by mating",
       subtitle = paste("Significance: adjusted p < 0.001 and |log2FC| >= 1. Mating additionally",
                        "requires the same direction against both virgin groups.\nAxes clipped at",
                        "|log2FC| 12 and -log10(q) 300; 38 of 148,968 points are drawn on the boundary.")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 8), legend.key.size = unit(3.5, "mm"),
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7, colour = "grey30"))
save_pdf(p25 + theme(plot.margin = SUPP_MARGIN), SUPP_OUT,
         "Figure_S25.pdf", SUPP_W_IN, SUPP_H_IN)

# =============================================================================
# Verification summary — compare against the header targets
# =============================================================================

# =============================================================================
# PANELS 4A / 4B — ALTERNATIVE VERSION, no GO colouring
#
# Requested 2026-09-17, matching the alternative 3A panels so the two figures
# can be assembled in the same visual language. Chemosensory genes red, every
# other significant gene light blue, non-significant pale grey.
#
# For 4B the strict definition still governs which points count as significant:
# genes that pass mated-vs-virgin-female but fail mated-vs-virgin-male are drawn
# as non-significant here rather than in a separate grey category.
#
# Written alongside the standard panels as Figure_4A_alt_*.pdf / Figure_4B_alt_*.pdf.
# Nothing is overwritten.
#
# Set ALT_BLUE_INCLUDES_NONSIG <- TRUE for literally every non-chemosensory
# point in light blue, including the non-significant cloud.
# =============================================================================
ALT_CHEMO_COL  <- "#A60000"
ALT_OTHER_COL  <- "#9ECAE1"   # light blue
ALT_NONSIG_COL <- "#E4E4E4"
ALT_BLUE_INCLUDES_NONSIG <- FALSE

volcano_alt <- function(df, keep, title_text, lab_neg, lab_pos, labels = TRUE) {
  df$neg_log10 <- cap_y(df$neg_log10)
  df$lfc       <- cap_x(df$lfc)
  sig <- df$is_sig & keep
  df$AltKey <- ifelse(df$is_chemo & sig, "Chemosensory",
               ifelse(sig | ALT_BLUE_INCLUDES_NONSIG,
                      "All other genes", "Not significant"))
  df$.labsig <- sig                      # store BEFORE reordering, see note
  df$AltKey <- factor(df$AltKey,
                      levels = c("Not significant", "All other genes", "Chemosensory"))
  df <- df[order(df$AltKey), ]

  # NOTE: `sig` was computed in the ORIGINAL row order; carry it as a column so
  # it reorders with df, otherwise labels attach to the wrong points.
  lab <- build_labels(df, df$lfc, "Gene", sig = df$.labsig)

  ggplot(df, aes(x = lfc, y = neg_log10, colour = AltKey)) +
    geom_point(size = 0.8 * SCALE, alpha = 0.65) +
    geom_hline(yintercept = -log10(PADJ_THR), linetype = "dashed",
               linewidth = 0.3 * SCALE) +
    geom_vline(xintercept = c(-LFC_THR, LFC_THR), linetype = "dashed",
               linewidth = 0.3 * SCALE) +
    (if (labels) volcano_labels(lab)) +
    scale_colour_manual(
      values = c(`Not significant` = ALT_NONSIG_COL,
                 `All other genes` = ALT_OTHER_COL,
                 Chemosensory      = ALT_CHEMO_COL),
      breaks = c("Chemosensory", "All other genes", "Not significant"),
      drop = FALSE, name = NULL) +
    volcano_y_scale() +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6,
             label = lab_neg, size = 2.6 * SCALE, colour = "grey30") +
    annotate("text", x =  Inf, y = Inf, hjust =  1.05, vjust = 1.6,
             label = lab_pos, size = 2.6 * SCALE, colour = "grey30") +
    labs(title = NULL,   # set titles in Illustrator; counts are in the CSV
         x = expression(Log[2](FoldChange)),
         y = expression(-log[10](italic(padj)))) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
          legend.position = "bottom", legend.text = element_text(size = 8),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    volcano_x_scale() + volcano_theme()
}

cat("\n--- Panels 4A / 4B (alternative: chemosensory vs light blue, no GO) ---\n")
for (a in APPENDAGES) {
  d <- DAT[[a]]; s <- d$sex
  n_f <- sum(s$is_sig & s$lfc > 0); n_m <- sum(s$is_sig & s$lfc < 0)
  ch_s <- sum(s$is_sig & s$is_chemo)
  p <- volcano_alt(s, rep(TRUE, nrow(s)),
                   paste0(a, " — sex bias in virgins : ", ch_s, "/", n_f + n_m),
                   paste0("male-biased\n", n_m), paste0("female-biased\n", n_f))
  print(p); save_volcano(p, MAIN_OUT, paste0("Figure_4A_alt_", gsub(" ", "_", a), ".pdf"))
  q <- volcano_alt(s, rep(TRUE, nrow(s)), "", paste0("male-biased\n", n_m),
                   paste0("female-biased\n", n_f), labels = FALSE)
  save_volcano(q, MAIN_OUT, paste0("Figure_4A_nolab_", gsub(" ", "_", a), ".pdf"))

  m <- d$mating; keep <- m$Gene %in% SETS[[a]]$mating_any
  n_up <- sum(keep & m$is_sig & m$lfc > 0); n_dn <- sum(keep & m$is_sig & m$lfc < 0)
  ch_m <- sum(keep & m$is_sig & m$is_chemo)
  p <- volcano_alt(m, keep,
                   paste0(a, " — mating-responsive : ", ch_m, "/", n_up + n_dn),
                   paste0("mating-suppressed\n", n_dn), paste0("mating-induced\n", n_up))
  print(p); save_volcano(p, MAIN_OUT, paste0("Figure_4B_alt_", gsub(" ", "_", a), ".pdf"))
  q <- volcano_alt(m, keep, "", paste0("mating-suppressed\n", n_dn),
                   paste0("mating-induced\n", n_up), labels = FALSE)
  save_volcano(q, MAIN_OUT, paste0("Figure_4B_nolab_", gsub(" ", "_", a), ".pdf"))
}

cat("\n================ CHECK AGAINST THE TARGETS ================\n")
cat("4A  sex DEGs         : ", paste(AB$sex_total, collapse = " / "),
    "   (target 535 / 555 / 263)\n", sep = "")
cat("4A  male-biased      : ", paste(AB$sex_male, collapse = " / "),
    "   (target 403 / 470 / 95)\n", sep = "")
cat("4B  mating strict    : ", paste(AB$mating_strict, collapse = " / "),
    "   (target 146 / 318 / 33)\n", sep = "")
cat("4B  up / down        : ", paste(AB$mating_up, AB$mating_down, sep = "/", collapse = "  "),
    "   (target 6/140  38/280  27/6)\n", sep = "")
cat("4D  sex-mating overlap: ", paste(sapply(APPENDAGES, function(a)
    length(intersect(SETS[[a]]$sex_any, SETS[[a]]$mating_any))), collapse = " / "),
    "   (target 41 / 168 / 3)\n", sep = "")
cat("4E  retained         : ", paste(ES$retained, collapse = " / "),
    "   (target 398 / 495 / 15)\n", sep = "")
cat("4E  reversals        : ", paste(ES$reversed, collapse = " / "),
    "   (target 0 / 8 / 3)\n", sep = "")
cat("\nIf any line disagrees with its target, stop and report it.\n")
cat("Done: Figure 4 and Supplementary Figures S20-S24.\n")
