################################################################################
# Figure 4, panel B — non-chemosensory genes responding to sex or mating
#
# SELECTION RULE (state this in the caption, it is not arbitrary):
#   1. non-chemosensory (no entry in the Name column of Additional files 4-9)
#   2. significant at padj < 0.001 and |log2FC| >= 1
#   3. must have a usable name (Additional file 11 description, else the NCBI
#      name in Additional file 27, else the D. melanogaster BLASTp match at
#      >= 50% query coverage, marked *); unnamed genes are dropped, not blank
#   4. transcript variants of one gene are collapsed to their best-scoring model
#   5. of those, the top TOP_N by pi-score |log2FC| x -log10(padj), computed
#      SEPARATELY within each contrast
#
# Mating uses the STRICT definition, matching every other Figure 4 panel.
# Box pinned to 65 mm so this sits beside panel A across the 180 mm text width.
#
# Output: 02_Main_Figures/panels_new/Figure_4B_other_lollipop.pdf
#         02_Main_Figures/panels_new/Figure_4B_name_key.csv   <- abbreviation key
################################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
PKG <- "C:/Users/dorpe/OneDrive/Desktop/PhD_Projects/BSF/Updated_GDrive_version/Submission_v5_2026-09"
AF  <- file.path(PKG,"04_Additional_Files"); OUT <- file.path(PKG,"02_Main_Figures","panels_new")
PADJ_THR <- 0.001; LFC_THR <- 1.0; TOP_N <- 20; MAXCHAR <- 24; BOX_W_MM <- 65

nm1 <- read.csv(file.path(AF,"Additional_file_11_GO_annotations_all_genes.csv"), stringsAsFactors=FALSE)
D1 <- setNames(trimws(sub("\\s*\\[Source:.*$","",nm1$description)), toupper(trimws(nm1$Gene)))
nm2 <- tryCatch(read.csv(file.path(AF,"Additional_file_27_Unknown_BLAST.csv"), stringsAsFactors=FALSE),
                error=function(e) NULL)
D2 <- if (is.null(nm2)) character(0) else
  setNames(trimws(nm2$NCBI_gene_name), toupper(trimws(nm2$Gene)))
# Third fallback: best D. melanogaster BLASTp match from Additional file 27, but
# ONLY at >= 50% query coverage (the rule the Methods state) and only where the
# Dmel protein has a real name - "uncharacterized protein Dmel_CG7214" and clone
# identifiers (ORF, AT10315p, IP06786p, E239.9-1) say no more than the accession
# does. These names are an INFERENCE from one BLAST hit, not a BSF annotation,
# so they are marked * and the caption must say so.
DMEL_MIN_COV <- 50
D3 <- character(0)
if (!is.null(nm2)) {
  cov <- suppressWarnings(as.numeric(nm2$Dmel_query_coverage_pct))
  dm  <- nm2[!is.na(cov) & cov >= DMEL_MIN_COV, , drop = FALSE]
  if (nrow(dm)) {
    # trimws AFTER stripping the pipe prefix: without it a leading space made
    # "^uncharacterized" fail to match and " uncharacterized protein
    # Dmel_CG9782*" was plotted as if it were a gene name.
    v <- trimws(sub("^.*\\|", "", as.character(dm$Dmel_gene_symbol)))
    v[grepl("^uncharacterized", v, ignore.case = TRUE)] <- ""
    v[grepl("^Dmel_CG[0-9]+$", v)] <- ""
    # clone / EST identifiers: two letters + digits + "p", plus the named forms
    v[grepl("^([A-Z]{2}[0-9]{4,}p|ORF|E[0-9.]+-?[0-9]*)$", v)] <- ""
    k <- toupper(trimws(dm$Gene)); keep <- nzchar(v) & nzchar(k)
    D3 <- setNames(v[keep], k[keep])[!duplicated(k[keep])]
  }
}

# ── Name shortening ── shared with Figure_4.R (panel D), so the two panels of
# one figure cannot disagree about what a gene is called.
source(file.path(PKG, "05_R_Scripts", "Figure_4_names.R"))

label_of <- function(g){
  k <- toupper(trimws(g))
  a <- unname(D1[k]); a[is.na(a)] <- ""
  a[grepl("^uncharacterized", a, ignore.case=TRUE)] <- ""
  b <- unname(D2[k]); b[is.na(b)] <- ""
  b[grepl("^uncharacterized", b, ignore.case=TRUE)] <- ""
  cc <- unname(D3[k]); cc[is.na(cc)] <- ""
  full  <- ifelse(nzchar(a), a, ifelse(nzchar(b), b, cc))
  star  <- !nzchar(a) & !nzchar(b) & nzchar(cc)
  short <- ifelse(nzchar(full), shorten(full), "")
  short <- ifelse(star & nzchar(short), paste0(short, "*"), short)
  short <- ifelse(nchar(short) > MAXCHAR, paste0(substr(short,1,MAXCHAR-1),"\u2026"), short)
  data.frame(lab = short, full = full, stringsAsFactors = FALSE)
}
rd <- function(f){
  d <- read.csv(file.path(AF,f), stringsAsFactors=FALSE)
  d$Gene <- trimws(d$Gene); d$Name <- trimws(ifelse(is.na(d$Name),"",d$Name))
  d$padj <- suppressWarnings(as.numeric(d$padj))
  d$lfc  <- -suppressWarnings(as.numeric(d$log2FoldChange))
  d[!nzchar(d$Name), ]                               # NON-chemosensory only
}
APP <- c("Antenna","Maxillary palp","Tarsi")
SEX <- c("Additional_file_4_Antenna_sex_VirginFemale_vs_VirginMale.csv",
         "Additional_file_6_MaxillaryPalp_sex_VirginFemale_vs_VirginMale.csv",
         "Additional_file_8_Tarsi_sex_VirginFemale_vs_VirginMale.csv")
MAT <- c("Additional_file_5_Antenna_mating_MatedFemale_vs_VirginFemale.csv",
         "Additional_file_7_MaxillaryPalp_mating_MatedFemale_vs_VirginFemale.csv",
         "Additional_file_9_Tarsi_mating_MatedFemale_vs_VirginFemale.csv")
NEW <- c("Additional_file_29_Antenna_matedsex_MatedFemale_vs_VirginMale.csv",
         "Additional_file_30_MaxillaryPalp_matedsex_MatedFemale_vs_VirginMale.csv",
         "Additional_file_31_Tarsi_matedsex_MatedFemale_vs_VirginMale.csv")
nc <- read.csv(file.path(AF,"Additional_file_10_normalized_counts_all_samples.csv"),
               stringsAsFactors=FALSE, check.names=FALSE)
pref <- c(Antenna="^Ant_", `Maxillary palp`="^P_", Tarsi="^Leg_")
EXPR <- lapply(pref, function(p) setNames(
  rowMeans(sapply(nc[grep(p,names(nc),value=TRUE)], function(x) suppressWarnings(as.numeric(x))), na.rm=TRUE),
  trimws(nc$Gene)))

rows <- list()
for (i in seq_along(APP)) {
  a <- APP[i]
  s <- rd(SEX[i]); s$sig <- !is.na(s$padj) & s$padj<PADJ_THR & abs(s$lfc)>=LFC_THR
  m <- rd(MAT[i]); n <- rd(NEW[i])
  mm <- merge(m[,c("Gene","lfc","padj")], n[,c("Gene","lfc","padj")], by="Gene", suffixes=c("","_n"))
  mm$sig <- !is.na(mm$padj)&mm$padj<PADJ_THR&abs(mm$lfc)>=LFC_THR&
            !is.na(mm$padj_n)&mm$padj_n<PADJ_THR&abs(mm$lfc_n)>=LFC_THR&sign(mm$lfc)==sign(mm$lfc_n)
  for (d in list(list(df=s,cn="Sex bias (VF vs Vm)"), list(df=mm,cn="Mating (strict)"))) {
    x <- d$df[d$df$sig, c("Gene","lfc","padj")]
    if (!nrow(x)) next
    x$Appendage <- a; x$Contrast <- d$cn
    x$expr <- unname(EXPR[[a]][x$Gene]); rows[[length(rows)+1]] <- x
  }
}
D <- bind_rows(rows)
D$score <- abs(D$lfc) * -log10(pmax(D$padj, 1e-300))

# rule 3: must have a name; unnamed genes are removed, not shown blank
.L <- label_of(D$Gene); D$lab <- .L$lab; D$full <- .L$full
D <- D[nzchar(D$lab), ]
# rule 4: collapse transcript variants to their best-scoring member. This now
# works, because the variant suffix was stripped BEFORE the names were compared.
bp <- D %>% group_by(lab, Gene) %>% summarise(s = max(score), .groups="drop") %>%
      group_by(lab) %>% slice_max(s, n = 1, with_ties = FALSE) %>% ungroup()
D <- D[D$Gene %in% bp$Gene, ]
# rule 5: top TOP_N SEPARATELY within each contrast, by pi-score
top <- D %>% group_by(Contrast, lab) %>% summarise(best = max(score), .groups="drop") %>%
       group_by(Contrast) %>% slice_max(best, n = TOP_N, with_ties = FALSE) %>% ungroup()
D <- merge(D, top[, c("Contrast","lab")], by = c("Contrast","lab"))
# merge drops the factor class; restore it before it is used for ordering
D$Contrast <- factor(as.character(D$Contrast),
                     levels = c("Sex bias (VF vs Vm)", "Mating (strict)"))
cat("panel B: ", nrow(top), " rows (", TOP_N, " per contrast), ",
    nrow(D), " cells\n", sep="")

# the abbreviation key, for the caption / an Additional file
key <- D %>% distinct(lab, full) %>% arrange(lab)
write.csv(key, file.path(OUT,"Figure_4B_name_key.csv"), row.names=FALSE)
cat("name key written:", nrow(key), "entries\n")

# each contrast has its own genes -> order within facet, free y scales
D <- D %>% group_by(Contrast, lab) %>% mutate(o = lfc[which.max(abs(lfc))]) %>% ungroup()
D$row <- paste(as.integer(D$Contrast), D$lab, sep="@@")
lv <- D %>% distinct(Contrast, lab, o, row) %>% arrange(Contrast, o) %>% pull(row) %>% unique()
D$row <- factor(D$row, levels = lv)
D$Appendage <- factor(D$Appendage, levels=APP)

D$expr[is.na(D$expr)] <- 1
# Two blues and a grey. The colour encodes the result: antenna and palp are the
# olfactory pair that respond to mating together (95 shared genes), so they take
# related blues; the tarsi respond independently and take a neutral grey.
APP_COL <- c(Antenna="#08306B", `Maxillary palp`="#4292C6", Tarsi="#737373")

# Say what each side of zero means, per facet - the two contrasts do not mean
# the same thing by "positive" and a reader should not have to find that in the
# caption.
ann <- data.frame(
  Contrast = factor(rep(levels(D$Contrast), each = 2), levels = levels(D$Contrast)),
  x   = c(-Inf, Inf, -Inf, Inf),
  hj  = c(-0.04, 1.04, -0.04, 1.04),
  # kept short: panel A gives each contrast only half of the 65 mm box, and the
  # longer wording left 0.9 mm between the two mating labels. The strip above
  # already says which contrast this is.
  lab = c("\u2190 male-biased", "female-biased \u2192",
          "\u2190 reduced", "induced \u2192"),
  stringsAsFactors = FALSE)

# ── Size scale, SHARED between panels A and B ───────────────────────────────
# Both panels map log10(mean counts + 1) onto the same 0.8-3.2 mm radius range,
# but each was taking its limits from its own data: the smallest dot meant 23
# counts in A and 113 counts in B, so a dot of a given size did not mean the
# same thing in the two halves of one figure. Fix the limits and the breaks so
# the panels are directly comparable. Keep both in sync if either changes.
SIZE_LIM    <- log10(c(20, 210000) + 1)
SIZE_BREAKS <- 2:5
SIZE_NAME   <- expression(log[10]*"(mean counts+1)")

p <- ggplot(D, aes(x=lfc, y=row)) +
  geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
  geom_segment(aes(x=0, xend=lfc, yend=row, colour=Appendage), linewidth=0.35) +
  geom_point(aes(colour=Appendage, size=log10(expr+1))) +
  geom_text(data=ann, aes(x=x, y=Inf, label=lab, hjust=hj), vjust=1.4,
            size=1.9, colour="grey30", inherit.aes=FALSE) +
  facet_wrap(~ Contrast, ncol=1, scales="free_y") +
  scale_y_discrete(labels=function(x) sub("^[0-9]+@@","",x),
                   expand=expansion(add=c(0.7, 2.4))) +
  scale_colour_manual(values=APP_COL, name=NULL) +
  scale_size_continuous(range=c(0.8,3.2), limits=SIZE_LIM, breaks=SIZE_BREAKS,
                        name=SIZE_NAME) +
  labs(x=expression(log[2]*" fold change"), y=NULL,
       title="Other genes responding to sex or mating",
       subtitle=paste0("Non-chemosensory genes, top ", TOP_N, " within each contrast by ",
                       "|log2FC| x -log10(padj).\nNames abbreviated; full names in the ",
                       "name key. Dot size = mean expression in that appendage.")) +
  theme_bw(base_size=8) +
  theme(panel.grid.major.y=element_line(colour="grey94", linewidth=0.25),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background=element_rect(fill="grey92", colour=NA),
        strip.text=element_text(face="bold", size=7),
        axis.text.y=element_text(size=5.8),
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=6.2, colour="grey30"),
        legend.position="bottom", legend.box="vertical",
        legend.key.size=unit(3,"mm"))


# Pin the WIDTH absolutely; leave the panel height a null unit and give ggsave a
# device height computed from the number of gene rows, so the panel stretches to
# fill it. Pinning both collapsed the panel to 42 mm for 37 rows.
ROW_MM <- 3.2
pin_and_save <- function(p, file, n_rows) {
  g <- ggplot2::ggplotGrob(p)
  pan <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  ncol_pan <- length(unique(pan$l))
  for (col in unique(pan$l)) g$widths[col] <- grid::unit(BOX_W_MM / ncol_pan, "mm")
  grDevices::pdf(NULL)
  w <- grid::convertWidth(sum(g$widths), "mm", valueOnly = TRUE)
  hs <- g$heights
  isnull <- vapply(seq_along(hs), function(i) grid::unitType(hs[i])[1] == "null", logical(1))
  h_fixed <- if (any(!isnull))
    grid::convertHeight(sum(hs[!isnull]), "mm", valueOnly = TRUE) else 0
  grDevices::dev.off()
  h <- h_fixed + n_rows * ROW_MM
  dev_fun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(file, g, width = w, height = h, units = "mm", device = dev_fun)
  cat(sprintf("written: %-42s box %.0f mm wide, page %.1f x %.1f mm (%d rows)\n",
              basename(file), BOX_W_MM, w, h, n_rows))
}

pin_and_save(p, file.path(OUT, "Figure_4B_other_lollipop.pdf"), nlevels(D$row))
