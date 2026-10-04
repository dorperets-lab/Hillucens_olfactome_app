################################################################################
# PROTOTYPE — Figure 4 as a ranked lollipop instead of six volcanoes
#
# After Houri-Zeevi et al. 2026 (bioRxiv 2026.09.23.753750) panel 4k: one row per
# unit, a stem from the null to the value, dot colour = category, dot size = how
# much data is behind it. Here: one row per chemosensory LOCUS, stem from 0 to
# the log2 fold change, colour = appendage, size = mean expression in that
# appendage (so a big shift in a barely-expressed gene does not read as loudly as
# a big shift in an abundant one).
#
# Same 37 loci as the matrix prototype, same thresholds, strict mating.
# Output: 02_Main_Figures/panels_new/Figure_4_lollipop_PROTOTYPE.pdf
################################################################################
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
PKG <- "C:/Users/dorpe/OneDrive/Desktop/PhD_Projects/BSF/Updated_GDrive_version/Submission_v5_2026-09"
AF  <- file.path(PKG, "04_Additional_Files")
OUT <- file.path(PKG, "02_Main_Figures", "panels_new")
PADJ_THR <- 0.001; LFC_THR <- 1.0

loc <- read.csv(file.path(AF,"Additional_file_37_transcript_model_to_gene_locus.csv"), stringsAsFactors=FALSE)
loc$sym <- trimws(loc$Hill_gene_symbol)
LOCUS_ID <- setNames(trimws(loc$gene_id), toupper(loc$sym))
.n <- suppressWarnings(as.numeric(gsub("^[^0-9]*","",loc$sym)))
.f <- loc[order(loc$gene_id, ifelse(is.na(.n),Inf,.n), loc$sym),]
LOCUS_REP <- setNames(.f$sym[!duplicated(.f$gene_id)], .f$gene_id[!duplicated(.f$gene_id)])
pretty_gene <- function(x){x<-trimws(x);p<-sub("^([A-Za-z]+).*$","\\1",x);r<-sub("^[A-Za-z]+","",x)
  paste0(toupper(substring(p,1,1)),tolower(substring(p,2)),r)}
fam_of <- function(n){n<-toupper(trimws(n));o<-rep(NA_character_,length(n))
  o[grepl("^TRP",n)]<-"Trp"; o[grepl("^ORCO",n)]<-"Or"
  o[is.na(o)&grepl("^OR",n)]<-"Or"; o[is.na(o)&grepl("^OBP",n)]<-"Obp"
  o[is.na(o)&grepl("^IR",n)]<-"Ir"; o[is.na(o)&grepl("^GR",n)]<-"Gr"
  o[is.na(o)&grepl("^PPK",n)]<-"Ppk"; o[is.na(o)&grepl("^CSP",n)]<-"Csp"; o}

# mean normalised expression per appendage, for dot size
nc <- read.csv(file.path(AF,"Additional_file_10_normalized_counts_all_samples.csv"),
               stringsAsFactors=FALSE, check.names=FALSE)
pref <- list(Antenna="^Ant_", `Maxillary palp`="^P_", Tarsi="^Leg_")
EXPR <- lapply(pref, function(p){
  cols <- grep(p, names(nc), value=TRUE)
  setNames(rowMeans(sapply(nc[cols], function(x) suppressWarnings(as.numeric(x))), na.rm=TRUE),
           trimws(nc$Gene))})

rd <- function(f){d<-read.csv(file.path(AF,f),stringsAsFactors=FALSE)
  d$Gene<-trimws(d$Gene); d$Name<-trimws(ifelse(is.na(d$Name),"",d$Name))
  d$padj<-suppressWarnings(as.numeric(d$padj))
  d$lfc <- -suppressWarnings(as.numeric(d$log2FoldChange)); d[nzchar(d$Name),]}
SEX <- list(Antenna="Additional_file_4_Antenna_sex_VirginFemale_vs_VirginMale.csv",
            `Maxillary palp`="Additional_file_6_MaxillaryPalp_sex_VirginFemale_vs_VirginMale.csv",
            Tarsi="Additional_file_8_Tarsi_sex_VirginFemale_vs_VirginMale.csv")
MAT <- list(Antenna="Additional_file_5_Antenna_mating_MatedFemale_vs_VirginFemale.csv",
            `Maxillary palp`="Additional_file_7_MaxillaryPalp_mating_MatedFemale_vs_VirginFemale.csv",
            Tarsi="Additional_file_9_Tarsi_mating_MatedFemale_vs_VirginFemale.csv")
NEW <- list(Antenna="Additional_file_29_Antenna_matedsex_MatedFemale_vs_VirginMale.csv",
            `Maxillary palp`="Additional_file_30_MaxillaryPalp_matedsex_MatedFemale_vs_VirginMale.csv",
            Tarsi="Additional_file_31_Tarsi_matedsex_MatedFemale_vs_VirginMale.csv")


# ── Pin the plot box, as for the volcanoes ───────────────────────────────────
# Two of these sit side by side across the 180 mm text width: 90 mm each, of
# which 65 mm is the plot box and the rest is gene names, axis and legend.
BOX_W_MM <- 65

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
  cat(sprintf("written: %-42s box %.0f mm wide, page %.1f x %.1f mm (%d rows)
",
              basename(file), BOX_W_MM, w, h, n_rows))
}

rows <- list()
for (a in names(SEX)) {
  s <- rd(SEX[[a]]); m <- rd(MAT[[a]]); n <- rd(NEW[[a]])
  s$sig <- !is.na(s$padj) & s$padj<PADJ_THR & abs(s$lfc)>=LFC_THR
  mm <- merge(m[,c("Gene","Name","lfc","padj")], n[,c("Gene","lfc","padj")], by="Gene", suffixes=c("","_n"))
  mm$sig <- !is.na(mm$padj)&mm$padj<PADJ_THR&abs(mm$lfc)>=LFC_THR&
            !is.na(mm$padj_n)&mm$padj_n<PADJ_THR&abs(mm$lfc_n)>=LFC_THR&sign(mm$lfc)==sign(mm$lfc_n)
  for (d in list(list(df=s,cn="Sex bias (VF vs Vm)"), list(df=mm,cn="Mating (strict)"))) {
    x <- d$df[d$df$sig, c("Gene","Name","lfc")]
    if(!nrow(x)) next
    x$locus <- unname(LOCUS_ID[toupper(x$Name)]); x$locus[is.na(x$locus)] <- toupper(x$Name[is.na(x$locus)])
    x <- x[order(-abs(x$lfc)),]; x <- x[!duplicated(x$locus),]
    x$expr <- unname(EXPR[[a]][x$Gene])
    rows[[length(rows)+1]] <- data.frame(locus=x$locus, lfc=x$lfc, expr=x$expr,
                                         Appendage=a, Contrast=d$cn, stringsAsFactors=FALSE)
  }
}
D <- bind_rows(rows)
D$gene   <- pretty_gene(ifelse(is.na(LOCUS_REP[D$locus]), D$locus, LOCUS_REP[D$locus]))
D$Family <- fam_of(D$gene)
D$expr[is.na(D$expr)] <- 1
cat("loci:", length(unique(D$locus)), " points:", nrow(D), "\n")

fam_order <- c("Or","Gr","Ir","Obp","Csp","Ppk","Trp")
D$Family <- factor(D$Family, levels=fam_order[fam_order %in% unique(D$Family)])
ord <- D %>% group_by(gene) %>% summarise(f=first(Family), v=lfc[which.max(abs(lfc))], .groups="drop") %>%
       arrange(f, v)
D$gene <- factor(D$gene, levels=ord$gene)
D$Appendage <- factor(D$Appendage, levels=names(SEX), labels=c("Antenna","Maxillary palp","Tarsi"))
D$Contrast  <- factor(D$Contrast, levels=c("Sex bias (VF vs Vm)","Mating (strict)"))
# Two blues and a grey. The colour encodes the result: antenna and palp are the
# olfactory pair that respond to mating together (95 shared genes), so they take
# related blues; the tarsi respond independently and take a neutral grey.
APP_COL <- c(Antenna="#08306B", `Maxillary palp`="#4292C6", Tarsi="#737373")

# Say what each side of zero means, as panel B does. facet_grid puts the two
# contrasts in COLUMNS, so the text goes in the top family row of each column
# and renders in the gap opened under the strip (clip = "off"). Pinning it to
# y = Inf rather than expanding the y scale matters: space = "free_y" would
# hand that expansion to every family panel, padding the one-gene families.
ann <- data.frame(
  Family   = factor(levels(D$Family)[1], levels = levels(D$Family)),
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

p <- ggplot(D, aes(x=lfc, y=gene)) +
  geom_vline(xintercept=0, colour="grey60", linewidth=0.3) +
  geom_segment(aes(x=0, xend=lfc, yend=gene, colour=Appendage), linewidth=0.35) +
  geom_point(aes(colour=Appendage, size=log10(expr+1))) +
  geom_text(data=ann, aes(x=x, y=Inf, label=lab, hjust=hj), vjust=-0.7,
            size=1.9, colour="grey30", inherit.aes=FALSE) +
  facet_grid(Family ~ Contrast, scales="free_y", space="free_y", switch="y") +
  coord_cartesian(clip="off") +
  scale_colour_manual(values=APP_COL, name=NULL) +
  scale_size_continuous(range=c(0.8,3.2), limits=SIZE_LIM, breaks=SIZE_BREAKS,
                        name=SIZE_NAME) +
  scale_x_continuous(limits=c(-9,9), breaks=seq(-8,8,4)) +
  labs(x=expression(log[2]*" fold change"), y=NULL,
       title="Chemosensory loci responding to sex or mating",
       subtitle=paste("One row per locus; stem from no change. Positive = female-biased or",
                      "mating-induced.\nDot size = mean expression in that appendage.",
                      "Only significant results shown (adj. p < 0.001, |log2FC| >= 1).")) +
  theme_bw(base_size=8) +
  theme(panel.grid.major.y=element_line(colour="grey94", linewidth=0.25),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        strip.background.y=element_rect(fill="grey92", colour=NA),
        strip.background.x=element_blank(),
        strip.text.y.left=element_text(angle=0, face="bold", size=7),
        strip.text.x=element_text(face="bold", size=7,
                                  margin=margin(t=1, b=9, unit="pt")),
        strip.placement="outside",
        axis.text.y=element_text(size=5.6),
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=6.2, colour="grey30"),
        legend.position="bottom", legend.box="vertical", legend.spacing.y=unit(0,"mm"),
        legend.key.size=unit(3,"mm"), panel.spacing.y=unit(0.6,"mm"))
pin_and_save(p, file.path(OUT, "Figure_4A_chemosensory_lollipop.pdf"), nlevels(D$gene))
