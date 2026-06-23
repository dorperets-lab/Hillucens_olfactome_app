suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(cluster)
  library(vegan)
  library(scales)
})

############################################################
## Settings
############################################################

INPUT_CSV <- "C:/Users/dorpe/Downloads/BSF2026/BSF2026/Supplemetary_Files/csv/BSF_normalized_counts_master.csv"

PERM_N <- 999
SEED <- 123
ELLIPSE_LEVEL <- 0.95
USE_PATCHWORK <- TRUE
OUTPUT_DIR <- "Figure1_R_outputs"

LIBS <- c(
  "Ant_MF1", "Ant_MF2", "Ant_MF3",
  "Ant_VF1", "Ant_VF2", "Ant_VF3",
  "Ant_Vm1", "Ant_Vm2", "Ant_Vm3",
  "Leg_MF1", "Leg_MF2", "Leg_MF3",
  "Leg_VF1", "Leg_VF2", "Leg_VF3",
  "Leg_Vm1", "Leg_Vm2", "Leg_Vm3",
  "P_MF1",   "P_MF2",   "P_MF3",
  "P_VF1",   "P_VF2",   "P_VF3",
  "P_Vm1",   "P_Vm2",   "P_Vm3"
)

GROUP_PATTERNS <- c(
  "Ant_Vm", "Ant_VF", "Ant_MF",
  "Leg_Vm", "Leg_VF", "Leg_MF",
  "P_Vm",   "P_VF",   "P_MF"
)

############################################################
## Report helpers
############################################################

# Force report printing to show complete tables.
# This prevents tibble truncation messages such as:
#   # i 17 more rows
#   # i 4 more variables
# Also explicitly define na.print as a character string.
# This avoids: invalid 'na.print' specification.
options(
  tibble.print_max = Inf,
  tibble.print_min = Inf,
  tibble.width = Inf,
  dplyr.print_max = Inf,
  dplyr.print_min = Inf,
  width = 240,
  max.print = 10000000,
  na.print = "NA"
)

report_lines <- character(0)

add_report_line <- function(...) {
  report_lines <<- c(report_lines, paste0(...))
}

add_report_blank <- function(n = 1) {
  report_lines <<- c(report_lines, rep("", n))
}

add_report_header <- function(title, char = "=") {
  line <- paste(rep(char, max(20, nchar(title))), collapse = "")
  report_lines <<- c(report_lines, line, title, line)
}

add_report_object <- function(title, obj) {
  add_report_header(title, "-")

  printed <- capture.output({
    if (inherits(obj, c("tbl_df", "tbl", "data.frame"))) {
      obj_df <- as.data.frame(obj, stringsAsFactors = FALSE)
      print(obj_df, row.names = FALSE, na.print = "NA", right = FALSE)
    } else if (inherits(obj, "matrix")) {
      obj_df <- as.data.frame(obj, stringsAsFactors = FALSE)
      obj_df <- tibble::rownames_to_column(obj_df, var = "row")
      print(obj_df, row.names = FALSE, na.print = "NA", right = FALSE)
    } else if (inherits(obj, "table")) {
      obj_df <- as.data.frame(obj, stringsAsFactors = FALSE)
      print(obj_df, row.names = FALSE, na.print = "NA", right = FALSE)
    } else if (is.atomic(obj) && !is.null(names(obj))) {
      obj_df <- data.frame(name = names(obj), value = as.vector(obj), stringsAsFactors = FALSE)
      print(obj_df, row.names = FALSE, na.print = "NA", right = FALSE)
    } else {
      print(obj)
    }
  })

  report_lines <<- c(report_lines, printed)
  add_report_blank()
}

add_report_csv_table <- function(title, csv_path) {
  add_report_header(title, "-")
  add_report_line("Source CSV: ", csv_path)
  add_report_blank()

  if (!file.exists(csv_path)) {
    add_report_line("CSV file was not found and could not be added to the report.")
    add_report_blank()
    return(invisible(NULL))
  }

  obj <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  add_report_object("FULL CSV CONTENT", obj)
}

append_exported_csvs_to_report <- function(out_dir = OUTPUT_DIR) {
  add_report_header("SECTION: EXPORTED CSV TABLES - FULL CONTENT")
  add_report_line("This section prints the complete content of every CSV exported by the script.")
  add_report_line("Therefore, the TXT report contains the same detailed tables as the CSV output files.")
  add_report_blank()

  csv_files <- c(
    "00_csv_index_read_me.csv",
    "01_pca_overview.csv",
    "02_pca_group_geometry_and_silhouette.csv",
    "03_permanova_summary.csv",
    "04a_pairwise_permanova_tests.csv",
    "04b_pairwise_mean_distances.csv",
    "05_pearson_correlation_long.csv",
    "06_pearson_correlation_summary.csv",
    "07_variation_partitioning.csv"
  )

  for (csv_file in csv_files) {
    add_report_csv_table(csv_file, file.path(out_dir, csv_file))
  }
}

write_report <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(report_lines, con = path, useBytes = TRUE)
}

############################################################
## Helper 0 - significance codes
############################################################

sig_code <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

############################################################
## Helper 0b - permutation p-value resolution
############################################################

perm_p_floor <- function(n_perm) 1 / (n_perm + 1)

format_perm_p <- function(p, n_perm, digits = 3) {
  floor_p <- perm_p_floor(n_perm)
  if (is.na(p)) return(NA_character_)
  if (p <= floor_p + .Machine$double.eps^0.5) {
    return(paste0("<= ", signif(floor_p, digits)))
  }
  as.character(signif(p, digits))
}

############################################################
## Helper 1 - mean within and between cluster distances
############################################################

compute_within_between <- function(dist_obj, group_factor) {
  dmat <- as.matrix(dist_obj)
  samples <- attr(dist_obj, "Labels")
  groups <- as.factor(group_factor)
  
  if (length(groups) != length(samples)) {
    stop("Group factor length does not match distance matrix size")
  }
  
  names(groups) <- samples
  
  within_vals <- c()
  between_vals <- c()
  
  for (i in seq_along(samples)) {
    for (j in seq_along(samples)) {
      if (j <= i) next
      si <- samples[i]
      sj <- samples[j]
      if (groups[si] == groups[sj]) {
        within_vals <- c(within_vals, dmat[si, sj])
      } else {
        between_vals <- c(between_vals, dmat[si, sj])
      }
    }
  }
  
  list(
    mean_within = mean(within_vals),
    mean_between = mean(between_vals)
  )
}

############################################################
## Helper 2 - pairwise mean distances between groups
############################################################

compute_pairwise_means <- function(dist_obj, group_factor) {
  dmat <- as.matrix(dist_obj)
  samples <- attr(dist_obj, "Labels")
  groups <- as.factor(group_factor)
  names(groups) <- samples
  
  levs <- levels(groups)
  res <- list()
  
  for (i in seq_along(levs)) {
    for (j in seq_along(levs)) {
      if (j <= i) next
      gi <- levs[i]
      gj <- levs[j]
      s_i <- names(groups)[groups == gi]
      s_j <- names(groups)[groups == gj]
      vals <- dmat[s_i, s_j, drop = FALSE]
      res[[length(res) + 1]] <- data.frame(
        group1 = gi,
        group2 = gj,
        mean_dist = mean(vals),
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, res)
}

############################################################
## Helper 2b - 2D area of clusters in PC space
############################################################

polygon_area <- function(x, y) {
  n <- length(x)
  if (n < 3) return(NA_real_)
  sum(x * c(y[-1], y[1]) - y * c(x[-1], x[1])) / 2
}

convex_hull_area <- function(df, xcol = "PC1", ycol = "PC2") {
  df <- df[complete.cases(df[, c(xcol, ycol)]), , drop = FALSE]
  if (nrow(df) < 3) return(NA_real_)
  hull_idx <- chull(df[[xcol]], df[[ycol]])
  hull <- df[hull_idx, , drop = FALSE]
  abs(polygon_area(hull[[xcol]], hull[[ycol]]))
}

ellipse_area <- function(df, xcol = "PC1", ycol = "PC2", level = 0.95) {
  df <- df[complete.cases(df[, c(xcol, ycol)]), , drop = FALSE]
  if (nrow(df) < 3) return(NA_real_)
  S <- stats::cov(df[, c(xcol, ycol), drop = FALSE])
  detS <- det(S)
  if (!is.finite(detS) || detS <= 0) return(NA_real_)
  q <- stats::qchisq(level, df = 2)
  pi * q * sqrt(detS)
}

compute_cluster_areas <- function(scores_df, group_col, level = 0.95) {
  group_col <- rlang::ensym(group_col)
  scores_df %>%
    dplyr::group_by(!!group_col) %>%
    dplyr::summarise(
      n = dplyr::n(),
      hull_area = convex_hull_area(dplyr::pick(PC1, PC2), xcol = "PC1", ycol = "PC2"),
      ellipse_area_95 = ellipse_area(dplyr::pick(PC1, PC2), xcol = "PC1", ycol = "PC2", level = level),
      .groups = "drop"
    )
}

############################################################
############################################################
## Helper 2d - publication-style silhouette plots
############################################################

make_silhouette_plot <- function(sil_df, scores_df, cluster_col, dataset_label, title) {
  cluster_col <- rlang::ensym(cluster_col)
  cluster_name <- rlang::as_name(cluster_col)

  sil_plot_df <- sil_df %>%
    dplyr::mutate(cluster = as.factor(.data[[cluster_name]])) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(dplyr::desc(sil_width), .by_group = TRUE) %>%
    dplyr::mutate(sample_order = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_id = dplyr::row_number())

  p_bar <- ggplot(sil_plot_df, aes(x = sample_id, y = sil_width, fill = cluster)) +
    geom_col(width = 1, color = NA) +
    geom_hline(yintercept = mean(sil_plot_df$sil_width, na.rm = TRUE),
               linetype = "dashed", linewidth = 0.6, color = "red3") +
    coord_flip() +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(
      title = title,
      subtitle = paste0(dataset_label, "; average silhouette width = ",
                        round(mean(sil_plot_df$sil_width, na.rm = TRUE), 3)),
      x = "Samples ordered within cluster",
      y = "Silhouette width",
      fill = cluster_name
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  p_scatter <- ggplot(scores_df, aes(x = PC1, y = PC2, color = as.factor(.data[[cluster_name]]))) +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_text(aes(label = Library), vjust = -0.8, size = 2.3, show.legend = FALSE) +
    labs(
      title = "PCA view of the same clusters",
      x = "PC1",
      y = "PC2",
      color = cluster_name
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  if (USE_PATCHWORK && requireNamespace("patchwork", quietly = TRUE)) {
    return(p_bar | p_scatter)
  }
  p_bar
}

## Helper 2c - PERMDISP diagnostics
############################################################

run_dispersion_test <- function(dist_obj, group_factor, permutations = 999) {
  grp <- as.factor(group_factor)
  
  bd <- vegan::betadisper(dist_obj, grp)
  pt <- vegan::permutest(bd, permutations = permutations)
  
  dist_to_centroid <- data.frame(
    sample = names(bd$distances),
    dist_to_centroid = unname(bd$distances),
    group = grp[names(bd$distances)],
    stringsAsFactors = FALSE
  )
  
  dist_summary <- dist_to_centroid %>%
    group_by(group) %>%
    summarise(
      n = n(),
      mean_dist = mean(dist_to_centroid),
      sd_dist = sd(dist_to_centroid),
      .groups = "drop"
    )
  
  list(
    betadisper = bd,
    permutest = pt,
    summary = data.frame(
      method = "betadisper_permutest",
      F = unname(pt$tab[1, "F"]),
      p = unname(pt$tab[1, "Pr(>F)"]),
      p_report = format_perm_p(unname(pt$tab[1, "Pr(>F)"]), permutations),
      permutations = permutations,
      stringsAsFactors = FALSE
    ),
    dist_to_centroid = dist_to_centroid,
    dist_summary = dist_summary
  )
}

############################################################
## Helper 3 - PERMANOVA interpretation and pairwise tests
############################################################

interpret_permanova <- function(perm_table, disp_summary, silhouette_mean = NA_real_, alpha = 0.05, context = "") {
  perm_p <- perm_table$`Pr(>F)`[1]
  perm_R2 <- perm_table$R2[1]
  perm_F <- perm_table$F[1]
  disp_p <- disp_summary$p[1]
  disp_F <- disp_summary$F[1]
  
  classification <- dplyr::case_when(
    !is.na(perm_p) && perm_p < alpha && !is.na(disp_p) && disp_p < alpha ~ "confounded_by_dispersion",
    !is.na(perm_p) && perm_p < alpha && (is.na(disp_p) || disp_p >= alpha) && !is.na(silhouette_mean) && silhouette_mean < 0.10 ~ "weak_centroid_separation",
    !is.na(perm_p) && perm_p < alpha && (is.na(disp_p) || disp_p >= alpha) ~ "centroid_separation_supported",
    (is.na(perm_p) || perm_p >= alpha) && !is.na(disp_p) && disp_p < alpha ~ "dispersion_difference_only",
    TRUE ~ "no_clear_multivariate_structure"
  )
  
  interpretation <- dplyr::case_when(
    classification == "confounded_by_dispersion" ~ "PERMANOVA is significant, but dispersion also differs among groups. Interpret centroid separation with caution because location and spread effects are confounded.",
    classification == "weak_centroid_separation" ~ "PERMANOVA supports centroid differences and dispersion is homogeneous, but silhouette support is weak. Treat the effect as statistically detectable but not strongly clustered.",
    classification == "centroid_separation_supported" ~ "PERMANOVA supports centroid differences and dispersion is homogeneous, so separation is more likely to reflect group location than spread.",
    classification == "dispersion_difference_only" ~ "Groups differ in dispersion without convincing PERMANOVA support for centroid separation.",
    TRUE ~ "Neither PERMANOVA nor dispersion diagnostics support strong multivariate separation."
  )
  
  data.frame(
    context = context,
    PERMANOVA_F = perm_F,
    PERMANOVA_R2 = perm_R2,
    PERMANOVA_p = perm_p,
    PERMANOVA_p_report = format_perm_p(perm_p, PERM_N),
    DISPERSION_F = disp_F,
    DISPERSION_p = disp_p,
    DISPERSION_p_report = ifelse(is.na(disp_p), NA_character_, format_perm_p(disp_p, PERM_N)),
    silhouette_mean = silhouette_mean,
    classification = classification,
    interpretation = interpretation,
    stringsAsFactors = FALSE
  )
}

pairwise_dispersion <- function(mat, group_factor, permutations = 999) {
  samples <- colnames(mat)
  groups <- as.factor(group_factor)
  names(groups) <- samples
  levs <- levels(groups)
  
  res_list <- list()
  
  for (i in seq_along(levs)) {
    for (j in seq_along(levs)) {
      if (j <= i) next
      gi <- levs[i]
      gj <- levs[j]
      keep_samp <- samples[groups %in% c(gi, gj)]
      grp_sub <- droplevels(groups[keep_samp])
      dist_sub <- dist(t(mat[, keep_samp, drop = FALSE]), method = "euclidean")
      disp_sub <- run_dispersion_test(dist_sub, grp_sub, permutations = permutations)
      
      res_list[[length(res_list) + 1]] <- data.frame(
        group1 = gi,
        group2 = gj,
        disp_F = disp_sub$summary$F[1],
        disp_p = disp_sub$summary$p[1],
        disp_p_report = disp_sub$summary$p_report[1],
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(res_list) == 0) return(NULL)
  out <- do.call(rbind, res_list)
  out$disp_p_adj_BH <- p.adjust(out$disp_p, method = "BH")
  out$disp_p_adj_report <- vapply(out$disp_p_adj_BH, function(x) as.character(signif(x, 3)), character(1))
  out
}

pairwise_permanova <- function(mat, group_factor, formula_var_name, permutations = 999, alpha = 0.05) {
  samples <- colnames(mat)
  groups <- as.factor(group_factor)
  names(groups) <- samples
  levs <- levels(groups)
  
  res_list <- list()
  
  for (i in seq_along(levs)) {
    for (j in seq_along(levs)) {
      if (j <= i) next
      gi <- levs[i]
      gj <- levs[j]
      
      keep_samp <- samples[groups %in% c(gi, gj)]
      mat_sub <- mat[, keep_samp, drop = FALSE]
      grp_sub <- droplevels(groups[keep_samp])
      
      meta_sub <- data.frame(
        sample = keep_samp,
        group = grp_sub,
        stringsAsFactors = FALSE
      )
      rownames(meta_sub) <- meta_sub$sample
      colnames(meta_sub)[2] <- formula_var_name
      
      set.seed(SEED)
      perm <- vegan::adonis2(
        t(mat_sub) ~ .,
        data = meta_sub[, formula_var_name, drop = FALSE],
        method = "euclidean",
        permutations = permutations
      )
      
      dist_sub <- dist(t(mat_sub), method = "euclidean")
      disp_sub <- run_dispersion_test(dist_sub, grp_sub, permutations = permutations)
      assessment <- interpret_permanova(
        perm_table = perm,
        disp_summary = disp_sub$summary,
        silhouette_mean = NA_real_,
        alpha = alpha,
        context = paste(gi, "vs", gj)
      )
      
      res_list[[length(res_list) + 1]] <- data.frame(
        group1 = gi,
        group2 = gj,
        R2 = perm$R2[1],
        F = perm$F[1],
        p = perm$`Pr(>F)`[1],
        p_report = format_perm_p(perm$`Pr(>F)`[1], permutations),
        sig = sig_code(perm$`Pr(>F)`[1]),
        disp_F = disp_sub$summary$F[1],
        disp_p = disp_sub$summary$p[1],
        disp_p_report = disp_sub$summary$p_report[1],
        classification = assessment$classification[1],
        interpretation = assessment$interpretation[1],
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- do.call(rbind, res_list)
  out$p_adj_BH <- p.adjust(out$p, method = "BH")
  out$p_adj_report <- vapply(out$p_adj_BH, function(x) as.character(signif(x, 3)), character(1))
  out$sig_adj <- vapply(out$p_adj_BH, sig_code, character(1))
  out$disp_p_adj_BH <- p.adjust(out$disp_p, method = "BH")
  out$disp_p_adj_report <- vapply(out$disp_p_adj_BH, function(x) as.character(signif(x, 3)), character(1))
  out
}

############################################################
## Helper 4 - verification utilities
############################################################

assert_no_na <- function(x, name) {
  if (any(is.na(x))) {
    stop(paste0("Verification failed: ", name, " contains NA values."))
  }
}

assert_same_order <- function(a, b, name_a, name_b) {
  if (!identical(a, b)) {
    stop(paste0("Verification failed: ", name_a, " and ", name_b, " are not identical in order."))
  }
}

warn_small_groups <- function(group_factor, context_label) {
  tab <- table(group_factor)
  if (any(tab < 3)) {
    warning(
      paste0(
        "Small group size detected in ", context_label, ": ",
        paste(names(tab), tab, sep = "=", collapse = ", "),
        ". PERMANOVA/pairwise tests may be unstable."
      )
    )
  }
  if (any(tab == 3)) {
    message(
      paste0(
        "Note: group size = 3 in ", context_label,
        ". Pairwise PERMANOVA will have coarse p-value resolution."
      )
    )
  }
}

############################################################
## Helper 5 - metadata builder
############################################################

build_sample_info <- function(library_names) {
  sample_info <- tibble(Library = library_names) %>%
    separate(Library, into = c("Appendage_raw", "Rest"), sep = "_", remove = FALSE) %>%
    mutate(
      Sex_Mating_status_raw = stringr::str_extract(Rest, "^[A-Za-z]+"),
      Replicate = stringr::str_extract(Rest, "[0-9]+"),
      Appendage = recode(
        Appendage_raw,
        "Ant" = "Antenna",
        "Leg" = "Tarsi",
        "P"   = "Maxillary palp",
        .default = Appendage_raw
      ),
      Sex = recode(
        Sex_Mating_status_raw,
        "MF" = "Female",
        "VF" = "Female",
        "Vm" = "Male",
        .default = Sex_Mating_status_raw
      ),
      Mating_status = recode(
        Sex_Mating_status_raw,
        "MF" = "Mated",
        "VF" = "Virgin",
        "Vm" = "Virgin",
        .default = Sex_Mating_status_raw
      ),
      Sex_Mating_status = paste(Sex, Mating_status, sep = " / "),
      Group = gsub("_\\d", "", Library)
    )
  
  sample_info$Appendage_raw <- factor(sample_info$Appendage_raw, levels = c("Ant", "P", "Leg"))
  sample_info$Sex_Mating_status_raw  <- factor(sample_info$Sex_Mating_status_raw,  levels = c("Vm", "VF", "MF"))
  sample_info$Sex <- factor(sample_info$Sex, levels = c("Male", "Female"))
  sample_info$Mating_status <- factor(sample_info$Mating_status, levels = c("Virgin", "Mated"))
  sample_info$Sex_Mating_status <- factor(sample_info$Sex_Mating_status, levels = c("Male / Virgin", "Female / Virgin", "Female / Mated"))
  
  sample_info
}

############################################################
## Helper 6 - variation partitioning
############################################################

compute_variation_partition <- function(log_counts_pca, sample_info) {
  X <- t(log_counts_pca)
  meta <- as.data.frame(sample_info)
  rownames(meta) <- meta$Library
  
  stopifnot(nrow(X) == nrow(meta))
  stopifnot(identical(rownames(X), rownames(meta)))
  
  meta$Appendage_raw <- factor(meta$Appendage_raw)
  meta$Sex_Mating_status_raw  <- factor(meta$Sex_Mating_status_raw)
  
  vp <- vegan::varpart(X, ~ Appendage_raw, ~ Sex_Mating_status_raw, data = meta)
  
  m_tissue_unique <- vegan::rda(X ~ Appendage_raw + Condition(Sex_Mating_status_raw), data = meta)
  m_state_unique  <- vegan::rda(X ~ Sex_Mating_status_raw + Condition(Appendage_raw), data = meta)
  
  tissue_unique <- vegan::RsquareAdj(m_tissue_unique)$adj.r.squared
  state_unique  <- vegan::RsquareAdj(m_state_unique)$adj.r.squared
  unexplained   <- 1 - (tissue_unique + state_unique)
  
  tissue_unique <- max(0, tissue_unique)
  state_unique  <- max(0, state_unique)
  unexplained   <- max(0, unexplained)
  
  fractions_df <- data.frame(
    Component = factor(
      c(
        "Appendage (unique effect)",
        "Sex / mating status (unique effect)",
        "Unexplained variance (other factors)"
      ),
      levels = c(
        "Appendage (unique effect)",
        "Sex / mating status (unique effect)",
        "Unexplained variance (other factors)"
      )
    ),
    Fraction = c(tissue_unique, state_unique, unexplained)
  )
  
  p_stack_h <- ggplot(
    fractions_df,
    aes(y = "Multivariate expression variance", x = Fraction, fill = Component)
  ) +
    geom_bar(stat = "identity", height = 0.4, color = "black") +
    geom_text(
      aes(label = sprintf("%.2f", Fraction),
          x = cumsum(Fraction) - Fraction / 2),
      y = "Multivariate expression variance",
      size = 3
    ) +
    scale_fill_manual(
      values = c(
        "Appendage (unique effect)"               = "grey20",
        "Sex / mating status (unique effect)"   = "grey60",
        "Unexplained variance (other factors)" = "white"
      )
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Fraction of total multivariate variance (Adj.R2-based)", y = "") +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      axis.text.y = element_text(size = 10)
    )
  
  list(
    vp = vp,
    tissue_unique = tissue_unique,
    state_unique = state_unique,
    unexplained = unexplained,
    fractions_df = fractions_df,
    plot = p_stack_h
  )
}

############################################################
## Helper 7 - correlation analysis
############################################################

run_correlation_analysis <- function(df_subset, libs, group_patterns, dataset_label) {
  counts_sel <- df_subset[, intersect(libs, colnames(df_subset)), drop = FALSE]
  counts_sel[is.na(counts_sel)] <- 0
  log_counts2 <- log2(counts_sel + 1)
  
  group_means <- lapply(group_patterns, function(pattern) {
    cols <- grep(paste0("^", pattern), colnames(log_counts2), value = TRUE)
    if (length(cols) == 0) return(NULL)
    rowMeans(log_counts2[, cols, drop = FALSE], na.rm = TRUE)
  })
  
  names(group_means) <- group_patterns
  group_means <- group_means[!sapply(group_means, is.null)]
  group_means_mat <- do.call(cbind, group_means)
  
  corr_mat <- cor(group_means_mat, method = "pearson", use = "pairwise.complete.obs")
  
  corr_df <- as.data.frame(corr_mat) %>%
    rownames_to_column("Row") %>%
    pivot_longer(cols = -Row, names_to = "Col", values_to = "corr") %>%
    mutate(
      Row = factor(Row, levels = colnames(corr_mat)),
      Col = factor(Col, levels = colnames(corr_mat))
    ) %>%
    filter(as.integer(Col) >= as.integer(Row)) %>%
    mutate(
      abs_corr = abs(corr),
      label = sprintf("%.2f", corr)
    )
  
  corr_df$Row <- factor(corr_df$Row, levels = rev(colnames(corr_mat)))
  corr_df$Col <- factor(corr_df$Col, levels = rev(colnames(corr_mat)))
  
  min_limit <- if (dataset_label == "ALL GENES") 0.86 else 0.50
  max_limit <- 1.00
  
  p_corr_bubble <- ggplot(corr_df, aes(x = Col, y = Row)) +
    geom_point(aes(size = abs_corr, fill = corr), shape = 21, color = "black") +
    scale_size(range = c(6, 12), name = "|Pearson r|") +
    scale_fill_gradientn(
      colours = c("white", "#ccffff", "#66ccff", "#3366ff", "#0033cc"),
      limits = c(min_limit, max_limit),
      oob = squish,
      name = "Pearson r"
    ) +
    scale_x_discrete(position = "top") +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      panel.grid = element_blank()
    ) +
    labs(
      x = "",
      y = "",
      title = paste0("Pearson correlation of mean expression - ", dataset_label)
    )
  
  p_corr_bubble_labeled <- p_corr_bubble +
    geom_text(aes(label = label), size = 5, color = "black") +
    labs(title = paste0("Pearson correlation of mean expression with labels - ", dataset_label))
  
  off_diag <- corr_df %>%
    filter(as.character(Row) != as.character(Col)) %>%
    pull(corr)
  
  summary_tbl <- tibble(
    n_pairs     = length(off_diag),
    min_corr    = min(off_diag, na.rm = TRUE),
    q1_corr     = quantile(off_diag, 0.25, na.rm = TRUE),
    median_corr = median(off_diag, na.rm = TRUE),
    mean_corr   = mean(off_diag, na.rm = TRUE),
    q3_corr     = quantile(off_diag, 0.75, na.rm = TRUE),
    max_corr    = max(off_diag, na.rm = TRUE)
  )
  
  top_pairs <- corr_df %>%
    filter(Row != Col) %>%
    arrange(desc(corr)) %>%
    head(30) %>%
    mutate(corr = round(corr, 3))
  
  bottom_pairs <- corr_df %>%
    filter(Row != Col) %>%
    arrange(corr) %>%
    head(30) %>%
    mutate(corr = round(corr, 3))
  
  list(
    group_means_mat = group_means_mat,
    corr_mat = corr_mat,
    corr_df = corr_df,
    summary_tbl = summary_tbl,
    top_pairs = top_pairs,
    bottom_pairs = bottom_pairs,
    plot = p_corr_bubble,
    plot_labeled = p_corr_bubble_labeled
  )
}

############################################################
## Helper 8 - appendage-level PCA
############################################################

run_appendage_pca <- function(appendage_raw_code, log_counts_mat, info_df, dataset_label) {
  libs_t <- info_df %>%
    filter(Appendage_raw == appendage_raw_code) %>%
    pull(Library)
  
  mat_t <- log_counts_mat[, libs_t, drop = FALSE]
  
  gene_sd_t <- apply(mat_t, 1, sd, na.rm = TRUE)
  zero_var_t <- sum(gene_sd_t == 0)
  
  mat_t_filt <- mat_t[gene_sd_t > 0, , drop = FALSE]
  
  pca_t <- prcomp(t(mat_t_filt), scale. = TRUE)
  
  scores_t <- data.frame(pca_t$x[, 1:2], stringsAsFactors = FALSE)
  colnames(scores_t) <- c("PC1", "PC2")
  scores_t$Library <- rownames(pca_t$x)
  
  scores_t <- scores_t %>% left_join(info_df, by = "Library")
  assert_no_na(scores_t$Sex_Mating_status_raw, paste0("scores_t$Sex_Mating_status_raw (", dataset_label, ", ", appendage_raw_code, ")"))
  
  var_expl_t <- pca_t$sdev^2 / sum(pca_t$sdev^2) * 100
  pc1_lab_t <- paste0("PC1 (", round(var_expl_t[1], 2), "%)")
  pc2_lab_t <- paste0("PC2 (", round(var_expl_t[2], 2), "%)")
  
  dist_t <- dist(t(mat_t_filt), method = "euclidean")
  cluster_t <- scores_t$Sex_Mating_status_raw
  
  warn_small_groups(cluster_t, paste0(dataset_label, " / sex/mating-status groups within appendage ", appendage_raw_code))
  
  sil_t <- silhouette(as.integer(cluster_t), dist_t)
  sil_t_df <- data.frame(
    Library = attr(dist_t, "Labels"),
    sil_width = sil_t[, "sil_width"],
    stringsAsFactors = FALSE
  ) %>%
    left_join(scores_t %>% select(Library, Sex_Mating_status_raw), by = "Library")
  
  assert_no_na(sil_t_df$Sex_Mating_status_raw, paste0("sil_t_df$Sex_Mating_status_raw (", dataset_label, ", ", appendage_raw_code, ")"))
  
  sil_t_summary <- sil_t_df %>%
    group_by(Sex_Mating_status_raw) %>%
    summarise(mean_sil_width = mean(sil_width), .groups = "drop")
  
  sil_t_mean <- mean(sil_t_df$sil_width)
  
  wb_t <- compute_within_between(dist_t, cluster_t)
  pairwise_state_meandist <- compute_pairwise_means(dist_t, cluster_t)
  
  meta_t <- as.data.frame(scores_t)
  rownames(meta_t) <- meta_t$Library
  
  set.seed(SEED)
  perm_state_t <- vegan::adonis2(
    t(mat_t_filt) ~ Sex_Mating_status_raw,
    data = meta_t,
    method = "euclidean",
    permutations = PERM_N
  )
  
  p_state_t <- perm_state_t$`Pr(>F)`[1]
  p_state_report <- format_perm_p(p_state_t, PERM_N)
  
  disp_t <- run_dispersion_test(dist_t, cluster_t, permutations = PERM_N)
  permanova_assessment_t <- interpret_permanova(
    perm_table = perm_state_t,
    disp_summary = disp_t$summary,
    silhouette_mean = sil_t_mean,
    alpha = 0.05,
    context = paste(dataset_label, unique(scores_t$Appendage), sep = " - ")
  )
  pairwise_perm_state <- pairwise_permanova(mat_t_filt, scores_t$Sex_Mating_status_raw, "Sex_Mating_status_raw", permutations = PERM_N)
  
  appendage_label <- unique(scores_t$Appendage)
  
  p_t <- ggplot(scores_t, aes(x = PC1, y = PC2)) +
    stat_ellipse(aes(group = Sex_Mating_status_raw), color = "grey30", linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(fill = Sex_Mating_status), shape = 21, size = 3, color = "black") +
    geom_text(aes(label = Library), vjust = -1, size = 2.5) +
    scale_fill_manual(
      values = c(
        "Male / Virgin"   = "lightblue",
        "Female / Virgin" = "lightpink",
        "Female / Mated"  = "purple"
      ),
      name = "Sex / mating status"
    ) +
    labs(
      title = paste("PCA -", appendage_label, "-", dataset_label),
      subtitle = paste0(
        "Clusters = sex / mating status; PERMANOVA p ", p_state_report,
        "; PERMDISP p ", disp_t$summary$p_report[1],
        "; mean silhouette = ", round(sil_t_mean, 3)
      ),
      x = pc1_lab_t,
      y = pc2_lab_t
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  centroids_state <- scores_t %>%
    group_by(Sex_Mating_status_raw, Sex_Mating_status) %>%
    summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2), .groups = "drop")
  
  radii_state <- scores_t %>%
    left_join(centroids_state, by = c("Sex_Mating_status_raw", "Sex_Mating_status")) %>%
    mutate(dist_to_centroid = sqrt((PC1 - mean_PC1)^2 + (PC2 - mean_PC2)^2)) %>%
    group_by(Sex_Mating_status_raw, Sex_Mating_status) %>%
    summarise(
      max_radius = max(dist_to_centroid),
      mean_radius = mean(dist_to_centroid),
      .groups = "drop"
    )
  
  areas_state <- compute_cluster_areas(scores_t, Sex_Mating_status_raw, level = ELLIPSE_LEVEL) %>%
    left_join(centroids_state %>% select(Sex_Mating_status_raw, Sex_Mating_status), by = "Sex_Mating_status_raw") %>%
    select(Sex_Mating_status_raw, Sex_Mating_status, n, hull_area, ellipse_area_95)
  
  centroid_state_mat <- as.matrix(centroids_state[, c("mean_PC1", "mean_PC2")])
  rownames(centroid_state_mat) <- centroids_state$Sex_Mating_status_raw
  d_state <- dist(centroid_state_mat)
  labs_state <- attr(d_state, "Labels")
  pairs_state <- t(combn(labs_state, 2))
  
  centroid_state_dist <- data.frame(
    Sex_Mating_status1 = pairs_state[, 1],
    Sex_Mating_status2 = pairs_state[, 2],
    centroid_distance = as.vector(d_state),
    stringsAsFactors = FALSE
  )
  
  list(
    plot = p_t,
    scores = scores_t,
    var_expl = var_expl_t,
    zero_var_removed = zero_var_t,
    sil_mean = sil_t_mean,
    sil_summary = sil_t_summary,
    sil_df = sil_t_df,
    silhouette_plot = make_silhouette_plot(sil_t_df, scores_t, Sex_Mating_status_raw, dataset_label, paste("Silhouette -", appendage_label)),
    wb = wb_t,
    pairwise_state_meandist = pairwise_state_meandist,
    perm = perm_state_t,
    perm_p_report = p_state_report,
    disp = disp_t,
    permanova_assessment = permanova_assessment_t,
    pairwise_perm = pairwise_perm_state,
    centroids = centroids_state,
    radii = radii_state,
    areas = areas_state,
    centroid_dist = centroid_state_dist
  )
}

############################################################
## Helper 9 - one unified analysis
############################################################

run_dataset_analysis <- function(df_input, dataset_label) {
  message("Running analysis for: ", dataset_label)
  
  missing <- setdiff(LIBS, colnames(df_input))
  if (length(missing) > 0) {
    warning("These expected library columns are missing: ", paste(missing, collapse = ", "))
  }
  
  df_selected_bsf <- df_input[, intersect(LIBS, colnames(df_input)), drop = FALSE]
  df_selected_bsf[is.na(df_selected_bsf)] <- 0
  
  log_counts <- log2(df_selected_bsf + 1)
  
  gene_sd <- apply(log_counts, 1, sd, na.rm = TRUE)
  sum_zero_var <- sum(gene_sd == 0)
  
  log_counts_pca <- log_counts[gene_sd > 0, , drop = FALSE]
  
  library_names <- colnames(log_counts_pca)
  sample_info <- build_sample_info(library_names)
  
  assert_no_na(sample_info$Appendage_raw, "Appendage")
  assert_no_na(sample_info$Sex_Mating_status_raw, "Sex / Mating status")
  assert_no_na(sample_info$Library, "Library")
  assert_same_order(colnames(log_counts_pca), sample_info$Library, "colnames(log_counts_pca)", "sample_info$Library")
  
  pca_global <- prcomp(t(log_counts_pca), scale. = TRUE)
  
  pca_global_df <- data.frame(pca_global$x[, 1:2], stringsAsFactors = FALSE)
  colnames(pca_global_df) <- c("PC1", "PC2")
  pca_global_df$Library <- rownames(pca_global$x)
  pca_global_df <- pca_global_df %>% left_join(sample_info, by = "Library")
  
  assert_no_na(pca_global_df$Appendage_raw, "pca_global_df$Appendage_raw")
  assert_no_na(pca_global_df$Sex_Mating_status_raw, "pca_global_df$Sex_Mating_status_raw")
  
  var_expl_g <- pca_global$sdev^2 / sum(pca_global$sdev^2) * 100
  pc1_lab_g <- paste0("PC1 (", round(var_expl_g[1], 2), "%)")
  pc2_lab_g <- paste0("PC2 (", round(var_expl_g[2], 2), "%)")
  
  global_centroids <- pca_global_df %>%
    group_by(Appendage_raw, Appendage) %>%
    summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2), .groups = "drop")
  
  global_radii <- pca_global_df %>%
    left_join(global_centroids, by = c("Appendage_raw", "Appendage")) %>%
    mutate(dist_to_centroid = sqrt((PC1 - mean_PC1)^2 + (PC2 - mean_PC2)^2)) %>%
    group_by(Appendage_raw, Appendage) %>%
    summarise(
      max_radius = max(dist_to_centroid),
      mean_radius = mean(dist_to_centroid),
      .groups = "drop"
    )
  
  global_areas <- compute_cluster_areas(pca_global_df, Appendage_raw, level = ELLIPSE_LEVEL) %>%
    left_join(global_centroids %>% select(Appendage_raw, Appendage), by = "Appendage_raw") %>%
    select(Appendage_raw, Appendage, n, hull_area, ellipse_area_95)
  
  centroid_mat <- as.matrix(global_centroids[, c("mean_PC1", "mean_PC2")])
  rownames(centroid_mat) <- global_centroids$Appendage_raw
  d_cent <- dist(centroid_mat)
  labs_cent <- attr(d_cent, "Labels")
  pairs_cent <- t(combn(labs_cent, 2))
  
  global_centroid_dist <- data.frame(
    Appendage1 = pairs_cent[, 1],
    Appendage2 = pairs_cent[, 2],
    centroid_distance = as.vector(d_cent),
    stringsAsFactors = FALSE
  )
  
  dist_global <- dist(t(log_counts_pca), method = "euclidean")
  cluster_global <- sample_info$Appendage_raw
  
  warn_small_groups(cluster_global, paste0(dataset_label, " / Global appendages"))
  
  sil_global <- silhouette(as.integer(cluster_global), dist_global)
  sil_global_df <- data.frame(
    Library = attr(dist_global, "Labels"),
    sil_width = sil_global[, "sil_width"],
    stringsAsFactors = FALSE
  ) %>%
    left_join(sample_info %>% select(Library, Appendage_raw), by = "Library")
  
  assert_no_na(sil_global_df$Appendage_raw, "sil_global_df$Appendage_raw")
  
  sil_global_summary <- sil_global_df %>%
    group_by(Appendage_raw) %>%
    summarise(mean_sil_width = mean(sil_width), .groups = "drop")
  
  sil_global_mean <- mean(sil_global_df$sil_width)
  
  wb_global <- compute_within_between(dist_global, cluster_global)
  pairwise_tissue_meandist <- compute_pairwise_means(dist_global, cluster_global)
  
  meta_for_perm <- as.data.frame(sample_info)
  rownames(meta_for_perm) <- meta_for_perm$Library
  
  set.seed(SEED)
  perm_tissue <- vegan::adonis2(
    t(log_counts_pca) ~ Appendage_raw,
    data = meta_for_perm,
    method = "euclidean",
    permutations = PERM_N
  )
  
  p_tissue <- perm_tissue$`Pr(>F)`[1]
  p_tissue_report <- format_perm_p(p_tissue, PERM_N)
  
  disp_global <- run_dispersion_test(dist_global, cluster_global, permutations = PERM_N)
  permanova_assessment_global <- interpret_permanova(
    perm_table = perm_tissue,
    disp_summary = disp_global$summary,
    silhouette_mean = sil_global_mean,
    alpha = 0.05,
    context = paste(dataset_label, "Global", sep = " - ")
  )
  pairwise_perm_tissue <- pairwise_permanova(log_counts_pca, sample_info$Appendage_raw, "Appendage_raw", permutations = PERM_N)
  
  p_global <- ggplot(pca_global_df, aes(x = PC1, y = PC2)) +
    stat_ellipse(aes(group = Appendage_raw), color = "grey30", linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(fill = Sex_Mating_status), shape = 21, size = 3, color = "black") +
    geom_text(aes(label = Library), vjust = -1, size = 2.5) +
    scale_fill_manual(
      values = c(
        "Male / Virgin"   = "lightblue",
        "Female / Virgin" = "lightpink",
        "Female / Mated"  = "purple"
      ),
      name = "Sex / mating status"
    ) +
    labs(
      title = paste("Global PCA -", dataset_label),
      subtitle = paste0("Clusters = appendages; PERMANOVA p ", p_tissue_report,
                        "; PERMDISP p ", disp_global$summary$p_report[1],
                        "; mean silhouette = ", round(sil_global_mean, 3)),
      x = pc1_lab_g,
      y = pc2_lab_g
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  res_ant <- run_appendage_pca("Ant", log_counts_pca, sample_info, dataset_label)
  res_p   <- run_appendage_pca("P",   log_counts_pca, sample_info, dataset_label)
  res_leg <- run_appendage_pca("Leg", log_counts_pca, sample_info, dataset_label)
  
  correlation <- run_correlation_analysis(df_input, LIBS, GROUP_PATTERNS, dataset_label)
  varpart_res <- compute_variation_partition(log_counts_pca, sample_info)
  
  list(
    dataset_label = dataset_label,
    n_genes_input = nrow(df_input),
    n_genes_used = nrow(log_counts_pca),
    zero_var_removed = sum_zero_var,
    sample_info = sample_info,
    log_counts_pca = log_counts_pca,
    pca_global = pca_global,
    pca_global_df = pca_global_df,
    var_expl_g = var_expl_g,
    sil_global_summary = sil_global_summary,
    sil_global_df = sil_global_df,
    sil_global_mean = sil_global_mean,
    silhouette_global_plot = make_silhouette_plot(sil_global_df, pca_global_df, Appendage_raw, dataset_label, "Silhouette - global appendage clusters"),
    wb_global = wb_global,
    pairwise_tissue_meandist = pairwise_tissue_meandist,
    global_centroids = global_centroids,
    global_radii = global_radii,
    global_areas = global_areas,
    global_centroid_dist = global_centroid_dist,
    perm_tissue = perm_tissue,
    p_tissue_report = p_tissue_report,
    disp_global = disp_global,
    permanova_assessment_global = permanova_assessment_global,
    pairwise_perm_tissue = pairwise_perm_tissue,
    plot_global = p_global,
    res_ant = res_ant,
    res_p = res_p,
    res_leg = res_leg,
    correlation = correlation,
    varpart = varpart_res
  )
}

############################################################
## Helper 10 - append one dataset to report
############################################################

append_dataset_report <- function(res) {
  add_report_header(paste0("SECTION: ", res$dataset_label))
  add_report_line("Input genes: ", res$n_genes_input)
  add_report_line("Genes used after zero-variance filtering: ", res$n_genes_used)
  add_report_line("Zero-variance genes removed: ", res$zero_var_removed)
  add_report_blank()
  
  add_report_object(paste0(res$dataset_label, " - Sample metadata"), res$sample_info)
  
  add_report_header(paste0(res$dataset_label, " - GLOBAL PCA SUMMARY"))
  add_report_line("PC1 variance percent: ", round(res$var_expl_g[1], 2))
  add_report_line("PC2 variance percent: ", round(res$var_expl_g[2], 2))
  add_report_line("Mean silhouette (global appendage PCA): ", round(res$sil_global_mean, 3))
  add_report_blank()
  
  add_report_object("Silhouette by appendage", res$sil_global_summary)
  add_report_object("Within vs between distances (global)", res$wb_global)
  add_report_object("Pairwise mean distances between appendages", res$pairwise_tissue_meandist)
  add_report_object("PCA centroids by appendage", res$global_centroids)
  add_report_object("PCA radii by appendage", res$global_radii)
  add_report_object("PCA cluster areas by appendage", res$global_areas)
  add_report_object("Pairwise centroid distances by appendage", res$global_centroid_dist)
  add_report_object("PERMANOVA appendage effect", res$perm_tissue)
  add_report_line("PERMANOVA p-value report: ", res$p_tissue_report)
  add_report_line("Permutation p-value floor: ", signif(perm_p_floor(PERM_N), 3))
  add_report_blank()
  add_report_object("Dispersion test for appendages", res$disp_global$summary)
  add_report_object("Integrated PERMANOVA interpretation - appendages", res$permanova_assessment_global)
  add_report_object("Dispersion summary for appendages", res$disp_global$dist_summary)
  add_report_object("Pairwise PERMANOVA between appendages", res$pairwise_perm_tissue)
  
  append_appendage <- function(res_obj, appendage_name) {
    add_report_header(paste0(res$dataset_label, " - ", toupper(appendage_name), " PCA SUMMARY"))
    add_report_line("PC1 variance percent: ", round(res_obj$var_expl[1], 2))
    add_report_line("PC2 variance percent: ", round(res_obj$var_expl[2], 2))
    add_report_line("Zero-variance genes removed within appendage: ", res_obj$zero_var_removed)
    add_report_line("Mean silhouette: ", round(res_obj$sil_mean, 3))
    add_report_blank()
    
    add_report_object("Silhouette by sex / mating status", res_obj$sil_summary)
    add_report_object("Within vs between distances", res_obj$wb)
    add_report_object("Pairwise mean distances between sex / mating-status groups", res_obj$pairwise_state_meandist)
    add_report_object("PCA centroids by sex / mating-status group", res_obj$centroids)
    add_report_object("PCA radii by sex / mating-status group", res_obj$radii)
    add_report_object("PCA cluster areas by sex / mating-status group", res_obj$areas)
    add_report_object("Pairwise centroid distances by sex / mating-status group", res_obj$centroid_dist)
    add_report_object("PERMANOVA sex / mating-status effect", res_obj$perm)
    add_report_line("PERMANOVA p-value report: ", res_obj$perm_p_report)
    add_report_line("Permutation p-value floor: ", signif(perm_p_floor(PERM_N), 3))
    add_report_blank()
    add_report_object("Dispersion test for sex / mating-status groups", res_obj$disp$summary)
    add_report_object("Integrated PERMANOVA interpretation - sex / mating-status groups", res_obj$permanova_assessment)
    add_report_object("Dispersion summary for sex / mating-status groups", res_obj$disp$dist_summary)
    add_report_object("Pairwise PERMANOVA between sex / mating-status groups", res_obj$pairwise_perm)
  }
  
  append_appendage(res$res_ant, "Antenna")
  append_appendage(res$res_p, "Maxillary palp")
  append_appendage(res$res_leg, "Tarsi")
  
  add_report_header(paste0(res$dataset_label, " - PEARSON CORRELATION ANALYSIS"))
  add_report_object("Group means matrix dimensions", dim(res$correlation$group_means_mat))
  add_report_object("Full Pearson correlation matrix", round(res$correlation$corr_mat, 3))
  add_report_object("Pearson correlation summary statistics", round(res$correlation$summary_tbl, 3))
  add_report_object("Top 30 highest correlation pairs", res$correlation$top_pairs %>% select(Row, Col, corr))
  add_report_object("Bottom 30 lowest correlation pairs", res$correlation$bottom_pairs %>% select(Row, Col, corr))
  
  add_report_header(paste0(res$dataset_label, " - VARIATION PARTITIONING"))
  add_report_object("varpart output", res$varpart$vp)
  add_report_line("Adj.R2 unique Appendage: ", signif(res$varpart$tissue_unique, 4))
  add_report_line("Adj.R2 unique Sex / mating status: ", signif(res$varpart$state_unique, 4))
  add_report_line("Unexplained         : ", signif(res$varpart$unexplained, 4))
  add_report_line(
    "Sum                 : ",
    signif(res$varpart$tissue_unique + res$varpart$state_unique + res$varpart$unexplained, 4)
  )
  add_report_blank()
}

############################################################

############################################################
## Helper 11 - compact, separated CSV outputs
############################################################

safe_dataset_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("(^_|_$)", "") %>%
    tolower()
}

adonis_first_row <- function(x, dataset, analysis_scope, effect, p_report) {
  tab <- as.data.frame(x)
  data.frame(
    dataset = dataset,
    analysis_scope = analysis_scope,
    effect = effect,
    df = tab$Df[1],
    sumsq = tab$SumOfSqs[1],
    r2 = tab$R2[1],
    f_statistic = tab$F[1],
    p_value = tab$`Pr(>F)`[1],
    p_report = p_report,
    stringsAsFactors = FALSE
  )
}


############################################################
## Helper 9 - 3-way variation partitioning
## Appendage | Sex | Mating
############################################################

compute_variation_partition_3way <- function(log_counts_pca, sample_info) {
  X <- t(log_counts_pca)
  meta <- as.data.frame(sample_info)
  rownames(meta) <- meta$Library
  stopifnot(identical(rownames(X), rownames(meta)))

  # Re-parameterise Sex_Mating_status_raw into two 1-df contrasts.
  meta$Sex <- as.integer(meta$Sex_Mating_status_raw == "Vm")
  meta$Mated <- as.integer(meta$Sex_Mating_status_raw == "MF")

  sex_mated_cor <- suppressWarnings(cor(meta$Sex, meta$Mated))

  vp <- vegan::varpart(X, ~ Appendage_raw, ~ Sex, ~ Mated, data = meta)

  m_tissue_unique <- vegan::rda(X ~ Appendage_raw + Condition(Sex) + Condition(Mated), data = meta)
  m_sex_unique <- vegan::rda(X ~ Sex + Condition(Appendage_raw) + Condition(Mated), data = meta)
  m_mated_unique <- vegan::rda(X ~ Mated + Condition(Appendage_raw) + Condition(Sex), data = meta)
  m_all <- vegan::rda(X ~ Appendage_raw + Sex + Mated, data = meta)

  tissue_unique <- max(0, vegan::RsquareAdj(m_tissue_unique)$adj.r.squared)
  sex_unique <- max(0, vegan::RsquareAdj(m_sex_unique)$adj.r.squared)
  mated_unique <- max(0, vegan::RsquareAdj(m_mated_unique)$adj.r.squared)
  explained_all <- max(0, vegan::RsquareAdj(m_all)$adj.r.squared)

  shared <- max(0, explained_all - (tissue_unique + sex_unique + mated_unique))
  unexplained <- max(0, 1 - (tissue_unique + sex_unique + mated_unique + shared))

  fractions_df <- tibble(
    Component = factor(
      c("Appendage unique", "Sex unique", "Mating unique", "Shared/confounded", "Unexplained"),
      levels = c("Appendage unique", "Sex unique", "Mating unique", "Shared/confounded", "Unexplained")
    ),
    Fraction = c(tissue_unique, sex_unique, mated_unique, shared, unexplained)
  )

  p_stack <- ggplot(fractions_df, aes(y = "Multivariate expression variance", x = Fraction, fill = Component)) +
    geom_col(width = 0.4, color = "black") +
    geom_text(aes(label = sprintf("%.3f", Fraction), x = cumsum(Fraction) - Fraction / 2),
              y = "Multivariate expression variance", size = 3) +
    scale_fill_manual(values = c(
      "Appendage unique" = "grey20",
      "Sex unique" = "#4f7fc2",
      "Mating unique" = "#c24f7f",
      "Shared/confounded" = "grey75",
      "Unexplained" = "white"
    )) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Fraction of total multivariate variance (Adj.R2-based)", y = "",
         title = "3-way variance partitioning: appendage | sex | mating") +
    theme_bw() +
    theme(legend.title = element_blank(), axis.text.y = element_text(size = 10))

  list(
    vp = vp,
    tissue_unique = tissue_unique,
    sex_unique = sex_unique,
    mated_unique = mated_unique,
    shared = shared,
    unexplained = unexplained,
    explained_all = explained_all,
    sex_mated_cor = sex_mated_cor,
    fractions_df = fractions_df,
    plot = p_stack
  )
}

make_variation_partition_csv <- function(res_all, res_chemo, vp3_all, vp3_chemo) {
  make_rows <- function(res, vp3) {
    dplyr::bind_rows(
      data.frame(
        dataset = res$dataset_label,
        model = "2-way: appendage + sex/mating-status group",
        component = c("Appendage unique", "Sex / mating-status unique", "Unexplained"),
        adj_r2_fraction = c(res$varpart$tissue_unique, res$varpart$state_unique, res$varpart$unexplained),
        sex_mated_correlation = NA_real_
      ),
      data.frame(
        dataset = res$dataset_label,
        model = "3-way: appendage + sex + mating",
        component = c("Appendage unique", "Sex unique", "Mating unique", "Shared/confounded", "Unexplained"),
        adj_r2_fraction = c(vp3$tissue_unique, vp3$sex_unique, vp3$mated_unique, vp3$shared, vp3$unexplained),
        sex_mated_correlation = c(vp3$sex_mated_cor, rep(NA_real_, 4))
      )
    )
  }
  dplyr::bind_rows(make_rows(res_all, vp3_all), make_rows(res_chemo, vp3_chemo))
}


append_3way_varpart_report <- function(vp3_all, vp3_chemo) {
  add_report_header("3-WAY VARIANCE PARTITIONING (Appendage | Sex | Mating)")
  add_report_blank()
  for (item in list(list(vp3 = vp3_all, label = "ALL GENES"), list(vp3 = vp3_chemo, label = "CHEMOSENSORY GENES"))) {
    add_report_header(paste0(item$label, " - 3-way partition"), "-")
    add_report_object("varpart output", item$vp3$vp)
    add_report_line("Adj.R2 Appendage unique  : ", signif(item$vp3$tissue_unique, 4))
    add_report_line("Adj.R2 Sex unique        : ", signif(item$vp3$sex_unique, 4))
    add_report_line("Adj.R2 Mating unique     : ", signif(item$vp3$mated_unique, 4))
    add_report_line("Shared/confounded        : ", signif(item$vp3$shared, 4))
    add_report_line("Unexplained              : ", signif(item$vp3$unexplained, 4))
    add_report_line("Sex-Mated correlation    : ", signif(item$vp3$sex_mated_cor, 4))
    add_report_blank()
  }
}

make_pca_summary_csv <- function(res) {
  dplyr::bind_rows(
    data.frame(dataset = res$dataset_label, analysis_scope = "global_appendage_pca", group_variable = "Appendage", pc1_variance = res$var_expl_g[1], pc2_variance = res$var_expl_g[2], mean_silhouette = res$sil_global_mean, n_genes_input = res$n_genes_input, n_genes_used = res$n_genes_used, zero_var_removed = res$zero_var_removed),
    data.frame(dataset = res$dataset_label, analysis_scope = "antenna_sex_mating_pca", group_variable = "Sex / Mating status", pc1_variance = res$res_ant$var_expl[1], pc2_variance = res$res_ant$var_expl[2], mean_silhouette = res$res_ant$sil_mean, n_genes_input = res$n_genes_input, n_genes_used = res$n_genes_used, zero_var_removed = res$res_ant$zero_var_removed),
    data.frame(dataset = res$dataset_label, analysis_scope = "palp_sex_mating_pca", group_variable = "Sex / Mating status", pc1_variance = res$res_p$var_expl[1], pc2_variance = res$res_p$var_expl[2], mean_silhouette = res$res_p$sil_mean, n_genes_input = res$n_genes_input, n_genes_used = res$n_genes_used, zero_var_removed = res$res_p$zero_var_removed),
    data.frame(dataset = res$dataset_label, analysis_scope = "tarsi_sex_mating_pca", group_variable = "Sex / Mating status", pc1_variance = res$res_leg$var_expl[1], pc2_variance = res$res_leg$var_expl[2], mean_silhouette = res$res_leg$sil_mean, n_genes_input = res$n_genes_input, n_genes_used = res$n_genes_used, zero_var_removed = res$res_leg$zero_var_removed)
  )
}

make_group_geometry_csv <- function(res) {
  global <- res$global_centroids %>%
    left_join(res$global_radii, by = c("Appendage_raw", "Appendage")) %>%
    left_join(res$global_areas, by = c("Appendage_raw", "Appendage")) %>%
    left_join(res$sil_global_summary, by = "Appendage_raw") %>%
    mutate(dataset = res$dataset_label, analysis_scope = "global_appendage_pca", group_variable = "Appendage") %>%
    rename(group_raw = Appendage_raw, group_label = Appendage)

  appendage_geom <- function(x, scope) {
    x$centroids %>%
      left_join(x$radii, by = c("Sex_Mating_status_raw", "Sex_Mating_status")) %>%
      left_join(x$areas, by = c("Sex_Mating_status_raw", "Sex_Mating_status")) %>%
      left_join(x$sil_summary, by = "Sex_Mating_status_raw") %>%
      mutate(dataset = res$dataset_label, analysis_scope = scope, group_variable = "Sex / Mating status") %>%
      rename(group_raw = Sex_Mating_status_raw, group_label = Sex_Mating_status)
  }

  dplyr::bind_rows(
    global,
    appendage_geom(res$res_ant, "antenna_sex_mating_pca"),
    appendage_geom(res$res_p, "palp_sex_mating_pca"),
    appendage_geom(res$res_leg, "tarsi_sex_mating_pca")
  )
}

make_permanova_csv <- function(res) {
  dplyr::bind_rows(
    adonis_first_row(res$perm_tissue, res$dataset_label, "global_appendage_pca", "Appendage", res$p_tissue_report),
    adonis_first_row(res$res_ant$perm, res$dataset_label, "antenna_sex_mating_pca", "Sex / Mating status", res$res_ant$perm_p_report),
    adonis_first_row(res$res_p$perm, res$dataset_label, "palp_sex_mating_pca", "Sex / Mating status", res$res_p$perm_p_report),
    adonis_first_row(res$res_leg$perm, res$dataset_label, "tarsi_sex_mating_pca", "Sex / Mating status", res$res_leg$perm_p_report)
  )
}

make_pairwise_permanova_csv <- function(res) {
  add_meta <- function(df, scope, group_variable) {
    df %>%
      mutate(
        dataset = res$dataset_label,
        analysis_scope = scope,
        group_variable = group_variable,
        test = "pairwise PERMANOVA",
        distance = "Euclidean distance on log2(x + 1) expression matrix",
        p_adjustment = "Benjamini-Hochberg within each analysis scope"
      ) %>%
      select(
        dataset, analysis_scope, group_variable, test, distance, group1, group2,
        R2, F, p, p_report, p_adj_BH, p_adj_report, sig_adj,
        disp_F, disp_p, disp_p_report, disp_p_adj_BH, disp_p_adj_report,
        classification, interpretation, p_adjustment, everything()
      )
  }

  dplyr::bind_rows(
    add_meta(res$pairwise_perm_tissue, "global_appendage_pca", "Appendage"),
    add_meta(res$res_ant$pairwise_perm, "antenna_sex_mating_pca", "Sex / Mating status within antenna"),
    add_meta(res$res_p$pairwise_perm, "palp_sex_mating_pca", "Sex / Mating status within maxillary palp"),
    add_meta(res$res_leg$pairwise_perm, "tarsi_sex_mating_pca", "Sex / Mating status within tarsi")
  )
}

make_pairwise_distance_csv <- function(res) {
  add_meta <- function(df, scope, group_variable) {
    df %>%
      mutate(
        dataset = res$dataset_label,
        analysis_scope = scope,
        group_variable = group_variable,
        metric = "mean pairwise Euclidean distance",
        value_definition = "Average distance between all samples in group1 and all samples in group2, calculated on log2(x + 1) expression values"
      ) %>%
      rename(mean_pairwise_distance = mean_dist) %>%
      select(dataset, analysis_scope, group_variable, metric, group1, group2, mean_pairwise_distance, value_definition)
  }

  dplyr::bind_rows(
    add_meta(res$pairwise_tissue_meandist, "global_appendage_pca", "Appendage"),
    add_meta(res$res_ant$pairwise_state_meandist, "antenna_sex_mating_pca", "Sex / Mating status within antenna"),
    add_meta(res$res_p$pairwise_state_meandist, "palp_sex_mating_pca", "Sex / Mating status within maxillary palp"),
    add_meta(res$res_leg$pairwise_state_meandist, "tarsi_sex_mating_pca", "Sex / Mating status within tarsi")
  )
}

make_csv_index <- function() {
  data.frame(
    file = c(
      "01_pca_overview.csv",
      "02_pca_group_geometry_and_silhouette.csv",
      "03_permanova_summary.csv",
      "04a_pairwise_permanova_tests.csv",
      "04b_pairwise_mean_distances.csv",
      "05_pearson_correlation_long.csv",
      "06_pearson_correlation_summary.csv",
      "07_variation_partitioning.csv"
    ),
    contains = c(
      "Variance explained by PC1/PC2 and mean silhouette for each PCA",
      "PCA centroid coordinates, centroid radius, convex-hull area, and silhouette summaries per group",
      "Global PERMANOVA summaries for appendage effects and within-appendage sex/mating-status effects",
      "Pairwise PERMANOVA results only: R2, F, p, BH-adjusted p, dispersion diagnostics, and interpretation",
      "Pairwise mean Euclidean distances only; descriptive distance metric, not a significance test",
      "Full long-format Pearson correlation matrix between appendage x sex/mating-status group means",
      "Compact Pearson correlation summaries",
      "2-way and 3-way adjusted-R2 variation partitioning: appendage, sex, mating status, shared/confounded, unexplained"
    ),
    stringsAsFactors = FALSE
  )
}

write_compact_csv_outputs <- function(res_all, res_chemo, vp3_all, vp3_chemo, out_dir = OUTPUT_DIR) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(make_csv_index(), file.path(out_dir, "00_csv_index_read_me.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(make_pca_summary_csv(res_all), make_pca_summary_csv(res_chemo)), file.path(out_dir, "01_pca_overview.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(make_group_geometry_csv(res_all), make_group_geometry_csv(res_chemo)), file.path(out_dir, "02_pca_group_geometry_and_silhouette.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(make_permanova_csv(res_all), make_permanova_csv(res_chemo)), file.path(out_dir, "03_permanova_summary.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(make_pairwise_permanova_csv(res_all), make_pairwise_permanova_csv(res_chemo)), file.path(out_dir, "04a_pairwise_permanova_tests.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(make_pairwise_distance_csv(res_all), make_pairwise_distance_csv(res_chemo)), file.path(out_dir, "04b_pairwise_mean_distances.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(res_all$correlation$corr_df %>% mutate(dataset = res_all$dataset_label), res_chemo$correlation$corr_df %>% mutate(dataset = res_chemo$dataset_label)), file.path(out_dir, "05_pearson_correlation_long.csv"), row.names = FALSE)
  write.csv(dplyr::bind_rows(res_all$correlation$summary_tbl %>% mutate(dataset = res_all$dataset_label), res_chemo$correlation$summary_tbl %>% mutate(dataset = res_chemo$dataset_label)), file.path(out_dir, "06_pearson_correlation_summary.csv"), row.names = FALSE)
  write.csv(make_variation_partition_csv(res_all, res_chemo, vp3_all, vp3_chemo), file.path(out_dir, "07_variation_partitioning.csv"), row.names = FALSE)
  invisible(out_dir)
}

## Main
############################################################

df_all <- read.csv(INPUT_CSV, row.names = 1, check.names = FALSE)

if (!("Name" %in% colnames(df_all))) {
  stop("Column 'Name' not found in input CSV. Please ensure there is a column named exactly 'Name'.")
}

df_chemo <- df_all[!is.na(df_all$Name) & trimws(df_all$Name) != "", , drop = FALSE]

cat("All genes table dimensions:\n")
print(dim(df_all))
cat("\nChemosensory subset dimensions:\n")
print(dim(df_chemo))

res_all   <- run_dataset_analysis(df_all,   "ALL GENES")
res_chemo <- run_dataset_analysis(df_chemo, "CHEMOSENSORY GENES (Name non-empty)")

cat("\n-- Running 3-way variance partitioning --\n")
vp3_all <- compute_variation_partition_3way(res_all$log_counts_pca, res_all$sample_info)
vp3_chemo <- compute_variation_partition_3way(res_chemo$log_counts_pca, res_chemo$sample_info)
print(vp3_all$fractions_df)
print(vp3_chemo$fractions_df)

write_compact_csv_outputs(res_all, res_chemo, vp3_all, vp3_chemo, OUTPUT_DIR)
cat("\nCSV tables written to: ", OUTPUT_DIR, "\n", sep = "")

append_dataset_report(res_all)
append_dataset_report(res_chemo)
append_3way_varpart_report(vp3_all, vp3_chemo)
append_exported_csvs_to_report(OUTPUT_DIR)
write_report(file.path(OUTPUT_DIR, "Figure1_full_text_report.txt"))
cat("Text report written to: ", file.path(OUTPUT_DIR, "Figure1_full_text_report.txt"), "\n", sep = "")

############################################################
## Show plots
############################################################

print(res_all$plot_global)
print(res_all$res_ant$plot)
print(res_all$res_p$plot)
print(res_all$res_leg$plot)
print(res_all$silhouette_global_plot)
print(res_all$res_ant$silhouette_plot)
print(res_all$res_p$silhouette_plot)
print(res_all$res_leg$silhouette_plot)
print(res_all$correlation$plot)
print(res_all$correlation$plot_labeled)
print(res_all$varpart$plot)
print(vp3_all$plot)

print(res_chemo$plot_global)
print(res_chemo$res_ant$plot)
print(res_chemo$res_p$plot)
print(res_chemo$res_leg$plot)
print(res_chemo$silhouette_global_plot)
print(res_chemo$res_ant$silhouette_plot)
print(res_chemo$res_p$silhouette_plot)
print(res_chemo$res_leg$silhouette_plot)
print(res_chemo$correlation$plot)
print(res_chemo$correlation$plot_labeled)
print(res_chemo$varpart$plot)
print(vp3_chemo$plot)

if (USE_PATCHWORK && requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  
  p_allgenes_panel <- (res_all$plot_global | res_all$res_ant$plot) /
    (res_all$res_p$plot | res_all$res_leg$plot)
  
  p_chemo_panel <- (res_chemo$plot_global | res_chemo$res_ant$plot) /
    (res_chemo$res_p$plot | res_chemo$res_leg$plot)
  
  print(p_allgenes_panel)
  print(p_chemo_panel)
}


# ── Session info (run to capture package versions for Methods reporting) ─────
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
