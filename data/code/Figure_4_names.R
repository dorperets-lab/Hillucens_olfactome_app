################################################################################
# Shared gene-name shortening for Figure 4.
#
# Sourced by Figure_4B_lollipop.R (panel B) and Figure_4.R (panel D scatter) so
# the two panels of one figure never disagree about what a gene is called.
# Depends on nothing but base R; the caller supplies the raw NCBI description.
################################################################################
# ── Name shortening ─────────────────────────────────────────────────────────
# NCBI descriptions run to 80 characters. Truncating them turned the strongest
# hit in the panel into "neuropeptide SIFamide rece..." and split single genes
# across two rows ("chitin deacetylase 1" and "chitin deacetylase 1, transcript
# variant X1"). Shorten in three stages, all of them declared:
#   (a) drop NCBI boilerplate - transcript-variant suffixes, "-like",
#       "probable/putative/predicted", "[UDP-forming]"
#   (b) an explicit abbreviation table for the names that recur here
#   (c) generic word abbreviations for anything the table misses
# The full name of every abbreviation is written to Figure_4B_name_key.csv so
# the caption can point at it.
strip_boiler <- function(x){
  x <- sub("\\s*\\[Source:.*$", "", x)
  x <- sub("\\s*\\[UDP-forming\\]", "", x)
  x <- sub(",?\\s*transcript variant\\s+[A-Za-z0-9]+\\s*$", "", x, ignore.case=TRUE)
  x <- sub("-?like\\s*$", "", x, ignore.case=TRUE)
  x <- sub("^(probable|putative|predicted)\\s+", "", x, ignore.case=TRUE)
  trimws(sub("[,\\-]\\s*$", "", x))
}
ABBR_EXACT <- c(
  "neuropeptide sifamide receptor"                            = "SIFamide receptor",
  "a-kinase anchor protein 14"                                = "Akap14",
  "alpha,alpha-trehalose-phosphate synthase"                  = "Tps1",
  "serine/threonine-protein kinase sik2"                      = "Sik2 kinase",
  "mediator of rna polymerase ii transcription subunit 12"    = "Med12",
  "cytochrome p450 4d2"                                       = "Cyp4d2",
  "female-specific protein transformer"                       = "Transformer (tra)",
  "atp-dependent dna helicase hfm1"                           = "Hfm1 helicase",
  "tyrosine-protein phosphatase non-receptor type 23"         = "Ptpn23",
  "gata zinc finger domain-containing protein 14"             = "Gatad14",
  "cdc42 effector protein 1"                                  = "Cdc42EP1",
  "protein decapentaplegic"                                   = "Dpp",
  "protein spaetzle"                                          = "Spaetzle",
  "protein obstructor-e"                                      = "Obst-E",
  "chitin deacetylase 1"                                      = "Cda1",
  "adult cuticle protein 1"                                   = "Acp1",
  "farnesol dehydrogenase"                                    = "Farnesol DH",
  "enhancer of split mbeta protein"                           = "E(spl)mbeta",
  "chemosensory protein a 87a"                                = "CheA87a",
  "vitellogenin-a1"                                           = "Vitellogenin-A1",
  "glycine-rich cell wall structural protein"                  = "Gly-rich cell wall prot"
)
# (c) generic substitutions, longest pattern first
ABBR_GEN <- c(
  "endocuticle structural glycoprotein" = "Ecg",
  "endocuticle structural protein"      = "Ecp",
  "larval cuticle protein"              = "Lcp",
  "pupal cuticle protein"               = "Pcp",
  "adult cuticle protein"               = "Acp",
  "cuticle protein"                     = "Cpr",
  "chitin deacetylase"                  = "Cda",
  "cytochrome p450"                     = "Cyp",
  "serine/threonine-protein kinase"     = "S/T kinase",
  "tyrosine-protein phosphatase"        = "Tyr phosphatase",
  "dehydrogenase"                       = "DH",
  "^protein "                           = ""
)
shorten <- function(x){
  b <- strip_boiler(x)
  k <- tolower(b)
  out <- unname(ABBR_EXACT[k])
  hit <- !is.na(out)
  for (i in which(!hit)) {
    s <- b[i]                      # one element at a time: sub() is vectorised
    for (j in seq_along(ABBR_GEN)) s <- sub(names(ABBR_GEN)[j], ABBR_GEN[j], s, ignore.case=TRUE)
    s <- trimws(s)
    # "cuticle protein-like" collapses to a bare "Cpr", which names nothing -
    # keep the spelt-out word when no identifier follows the abbreviation
    if (s %in% ABBR_GEN) s <- b[i]
    s <- sub("^(Cpr|Lcp|Pcp|Acp|Ecg|Ecp|Cda|Cyp) +([0-9])", "\\1\\2", s)
    out[i] <- if (nzchar(s)) paste0(toupper(substring(s,1,1)), substring(s,2)) else ""
  }
  out
}
