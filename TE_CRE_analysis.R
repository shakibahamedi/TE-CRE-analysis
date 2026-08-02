#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

# ============================================================
# CRE / transposable-element (TE) overlap analysis
#
# Takes ONLY two inputs and writes all outputs automatically:
#   1) a CRE .tsv file   (columns: CREID, promoter_type, ...)
#   2) a UCSC RepeatMasker rmsk.txt / rmsk.txt.gz file
#
# Usage:
#   Rscript TE_CRE_analysis.R  <CRE_file.tsv>  <rmsk_file.txt.gz>
#
# Example:
#   Rscript TE_CRE_analysis.R THP1_Series.CRE.info.p.e.se.tsv rmsk.txt.gz
#
# The sample name is taken automatically from the CRE file name,
# and results are written to a new folder named "<sample>_CRE_TE_output"
# next to where you run the script. Nothing else needs to be set --
# run the same command for any similar dataset (THP1, Neuron, HEK, ...).
#
# CRE groups analyzed together:
#   promoter-like, enhancer-like, CTCF-alone, unclassed
#   (any other promoter_type label is excluded and saved separately)
#
# Repeat scope:
#   DNA, LINE, SINE, LTR, RC, Retroposon (including "?" uncertain labels)
#
# Coordinates: CRE and rmsk inputs are assumed 0-based, half-open (BED/UCSC).
# CREID format: chr_start_end_strand (parsed from the right, so chromosome
# names that contain underscores, e.g. chr1_KI270706v1_random, still work).
# ============================================================

STRAND_MODE      <- "same"   # same | opposite | both
MIN_OVERLAP_BP   <- 25L
MIN_FRAC_CRE     <- 0.10
TE_CLASS_PATTERN <- "^(DNA|LINE|SINE|LTR|RC|Retroposon)\\??$"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript TE_CRE_analysis.R <CRE_file.tsv> <rmsk_file.txt.gz>",
    call. = FALSE
  )
}
cre_file  <- args[[1L]]
rmsk_file <- args[[2L]]

for (input_file in c(cre_file, rmsk_file)) {
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)
}

# Derive a clean sample name from the CRE file name,
# e.g. "THP1_Series.CRE.info.p.e.se.tsv" -> "THP1"
sample_name <- basename(cre_file)
sample_name <- sub("\\.(tsv|txt|csv)(\\.gz)?$", "", sample_name, ignore.case = TRUE)
sample_name <- sub("[._-]?(Series)?[._-]?CRE.*$", "", sample_name, ignore.case = TRUE)
sample_name <- sub("[._-]+$", "", sample_name)
if (sample_name == "" || is.na(sample_name)) sample_name <- "sample"

out_dir <- paste0(sample_name, "_CRE_TE_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- function(suffix) file.path(out_dir, paste0(sample_name, "_", suffix))

message("============================================================")
message("[", sample_name, "] Starting CRE/TE overlap analysis")
message("[", sample_name, "] CRE file:  ", cre_file)
message("[", sample_name, "] rmsk file: ", rmsk_file)
message("[", sample_name, "] Output:    ", normalizePath(out_dir, mustWork = FALSE))
message("============================================================")

# ------------------------------------------------------------
# 1. Read and validate CRE data
# ------------------------------------------------------------
message("[", sample_name, "] Reading CRE data...")
cre <- fread(cre_file, header = TRUE, sep = "\t", na.strings = c("", "NA"))

required_cre_cols <- c("CREID", "promoter_type")
missing_cre_cols <- setdiff(required_cre_cols, names(cre))
if (length(missing_cre_cols) > 0L) {
  stop("Missing required CRE columns: ", paste(missing_cre_cols, collapse = ", "))
}

cre[, CREID := trimws(as.character(CREID))]
if (anyNA(cre$CREID) || any(cre$CREID == "")) stop("CREID contains missing or empty values.")
if (anyDuplicated(cre$CREID)) stop("CREID contains duplicated values.")

# CREID format: chr_start_end_strand, parsed from the right so
# chromosome labels containing underscores are supported.
cre[, cre_strand := sub("^.*_([^_]*)$", "\\1", CREID)]
cre[, id_without_strand := sub("_[^_]*$", "", CREID)]
cre[, end := suppressWarnings(as.integer(sub("^.*_([^_]*)$", "\\1", id_without_strand)))]
cre[, id_without_end := sub("_[^_]*$", "", id_without_strand)]
cre[, start := suppressWarnings(as.integer(sub("^.*_([^_]*)$", "\\1", id_without_end)))]
cre[, chr := sub("_[^_]*$", "", id_without_end)]
cre[, c("id_without_strand", "id_without_end") := NULL]

invalid_cre <- cre[
  is.na(chr) | chr == "" | is.na(start) | is.na(end) |
    end <= start | !cre_strand %chin% c("+", "-", ".")
]
if (nrow(invalid_cre) > 0L) {
  fwrite(invalid_cre, output_file("CRE_invalid_coordinates.tsv"), sep = "\t")
  stop(
    nrow(invalid_cre),
    " CRE rows have invalid coordinates or strand. See ",
    output_file("CRE_invalid_coordinates.tsv")
  )
}

# ------------------------------------------------------------
# 2. Standardize CRE classes
# ------------------------------------------------------------
cre[, promoter_type_clean := tolower(trimws(as.character(promoter_type)))]
cre[, CRE_type := fcase(
  promoter_type_clean == "promoter-like", "promoter",
  promoter_type_clean == "enhancer-like", "enhancer",
  promoter_type_clean == "ctcf-alone",    "CTCF-alone",
  promoter_type_clean == "unclassed",     "unclassed",
  default = "other"
)]

cre_other <- cre[CRE_type == "other"]
if (nrow(cre_other) > 0L) {
  fwrite(cre_other, output_file("other_CRE_labels_excluded.tsv.gz"), sep = "\t")
}

cre <- cre[CRE_type %chin% c("promoter", "enhancer", "CTCF-alone", "unclassed")]

cre_unknown_strand <- cre[!cre_strand %chin% c("+", "-")]
if (nrow(cre_unknown_strand) > 0L) {
  fwrite(cre_unknown_strand, output_file("CRE_unknown_strand_excluded.tsv.gz"), sep = "\t")
}

cre <- cre[cre_strand %chin% c("+", "-")]
if (nrow(cre) == 0L) stop("No eligible CREs with known strand remain.")
cre[, cre_index := .I]

message("[", sample_name, "] CRE counts used in analysis:")
print(cre[, .N, by = .(CRE_type, cre_strand)][order(CRE_type, cre_strand)])

# ------------------------------------------------------------
# 3. Read and filter UCSC RepeatMasker data
# ------------------------------------------------------------
message("[", sample_name, "] Reading RepeatMasker data and retaining TE classes...")

# UCSC rmsk fields:
# 6 genoName, 7 genoStart, 8 genoEnd, 10 strand,
# 11 repName, 12 repClass, 13 repFamily, 17 id
rmsk <- fread(
  rmsk_file,
  header = FALSE,
  sep = "\t",
  select = c(6L, 7L, 8L, 10L, 11L, 12L, 13L, 17L),
  col.names = c(
    "chr", "start", "end", "rmsk_strand",
    "repName", "repClass", "repFamily", "rmsk_id"
  ),
  na.strings = c("", "NA")
)

rmsk[, `:=`(
  chr         = trimws(as.character(chr)),
  start       = suppressWarnings(as.integer(start)),
  end         = suppressWarnings(as.integer(end)),
  rmsk_strand = trimws(as.character(rmsk_strand)),
  repName     = trimws(as.character(repName)),
  repClass    = trimws(as.character(repClass)),
  repFamily   = trimws(as.character(repFamily))
)]

rmsk <- rmsk[
  chr %chin% unique(cre$chr) &
    !is.na(start) & !is.na(end) & end > start &
    rmsk_strand %chin% c("+", "-") &
    grepl(TE_CLASS_PATTERN, repClass, ignore.case = TRUE)
]

if (nrow(rmsk) == 0L) stop("No valid transposable-element records remain.")
rmsk[, te_index := .I]

message("[", sample_name, "] TE rows used: ", format(nrow(rmsk), big.mark = ","))
print(rmsk[, .N, by = repClass][order(-N)])

# ------------------------------------------------------------
# 4. Build GRanges objects
# ------------------------------------------------------------
optional_cre_col <- function(column_name) {
  if (column_name %in% names(cre)) cre[[column_name]] else rep(NA_character_, nrow(cre))
}

cre_gr <- GRanges(
  seqnames = cre$chr,
  ranges = IRanges(start = cre$start + 1L, end = cre$end),
  strand = cre$cre_strand,
  CRE_ID = cre$CREID,
  CREName = optional_cre_col("CREName"),
  CRE_type = cre$CRE_type,
  promoter_type = cre$promoter_type,
  proximity = optional_cre_col("proximity"),
  typeStr = optional_cre_col("typeStr"),
  gene_promoter = optional_cre_col("gene_promoter"),
  geneNameStr = optional_cre_col("geneNameStr"),
  cre_index = cre$cre_index
)

rmsk_gr <- GRanges(
  seqnames = rmsk$chr,
  ranges = IRanges(start = rmsk$start + 1L, end = rmsk$end),
  strand = rmsk$rmsk_strand,
  rmsk_id = rmsk$rmsk_id,
  repName = rmsk$repName,
  repClass = rmsk$repClass,
  repFamily = rmsk$repFamily,
  te_index = rmsk$te_index
)

# ------------------------------------------------------------
# 5. Find overlaps and assign strand relationship
# ------------------------------------------------------------
message("[", sample_name, "] Finding genomic overlaps...")
hits <- findOverlaps(cre_gr, rmsk_gr, ignore.strand = TRUE)
if (length(hits) == 0L) stop("No positional CRE/TE overlaps were found.")

q <- queryHits(hits)
s <- subjectHits(hits)
strand_relation <- fifelse(
  as.character(strand(cre_gr[q])) == as.character(strand(rmsk_gr[s])),
  "same",
  "opposite"
)

keep_hits <- switch(
  STRAND_MODE,
  same = strand_relation == "same",
  opposite = strand_relation == "opposite",
  both = rep(TRUE, length(strand_relation))
)

q <- q[keep_hits]
s <- s[keep_hits]
strand_relation <- strand_relation[keep_hits]
if (length(q) == 0L) stop("No overlaps remain for strand_mode = ", STRAND_MODE)

# ------------------------------------------------------------
# 6. Construct overlap table
# ------------------------------------------------------------
overlap_start <- pmax(start(cre_gr[q]), start(rmsk_gr[s]))
overlap_end <- pmin(end(cre_gr[q]), end(rmsk_gr[s]))
overlap_len <- overlap_end - overlap_start + 1L

cre_meta <- mcols(cre_gr)
rmsk_meta <- mcols(rmsk_gr)

overlap_all <- data.table(
  CRE_chr = as.character(seqnames(cre_gr[q])),
  CRE_start = start(cre_gr[q]) - 1L,
  CRE_end = end(cre_gr[q]),
  CRE_strand = as.character(strand(cre_gr[q])),
  CRE_ID = cre_meta$CRE_ID[q],
  CREName = cre_meta$CREName[q],
  CRE_type = cre_meta$CRE_type[q],
  promoter_type = cre_meta$promoter_type[q],
  proximity = cre_meta$proximity[q],
  typeStr = cre_meta$typeStr[q],
  gene_promoter = cre_meta$gene_promoter[q],
  geneNameStr = cre_meta$geneNameStr[q],
  rmsk_chr = as.character(seqnames(rmsk_gr[s])),
  rmsk_start = start(rmsk_gr[s]) - 1L,
  rmsk_end = end(rmsk_gr[s]),
  rmsk_strand = as.character(strand(rmsk_gr[s])),
  rmsk_id = rmsk_meta$rmsk_id[s],
  repName = rmsk_meta$repName[s],
  repClass = rmsk_meta$repClass[s],
  repFamily = rmsk_meta$repFamily[s],
  strand_relation = strand_relation,
  overlap_start = overlap_start - 1L,
  overlap_end = overlap_end,
  overlap_len = as.numeric(overlap_len),
  CRE_len = as.numeric(width(cre_gr[q])),
  rmsk_len = as.numeric(width(rmsk_gr[s])),
  cre_index = cre_meta$cre_index[q],
  te_index = rmsk_meta$te_index[s]
)

overlap_all[, `:=`(
  frac_CRE_overlap = overlap_len / CRE_len,
  frac_rmsk_overlap = overlap_len / rmsk_len
)]
setorder(overlap_all, CRE_chr, CRE_start, rmsk_start, repClass, repName)

overlap_filtered <- overlap_all[
  overlap_len >= MIN_OVERLAP_BP & frac_CRE_overlap >= MIN_FRAC_CRE
]
message("[", sample_name, "] Filtered overlap pairs: ", format(nrow(overlap_filtered), big.mark = ","))

if (nrow(overlap_filtered) == 0L) {
  warning("No overlaps passed the significant-overlap thresholds; empty summaries will be written.")
}

# ------------------------------------------------------------
# 7. Summary tables
# ------------------------------------------------------------
total_by_type <- cre[, .(n_CRE_total = .N), by = CRE_type]

raw_by_type <- overlap_all[, .(
  n_CRE_with_any_overlap = uniqueN(cre_index),
  n_raw_overlap_pairs = .N
), by = CRE_type]

filtered_by_type <- overlap_filtered[, .(
  n_CRE_with_significant_overlap = uniqueN(cre_index),
  n_filtered_overlap_pairs = .N,
  total_overlap_bp = sum(overlap_len),
  mean_overlap_bp = mean(overlap_len),
  median_overlap_bp = median(overlap_len),
  mean_frac_CRE_overlap = mean(frac_CRE_overlap),
  median_frac_CRE_overlap = median(frac_CRE_overlap),
  mean_frac_rmsk_overlap = mean(frac_rmsk_overlap),
  median_frac_rmsk_overlap = median(frac_rmsk_overlap)
), by = CRE_type]

summary_by_type <- merge(total_by_type, raw_by_type, by = "CRE_type", all.x = TRUE)
summary_by_type <- merge(summary_by_type, filtered_by_type, by = "CRE_type", all.x = TRUE)

count_cols <- c(
  "n_CRE_with_any_overlap", "n_raw_overlap_pairs",
  "n_CRE_with_significant_overlap", "n_filtered_overlap_pairs",
  "total_overlap_bp"
)
for (column_name in intersect(count_cols, names(summary_by_type))) {
  set(summary_by_type, which(is.na(summary_by_type[[column_name]])), column_name, 0)
}

summary_by_type[, `:=`(
  frac_CRE_with_any_overlap = n_CRE_with_any_overlap / n_CRE_total,
  percent_CRE_with_any_overlap = 100 * n_CRE_with_any_overlap / n_CRE_total,
  frac_CRE_with_significant_overlap = n_CRE_with_significant_overlap / n_CRE_total,
  percent_CRE_with_significant_overlap = 100 * n_CRE_with_significant_overlap / n_CRE_total
)]
setorder(summary_by_type, CRE_type)

summary_overall <- data.table(
  sample = sample_name,
  strand_mode = STRAND_MODE,
  min_overlap_bp = MIN_OVERLAP_BP,
  min_frac_CRE_overlap = MIN_FRAC_CRE,
  retained_TE_class_pattern = TE_CLASS_PATTERN,
  n_CRE_total = nrow(cre),
  n_promoter = cre[CRE_type == "promoter", .N],
  n_enhancer = cre[CRE_type == "enhancer", .N],
  n_CTCF = cre[CRE_type == "CTCF-alone", .N],
  n_unclassed = cre[CRE_type == "unclassed", .N],
  n_other_labels_excluded = nrow(cre_other),
  n_unknown_strand_excluded = nrow(cre_unknown_strand),
  n_CRE_with_any_overlap = uniqueN(overlap_all$cre_index),
  n_CRE_with_significant_overlap = uniqueN(overlap_filtered$cre_index),
  n_raw_overlap_pairs = nrow(overlap_all),
  n_filtered_overlap_pairs = nrow(overlap_filtered)
)
summary_overall[, `:=`(
  frac_CRE_with_any_overlap = n_CRE_with_any_overlap / n_CRE_total,
  frac_CRE_with_significant_overlap = n_CRE_with_significant_overlap / n_CRE_total
)]

summarize_repeats <- function(group_columns) {
  overlap_filtered[, .(
    n_overlap_pairs = .N,
    n_CRE = uniqueN(cre_index),
    total_overlap_bp = sum(overlap_len),
    mean_overlap_bp = mean(overlap_len),
    median_overlap_bp = median(overlap_len),
    mean_frac_CRE_overlap = mean(frac_CRE_overlap),
    median_frac_CRE_overlap = median(frac_CRE_overlap),
    mean_frac_rmsk_overlap = mean(frac_rmsk_overlap),
    median_frac_rmsk_overlap = median(frac_rmsk_overlap)
  ), by = group_columns]
}

summary_repClass <- summarize_repeats(c("CRE_type", "strand_relation", "repClass"))
summary_repFamily <- summarize_repeats(c("CRE_type", "strand_relation", "repClass", "repFamily"))
summary_repName <- summarize_repeats(c("CRE_type", "strand_relation", "repClass", "repFamily", "repName"))

setorderv(summary_repClass, c("CRE_type", "strand_relation", "n_CRE", "n_overlap_pairs"), c(1L, 1L, -1L, -1L))
setorderv(summary_repFamily, c("CRE_type", "strand_relation", "n_CRE", "n_overlap_pairs"), c(1L, 1L, -1L, -1L))
setorderv(summary_repName, c("CRE_type", "strand_relation", "n_CRE", "n_overlap_pairs"), c(1L, 1L, -1L, -1L))

hit_by_type_class <- overlap_filtered[, .(
  n_CRE_hit = uniqueN(cre_index),
  n_overlap_pairs = .N,
  total_overlap_bp = sum(overlap_len),
  mean_frac_CRE_overlap = mean(frac_CRE_overlap),
  median_frac_CRE_overlap = median(frac_CRE_overlap),
  mean_frac_rmsk_overlap = mean(frac_rmsk_overlap),
  median_frac_rmsk_overlap = median(frac_rmsk_overlap)
), by = .(CRE_type, strand_relation, repClass)]

hit_rate_by_type_class <- merge(hit_by_type_class, total_by_type, by = "CRE_type", all.x = TRUE)
hit_rate_by_type_class[, `:=`(
  frac_CRE_hit = n_CRE_hit / n_CRE_total,
  percent_CRE_hit = 100 * n_CRE_hit / n_CRE_total
)]
setorderv(hit_rate_by_type_class, c("CRE_type", "strand_relation", "percent_CRE_hit"), c(1L, 1L, -1L))

# ------------------------------------------------------------
# 8. Write results
# ------------------------------------------------------------
message("[", sample_name, "] Writing output files...")

fwrite(overlap_all, output_file("CRE_TE_overlaps_all_raw.tsv.gz"), sep = "\t")
fwrite(overlap_filtered, output_file("CRE_TE_overlaps_filtered.tsv.gz"), sep = "\t")

for (cre_group in c("promoter", "enhancer", "CTCF-alone", "unclassed")) {
  file_label <- if (cre_group == "CTCF-alone") "CTCF" else cre_group
  fwrite(
    overlap_filtered[CRE_type == cre_group],
    output_file(paste0(file_label, "_TE_overlaps_filtered.tsv.gz")),
    sep = "\t"
  )
}

fwrite(summary_overall, output_file("CRE_TE_overall_summary.tsv"), sep = "\t")
fwrite(summary_by_type, output_file("CRE_TE_summary_by_CRE_type.tsv"), sep = "\t")
fwrite(summary_repClass, output_file("CRE_TE_summary_by_repClass.tsv"), sep = "\t")
fwrite(summary_repFamily, output_file("CRE_TE_summary_by_repFamily.tsv"), sep = "\t")
fwrite(summary_repName, output_file("CRE_TE_summary_by_repName.tsv"), sep = "\t")
fwrite(hit_rate_by_type_class, output_file("CRE_TE_percent_hit_by_repClass.tsv"), sep = "\t")

# ------------------------------------------------------------
# 9. Console report
# ------------------------------------------------------------
cat("\n============================================================\n")
cat("[", sample_name, "] CRE/TE analysis completed successfully\n", sep = "")
cat("============================================================\n")
cat("Strand mode:                      ", STRAND_MODE, "\n")
cat("Minimum overlap:                  ", MIN_OVERLAP_BP, " bp\n")
cat("Minimum CRE fraction:             ", MIN_FRAC_CRE, "\n")
cat("CREs analyzed:                    ", format(nrow(cre), big.mark = ","), "\n")
cat("Promoters:                        ", format(cre[CRE_type == "promoter", .N], big.mark = ","), "\n")
cat("Enhancers:                        ", format(cre[CRE_type == "enhancer", .N], big.mark = ","), "\n")
cat("CTCF-alone:                       ", format(cre[CRE_type == "CTCF-alone", .N], big.mark = ","), "\n")
cat("Unclassed:                        ", format(cre[CRE_type == "unclassed", .N], big.mark = ","), "\n")
cat("Unexpected labels excluded:       ", format(nrow(cre_other), big.mark = ","), "\n")
cat("Unknown-strand CREs excluded:     ", format(nrow(cre_unknown_strand), big.mark = ","), "\n")
cat("Raw overlap pairs:                ", format(nrow(overlap_all), big.mark = ","), "\n")
cat("Filtered overlap pairs:           ", format(nrow(overlap_filtered), big.mark = ","), "\n")
cat("CREs with significant overlap:    ", format(uniqueN(overlap_filtered$cre_index), big.mark = ","), "\n")
cat("\nSummary by CRE type:\n")
print(summary_by_type)
cat("\nOutput directory: ", normalizePath(out_dir), "\n", sep = "")
