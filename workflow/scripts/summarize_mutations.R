#!/usr/bin/env Rscript

# Redirect all stdout/stderr (including package-attach messages, warnings,
# and message() calls) into the Snakemake log file instead of the console.
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

#############################
## INPUT / OUTPUT (SNAKEMAKE)
#############################

input_dir <- dirname(snakemake@input[["done"]])
gene_map_file <- snakemake@input[["gene_map"]]  # NULL if not declared in the rule

output_long <- snakemake@output[["long"]]
output_wide <- snakemake@output[["wide"]]

input_files <- list.files(
  input_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

# Label used when a mapped gene_id was never seen anywhere in the input
# data at all (as opposed to "Wildtype", which means the gene WAS seen but
# this particular sample had no fixed-alt call there).
NOT_FOUND_LABEL <- "ERROR_NOT_FOUND"

#############################
## GENE MAPPING (OPTIONAL)
#############################

# Gene mapping is optional. It's considered "not in use" if:
#   - the gene_map input wasn't declared at all (config key absent/commented out), or
#   - the file exists but has zero data rows (header-only)
# In that case we run in pass-through mode: no gene filtering, and
# gene_name is just set equal to gene_id.

use_gene_mapping <- !is.null(gene_map_file)

gene_map <- NULL

if (use_gene_mapping) {

  gene_map <- read_tsv(gene_map_file, show_col_types = FALSE) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), str_trim))

  required_map_cols <- c("gene_id", "gene_name")
  if (!all(required_map_cols %in% names(gene_map))) {
    stop(
      "Gene mapping file '", gene_map_file, "' must contain columns: ",
      paste(required_map_cols, collapse = ", "), ". Found: ",
      paste(names(gene_map), collapse = ", ")
    )
  }

  if (nrow(gene_map) == 0) {
    message(
      "Gene mapping file '", gene_map_file, "' has no rows - ",
      "running without gene filtering (pass-through mode)."
    )
    use_gene_mapping <- FALSE
    gene_map <- NULL
  } else {

    if (any(duplicated(gene_map$gene_id))) {
      dupes <- unique(gene_map$gene_id[duplicated(gene_map$gene_id)])
      stop("Duplicate gene_id values in gene mapping file: ", paste(dupes, collapse = ", "))
    }

    if (any(duplicated(gene_map$gene_name))) {
      dupes <- unique(gene_map$gene_name[duplicated(gene_map$gene_name)])
      stop(
        "Duplicate gene_name values in gene mapping file: ", paste(dupes, collapse = ", "),
        ". Each gene_name must be unique - it becomes a column header in the wide table."
      )
    }
  }
}

#############################
## HELPERS
#############################

is_silent_aa <- function(x) {
  !is.na(x) & grepl("^([A-Z])\\d+\\1$", x)
}

extract_position <- function(mutation) {
  suppressWarnings(
    as.numeric(gsub("^[A-Z]+([0-9]+).*", "\\1", mutation))
  )
}

# A genotype counts as "fixed" (i.e. the sample unambiguously carries the
# alt allele, not just a heterozygous call) if:
#   - it's haploid: rows only exist here for non-ref calls to begin with
#     (the upstream awk script drops 0/. calls), so any haploid row is fixed
#   - it's diploid/polyploid: all alleles in GT are identical and non-ref,
#     e.g. "1/1", "1|1", "2/2", "2|2" - not just the hardcoded "1/1"/"1|1"
is_fixed_alt <- function(gt, ploidy) {
  diploid_alleles <- str_split(gt, "[/|]")
  is_diploid_homozygous <- map_lgl(diploid_alleles, function(a) {
    length(a) > 1 && length(unique(a)) == 1 && a[1] != "0"
  })
  ploidy == "haploid" | is_diploid_homozygous
}

# Collapse all mutations observed for one (Sample, Gene) group into a single
# label. LOF/NMD take priority over AA changes. Multiple AA changes are
# ordered by residue position (ascending) rather than table row order
summarize_mutations <- function(mutation_label, position) {

  if (any(mutation_label == "LOF_NMD", na.rm = TRUE)) return("LOF_NMD")
  if (any(mutation_label == "LOF", na.rm = TRUE)) return("LOF")
  if (any(mutation_label == "NMD", na.rm = TRUE)) return("NMD")

  aa_idx <- !is.na(mutation_label) & !(mutation_label %in% c("LOF", "NMD", "LOF_NMD"))

  if (!any(aa_idx)) return("Wildtype")

  labs <- mutation_label[aa_idx]
  pos <- position[aa_idx]

  ord <- order(pos, labs, na.last = TRUE)
  ordered_labs <- labs[ord]
  ordered_labs <- ordered_labs[!duplicated(ordered_labs)]

  paste(ordered_labs, collapse = "/")
}

#############################
## CORE FUNCTION
#############################

process_batch <- function(annotated_table, use_gene_mapping, target_gene_ids) {

  if (!"Ploidy" %in% names(annotated_table)) {
    stop(
      "Input table is missing the 'Ploidy' column. ",
      "Regenerate it with the updated vcf2table.sh script."
    )
  }

  mutations2 <- annotated_table %>%
    mutate(
      AA_Change = na_if(AA_Change, ""),
      AA_Change = na_if(AA_Change, " "),
      LOF = as.integer(LOF),
      NMD = as.integer(NMD)
    ) %>%
    mutate(
      Is_silent = ifelse(!is.na(AA_Change), is_silent_aa(AA_Change), FALSE),
      Mutation_label = case_when(
        LOF == 1 & NMD == 1 ~ "LOF_NMD",
        LOF == 1 ~ "LOF",
        NMD == 1 ~ "NMD",
        !is.na(AA_Change) & !Is_silent ~ AA_Change,
        TRUE ~ NA_character_
      )
    )

  # Capture the full sample list BEFORE any gene filtering, so samples that
  # are wildtype across the whole target gene panel are still known to exist
  # (and end up as explicit "Wildtype" rows/cells) instead of silently
  # disappearing from the output.
  all_samples_batch <- unique(mutations2$Sample)

  # Capture every raw Gene ID that appears in this batch BEFORE the
  # gene-of-interest filter is applied. Used both for the FOUND/NOT FOUND
  # log report and to distinguish "confirmed Wildtype" from "gene_id never
  # seen at all" in the final table.
  raw_genes_batch <- unique(mutations2$Gene)

  mutations2 <- mutations2 %>%
    mutate(
      Mutation_for_order = ifelse(!is.na(AA_Change), AA_Change, "Wildtype"),
      Position = extract_position(Mutation_for_order)
    ) %>%
    filter(is_fixed_alt(GT, Ploidy))

  if (use_gene_mapping) {
    mutations2 <- mutations2 %>% filter(Gene %in% target_gene_ids)
  }

  genes_observed_batch <- unique(mutations2$Gene)

  mutation_summary <- mutations2 %>%
    group_by(Sample, Gene) %>%
    summarise(
      Final_mutation = summarize_mutations(Mutation_label, Position),
      .groups = "drop"
    ) %>%
    arrange(Sample)

  list(
    summary = mutation_summary,
    samples = all_samples_batch,
    genes = genes_observed_batch,
    raw_genes = raw_genes_batch
  )
}

#############################
## BATCH PROCESSING
#############################

target_gene_ids <- if (use_gene_mapping) gene_map$gene_id else NULL

batch_results <- map(input_files, function(f) {

  message("Processing ", basename(f))

  annotated_table <- read_tsv(f, show_col_types = FALSE) %>%
    mutate(across(everything(), as.character))

  process_batch(annotated_table, use_gene_mapping, target_gene_ids)
})

long_summary <- map(batch_results, "summary") %>%
  bind_rows() %>%
  distinct()

all_samples <- map(batch_results, "samples") %>%
  unlist() %>%
  unique()

# Every raw Gene ID seen anywhere in the input data, regardless of the
# fixed-alt/gene-of-interest filters. Used to tell "confirmed Wildtype"
# apart from "gene_id never found at all".
raw_genes_all <- map(batch_results, "raw_genes") %>% unlist() %>% unique()

#############################
## FINAL GENE PANEL
#############################

if (use_gene_mapping) {

  message("---- Gene mapping match status (", nrow(gene_map), " genes in ", gene_map_file, ") ----")
  for (i in seq_len(nrow(gene_map))) {
    gid <- gene_map$gene_id[i]
    gname <- gene_map$gene_name[i]
    if (gid %in% raw_genes_all) {
      message(sprintf("  FOUND     %-20s (%s)", gid, gname))
    } else {
      message(sprintf(
        "  NOT FOUND %-20s (%s) - check for typos, a 'gene-' prefix, case mismatch, or genuinely zero variants",
        gid, gname
      ))
    }
  }
  message("----------------------------------------------------------------")

  # Predefined panel: every gene in the mapping file, even ones with zero
  # observed variants anywhere, so they still show up in the output.
  final_gene_map <- gene_map

} else {
  # No mapping supplied: the "panel" is just whatever genes were actually
  # observed across all batches. gene_name = gene_id (identity mapping).
  # Every gene here was by definition found in the data, so the
  # NOT_FOUND_LABEL case never triggers in this branch.
  all_genes_observed <- map(batch_results, "genes") %>% unlist() %>% unique()
  final_gene_map <- tibble(gene_id = all_genes_observed, gene_name = all_genes_observed)
}

#############################
## BUILD FULL SAMPLE x GENE GRID
#############################

full_grid <- expand_grid(
  Sample = all_samples,
  gene_id = final_gene_map$gene_id
) %>%
  left_join(final_gene_map, by = "gene_id") %>%
  mutate(Gene_found_in_data = gene_id %in% raw_genes_all)

final_long <- full_grid %>%
  left_join(
    long_summary %>% rename(gene_id = Gene),
    by = c("Sample", "gene_id")
  ) %>%
  mutate(
    Final_mutation = case_when(
      !is.na(Final_mutation) ~ Final_mutation,
      Gene_found_in_data ~ "Wildtype",
      TRUE ~ NOT_FOUND_LABEL
    )
  ) %>%
  select(Sample, Gene_id = gene_id, Gene_name = gene_name, Final_mutation) %>%
  arrange(Sample, Gene_name)

final_wide <- final_long %>%
  select(Sample, Gene_name, Final_mutation) %>%
  pivot_wider(
    names_from = Gene_name,
    values_from = Final_mutation,
    values_fill = "Wildtype"
  ) %>%
  arrange(Sample)

write_tsv(final_long, output_long)
write_tsv(final_wide, output_wide)