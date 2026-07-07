#!/usr/bin/env Rscript

# Redirect all stdout/stderr (including package-attach messages, warnings,
# and message() calls) into the Snakemake log file instead of the console.
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file)
sink(log_file, type = "message")

# Load required libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(tidyverse)
# Get command line arguments
# args <- commandArgs(trailingOnly = TRUE)

input_files <- list(snakemake@input[["snpeff"]])
genes_file <- snakemake@input[["genes"]]
output_file <- snakemake@output[["wide"]]

base_output <- sub("\\.tsv$", "", output_file)

#############################
## READ DATA
#############################

annotated_table <- map_dfr(
  input_files,
  ~ read_tsv(.x, show_col_types = FALSE) %>%
      mutate(across(everything(), as.character))
)

gene_map <- read_tsv(
  genes_file,
  show_col_types = FALSE
)
colnames(gene_map) <- c("Gene", "Name")

genes_of_interest <- gene_map$Name

#############################
## INITIAL FILTERING on functionality
#############################
is_silent_aa <- function(x) {
  grepl("^([A-Z])\\d+\\1$", x) & !is.na(x)
}

mutations2 <- annotated_table %>%
  mutate(
    AA_Change = ifelse(AA_Change %in% c(".", " ", ""), NA, AA_Change),
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

#############################
## HELPER FUNCTION
#############################

extract_position <- function(mutation) {
  suppressWarnings(as.numeric(case_when(
    grepl("^TR[0-9]+$", mutation) ~ 0,
    TRUE ~ as.numeric(gsub("^[A-Z]+([0-9]+).*", "\\1", mutation))
  )))
}

#############################
## GENE NAME HANDLING CHANGE
#############################

final_mapped <- mutations2 %>%
  left_join(gene_map, by = "Gene")

# replace Gene_final with Name and otherwise Gene
mutations3 <- final_mapped %>%
  mutate(
    Gene_final = ifelse(!is.na(Name) & Name != "", Name, Gene)
  )

all_samples <- unique(mutations3$Sample)

#############################
## TR LOGIC FOR cyp51A
#############################

mutations3 <- mutations3 %>%
  mutate(
    start = as.integer(str_match(HGVS.c, "c\\.(-?\\d+)_")[, 2]),
    ins_seq = str_match(HGVS.c, "ins([ACGT]+)")[, 2],
    Mutation_label = if_else(
      Gene_final == "cyp51A" &
      !is.na(ins_seq) &
      !is.na(start) &
      abs(start) >= 250 &
      abs(start) <= 350,
      paste0("TR", nchar(ins_seq)),
      Mutation_label
    ),
    AA_Change = if_else(
      Gene_final == "cyp51A" &
      !is.na(ins_seq) &
      !is.na(start) &
      abs(start) >= 250 &
      abs(start) <= 350,
      paste0("TR", nchar(ins_seq)),
      AA_Change
    )
  ) %>%
  select(-start, -ins_seq)

mutations3 <- mutations3 %>%
  mutate(
    Mutation_for_order = ifelse(
      !is.na(AA_Change),
      AA_Change,
      "Wildtype"
    ),
    Position = extract_position(Mutation_for_order)
  )

### Filter on AA_Change/HGVS.p for relevant genes ###


mutations3 <- mutations3 %>%
  filter(Gene_final %in% genes_of_interest)

## message about heterozygous mutations
het_mutations <- mutations3 %>%
  filter(!GT %in% c("1/1", "1|1")) %>% filter(!is.na(Mutation_label) & Mutation_label != "Wildtype")
if (nrow(het_mutations) > 0) {
  message("Heterozygous mutations detected in samples, see tsv output for details: see heterozygous_mutations.tsv")
}

write_tsv(
  het_mutations,
  paste0(output_file, "_heterozygous_mutations.tsv")
)

############################
## MUTATION IDS
#############################

########## REMOVE non 1/1 #######
# mutations3 <- mutations3 %>%
#   filter(GT %in% c("1/1", "1|1"))

mut_rows <- mutations3 %>%
  select(Sample, Gene_final, AA_Change, GT)

mutations_combined <- mutations3 %>%
  distinct() %>%
  mutate(
    Mutation_ID = paste0(Gene_final, ".", AA_Change)
  ) #%>%
#mutate(across(everything(), as.character))

## Add wildtype rows for samples without mutations ###
samples_missing <- setdiff(
  all_samples,
  unique(mutations_combined$Sample)
)

wildtype_rows <- data.frame(
  Sample = samples_missing,
  Gene_final = rep(NA_character_, length(samples_missing)),
  AA_Change = rep(NA_character_, length(samples_missing)),
  Mutation_ID = rep(NA_character_, length(samples_missing)),
  stringsAsFactors = FALSE
)

mutations_combined <- bind_rows(mutations_combined, wildtype_rows) %>%  distinct() 

#############################
## MUTATION SUMMARY
#############################

mutation_summary <- mutations3 %>%
  group_by(Sample, Gene_final) %>%
  arrange(Position, .by_group = TRUE) %>%
  summarise(
    Final_mutation = case_when(
      any(Mutation_label == "LOF_NMD") ~ "LOF_NMD",
      any(Mutation_label == "LOF") ~ "LOF",
      any(Mutation_label == "NMD") ~ "NMD",
      any(!is.na(Mutation_label)) ~
        paste(unique(Mutation_label[!is.na(Mutation_label)]), collapse = "/"),
      TRUE ~ "Wildtype"
    ),
    .groups = "drop"
  ) %>%
  complete(
    Sample = all_samples,
    Gene_final = genes_of_interest,
    fill = list(Final_mutation = "Wildtype")
  ) %>%
  pivot_wider(
    names_from = Gene_final,
    values_from = Final_mutation
  )

final_output <- mutation_summary %>%
  arrange(Sample)

write_tsv(
  final_output,
  output_file
)

#############################
## PRESENCE / ABSENCE MATRIX
#############################

mutations3_matrix <- mutations_combined %>%
  select(Sample, Mutation_ID)

mutation_matrix <- mutations3_matrix %>%
  mutate(has_mut = 1) %>%
  distinct(Sample, Mutation_ID, .keep_all = TRUE) %>%
  pivot_wider(
    names_from = Mutation_ID,
    values_from = has_mut,
    values_fill = 0
  ) %>%
  select(-any_of("NA"))

# ---- filter mutation columns (same logic as before) ----
mut_cols <- colnames(mutation_matrix)[-1]

cols_to_keep <- c(
  "Sample",
  mut_cols[sapply(mutation_matrix[mut_cols], function(x) {
    sum(as.numeric(x), na.rm = TRUE) >= 1
  })]
)

mutation_matrix <- mutation_matrix %>%
  select(all_of(cols_to_keep))

write_tsv(
  mutation_matrix,
  paste0(base_output, "_mutation_presence_absence_matrix_all.tsv")
)
