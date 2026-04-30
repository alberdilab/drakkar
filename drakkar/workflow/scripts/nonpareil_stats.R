#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    stop(paste("Missing required argument", flag), call. = FALSE)
  }
  args[[idx + 1]]
}

sample_name <- get_arg("--sample")
npo_file <- get_arg("--npo")
output_file <- get_arg("--output")

fields <- c("kappa", "C", "LR", "modelR", "LRstar", "diversity")

write_stats <- function(values) {
  row <- as.data.frame(as.list(values), check.names = FALSE)
  row <- cbind(sample = sample_name, row)
  write.table(row, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!file.exists(npo_file) || file.info(npo_file)$size == 0) {
  values <- setNames(rep("0", length(fields)), fields)
  write_stats(values)
  quit(status = 0)
}

suppressPackageStartupMessages(library(Nonpareil))

normalise_name <- function(value) {
  gsub("[^a-z0-9]", "", tolower(value))
}

extract_value <- function(stats, candidates) {
  candidate_names <- normalise_name(candidates)

  if (is.null(dim(stats))) {
    stat_names <- names(stats)
    if (!is.null(stat_names)) {
      idx <- match(candidate_names, normalise_name(stat_names))
      idx <- idx[!is.na(idx)]
      if (length(idx) > 0) {
        return(unname(as.character(stats[[idx[[1]]]])))
      }
    }
    return("NA")
  }

  col_names <- colnames(stats)
  if (!is.null(col_names)) {
    idx <- match(candidate_names, normalise_name(col_names))
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      return(unname(as.character(stats[1, idx[[1]]])))
    }
  }

  row_names <- rownames(stats)
  if (!is.null(row_names)) {
    idx <- match(candidate_names, normalise_name(row_names))
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) {
      return(unname(as.character(stats[idx[[1]], 1])))
    }
  }

  "NA"
}

curve <- Nonpareil.curve(
  npo_file,
  plot = FALSE,
  label = sample_name,
  enforce.consistency = FALSE
)
stats <- summary(curve)

values <- c(
  kappa = extract_value(stats, c("kappa")),
  C = extract_value(stats, c("C")),
  LR = extract_value(stats, c("LR")),
  modelR = extract_value(stats, c("modelR", "modelRt")),
  LRstar = extract_value(stats, c("LRstar", "LR*")),
  diversity = extract_value(stats, c("diversity"))
)

write_stats(values)
