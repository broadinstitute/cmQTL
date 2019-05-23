#!/usr/bin/env Rscript

'normalize

Usage:
  normalize.R -p <id> -s <query> [-c] [-r <op>] [-t <dir>]

Options:
  -h --help                     Show this screen.
  -p <id> --plate_id=<id>       Plate ID.
  -r <op> --operation=<op>      Normalization operation [default: robustize].
  -s <query> --subset=<query>   Query to specify the sample for estimating normalization parameters.
  -c <sc> --sc_type=<sc>        What type of sc data to load [default: colony]
  -t <dir> --tmpdir=<dir>       Temporary directory [default: /tmp]' -> doc

suppressWarnings(suppressMessages(library(docopt)))

suppressWarnings(suppressMessages(library(dplyr)))

suppressWarnings(suppressMessages(library(magrittr)))

suppressWarnings(suppressMessages(library(stringr)))

opts <- docopt(doc)

plate_id <- opts[["plate_id"]]
operation <- opts[["operation"]]
subset <- opts[["subset"]] # e.g. "Metadata_broad_sample_type == '''control'''"
sc_type <- opts[["sc_type"]]

# load profiles
profiles <- suppressMessages(
  readr::read_tsv(
    file.path(
      "data",
      paste0(plate_id, "_single_cell_", sc_type, "_profiles.tsv.gz")
    )
  )
)

# compartment tag converts nuclei to ^Nuclei_
compartment_tag <- function(compartment) {
  str_c("^", str_sub(compartment, 1, 1) %>% str_to_upper(), str_sub(compartment, 2), "_")

}

load_profiles <- function(compartment) {
  profiles %>%
    select(matches(str_c("Metadata_", "|", compartment_tag(compartment))))

}

normalize_profiles <- function(compartment) {

  sample <- load_profiles(compartment = compartment)

  variables <- colnames(sample) %>% str_subset(compartment_tag(compartment))

  sample %<>%
    filter_(subset) %>%
    collect(n=Inf) %>%
    mutate_at(variables, as.double)

  normalized <-
    cytominer::normalize(
      population = load_profiles(compartment = compartment),
      variables = variables,
      strata =  c("Metadata_Plate"),
      sample = sample,
      operation = operation
    )

}

metadata <-
  colnames(profiles) %>% str_subset("^Metadata_")

normalized <-
  normalize_profiles("cells") %>%
  inner_join(normalize_profiles("cytoplasm"),
             by = metadata) %>%
  inner_join(normalize_profiles("nuclei"),
             by = metadata)

normalized %>%
  readr::write_csv(
    file.path(
      "data",
      paste0(plate_id, "single_cell_", sc_type, "_normalized.csv")
    ),
    sep = "/"
  )
