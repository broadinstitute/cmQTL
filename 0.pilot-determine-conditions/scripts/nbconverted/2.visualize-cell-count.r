
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

set.seed(123)

# Read in well counts
sc_col_types <- readr::cols(
    Metadata_Well = readr::col_character(),
    cell_count = readr::col_integer(),
    batch_id = readr::col_character(),
    sc_type = readr::col_character()
)

file <- file.path("results", "well_cell_counts.tsv")
num_sc_df <- readr::read_tsv(file, col_types = sc_col_types)

dim(num_sc_df)
head(num_sc_df, 2)

# Read in example file and merge with well
profile_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_Assay_Plate_Barcode = readr::col_character(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_well_position = readr::col_character(),
    Metadata_plating_density = readr::col_integer(),
    Metadata_line_ID = readr::col_character(),
    Metadata_timepoint = readr::col_integer()
)

batches <- c("BR00103267", "BR00103268")

eg_list <- list()
for (batch in batches) {
    file <- file.path("data", paste0(batch, "_normalized_variable_selected.csv"))
    eg_list[[batch]] <- readr::read_csv(file, col_types = profile_cols)
}

eg_df <- do.call(rbind, eg_list)
eg_df <- eg_df %>% dplyr::select(starts_with("Metadata_"))

dim(eg_df)
head(eg_df, 2)

# Merge with well
num_sc_df <- num_sc_df %>% 
    dplyr::left_join(eg_df,
                     by = c("Metadata_Well" = "Metadata_Well",
                            "batch_id" = "Metadata_Assay_Plate_Barcode")) %>%
    dplyr::select(-X1)

head(num_sc_df, 2)

append_timepoint <- function(string) paste("Time:", string)
append_density <- function(string) paste("Density:", string)

ggplot(num_sc_df,
       aes(x = sc_type,
           y = cell_count,
           color = Metadata_line_ID)) +
        geom_jitter(size = 0.6,
                    alpha = 0.6,
                    height = 0,
                    width = 0.3) +
        facet_wrap(Metadata_timepoint ~ Metadata_plating_density,
                   ncol = 4, 
                   scales = "free_y",
                   labeller = labeller(Metadata_timepoint = as_labeller(append_timepoint),
                                       Metadata_plating_density = as_labeller(append_density))) +
        scale_color_manual(name = "Cell ID",
                           values = c("#a6cee3",
                                      "#d95f02",
                                      "#7570b3",
                                      "#e7298a",
                                      "#66a61e",
                                      "#e6ab02")) +
        ylab("Total Number of Cells Per Well") +
        xlab("Cell Context") +
        theme_bw() + 
        theme(axis.text.y = element_text(size = 11),
              axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5),
              axis.title = element_text(size = 11),
              legend.text = element_text(size = 10),
              strip.text = element_text(size = 8),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))

file_base <- file.path("figures", paste0("sc_num_cells_pilot"))
for (extension in c('.png', '.pdf')) {
    ggsave(filename = paste0(file_base, extension),
           height = 5,
           width = 8)
}

file <- file.path("results", "well_cell_counts_with_metadata.tsv")
readr::write_tsv(num_sc_df, file)
