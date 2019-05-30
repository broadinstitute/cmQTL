
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(corrplot))

set.seed(123)

# Read in both btches of data
profile_dir <- file.path("..", "..", "..", "backend", "2019_05_13_Batch2")
profile_files <- list.files(profile_dir,
                            recursive = TRUE,
                            full.names = TRUE,
                            pattern = "variable_selected.csv")

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

profile_list <- list()
for (norm_file in profile_files) {
    if (grepl("colony", norm_file)) {
        profile_type <- "colony"
    } else if (grepl("isolated", norm_file)) {
        profile_type <- "isolated"
    } else {
        profile_type <- "all"
    }
    
    profile_list[[norm_file]] <-
    readr::read_csv(norm_file, col_types = profile_cols) %>%
        dplyr::mutate(Metadata_profile_type = profile_type)
    
}

profile_df <- do.call(rbind, profile_list) %>% tibble::remove_rownames()

# Note that with the addition of single cell data, some variables recieve NA
head(sort(colSums(is.na(profile_df)), decreasing = TRUE))

# We need to drop these
profile_df <- profile_df %>% dplyr::select_if(~ !any(is.na(.)))

dim(profile_df)
table(profile_df$Metadata_profile_type)
head(profile_df, 2)

# Separate different cell profiler data
cp_features <- colnames(profile_df) %>%
    stringr::str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

length(cp_features)

cp_metadata <- colnames(profile_df) %>%
    stringr::str_subset("^Metadata_")

length(cp_metadata)

# Create a metadata dictionary and dummy variable "group_id"
# "group_id" distinguishes each separate condition including cell line
# "condition_group_id" distinguishes separate conditions ignoring cell line
metadata_df <- profile_df %>%
    dplyr::select(cp_metadata) %>%
    dplyr::mutate(dictionary_id = paste0("id_", dplyr::row_number()),
                  group_id = group_indices(.,
                                           Metadata_profile_type,
                                           Metadata_line_ID,
                                           Metadata_plating_density,
                                           Metadata_timepoint,
                                           Metadata_Plate),
                  condition_group_id = group_indices(.,
                                                     Metadata_profile_type,
                                                     Metadata_plating_density,
                                                     Metadata_timepoint,
                                                     Metadata_Plate))

tail(metadata_df)

# Create a dataframe of variables for each group
group_id_df <- metadata_df %>%
    dplyr::select(group_id,
                  Metadata_profile_type,
                  Metadata_line_ID,
                  Metadata_timepoint,
                  Metadata_plating_density) %>%
    dplyr::distinct() %>%
    dplyr::arrange(group_id)

dim(group_id_df)
tail(group_id_df)

table(
    profile_df$Metadata_profile_type,
    profile_df$Metadata_line_ID,
    profile_df$Metadata_plating_density,
    profile_df$Metadata_timepoint,
    profile_df$Metadata_Plate
)

cor_df <- profile_df %>%
    dplyr::select(cp_features) %>%
    t() %>%
    cor() %>%
    dplyr::as_tibble() %>%
    magrittr::set_colnames(metadata_df$dictionary_id)

cor_melt_df <- metadata_df %>%
    dplyr::select(-group_id,
                  -condition_group_id) %>%
    dplyr::bind_cols(
        replace(cor_df,
                lower.tri(cor_df, TRUE), NA)
    ) %>%
    dplyr::select(-cp_metadata) %>%
    reshape2::melt(id.vars = 'dictionary_id',
                   variable.name = 'correlation_id', 
                   value.name = "pearson_cor",
                   na.rm = TRUE) %>%
    tibble::remove_rownames()

dim(cor_melt_df)
head(cor_melt_df)

# Map group IDs and condition IDs onto the correlation dataframe
# We are interested in correlations between specific groups
cor_group_df <- cor_melt_df %>%
    dplyr::inner_join(
        metadata_df %>%
        select(dictionary_id,
               group_id,
               condition_group_id),
        by = 'dictionary_id'
    ) %>%
    dplyr::rename(pair_a = group_id,
                  pair_a_condition = condition_group_id) %>%
    dplyr::inner_join(
        metadata_df %>%
        select(dictionary_id,
               group_id,
               condition_group_id),
        by = c('correlation_id' = 'dictionary_id')
    ) %>%
    dplyr::rename(pair_b = group_id,
                  pair_a_id = dictionary_id,
                  pair_b_id = correlation_id,
                  pair_b_condition = condition_group_id)

dim(cor_group_df)
head(cor_group_df)

# Remove self correlations and determine median correlation between all groups
# Also create a variable that represents correlations across cell lines within
# the same condition. This variable will be used as the null distribution.
cor_group_df <- cor_group_df %>%
    dplyr::mutate(
        within_group_cor =
            as.numeric(cor_group_df$pair_a == cor_group_df$pair_b),
        within_group_across_cell_line_cor =
            as.numeric(cor_group_df$pair_a_condition == cor_group_df$pair_b_condition &
                       cor_group_df$pair_a != cor_group_df$pair_b)
    ) %>%
    dplyr::filter(cor_group_df$pair_a_id != cor_group_df$pair_b_id) %>%
    dplyr::group_by(
        pair_a,
        pair_b
    ) %>%
    dplyr::mutate(median_cor = median(pearson_cor)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(median_cor))

dim(cor_group_df)
head(cor_group_df)

# Join Replicate Correlations and Null Distribution Correlations
within_group_cor_df <- cor_group_df %>%
    dplyr::filter(within_group_cor == 1) %>%
    dplyr::group_by(pair_b) %>%
    dplyr::mutate(pair_b_median_cor = median(pearson_cor),
                  null_data = "Replicate Correlation") %>%
    dplyr::arrange(desc(pair_b_median_cor)) %>%
    dplyr::ungroup()

null_group_cor_df <- cor_group_df %>%
    dplyr::filter(within_group_across_cell_line_cor == 1) %>%
    dplyr::group_by(pair_b) %>%
    dplyr::mutate(pair_b_median_cor = median(pearson_cor),
                  null_data = "Matched Cell Lines") %>%
    dplyr::arrange(desc(pair_b_median_cor)) %>%
    dplyr::ungroup()

full_plot_ready <- within_group_cor_df %>%
    dplyr::bind_rows(
        null_group_cor_df
    )

dim(full_plot_ready)
head(full_plot_ready)

# Perform KS tests between real and null distributions
all_results <- list()
for (group_condition in unique(full_plot_ready$pair_b)) {
    full_plot_group_df = full_plot_ready %>%
        dplyr::filter(pair_b == group_condition)

    replicate_corr_df <- full_plot_group_df %>%
        dplyr::filter(within_group_across_cell_line_cor == 0)
    null_corr_df <- full_plot_group_df %>%
        dplyr::filter(within_group_across_cell_line_cor == 1)

    ks_result = ks.test(x = replicate_corr_df$pearson_cor,
                        y = null_corr_df$pearson_cor,
                        alternative = "less")

    k_stat = as.numeric(ks_result$statistic)
    k_p = as.numeric(ks_result$p.value)
    all_results[[group_condition]] <- c(group_condition, k_stat, k_p, -log10(k_p)) 
}

ks_result_df <- dplyr::as_tibble(do.call(rbind, all_results))
colnames(ks_result_df) <- c("group_condition_id", "ks_stat", "ks_p_value", "ks_log_10_p")

ks_result_df <- ks_result_df %>% dplyr::arrange(desc(as.numeric(paste(ks_log_10_p))))

dim(ks_result_df)
head(ks_result_df)

# Merge plot ready data with info on group ID
full_plot_ready <- full_plot_ready %>%
    dplyr::left_join(group_id_df, by = c("pair_b" = "group_id"))

# Sort the plot ready dataframe in order of KS significance
full_plot_ready$pair_b <- factor(full_plot_ready$pair_b,
                                 levels = ks_result_df$group_condition_id)

full_plot_ready <- full_plot_ready %>%
    dplyr::mutate(grouped_plot = paste0(Metadata_timepoint,
                                        null_data))

for (sc_type in unique(full_plot_ready$Metadata_profile_type)) {
    sc_full_plot_ready <- full_plot_ready %>%
        dplyr::filter(Metadata_profile_type == sc_type)
    
    condition_gg <-
    ggplot(sc_full_plot_ready,
           aes(x = pearson_cor,
               y = as.factor(Metadata_timepoint),
               fill = null_data)) +
        geom_density_ridges(alpha = 0.8) +
        ylab("Conditions") +
        xlab("Profile Correlation (Pearson)") +
        theme_bw() +
        ggtitle(paste("Correlations within", sc_type)) +
        facet_grid(Metadata_line_ID~Metadata_plating_density,
                   scales="free_y") +
        scale_fill_manual(name = "",
                          values = c("#FFC107", "#004D40")) +
        theme(axis.text.y = element_text(size = 12),
              axis.text.x = element_text(size = 8),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              strip.text = element_text(size = 10),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))

    print(condition_gg)
    
    file_base <- file.path("figures", paste0("condition_pilot_", sc_type))
    for (extension in c('.png', '.pdf')) {
        ggsave(condition_gg,
               filename = paste0(file_base, extension),
               height = 8,
               width = 10)
    }
}

ks_result_df$group_condition_id <- as.factor(ks_result_df$group_condition_id)

final_results_df <- full_plot_ready %>%
    dplyr::left_join(ks_result_df,
                     by = c("pair_b" = "group_condition_id")) %>%
    dplyr::group_by(
        Metadata_profile_type,
        Metadata_timepoint,
        Metadata_plating_density,
        Metadata_line_ID
    ) %>%
    dplyr::mutate(
        median_cell_ks = median(ks_stat)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(
        Metadata_profile_type,
        Metadata_timepoint,
        Metadata_plating_density
    ) %>%
    dplyr::mutate(
        median_condition_ks = median(ks_stat)
    ) %>%
    dplyr::select(
        Metadata_profile_type,
        Metadata_timepoint,
        Metadata_plating_density,
        Metadata_line_ID,
        median_cell_ks,
        median_condition_ks
    ) %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

append_timepoint <- function(string) paste("Timepoint:", string)

ks_test_gg <- ggplot(final_results_df) +
    geom_jitter(aes(y = median_cell_ks,
                    x = as.factor(Metadata_plating_density),
                    color = Metadata_line_ID),
                size = 1.25,
                height = 0,
                width = 0.2,
                alpha = 0.8) +
    geom_point(aes(y = median_condition_ks,
                   x = as.factor(Metadata_plating_density)),
               color = "red",
               shape = 23,
               size = 3) +
    scale_color_manual(name = "Cell ID",
                       values = c("#a6cee3",
                                  "#d95f02",
                                  "#7570b3",
                                  "#e7298a",
                                  "#66a61e",
                                  "#e6ab02")) +
    facet_grid(Metadata_profile_type~Metadata_timepoint,
               scales = "free_x",
               labeller = labeller(Metadata_timepoint = as_labeller(append_timepoint))) +
    xlab("Plating Density") +
    ylab("KS Statistic") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 8),
          axis.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 9),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

ks_test_gg

file_base <- file.path("figures", "condition_pilot_kstest")
for (extension in c('.png', '.pdf')) {
    ggsave(ks_test_gg,
           filename = paste0(file_base, extension),
           height = 6,
           width = 6)
}

final_results_df <- full_plot_ready %>%
    dplyr::left_join(ks_result_df,
                     by = c("pair_b" = "group_condition_id")) %>%
    dplyr::group_by(
        Metadata_profile_type,
        Metadata_timepoint,
        Metadata_plating_density
    ) %>%
    dplyr::mutate(
        median_ks = median(ks_stat),
        low_ks = min(ks_stat),
        high_ks = max(ks_stat)
    ) %>%
    dplyr::select(
        Metadata_profile_type,
        Metadata_timepoint,
        Metadata_plating_density,
        median_ks,
        low_ks,
        high_ks
    ) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(median_ks))

file <- file.path("results", "pilot_condition_results.tsv")
readr::write_tsv(final_results_df, file)

final_results_df
