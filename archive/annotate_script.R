# This a hack to annotate this plate, because cytotools annotate would need us to overspecify the metadata
# We only need cell line and dose information

library(tidyverse)
library(magrittr)

df <- 
  read_csv("../../backend/PILOT_1/BR00098071/BR00098071.csv") %>%
  select(matches("Metadata")) %>%
  mutate(well_col = as.integer(str_sub(Metadata_Well, 2, 3)))

df %<>% 
  mutate(
    Metadata_cell_line = 
      case_when(
        well_col >= 1  & well_col <= 4  ~ "A",
        well_col >= 5  & well_col <= 8  ~ "B",
        well_col >= 9  & well_col <= 12 ~ "C",
        well_col >= 13 & well_col <= 16 ~ "D",
        well_col >= 17 & well_col <= 20 ~ "E",
        well_col >= 21 & well_col <= 24 ~ "F"
      )
  )

seeding_density_level <-
  tribble(~well_col, ~Metadata_seeding_density_level,
          1, 1,
          2, 2,
          3, 3,
          4, 4,
          5, 4,
          6, 3,
          7, 2,
          8, 1,
          9, 1,
          10, 2,
          11, 3,
          12, 4,
          13, 4,
          14, 3,
          15, 2,
          16, 1,
          17, 1,
          18, 2,
          19, 3,
          20, 4,
          21, 4,
          22, 3,
          23, 2,
          24, 1
  ) %>%
    inner_join(
      tribble(~Metadata_seeding_density_level, ~Metadata_confluence,
              1, 40,
              2, 60,
              3, 80,
              4, 100
              )
    )

df %<>%
  inner_join(seeding_density_level) %>%
  select(-well_col)


read_csv("../../backend/PILOT_1/BR00098071/BR00098071.csv") %>%
  inner_join(df) %>%
  write_csv("../../backend/PILOT_1/BR00098071/BR00098071_augmented.csv")
  
library(viridis)
library(platetools)

p <- 
  raw_map(data = df$Metadata_cell_line,
        well = df$Metadata_Well,
        plate = 384) +
  ggtitle("PILOT_1") +
  theme_dark() +
  scale_fill_discrete()

p

ggsave("cell_line_layout.pdf", width = 8)

p <- 
  raw_map(data = df$Metadata_confluence,
        well = df$Metadata_Well,
        plate = 384) +
  ggtitle("PILOT_1") +
  theme_dark() +
  scale_fill_continuous()

p

ggsave("density_layout.pdf", width = 8)
