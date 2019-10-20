# Profile Stem Cell Lines

Our goal with this module is to generate profiles for all the 267 stem cell lines. 

There are 10 profiles files for each plate, stored in [profiles](https://github.com/broadinstitute/cmQTL/tree/master/1.profile-cell-lines/profiles), corresponding to different transformations of the data:

| file | description |
|-------------|---|
|augmented | mean profiles |
|colony_augmented | mean profiles of cells growing in colonies |
|isolated_augmented | mean profiles of isolated cells |
|count | cell counts |
|normalized | z-scored profiles |
|normalized_variable_selected | z-profiles, with redundant variables removed |
|colony_normalized | (similar to normalized) |
|colony_normalized_variable_selected | (similar to normalized_variable_selected) |
|isolated_normalized | (similar to normalized) |
|isolated_normalized_variable_selected | (similar to normalized_variable_selected) |

These are the counts of cell lines across the full dataset (n=267 in total)

|Metadata_Plate        |  n|
|:---------------------|--:|
|BR00106708            | 48|
|BR00106709            | 45|
|BR00107338            | 40|
|BR00107339            | 48|
|cmqtlpl1.5-31-2019-mt | 48|
|cmqtlpl261-2019-mt    | 48|
