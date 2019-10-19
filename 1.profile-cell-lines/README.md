# Profile Stem Cell Lines

Our goal with this module is to generate profiles for all the 480 stem cell lines. 

We will process data in batches and update this readme with notes as we go.

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

