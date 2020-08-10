# Profile Stem Cell Lines

Our goal with this module is to generate profiles for 308 stem cell lines. 

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

The single cell SQLite files are available on S3:

|Metadata_Plate|
|:-------------|
| [cmqtlpl1.5-31-2019-mt](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_06_10_Batch3/cmqtlpl1.5-31-2019-mt/cmqtlpl1.5-31-2019-mt.sqlite) |
| [cmqtlpl261-2019-mt](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_06_10_Batch3/cmqtlpl261-2019-mt/cmqtlpl261-2019-mt.sqlite) |
| [BR00106708](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_08_15_Batch4/BR00106708/BR00106708.sqlite) |
| [BR00106709](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_08_15_Batch4/BR00106709/BR00106709.sqlite) |
| [BR00107338](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_09_06_Batch5/BR00107338/BR00107338.sqlite) |
| [BR00107339](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2019_09_06_Batch5/BR00107339/BR00107339.sqlite) |
| [cmQTLplate7-7-22-20](https://imaging-platform.s3.amazonaws.com/projects/2018_06_05_cmQTL/workspace/backend/2020_07_22_Batch7/cmQTLplate7-7-22-20/cmQTLplate7-7-22-20.sqlite) |


These are the counts of cell lines across the full dataset (n=308 in total)

|Metadata_Plate        |  n|
|:---------------------|--:|
|cmqtlpl1.5-31-2019-mt | 48|
|cmqtlpl261-2019-mt    | 48|
|BR00106708            | 48|
|BR00106709            | 45|
|BR00107338            | 40|
|BR00107339            | 48|
|cmQTLplate7-7-22-20   | 48|

These are the details of batches – when they were imaged and the plate IDs of the constituent plates (a.k.a. `Metadata_Plate`)

| Fixation date | Batch ID | Metadata_Plate |
|:------------|----------| :-------|
| 2020-07-22 | 2020_07_22_Batch7 | cmQTLplate7-7-22-20 |
| 2019-08-27 | 2019_09_06_Batch5 | BR00107338 |
| 2019-08-20 | 2019_09_06_Batch5 | BR00107339 |
| 2019-08-09 | 2019_08_15_Batch4 | BR00106708 |
| 2019-08-05 | 2019_08_15_Batch4 | BR00106709 |
| 2019-06-03 | 2019_06_10_Batch3 | cmqtlpl261-2019-mt |
| 2019-05-31 | 2019_06_10_Batch3 | cmqtlpl1.5-31-2019-mt |


## Downloading sample images

The notebook `⁨1.profile-cell-lines⁩/7.select_images_to_print.Rmd` shows how to download sample images.

