#!/bin/bash
#
# cmQTL processing pipeline
# Shantanu Singh, 2019
#
# Here, we copy the profiles generated from the backend into this repo

echo "batch_id,plate_id
2019_06_10_Batch3,cmqtlpl1.5-31-2019-mt
2019_06_10_Batch3,cmqtlpl261-2019-mt
2019_08_15_Batch4,BR00106709
2019_08_15_Batch4,BR00106708
2019_09_06_Batch5,BR00107338
2019_09_06_Batch5,BR00107339" > /tmp/batch_plate.txt

parallel \
    -a batch_plate.txt \
    --keep-order \
    --max-procs 2 \
    --header ".*\n" \
    -C "," \
    cp ../../../backend/{1}/{2}/{2}{3}.csv profiles/ ::: \
    _augmented _colony _colony_augmented _colony_normalized _colony_normalized_variable_selected _count _isolated _isolated_augmented _isolated_normalized _isolated_normalized_variable_selected _normalized _normalized_variable_selected
