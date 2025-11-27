#!/bin/bash

# Run metric calculation for a single generator
python3 ../src/metrics_define_path.py \
    --type_cluster sim \
    --type_scaffold csk \
    --generator Molpher_62.5k \
    --receptor Glucocorticoid_receptor \
    --recall_set_path /home/filv/phd_projects/iga_2023/git_reccal/new/data/input_recall_sets/Glucocorticoid_receptor \
    --output_set_path /home/filv/phd_projects/iga_2023/git_reccal/new/data/output_sets/Glucocorticoid_receptor/Molpher_62.5k \
    --data_folder ../../ \
    --ncpus 1


