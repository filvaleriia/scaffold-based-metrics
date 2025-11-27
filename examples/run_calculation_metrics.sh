#!/bin/bash

# Run metric calculation for a single generator
python3 ../src/metrics.py \
    --type_cluster sim \
    --type_scaffold csk \
    --generator Molpher_62.5k \
    --receptor Glucocorticoid_receptor \
    --data_folder ../../ \
    --ncpus 1


