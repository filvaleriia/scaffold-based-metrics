#!/bin/bash

# Posibility to add all mean values Generator to one DataFrame
python3 ../src/metrics_connection.py \
    --type_cluster sim \
    --type_scaffold csk \
    --generator_list Molpher DrugEx_GT_epsilon_0.6 REINVENT DrugEx_GT_epsilon_0.1 \
    --receptor Glucocorticoid_receptor \
    --subset _62.5k \
    --data_folder ../../