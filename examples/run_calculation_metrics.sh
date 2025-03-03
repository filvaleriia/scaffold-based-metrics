# Run metric calculation for a single generator
python3 ../src/metrics.py \
    --type_cluster sim \
    --type_scaffold csk \
    --generator Molpher_125k \
    --receptor Glucocorticoid_receptor \
    --save_option True \
    --ncpus 1


