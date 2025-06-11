#!/bin/bash
#SBATCH --job-name=create_matrix
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "=== Step 1: Creating CSV Matrices ==="
for threshold in 0.5 0.7 0.8 0.9 1; do
    echo "Processing threshold: $threshold"
    if python create_csv/process_bam_folder.py bam . --threshold $threshold > res$threshold.txt ; then
        echo "✅ CSV creation completed successfully for $threshold!"
    else
        echo "❌ CSV creation failed for threshold $threshold"
        exit 1
    fi
done