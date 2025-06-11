#!/bin/bash
#SBATCH --job-name=create_matrix
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "=== Step 1: Creating CSV Matrices ==="
echo "Working directory: $(pwd)"
echo "Available BAM files:"
ls -la bam/

for threshold in 0.5 0.7 0.8 0.9 1; do
    echo "Processing threshold: $threshold"
    echo "Starting at: $(date)"
    
    # Add timeout to prevent hanging
    timeout 1800 python create_csv/process_bam_folder.py bam . --threshold $threshold
    
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "✅ CSV creation completed successfully for $threshold!"
    elif [ $exit_code -eq 124 ]; then
        echo "❌ CSV creation timed out after 30 minutes for threshold $threshold"
        exit 1
    else
        echo "❌ CSV creation failed for threshold $threshold (exit code: $exit_code)"
        exit 1
    fi
    echo "Finished at: $(date)"
    echo ""
done