#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "🧬 StrainMiner - KNN Complete Pipeline"
echo "======================================"

# Create CSV
echo "=== Creating CSV Matrices ==="
python create_csv/process_bam_folder.py bam .

# KNN Experiment
echo "=== KNN Experiment Mode ==="
python Knn/impute_matrices.py matrices_no_binarize --experiment

# KNN Normal
echo "=== KNN Normal Mode ==="
python Knn/impute_matrices.py matrices_no_binarize

# ilphaplo
echo "=== ilphaplo Analysis ==="
bash ilphaplo/run_batch.sh

echo "🎉 Complete KNN pipeline finished!"