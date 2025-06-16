#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "🚀 StrainMiner - Sequential Pipeline"

# Step 2: KNN
echo "=== KNN Imputation ==="
python Knn/impute_matrices.py matrices_no_binarize

# Step 3: Run solver comparison
echo "=== Solver Comparison ==="
bash run_solver_comparison.sh

echo "🎉 Done!"