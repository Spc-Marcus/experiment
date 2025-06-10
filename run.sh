#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "🚀 StrainMiner - Sequential Pipeline"
echo "==================================="

# Step 1: Create CSV
echo "=== Step 1: Creating CSV Matrices ==="
if python create_csv/process_bam_folder.py bam .; then
    echo "✅ CSV creation completed successfully!"
else
    echo "❌ CSV creation failed"
    exit 1
fi

# Step 2: KNN
echo "=== Step 2: KNN Imputation ==="
if python Knn/impute_matrices.py matrices_no_binarize; then
    echo "✅ KNN imputation completed successfully!"
    echo "Matrices saved in: matrices/"
else
    echo "❌ KNN imputation failed"
    exit 1
fi

# Step 3: ilphaplo analysis
echo "=== Step 3: ilphaplo Analysis ==="
if bash ilphaplo/run_batch.sh; then
    echo "✅ ilphaplo analysis completed successfully!"
else
    echo "❌ ilphaplo analysis failed"
    echo "Note: This may be normal if no suitable matrices were found"
fi

echo "🎉 Sequential pipeline finished!"