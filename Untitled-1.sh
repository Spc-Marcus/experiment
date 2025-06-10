#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "🧬 StrainMiner Pipeline Scripts"
echo "================================"
echo ""
echo "Available scripts:"
echo "  📋 run_pipeline.sh: Complete pipeline (CSV → KNN → ilphaplo)"
echo "  📊 run_create.sh: Create CSV matrices only"
echo "  🤖 run_knn.sh: KNN imputation (both experiment and normal modes)"
echo "  🧩 run.sh: Full workflow with BAM indexing"
echo ""
echo "Usage: sbatch <script_name>"
echo ""
echo "Current directory: $(pwd)"
echo "Available directories:"
ls -la | grep "^d"