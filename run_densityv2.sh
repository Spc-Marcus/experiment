#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh
echo "📊 StrainMiner - Density Analysis
=================================="

bash densityv2/run_batch.sh

if [ $? -eq 0 ]; then
    echo "✅ Density analysis completed successfully!"
else
    echo "❌ Density analysis failed"
    exit 1
fi