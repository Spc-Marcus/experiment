#!/bin/bash
#SBATCH --job-name=ClusteringH
#SBATCH --output=res.txt
#SBATCH --mem=15G

. /local/env/envconda.sh
conda activate strainminer

if python test_pre_post.py; then
    echo "✅ CSV creation completed successfully!"
else
    echo "❌ CSV creation failed"
    exit 1
fi