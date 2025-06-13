#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh





# Step 3: ilpnoseed analysis
echo "=== Step 3: ilpnoseed Analysis ==="
if bash ilpnoseed/run_batch.sh; then
    echo "✅ ilpnoseed analysis completed successfully!"
else
    echo "❌ ilpnoseed analysis failed"
    echo "Note: This may be normal if no suitable matrices were found"
fi

echo "🎉 Sequential pipeline finished!"
