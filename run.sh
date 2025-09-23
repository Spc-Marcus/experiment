#!/bin/bash
#SBATCH --job-name=Test_max_e_wr
#SBATCH --output=test.txt
#SBATCH --ntasks=1
#SBATCH --mem=12G

. /local/env/envconda.sh
conda activate strainminer


echo "Test run max e wr (avec warm start)"
echo "==================================="

if python ilp_time/test_all.py; then
    echo "✅ ILP optimization completed successfully!"
else
    echo "❌ ILP optimization failed"
    exit 1
fi