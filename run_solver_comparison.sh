#!/bin/bash
#SBATCH --job-name=strainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=10G

. /local/env/envconda.sh
conda activate strainminer
. /local/env/envsamtools-1.15.sh

echo "=== Running Solver Comparison Experiments ==="
echo "This will test all three solvers: gurobi, max_ones, max_ones_comp"

# Set up environment
cd "$(dirname "$0")"
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Create results directory
mkdir -p comparison_results
timestamp=$(date +%Y%m%d_%H%M%S)

# Array of solvers to test
solvers=("max_ones" "max_ones_comp")

echo ""
echo "Starting solver comparison at $(date)"
echo "Results will be saved with timestamp: $timestamp"

# Run experiments for each solver
for solver in "${solvers[@]}"; do
    echo ""
    echo "=== Testing $solver solver ==="
    
    if bash ILPV2/run_batch.sh --solver "$solver"; then
        echo "✅ $solver experiments completed successfully"
        
        # Move results to comparison directory
        if ls experiment_results_${solver}_*.csv 1> /dev/null 2>&1; then
            latest_result=$(ls -t experiment_results_${solver}_*.csv | head -n1)
            cp "$latest_result" "comparison_results/results_${solver}_${timestamp}.csv"
            echo "Results copied to: comparison_results/results_${solver}_${timestamp}.csv"
        fi
    else
        echo "❌ $solver experiments failed"
    fi
done
echo "✅ All solver comparisons completed"