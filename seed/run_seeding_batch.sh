#!/bin/bash

echo "Starting ILP Seeding Parameter Experiments..."
echo "Working directory: $(pwd)"

# Set Python path to include experiment directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Check for existing matrix directories
echo "Searching for matrix data directories..."

found_dir=""
csv_count=0

# Check for matrix directories
for matrices_path in "../matrices" "matrices" "./matrices"; do
    if [ -d "$matrices_path" ]; then
        echo "Found matrices directory: $matrices_path"
        
        # Count CSV files
        local_csv_count=$(find "$matrices_path" -name "*.csv" | wc -l)
        
        if [ $local_csv_count -gt 0 ]; then
            found_dir="$matrices_path"
            csv_count=$local_csv_count
            echo "  Total CSV files: $csv_count"
            break
        fi
    fi
done

if [ -z "$found_dir" ]; then
    echo "No matrix data directories found"
    exit 1
fi

if [ $csv_count -eq 0 ]; then
    echo "Warning: No CSV files found in matrix directories"
    exit 1
fi

# Run seeding experiments
echo ""
echo "Running seeding experiments..."
echo "Parameters to test:"
echo "  X_factor: [2, 3, 4, 5, 6]"
echo "  step_n: [2, 4, 6, 8, 10, 12, 14]"
echo "  Condition: Save only if ILP calls > 1"
echo ""

if python seeding_experiment.py "$found_dir"; then
    echo "✅ Seeding experiments completed successfully!"
else
    echo "❌ Seeding experiments failed"
    exit 1
fi

# Check if results were generated
if ls seeding_experiment_results_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t seeding_experiment_results_*.csv | head -n1)
    echo "Seeding experiments completed. Results saved to: $latest_result"
    
    # Show summary
    echo ""
    echo "=== EXPERIMENT SUMMARY ==="
    python -c "
import pandas as pd
try:
    df = pd.read_csv('$latest_result')
    print(f'Total experiments with ILP > 1: {len(df)}')
    print(f'X_factor values: {sorted(df[\"X_factor\"].unique())}')
    print(f'step_n values: {sorted(df[\"step_n\"].unique())}')
    print(f'Average ILP calls: {df[\"nb_ilp_calculated\"].mean():.2f}')
    print(f'Average total time: {df[\"total_time\"].mean():.3f}s')
except Exception as e:
    print(f'Error reading results: {e}')
"
else
    echo "Error: No seeding results file found"
    exit 1
fi
