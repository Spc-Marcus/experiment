#!/bin/bash

echo "=== Running Solver Comparison Experiments ==="

# Set up environment
if [ -n "$SLURM_JOB_ID" ]; then
    if [ -n "$SLURM_SUBMIT_DIR" ]; then
        cd "$SLURM_SUBMIT_DIR"
    fi
else
    cd "$(dirname "$0")"
fi

export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Find matrices directory
matrices_found=false
for matrices_path in "matrices" "../matrices" "./matrices" "matrice" "../matrice" "./matrice"; do
    if [ -d "$matrices_path" ]; then
        csv_count=0
        for hap in {2..10}; do
            if [ -d "$matrices_path/$hap" ]; then
                hap_csv_count=$(find "$matrices_path/$hap" -name "*.csv" 2>/dev/null | wc -l)
                csv_count=$((csv_count + hap_csv_count))
            fi
        done
        direct_csv_count=$(find "$matrices_path" -maxdepth 1 -name "*.csv" 2>/dev/null | wc -l)
        csv_count=$((csv_count + direct_csv_count))
        
        if [ $csv_count -gt 0 ]; then
            matrices_found=true
            matrices_dir="$matrices_path"
            break
        fi
    fi
done

if [ "$matrices_found" = false ]; then
    echo "Error: No matrices directory found"
    exit 1
fi

# Run experiments for each solver
solvers=("gurobi" "max_ones" "max_ones_comp")

for solver in "${solvers[@]}"; do
    echo "Testing $solver solver..."
    bash ILPV2/run_batch.sh --solver "$solver" "$matrices_dir"
done

echo "Done!"
        echo "✅ $solver experiments completed successfully"
    else
        echo "❌ $solver experiments failed"
    fi
done

