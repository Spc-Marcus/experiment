#!/bin/bash
#SBATCH --job-name=SeedStrainminer
#SBATCH --output=res.txt
#SBATCH --ntasks=1
#SBATCH --mem=12G

echo "🚀 Starting ILP Seeding Parameter Experiments (Haplotypes >= 4)..."
echo "Working directory: $(pwd)"

# Set Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd):$(pwd)/.."

# Step 1: Test seeding parameter support
echo ""
echo "Step 1: Testing seeding parameter support..."
if python test_seeding_support.py; then
    echo "✅ Seeding support test completed!"
else
    echo "⚠️  Seeding support test detected issues!"
    echo "Experiments will run with default clustering parameters."
fi

# Step 2: Find matrix directories with haplotype >= 4
echo ""
echo "Step 2: Searching for matrix data directories (haplotype >= 4)..."

found_dir=""
csv_count=0

for matrices_path in "../matrices" "matrices" "./matrices"; do
    if [ -d "$matrices_path" ]; then
        echo "Found matrices directory: $matrices_path"
        local_csv_count=0
        
        # Check for numbered subdirectories >= 4
        for hap in {4..10}; do
            if [ -d "$matrices_path/$hap" ]; then
                hap_csv_count=$(find "$matrices_path/$hap" -name "*.csv" | wc -l)
                local_csv_count=$((local_csv_count + hap_csv_count))
                echo "  Haplotype $hap: $hap_csv_count CSV files"
            fi
        done
        
        # Also check for direct CSV files (will be filtered by haplotype count)
        direct_csv_count=$(find "$matrices_path" -maxdepth 1 -name "*.csv" | wc -l)
        if [ $direct_csv_count -gt 0 ]; then
            echo "  Direct CSV files: $direct_csv_count (will be filtered)"
            local_csv_count=$((local_csv_count + direct_csv_count))
        fi
        
        if [ $local_csv_count -gt 0 ]; then
            found_dir="$matrices_path"
            csv_count=$local_csv_count
            echo "  Total CSV files: $csv_count"
            break
        fi
    fi
done

if [ -z "$found_dir" ]; then
    echo "❌ No matrix data directories found with haplotype >= 4"
    exit 1
fi

if [ $csv_count -eq 0 ]; then
    echo "❌ No CSV files found in matrix directories with haplotype >= 4"
    exit 1
fi

# Step 3: Run seeding experiments
echo ""
echo "Step 3: Running seeding experiments (haplotypes >= 4 only)..."
echo "Parameters to test:"
echo "  X_factor: [2, 3, 4, 5, 6]"
echo "  step_n: [2, 4, 6, 8, 10, 12, 14]"
echo "  Condition: Save only if ILP calls > 1"
echo "  Filter: Only haplotype count >= 4"
echo "  Target matrices: Size >= 50x20"
echo ""

if python seeding_experiment.py "$found_dir"; then
    echo "✅ Seeding experiments completed successfully!"
else
    echo "❌ Seeding experiments failed"
    exit 1
fi

# Step 4: Check and analyze results
echo ""
echo "Step 4: Checking results..."

if ls seeding_experiment_results_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t seeding_experiment_results_*.csv | head -n1)
    echo "✅ Results saved to: $latest_result"
    
    # Show summary
    echo ""
    echo "=== QUICK SUMMARY (HAPLOTYPES >= 4) ==="
    python -c "
import pandas as pd
try:
    df = pd.read_csv('$latest_result')
    print(f'✅ Total experiments with ILP > 1: {len(df)}')
    if len(df) > 0:
        print(f'📊 Haplotype counts: {sorted(df[\"nb_haplotypes\"].unique())}')
        print(f'📊 X_factor values tested: {sorted(df[\"X_factor\"].unique())}')
        print(f'📊 step_n values tested: {sorted(df[\"step_n\"].unique())}')
        print(f'📊 Average ILP calls: {df[\"nb_ilp_calculated\"].mean():.2f}')
        print(f'📊 Average total time: {df[\"total_time\"].mean():.3f}s')
        print(f'📊 Files processed: {df[\"read_name\"].nunique()}')
        
        # Check seeding usage
        if 'seeding_used' in df.columns:
            seeding_used = df['seeding_used'].sum()
            print(f'📊 Experiments using seeding: {seeding_used} out of {len(df)}')
            if seeding_used == 0:
                print('⚠️  WARNING: No experiments used seeding parameters!')
                print('   All results represent baseline performance only.')
        
        # Best performing combination
        if len(df) > 5:
            best_combo = df.groupby(['X_factor', 'step_n'])['total_time'].mean().idxmin()
            best_time = df.groupby(['X_factor', 'step_n'])['total_time'].mean().min()
            print(f'🏆 Best combination: X_factor={best_combo[0]}, step_n={best_combo[1]} (avg time: {best_time:.3f}s)')
    else:
        print('⚠️  No experiments met the ILP > 1 criteria')
        print('   This means all matrices (haplotypes >= 4) were solved efficiently without ILP')
except Exception as e:
    print(f'❌ Error reading results: {e}')
"
else
    echo "❌ No results file found"
    exit 1
fi

echo ""
echo "🎉 Seeding experiments completed (haplotypes >= 4)!"
echo ""
echo "Next steps:"
echo "1. Analyze the detailed results in: $latest_result"
echo "2. If seeding wasn't used, check clustering.py parameter support"
echo "3. Generate plots for visualization"
