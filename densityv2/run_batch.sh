#!/bin/bash

echo "🔍 Starting batch processing of target matrices for density constraint validation..."
echo "Working directory: $(pwd)"

# Set Python path to include experiment directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Check for target files list
target_file="filepaths_with_ilp_calls.txt"
if [ ! -f "$target_file" ]; then
    echo "❌ Error: Target files list not found at $target_file"
    echo "Please ensure filepaths_with_ilp_calls.txt exists in the current directory"
    exit 1
fi

target_count=$(wc -l < "$target_file")
echo "📋 Found target files list with $target_count files to process"

# Check for matrices directory
found_dir=""
for matrices_path in "../matrices" "matrices" "./matrices"; do
    if [ -d "$matrices_path" ]; then
        echo "📁 Found matrices directory: $matrices_path"
        found_dir="$matrices_path"
        break
    fi
done

if [ -z "$found_dir" ]; then
    echo "❌ No matrices directory found"
    echo "Expected locations: ../matrices, matrices, ./matrices"
    exit 1
fi

# Verify some target files exist
sample_files=($(head -5 "$target_file"))
found_samples=0
for file in "${sample_files[@]}"; do
    if [ -f "$found_dir/$file" ]; then
        found_samples=$((found_samples + 1))
    fi
done

echo "✅ Verified $found_samples out of 5 sample target files exist in $found_dir"

if [ $found_samples -eq 0 ]; then
    echo "⚠️  Warning: No sample target files found in matrices directory"
    echo "Please check that the file paths in filepaths_with_ilp_calls.txt match the matrices directory structure"
fi

# Run the density constraint validation experiment
echo ""
echo "🚀 Running density constraint validation experiment..."
echo "Processing $target_count matrices to check for density violations..."

if python densityv2/run_experiments.py "$found_dir"; then
    echo "✅ Density validation experiment completed successfully!"
else
    echo "❌ Density validation experiment failed"
    exit 1
fi

# Check if constraint violations were found
constraint_files=$(ls constraint_violations_*.csv 2>/dev/null || true)
if [ -n "$constraint_files" ]; then
    latest_violations=$(ls -t constraint_violations_*.csv | head -n1)
    echo ""
    echo "⚠️  CONSTRAINT VIOLATIONS DETECTED!"
    echo "📋 Violations logged to: $latest_violations"
    
    violation_count=$(tail -n +2 "$latest_violations" | wc -l)
    echo "🔢 Total constraint violations found: $violation_count"
    
    # Count unique problematic matrices
    matrices_list=$(ls constraint_violations_*_matrices_list.txt 2>/dev/null | head -n1)
    if [ -f "$matrices_list" ]; then
        matrix_count=$(wc -l < "$matrices_list")
        echo "🔴 Matrices with violations: $matrix_count out of $target_count"
        echo "📄 List of problematic matrices: $matrices_list"
        
        # Show first few problematic matrices
        echo ""
        echo "🔍 First 5 problematic matrices:"
        head -5 "$matrices_list"
        if [ $matrix_count -gt 5 ]; then
            echo "... and $((matrix_count - 5)) more"
        fi
    fi
    
    # Check if experiment results were generated (only for problematic matrices)
    if ls experiment_results_*.csv 1> /dev/null 2>&1; then
        latest_result=$(ls -t experiment_results_*.csv | head -n1)
        echo "📊 Detailed analysis of problematic matrices: $latest_result"
    fi
else
    echo ""
    echo "🎉 PERFECT RESULT!"
    echo "✅ No constraint violations found - all zones satisfy density constraints"
    echo "🎯 All $target_count matrices passed density validation"
    echo "ℹ️  No results file created (no problematic matrices to report)"
fi

echo ""
echo "📝 Summary:"
echo "   - Total matrices processed: $target_count"
if [ -n "$constraint_files" ]; then
    violation_matrices=$(wc -l < "$matrices_list" 2>/dev/null || echo "0")
    clean_matrices=$((target_count - violation_matrices))
    echo "   - Clean matrices (no violations): $clean_matrices"
    echo "   - Problematic matrices: $violation_matrices"
    echo "   - Success rate: $(( clean_matrices * 100 / target_count ))%"
else
    echo "   - Clean matrices (no violations): $target_count"
    echo "   - Problematic matrices: 0"
    echo "   - Success rate: 100%"
fi

echo ""
echo "🏁 Density constraint validation completed!"
