#!/bin/bash

echo "Starting batch processing of target matrices..."
echo "Working directory: $(pwd)"

# Set Python path to include experiment directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Check for target files list
target_file="filepaths_with_ilp_calls.txt"
if [ ! -f "$target_file" ]; then
    echo "Error: Target files list not found at $target_file"
    echo "Please ensure filepaths_with_ilp_calls.txt exists in the ilphaplo directory"
    exit 1
fi

target_count=$(wc -l < "$target_file")
echo "Found target files list with $target_count files to process"

# Check for matrices directory
found_dir=""
for matrices_path in "../matrices" "matrices" "./matrices"; do
    if [ -d "$matrices_path" ]; then
        echo "Found matrices directory: $matrices_path"
        found_dir="$matrices_path"
        break
    fi
done

if [ -z "$found_dir" ]; then
    echo "No matrices directory found"
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

echo "Verified $found_samples out of 5 sample target files exist in $found_dir"

if [ $found_samples -eq 0 ]; then
    echo "Warning: No sample target files found in matrices directory"
    echo "Please check that the file paths in filepaths_with_ilp_calls.txt match the matrices directory structure"
fi

# Run the experiment - use correct path to Python script
echo ""
echo "Running density constraint validation experiment..."
if python ilphaplo/run_experiments.py "$found_dir"; then
    echo "✅ Experiment script completed successfully!"
else
    echo "❌ Experiment script failed"
    exit 1
fi

# Check if results were generated
if ls experiment_results_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t experiment_results_*.csv | head -n1)
    echo "Experiments completed. Results saved to: $latest_result"
    
    # Show constraint validation summary
    constraint_files=$(ls constraint_violations_*.csv 2>/dev/null || true)
    if [ -n "$constraint_files" ]; then
        latest_violations=$(ls -t constraint_violations_*.csv | head -n1)
        echo "Constraint violations logged to: $latest_violations"
        violation_count=$(tail -n +2 "$latest_violations" | wc -l)
        echo "Total constraint violations found: $violation_count"
    else
        echo "✅ No constraint violations found - all zones satisfy density constraints"
    fi
    
    # Run analysis if analyze_results.py exists
    if [ -f "ilphaplo/analyze_results.py" ]; then
        echo "Running analysis..."
        python ilphaplo/analyze_results.py "$latest_result" --plots
        echo "Analysis completed. Check the plots/ directory for visualizations."
    else
        echo "Analysis script not found, skipping visualization step."
        echo "You can run analysis manually with:"
        echo "python ilphaplo/analyze_results.py $latest_result --plots"
    fi
else
    echo "Error: No results file found"
    exit 1
fi
    echo ""
    echo "Please ensure matrices are available in numbered subdirectories"
    exit 1
fi

if [ $csv_count -eq 0 ]; then
    echo "Warning: No CSV files found in matrix directories"
    echo "Please ensure your data files have .csv extension"
    exit 1
fi

# Run the experiment - use correct path to Python script
echo ""
echo "Running experiment script..."
if python ilphaplo/run_experiments.py "$found_dir"; then
    echo "✅ Experiment script completed successfully!"
else
    echo "❌ Experiment script failed"
    exit 1
fi

# Check if results were generated
if ls experiment_results_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t experiment_results_*.csv | head -n1)
    echo "Experiments completed. Results saved to: $latest_result"
    
    # Show first few lines of results for verification
    echo ""
    echo "=== FIRST FEW RESULTS ==="
    head -n 3 "$latest_result"
    echo ""
    
    # Run analysis if analyze_results.py exists
    if [ -f "ilphaplo/analyze_results.py" ]; then
        echo "Running analysis..."
        python ilphaplo/analyze_results.py "$latest_result" --plots
        echo "Analysis completed. Check the plots/ directory for visualizations."
    else
        echo "Analysis script not found, skipping visualization step."
        echo "You can run analysis manually with:"
        echo "python ilphaplo/analyze_results.py $latest_result --plots"
    fi
else
    echo "Error: No results file found"
    echo "Check the console output above for error details"
    exit 1
fi
