#!/bin/bash

echo "Starting batch processing of matrices..."
echo "Working directory: $(pwd)"

# Set Python path to include experiment directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Check for existing matrix directories
echo "Searching for matrix data directories..."

found_dir=""
csv_count=0

# Check KNN output directory first - look for numbered subdirectories
for matrices_path in "../matrices" "matrices" "./matrices"; do
    if [ -d "$matrices_path" ]; then
        echo "Found KNN matrices directory: $matrices_path"
        
        # Check for numbered subdirectories (haplotype counts)
        subdirs=""
        local_csv_count=0
        for hap in {2..10}; do
            if [ -d "$matrices_path/$hap" ]; then
                subdirs="$subdirs $hap"
                hap_csv_count=$(find "$matrices_path/$hap" -name "*.csv" | wc -l)
                local_csv_count=$((local_csv_count + hap_csv_count))
                echo "  Haplotype $hap: $hap_csv_count CSV files"
            fi
        done
        
        # Also check for direct CSV files in matrices directory
        direct_csv_count=$(find "$matrices_path" -maxdepth 1 -name "*.csv" | wc -l)
        if [ $direct_csv_count -gt 0 ]; then
            echo "  Direct CSV files: $direct_csv_count"
            local_csv_count=$((local_csv_count + direct_csv_count))
        fi
        
        if [ $local_csv_count -gt 0 ]; then
            found_dir="$matrices_path"
            csv_count=$local_csv_count
            if [ ! -z "$subdirs" ]; then
                echo "  Haplotype subdirectories found:$subdirs"
            fi
            echo "  Total CSV files: $csv_count"
            break
        fi
    fi
done

# If no KNN matrices, check for organized haplotype directories
if [ -z "$found_dir" ]; then
    for base_path in ".." "."; do
        for dir in matrice matrices data datasets matrix_data haplotype_data; do
            full_path="$base_path/$dir"
            if [ -d "$full_path" ]; then
                echo "Found directory: $full_path"
                
                # Check for numbered subdirectories (2, 3, 4, 5, 6, etc.)
                subdirs=""
                local_csv_count=0
                for hap in {2..10}; do
                    if [ -d "$full_path/$hap" ]; then
                        subdirs="$subdirs $hap"
                        hap_csv_count=$(find "$full_path/$hap" -name "*.csv" | wc -l)
                        local_csv_count=$((local_csv_count + hap_csv_count))
                        echo "  Haplotype $hap: $hap_csv_count CSV files"
                    fi
                done
                
                if [ ! -z "$subdirs" ]; then
                    found_dir="$full_path"
                    csv_count=$local_csv_count
                    echo "  Haplotype subdirectories found:$subdirs"
                    echo "  Total CSV files: $csv_count"
                    break 2
                fi
            fi
        done
    done
fi

if [ -z "$found_dir" ]; then
    echo "No matrix data directories found"
    echo ""
    echo "Expected locations:"
    echo "1. KNN output: matrices/2/, matrices/3/, etc. (from KNN imputation)"
    echo "2. Organized data: matrice/2/, matrice/3/, etc."
    echo ""
    echo "Current directory contents:"
    ls -la
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
if python ilpnoseed/run_experiments.py "$found_dir"; then
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
    if [ -f "ilpnoseed/analyze_results.py" ]; then
        echo "Running analysis..."
        python ilpnoseed/analyze_results.py "$latest_result" --plots
        echo "Analysis completed. Check the plots/ directory for visualizations."
    else
        echo "Analysis script not found, skipping visualization step."
        echo "You can run analysis manually with:"
        echo "python ilpnoseed/analyze_results.py $latest_result --plots"
    fi
else
    echo "Error: No results file found"
    echo "Check the console output above for error details"
    exit 1
fi
