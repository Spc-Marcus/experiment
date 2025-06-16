#!/bin/bash

echo "Starting batch processing of matrices..."

# Parse arguments
SOLVER="gurobi"
while [[ $# -gt 0 ]]; do
    case $1 in
        --solver)
            SOLVER="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--solver gurobi|max_ones|max_ones_comp] [matrix_directory]"
            exit 0
            ;;
        *)
            MATRIX_DIR="$1"
            shift
            ;;
    esac
done

# Validate solver
case $SOLVER in
    gurobi|max_ones|max_ones_comp)
        echo "Using solver: $SOLVER"
        ;;
    *)
        echo "Error: Invalid solver '$SOLVER'"
        exit 1
        ;;
esac

export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Find matrices directory if not provided
if [ -z "$MATRIX_DIR" ]; then
    for matrices_path in "../matrices" "matrices" "./matrices"; do
        if [ -d "$matrices_path" ]; then
            local_csv_count=0
            for hap in {2..10}; do
                if [ -d "$matrices_path/$hap" ]; then
                    hap_csv_count=$(find "$matrices_path/$hap" -name "*.csv" | wc -l)
                    local_csv_count=$((local_csv_count + hap_csv_count))
                fi
            done
            direct_csv_count=$(find "$matrices_path" -maxdepth 1 -name "*.csv" | wc -l)
            local_csv_count=$((local_csv_count + direct_csv_count))
            
            if [ $local_csv_count -gt 0 ]; then
                MATRIX_DIR="$matrices_path"
                break
            fi
        fi
    done
fi

if [ -z "$MATRIX_DIR" ]; then
    echo "No matrix directory found"
    exit 1
fi

# Run experiment
python ilphaplo/run_experiments.py --solver "$SOLVER" "$MATRIX_DIR"
    echo "No matrix data directories found"
    exit 1
fi

if [ $csv_count -eq 0 ]; then
    echo "No CSV files found"
    exit 1
fi

# Run experiment
echo "Running experiment with $SOLVER solver..."
if [ -n "$MATRIX_DIR" ]; then
    python ilphaplo/run_experiments.py --solver "$SOLVER" "$MATRIX_DIR"
else
    python ilphaplo/run_experiments.py --solver "$SOLVER" "$found_dir"
fi

# Check results
if ls experiment_results_${SOLVER}_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t experiment_results_${SOLVER}_*.csv | head -n1)
    echo "✅ Results saved to: $latest_result"
else
    echo "❌ No results file found"
    exit 1
fi
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

# Run the experiment with selected solver
echo ""
echo "Running experiment script with $SOLVER solver..."
if [ -n "$MATRIX_DIR" ]; then
    # Use specified directory
    if python ilphaplo/run_experiments.py --solver "$SOLVER" "$MATRIX_DIR"; then
        echo "✅ Experiment script completed successfully!"
    else
        echo "❌ Experiment script failed"
        exit 1
    fi
else
    # Auto-detect directory
    if python ilphaplo/run_experiments.py --solver "$SOLVER" "$found_dir"; then
        echo "✅ Experiment script completed successfully!"
    else
        echo "❌ Experiment script failed"
        exit 1
    fi
fi

# Check if results were generated
if ls experiment_results_${SOLVER}_*.csv 1> /dev/null 2>&1; then
    latest_result=$(ls -t experiment_results_${SOLVER}_*.csv | head -n1)
    echo "Experiments completed with $SOLVER solver. Results saved to: $latest_result"
    
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
    echo "Error: No results file found for $SOLVER solver"
    echo "Check the console output above for error details"
    exit 1
fi
