#!/bin/bash

echo "Starting batch processing of matrices..."
echo "Working directory: $(pwd)"

# Set Python path to include experiment directory
export PYTHONPATH="${PYTHONPATH}:$(pwd)/.."

# Check for existing matrix directories
echo "Searching for matrix data directories..."

found_dir=""
csv_count=0

# Check for various possible directory names
for dir in matrice matrices data datasets matrix_data haplotype_data; do
    if [ -d "../$dir" ]; then
        echo "Found directory: ../$dir"
        
        # Check for numbered subdirectories (2, 3, 4, 5, 6)
        subdirs=""
        for hap in {2..8}; do
            if [ -d "../$dir/$hap" ]; then
                subdirs="$subdirs $hap"
                hap_csv_count=$(find "../$dir/$hap" -name "*.csv" | wc -l)
                csv_count=$((csv_count + hap_csv_count))
                echo "  Haplotype $hap: $hap_csv_count CSV files"
            fi
        done
        
        if [ ! -z "$subdirs" ]; then
            found_dir="../$dir"
            echo "  Haplotype subdirectories found:$subdirs"
            echo "  Total CSV files: $csv_count"
            break
        fi
    fi
done

if [ -z "$found_dir" ]; then
    echo "No matrix data directories found with numbered subdirectories (2, 3, 4, 5, 6)"
    echo ""
    echo "Please organize your data as follows:"
    echo "your_data_directory/"
    echo "├── 2/"
    echo "│   ├── matrix1.csv"
    echo "│   └── matrix2.csv"
    echo "├── 3/"
    echo "│   └── matrix3.csv"
    echo "├── 4/"
    echo "├── 5/"
    echo "└── 6/"
    exit 1
fi

if [ $csv_count -eq 0 ]; then
    echo "Warning: No CSV files found in matrix directories"
    echo "Please ensure your data files have .csv extension"
    exit 1
fi

# Run the experiment
echo ""
echo "Running experiment script..."
python run_experiments.py

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
    if [ -f "analyze_results.py" ]; then
        echo "Running analysis..."
        python analyze_results.py "$latest_result" --plots
        echo "Analysis completed. Check the plots/ directory for visualizations."
    else
        echo "Analysis script not found, skipping visualization step."
        echo "You can run analysis manually with:"
        echo "python analyze_results.py $latest_result --plots"
    fi
else
    echo "Error: No results file found"
    echo "Check the console output above for error details"
    echo ""
    echo "Make sure your data is organized as follows:"
    echo "matrice/"
    echo "├── 2/"
    echo "│   ├── matrix1.csv"
    echo "│   └── matrix2.csv"
    echo "├── 3/"
    echo "│   └── matrix3.csv"
    echo "└── ..."
    exit 1
fi
