# ILP Haplotype Analysis Experiment Framework

This framework provides a comprehensive benchmarking system for analyzing the performance of ILP (Integer Linear Programming) based haplotype reconstruction algorithms. It processes binary matrices representing genetic variant data and measures algorithm efficiency, particularly focusing on cases where optimization can be avoided.

## 📋 Overview

// ...existing code...

## 📁 Data Organization

Your matrix data must be organized with numbered haplotype subdirectories:

### Required Structure
```
matrices_bin/                   # Your input directory
├── 1/                         # 1-haplotype matrices  
│   ├── matrix_chr1_0.csv
│   ├── matrix_chr1_5000.csv
│   ├── matrix_chr2_0.csv
│   └── matrix_chr2_5000.csv
├── 2/                         # 2-haplotype matrices
│   ├── matrix_chr1_0.csv
│   ├── matrix_chr1_5000.csv
│   └── matrix_chr3_0.csv
└── 3/                         # 3-haplotype matrices
    ├── matrix_chr1_0.csv
    └── matrix_chr2_0.csv
```

### Matrix Format
- CSV files with binary values (0, 1)
- Rows represent reads/sequences  
- Columns represent variant positions
- First row/column can contain headers

## 🚀 Quick Start

### 1. Basic Usage (Required Arguments)

```bash
# Navigate to experiment directory
cd /udd/mfoin/Dev/refactor/experiment/ILP_call

# Run with required matrix and output directories
./run_batch.sh -m /path/to/matrices_bin -o /path/to/experiment

# Example with specific paths
./run_batch.sh -m /data/matrices_bin -o /results/experiment_2024
```

### 2. Python Usage

```bash
# Run experiments directly
python run_experiments.py -m /path/to/matrices_bin -o /path/to/experiment

# Analyze results (separate step)
cd ../analyze
python analyze_results.py /path/to/experiment/experiment_results_*.csv --plots
```

## 📊 Understanding Results

// ...existing code...

## 🔧 Advanced Usage

### Command Line Options

#### run_batch.sh (Required Arguments)
```bash
./run_batch.sh -m MATRIX_DIR -o OUTPUT_DIR

Required Arguments:
  -m, --matrix-dir DIR    Directory containing haplotype subdirectories (1/, 2/, 3/, etc.)
  -o, --output-dir DIR    Directory to save analysis results

Optional Arguments:  
  -h, --help             Show help message

Examples:
  ./run_batch.sh -m /data/matrices_bin -o /results/experiment
  ./run_batch.sh --matrix-dir ./matrices_bin --output-dir ./output
```

#### run_experiments.py (Required Arguments)
```bash
python run_experiments.py -m MATRIX_DIR -o OUTPUT_DIR

Required Arguments:
  --matrix-dir, -m DIR    Matrix directory path (must contain 1/, 2/, 3/ subdirs)
  --output-dir, -o DIR    Output directory path

Examples:
  python run_experiments.py -m ./matrices_bin -o ./results
  python run_experiments.py --matrix-dir /data/matrices_bin --output-dir /output/exp1
```

#### analyze_results.py (Separate Analysis Step)
```bash
python analyze_results.py RESULTS_CSV [OPTIONS]

Arguments:
  RESULTS_CSV             Path to experiment results CSV file

Options:
  --plots                 Generate visualization plots  
  --output-dir, -o DIR    Directory to save plots (default: plots)

Examples:
  python analyze_results.py /results/experiment_results_20241201_143022.csv --plots
  python analyze_results.py results.csv --plots -o /output/plots
```

// ...existing code...

## 📝 Example Workflow

```bash
# 1. Prepare your data with required structure
mkdir -p /data/matrices_bin/{1,2,3}
# ... organize your CSV files in numbered subdirectories ...

# 2. Run the complete analysis with required arguments
cd /udd/mfoin/Dev/refactor/experiment/ILP_call
./run_batch.sh -m /data/matrices_bin -o /results/my_experiment

# 3. Analyze results with visualizations (separate step)
cd ../analyze  
python analyze_results.py /results/my_experiment/experiment_results_*.csv --plots -o /results/my_experiment/plots

# 4. Review results
ls /results/my_experiment/                    # CSV file with statistics
ls /results/my_experiment/plots/             # PNG files with visualizations
```

## 🔧 Troubleshooting

### Common Issues

#### "Matrix directory does not exist"
- Ensure the provided path exists and is accessible
- Use absolute paths to avoid confusion
- Check directory permissions

#### "No numbered subdirectories found"  
- Ensure your data has subdirectories named 1/, 2/, 3/, etc.
- Check that subdirectories contain CSV files
- Verify the directory structure matches the required format

#### "No CSV files found"
- Check that CSV files exist in the numbered subdirectories
- Ensure files have .csv extension
- Verify file permissions

#### Arguments missing
- Both -m (matrix directory) and -o (output directory) are required
- Use -h flag to see usage help
- Check command syntax

// ...existing code...