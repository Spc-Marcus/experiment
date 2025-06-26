# ILP Haplotype Analysis Experiment Package

This package provides tools for analyzing binary haplotype matrices using Integer Linear Programming (ILP) optimization. It performs preprocessing, clustering, and pattern detection to identify haplotype structures in genomic data.

## Overview

The package implements a three-phase pipeline:
1. **Preprocessing**: Identifies inhomogeneous regions using hierarchical clustering
2. **Clustering**: Applies iterative biclustering with ILP optimization 
3. **Analysis**: Collects comprehensive performance statistics and results

## Installation

### Prerequisites

- Python 3.8+
- Gurobi Optimizer (license required)
- Required Python packages:
  ```bash
  pip install numpy pandas scikit-learn gurobipy pathlib
  ```

### Gurobi Setup

This package requires Gurobi with valid licensing. The code includes WLS credentials for academic use, but you may need to configure your own license.

## Usage

### Basic Command

```bash
python run_experiments.py -m <matrix_directory> -o <output_directory>
```

### Example Commands

```bash
# Basic usage with relative paths
python run_experiments.py -m ../data/imputed_matrices -o ../data/experiments

# Full path example
python run_experiments.py --matrix-dir /path/to/matrices_bin --output-dir /path/to/results

# With specific parameters (modify source code for custom parameters)
python run_experiments.py -m ./data/matrices -o ./results
```

### Input Directory Structure

Your matrix directory should be organized as follows:
```
matrices_bin/
├── 1/           # Haplotype count 1
│   ├── matrix1.csv
│   ├── matrix2.csv
│   └── ...
├── 2/           # Haplotype count 2
│   ├── matrix1.csv
│   └── ...
├── 3/           # Haplotype count 3
└── ...
```

Each numbered subdirectory represents a different haplotype count, containing CSV files with binary matrices (0/1 values).

### CSV Matrix Format

Input matrices should be CSV files with:
- First row: column headers (feature names)
- First column: row indices (sample names)  
- Values: binary (0 or 1) indicating presence/absence
- Missing values automatically converted to 0

Example:
```csv
,feature1,feature2,feature3,feature4
sample1,1,0,1,0
sample2,0,1,1,1
sample3,1,1,0,0
```

## Output

### Results File

The analysis generates `ilp_call_stats.csv` in the output directory containing:

#### File Information
- `filename`: Original CSV filename
- `filepath`: Relative path from matrix directory
- `haplotype_count`: Number of haplotypes (from directory structure)
- `load_time`: Time to load and process input

#### Matrix Characteristics  
- `matrix_rows`, `matrix_cols`: Matrix dimensions
- `matrix_size`: Total number of cells
- `matrix_density`: Fraction of 1s in matrix
- `ones_count`, `zeros_count`: Count of 1s and 0s

#### Statistical Summaries
- `row_min_sum`, `row_max_sum`, `row_mean_sum`: Row sum statistics
- `col_min_sum`, `col_max_sum`, `col_mean_sum`: Column sum statistics

#### Performance Metrics
- `ilp_calls_total`: Total ILP optimization calls
- `preprocessing_time`: Time spent in preprocessing phase
- `clustering_time`: Time spent in clustering phase  
- `total_time`: Total execution time
- `ilp_time_total`: Time spent in ILP optimization
- `patterns_found`: Number of significant patterns detected

#### Algorithm Results
- `regions_found`: Number of inhomogeneous regions identified
- `clustering_steps`: Number of clustering iterations
- `avg_region_size`: Average size of processed regions

### Summary Statistics

The tool prints comprehensive statistics including:
- Algorithm efficiency rate (percentage of matrices requiring no ILP)
- Performance breakdown by haplotype count
- Total processing time and ILP calls across all matrices

## File Descriptions

### Core Modules

- **`run_experiments.py`**: Main script for batch processing experiments
- **`pipeline.py`**: Orchestrates the complete analysis pipeline
- **`preprocessing.py`**: Handles matrix preprocessing and region identification
- **`clustering.py`**: Implements ILP-based biclustering algorithms
- **`get_data.py`**: Handles CSV loading and data preparation

### Key Functions

#### `run_experiments(matrix_dir, output_dir)`
Processes all matrices in directory and generates comprehensive statistics.

#### `run_pipeline(binary_matrix, **params)` 
Executes complete preprocessing + clustering pipeline on single matrix.

#### `clustering_full_matrix(matrix, regions, steps, **params)`
Performs exhaustive iterative biclustering using ILP optimization.

#### `find_quasi_biclique(matrix, error_rate)`
Core ILP optimization for dense pattern detection.

## Configuration Parameters

Key parameters (modify in source code):

```python
min_col_quality = 3      # Minimum columns for valid regions
min_row_quality = 5      # Minimum rows for valid clusters  
error_rate = 0.025       # Noise tolerance (2.5%)
```

## Examples

### Single Matrix Analysis

```python
from pipeline import run_pipeline
import numpy as np

# Load your matrix
matrix = np.loadtxt('matrix.csv', delimiter=',', dtype=int)

# Run analysis
results = run_pipeline(
    matrix, 
    min_col_quality=3,
    min_row_quality=5, 
    error_rate=0.025
)

print(f"ILP calls: {results['ilp_calls_total']}")
print(f"Patterns found: {results['patterns_found']}")
```

### Batch Processing

```python
from run_experiments import run_experiments

# Process all matrices
results_file = run_experiments(
    matrix_dir="../data/imputed_matrices",
    output_dir="../data/experiments" 
)
```

## Expected Output Example

```
Starting ILP Haplotype Analysis Experiments
==================================================
Matrix directory: ../data/imputed_matrices
Output directory: ../data/experiments

Found matrix directory: /path/to/imputed_matrices
Haplotype subdirectories: ['1', '2', '3', '4', '5']
Total CSV files: 150

Processing 1/150: 1/matrix_001.csv
  -> Completed: 0 ILP calls, 5 patterns, 0.123s

...

=== EXPERIMENT SUMMARY ===
Total matrices processed: 150
Efficient runs (no ILP needed): 89
Complex runs (ILP required): 61
Algorithm efficiency rate: 59.3%

Haplotype counts: [1, 2, 3, 4, 5]
Matrix sizes range: 1000 - 50000
Total ILP calls across all runs: 234
Total processing time: 45.67s
Average ILP calls per matrix: 1.6
Average patterns found per matrix: 8.3

Results saved to: /path/to/experiments/ilp_call_stats.csv
```

## Troubleshooting

### Common Issues

1. **Gurobi License Error**
   ```
   Error: No Gurobi license found
   ```
   - Ensure Gurobi is properly licensed
   - Check WLS credentials in clustering.py

2. **No CSV Files Found**
   ```
   Error: No CSV files found in numbered subdirectories
   ```
   - Verify directory structure (numbered subdirectories)
   - Check file extensions (.csv)

3. **Memory Issues**
   - Reduce matrix size or adjust Gurobi memory settings
   - Process matrices in smaller batches

4. **Import Errors**
   ```
   ImportError: No module named 'gurobipy'
   ```
   - Install Gurobi: `pip install gurobipy`
   - Verify license activation

### Performance Tips

- Use SSD storage for large matrix datasets
- Increase `TimeLimit` in clustering.py for complex matrices
- Adjust `error_rate` based on data quality
- Monitor memory usage for very large matrices

## Algorithm Details

The package implements a novel haplotype reconstruction algorithm:

1. **Hierarchical Preprocessing**: Uses FeatureAgglomeration to identify coherent column regions
2. **Iterative Biclustering**: Alternates between positive and negative pattern detection
3. **ILP Optimization**: Uses Gurobi to solve quasi-biclique detection problems
4. **Pattern Accumulation**: Combines multiple patterns to reconstruct haplotype structure

## Citation

If you use this package in research, please cite:
[Add appropriate citation information]

## License

[Add license information]

## Support

For issues and questions:
- Check troubleshooting section above
- Review log output for specific error details
- Ensure input data format matches requirements
