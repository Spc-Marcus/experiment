# Density Clustering Analysis Package

This package provides tools for analyzing binary matrices using density-based clustering with Integer Linear Programming (ILP) optimization. It specializes in detecting quasi-bicliques and pattern separation in dense binary matrices.

## Overview

The package implements a density-focused clustering pipeline:
1. **Matrix Loading**: Supports CSV files with binary data (0s and 1s)
2. **Quasi-Biclique Detection**: Uses ILP optimization to find dense patterns
3. **Iterative Clustering**: Alternates between positive and negative pattern detection
4. **Comprehensive Logging**: Tracks all intermediate results and quasi-biclique discoveries

## Installation

### Prerequisites

- Python 3.8+
- Gurobi Optimizer (license required)
- Required Python packages:
  ```bash
  pip install numpy pandas gurobipy pathlib
  ```

### Gurobi Setup

This package requires Gurobi with valid licensing. The code includes WLS credentials for academic use.

## Usage

### Basic Matrix Processing

```bash
python matrix_processor.py matrix_edge_17_10000.csv --output ../results/clustering_results.json
```

### Command Line Options

```bash
python matrix_processor.py <csv_file> [OPTIONS]

Options:
  --output FILE          Output JSON file (default: clustering_results_base.json)
  --error-rate FLOAT     Error rate for clustering (default: 0.025)
  --min-row-quality INT  Minimum row quality (default: 5)
  --min-col-quality INT  Minimum column quality (default: 3)
```

### Example Commands

```bash
# Basic processing with default parameters
python matrix_processor.py matrix_edge_17_10000.csv

# Custom output file and parameters
python matrix_processor.py matrix_edge_17_10000.csv \
  --output ../results/my_results.json \
  --error-rate 0.05 \
  --min-row-quality 3 \
  --min-col-quality 2

# High sensitivity analysis
python matrix_processor.py matrix_edge_17_10000.csv \
  --output ../results/high_sens_results.json \
  --error-rate 0.01 \
  --min-row-quality 3
```

## Input Format

### CSV Matrix Requirements

Input matrices should be CSV files with:
- First row: column headers (C0, C1, C2, ...)
- First column: row indices (L0, L1, L2, ...)
- Values: binary (0 or 1) only
- Missing values automatically converted to 0

### Example Input Matrix

```csv
,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16
L0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
L1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
L2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
L3,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1
L4,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1
L5,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1
```

## Output Results

### JSON Results Structure

The analysis generates a comprehensive JSON file containing:

```json
{
  "timestamp": "2025-06-27 13:48:11",
  "matrix_info": {
    "shape": [50, 26],
    "total_elements": 1300,
    "ones_count": 1130,
    "zeros_count": 170,
    "density": 0.8692
  },
  "quasi_biclique_results": [...],
  "clustering_steps": [...],
  "final_results": {...}
}
```

### Matrix Information

Basic statistics about the input matrix:
- **shape**: Matrix dimensions [rows, columns]
- **density**: Fraction of 1s in the matrix
- **ones_count/zeros_count**: Distribution of values

### Quasi-Biclique Results

Detailed information about each quasi-biclique detected:

```json
{
  "phase": "iteration_1_positive",
  "input_shape": [50, 26],
  "selected_rows": [0, 1, 2, 3, ...],
  "selected_cols": [4, 5, 6, 7, ...],
  "success": true,
  "selected_density": 0.9824,
  "zeros_positions": [[7, 3], [10, 5], ...],
  "zeros_positions_count": 16
}
```

### Clustering Steps

Summary of each clustering iteration:

```json
{
  "step_number": 1,
  "reads_group1": [0, 1, 2, ...],
  "reads_group0": [37, 38, 39, ...],
  "columns": [0, 1, 4, 5, ...],
  "ilp_calls": 2,
  "group1_density": {
    "ones_density": 0.9579,
    "zeros_density": 0.0421
  },
  "group0_density": {
    "ones_density": 0.0,
    "zeros_density": 1.0
  },
  "density_contrast": {
    "ones_difference": 0.9579,
    "zeros_difference": -0.9579
  }
}
```

## Examples

### Example 1: Processing Dense Matrix

Given the provided matrix `matrix_edge_17_10000.csv`:

```bash
python matrix_processor.py matrix_edge_17_10000.csv --output dense_analysis.json
```

**Expected Output:**
```
Loaded matrix from matrix_edge_17_10000.csv
Shape: (50, 26)
Density: 0.869

Processing matrix of shape (50, 26)
Using 1 regions
Error rate: 0.025
Min row quality: 5
Min column quality: 3

Processing completed in 2.45 seconds
Found 1 valid clustering steps
Results saved to dense_analysis.json

Final results:
  Step 1: 38 vs 12 reads, 10 columns
```

### Example 2: High Sensitivity Analysis

For detecting more subtle patterns:

```bash
python matrix_processor.py matrix_edge_17_10000.csv \
  --output sensitive_analysis.json \
  --error-rate 0.01 \
  --min-row-quality 3 \
  --min-col-quality 2
```

### Example 3: Programmatic Usage

```python
from matrix_processor import load_csv_matrix, process_matrix_with_logging
import numpy as np

# Load matrix
matrix = load_csv_matrix("matrix_edge_17_10000.csv")

# Process with custom parameters
results = process_matrix_with_logging(
    matrix,
    output_file="../results/my_analysis.json",
    error_rate=0.025,
    min_row_quality=5,
    min_col_quality=3
)

print(f"Found {len(results)} clustering steps")
for i, (group1, group0, cols) in enumerate(results):
    print(f"Step {i+1}: {len(group1)} vs {len(group0)} reads, {len(cols)} columns")
```

### Example 4: Using the Clustering Logger

```python
from matrix_clustering_logger import process_matrix_with_logging
import numpy as np

# Create example matrix
matrix = np.random.choice([0, 1], size=(30, 20), p=[0.3, 0.7])

# Process with comprehensive logging
results = process_matrix_with_logging(
    matrix,
    output_file="../results/random_matrix_analysis.json",
    error_rate=0.05,
    min_row_quality=4
)
```

## Results Storage

All results are stored in the `../results/` directory:

```
../results/
├── clustering_results.json      # Default output
├── dense_analysis.json         # Example 1 output
├── sensitive_analysis.json     # Example 2 output
├── my_analysis.json           # Example 3 output
└── random_matrix_analysis.json # Example 4 output
```

## Configuration Parameters

### Key Parameters

- **error_rate** (default: 0.025): Noise tolerance for quasi-biclique detection
  - Lower values = more stringent pattern detection
  - Higher values = more tolerant of noise

- **min_row_quality** (default: 5): Minimum rows required for valid clusters
  - Higher values = larger, more stable patterns
  - Lower values = detect smaller patterns

- **min_col_quality** (default: 3): Minimum columns for valid regions
  - Controls granularity of pattern detection

### Optimization Settings

The ILP optimization uses these Gurobi parameters:
- **MIPGAP**: 0.05 (5% optimality gap)
- **TimeLimit**: 200 seconds per optimization
- **OutputFlag**: 0 (suppressed output)

## Algorithm Details

### Quasi-Biclique Detection

The core algorithm uses a three-phase approach:

1. **Seeding**: Find dense sub-regions with >99% density
2. **Row Extension**: Add compatible rows (>50% compatibility)
3. **Column Extension**: Add compatible columns (>90% compatibility)

### Mathematical Formulation

**Objective**: Maximize Σ M[i,j] × x[i,j] (sum of 1s in selected cells)

**Constraints**:
- x[i,j] ≤ row[i] (cell selected only if row selected)
- x[i,j] ≤ col[j] (cell selected only if column selected)
- Σ x[i,j] × (1-M[i,j]) ≤ error_rate × Σ x[i,j] (density constraint)

### Iterative Clustering

The algorithm alternates between:
- **Positive patterns**: Dense regions of 1s
- **Negative patterns**: Dense regions of 0s (inverted matrix)


### Debug Mode

Enable detailed logging by modifying clustering.py:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Advanced Usage

### Custom Region Processing

```python
from clustering import clustering_full_matrix
import numpy as np

matrix = np.loadtxt("matrix.csv", delimiter=",", dtype=int)

# Define custom column regions
regions = [
    [0, 1, 2, 3, 4],      # First region: columns 0-4
    [5, 6, 7, 8, 9],      # Second region: columns 5-9
    [10, 11, 12, 13, 14]  # Third region: columns 10-14
]

# Process with custom regions
results = clustering_full_matrix(
    matrix,
    regions=regions,
    steps=[],
    error_rate=0.025,
    log_results=True
)
```

### Batch Processing Multiple Files

```python
import os
from pathlib import Path
from matrix_processor import process_matrix_with_logging, load_csv_matrix

input_dir = "input_matrices/"
output_dir = "../results/"

for csv_file in Path(input_dir).glob("*.csv"):
    print(f"Processing {csv_file.name}")
    
    matrix = load_csv_matrix(str(csv_file))
    output_file = f"{output_dir}/{csv_file.stem}_results.json"
    
    results = process_matrix_with_logging(
        matrix,
        output_file=output_file,
        error_rate=0.025
    )
    
    print(f"  -> {len(results)} patterns found")
```

