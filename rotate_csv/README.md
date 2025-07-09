# CSV Matrix Rotation Tool

This tool processes CSV matrices by transposing (rotating) them, effectively swapping rows and columns. This is useful for data analysis workflows where you need to change the orientation of your data matrix for different types of analysis.

## Overview

The rotation tool consists of two main components:

1. **`matrix_rotator.py`**: Core functionality for transposing CSV matrices
2. **`process_csv_folder.py`**: Main script that orchestrates the processing of multiple CSV files

## What the Program Does

### Matrix Transposition Process

1. **CSV Reading**: Loads CSV matrices using pandas for efficient data handling
2. **Transposition**: Swaps rows and columns using pandas transpose functionality
3. **Data Preservation**: Maintains all original data values and structure
4. **Output Generation**: Saves transposed matrices as new CSV files

### Transformation Example

**Original Matrix:**
```csv
,pos1,pos2,pos3,pos4
read_001,1,0,1,1
read_002,1,0,1,0
read_003,0,1,1,1
read_004,1,0,0,1
```

**After Transposition:**
```csv
,read_001,read_002,read_003,read_004
pos1,1,1,0,1
pos2,0,0,1,0
pos3,1,1,1,0
pos4,1,0,1,1
```

## Requirements

### Dependencies
```bash
pip install pandas numpy
```

### Input Files
- **CSV files**: Properly formatted CSV matrices with headers
- **File structure**: CSV files should have row and column headers

## Usage Examples

### Basic Usage

From the repository root (`/udd/mfoin/Dev/refactor/experiment/`):

```bash
# Rotate all CSV files in a directory
python rotate_csv/process_csv_folder.py /path/to/csv_files /path/to/output

# Example with actual paths
python rotate_csv/process_csv_folder.py ./data/matrices ./data/rotated
```

### Advanced Usage with Parameters

```bash
# Process specific directory with custom settings
python rotate_csv/process_csv_folder.py \
    ./data/matrices_no_binarize \
    ./data/rotated_matrices \
    --suffix "_transposed"

# Process with verbose output
python rotate_csv/process_csv_folder.py \
    /XXXX/datasets/matrices \
    /XXXX/results/rotated \
    --verbose
```

## Command Line Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_folder` | Path | Required | Directory containing CSV files |
| `output_folder` | Path | Required | Output directory for rotated matrices |
| `--suffix` | str | "_rotated" | Suffix to add to output filenames |
| `--verbose` | flag | False | Enable verbose output with progress information |
| `--recursive` | flag | False | Process subdirectories recursively |

### Output Naming Convention

- **Input file**: `matrix_chr1_0.csv`
- **Output file**: `matrix_chr1_0_rotated.csv` (with default suffix)

## Real-World Examples

### Example 1: Basic Matrix Rotation
```bash
# From repository root
python rotate_csv/process_csv_folder.py \
    ./data/matrices_no_binarize \
    ./data/rotated_matrices
```

### Example 2: Custom Suffix for Transposed Matrices
```bash
# Use "_transposed" suffix instead of default
python rotate_csv/process_csv_folder.py \
    ./data/variant_matrices \
    ./data/transposed_matrices \
    --suffix "_transposed"
```

### Example 3: Verbose Processing with Progress
```bash
# Show detailed progress information
python rotate_csv/process_csv_folder.py \
    /data/analysis_matrices \
    /results/rotated_analysis \
    --verbose
```

### Example 4: Recursive Processing
```bash
# Process all CSV files in subdirectories
python rotate_csv/process_csv_folder.py \
    ./data/nested_matrices \
    ./data/rotated_nested \
    --recursive \
    --verbose
```

### Example 5: Complete Parameter Specification
```bash
# Full parameter control
python rotate_csv/process_csv_folder.py \
    /data/bam_matrices \
    /results/rotated_bam \
    --suffix "_transposed" \
    --verbose \
    --recursive
```

## Expected Folder Structure

### Input Structure
```
project_root/
├── data/
│   └── matrices_no_binarize/
│       ├── sample1/
│       │   ├── matrix_chr1_0.csv
│       │   ├── matrix_chr1_5000.csv
│       │   └── matrix_chr2_0.csv
│       ├── sample2/
│       │   ├── matrix_chr1_0.csv
│       │   └── matrix_chr2_0.csv
│       └── sample3/
│           └── matrix_chr1_0.csv
└── rotate_csv/
    ├── process_csv_folder.py
    ├── matrix_rotator.py
    └── README.md
```

### Output Structure
```
rotated_matrices/
    ├── sample1/
    │   ├── matrix_chr1_0_rotated.csv
    │   ├── matrix_chr1_5000_rotated.csv
    │   └── matrix_chr2_0_rotated.csv
    ├── sample2/
    │   ├── matrix_chr1_0_rotated.csv
    │   └── matrix_chr2_0_rotated.csv
    └── sample3/
        └── matrix_chr1_0_rotated.csv
```

## Output Files

### Transposed CSV Format
Each rotated CSV file contains:
- **Rows**: Original column headers (e.g., genomic positions)
- **Columns**: Original row headers (e.g., read names)
- **Values**: Same data values, but with rows and columns swapped

### Example Transposed CSV Content
```csv
,read_001,read_002,read_003,read_004
1000,1,1,0,1
1001,0,0,1,0
1002,1,1,1,0
1003,1,0,1,1
```

## Use Cases

### 1. Analysis Perspective Change
- **Original**: Rows = reads, Columns = positions
- **Rotated**: Rows = positions, Columns = reads
- **Benefit**: Analyze position-wise patterns instead of read-wise patterns

### 2. Statistical Analysis
- **Original**: Each row represents a read's variant profile
- **Rotated**: Each row represents a position's read profile
- **Benefit**: Calculate position-based statistics and correlations

### 3. Visualization
- **Original**: Heatmaps with reads on y-axis, positions on x-axis
- **Rotated**: Heatmaps with positions on y-axis, reads on x-axis
- **Benefit**: Different visual perspectives for pattern identification

### 4. Machine Learning
- **Original**: Features = positions, Samples = reads
- **Rotated**: Features = reads, Samples = positions
- **Benefit**: Different feature space for clustering and classification

## Performance Notes

- **Memory Usage**: Processes one file at a time to manage memory efficiently
- **Processing Time**: Depends on CSV file size and number of files
- **File Size**: Transposed files maintain the same data volume
- **Compatibility**: Works with any CSV format that pandas can read

## Troubleshooting

### Common Issues

1. **Missing Input Directory**
   ```bash
   Error: Input directory /path/to/csv_files not found
   Solution: Check the path and ensure the directory exists
   ```

2. **Permission Errors**
   ```bash
   Error: Permission denied to create/write to /output
   Solution: Use a directory in your home folder or check permissions
   ```

3. **No CSV Files Found**
   ```bash
   No CSV files found in /path/to/folder
   Solution: Check path and ensure .csv files are present
   ```

4. **Empty CSV Files**
   ```bash
   Warning: File matrix.csv is empty or contains no data
   Solution: Check if the input CSV files contain valid data
   ```

5. **Memory Issues with Large Files**
   ```bash
   Error: Memory error when processing large CSV file
   Solution: Process files individually or increase system memory
   ```

6. **Invalid CSV Format**
   ```bash
   Error: Unable to parse CSV file matrix.csv
   Solution: Check CSV format and ensure proper headers
   ```

## Data Integrity

The rotation process preserves:
- **All data values**: No data loss during transposition
- **Data types**: Numeric values remain numeric
- **Missing values**: NaN or empty cells are preserved
- **File structure**: Maintains CSV format compatibility

## Integration with Other Tools

This rotation tool is designed to work seamlessly with:
- **create_csv**: Rotate matrices generated from BAM processing
- **ILP_call**: Prepare rotated matrices for integer linear programming
- **KNN**: Rotate matrices for k-nearest neighbors analysis
- **density_contre_exemple**: Rotate matrices for density-based clustering 