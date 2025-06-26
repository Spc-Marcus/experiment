# BAM File Processing Pipeline for Variant Matrix Generation

This pipeline processes BAM (Binary Alignment Map) files to identify genomic variants and generate CSV matrices for downstream strain separation analysis. The tool is designed to detect heterozygous positions in genomic sequences and create variant matrices suitable for clustering and phylogenetic analysis.

## Overview

The pipeline consists of three main components:

1. **`core.py`**: Extracts suspicious genomic positions with variant alleles from BAM files
2. **`matrix.py`**: Converts variant data into filtered matrices and saves them as CSV files
3. **`process_bam_folder.py`**: Main script that orchestrates the processing of multiple BAM files

## What the Program Does

### Variant Detection Process

1. **Pileup Analysis**: For each genomic position, the program analyzes all reads covering that position
2. **Allele Frequency Calculation**: Identifies positions where multiple alleles are present with sufficient frequency
3. **Variant Calling**: Flags positions where the major allele frequency is below 95% (configurable)
4. **Binary Matrix Generation**: Creates matrices where:
   - `1` = read carries the major (reference-like) allele
   - `0` = read carries the minor (variant) allele
   - Missing values = read doesn't cover the position or was filtered out

### Filtering Criteria

- **Minimum Coverage**: Positions must be covered by at least 5 reads (default)
- **Quality Filtering**: Only bases with quality scores ≥ 10 are considered
- **Window Spanning**: Reads must have sufficient coverage in both the beginning and end thirds of variant positions (controlled by `--cover`)
- **Position Coverage**: Positions must have sufficient overall coverage across reads (controlled by `--threshold`)
- **Matrix Size**: Only matrices with ≥`--rows` reads and ≥`--cols` variant positions are saved (default: 20x20)

## Requirements

### Dependencies
```bash
pip install pandas numpy scikit-learn pysam
```

### Input Files
- **BAM files**: Properly mapped and sorted alignment files
- **BAI index files**: Each BAM file must have a corresponding `.bai` index file
  ```bash
  samtools index your_file.bam
  ```

## Usage Examples

### Basic Usage

From the repository root (`/udd/mfoin/Dev/refactor/experiment/`):

```bash
# Process all BAM files in a directory
python create_csv/process_bam_folder.py /path/to/bam_files /path/to/output

# Example with actual paths
python create_csv/process_bam_folder.py ./data/bam ./data
```

### Advanced Usage with Parameters

```bash
# Custom window size and coverage threshold
python create_csv/process_bam_folder.py \
    ./data/bam_samples \
    . \
    --window-size 10000 \
    --threshold 0.7

# Process specific sample directory
python create_csv/process_bam_folder.py \
    /XXXX/datasets/illumina_samples \
    /XXXX/results/variant_analysis \
    --window-size 5000 \
    --threshold 0.6
```

## Command Line Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bam_folder` | Path | Required | Directory containing BAM files |
| `output_folder` | Path | Required | Output directory for results |
| `--window-size` | int | 5000 | Size of genomic windows in base pairs |
| `--overlap` | int | 10000 | Overlap between consecutive windows |
| `--threshold` | float | 0.6 | Minimum read coverage threshold (0.0-1.0) |
| `--cover` | float | 0.6 | Read spanning coverage threshold for beginning and end positions (0.0-1.0) |
| `--rows` | int | 20 | Minimum number of rows (reads) required for saving matrix |
| `--cols` | int | 20 | Minimum number of columns (positions) required for saving matrix |

### Coverage Parameters Explained

- **`--threshold`**: Controls column (position) filtering. Positions with coverage below this threshold across all reads are removed from the matrix.

- **`--cover`**: Controls read (row) filtering based on spanning coverage. Reads must have data in both the beginning third AND end third of variant positions at this coverage level. Gaps in the middle positions are ignored.

**Example**: With `--cover 0.6`:
- If there are 30 variant positions, reads must cover at least 60% of the first 10 positions AND 60% of the last 10 positions
- Coverage in the middle 10 positions doesn't affect this filtering

### Matrix Size Parameters

- **`--rows`**: Minimum number of reads (rows) required for a matrix to be saved as CSV. Smaller matrices are discarded.

- **`--cols`**: Minimum number of variant positions (columns) required for a matrix to be saved as CSV. Matrices with fewer variants are discarded.

**Example**: With `--rows 30 --cols 25`:
- Only matrices with at least 30 reads AND at least 25 variant positions will be saved
- This helps ensure matrices have sufficient data for downstream analysis

## Real-World Examples

### Example 1: Processing Illumina Sequencing Data
```bash
# From repository root
python create_csv/process_bam_folder.py \
    /data/illumina_project/aligned_reads \
    /results/strain_analysis \
    --threshold 0.8
```

### Example 2: High-Resolution Analysis with Strict Coverage
```bash
# Smaller windows with strict read spanning requirements
python create_csv/process_bam_folder.py \
    ./samples/pacbio_data \
    ./output/high_res \
    --window-size 2000 \
    --threshold 0.7 \
    --cover 0.8 \
    --rows 30 \
    --cols 25
```

### Example 3: Lenient Coverage for Sparse Data
```bash
# More permissive settings for low-coverage samples
python create_csv/process_bam_folder.py \
    ~/datasets/viral_samples \
    ~/analysis_results/variants \
    --threshold 0.5 \
    --cover 0.4 \
    --rows 10 \
    --cols 15
```

### Example 4: Complete Parameter Specification
```bash
# Full parameter control
python create_csv/process_bam_folder.py \
    /data/bam_files \
    /results/matrices \
    --window-size 10000 \
    --threshold 0.6 \
    --cover 0.7 \
    --rows 25 \
    --cols 30
```

### Example 5: High-Quality Analysis
```bash
# Strict requirements for high-quality matrices
python create_csv/process_bam_folder.py \
    /data/high_coverage_samples \
    /results/strict_analysis \
    --threshold 0.8 \
    --cover 0.8 \
    --rows 50 \
    --cols 40
```

## Expected Folder Structure

### Input Structure
```
project_root/
├── data/
│   └── bam_samples/
│       ├── sample1.bam
│       ├── sample1.bam.bai
│       ├── sample2.bam
│       ├── sample2.bam.bai
│       ├── sample3.bam
│       └── sample3.bam.bai
└── create_csv/
    ├── process_bam_folder.py
    ├── core.py
    ├── matrix.py
    └── README.md
```

### Output Structure
```
matrices_no_binarize/
    ├── sample1/
    │   ├── matrix_chr1_0.csv
    │   ├── matrix_chr1_5000.csv
    │   ├── matrix_chr2_0.csv
    │   └── matrix_chr2_5000.csv
    ├── sample2/
    │   ├── matrix_chr1_0.csv
    │   ├── matrix_chr1_5000.csv
    │   └── matrix_chr3_0.csv
    └── sample3/
        ├── matrix_chr1_0.csv
        └── matrix_chr2_0.csv
```

## Output Files

### CSV Matrix Format
Each CSV file contains:
- **Rows**: Read names (unique identifiers for sequencing reads)
- **Columns**: Genomic positions where variants were detected
- **Values**: 
  - `1`: Read matches the major allele (reference-like)
  - `0`: Read matches the minor allele (variant)
  - Empty cells: Read doesn't cover this position

### Example CSV Content
```csv
,1000,1001,1002,1003
read_001,1,0,1,1
read_002,1,0,1,0
read_003,0,1,1,1
read_004,1,0,0,1
```

## Performance Notes

- **Large Contigs**: Automatically uses 5000 bp windows for contigs >90,000 bp
- **Small Contigs**: Contigs ≤90,000 bp are skipped entirely
- **Memory Usage**: Processes one window at a time to manage memory efficiently
- **Processing Time**: Depends on BAM file size and number of variants

## Troubleshooting

### Common Issues

1. **Missing Index Files**
   ```bash
   Error: Index file sample.bam.bai not found
   Solution: samtools index sample.bam
   ```

2. **Permission Errors**
   ```bash
   Error: Permission denied to create/write to /results
   Solution: Use a directory in your home folder or check permissions
   ```

3. **No BAM Files Found**
   ```bash
   No BAM files found in /path/to/folder
   Solution: Check path and ensure .bam files are present
   ```

4. **No Matrices Saved**
   ```bash
   Problem: Matrices are created but no CSV files are saved
   Possible causes:
   - Matrix dimensions below --rows/--cols thresholds
   - Contig size ≤90,000 bp (small contigs are skipped)
   - Insufficient variant positions after filtering
   Solution: Lower --rows and --cols values, or adjust --threshold and --cover
   ```

5. **Parameter Validation Errors**
   ```bash
   Error: Minimum rows/columns must be positive
   Solution: Ensure --rows and --cols are > 0
   
   Error: Cover threshold must be between 0.0 and 1.0
   Solution: Use decimal values like 0.6, not percentages like 60
   ```

