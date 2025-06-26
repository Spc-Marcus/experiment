# Experiment Results Analysis

This folder contains analysis tools for experimentation results of different imputation algorithms.

## Available Analysis Tools

- [`KNN.py`](#knn-analysis) - K-Nearest Neighbors imputation analysis
- `README.md` - This documentation file

## KNN Analysis

### Description

The `KNN.py` script analyzes K-Nearest Neighbors algorithm experimentation results for missing value imputation in phylogenetic matrices.

### Prerequisites

```bash
pip install pandas matplotlib seaborn numpy
```

### Usage

#### Syntax
```bash
python KNN.py [CSV_file] [options]
```

#### Options
- `--plots` : Generate visualization plots
- `--output`, `-o` : Output directory (default: ./results)

#### Examples

**Basic analysis:**
```bash
python KNN.py ../data/experiments/KNN_exp.csv
```

**Complete analysis with plots:**
```bash
python KNN.py ../data/experiments/KNN_exp.csv --plots -o ../results
```

### Input Format

The CSV file must contain:
- `sous_dossier` : Number of haplotypes
- `fichier` : Matrix filename
- `k_value` : K value tested
- `taille_matrice` : Dimensions (format "NxM")
- `nan_avant` : Missing values before imputation
- `nan_après` : Missing values after imputation

### Output Structure

```
[output_directory]/
└── KNN/
    ├── analysis_summary.txt           # Complete summary
    ├── correlation_analysis.png       # Correlation analysis (if --plots)
    └── performance_analysis.png       # Performance analysis (if --plots)
```

### Generated Analyses

1. **Performance by K value**
2. **Perfect imputation statistics**
3. **Key correlations**
4. **Performance by haplotype count**
5. **Optimal K distribution**
6. **Executive summary**

### Calculated Metrics

- **Reduction rate** : `((nan_before - nan_after) / nan_before) * 100`
- **NaN density** : `(nan / total_size) * 100`


[🔝 Back to top](#experiment-results-analysis)

