# Experiment Results Analysis

This folder contains analysis tools for experimentation results of different imputation algorithms.

## Available Analysis Tools

- [`KNN.py`](#knn-analysis) - K-Nearest Neighbors imputation analysis
- [`ILP_call.py`](#ilp-call-analysis) - ILP optimization usage analysis
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

## ILP Call Analysis

### Description

The `ILP_call.py` script analyzes Integer Linear Programming (ILP) optimization usage in haplotype analysis experiments. It evaluates algorithm efficiency by examining when ILP optimization is required versus when problems can be solved efficiently without it.

### Prerequisites

```bash
pip install pandas matplotlib seaborn numpy scipy
```

### Usage

#### Syntax
```bash
python ILP_call.py [CSV_file] [options]
```

#### Options
- `--plots` : Generate visualization plots
- `--output`, `-o` : Output directory (default: ./results)
- `--subdir` : Subdirectory name for ILP analysis results (default: ilp_call)

#### Examples

**Basic analysis:**
```bash
python ILP_call.py ../data/experiments/ilp_call_stats.csv
```

**Complete analysis with plots:**
```bash
python ILP_call.py ./results/experiment_results.csv --plots
```

**Custom output directory:**
```bash
python ILP_call.py /path/to/results.csv --plots --output ./analysis_plots
```

**Custom subdirectory name:**
```bash
python ILP_call.py data.csv --plots --subdir my_ilp_analysis
```

**Complete custom paths:**
```bash
python ILP_call.py data.csv --plots --output ./custom_results --subdir experiment_1
```

### Input Format

The CSV file must contain:
- `haplotype_count` : Number of haplotypes
- `ilp_calls_total` : Total ILP optimization calls
- `matrix_size` : Matrix size (rows × cols)
- `total_time` : Total execution time
- `matrix_rows`, `matrix_cols` : Matrix dimensions
- `matrix_density` : Matrix density (optional)
- `patterns_found` : Number of patterns found (optional)

### Output Structure

```
[output_directory]/
└── [subdir_name]/                              # Default: ilp_call/
    ├── analysis_summary.txt                    # Complete analysis summary
    ├── time_vs_size_all_data.png              # Execution time vs matrix size (if --plots)
    ├── complex_vs_efficient_comparison.png    # ILP vs efficient runs comparison (if --plots)
    ├── ilp_calls_distribution.png             # ILP usage distribution (if --plots)
    └── ilp_utilization_analysis.png           # ILP utilization patterns (if --plots)
```

### Generated Analyses

1. **Performance by haplotype count**
2. **ILP utilization rates**
3. **Scalability analysis**
4. **Algorithm efficiency patterns**
5. **Statistical comparisons**
6. **Optimal characteristics identification**

### Calculated Metrics

- **ILP utilization rate** : `(runs_with_ILP / total_runs) * 100`
- **Matrix complexity** : `matrix_size × matrix_density`
- **Algorithm efficiency** : Runs solved without ILP optimization
- **Performance correlations** : Matrix characteristics vs ILP requirements

### Key Insights

- **Efficiency classification** : Distinguishes between runs requiring ILP optimization vs efficient runs
- **Scalability patterns** : Identifies matrix characteristics that lead to complex optimization requirements
- **Sweet spot identification** : Determines optimal matrix ranges for efficient processing


[🔝 Back to top](#experiment-results-analysis)

