# Experiment Results Analysis

This folder contains analysis tools for experimentation results of different imputation algorithms.

## Available Analysis Tools

- [`KNN.py`](#knn-analysis) - K-Nearest Neighbors imputation analysis
- [`ILP_call.py`](#ilp-call-analysis) - ILP optimization usage analysis
- [`Efficience.py`](#efficiency-analysis) - Efficiency heatmap analysis
- [`efficience_ilp.py`](#ilp-efficiency-analysis) - ILP vs Pre-ILP efficiency analysis
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

## Efficiency Analysis

### Description

The `Efficience.py` script analyzes clustering efficiency by evaluating how often the number of final clusters matches the number of haplotypes across different certainty thresholds and error rates. It produces a heatmap of success rates per haplotype count.

### Prerequisites

```bash
pip install pandas matplotlib seaborn numpy
```

### Usage

#### Syntax
```bash
python Efficience.py
```

Note: This script currently uses a default CSV path defined inside the file. To use a custom CSV or output directory, call the analyzer programmatically (see examples).

#### Options
- No CLI options. Edit `csv_file` in `main()` or use the programmatic examples below.

#### Examples

**Basic analysis (uses the default CSV path inside the script):**
```bash
python Efficience.py
```

**Programmatic usage with a custom CSV and default output path:**
```bash
python -c "from Efficience import EfficiencyAnalyzer; a=EfficiencyAnalyzer('./results.csv'); a.generate_all_plots()"
```

**Programmatic usage with a custom CSV and custom output directory:**
```bash
python -c "from Efficience import EfficiencyAnalyzer; a=EfficiencyAnalyzer('./results.csv'); a.generate_all_plots('./plots')"
```

### Input Format

The CSV file should contain:
- `Threshold` : Certainty threshold
- `Error-Rate` : Error rate used
- `Strip` : Strip size or related parameter
- `Haplotype` : Expected number of haplotypes
- `Steps-Count` : Number of steps
- `Unused-Cols` : Count of unused columns
- `Final-Clusters` : Number of clusters produced
- `Orphan-Reads` : Number of orphan reads
- `Time-Pre-Processing` : Pre-processing time
- `Time-Post-Processing` : Post-processing time
- `Matrix-Size` : Tuple string "(rows, cols)" (unquoted tuples are auto-fixed)

### Output Structure

```
[output_directory]/                         # Default: /home/mafoin/stage/experiment/analyze/plots
└── heatmap_certitude_error_by_haplotype.png
```

### Generated Analyses

1. Success-rate heatmap per haplotype count across Threshold × Error-Rate
2. Console summary (dataset size, ranges, and unique values)

### Calculated Metrics

- Success rate: `1` if `Final-Clusters == Haplotype`, else `0`; averaged per (Threshold, Error-Rate) grid.

### Key Insights

- Identifies parameter regions (certainty/error rate) yielding perfect clustering.
- Highlights sensitivity of results to error and certainty settings across haplotype counts.

## ILP Efficiency Analysis

### Description

The `efficience_ilp.py` script compares clustering efficiency before and after ILP. It computes success metrics (clusters match haplotypes), orphan/steps reductions, and generates a comparison plot with optional parameter filters. It also exports per-parameter sweep plots.

### Prerequisites

```bash
pip install pandas numpy matplotlib
```

### Usage

#### Syntax
```bash
python efficience_ilp.py
```

Note: This script uses a default CSV path defined inside the file. To use a custom CSV or apply parameter filters, call the analyzer programmatically.

#### Options
- No CLI options. Edit `csv_file` in `main()` or use programmatic calls.

#### Examples

**Basic analysis (default CSV path):**
```bash
python efficience_ilp.py
```

**Programmatic usage (custom CSV and default outputs):**
```bash
python -c "from efficience_ilp import ILPEfficiencyAnalyzer; a=ILPEfficiencyAnalyzer('./results.csv'); a.generate_all_plots()"
```

**Programmatic usage with filters (saved to a custom path):**
```bash
python -c "from efficience_ilp import ILPEfficiencyAnalyzer; a=ILPEfficiencyAnalyzer('./results.csv'); a.plot_success_rate_comparison('./plots/success_filtered.png', threshold=0.9, error_rate=0.02, distance=0.1, show=False)"
```

### Input Format

The CSV file should contain:
- `Threshold`
- `Error-Rate`
- `Strip`
- `Haplotype`
- `Matrix-Size-cols`
- `Matrix-Size-rows`
- `Time-Pre-Processing`
- `Steps-Count-Pre`
- `Unused-Cols-Pre`
- `Final-Clusters-Pre`
- `Orphan-Reads-Pre`
- `Steps-Count-ILP`
- `Final-Clusters-ILP`
- `Orphan-Reads-ILP`
- `Distance`

### Output Structure

```
[output_directory]/                      # Default: /home/mafoin/stage/experiment/analyze/plots_ilp
├── success_rate_comparison.png
└── success_rate_by_param/
    ├── threshold/threshold_<val>.png
    ├── error_rate/error_rate_<val>.png
    └── distance/distance_<val>.png
```

### Generated Analyses

1. Success rate comparison (Pre-ILP vs Post-ILP) by haplotype count
2. Exported sweeps for Threshold, Error-Rate, and Distance
3. Summary statistics printed to console

### Calculated Metrics

- `Success_Rate_Pre`, `Success_Rate_ILP`
- `Clusters_Improvement` = `Final-Clusters-ILP - Final-Clusters-Pre`
- `Orphan_Reduction` and `Orphan_Reduction_Rate`
- `Steps_Reduction`

### Key Insights

- Quantifies ILP’s contribution to success rate and orphan reduction
- Highlights sensitivity to Threshold, Error-Rate, and Distance parameters


[🔝 Back to top](#experiment-results-analysis)

