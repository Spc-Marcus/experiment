# KNN Imputation for Genomic Matrices

This tool performs K-Nearest Neighbors (KNN) imputation on genomic variant matrices to fill missing values and binarize uncertain values.

## Features

- **KNN Imputation**: Fill missing values using K-nearest neighbors algorithm
- **Binarization**: Convert uncertain values to binary (0/1) based on confidence thresholds
- **Batch Processing**: Process single files or entire directory structures
- **Experiment Mode**: Test multiple K values to find optimal parameters
- **Dual Output**: Can perform both matrix imputation and experiments simultaneously

## Files

- `matrix.py`: Core imputation functions
- `impute_matrices.py`: Main script for processing matrices
- `Knn/`: Directory containing the KNN implementation

## Data Structure

```
data/
├── matrices_no_binarize/   
│   ├── subdirectory1/
│   │   ├── matrix1.csv
│   │   └── matrix2.csv
│   └── subdirectory2/
│       └── matrix3.csv
├── imputed_matrices/    # KNN imputed matrices
│   ├── subdirectory1/
│   │   ├── matrix1.csv
│   │   └── matrix2.csv
│   └── subdirectory2/
│       └── matrix3.csv
└── experiments/         # Experiment results
    └── KNN_exp.csv
```

## Usage Examples

### 1. Process a Single Matrix File

```bash
# Only matrix imputation with K=10
python Knn/impute_matrices.py data/matrix1.csv --matrix-output results/imputed/

# Only experiment mode (K values 5-25)
python Knn/impute_matrices.py data/matrix1.csv --experiment --exp-output results/experiments/

# Both matrix imputation AND experiment
python Knn/impute_matrices.py data/matrix1.csv --k 15 --matrix-output results/imputed/ --experiment --exp-output results/experiments/
```

### 2. Process Directory Structure

```bash
# Only process matrices with K=10
python Knn/impute_matrices.py data/matrices_no_binarize --k 10 --matrix-output data/imputed_matrices/

# Only run experiments
python Knn/impute_matrices.py data/matrices_no_binarize --experiment --exp-output data/experiments/

# Both operations simultaneously
python Knn/impute_matrices.py data/matrices_no_binarize --k 15 --matrix-output data/imputed_matrices/ --experiment --exp-output data/experiments/
```

### 3. Complete Workflow Examples

#### Research Pipeline
```bash
# Step 1: Run experiments to find optimal K
python Knn/impute_matrices.py data/ --experiment --exp-output experiments/

# Step 2: Analyze KNN_exp.csv to determine best K value

# Step 3: Process all matrices with optimal K
python Knn/impute_matrices.py data/ --k [optimal_k] --matrix-output final_matrices/
```

## Parameters

- `input`: Input CSV file or directory (required)
- `--k`: Number of neighbors for matrix imputation (default: 10)
- `--matrix-output`: Output directory for KNN imputed matrices
- `--exp-output`: Output directory for experiment results
- `--experiment`: Enable experiment mode (tests K values 5-25)

**Requirements:**
- At least one of `--matrix-output` or `--experiment` must be specified
- If `--experiment` is used, `--exp-output` is required
- Both operations can run simultaneously if both outputs are specified

## Input Matrix Format

CSV files with:
- **Rows**: Genomic reads
- **Columns**: Variant positions  
- **Values**: 
  - `1`: Reference allele
  - `0`: Alternative allele
  - `NaN/empty`: Missing values (to be imputed)
  - `0.0-1.0`: Probability scores

Example:
```csv
,pos_1,pos_2,pos_3,pos_4
read_1,1,0,,0.7
read_2,0,1,1,
read_3,,0,1,0.2
```

## Output Structure

### Matrix Imputation Output
- **Location**: Specified `--matrix-output` directory
- **Structure**: Maintains original subdirectory organization
- **Files**: 
  - Single files: `original_name_knn_imputed.csv`
  - Directories: Original filenames preserved in subdirectories
- **Content**: Binarized matrices (0/1 values only)

### Experiment Output
- **Location**: Specified `--exp-output` directory  
- **File**: Always named `KNN_exp.csv`
- **Content**: Results for K values 5-25

**KNN_exp.csv columns:**
- `sous_dossier`: Subdirectory name (directory processing)
- `fichier`: Original filename
- `k_value`: K parameter tested (5-25)
- `taille_matrice`: Matrix dimensions (reads×positions)
- `nan_avant`: Missing values before imputation
- `nan_après`: Uncertain values after binarization

## Algorithm Details

### KNN Imputation Process
1. **Load matrix**: Parse CSV with proper NaN handling
2. **Apply KNN**: Use scikit-learn's KNNImputer with specified K
3. **Binarize**: Apply 0.3 confidence threshold
   - Values ≤ 0.3 → 0
   - Values ≥ 0.7 → 1  
   - Values 0.3-0.7 → 0 (default)
4. **Save**: Preserve original structure and indexing

### Experiment Process
- Tests K values from 5 to 25
- For each K: performs full imputation + binarization
- Counts uncertain values (0.3 < value < 0.7) after binarization
- Aggregates results across all files and K values

## Performance Notes

- **Experiment mode**: Tests 21 different K values (5-25)
- **Directory limit**: Max 150 files per subdirectory in experiment mode
- **Memory**: Large matrices may require sufficient RAM for KNN computation
- **Time**: Experiment mode takes ~21x longer than single K processing

## Requirements

```bash
pip install pandas numpy scikit-learn
```
