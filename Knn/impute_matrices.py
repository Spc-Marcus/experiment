import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
import glob
from matrix import impute_missing_values

def load_matrix(csv_path):
    """Load matrix from CSV file with proper indexing."""
    # Force empty strings and whitespace to be treated as NaN
    df = pd.read_csv(csv_path, index_col=0, na_values=['', ' ', 'NA', 'nan', 'NaN'])
    return df

def save_matrix(df, output_path):
    """Save matrix to CSV file preserving original format."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path)

def count_nan(matrix):
    """Count NaN values in matrix."""
    return np.isnan(matrix).sum()

def binarize_matrix(matrix: np.ndarray, certitude: float = 0.3, default: int = 0) -> tuple:
    """
    Binarize a matrix and count uncertain values.
    
    Returns:
    --------
    tuple: (binarized_matrix, uncertain_count)
    """
    if matrix.size == 0:
        return matrix.astype(int), 0
    
    # Calculate thresholds
    lower_threshold = certitude
    upper_threshold = 1.0 - certitude
    
    # Count uncertain values (between thresholds, exclusive)
    uncertain_mask = (matrix > lower_threshold) & (matrix < upper_threshold)
    uncertain_count = uncertain_mask.sum()
    
    # Apply binarization
    result = np.full_like(matrix, default, dtype=int)
    result[matrix <= lower_threshold] = 0
    result[matrix >= upper_threshold] = 1
    
    return result, uncertain_count

def process_matrix_normal(csv_path, output_path, k):
    """Process single matrix in normal mode."""
    df = load_matrix(csv_path)
    matrix = df.values
    
    # Apply KNN imputation
    imputed_matrix = impute_missing_values(matrix, n_neighbors=k)
    
    # Binariser la matrice imputée
    binarized_matrix, _ = binarize_matrix(imputed_matrix, certitude=0.3, default=0)
    
    # Create new dataframe with binarized values
    imputed_df = pd.DataFrame(binarized_matrix, index=df.index, columns=df.columns)
    
    # Save imputed matrix
    save_matrix(imputed_df, output_path)

def process_matrix_experiment(csv_path, k_values):
    """Process single matrix in experiment mode."""
    df = load_matrix(csv_path)
    matrix = df.values
    
    results = []
    for k in k_values:
        nan_before = count_nan(matrix)
        imputed_matrix = impute_missing_values(matrix, n_neighbors=k)
        
        # Binariser la matrice imputée et compter les valeurs incertaines
        binarized_matrix, uncertain_count = binarize_matrix(imputed_matrix, certitude=0.3, default=0)
        
        results.append({
            'k_value': k,
            'taille_matrice': f"{matrix.shape[0]}×{matrix.shape[1]}",
            'nan_avant': nan_before,
            'nan_après': uncertain_count  # Nombre de valeurs incertaines après binarisation
        })
    
    return results

def get_csv_files(input_dir, max_files=None):
    """Get all CSV files organized by subdirectory."""
    csv_files = {}
    
    for subdir in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, subdir)
        if os.path.isdir(subdir_path):
            pattern = os.path.join(subdir_path, "*.csv")
            files = glob.glob(pattern)
            
            if max_files:
                files = files[:max_files]
            
            csv_files[subdir] = files
    
    return csv_files

def main():
    parser = argparse.ArgumentParser(description="Impute missing values in genomic matrices using KNN")
    parser.add_argument("input_dir", help="Input directory containing subdirectories with CSV matrices")
    parser.add_argument("--k", type=int, default=10, help="KNN value (default: 10)")
    parser.add_argument("--experiment", action="store_true", help="Run in experiment mode")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory {args.input_dir} does not exist")
        return
    
    if args.experiment:
        print("Running in experiment mode...")
        k_values = list(range(5, 26))
        csv_files = get_csv_files(args.input_dir, max_files=150)
        
        all_results = []
        
        for subdir, files in csv_files.items():
            print(f"Processing subdirectory: {subdir}")
            
            for csv_path in files:
                filename = os.path.basename(csv_path)
                print(f"  Processing file: {filename}")
                
                try:
                    results = process_matrix_experiment(csv_path, k_values)
                    
                    for result in results:
                        result['sous_dossier'] = subdir
                        result['fichier'] = filename
                        all_results.append(result)
                        
                except Exception as e:
                    print(f"    Error processing {filename}: {e}")
        
        # Save experiment results
        results_df = pd.DataFrame(all_results)
        results_df = results_df[['sous_dossier', 'fichier', 'k_value', 'taille_matrice', 'nan_avant', 'nan_après']]
        results_df.to_csv("experiment_results.csv", index=False)
        print(f"Experiment results saved to experiment_results.csv")
        
    else:
        print(f"Running in normal mode with k={args.k}...")
        csv_files = get_csv_files(args.input_dir)
        
        output_base = "matrices"
        
        for subdir, files in csv_files.items():
            print(f"Processing subdirectory: {subdir}")
            
            for csv_path in files:
                filename = os.path.basename(csv_path)
                print(f"  Processing file: {filename}")
                
                # Create output path maintaining directory structure
                output_path = os.path.join(output_base, subdir, filename)
                
                try:
                    process_matrix_normal(csv_path, output_path, args.k)
                    print(f"    Saved to: {output_path}")
                    
                except Exception as e:
                    print(f"    Error processing {filename}: {e}")
        
        print(f"Imputed matrices saved in {output_base}/ directory")

if __name__ == "__main__":
    main()
