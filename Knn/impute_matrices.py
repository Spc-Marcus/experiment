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
    df = pd.read_csv(csv_path, index_col=0, na_values=['', ' ', 'NA', 'nan', 'NaN', '-1'])
    return df

def save_matrix(df, output_path):
    """Save matrix to CSV file preserving original format."""
    # Only create directory if output_path contains a directory part
    dir_name = os.path.dirname(output_path)
    if dir_name and dir_name != '':
        os.makedirs(dir_name, exist_ok=True)
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

def process_single_file(input_path, k, matrix_output_dir=None, exp_output_dir=None, experiment=False):
    """Process a single CSV file."""
    if not os.path.exists(input_path):
        print(f"Error: File {input_path} does not exist")
        return
    
    if not input_path.endswith('.csv'):
        print(f"Error: File must be a CSV file")
        return
    
    print(f"Processing single file: {input_path}")
    filename = os.path.basename(input_path)
    
    # Always perform normal processing if matrix_output_dir is provided
    if matrix_output_dir:
        print(f"Running in normal mode with k={k}...")
        
        # Create output filename
        os.makedirs(matrix_output_dir, exist_ok=True)
        output_name = os.path.join(matrix_output_dir, f"{filename[:-4]}_knn_imputed.csv")
        
        try:
            process_matrix_normal(input_path, output_name, k)
            print(f"Imputed matrix saved to: {output_name}")
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    # Additionally perform experiment if requested
    if experiment:
        if not exp_output_dir:
            print("Error: --exp-output is required in experiment mode")
            return
            
        print("Running in experiment mode...")
        k_values = list(range(5, 26))
        
        try:
            results = process_matrix_experiment(input_path, k_values)
            
            # Add file info to results
            for result in results:
                result['fichier'] = filename
            
            # Save experiment results
            results_df = pd.DataFrame(results)
            results_df = results_df[['fichier', 'k_value', 'taille_matrice', 'nan_avant', 'nan_après']]
            
            os.makedirs(exp_output_dir, exist_ok=True)
            output_name = os.path.join(exp_output_dir, "KNN_exp.csv")
            
            results_df.to_csv(output_name, index=False)
            print(f"Experiment results saved to {output_name}")
            
        except Exception as e:
            print(f"Error processing experiment for {filename}: {e}")

def process_directory(input_dir, k, matrix_output_dir=None, exp_output_dir=None, experiment=False):
    """Process all CSV files in directory structure."""
    csv_files = get_csv_files(input_dir)
    
    # Always perform normal processing if matrix_output_dir is provided
    if matrix_output_dir:
        print(f"Running in normal mode with k={k}...")
        
        for subdir, files in csv_files.items():
            print(f"Processing subdirectory: {subdir}")
            
            for csv_path in files:
                filename = os.path.basename(csv_path)
                print(f"  Processing file: {filename}")
                
                # Create output path maintaining directory structure
                output_path = os.path.join(matrix_output_dir, subdir, filename)
                
                try:
                    process_matrix_normal(csv_path, output_path, k)
                    print(f"    Saved to: {output_path}")
                    
                except Exception as e:
                    print(f"    Error processing {filename}: {e}")
        
        print(f"Imputed matrices saved in {matrix_output_dir}/ directory")
    
    # Additionally perform experiment if requested
    if experiment:
        if not exp_output_dir:
            print("Error: --exp-output is required in experiment mode")
            return
            
        print("Running in experiment mode...")
        k_values = list(range(5, 26))
        csv_files_exp = get_csv_files(input_dir, max_files=150)
        
        all_results = []
        
        for subdir, files in csv_files_exp.items():
            print(f"Processing experiment for subdirectory: {subdir}")
            
            for csv_path in files:
                filename = os.path.basename(csv_path)
                print(f"  Experimenting with file: {filename}")
                
                try:
                    results = process_matrix_experiment(csv_path, k_values)
                    
                    for result in results:
                        result['sous_dossier'] = subdir
                        result['fichier'] = filename
                        all_results.append(result)
                        
                except Exception as e:
                    print(f"    Error in experiment for {filename}: {e}")
        
        # Save experiment results
        results_df = pd.DataFrame(all_results)
        results_df = results_df[['sous_dossier', 'fichier', 'k_value', 'taille_matrice', 'nan_avant', 'nan_après']]
        
        os.makedirs(exp_output_dir, exist_ok=True)
        results_path = os.path.join(exp_output_dir, "KNN_exp.csv")
        
        results_df.to_csv(results_path, index=False)
        print(f"Experiment results saved to {results_path}")

def main():
    parser = argparse.ArgumentParser(description="Impute missing values in genomic matrices using KNN")
    parser.add_argument("input", help="Input directory containing subdirectories with CSV matrices OR single CSV file")
    parser.add_argument("--k", type=int, default=10, help="KNN value (default: 10)")
    parser.add_argument("--matrix-output", help="Output directory for KNN imputed matrices")
    parser.add_argument("--exp-output", help="Output directory for experiment results (required if --experiment is used)")
    parser.add_argument("--experiment", action="store_true", help="Run in experiment mode (in addition to normal processing)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input path {args.input} does not exist")
        return
    
    # Validate required arguments
    if not args.matrix_output and not args.experiment:
        print("Error: Either --matrix-output or --experiment (or both) must be specified")
        return
    
    if args.experiment and not args.exp_output:
        print("Error: --exp-output is required when --experiment is used")
        return
    
    # Check if input is a file or directory
    if os.path.isfile(args.input):
        process_single_file(args.input, args.k, args.matrix_output, args.exp_output, args.experiment)
    elif os.path.isdir(args.input):
        process_directory(args.input, args.k, args.matrix_output, args.exp_output, args.experiment)
    else:
        print(f"Error: {args.input} is neither a file nor a directory")

if __name__ == "__main__":
    main()