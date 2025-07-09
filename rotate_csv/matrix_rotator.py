"""
CSV Matrix Rotation Tool

This module provides functionality to transpose (rotate) CSV matrices,
effectively swapping rows and columns while preserving all data values.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path
from typing import Union, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MatrixRotator:
    """
    A class to handle CSV matrix transposition operations.
    
    This class provides methods to load, transpose, and save CSV matrices
    while maintaining data integrity and providing detailed logging.
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the MatrixRotator.
        
        Args:
            verbose (bool): Enable verbose logging output
        """
        self.verbose = verbose
        if verbose:
            logging.getLogger().setLevel(logging.DEBUG)
    
    def load_matrix(self, file_path: Union[str, Path]) -> pd.DataFrame:
        """
        Load a CSV matrix from file.
        
        Args:
            file_path (Union[str, Path]): Path to the CSV file
            
        Returns:
            pd.DataFrame: Loaded matrix as pandas DataFrame
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            pd.errors.EmptyDataError: If the file is empty
            pd.errors.ParserError: If the CSV format is invalid
        """
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        if self.verbose:
            logger.info(f"Loading matrix from: {file_path}")
        
        try:
            # Load CSV with index_col=0 to preserve the first column as index
            matrix = pd.read_csv(file_path, index_col=0)
            
            if matrix.empty:
                logger.warning(f"File {file_path} is empty or contains no data")
                return matrix
            
            if self.verbose:
                logger.info(f"Loaded matrix with shape: {matrix.shape}")
                logger.info(f"Matrix columns: {list(matrix.columns)}")
                logger.info(f"Matrix index: {list(matrix.index)}")
            
            return matrix
            
        except pd.errors.EmptyDataError:
            logger.error(f"File {file_path} is empty")
            raise
        except pd.errors.ParserError as e:
            logger.error(f"Error parsing CSV file {file_path}: {e}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error loading {file_path}: {e}")
            raise
    
    def transpose_matrix(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """
        Transpose a matrix (swap rows and columns).
        
        Args:
            matrix (pd.DataFrame): Input matrix to transpose
            
        Returns:
            pd.DataFrame: Transposed matrix
            
        Raises:
            ValueError: If the matrix is empty or invalid
        """
        if matrix.empty:
            raise ValueError("Cannot transpose an empty matrix")
        
        if self.verbose:
            logger.info(f"Transposing matrix with shape: {matrix.shape}")
        
        try:
            # Transpose the matrix
            transposed = matrix.transpose()
            
            if self.verbose:
                logger.info(f"Transposed matrix shape: {transposed.shape}")
                logger.info(f"Original shape: {matrix.shape} -> Transposed shape: {transposed.shape}")
            
            return transposed
            
        except Exception as e:
            logger.error(f"Error transposing matrix: {e}")
            raise
    
    def save_matrix(self, matrix: pd.DataFrame, output_path: Union[str, Path]) -> None:
        """
        Save a matrix to CSV file.
        
        Args:
            matrix (pd.DataFrame): Matrix to save
            output_path (Union[str, Path]): Output file path
            
        Raises:
            PermissionError: If unable to write to the output location
            Exception: For other file writing errors
        """
        output_path = Path(output_path)
        
        # Create output directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if self.verbose:
            logger.info(f"Saving transposed matrix to: {output_path}")
        
        try:
            # Save with index=True to preserve row names
            matrix.to_csv(output_path, index=True)
            
            if self.verbose:
                logger.info(f"Successfully saved matrix to: {output_path}")
                
        except PermissionError:
            logger.error(f"Permission denied writing to: {output_path}")
            raise
        except Exception as e:
            logger.error(f"Error saving matrix to {output_path}: {e}")
            raise
    
    def process_single_file(self, input_path: Union[str, Path], 
                          output_path: Union[str, Path]) -> bool:
        """
        Process a single CSV file: load, transpose, and save.
        
        Args:
            input_path (Union[str, Path]): Input CSV file path
            output_path (Union[str, Path]): Output CSV file path
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Load the matrix
            matrix = self.load_matrix(input_path)
            
            if matrix.empty:
                logger.warning(f"Skipping empty file: {input_path}")
                return False
            
            # Transpose the matrix
            transposed = self.transpose_matrix(matrix)
            
            # Save the transposed matrix
            self.save_matrix(transposed, output_path)
            
            return True
            
        except Exception as e:
            logger.error(f"Error processing {input_path}: {e}")
            return False
    
    def get_output_filename(self, input_path: Union[str, Path], 
                          suffix: str = "_rotated") -> str:
        """
        Generate output filename with suffix.
        
        Args:
            input_path (Union[str, Path]): Input file path
            suffix (str): Suffix to add to filename
            
        Returns:
            str: Output filename
        """
        input_path = Path(input_path)
        stem = input_path.stem
        suffix_with_extension = suffix + input_path.suffix
        return stem + suffix_with_extension


def rotate_single_matrix(input_file: Union[str, Path], 
                        output_file: Union[str, Path], 
                        verbose: bool = False) -> bool:
    """
    Convenience function to rotate a single matrix file.
    
    Args:
        input_file (Union[str, Path]): Input CSV file path
        output_file (Union[str, Path]): Output CSV file path
        verbose (bool): Enable verbose logging
        
    Returns:
        bool: True if successful, False otherwise
    """
    rotator = MatrixRotator(verbose=verbose)
    return rotator.process_single_file(input_file, output_file)


def validate_matrix_format(matrix: pd.DataFrame) -> bool:
    """
    Validate that a matrix has the expected format for rotation.
    
    Args:
        matrix (pd.DataFrame): Matrix to validate
        
    Returns:
        bool: True if valid, False otherwise
    """
    if matrix.empty:
        return False
    
    # Check if matrix has at least 2 rows and 2 columns
    if matrix.shape[0] < 2 or matrix.shape[1] < 2:
        logger.warning("Matrix is too small for meaningful transposition")
        return False
    
    # Check for non-numeric data (warn but don't fail)
    non_numeric_cols = []
    for col in matrix.columns:
        if not pd.api.types.is_numeric_dtype(matrix[col]):
            non_numeric_cols.append(col)
    
    if non_numeric_cols:
        logger.warning(f"Non-numeric columns found: {non_numeric_cols}")
    
    return True


if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python matrix_rotator.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    success = rotate_single_matrix(input_file, output_file, verbose=True)
    
    if success:
        print(f"Successfully rotated {input_file} -> {output_file}")
    else:
        print(f"Failed to rotate {input_file}")
        sys.exit(1) 