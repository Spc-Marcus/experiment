#!/usr/bin/env python3
"""
CSV Matrix Rotation Pipeline

This script processes multiple CSV files in a directory structure,
transposing each matrix and saving the results with a configurable suffix.

Usage:
    python process_csv_folder.py <input_folder> <output_folder> [options]
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple
import time

# Import our matrix rotator module
from matrix_rotator import MatrixRotator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def find_csv_files(input_folder: Path, recursive: bool = False) -> List[Path]:
    """
    Find all CSV files in the input folder.
    
    Args:
        input_folder (Path): Directory to search for CSV files
        recursive (bool): Whether to search subdirectories recursively
        
    Returns:
        List[Path]: List of CSV file paths found
    """
    csv_files = []
    
    if recursive:
        # Search recursively in all subdirectories
        pattern = "**/*.csv"
        csv_files = list(input_folder.glob(pattern))
    else:
        # Search only in the input folder
        csv_files = list(input_folder.glob("*.csv"))
    
    # Filter out any directories that might match the pattern
    csv_files = [f for f in csv_files if f.is_file()]
    
    return sorted(csv_files)


def create_output_path(input_file: Path, output_folder: Path, 
                      suffix: str = "_rotated") -> Path:
    """
    Create the output path for a transposed matrix.
    
    Args:
        input_file (Path): Input CSV file path
        output_folder (Path): Base output directory
        suffix (str): Suffix to add to filename
        
    Returns:
        Path: Output file path
    """
    # Get the relative path from input folder to preserve directory structure
    try:
        # If input_file is inside input_folder, preserve the relative structure
        relative_path = input_file.relative_to(input_file.parents[len(input_file.parts) - len(output_folder.parts) - 1])
    except ValueError:
        # If not, just use the filename
        relative_path = Path(input_file.name)
    
    # Create output filename with suffix
    output_filename = input_file.stem + suffix + input_file.suffix
    output_path = output_folder / relative_path.parent / output_filename
    
    return output_path


def process_csv_folder(input_folder: str, output_folder: str, 
                      suffix: str = "_rotated", verbose: bool = False,
                      recursive: bool = False) -> Tuple[int, int]:
    """
    Process all CSV files in a folder and its subdirectories.
    
    Args:
        input_folder (str): Input directory containing CSV files
        output_folder (str): Output directory for rotated matrices
        suffix (str): Suffix to add to output filenames
        verbose (bool): Enable verbose logging
        recursive (bool): Process subdirectories recursively
        
    Returns:
        Tuple[int, int]: (total_files, successful_files)
    """
    input_path = Path(input_folder)
    output_path = Path(output_folder)
    
    # Validate input directory
    if not input_path.exists():
        logger.error(f"Input directory does not exist: {input_folder}")
        return 0, 0
    
    if not input_path.is_dir():
        logger.error(f"Input path is not a directory: {input_folder}")
        return 0, 0
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all CSV files
    csv_files = find_csv_files(input_path, recursive=recursive)
    
    if not csv_files:
        logger.warning(f"No CSV files found in {input_folder}")
        if recursive:
            logger.info("Searched recursively in all subdirectories")
        return 0, 0
    
    logger.info(f"Found {len(csv_files)} CSV files to process")
    if verbose:
        for file in csv_files:
            logger.debug(f"  - {file}")
    
    # Initialize the matrix rotator
    rotator = MatrixRotator(verbose=verbose)
    
    # Process each file
    successful_files = 0
    start_time = time.time()
    
    for i, csv_file in enumerate(csv_files, 1):
        try:
            # Create output path
            output_file = create_output_path(csv_file, output_path, suffix)
            
            if verbose:
                logger.info(f"Processing file {i}/{len(csv_files)}: {csv_file.name}")
                logger.info(f"  Input:  {csv_file}")
                logger.info(f"  Output: {output_file}")
            
            # Process the file
            success = rotator.process_single_file(csv_file, output_file)
            
            if success:
                successful_files += 1
                if verbose:
                    logger.info(f"  ✓ Successfully rotated")
                else:
                    logger.info(f"✓ {csv_file.name}")
            else:
                logger.error(f"  ✗ Failed to process {csv_file.name}")
                
        except Exception as e:
            logger.error(f"Error processing {csv_file.name}: {e}")
    
    # Summary
    elapsed_time = time.time() - start_time
    logger.info(f"\nProcessing complete!")
    logger.info(f"Total files: {len(csv_files)}")
    logger.info(f"Successful: {successful_files}")
    logger.info(f"Failed: {len(csv_files) - successful_files}")
    logger.info(f"Time elapsed: {elapsed_time:.2f} seconds")
    
    return len(csv_files), successful_files


def main():
    """Main function to handle command line arguments and execute the pipeline."""
    parser = argparse.ArgumentParser(
        description="Rotate (transpose) CSV matrices in a folder structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python process_csv_folder.py ./data/matrices ./data/rotated
  
  # With custom suffix
  python process_csv_folder.py ./data/matrices ./data/rotated --suffix "_transposed"
  
  # Verbose output with recursive processing
  python process_csv_folder.py ./data/matrices ./data/rotated --verbose --recursive
  
  # Process specific directory
  python process_csv_folder.py /path/to/matrices /path/to/output --suffix "_rotated"
        """
    )
    
    parser.add_argument(
        "input_folder",
        help="Directory containing CSV files to rotate"
    )
    
    parser.add_argument(
        "output_folder", 
        help="Output directory for rotated matrices"
    )
    
    parser.add_argument(
        "--suffix",
        default="_rotated",
        help="Suffix to add to output filenames (default: '_rotated')"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output with detailed progress information"
    )
    
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Process subdirectories recursively"
    )
    
    args = parser.parse_args()
    
    # Set logging level based on verbose flag
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)
    
    # Display configuration
    logger.info("CSV Matrix Rotation Pipeline")
    logger.info("=" * 40)
    logger.info(f"Input folder:  {args.input_folder}")
    logger.info(f"Output folder: {args.output_folder}")
    logger.info(f"Suffix:        {args.suffix}")
    logger.info(f"Recursive:     {args.recursive}")
    logger.info(f"Verbose:       {args.verbose}")
    logger.info("=" * 40)
    
    try:
        # Process the folder
        total_files, successful_files = process_csv_folder(
            input_folder=args.input_folder,
            output_folder=args.output_folder,
            suffix=args.suffix,
            verbose=args.verbose,
            recursive=args.recursive
        )
        
        # Exit with appropriate code
        if total_files == 0:
            logger.warning("No files were processed")
            sys.exit(0)
        elif successful_files == total_files:
            logger.info("All files processed successfully!")
            sys.exit(0)
        else:
            logger.error(f"Some files failed to process ({successful_files}/{total_files} successful)")
            sys.exit(1)
            
    except KeyboardInterrupt:
        logger.info("\nProcessing interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 