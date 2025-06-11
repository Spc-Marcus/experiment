import argparse
import sys
from pathlib import Path
import pysam as ps
import pandas as pd
import numpy as np

# Ensure we can import from the same directory as this script
script_dir = Path(__file__).parent.absolute()
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

try:
    from core import get_data
    from matrix import create_matrix
except ImportError as e:
    print(f"Error importing modules: {e}")
    print(f"Make sure core.py and matrix.py are in: {script_dir}")
    sys.exit(1)

def process_single_bam(bam_file_path: Path, output_base_dir: Path, 
                      window_size: int = 100000, overlap: int = 10000,
                      filtered_col_threshold: float = 0.6):
    """
    Process a single BAM file and save matrices for all contigs.
    
    Parameters
    ----------
    bam_file_path : Path
        Path to the BAM file
    output_base_dir : Path
        Base output directory (matrices/)
    window_size : int
        Size of genomic windows to process
    overlap : int
        Overlap between consecutive windows
    filtered_col_threshold : float
        Coverage threshold for columns filtering
    """
    # Get BAM filename without extension for output folder
    bam_name = bam_file_path.stem  # removes .bam extension
    
    # Create output directory: matrices/X/ where X is BAM name
    output_dir = output_base_dir / bam_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing {bam_file_path} -> {output_dir}")
    
    total_matrices_saved = 0
    total_variants = 0
    
    try:
        with ps.AlignmentFile(str(bam_file_path), "rb") as bam_file:
            # Get all contigs from BAM file
            contigs = bam_file.references
            contig_lengths = bam_file.lengths
            
            print(f"Found {len(contigs)} contigs in {bam_name}")
            
            for contig_name, contig_length in zip(contigs, contig_lengths):
                print(f"  Processing contig {contig_name} (length: {contig_length:,} bp)")
                
                # Use different window sizes based on contig length
                if contig_length > 90000:
                    # For large contigs, use small windows of 5000 bp
                    contig_window_size = 5000
                    contig_overlap = 0  # No overlap for small windows
                    print(f"    Large contig detected - using {contig_window_size} bp windows")
                else:
                    # For small contigs, use provided window size
                    contig_window_size = window_size
                    contig_overlap = overlap
                
                # Process contig in sliding windows
                start_pos = 0
                window_count = 0
                
                while start_pos < contig_length:
                    end_pos = min(start_pos + contig_window_size, contig_length)
                    
                    try:
                        # Extract variant data
                        dict_of_sus_pos, contig_size = get_data(
                            bam_file, contig_name, start_pos, end_pos
                        )
                        
                        total_variants += len(dict_of_sus_pos)
                        
                        if dict_of_sus_pos:
                            # Create and filter matrix
                            X_matrix, reads = create_matrix(
                                dict_of_sus_pos,
                                filtered_col_threshold,
                                contig_size=contig_size,
                                contig_name=contig_name,
                                start_pos=start_pos,
                                csv_output_dir=output_dir
                            )
                            
                            # Check if CSV was saved
                            if (contig_size > 90000 and 
                                X_matrix.size > 0 and 
                                X_matrix.shape[0] >= 20 and 
                                X_matrix.shape[1] >= 20):
                                total_matrices_saved += 1
                                print(f"    Window {start_pos}-{end_pos}: CSV saved (matrix: {X_matrix.shape})")
                            else:
                                print(f"    Window {start_pos}-{end_pos}: {len(dict_of_sus_pos)} variants, matrix: {X_matrix.shape if X_matrix.size > 0 else 'empty'}")
                        
                        window_count += 1
                        
                    except Exception as e:
                        print(f"    Error processing window {start_pos}-{end_pos}: {e}")
                    
                    # Move to next window
                    start_pos += contig_window_size - contig_overlap
                
                print(f"    Completed {window_count} windows for {contig_name}")
    
    except Exception as e:
        print(f"Error processing {bam_file_path}: {e}")
        return False
    
    print(f"Completed {bam_name}: {total_variants} total variants, {total_matrices_saved} CSV matrices saved")
    return True

def process_bam_folder(bam_folder: Path, output_folder: Path, 
                      window_size: int = 100000, overlap: int = 10000,
                      filtered_col_threshold: float = 0.6):
    """
    Process all BAM files in a folder.
    
    Parameters
    ----------
    bam_folder : Path
        Directory containing BAM files
    output_folder : Path
        Output directory (will create matrices_no_binarize/ subdirectory)
    window_size : int
        Size of genomic windows to process
    overlap : int
        Overlap between consecutive windows
    filtered_col_threshold : float
        Coverage threshold for columns filtering
    """
    # Find all BAM files
    bam_files = list(bam_folder.glob("*.bam"))
    
    if not bam_files:
        print(f"No BAM files found in {bam_folder}")
        return
    
    print(f"Found {len(bam_files)} BAM files to process")
    
    # Create output directory structure
    matrices_dir = output_folder / f"matrices_no_binarize_{filtered_col_threshold}"
    matrices_dir.mkdir(parents=True, exist_ok=True)
    
    successful = 0
    failed = 0
    
    for bam_file in bam_files:
        print(f"\n{'='*60}")
        print(f"Processing BAM file {successful + failed + 1}/{len(bam_files)}: {bam_file.name}")
        
        # Check if BAM file is indexed
        bai_file = Path(str(bam_file) + ".bai")
        if not bai_file.exists():
            print(f"Error: Index file {bai_file.name} not found.")
            print(f"Please create index with: samtools index {bam_file}")
            failed += 1
            continue
        
        try:
            success = process_single_bam(
                bam_file, 
                matrices_dir, 
                window_size, 
                overlap, 
                filtered_col_threshold
            )
            
            if success:
                successful += 1
            else:
                failed += 1
                
        except Exception as e:
            print(f"Failed to process {bam_file.name}: {e}")
            failed += 1
    
    print(f"\n{'='*60}")
    print(f"Processing complete: {successful} successful, {failed} failed")
    print(f"Output saved to: {matrices_dir}")

def main():
    parser = argparse.ArgumentParser(
        description='Process all BAM files in a directory and create CSV matrices',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python process_bam_folder.py bam ./output
  python process_bam_folder.py /path/to/bam_files /path/to/output
  python process_bam_folder.py /data/bam /results --window-size 50000 --threshold 0.7
        
Output structure:
  output_folder/
  └── matrices_no_binarize/
      ├── sample1/
      │   ├── matrix_chr1_0.csv
      │   ├── matrix_chr1_100000.csv
      │   └── ...
      ├── sample2/
      │   ├── matrix_chr1_0.csv
      │   └── ...
      └── sample3/
          └── ...
        """
    )
    
    parser.add_argument('bam_folder', type=Path, help='Directory containing BAM files')
    parser.add_argument('output_folder', type=Path, help='Output directory (will create matrices/ subdirectory)')
    parser.add_argument('--window-size', type=int, default=5000, 
                       help='Size of genomic windows (default: 5000)')
    parser.add_argument('--overlap', type=int, default=10000,
                       help='Overlap between windows (default: 10000)')
    parser.add_argument('--threshold', type=float, default=0.6,
                       help='Minimum read coverage threshold: reads must cover at least this fraction of variant positions (default: 0.6)')
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    bam_folder = args.bam_folder.resolve()
    output_folder = args.output_folder.resolve()
    
    # Validate arguments
    if not bam_folder.exists():
        print(f"Error: BAM folder {bam_folder} does not exist")
        sys.exit(1)
    
    if not bam_folder.is_dir():
        print(f"Error: {bam_folder} is not a directory")
        sys.exit(1)
    
    # Check if we can write to output directory
    try:
        output_folder.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        print(f"Error: Permission denied to create/write to {output_folder}")
        print("Try using a directory in your home folder or current working directory")
        sys.exit(1)
    except Exception as e:
        print(f"Error creating output directory {output_folder}: {e}")
        sys.exit(1)
    
    if args.window_size <= 0:
        print("Error: Window size must be positive")
        sys.exit(1)
    
    if not 0.0 <= args.threshold <= 1.0:
        print("Error: Threshold must be between 0.0 and 1.0")
        sys.exit(1)
    
    print(f"Input BAM folder: {bam_folder}")
    print(f"Output folder: {output_folder}")
    
    # Process all BAM files
    try:
        process_bam_folder(
            bam_folder,
            output_folder, 
            args.window_size,
            args.overlap,
            args.threshold
        )
    except KeyboardInterrupt:
        print("\nProcessing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
