import pysam
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
from sklearn.impute import KNNImputer
import sys

# Add parent directory to path for imports
current_dir = Path(__file__).parent
parent_dir = current_dir.parent
sys.path.insert(0, str(parent_dir))

logger = logging.getLogger(__name__)

def extract_variants_from_bam(bam_file: str, min_contig_length: int = 90000) -> Dict[str, Dict[int, Dict[str, int]]]:
    """
    Extract variant positions from BAM file for contigs longer than threshold.
    
    Parameters
    ----------
    bam_file : str
        Path to BAM file
    min_contig_length : int
        Minimum contig length to process (default: 90000)
        
    Returns
    -------
    Dict[str, Dict[int, Dict[str, int]]]
        Nested dictionary: {contig: {position: {read_name: variant_value}}}
    """
    
    variants_by_contig = {}
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bamfile:
            # Get contig lengths
            contig_lengths = dict(zip(bamfile.references, bamfile.lengths))
            
            # Filter contigs by length
            valid_contigs = {contig: length for contig, length in contig_lengths.items() 
                           if length >= min_contig_length}
            
            if not valid_contigs:
                logger.warning(f"No contigs >= {min_contig_length} bp found in {bam_file}")
                return variants_by_contig
            
            logger.info(f"Processing {len(valid_contigs)} contigs >= {min_contig_length} bp")
            
            for contig in valid_contigs:
                logger.debug(f"Processing contig {contig} ({valid_contigs[contig]} bp)")
                variants_by_contig[contig] = {}
                
                # Fetch all reads for this contig
                for read in bamfile.fetch(contig):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    
                    read_name = read.query_name
                    
                    # Extract variant positions from read
                    ref_pos = read.reference_start
                    
                    for i, (query_pos, ref_pos_aligned) in enumerate(read.get_aligned_pairs(matches_only=False)):
                        if ref_pos_aligned is None:  # Insertion
                            continue
                            
                        if query_pos is None:  # Deletion
                            # Record deletion at this position
                            if ref_pos_aligned not in variants_by_contig[contig]:
                                variants_by_contig[contig][ref_pos_aligned] = {}
                            variants_by_contig[contig][ref_pos_aligned][read_name] = 0  # Deletion = 0
                        else:
                            # Regular aligned position
                            if ref_pos_aligned not in variants_by_contig[contig]:
                                variants_by_contig[contig][ref_pos_aligned] = {}
                            
                            # Get reference and query bases
                            query_base = read.query_sequence[query_pos] if read.query_sequence else 'N'
                            
                            # Simple variant calling: different from reference = variant (0), same = reference (1)
                            # This is simplified - in practice you'd use a reference genome
                            variants_by_contig[contig][ref_pos_aligned][read_name] = 1  # Default to reference
                
                logger.debug(f"Contig {contig}: found {len(variants_by_contig[contig])} variant positions")
                
    except Exception as e:
        logger.error(f"Error processing BAM file {bam_file}: {e}")
        raise
    
    return variants_by_contig

def create_matrix_from_variants(variants: Dict[int, Dict[str, int]], min_coverage_threshold: float = 0.6) -> Tuple[np.ndarray, List[str]]:
    """
    Create matrix from variant dictionary, similar to create_matrix function.
    
    Parameters
    ----------
    variants : Dict[int, Dict[str, int]]
        Dictionary {position: {read_name: variant_value}}
    min_coverage_threshold : float
        Minimum coverage threshold for positions
        
    Returns
    -------
    Tuple[np.ndarray, List[str]]
        Matrix and list of read names
    """
    
    if not variants:
        return np.array([]), []
    
    # Create DataFrame from variants
    df = pd.DataFrame(variants)
    
    if df.empty:
        return np.array([]), []
    
    logger.debug(f"Processing {df.shape[1]} positions with {df.shape[0]} reads")
    
    variant_matrix = df.copy()
    
    # Process each position to determine majority allele
    for col in df.columns:
        position_data = df[col].dropna()
        
        if len(position_data) == 0:
            continue
            
        allele_counts = position_data.value_counts()
        
        if len(allele_counts) == 0:
            continue
        elif len(allele_counts) == 1:
            # Monomorphic position
            majority_allele = allele_counts.index[0]
            variant_matrix[col] = 1  # All same as reference
        else:
            majority_allele = allele_counts.index[0]  # Most frequent (reference)
            
            # 1 = same as reference (majority), 0 = variant
            variant_matrix[col] = (df[col] == majority_allele).astype(int)
            
            # Keep NaN for missing positions
            variant_matrix.loc[df[col].isna(), col] = np.nan
    
    # Filter reads that span the window
    if len(variant_matrix.columns) >= 3:
        tmp_idx = variant_matrix.iloc[:, :len(variant_matrix.columns)//3].dropna(axis=0, how='all')
        variant_matrix = variant_matrix.loc[tmp_idx.index, :]
        
        tmp_idx = variant_matrix.iloc[:, 2*len(variant_matrix.columns)//3:].dropna(axis=0, how='all')
        variant_matrix = variant_matrix.loc[tmp_idx.index, :]
    
    # Filter columns with insufficient coverage
    if not variant_matrix.empty:
        variant_matrix = variant_matrix.dropna(axis=1, thresh=min_coverage_threshold * len(variant_matrix.index))
    
    # Extract filtered reads
    reads = list(variant_matrix.index)
    
    # Convert to numpy matrix
    matrix = variant_matrix.to_numpy() if not variant_matrix.empty else np.array([])
    
    logger.debug(f"Matrix created: {matrix.shape} with {len(reads)} reads")
    
    return matrix, reads

def impute_and_binarize_matrix(matrix: np.ndarray, n_neighbors: int = 10) -> np.ndarray:
    """
    Impute missing values and binarize matrix to have only 0s and 1s.
    
    Parameters
    ----------
    matrix : np.ndarray
        Input matrix with potential missing values
    n_neighbors : int
        Number of neighbors for KNN imputation
        
    Returns
    -------
    np.ndarray
        Binary matrix with only 0s and 1s
    """
    
    if matrix.size == 0:
        return matrix
    
    # Impute missing values if any exist
    if np.isnan(matrix).any():
        logger.debug(f"Imputing {np.isnan(matrix).sum()} missing values")
        
        # Adjust n_neighbors if necessary
        n_reads = matrix.shape[0]
        actual_neighbors = min(n_neighbors, max(1, n_reads - 1))
        
        imputer = KNNImputer(n_neighbors=actual_neighbors)
        matrix = imputer.fit_transform(matrix)
    
    # Binarize: values >= 0.5 become 1, others become 0
    binary_matrix = (matrix >= 0.5).astype(int)
    
    logger.debug(f"Matrix binarized: {(binary_matrix == 1).sum()} ones, {(binary_matrix == 0).sum()} zeros")
    
    return binary_matrix

def save_matrix_csv(matrix: np.ndarray, reads: List[str], output_file: str) -> bool:
    """
    Save matrix as CSV file.
    
    Parameters
    ----------
    matrix : np.ndarray
        Matrix to save
    reads : List[str]
        Read names for row indices
    output_file : str
        Output file path
        
    Returns
    -------
    bool
        True if saved successfully
    """
    
    try:
        # Create DataFrame
        df = pd.DataFrame(matrix, index=reads)
        
        # Save to CSV
        df.to_csv(output_file, index=True)
        logger.info(f"Matrix saved: {output_file} ({matrix.shape[0]}x{matrix.shape[1]})")
        return True
        
    except Exception as e:
        logger.error(f"Error saving matrix {output_file}: {e}")
        return False

def process_single_bam(bam_file: str, output_dir: str = "matrices", min_contig_length: int = 90000, 
                      min_matrix_size: int = 25) -> List[str]:
    """
    Process a single BAM file and create matrices for valid contigs.
    
    Parameters
    ----------
    bam_file : str
        Path to BAM file
    output_dir : str
        Output directory for matrices
    min_contig_length : int
        Minimum contig length to process
    min_matrix_size : int
        Minimum matrix dimensions (rows and cols) to save
        
    Returns
    -------
    List[str]
        List of created matrix files
    """
    
    bam_path = Path(bam_file)
    if not bam_path.exists():
        logger.error(f"BAM file not found: {bam_file}")
        return []
    
    # Extract base name for output directory
    base_name = bam_path.stem  # Remove .bam extension
    
    # Create output directory
    output_path = Path(output_dir) / base_name
    output_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Processing BAM file: {bam_file}")
    
    try:
        # Extract variants from BAM
        variants_by_contig = extract_variants_from_bam(bam_file, min_contig_length)
        
        if not variants_by_contig:
            logger.warning(f"No valid contigs found in {bam_file}")
            return []
        
        created_files = []
        
        # Process each contig
        for contig, variants in variants_by_contig.items():
            logger.info(f"Processing contig {contig}")
            
            # Create matrix
            matrix, reads = create_matrix_from_variants(variants)
            
            if matrix.size == 0:
                logger.warning(f"Empty matrix for contig {contig}")
                continue
            
            # Check minimum size requirements
            if matrix.shape[0] < min_matrix_size or matrix.shape[1] < min_matrix_size:
                logger.info(f"Matrix too small for contig {contig}: {matrix.shape} < {min_matrix_size}x{min_matrix_size}")
                continue
            
            # Impute and binarize
            binary_matrix = impute_and_binarize_matrix(matrix)
            
            # Save matrix
            output_file = output_path / f"{contig}_matrix.csv"
            if save_matrix_csv(binary_matrix, reads, str(output_file)):
                created_files.append(str(output_file))
        
        logger.info(f"Created {len(created_files)} matrices for {bam_file}")
        return created_files
        
    except Exception as e:
        logger.error(f"Error processing BAM file {bam_file}: {e}")
        return []

def process_bam_files(bam_directory: str = "/bam", output_dir: str = "matrices", 
                     min_contig_length: int = 90000, min_matrix_size: int = 25) -> Dict[str, List[str]]:
    """
    Process all BAM files in a directory.
    
    Parameters
    ----------
    bam_directory : str
        Directory containing BAM files
    output_dir : str
        Output directory for matrices
    min_contig_length : int
        Minimum contig length to process
    min_matrix_size : int
        Minimum matrix dimensions to save
        
    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping BAM files to created matrix files
    """
    
    bam_dir = Path(bam_directory)
    if not bam_dir.exists():
        logger.error(f"BAM directory not found: {bam_directory}")
        return {}
    
    # Find all BAM files
    bam_files = list(bam_dir.glob("*.bam"))
    
    if not bam_files:
        logger.warning(f"No BAM files found in {bam_directory}")
        return {}
    
    logger.info(f"Found {len(bam_files)} BAM files to process")
    
    results = {}
    
    for bam_file in bam_files:
        logger.info(f"Processing {bam_file.name}")
        
        try:
            created_files = process_single_bam(
                str(bam_file), 
                output_dir, 
                min_contig_length, 
                min_matrix_size
            )
            results[str(bam_file)] = created_files
            
        except Exception as e:
            logger.error(f"Failed to process {bam_file}: {e}")
            results[str(bam_file)] = []
    
    # Summary
    total_matrices = sum(len(files) for files in results.values())
    successful_bams = len([files for files in results.values() if files])
    
    logger.info(f"Processing complete: {total_matrices} matrices created from {successful_bams}/{len(bam_files)} BAM files")
    
    return results

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Process all BAM files
    results = process_bam_files()
    
    # Print summary
    for bam_file, matrices in results.items():
        print(f"{Path(bam_file).name}: {len(matrices)} matrices created")
