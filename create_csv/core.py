import os
import time
from pathlib import Path
from typing import Dict, Tuple, List
import numpy as np
import pysam as ps




def get_data(input_file: ps.AlignmentFile, 
             contig_name: str, 
             start_pos: int, 
             end_pos: int,
             min_coverage: int = 5,
             min_base_quality: int = 10,
             max_major_allele_freq: float = 0.95) -> Tuple[Dict[int, Dict[str, int]], int]:
    """
    Extract suspicious genomic positions with variant alleles from alignment data.
    
    This function performs pileup analysis on a specified genomic region to identify
    positions where multiple alleles are present with sufficient frequency to suggest
    genuine variants (as opposed to sequencing errors). It's the primary data extraction
    function for downstream strain separation analysis.
    
    The algorithm identifies heterozygous positions by examining allele frequencies
    at each covered position, filtering for positions where the major allele frequency
    is below a specified threshold, indicating potential strain variants or SNPs.
    
    Parameters
    ----------
    input_file : ps.AlignmentFile
        Opened pysam AlignmentFile object (BAM/SAM) containing mapped sequencing reads.
        File must be properly indexed for efficient region-based access.
    contig_name : str
        Name of the reference contig/chromosome to analyze. Must match header names
        in the alignment file exactly (case-sensitive).
    start_pos : int
        Starting genomic position for analysis (0-based coordinate system).
        Must be non-negative and less than contig length.
    end_pos : int
        Ending genomic position for analysis (0-based, exclusive).
        Must be greater than start_pos and within contig bounds.
    min_coverage : int, optional
        Minimum number of reads required to cover a position for variant calling.
        Positions with fewer reads are skipped to avoid low-confidence calls.
        Default is 5 reads.
    min_base_quality : int, optional
        Minimum base quality score (Phred scale) required for reads to contribute
        to variant calling. Low-quality bases are filtered out. Default is 10.
    max_major_allele_freq : float, optional
        Maximum frequency allowed for the most common allele at a position.
        Positions where the major allele exceeds this frequency are considered
        monomorphic and excluded. Default is 0.95 (95%).
        
    Returns
    -------
    Tuple[Dict[int, Dict[str, int]], int]
        - Nested dictionary mapping genomic positions to read variant calls:
          - Outer key: genomic position (int)
          - Inner key: read name (str)  
          - Inner value: allele call (int)
            - 0: read has minor allele (variant)
            - 1: read has major allele (reference-like)
            - Missing: read doesn't cover position or filtered out
          
          Only positions with evidence for multiple alleles are included.
        - int
          - The size of the contig in base pairs
          
    Algorithm Details
    -----------------
    1. **Pileup Generation**: 
       - Use pysam pileup with specified quality and coverage filters
       - Process each position in the specified genomic window
       
    2. **Allele Frequency Analysis**:
       - Extract all base calls at each position
       - Calculate frequency of each unique base
       - Sort alleles by decreasing frequency
       
    3. **Variant Position Detection**:
       - Skip positions where major allele frequency ≥ max_major_allele_freq
       - Require at least 2 distinct alleles with sufficient reads
       
    4. **Read Assignment**:
       - Assign binary codes based on allele identity:
         * 1: read carries the least frequent (minor) allele
         * 0: read carries the second most frequent allele  
         * Filter: reads with other alleles or quality issues
    
    Raises
    ------
    ValueError
        If start_pos >= end_pos or coordinates are invalid
    KeyError  
        If contig_name is not found in alignment file header
    IOError
        If alignment file is not properly indexed or accessible
    MemoryError
        If region is too large for available system memory
        
    See Also
    --------
    create_matrix : Converts variant dictionary to numerical matrix format
    pre_processing : Downstream processing of variant matrices  
    pysam.AlignmentFile.pileup : Underlying pileup generation function
    
    """
    # Input validation
    if start_pos >= end_pos:
        raise ValueError(f"start_pos ({start_pos}) must be less than end_pos ({end_pos})")
    
    if start_pos < 0:
        raise ValueError(f"start_pos ({start_pos}) must be non-negative")
        
    # Check if contig exists in alignment file
    if not hasattr(input_file, 'references') or input_file.references is None:
        raise ValueError("Alignment file does not contain reference information")
        
    if contig_name not in input_file.references:
        available_contigs = list(input_file.references)
        if len(available_contigs) == 0:
            raise ValueError("No contigs found in alignment file")
        
        # Show first 5 contigs in error message
        contig_preview = ', '.join(available_contigs[:5])
        if len(available_contigs) > 5:
            contig_preview += f"... (and {len(available_contigs) - 5} more)"
            
        raise ValueError(f"Contig '{contig_name}' not found in alignment file. "
                        f"Available contigs: {contig_preview}")
    
    
    # Initialize results dictionary
    suspicious_positions = {}
    
    
    try:
        # Generate pileup over specified region with quality filters
        pileup_iter = input_file.pileup(
            contig=contig_name,
            start=start_pos, 
            stop=end_pos,
            truncate=True,
            min_base_quality=min_base_quality
        )
        
        # Process each position in the pileup
        for pileupcolumn in pileup_iter:
            if pileupcolumn.nsegments >=5:
                
                tmp_dict = {}
                sequence = np.char.upper(np.array(pileupcolumn.get_query_sequences()))
                bases, freq = np.unique(np.array(sequence),return_counts=True)
        
                if bases[0] =='':    
                    ratio = freq[1:]/sum(freq[1:])
                    bases = bases[1:]
                else:
                    ratio = freq/sum(freq)
                
                if len(ratio) >0:
                    idx_sort = np.argsort(ratio)
                    
                    if ratio[idx_sort[-1]]<0.95:
                        for pileupread in pileupcolumn.pileups:        
                            if not pileupread.is_del and not pileupread.is_refskip:
                                # query position is None if is_del or is_refskip is set.
                                if pileupread.alignment.query_sequence[pileupread.query_position] == bases[idx_sort[-1]]:
                                    tmp_dict[pileupread.alignment.query_name] =  1
                                elif pileupread.alignment.query_sequence[pileupread.query_position] == bases[idx_sort[-2]]:
                                    tmp_dict[pileupread.alignment.query_name] =  0
                                else :
                                    tmp_dict[pileupread.alignment.query_name] = np.nan
                    
                        suspicious_positions[pileupcolumn.reference_pos] = tmp_dict        

        
                   
    except Exception as e:
        raise RuntimeError(f"Pileup processing failed: {e}") from e
    
    contig_size = get_contig_size(input_file, contig_name)
   
    return suspicious_positions, contig_size

def get_contig_size(input_file: ps.AlignmentFile, contig_name: str) -> int:
    """
    Get the size of a specific contig from the alignment file.
    
    Parameters
    ----------
    input_file : ps.AlignmentFile
        Opened pysam AlignmentFile object
    contig_name : str
        Name of the contig
        
    Returns
    -------
    int
        Size of the contig in base pairs
    """
    try:
        contig_index = input_file.references.index(contig_name)
        return input_file.lengths[contig_index]
    except (ValueError, IndexError):
        return 0