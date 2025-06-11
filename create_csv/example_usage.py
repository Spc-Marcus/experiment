import pysam as ps
from pathlib import Path
import json
from core import get_data
from matrix import create_matrix

def process_genomic_data(bam_file_path: str, output_dir: str, 
                        contig_name: str, start_pos: int, end_pos: int,
                        filtered_row_threshold: float = 0.6):
    """
    Example function showing how to use the updated get_data and create_matrix functions.
    """
    # Create output directory
    tmp_dir = Path(output_dir)
    tmp_dir.mkdir(exist_ok=True)
    
    # Open BAM file
    with ps.AlignmentFile(bam_file_path, "rb") as file:
        # Extract variant data (now returns contig size as well)
        dict_of_sus_pos, contig_size = get_data(file, contig_name, start_pos, end_pos)
        total_variants = len(dict_of_sus_pos)
        
        print(f"Found {total_variants} variant positions")
        print(f"Contig {contig_name} size: {contig_size} bp")
        
        # Save intermediate data
        output_file = tmp_dir / f"dict_of_sus_pos_{contig_name}_{start_pos}.json"
        with open(output_file, "w", encoding="utf-8") as f:
            # Convert numpy types to native Python types for JSON serialization
            serializable_dict = {}
            for pos, reads_dict in dict_of_sus_pos.items():
                serializable_dict[str(pos)] = {}
                for read_name, value in reads_dict.items():
                    if pd.isna(value):
                        serializable_dict[str(pos)][read_name] = None
                    else:
                        serializable_dict[str(pos)][read_name] = int(value)
            json.dump(serializable_dict, f, ensure_ascii=False, indent=2)
        
        haplotypes_here = {}
        
        if dict_of_sus_pos:
            # Create and filter matrix using create_matrix function with CSV saving
            X_matrix, reads = create_matrix(
                dict_of_sus_pos, 
                filtered_row_threshold,
                contig_size=contig_size,
                contig_name=contig_name,
                start_pos=start_pos,
                csv_output_dir=tmp_dir
            )
            
            print(f"Matrix shape: {X_matrix.shape}")
            print(f"Number of reads: {len(reads)}")
            
            # Check if CSV was saved
            if (contig_size > 90000 and 
                X_matrix.size > 0 and 
                X_matrix.shape[0] >= 25 and 
                X_matrix.shape[1] >= 25):
                print(f"CSV matrix saved for large contig ({contig_size} bp) with matrix size {X_matrix.shape}")
            else:
                print("CSV not saved - conditions not met")
                print(f"  - Contig size > 90k: {contig_size > 90000}")
                if X_matrix.size > 0:
                    print(f"  - Matrix >= 25x25: {X_matrix.shape[0] >= 25 and X_matrix.shape[1] >= 25}")
                else:
                    print("  - Matrix is empty")
        
        return X_matrix, reads, contig_size

# Example usage
if __name__ == "__main__":
    # Replace with your actual file paths and parameters
    bam_file = "/path/to/your/file.bam"
    output_directory = "/udd/mfoin/Dev/experiment/create_csv/output"
    contig = "chr1"  # or whatever contig you want to analyze
    start = 1000000
    end = 1100000
    
    X_matrix, reads, contig_size = process_genomic_data(
        bam_file, 
        output_directory, 
        contig, 
        start, 
        end
    )
