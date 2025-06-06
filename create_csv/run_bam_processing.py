import argparse
import logging
from pathlib import Path
from bam_to_matrix import process_bam_files, process_single_bam

def main():
    parser = argparse.ArgumentParser(description='Process BAM files to create CSV matrices')
    parser.add_argument('--bam-dir', default='/bam', 
                       help='Directory containing BAM files (default: /bam)')
    parser.add_argument('--output-dir', default='matrices',
                       help='Output directory for matrices (default: matrices)')
    parser.add_argument('--min-contig-length', type=int, default=90000,
                       help='Minimum contig length to process (default: 90000)')
    parser.add_argument('--min-matrix-size', type=int, default=25,
                       help='Minimum matrix size (rows and cols) to save (default: 25)')
    parser.add_argument('--single-bam', type=str,
                       help='Process a single BAM file instead of directory')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    logger = logging.getLogger(__name__)
    
    try:
        if args.single_bam:
            # Process single BAM file
            logger.info(f"Processing single BAM file: {args.single_bam}")
            created_files = process_single_bam(
                args.single_bam,
                args.output_dir,
                args.min_contig_length,
                args.min_matrix_size
            )
            
            if created_files:
                logger.info(f"Successfully created {len(created_files)} matrices:")
                for file in created_files:
                    logger.info(f"  - {file}")
            else:
                logger.warning("No matrices were created")
                
        else:
            # Process all BAM files in directory
            logger.info(f"Processing BAM files in directory: {args.bam_dir}")
            results = process_bam_files(
                args.bam_dir,
                args.output_dir,
                args.min_contig_length,
                args.min_matrix_size
            )
            
            # Print detailed summary
            total_matrices = 0
            successful_bams = 0
            
            for bam_file, matrices in results.items():
                if matrices:
                    successful_bams += 1
                    total_matrices += len(matrices)
                    logger.info(f"{Path(bam_file).name}: {len(matrices)} matrices created")
                else:
                    logger.warning(f"{Path(bam_file).name}: No matrices created")
            
            logger.info(f"Summary: {total_matrices} total matrices from {successful_bams}/{len(results)} BAM files")
            
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
