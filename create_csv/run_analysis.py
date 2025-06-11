import argparse
import sys
from example_usage import process_genomic_data

def main():
    parser = argparse.ArgumentParser(description='Process genomic data and create CSV matrices')
    parser.add_argument('bam_file', help='Path to BAM file')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('contig', help='Contig name')
    parser.add_argument('start', type=int, help='Start position')
    parser.add_argument('end', type=int, help='End position')
    parser.add_argument('--threshold', type=float, default=0.6, help='Minimum read coverage threshold: reads must cover at least this fraction of variant positions (default: 0.6)')
    
    args = parser.parse_args()
    
    try:
        X_matrix, reads, contig_size = process_genomic_data(
            args.bam_file,
            args.output_dir,
            args.contig,
            args.start,
            args.end,
            args.threshold
        )
        print("Analysis completed successfully!")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
