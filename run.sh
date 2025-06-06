#!/bin/bash
# filepath: /udd/mfoin/Dev/experiment/index_bam_files.sh

# Script to index BAM files that don't have corresponding .bai index files
# Usage: ./index_bam_files.sh

# Change to the bam directory
cd "$(dirname "$0")/bam" || {
    echo "Error: Could not change to bam directory"
    exit 1
}

echo "Checking BAM files in: $(pwd)"

# Initialize counters
indexed_count=0
skipped_count=0
error_count=0

# Process all BAM files
for bam_file in *.bam; do
    # Check if any BAM files exist
    if [ "$bam_file" = "*.bam" ]; then
        echo "No BAM files found in $(pwd)"
        exit 0
    fi
    
    # Check if corresponding .bai file exists
    bai_file="${bam_file}.bai"
    
    if [ -f "$bai_file" ]; then
        echo "✓ $bam_file already indexed (${bai_file} exists)"
        ((skipped_count++))
    else
        echo "⏳ Indexing $bam_file..."
        
        # Run samtools index
        if samtools index "$bam_file"; then
            echo "✅ Successfully indexed $bam_file"
            ((indexed_count++))
        else
            echo "❌ Failed to index $bam_file"
            ((error_count++))
        fi
    fi
done

# Print summary
echo ""
echo "=== Indexing Summary ==="
echo "Files indexed: $indexed_count"
echo "Files skipped (already indexed): $skipped_count"
echo "Errors: $error_count"

if [ $error_count -eq 0 ]; then
    echo "✅ All BAM files are now indexed!"
    cd ..
    # Change back to parent directory for Python script
    cd "$(dirname "$0")" || {
        echo "Error: Could not change to parent directory"
        exit 1
    }
    
    echo ""
    echo "=== Running Matrix Generation ==="
    echo "⏳ Processing BAM files to create CSV matrices..."
    
    # Run the Python script to process BAM files
    if python create_csv/process_bam_folder.py bam .; then
        echo "✅ Matrix generation completed successfully!"
        
        echo ""
        echo "=== Running ilphaplo ==="
        echo "⏳ Running ilphaplo analysis..."
        
        # Run ilphaplo batch script
        if cd ilphaplo && bash run_batch.sh; then
            echo "✅ ilphaplo analysis completed successfully!"
            echo ""
            echo "🎉 All processing steps completed successfully!"
            exit 0
        else
            echo "❌ ilphaplo analysis failed"
            echo "Check if ilphaplo directory exists and contains run_batch.sh"
            exit 1
        fi
    else
        echo "❌ Matrix generation failed"
        exit 1
    fi
else
    echo "⚠️  Some files failed to index. Check samtools installation and file permissions."
    exit 1
fi