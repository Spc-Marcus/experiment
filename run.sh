#!/bin/bash
# filepath: /udd/mfoin/Dev/experiment/index_bam_files.sh

# Script to index BAM files that don't have corresponding .bai index files
# Usage: ./index_bam_files.sh (run from project root)

echo "Checking BAM files in: bam/"

# Check if bam directory exists
if [ ! -d "bam" ]; then
    echo "Error: BAM directory 'bam' not found in current directory"
    exit 1
fi

# Initialize counters
indexed_count=0
skipped_count=0
error_count=0

# Process all BAM files
for bam_file in bam/*.bam; do
    # Check if any BAM files exist
    if [ ! -e "$bam_file" ]; then
        echo "No BAM files found in bam/"
        exit 0
    fi
    
    # Get just the filename for display
    bam_filename=$(basename "$bam_file")
    
    # Check if corresponding .bai file exists
    bai_file="${bam_file}.bai"
    
    if [ -f "$bai_file" ]; then
        echo "✓ $bam_filename already indexed ($(basename "$bai_file") exists)"
        ((skipped_count++))
    else
        echo "⏳ Indexing $bam_filename..."
        
        # Run samtools index
        if samtools index "$bam_file"; then
            echo "✅ Successfully indexed $bam_filename"
            ((indexed_count++))
        else
            echo "❌ Failed to index $bam_filename"
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
        if [ -d "ilphaplo" ] && [ -f "ilphaplo/run_batch.sh" ]; then
            if bash ilphaplo/run_batch.sh; then
                echo "✅ ilphaplo analysis completed successfully!"
                echo ""
                echo "🎉 All processing steps completed successfully!"
                exit 0
            else
                echo "❌ ilphaplo analysis failed"
                exit 1
            fi
        else
            echo "❌ ilphaplo directory or run_batch.sh not found"
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