# Contig Extraction and Reference Generation

This package provides automated tools for extracting contigs from BAM files and generating reference sequences using majority allele consensus. It's designed for rapid processing of alignment data to create target references for downstream phylogenetic analysis.

## Overview

The package implements a streamlined pipeline for:
1. **Automatic Contig Selection**: Extracts the first N contigs from BAM headers
2. **Fast Consensus Generation**: Uses majority allele calling for reference creation
3. **Multi-format Output**: Generates BAM, FASTQ, FASTA, and GFA files
4. **Quality Control**: Validates sequence lengths and provides detailed statistics
5. **Performance Optimized**: Uses native samtools commands for maximum speed

## Key Features

- **Rapid Processing**: Optimized for speed using samtools consensus when available
- **Majority Allele Consensus**: Generates reference based on most frequent alleles
- **Multiple Output Formats**: BAM, FASTQ, FASTA, and GFA compatible with downstream tools
- **Automatic Validation**: Checks sequence lengths against expected contig sizes
- **Comprehensive Logging**: Detailed statistics and quality metrics
- **Fallback Methods**: Multiple algorithms ensure reliable sequence generation

## Installation

### Prerequisites

- Python 3.8+
- samtools (latest version recommended)
- Required Python packages:
  ```bash
  pip install pysam numpy
  ```

### Dependencies

```bash
# Install samtools (Ubuntu/Debian)
sudo apt-get install samtools

# Install samtools (macOS with Homebrew)
brew install samtools

# Install Python dependencies
pip install pysam numpy
```

## Usage

### Basic Contig Extraction

```bash
./extract_contigs.sh
```

### Input Requirements

The script expects:
- **Input BAM file**: `6.bam` (configurable in script)
- **BAM index**: `6.bam.bai` (auto-generated if missing)
- **Proper headers**: BAM must contain @SQ entries with contig information

### Command Line Execution

```bash
# Basic extraction (5 contigs by default)
./extract_contigs.sh

# Make script executable if needed
chmod +x extract_contigs.sh
```

### Configuration

Edit the script to modify default parameters:

```bash
# Input/Output configuration
INPUT_BAM="6.bam"                    # Input BAM file
OUTPUT_DIR="output_5_contigs"        # Output directory

# Processing parameters (modify in script)
NUM_CONTIGS=5                        # Number of contigs to extract
MIN_COVERAGE=3                       # Minimum read depth for consensus
MIN_BASE_QUALITY=10                  # Minimum base quality score
```

## Input Format

### BAM File Requirements

Input BAM files must contain:
- **Proper headers**: @SQ entries with sequence names and lengths
- **Aligned reads**: Mapped sequencing data
- **Index file**: .bai index (auto-generated if missing)

### Example BAM Header

```
@SQ	SN:edge_1	LN:619
@SQ	SN:edge_2	LN:32764
@SQ	SN:edge_3	LN:15423
@SQ	SN:edge_4	LN:8901
@SQ	SN:edge_5	LN:12056
```

## Output Results

### Generated Files

```
output_5_contigs/
├── filtered.bam          # Filtered BAM with selected contigs
├── filtered.bam.bai      # BAM index file
├── reads.fastq           # Extracted reads in FASTQ format
├── target.fasta          # Reference sequences (majority consensus)
└── target.gfa            # Graph format for assembly tools
```

### Example Output

```
=== Vérifications ===
Index BAI présent: OUI
Taille du BAM filtré: 8.7M
Nombre de reads dans reads.fastq: 462
Nombre de contigs dans target.fasta: 5
Lignes dans le GFA: 6

=== Aperçu des séquences générées ===
Contig edge_1 (premiers 100 caractères):
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT...
  Composition: A:156 T:154 C:155 G:154 N:0

Contig edge_2 (premiers 100 caractères):
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG...
  Composition: A:8190 T:8191 C:8192 G:8191 N:0
```

### File Descriptions

#### 1. filtered.bam
- **Content**: Reads mapping to selected contigs only
- **Format**: Standard BAM format
- **Usage**: Input for variant calling, visualization

#### 2. reads.fastq
- **Content**: All reads from selected contigs in FASTQ format
- **Format**: Standard FASTQ with quality scores
- **Usage**: Read-based analysis, quality assessment

#### 3. target.fasta
- **Content**: Consensus reference sequences using majority alleles
- **Format**: Multi-FASTA with contig names as headers
- **Usage**: Reference for alignment, phylogenetic analysis

#### 4. target.gfa
- **Content**: Graph assembly format compatible with assemblers
- **Format**: GFA v1.0 with sequence entries
- **Usage**: Assembly tools, graph visualization

## Algorithm Details

### Consensus Generation Methods

The script uses multiple methods for optimal performance:

#### Method 1: samtools consensus (Preferred)
```bash
samtools consensus -f fasta "$INPUT_BAM" -r "$contig"
```
- **Speed**: Fastest available method
- **Quality**: Native C implementation
- **Requirements**: Recent samtools version

#### Method 2: mpileup + majority calling
```bash
samtools mpileup -r "$contig" "$INPUT_BAM" | awk [majority_script]
```
- **Speed**: Fast, reliable fallback
- **Quality**: Custom majority allele detection
- **Requirements**: Standard samtools + awk

#### Method 3: Position-by-position analysis
- **Speed**: Slower but comprehensive
- **Quality**: Full coverage validation
- **Requirements**: Python + pysam

### Majority Allele Algorithm

```
For each position:
1. Extract all bases covering the position
2. Filter by minimum base quality (≥10)
3. Count frequency of each base (A, T, C, G)
4. Select most frequent base as consensus
5. Use 'N' if no base meets minimum coverage (≥3 reads)
```

### Quality Control

- **Length validation**: Ensures generated sequences match expected contig lengths
- **Coverage assessment**: Reports positions with insufficient coverage
- **Base composition**: Analyzes A/T/C/G/N distribution
- **Density metrics**: Calculates information content

## Examples

### Example : Basic Processing

Process default BAM file with 5 contigs:

```bash
./extract_contigs.sh
```

**Expected Output:**
```
=== Vérification des fichiers d'entrée ===
Fichier 6.bam trouvé

=== Extraction des 5 premiers contigs du BAM ===
Contigs sélectionnés:
edge_1
edge_2
edge_3
edge_4
edge_5

=== Création du fichier target.fasta (méthode rapide) ===
Utilisation de samtools consensus...
Fichier target.fasta créé rapidement!

=== Extraction terminée avec succès ===
Tous les fichiers ont été créés dans: output_5_contigs
```

## Configuration Parameters

### Script Variables

```bash
# Core configuration
INPUT_BAM="6.bam"                    # Input BAM file path
OUTPUT_DIR="output_5_contigs"        # Output directory name
```

### Algorithm Selection

The script automatically selects the best available method:

1. **samtools consensus** (if available) - Fastest
2. **mpileup + awk** (fallback) - Reliable 
3. **Python pysam** (emergency) - Comprehensive

## File Format Specifications

### target.fasta Format

```
>edge_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>edge_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

### target.gfa Format

```
H	VN:Z:1.0
S	edge_1	ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
S	edge_2	GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

### reads.fastq Format

```
@read_1
ATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIII
@read_2
GCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIII
```
