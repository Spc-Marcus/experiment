#!/bin/bash
# Extraire automatiquement les premiers contigs du BAM

INPUT_BAM="6.bam"
OUTPUT_DIR="output_${NB_CONTIGS}_contigs"
OUTPUT_BAM="$OUTPUT_DIR/filtered.bam"
OUTPUT_BAI="$OUTPUT_DIR/filtered.bam.bai"
OUTPUT_READS_FASTQ="$OUTPUT_DIR/reads.fastq"
OUTPUT_TARGET_FASTA="$OUTPUT_DIR/target.fasta"
OUTPUT_GFA="$OUTPUT_DIR/target.gfa"

NB_CONTIGS=10
# Nettoyer et créer le répertoire de sortie
rm -rf $OUTPUT_DIR
mkdir -p $OUTPUT_DIR

echo "=== Vérification des fichiers d'entrée ==="

# Vérifier si le fichier BAM existe
if [ ! -f "$INPUT_BAM" ]; then
    echo "Erreur: Le fichier $INPUT_BAM n'existe pas"
    exit 1
fi

# Vérifier si l'index existe, sinon le créer
if [ ! -f "${INPUT_BAM}.bai" ]; then
    echo "Création de l'index pour $INPUT_BAM..."
    samtools index $INPUT_BAM
fi

echo "=== Visualisation des contigs disponibles ==="
echo "Contigs disponibles dans le BAM:"
samtools view -H $INPUT_BAM | grep "^@SQ" | head -10 | cut -f2 | sed 's/SN://'

echo "=== Extraction des $NB_CONTIGS premiers contigs du BAM ==="

# Récupérer les premiers noms de contigs directement du header BAM
SELECTED_CONTIGS=$(samtools view -H $INPUT_BAM | grep "^@SQ" | head -$NB_CONTIGS | cut -f2 | sed 's/SN://')

echo "Contigs sélectionnés:"
echo "$SELECTED_CONTIGS"

# Convertir en array pour samtools
CONTIG_ARRAY=($SELECTED_CONTIGS)

# 1. Extraire les contigs du fichier BAM
echo "=== Création du fichier BAM filtré ==="

# Créer une liste de régions pour samtools
REGIONS=""
for contig in "${CONTIG_ARRAY[@]}"; do
    if [ -z "$REGIONS" ]; then
        REGIONS="$contig"
    else
        REGIONS="$REGIONS $contig"
    fi
done

echo "Régions à extraire: $REGIONS"
samtools view -b $INPUT_BAM $REGIONS > $OUTPUT_BAM


# 2. Créer l'index BAI
echo "=== Création de l'index BAI ==="
if [ -s "$OUTPUT_BAM" ]; then
    samtools index $OUTPUT_BAM
    echo "Index BAI créé: $OUTPUT_BAI"
else
    echo "Impossible de créer l'index: le fichier BAM est vide"
fi

# 3. Extraire les reads en format FASTQ
echo "=== Création du fichier reads.fastq ==="
if [ -s "$OUTPUT_BAM" ]; then
    samtools fastq $OUTPUT_BAM > $OUTPUT_READS_FASTQ
else
    touch $OUTPUT_READS_FASTQ
fi

# 4. Créer le fichier target.fasta avec l'allèle majoritaire
echo "=== Diagnostic de couverture par contig ==="
for contig in "${CONTIG_ARRAY[@]}"; do
    read_count=$(samtools view -c "$INPUT_BAM" "$contig" 2>/dev/null || echo "0")
    if [ "$read_count" -eq 0 ]; then
        echo "  ⚠️  $contig: AUCUN READ mappé"
    else
        echo "  ✅ $contig: $read_count reads mappés"
    fi
done

echo "=== Création du fichier target.fasta (méthode rapide) ==="

# Méthode 1: Essayer samtools consensus (le plus rapide)
if samtools consensus --help >/dev/null 2>&1; then
    echo "Utilisation de samtools consensus..."
    {
        for contig in "${CONTIG_ARRAY[@]}"; do
            echo ">$contig"
            samtools consensus -f fasta "$INPUT_BAM" -r "$contig" 2>/dev/null | grep -v "^>" | tr -d '\n' || {
                # Si échec, utiliser une séquence de N basée sur la longueur du contig
                length=$(samtools view -H $INPUT_BAM | grep "^@SQ" | grep "SN:$contig" | cut -f3 | sed 's/LN://')
                python3 -c "print('N' * int('$length' if '$length' else '1000'))"
            }
            echo ""
        done
    } > $OUTPUT_TARGET_FASTA
    
# Méthode 2: Approche simplifiée avec mpileup
else
    echo "samtools consensus non disponible, utilisation de mpileup simplifié..."
    > $OUTPUT_TARGET_FASTA  # Vider le fichier
    
    for contig in "${CONTIG_ARRAY[@]}"; do
        echo "Traitement du contig: $contig"
        echo ">$contig" >> $OUTPUT_TARGET_FASTA
        
        # Vérifier d'abord si des reads mappent sur ce contig
        read_count=$(samtools view -c "$INPUT_BAM" "$contig" 2>/dev/null || echo "0")
        
        # Obtenir la longueur du contig
        length=$(samtools view -H $INPUT_BAM | grep "^@SQ" | grep "SN:$contig" | cut -f3 | sed 's/LN://')
        echo "  Longueur attendue: $length bases"
        echo "  Reads mappés: $read_count"
        
        if [ -z "$length" ] || [ "$length" -eq 0 ]; then
            echo "  ERREUR: Impossible de déterminer la longueur du contig"
            echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> $OUTPUT_TARGET_FASTA
            continue
        fi
        
        if [ "$read_count" -eq 0 ]; then
            echo "  ATTENTION: Aucun read mappé, génération d'une séquence de N"
            python3 -c "
seq = 'N' * $length
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" >> $OUTPUT_TARGET_FASTA
            continue
        fi
        
        # Créer un fichier temporaire pour stocker les positions couvertes
        temp_coverage="/tmp/coverage_${contig}.txt"
        
        # Extraire les positions couvertes avec mpileup
        samtools mpileup -r "$contig" "$INPUT_BAM" 2>/dev/null | \
        awk -v contig_length="$length" '
        {
            pos = $2
            depth = $4
            bases = toupper($5)
            
            if(depth >= 1) {
                # Nettoyer la colonne des bases
                gsub(/\^./, "", bases)  # Supprimer les marqueurs de début
                gsub(/\$/, "", bases)   # Supprimer les marqueurs de fin
                gsub(/[+-][0-9]+[ACGTNacgtn]*/, "", bases)  # Supprimer les indels
                
                # Compter les bases valides
                a_count = gsub(/A/, "", bases)
                t_count = gsub(/T/, "", bases) 
                c_count = gsub(/C/, "", bases)
                g_count = gsub(/G/, "", bases)
                
                # Trouver la base majoritaire
                max_count = 0
                majority_base = "N"
                
                if(a_count > max_count) { max_count = a_count; majority_base = "A" }
                if(t_count > max_count) { max_count = t_count; majority_base = "T" }
                if(c_count > max_count) { max_count = c_count; majority_base = "C" }
                if(g_count > max_count) { max_count = g_count; majority_base = "G" }
                
                print pos "\t" majority_base
            }
        }' > $temp_coverage
        
        # Générer la séquence complète en Python
        python3 << EOF >> $OUTPUT_TARGET_FASTA
import sys

# Lire les positions couvertes
covered_positions = {}
try:
    with open('$temp_coverage', 'r') as f:
        for line in f:
            if line.strip():
                pos, base = line.strip().split('\t')
                covered_positions[int(pos)] = base
except FileNotFoundError:
    pass

# Créer la séquence de longueur complète
length = $length
sequence = []

for i in range(1, length + 1):
    if i in covered_positions:
        sequence.append(covered_positions[i])
    else:
        sequence.append('N')

# Imprimer la séquence par lignes de 80 caractères
full_sequence = ''.join(sequence)
for i in range(0, len(full_sequence), 80):
    print(full_sequence[i:i+80])
EOF
        
        # Nettoyer le fichier temporaire
        rm -f $temp_coverage
        
        # Vérifier la longueur générée
        actual_length=$(tail -n +2 $OUTPUT_TARGET_FASTA | grep -A 1000 "^>$contig" | head -1|grep -v "^>" | tr -d '\n' | wc -c)
        echo "  Longueur générée: $actual_length bases"
        
        if [ "$actual_length" -ne "$length" ]; then
            echo "  ATTENTION: Longueur incorrecte, régénération avec N..."
            # Supprimer les lignes incorrectes et ajouter la bonne séquence
            head -n -1 $OUTPUT_TARGET_FASTA > temp_fasta
            mv temp_fasta $OUTPUT_TARGET_FASTA
            python3 -c "
seq = 'N' * $length
for i in range(0, len(seq), 80):
    print(seq[i:i+80])
" >> $OUTPUT_TARGET_FASTA
        fi
    done
fi

echo "Fichier target.fasta créé!"

# Filtrage des contigs avec plus de 1% de N
echo "=== Filtrage des contigs avec >1% de N ==="
TEMP_FASTA="$OUTPUT_DIR/target_filtered.fasta"
> $TEMP_FASTA

FILTERED_CONTIGS=()
for contig in "${CONTIG_ARRAY[@]}"; do
    # Extraire la séquence du contig
    sequence=$(awk -v contig="$contig" '
        /^>/ { 
            if($0 == ">"contig) found=1; else found=0; next 
        } 
        found { 
            gsub(/\n/, ""); seq=seq$0 
        } 
        END { print seq }
    ' $OUTPUT_TARGET_FASTA)
    
    if [ -n "$sequence" ] && [ "$sequence" != "" ]; then
        # Calculer le pourcentage de N
        total_length=${#sequence}
        n_count=$(echo "$sequence" | grep -o "N" | wc -l)
        
        if [ "$total_length" -gt 0 ]; then
            n_percentage=$(echo "scale=2; $n_count * 100 / $total_length" | bc -l)
            
            echo "  $contig: ${n_count}N/${total_length}bp = ${n_percentage}%"
            
            # Comparer directement avec bc (plus précis que la conversion en entier)
            is_acceptable=$(echo "$n_percentage <= 1" | bc -l)
            
            if [ "$is_acceptable" -eq 1 ]; then
                echo "    ✅ Contig conservé"
                echo ">$contig" >> $TEMP_FASTA
                echo "$sequence" | fold -w 80 >> $TEMP_FASTA
                FILTERED_CONTIGS+=("$contig")
            else
                echo "    ❌ Contig supprimé (${n_percentage}% N > 1%)"
            fi
        else
            echo "  $contig: séquence vide - supprimé"
        fi
    else
        echo "  $contig: pas de séquence trouvée - supprimé"
    fi
done

# Remplacer le fichier original par le filtré
mv $TEMP_FASTA $OUTPUT_TARGET_FASTA

# Mettre à jour la liste des contigs
CONTIG_ARRAY=("${FILTERED_CONTIGS[@]}")

echo ""
echo "Contigs conservés après filtrage: ${#CONTIG_ARRAY[@]}"
for contig in "${CONTIG_ARRAY[@]}"; do
    echo "  - $contig"
done

# Régénérer le GFA avec les contigs filtrés
echo "=== Régénération du fichier GFA avec contigs filtrés ==="
{
    echo "H	VN:Z:1.0"
    for contig in "${CONTIG_ARRAY[@]}"; do
        # Extraire directement la séquence du fichier target.fasta filtré
        sequence=$(awk -v contig="$contig" '
            /^>/ { 
                if($0 == ">"contig) found=1; else found=0; next 
            } 
            found { 
                gsub(/\n/, ""); seq=seq$0 
            } 
            END { print seq }
        ' $OUTPUT_TARGET_FASTA)
        
        if [ -n "$sequence" ] && [ "$sequence" != "" ]; then
            echo "S	$contig	$sequence"
        fi
    done
} > $OUTPUT_GFA

# 6. Afficher les statistiques
echo "=== Statistiques des fichiers générés ==="
echo "Contigs conservés après filtrage:"
for i in "${!CONTIG_ARRAY[@]}"; do
    echo "  $((i+1)). ${CONTIG_ARRAY[$i]}"
done

echo ""
if [ -s "$OUTPUT_BAM" ]; then
    echo "BAM filtré:"
    samtools flagstat $OUTPUT_BAM
    
    echo ""
    echo "Détail des contigs dans le BAM filtré:"
    samtools idxstats $OUTPUT_BAM
else
    echo "Le fichier BAM est vide - aucun read trouvé pour ces contigs"
fi

echo ""
echo "Fichiers générés dans $OUTPUT_DIR:"
ls -lh $OUTPUT_DIR/

# 7. Vérifications
echo "=== Vérifications ==="
echo "Index BAI présent: $([ -f "$OUTPUT_BAI" ] && echo "OUI" || echo "NON")"
echo "Taille du BAM filtré: $(du -h $OUTPUT_BAM 2>/dev/null | cut -f1 || echo "0")"
echo "Nombre de reads dans reads.fastq: $(grep -c "^@" $OUTPUT_READS_FASTQ 2>/dev/null || echo "0")"
echo "Nombre de contigs dans target.fasta: $(grep -c "^>" $OUTPUT_TARGET_FASTA 2>/dev/null || echo "0")"
echo "Lignes dans le GFA: $(wc -l < $OUTPUT_GFA 2>/dev/null || echo "0")"

echo ""
echo "=== Aperçu des séquences générées ==="
for contig in "${CONTIG_ARRAY[@]}"; do
    echo "Contig $contig (premiers 100 caractères):"
    sequence_preview=$(awk -v contig="$contig" '
        /^>/ { 
            if($0 == ">"contig) found=1; else found=0; next 
        } 
        found { 
            printf "%s", $0 
        }
    ' $OUTPUT_TARGET_FASTA | head -c 100)
    
    if [ -n "$sequence_preview" ]; then
        echo "$sequence_preview..."
    else
        echo "Pas de séquence trouvée"
    fi
    
    # Afficher aussi les statistiques de bases
    base_stats=$(awk -v contig="$contig" '
        /^>/ { 
            if($0 == ">"contig) found=1; else found=0; next 
        } 
        found { 
            seq=seq$0 
        } 
        END { 
            a=gsub(/A/,"",seq); t=gsub(/T/,"",seq); c=gsub(/C/,"",seq); g=gsub(/G/,"",seq); n=gsub(/N/,"",seq)
            printf "A:%d T:%d C:%d G:%d N:%d", a, t, c, g, n
        }
    ' $OUTPUT_TARGET_FASTA)
    echo "  Composition: $base_stats"
    echo ""
done

echo ""
echo "=== Extraction terminée avec succès ==="
echo "Tous les fichiers ont été créés dans: $OUTPUT_DIR"
