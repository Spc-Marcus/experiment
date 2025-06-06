# BAM to Matrix Converter

Ce module convertit les fichiers BAM en matrices CSV pour l'analyse d'haplotypes.

## Fonctionnalités

- Lit tous les fichiers `.bam` dans un dossier spécifié
- Ne traite que les contigs avec plus de 90k bases (configurable)
- Crée des matrices binaires (0/1) avec imputation des valeurs manquantes
- Sauvegarde dans `matrices/X/` où X est le nom du fichier BAM
- Ne sauvegarde que les matrices d'au moins 25x25 (configurable)

## Utilisation

### Traitement de tous les BAM d'un dossier

```bash
python run_bam_processing.py --bam-dir /path/to/bam/files
```

### Traitement d'un seul fichier BAM

```bash
python run_bam_processing.py --single-bam /path/to/file.bam
```

### Options avancées

```bash
python run_bam_processing.py \
    --bam-dir /custom/bam/path \
    --output-dir custom_matrices \
    --min-contig-length 50000 \
    --min-matrix-size 30 \
    --verbose
```

## Structure de sortie

```
matrices/
├── sample1/
│   ├── contig1_matrix.csv
│   ├── contig2_matrix.csv
│   └── ...
├── sample2/
│   ├── contig1_matrix.csv
│   └── ...
└── ...
```

## Dépendances

- pysam
- pandas
- numpy
- scikit-learn

```bash
pip install pysam pandas numpy scikit-learn
```

## Algorithme

1. **Lecture BAM** : Extraction des variants pour chaque contig ≥ 90k bp
2. **Création matrice** : Conversion en matrice read×position 
3. **Imputation** : Remplissage des valeurs manquantes (KNN)
4. **Binarisation** : Conversion en valeurs 0/1 uniquement
5. **Filtrage** : Sauvegarde seulement si ≥ 25×25
6. **Sauvegarde** : CSV dans `matrices/nom_bam/contig_matrix.csv`
