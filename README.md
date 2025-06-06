# Experiment Runner

## Structure des Données

Organisez vos matrices CSV dans la structure suivante :
```
matrice/
├── 2/
│   ├── small_matrix_2_1.csv      (20x10)
│   ├── medium_matrix_2_2.csv     (50x25)
│   └── large_matrix_2_3.csv      (80x40)
├── 3/
│   ├── matrix_3_1.csv            (30x15)
│   └── matrix_3_2.csv            (60x30)
├── 4/
│   ├── matrix_4_1.csv            (40x20)
│   └── matrix_4_2.csv            (70x35)
├── 5/
│   ├── matrix_5_1.csv            (50x25)
│   └── matrix_5_2.csv            (90x45)
└── 6/
    ├── matrix_6_1.csv            (60x30)
    └── matrix_6_2.csv            (100x50)
```

Où le nombre (2, 3, 4, 5, 6) représente le nombre d'haplotypes.

## Caractéristiques des Matrices d'Échantillon

- **2 haplotypes** : 3 matrices de complexité croissante (20x10 à 80x40)
- **3 haplotypes** : 2 matrices de taille moyenne (30x15, 60x30)  
- **4 haplotypes** : 2 matrices de complexité modérée (40x20, 70x35)
- **5 haplotypes** : 2 matrices de haute complexité (50x25, 90x45)
- **6 haplotypes** : 2 matrices très complexes (60x30, 100x50)

Les densités varient de 0.3 à 0.5 pour simuler différents niveaux de variabilité génétique.

## Comment Lancer les Expériences

### Option 1: Script Python Direct

```bash
cd /udd/mfoin/Dev/experiment/ilphaplo
python run_experiments.py
```

### Option 2: Script Bash (recommandé)

```bash
cd /udd/mfoin/Dev/experiment/ilphaplo
chmod +x run_batch.sh
./run_batch.sh
```

## Métriques d'Analyse

Les expériences analysent :

- **Scalabilité** : Performance vs taille de matrice
- **Complexité** : Impact du nombre d'haplotypes  
- **Efficacité** : Ratio appels ILP / patterns trouvés
- **Temps d'exécution** : Décomposition par phase (preprocessing, clustering)
- **Qualité des résultats** : Nombre de patterns et régions découverts

## Résultats Attendus

Pour chaque matrice, vous obtiendrez :

1. **Métriques de performance** :
   - Nombre d'appels ILP (augmente avec la complexité)
   - Temps d'exécution total et par phase
   - Ratio temps ILP / temps total

2. **Qualité de l'analyse** :
   - Nombre de patterns découverts
   - Nombre de régions identifiées
   - Taille moyenne des régions

3. **Caractéristiques de la matrice** :
   - Dimensions et densité
   - Statistiques sur les lignes/colonnes

## Exemple de Sortie Attendu

```
Found 12 CSV files to process
Processing 1/12: matrice/2/small_matrix_2_1.csv
  -> Completed: 5 ILP calls, 2 patterns, 0.234s
Processing 2/12: matrice/2/medium_matrix_2_2.csv  
  -> Completed: 12 ILP calls, 3 patterns, 0.567s
...
Processing 12/12: matrice/6/matrix_6_2.csv
  -> Completed: 45 ILP calls, 8 patterns, 2.134s

=== EXPERIMENT SUMMARY ===
Total matrices processed: 12
Successful runs: 12
Failed runs: 0

=== BY HAPLOTYPE COUNT ===
Haplotypes 2: 3 matrices, avg ILP calls: 8.7, avg time: 0.45s
Haplotypes 3: 2 matrices, avg ILP calls: 15.5, avg time: 0.89s
Haplotypes 4: 2 matrices, avg ILP calls: 22.0, avg time: 1.23s
Haplotypes 5: 2 matrices, avg ILP calls: 31.5, avg time: 1.67s
Haplotypes 6: 2 matrices, avg ILP calls: 42.0, avg time: 2.01s
```
