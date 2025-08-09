# Efficience_ALL – Pré/ILP/Post-Processing sur matrices binaires

Suite d’expérimentations pour évaluer l’impact combiné:
- du prétraitement (biclustering/steps),
- d’un enrichissement ILP (détection de quasi-bicliques),
- et du post-traitement (fusion/filtrage/réaffectation)
sur des matrices binaires synthétiques (0/1).

Le module génère des matrices, exécute le prétraitement, ajoute des steps via ILP, applique le post-traitement avec plusieurs seuils de distance et journalise les métriques dans un CSV.

## Fonctionnalités

- Génération de matrices binaires synthétiques avec:
  - Extension aléatoire des lignes/colonnes
  - Injection de bruit (taux d’erreur)
  - Permutation aléatoire des lignes/colonnes
- Prétraitement: calcule des étapes (steps) de séparation de lectures
- ILP (Gurobi): extraction de patterns denses pour proposer des steps additionnels
- Post-traitement:
  - Application séquentielle des steps
  - Fusion (distance de Hamming, AgglomerativeClustering)
  - Réaffectation des orphelins
  - Filtrage par taille minimale
  - Construction d’une matrice réduite (profils consensus par step)
- Évaluation multi-seuils de distance de fusion
- Export CSV (comparatif pré seul vs pré+ILP)

## Arborescence

- `test_all.py` : script principal (grille + CSV comparatif)
- `create_matrix.py` : création/extension/bruit/mélange des matrices
- `pre_processing.py` : prétraitement (externe, requis)
- `post_processing.py` : post-traitement (fusion/réaffectation/filtrage)
- `ilp.py` : génération de steps via ILP (Gurobi)
- `results.csv` : sortie des expériences

## Installation

Dépendances Python:
```bash
pip install numpy scikit-learn gurobipy
```

Prétraitement attendu (fichier `pre_processing.py` au même niveau) doit exposer:
```python
inhomogeneous_regions, steps = pre_processing(
    matrix, min_col_quality: int, certitude: float, error_rate: float
)
```

ILP:
- `ilp.py` appelle `ilp_grb.find_quasi_dens_matrix_max_ones` (Gurobi requis).

## Utilisation

### 1) Exécution rapide (valeurs par défaut)
Lance la grille et écrit `results.csv`:
```bash
python test_all.py
```

Paramètres par défaut (dans le main):
- nb_matrix_permutations = 25
- min_rows = 3, min_cols = 3, max_rows = 20, max_cols = 12
- thresholds = [0, 0.01, ..., 0.2]
- error_rates = [0.0, 0.005, ..., 0.1]
- strips = [3..9]
- haplotypes = [4..10]
- distance_thresh = [0.0, 0.05, 0.1, 0.15]
- csv_file = "results.csv"

Taille de la grille:
len(thresholds) × len(error_rates) × len(strips) × len(haplotypes) × nb_matrix_permutations × len(distance_thresh)

### 2) Appel programmatique
```python
from test_all import test_all

test_all(
    thresholds=[0.02, 0.05],
    error_rates=[0.0, 0.02],
    strips=[4, 5],
    haplotypes=[4, 6],
    nb_matrix_permutations=5,
    min_rows=3, min_cols=3, max_rows=20, max_cols=12,
    csv_file="results.csv",
    distance_thresh=[0.05, 0.1]
)
```

### 3) Reproductibilité
```python
import random, numpy as np
random.seed(42)
np.random.seed(42)
```

## Détails du pipeline

- Génération: `create_simple_matrix` → `extend_matrix` → `add_noise_to_matrix` → `mix_matrices`
- Prétraitement: `pre_processing(matrix_data, min_col_quality=min_cols, certitude=threshold, error_rate=error_rate)`
- ILP: `clustering_full_matrix(matrix_data, inhomogeneous_regions, min_rows, min_cols, error_rate)`
  - Retourne des steps ILP s’ils respectent les seuils (colonnes/qualité)
- Post-traitement: `post_processing(matrix, steps, read_names, distance_thresh, min_reads_per_cluster=min_rows)`
  - Exécuté pour:
    - steps prétraitement seul
    - steps prétraitement + ILP
  - Pour chaque distance dans `distance_thresh`

## Sortie CSV

En-tête (dans `test_all.py`):
```
Threshold,Error-Rate,Strip,Haplotype,Matrix-Size-cols,Matrix-Size-rows,Time-Pre-Processing,Steps-Count-Pre,Unused-Cols-Pre,Final-Clusters-Pre,Orphan-Reads-Pre,Steps-Count-ILP,Final-Clusters-ILP,Orphan-Reads-ILP,Distance
```

- Threshold: seuil de certitude (prétraitement)
- Error-Rate: bruit injecté
- Strip / Haplotype: colonnes/lignes de base avant extension
- Matrix-Size-cols / rows: dimensions finales après extension/bruit/mélange
- Time-Pre-Processing: temps du prétraitement
- Steps-Count-Pre: nombre d’étapes prétraitement
- Unused-Cols-Pre: colonnes inutilisées (comptées après post-traitement; voir code)
- Final-Clusters-Pre / Orphan-Reads-Pre: résultats post-traitement (pré seul)
- Steps-Count-ILP: nombre d’étapes ILP ajoutées
- Final-Clusters-ILP / Orphan-Reads-ILP: résultats post-traitement (pré+ILP)
- Distance: seuil de fusion (Hamming) utilisé

## Paramètres clés

- distance_thresh (liste): seuils de fusion, ex. [0.0, 0.05, 0.1, 0.15]
- min_reads_per_cluster: fixé à `min_rows` dans `post_processing`
- ILP:
  - `min_row_quality` = `min_rows`
  - `min_col_quality` = `min_cols`
  - `error_rate`: identique au taux de bruit

## Dépannage

- Erreur "Limite fondamentale dépassée":
  - Cause: `haplotypes > 2^strips` dans `create_simple_matrix`
  - Solutions: augmenter `strip` ou diminuer `haplotype`
- Gurobi/licence:
  - Installer `gurobipy` et configurer la licence (académique si possible)
- Exécutions longues:
  - Réduire les listes de la grille, `nb_matrix_permutations` et/ou `distance_thresh`

