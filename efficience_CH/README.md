# Efficience_CH – Test Pré/Post-Processing sur matrices binaires

Suite d’expérimentation pour évaluer l’impact des paramètres de prétraitement et post-traitement sur des matrices binaires synthétiques (0/1). Le module génère des matrices, applique un prétraitement (biclustering/étapes), réalise un post-traitement (fusion/filtrage/réaffectation) et journalise les métriques de performance dans un CSV.

## Fonctionnalités

- Génération de matrices binaires synthétiques avec:
  - Extension aléatoire des lignes/colonnes
  - Injection de bruit (taux d’erreur)
  - Permutation aléatoire des lignes/colonnes
- Boucle d’expériences par grille de paramètres:
  - Seuils de certitude (thresholds)
  - Taux d’erreur (error_rates)
  - Strips (colonnes de base)
  - Haplotypes (lignes de base)
- Prétraitement: calcule des étapes (steps) de séparation de lectures
- Post-traitement:
  - Fusion de clusters similaires (distance de Hamming)
  - Réaffectation des reads orphelins
  - Filtrage par taille minimale de cluster
  - Construction d’une matrice réduite par profils de cluster
- Export CSV détaillé (temps, tailles, compte d’étapes, colonnes inutilisées, clusters finaux, orphelins)

## Arborescence

- `test_pre_post.py` : script principal d’expérimentations (grille + CSV)
- `create_matrix.py` : utilitaires de création/extension/bruit/mélange de matrices
- `pre_processing.py` : prétraitement (non inclus ici, mais requis)
- `post_processing.py` : post-traitement des étapes en clusters finaux
- `results.csv` : sortie des expériences (si lancé avec les paramètres par défaut)

## Installation

Dépendances Python:
```bash
pip install numpy scikit-learn
```

Le fichier `pre_processing.py` doit être présent dans ce répertoire et exposer une fonction:
```python
inhomogeneous_regions, steps = pre_processing(
    matrix, min_col_quality: int, certitude: float, error_rate: float
)
```

## Utilisation

### 1) Exécution rapide (valeurs par défaut)
Lance une grille d’expériences et écrit `results.csv`:
```bash
python test_pre_post.py
```

Paramètres par défaut (dans le bloc main):
- nb_matrix_permutations = 25
- min_rows = 3, min_cols = 3, max_rows = 20, max_cols = 12
- thresholds = [0, 0.01, ..., 0.2]
- error_rates = [0.0, 0.005, ..., 0.1]
- strips = [3..9]
- haplotypes = [4..10]
- csv_file = "results.csv"

Attention: la grille peut être volumineuse. Le temps total est proportionnel à:
