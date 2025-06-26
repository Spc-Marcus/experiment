# ILP Tracker Decorators

Ce module fournit des décorateurs pour surveiller et analyser les performances des opérations de programmation linéaire en nombres entiers (ILP) et des pipelines de traitement associés.

## Pourquoi utiliser ce module ?

### Problèmes résolus
- **Surveillance des performances** : Mesure automatique du temps d'exécution des fonctions utilisant des solveurs ILP
- **Debugging** : Identification des goulots d'étranglement dans les pipelines de traitement
- **Optimisation** : Données quantitatives pour améliorer l'efficacité des algorithmes
- **Monitoring** : Suivi des taux de succès/échec des appels aux solveurs

### Avantages
- Intégration transparente via décorateurs Python
- Thread-safe pour applications multi-threadées
- Métriques détaillées sans modification du code existant
- Rapport de performance complet

## Comment utiliser

### 1. Décorateurs de base

```python
from utils.decorators.ilp_tracker import track_ilp_call, track_preprocessing, track_clustering

@track_ilp_call(phase='optimization')
def solve_ilp_problem(matrix, constraints):
    # Votre code de résolution ILP
    # Note: Ceci compte 1 appel même si le solveur fait plusieurs optimisations internes
    return rows, cols, success_flag

@track_preprocessing()
def preprocess_data(raw_matrix):
    # Prétraitement des données
    return processed_matrix, regions, steps

@track_clustering()
def cluster_regions(matrix, regions):
    # Algorithme de clustering
    return clusters
```

### 2. Suivi de pipeline complet

```python
from utils.decorators.ilp_tracker import PipelineTracker

with PipelineTracker() as tracker:
    # Exécution du pipeline complet
    processed_data = preprocess_data(raw_data)
    clusters = cluster_regions(processed_data)
    results = solve_ilp_problem(clusters)
    
    # Récupération des métriques
    performance_report = tracker.get_results()
    print(f"Temps total: {performance_report['pipeline_total_time']:.2f}s")
```

### 3. Analyse des performances

```python
from utils.decorators.ilp_tracker import get_performance_summary, reset_tracker

# Réinitialiser les compteurs
reset_tracker()

# Après exécution de vos fonctions décorées
summary = get_performance_summary()

print(f"Appels de fonctions ILP: {summary['ilp_calls_total']}")  # Nombre d'appels de fonctions
print(f"Temps ILP moyen par appel: {summary['avg_ilp_time']:.3f}s")
print(f"Taux de succès des appels: {summary['solver_status_counts']}")
```

## API Reference

### Décorateurs

#### `@track_ilp_call(phase='unknown')`
Surveille les appels aux fonctions utilisant des solveurs ILP.
- **phase** : Nom de la phase pour identification
- **Mesure** : Temps d'exécution de la fonction, statut de succès du retour, patterns trouvés
- **Important** : Compte les appels à la fonction décorée, pas les optimisations internes du solveur

#### `@track_preprocessing()`
Surveille les opérations de prétraitement.
- **Mesure** : Temps de prétraitement, forme des matrices, régions traitées

#### `@track_clustering()`
Surveille les opérations de clustering.
- **Mesure** : Temps de clustering, nombre d'étapes

### Fonctions utilitaires

#### `get_performance_summary() -> Dict[str, Any]`
Retourne un résumé complet des performances :
```python
{
    'ilp_calls_total': int,        # Nombre d'appels aux fonctions @track_ilp_call
    'ilp_time_total': float,       # Temps total passé dans ces fonctions
    'preprocessing_time': float,
    'clustering_time': float,
    'patterns_found': int,         # Basé sur l'analyse des retours de fonction
    'regions_processed': int,
    'solver_status_counts': dict,  # Statuts extraits des valeurs de retour
    'avg_ilp_time': float,         # Temps moyen par appel de fonction
    'ilp_time_ratio': float
}
```

#### `reset_tracker()`
Remet à zéro tous les compteurs de performance.

### Context Manager

#### `PipelineTracker`
Gestionnaire de contexte pour le suivi complet d'un pipeline :
```python
with PipelineTracker() as tracker:
    # Votre code
    results = tracker.get_results()
```

## Exemple d'utilisation complète

```python
import numpy as np
from utils.decorators.ilp_tracker import (
    track_ilp_call, track_preprocessing, track_clustering, 
    PipelineTracker, get_performance_summary
)

@track_preprocessing()
def load_and_preprocess(data_path):
    # Simulation de chargement et prétraitement
    matrix = np.random.choice([0, 1], size=(100, 50))
    regions = [[0, 10], [10, 20], [20, 30]]
    return matrix, regions

@track_clustering()
def perform_clustering(matrix, regions):
    # Simulation de clustering avec appels ILP
    clusters = []
    for region in regions:
        result = solve_region(matrix, region)  # 3 appels seront comptés
        clusters.append(result)
    return clusters

@track_ilp_call(phase='region_optimization')
def solve_region(matrix, region):
    # Simulation de résolution ILP
    # Même si le solveur interne fait 100 optimisations, ceci compte comme 1 appel
    import time
    time.sleep(0.1)  # Simule le temps de calcul
    rows, cols = [1, 2, 3], [0, 1]
    success = True
    return rows, cols, success

# Utilisation
def main():
    with PipelineTracker() as pipeline:
        # Exécution du pipeline
        matrix, regions = load_and_preprocess("data.csv")
        clusters = perform_clustering(matrix, regions)  # 3 régions = 3 appels ILP
        
        # Analyse des résultats
        performance = pipeline.get_results()
        
        print("=== Rapport de Performance ===")
        print(f"Temps total pipeline: {performance['pipeline_total_time']:.2f}s")
        print(f"Appels aux fonctions ILP: {performance['ilp_calls_total']}")  # Affichera 3
        print(f"Temps moyen par appel de fonction: {performance['avg_ilp_time']:.3f}s")
        print(f"Régions traitées: {performance['regions_processed']}")
        print(f"Patterns trouvés: {performance['patterns_found']}")
        
        if performance['solver_status_counts']:
            print(f"Statuts des appels: {performance['solver_status_counts']}")

if __name__ == "__main__":
    main()
```

## Métriques collectées

- **Timing** : Temps d'exécution par phase et total
- **Function calls** : Nombre d'appels aux fonctions décorées (pas les optimisations internes)
- **Success rates** : Taux de succès/échec basé sur les valeurs de retour des fonctions
- **Matrix info** : Taille et forme des matrices traitées
- **Ratios** : Proportion du temps passé dans les fonctions ILP vs autres opérations

## Thread Safety

Le module utilise `threading.local()` pour assurer la sécurité des threads, permettant un suivi indépendant par thread dans les applications multi-threadées.

## Cas d'usage typiques

1. **Benchmarking** : Comparer les performances de différents algorithmes au niveau fonction
2. **Debugging** : Identifier les fonctions qui prennent le plus de temps
3. **Monitoring production** : Surveiller les performances en temps réel
4. **Optimisation** : Mesurer l'impact des améliorations de code
5. **Reporting** : Générer des rapports de performance pour les stakeholders

## Notes importantes

⚠️ **Limitation** : Ce tracker mesure les appels aux fonctions Python décorées, pas les optimisations internes du solveur ILP. Si le solveur fait plusieurs passes d'optimisation dans une seule fonction, cela sera compté comme un seul appel.
