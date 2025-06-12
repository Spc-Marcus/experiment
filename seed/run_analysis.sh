#!/bin/bash

echo "🔍 Lancement de l'analyse des résultats de seeding..."

# Trouver le fichier de résultats le plus récent
latest_result=$(ls -t seeding_experiment_results_*.csv 2>/dev/null | head -n1)

if [ -z "$latest_result" ]; then
    echo "❌ Aucun fichier de résultats trouvé"
    echo "Veuillez d'abord exécuter les expériences de seeding"
    exit 1
fi

echo "📁 Fichier de résultats trouvé: $latest_result"

# Vérifier si matplotlib est disponible
python -c "import matplotlib.pyplot as plt" 2>/dev/null
if [ $? -eq 0 ]; then
    echo "📊 Matplotlib disponible - Génération de l'analyse complète..."
    python analyze_seeding_results.py "$latest_result" --output-dir "seeding_analysis_$(date +%Y%m%d_%H%M%S)"
else
    echo "⚠️  Matplotlib non disponible - Analyse rapide seulement..."
    python quick_analysis.py "$latest_result"
fi

echo "✅ Analyse terminée!"
