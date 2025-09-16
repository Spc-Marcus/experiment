import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

print(f"Recherche dans le répertoire: {DATA_DIR}")
print("Fichiers CSV disponibles:")
all_csv = glob.glob(os.path.join(DATA_DIR, '*.csv'))
for f in all_csv:
    print(f"  {os.path.basename(f)}")

# Charger les fichiers de meilleurs paramètres
best_params_files = glob.glob(os.path.join(DATA_DIR, 'best_by_error_rate_*.csv'))
# Charger le fichier de résultats
results_file = os.path.join(DATA_DIR, 'results_max_one.csv')

if not best_params_files:
    print(f"Aucun fichier best_by_error_rate trouvé dans {DATA_DIR}")
    raise SystemExit(0)

if not os.path.exists(results_file):
    print(f"Fichier {results_file} non trouvé")
    raise SystemExit(0)

def charger_csv_robuste(path: str) -> pd.DataFrame:
    # Détection du séparateur
    with open(path, 'r', encoding='utf-8') as f:
        first = f.readline().rstrip('\n')
    
    if '\t' in first and ',' not in first:
        sep = '\t'
    elif ',' in first and '\t' not in first:
        sep = ','
    else:
        sep = None

    if sep:
        df = pd.read_csv(path, sep=sep)
    else:
        for s in ('\t', ',', ';'):
            df = pd.read_csv(path, sep=s)
            if len(df.columns) > 1:
                break

    if len(df.columns) == 1 and ',' in df.columns[0]:
        df = pd.read_csv(path, sep=',')

    # Nettoyage des noms de colonnes
    df.columns = [c.strip() for c in df.columns]
    return df

# Créer les répertoires de sortie
output_dir = os.path.join(os.path.dirname(__file__), 'figures')
by_haplo_dir = os.path.join(output_dir, 'by_haplo')
by_stripe_dir = os.path.join(output_dir, 'by_stripe')
os.makedirs(by_haplo_dir, exist_ok=True)
os.makedirs(by_stripe_dir, exist_ok=True)

# Charger les meilleurs paramètres
best_params = {}
for best_file in best_params_files:
    base_name = os.path.basename(best_file).replace('best_by_error_rate_', '').replace('.csv', '')
    df_best = charger_csv_robuste(best_file)
    
    # Normaliser les noms de colonnes
    column_mapping = {}
    for col in df_best.columns:
        canon = col.lower().replace('-', '').replace('_', '')
        if 'error' in canon and 'rate' in canon:
            column_mapping[col] = 'ErrorRate'
        elif 'threshold' in canon:
            column_mapping[col] = 'Threshold'
        elif 'distance' in canon:
            column_mapping[col] = 'Distance'
    
    df_best = df_best.rename(columns=column_mapping)
    
    if 'ErrorRate' in df_best.columns and 'Threshold' in df_best.columns and 'Distance' in df_best.columns:
        best_params[base_name] = df_best
        print(f"Chargé {len(df_best)} meilleurs paramètres pour {base_name}")

# Charger les résultats détaillés
df_results = charger_csv_robuste(results_file)
print(f"Colonnes dans results_max_one.csv: {df_results.columns.tolist()}")

# Normaliser les noms de colonnes pour les résultats
column_mapping = {}
for col in df_results.columns:
    canon = col.lower().replace('-', '').replace('_', '')
    if 'error' in canon and 'rate' in canon:
        column_mapping[col] = 'ErrorRate'
    elif 'threshold' in canon:
        column_mapping[col] = 'Threshold'
    elif 'distance' in canon:
        column_mapping[col] = 'Distance'
    elif 'haplotype' in canon:
        column_mapping[col] = 'Haplotype'
    elif 'strip' in canon and canon == 'strip':  # Nombre de stripes attendues
        column_mapping[col] = 'ExpectedStripes'
    elif 'final' in canon and 'cluster' in canon and 'ilp' in canon:
        column_mapping[col] = 'FinalClustersILP'
    elif 'steps' in canon and 'count' in canon and 'pre' in canon:
        column_mapping[col] = 'StepsCountPre'
    elif 'steps' in canon and 'count' in canon and 'ilp' in canon:
        column_mapping[col] = 'StepsCountILP'

df_results = df_results.rename(columns=column_mapping)

required_cols = ['ErrorRate', 'Threshold', 'Distance', 'Haplotype', 'ExpectedStripes', 'FinalClustersILP', 'StepsCountPre', 'StepsCountILP']
missing_cols = [col for col in required_cols if col not in df_results.columns]
if missing_cols:
    print(f"Colonnes manquantes dans results_max_one.csv: {missing_cols}")
    print(f"Colonnes disponibles après mapping: {df_results.columns.tolist()}")
    raise SystemExit(0)

# Convertir en numérique
for col in required_cols:
    df_results[col] = pd.to_numeric(df_results[col], errors='coerce')

# Calculer le nombre total de stripes trouvées
df_results['FoundStripes'] = df_results['StepsCountPre'] + df_results['StepsCountILP']

# Utiliser les meilleurs paramètres (prendre le premier dataset disponible)
if not best_params:
    print("Aucun paramètre optimal trouvé")
    raise SystemExit(0)

dataset_name = list(best_params.keys())[0]
df_best = best_params[dataset_name]
print(f"Utilisation des paramètres de: {dataset_name}")

# Filtrer selon les meilleurs paramètres
filtered_results = []

for _, best_row in df_best.iterrows():
    error_rate = best_row['ErrorRate']
    threshold = best_row['Threshold']
    distance = best_row['Distance']
    
    # Filtrer les résultats pour ces paramètres exacts
    mask = (
        (df_results['ErrorRate'] == error_rate) &
        (df_results['Threshold'] == threshold) &
        (df_results['Distance'] == distance)
    )
    
    filtered = df_results[mask].copy()
    if not filtered.empty:
        # Calculer les différences
        filtered['ClusterDifference'] = filtered['FinalClustersILP'] - filtered['Haplotype']
        filtered['StripeDifference'] = filtered['FoundStripes'] - filtered['ExpectedStripes']
        filtered_results.append(filtered)
        print(f"Trouvé {len(filtered)} lignes pour error_rate={error_rate}, threshold={threshold}, distance={distance}")

if not filtered_results:
    print("Aucun résultat filtré trouvé.")
    raise SystemExit(0)

# Combiner tous les résultats filtrés
combined_df = pd.concat(filtered_results, ignore_index=True)

# === GRAPHIQUES PAR NOMBRE D'HAPLOTYPES ===
unique_haplotypes = sorted(combined_df['Haplotype'].unique())
print(f"Nombres d'haplotypes trouvés: {unique_haplotypes}")

for haplotype_count in unique_haplotypes:
    haplotype_data = combined_df[combined_df['Haplotype'] == haplotype_count]
    
    if haplotype_data.empty:
        continue
    
    # Compter le nombre de matrices pour chaque taux d'erreur
    error_rates = sorted(haplotype_data['ErrorRate'].unique())
    matrix_counts_by_error = {}
    
    for er in error_rates:
        count = len(haplotype_data[haplotype_data['ErrorRate'] == er])
        matrix_counts_by_error[er] = count
    
    # Vérifier si le nombre de matrices est constant
    unique_counts = set(matrix_counts_by_error.values())
    if len(unique_counts) == 1:
        matrix_count_text = f" ({list(unique_counts)[0]} matrices)"
    else:
        matrix_count_text = f" (matrices: {dict(matrix_counts_by_error)})"
    
    plt.figure(figsize=(14, 8))
    
    boxplot_data = []
    labels = []
    
    for er in error_rates:
        data_for_er = haplotype_data[haplotype_data['ErrorRate'] == er]['ClusterDifference']
        if len(data_for_er) > 0:
            boxplot_data.append(data_for_er)
            labels.append(f'{er:.3f}')
    
    if boxplot_data:
        plt.boxplot(boxplot_data, labels=labels)
        plt.xlabel('Taux d\'erreur')
        plt.ylabel('Différence (Clusters trouvés - Clusters attendus)')
        plt.title(f'Distribution de la différence de clusters par taux d\'erreur\n({int(haplotype_count)} haplotypes{matrix_count_text})')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)
        plt.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Différence nulle')
        plt.legend()
        plt.tight_layout()
        
        # Sauvegarder
        boxplot_path = os.path.join(by_haplo_dir, f"cluster_difference_{int(haplotype_count)}_haplotypes.png")
        plt.savefig(boxplot_path, dpi=150)
        print(f"Boîtes à moustaches pour {int(haplotype_count)} haplotypes sauvegardées: {boxplot_path}")
        plt.close()

# === GRAPHIQUES PAR NOMBRE DE STRIPES ATTENDUES ===
unique_expected_stripes = sorted(combined_df['ExpectedStripes'].unique())
print(f"Nombres de stripes attendues trouvés: {unique_expected_stripes}")

for stripe_count in unique_expected_stripes:
    stripe_data = combined_df[combined_df['ExpectedStripes'] == stripe_count]
    
    if stripe_data.empty:
        continue
    
    # Compter le nombre de matrices pour chaque taux d'erreur
    error_rates = sorted(stripe_data['ErrorRate'].unique())
    matrix_counts_by_error = {}
    
    for er in error_rates:
        count = len(stripe_data[stripe_data['ErrorRate'] == er])
        matrix_counts_by_error[er] = count
    
    # Vérifier si le nombre de matrices est constant
    unique_counts = set(matrix_counts_by_error.values())
    if len(unique_counts) == 1:
        matrix_count_text = f" ({list(unique_counts)[0]} matrices)"
    else:
        matrix_count_text = f" (matrices: {dict(matrix_counts_by_error)})"
    
    plt.figure(figsize=(14, 8))
    
    boxplot_data = []
    labels = []
    
    for er in error_rates:
        data_for_er = stripe_data[stripe_data['ErrorRate'] == er]['StripeDifference']
        if len(data_for_er) > 0:
            boxplot_data.append(data_for_er)
            labels.append(f'{er:.3f}')
    
    if boxplot_data:
        plt.boxplot(boxplot_data, labels=labels)
        plt.xlabel('Taux d\'erreur')
        plt.ylabel('Différence (Stripes trouvées - Stripes attendues)')
        plt.title(f'Distribution de la différence de stripes par taux d\'erreur\n({int(stripe_count)} stripes attendues{matrix_count_text})')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)
        plt.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Différence nulle')
        plt.legend()
        plt.tight_layout()
        
        # Sauvegarder
        boxplot_path = os.path.join(by_stripe_dir, f"stripe_difference_{int(stripe_count)}_expected_stripes.png")
        plt.savefig(boxplot_path, dpi=150)
        print(f"Boîtes à moustaches pour {int(stripe_count)} stripes attendues sauvegardées: {boxplot_path}")
        plt.close()

# === GRAPHIQUE COMBINÉ PAR HAPLOTYPES ===
if len(unique_haplotypes) > 0:
    fig, axes = plt.subplots(len(unique_haplotypes), 1, figsize=(16, 6*len(unique_haplotypes)))
    
    if len(unique_haplotypes) == 1:
        axes = [axes]

    for i, haplotype_count in enumerate(unique_haplotypes):
        haplotype_data = combined_df[combined_df['Haplotype'] == haplotype_count]
        
        # Compter le nombre de matrices
        error_rates = sorted(haplotype_data['ErrorRate'].unique())
        matrix_counts_by_error = {}
        
        for er in error_rates:
            count = len(haplotype_data[haplotype_data['ErrorRate'] == er])
            matrix_counts_by_error[er] = count
        
        unique_counts = set(matrix_counts_by_error.values())
        if len(unique_counts) == 1:
            matrix_count_text = f" ({list(unique_counts)[0]} matrices)"
        else:
            matrix_count_text = f""
        
        boxplot_data = []
        labels = []
        
        for er in error_rates:
            data_for_er = haplotype_data[haplotype_data['ErrorRate'] == er]['ClusterDifference']
            if len(data_for_er) > 0:
                boxplot_data.append(data_for_er)
                labels.append(f'{er:.3f}')
        
        if boxplot_data:
            axes[i].boxplot(boxplot_data, labels=labels)
            axes[i].set_xlabel('Taux d\'erreur')
            axes[i].set_ylabel('Différence clusters')
            axes[i].set_title(f'{int(haplotype_count)} haplotypes{matrix_count_text}')
            axes[i].tick_params(axis='x', rotation=45)
            axes[i].grid(True, alpha=0.3)
            axes[i].axhline(y=0, color='red', linestyle='--', alpha=0.7)

    plt.tight_layout()
    combined_haplo_path = os.path.join(by_haplo_dir, "cluster_difference_all_haplotypes.png")
    plt.savefig(combined_haplo_path, dpi=150)
    print(f"Graphique combiné haplotypes sauvegardé: {combined_haplo_path}")
    plt.close()

# === GRAPHIQUE COMBINÉ PAR STRIPES ===
if len(unique_expected_stripes) > 0:
    fig, axes = plt.subplots(len(unique_expected_stripes), 1, figsize=(16, 6*len(unique_expected_stripes)))
    
    if len(unique_expected_stripes) == 1:
        axes = [axes]

    for i, stripe_count in enumerate(unique_expected_stripes):
        stripe_data = combined_df[combined_df['ExpectedStripes'] == stripe_count]
        
        # Compter le nombre de matrices
        error_rates = sorted(stripe_data['ErrorRate'].unique())
        matrix_counts_by_error = {}
        
        for er in error_rates:
            count = len(stripe_data[stripe_data['ErrorRate'] == er])
            matrix_counts_by_error[er] = count
        
        unique_counts = set(matrix_counts_by_error.values())
        if len(unique_counts) == 1:
            matrix_count_text = f" ({list(unique_counts)[0]} matrices)"
        else:
            matrix_count_text = f""
        
        boxplot_data = []
        labels = []
        
        for er in error_rates:
            data_for_er = stripe_data[stripe_data['ErrorRate'] == er]['StripeDifference']
            if len(data_for_er) > 0:
                boxplot_data.append(data_for_er)
                labels.append(f'{er:.3f}')
        
        if boxplot_data:
            axes[i].boxplot(boxplot_data, labels=labels)
            axes[i].set_xlabel('Taux d\'erreur')
            axes[i].set_ylabel('Différence stripes')
            axes[i].set_title(f'{int(stripe_count)} stripes attendues{matrix_count_text}')
            axes[i].tick_params(axis='x', rotation=45)
            axes[i].grid(True, alpha=0.3)
            axes[i].axhline(y=0, color='red', linestyle='--', alpha=0.7)

    plt.tight_layout()
    combined_stripe_path = os.path.join(by_stripe_dir, "stripe_difference_all_expected_stripes.png")
    plt.savefig(combined_stripe_path, dpi=150)
    print(f"Graphique combiné stripes sauvegardé: {combined_stripe_path}")
    plt.close()

print(f"\nTous les graphiques ont été sauvegardés dans:")
print(f"- Haplotypes: {by_haplo_dir}")
print(f"- Stripes: {by_stripe_dir}")