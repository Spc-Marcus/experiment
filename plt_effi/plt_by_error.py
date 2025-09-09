import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
csv_files = glob.glob(os.path.join(DATA_DIR, 'best_by_error_rate_*.csv'))

if not csv_files:
    print(f"Aucun fichier trouvé dans {DATA_DIR}")
    raise SystemExit(0)

plt.figure(figsize=(10, 6))
courbes_tracees = 0

def charger_csv_robuste(path: str) -> pd.DataFrame:
    # Détection simple du séparateur à partir de la première ligne
    with open(path, 'r', encoding='utf-8') as f:
        first = f.readline().rstrip('\n')
    if '\t' in first and ',' not in first:
        sep = '\t'
    elif ',' in first and '\t' not in first:
        sep = ','
    else:
        # Ambigu : on tente tab puis virgule
        sep = None

    if sep:
        df = pd.read_csv(path, sep=sep)
    else:
        # Essais successifs
        for s in ('\t', ',', ';'):
            df = pd.read_csv(path, sep=s)
            if len(df.columns) > 1:
                break

    # Si encore une seule colonne contenant des virgules => reparse en ','
    if len(df.columns) == 1 and ',' in df.columns[0]:
        df = pd.read_csv(path, sep=',')

    # Nettoyage noms de colonnes
    df.columns = [c.strip() for c in df.columns]

    # Harmonisation : créer un mapping standard
    mapping = {}
    for c in df.columns:
        canon = c.lower().replace('-', '').replace('_', '')
        if canon == 'errorrate':
            mapping[c] = 'ErrorRate'
        elif canon == 'threshold':
            mapping[c] = 'Threshold'
        elif canon == 'distance':
            mapping[c] = 'Distance'
        elif canon == 'meansuccess':
            mapping[c] = 'MeanSuccess'
        elif canon == 'support':
            mapping[c] = 'Support'
    df = df.rename(columns=mapping)

    # Filtrer uniquement si colonnes essentielles présentes
    needed = {'ErrorRate', 'MeanSuccess'}
    if not needed.issubset(df.columns):
        print(f"[AVERTISSEMENT] Colonnes manquantes dans {path}: {df.columns.tolist()}")
        return pd.DataFrame()

    # Conversion numérique sûre
    for col in ['ErrorRate', 'MeanSuccess']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.dropna(subset=['ErrorRate', 'MeanSuccess'])
    df = df.sort_values('ErrorRate')
    return df

output_dir = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(output_dir, exist_ok=True)

for csv_path in csv_files:
    base = os.path.basename(csv_path)
    label = base.replace('best_by_error_rate_', '').replace('.csv', '')
    df = charger_csv_robuste(csv_path)
    if df.empty:
        continue

    print(f"{csv_path} -> colonnes: {df.columns.tolist()}  lignes: {len(df)}")

    # Tracé sur la figure combinée
    plt.plot(df['ErrorRate'], df['MeanSuccess'], marker='o', label=label)
    courbes_tracees += 1

    # Figure individuelle
    fig_indiv, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df['ErrorRate'], df['MeanSuccess'], marker='o')
    ax.set_xlabel('Error Rate')
    ax.set_ylabel('Mean Success')
    ax.set_title(f'Succès vs Error Rate - {label}')
    ax.grid(True)
    fig_indiv.tight_layout()
    indiv_path_png = os.path.join(output_dir, f"{label}.png")
    fig_indiv.savefig(indiv_path_png, dpi=150)
    plt.close(fig_indiv)

if courbes_tracees == 0:
    print("Aucune courbe tracée (colonnes introuvables).")
else:
    plt.xlabel('Error Rate')
    plt.ylabel('Mean Success')
    plt.title('Succès moyen en fonction du taux d\'erreur')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    combined_path = os.path.join(output_dir, "combined.png")
    plt.savefig(combined_path, dpi=150)
    print(f"Figure combinée sauvegardée: {combined_path}")
    plt.show()