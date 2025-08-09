from pathlib import Path
import pandas as pd
import numpy as np

def load_ilp_data(csv_path: str) -> pd.DataFrame:
    """Charge les données et calcule Success_Rate_Pre / Success_Rate_ILP."""
    dtype_dict = {
        'Threshold': 'float64',
        'Error-Rate': 'float64',
        'Strip': 'int64',
        'Haplotype': 'int64',
        'Matrix-Size-cols': 'int64',
        'Matrix-Size-rows': 'int64',
        'Time-Pre-Processing': 'float64',
        'Steps-Count-Pre': 'int64',
        'Unused-Cols-Pre': 'int64',
        'Final-Clusters-Pre': 'int64',
        'Orphan-Reads-Pre': 'int64',
        'Steps-Count-ILP': 'int64',
        'Final-Clusters-ILP': 'int64',
        'Orphan-Reads-ILP': 'int64',
        'Distance': 'float64'
    }
    df = pd.read_csv(csv_path, dtype=dtype_dict)
    # Colonnes de succès (cluster == haplotype)
    df['Success_Rate_Pre'] = (df['Final-Clusters-Pre'] == df['Haplotype']).astype(int)
    df['Success_Rate_ILP'] = (df['Final-Clusters-ILP'] == df['Haplotype']).astype(int)
    return df

def _select_best_by_key(df: pd.DataFrame, key_col: str, metric_col: str) -> pd.DataFrame:
    """
    Pour chaque valeur de key_col (Error-Rate ou Haplotype), choisir le couple (Threshold, Distance)
    qui maximise la moyenne de metric_col. Strip est ignoré (moyenné).
    Ex aequo: succès desc, support desc, Distance asc, Threshold desc.
    """
    results = []
    params = ['Threshold', 'Distance']  # Strip retiré (moyenné)

    for key_val, sub in df.groupby(key_col):
        # Agrégation sur (Threshold, Distance) uniquement
        agg = (
            sub.groupby(params)[metric_col]
            .agg(['mean', 'count'])
            .reset_index()
            .rename(columns={'mean': 'Mean_Success', 'count': 'Support'})
        )
        if agg.empty:
            continue

        # Triage pour tie-breaks (sans Strip)
        agg = agg.sort_values(
            by=['Mean_Success', 'Support', 'Distance', 'Threshold'],
            ascending=[False, False, True, False]
        ).reset_index(drop=True)

        best = agg.iloc[0].to_dict()
        best[key_col] = key_val
        results.append(best)

    if not results:
        return pd.DataFrame(columns=[key_col, *params, 'Mean_Success', 'Support'])

    cols = [key_col, *params, 'Mean_Success', 'Support']
    return pd.DataFrame(results)[cols]

def best_params_when_known_error_rate(df: pd.DataFrame, prefer: str = 'ILP') -> pd.DataFrame:
    """
    Maximise le succès quand seul Error-Rate est connu.
    - On cherche le meilleur couple (Threshold, Distance) pour chaque Error-Rate.
    - Strip (et le reste) sont moyennés.
    prefer: 'ILP' (par défaut) ou 'PRE' pour choisir la métrique.
    """
    metric_col = 'Success_Rate_ILP' if prefer.upper() == 'ILP' else 'Success_Rate_Pre'
    return _select_best_by_key(df, key_col='Error-Rate', metric_col=metric_col)

def export_best_tables(
    csv_path="/home/mafoin/stage/experiment/efficience_ALL/results.csv",
    out_dir="/home/mafoin/stage/experiment/Param_best",
    prefer='ILP'
):
    """Charge les données, calcule les meilleurs paramètres et exporte en CSV."""
    df = load_ilp_data(csv_path)

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    best_err = best_params_when_known_error_rate(df, prefer=prefer)

    err_file = out_path / f"best_by_error_rate_{prefer.lower()}.csv"
    best_err.to_csv(err_file, index=False)

    print(f"Export OK:")
    print(f"- {err_file}")

if __name__ == "__main__":
    # Exemple d'utilisation par défaut (maximisation sur ILP)
    export_best_tables(prefer='ILP')
