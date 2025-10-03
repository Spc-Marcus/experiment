import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_results(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    # Normalize expected columns
    expected = {
        'Error-Rate': 'error_rate',
        'Haplotype': 'haplotype',
        'Time-ILP-V1': 'time_v1',
        'Time-ILP-V2': 'time_v2',
        'ILP-Steps-V1': 'steps_v1',
        'ILP-Steps-V2': 'steps_v2',
        'Matrix-Rows': 'rows',
        'Matrix-Cols': 'cols',
        'Clusters-Equivalent': 'clusters_equivalent',
    }
    missing = [c for c in expected.keys() if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")
    df = df.rename(columns=expected)
    # Types
    df['error_rate'] = df['error_rate'].astype(float)
    df['haplotype'] = df['haplotype'].astype(int)
    df['time_v1'] = df['time_v1'].astype(float)
    df['time_v2'] = df['time_v2'].astype(float)
    df['steps_v1'] = df['steps_v1'].astype(int)
    df['steps_v2'] = df['steps_v2'].astype(int)
    if 'speedup_v2_over_v1' not in df.columns:
        df['speedup_v2_over_v1'] = df['time_v1'] / df['time_v2']
    if 'speedup_v1_over_v2' not in df.columns:
        df['speedup_v1_over_v2'] = df['time_v2'] / df['time_v1']
    return df


def remove_outliers(df: pd.DataFrame, columns: list = None, percentile: float = 1.0) -> pd.DataFrame:
    """Remove top percentile outliers for specified columns"""
    if columns is None:
        columns = ['time_v1', 'time_v2', 'steps_v1', 'steps_v2']
    
    df_clean = df.copy()
    initial_count = len(df_clean)
    
    for col in columns:
        if col in df_clean.columns:
            # Supprimer seulement les valeurs dans le top percentile%
            threshold = df_clean[col].quantile(1 - percentile/100)
            df_clean = df_clean[df_clean[col] <= threshold]
    
    removed_count = initial_count - len(df_clean)
    if removed_count > 0:
        print(f"Removed {removed_count} outliers ({removed_count/initial_count*100:.1f}% of data, target: {percentile}%)")
    
    return df_clean


def aggregate(df: pd.DataFrame, agg: str = 'median') -> pd.DataFrame:
    agg_fn = {'mean': 'mean', 'median': 'median'}[agg]
    grouped = (df
               .groupby(['error_rate', 'haplotype'], as_index=False)
               .agg(time_v1=("time_v1", agg_fn), time_v2=("time_v2", agg_fn),
                    steps_v1=("steps_v1", agg_fn), steps_v2=("steps_v2", agg_fn),
                    speedup_v2_over_v1=("speedup_v2_over_v1", agg_fn)))
    counts = df.groupby(['error_rate', 'haplotype'], as_index=False).size().rename(columns={'size': 'n'})
    grouped = grouped.merge(counts, on=['error_rate', 'haplotype'], how='left')
    return grouped


def plot_times_by_error_rate(df_ag: pd.DataFrame, outdir: str, v1_name: str = "V1", v2_name: str = "V2"):
    # One subplot per haplotype
    haplos = sorted(df_ag['haplotype'].unique())
    n = len(haplos)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3.2*rows), squeeze=False, sharex=True)
    for i, h in enumerate(haplos):
        ax = axes[i//cols][i%cols]
        sub = df_ag[df_ag['haplotype'] == h].sort_values('error_rate')
        ax.plot(sub['error_rate'], sub['time_v1'], marker='o', label=v1_name)
        ax.plot(sub['error_rate'], sub['time_v2'], marker='s', label=v2_name)
        # X tick labels with sample size per point
        ax.set_xticks(sub['error_rate'].tolist())
        xticklabels = [f"{er:.3f} (n={int(nv)})" for er, nv in zip(sub['error_rate'], sub.get('n', pd.Series([1]*len(sub))))]
        ax.set_xticklabels(xticklabels, rotation=30, ha='right')
        # Title with total n for this haplotype
        ax.set_title(f'Haplotype={h} (n={int(sub["n"].sum()) if "n" in sub.columns else len(sub)})')
        ax.set_xlabel('Error rate')
        ax.set_ylabel('Time (s)')
        ax.grid(True, alpha=0.3)
        ax.legend()
    # Hide unused axes
    for j in range(n, rows*cols):
        fig.delaxes(axes[j//cols][j%cols])
    fig.tight_layout()
    out = os.path.join(outdir, 'times_by_error_rate_per_haplotype.png')
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_speedup_heatmap(df_ag: pd.DataFrame, outdir: str, v1_name: str = "V1", v2_name: str = "V2"):
    # Pivot: rows=error_rate, cols=haplotype
    pivot = df_ag.pivot(index='error_rate', columns='haplotype', values='speedup_v2_over_v1')
    counts = df_ag.pivot(index='error_rate', columns='haplotype', values='n') if 'n' in df_ag.columns else None
    fig, ax = plt.subplots(figsize=(max(6, 1.5*len(pivot.columns)), max(4, 0.8*len(pivot.index))))
    im = ax.imshow(pivot.values, aspect='auto', cmap='RdYlGn', vmin=0, vmax=max(2, float(pivot.max().max())))
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"{er:.3f}" for er in pivot.index])
    ax.set_xlabel('Haplotype')
    ax.set_ylabel('Error rate')
    ax.set_title(f'Speedup ({v2_name} faster when > 1)')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f'time_{v1_name} / time_{v2_name}')
    # Annotate values
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            if counts is not None:
                n = counts.values[i, j] if not np.isnan(counts.values[i, j]) else 0
                ax.text(j, i, f"{val:.2f} (n={int(n)})", ha='center', va='center', fontsize=8, color='black')
            else:
                ax.text(j, i, f"{val:.2f}", ha='center', va='center', fontsize=8, color='black')
    fig.tight_layout()
    out = os.path.join(outdir, 'speedup_heatmap.png')
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_v1_vs_v2_scatter(df: pd.DataFrame, outdir: str, v1_name: str = "V1", v2_name: str = "V2"):
    fig, ax = plt.subplots(figsize=(6, 5))
    # Size or color by haplotype
    haplos = sorted(df['haplotype'].unique())
    cmap = plt.get_cmap('tab10')
    for idx, h in enumerate(haplos):
        sub = df[df['haplotype'] == h]
        ax.scatter(sub['time_v1'], sub['time_v2'], alpha=0.6, label=f'h={h} (n={len(sub)})', color=cmap(idx % 10))
    limit = max(df['time_v1'].max(), df['time_v2'].max())
    ax.plot([0, limit], [0, limit], 'k--', linewidth=1, label='y=x')
    ax.set_xlabel(f'Time {v1_name} (s)')
    ax.set_ylabel(f'Time {v2_name} (s)')
    ax.set_title(f'{v1_name} vs {v2_name} runtime')
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    out = os.path.join(outdir, 'scatter_time_v1_vs_v2.png')
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Plot graphs from ILP results CSV (like results.csv).')
    ap.add_argument('--csv', required=True, help='Path to results CSV')
    ap.add_argument('--outdir', default=None, help='Directory to save plots (default: benchmarks/figures/YYYYmmdd_HHMMSS)')
    ap.add_argument('--agg', choices=['mean', 'median'], default='median', help='Aggregation for repeated points')
    ap.add_argument('--v1-name', default='Max_one_V1', help='Name for V1 data in plots (default: Max_one_V1)')
    ap.add_argument('--v2-name', default='Max_one_V2', help='Name for V2 data in plots (default: Max_one_V2)')
    ap.add_argument('--remove-outliers', action='store_true', help='Remove top 1%% outliers')
    ap.add_argument('--outlier-percentile', type=float, default=1.0, help='Percentage of top outliers to remove (default: 1.0)')
    args = ap.parse_args()

    df = load_results(args.csv)
    
    # Remove outliers if requested
    if args.remove_outliers:
        df = remove_outliers(df, percentile=args.outlier_percentile)

    # Prepare output dir
    if args.outdir is None:
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        outdir = os.path.join(os.path.dirname(args.csv) or '.', 'figures', ts)
    else:
        outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Aggregated plots by (error_rate, haplotype)
    df_ag = aggregate(df, agg=args.agg)

    plot_times_by_error_rate(df_ag, outdir, args.v1_name, args.v2_name)
    plot_speedup_heatmap(df_ag, outdir, args.v1_name, args.v2_name)
    plot_v1_vs_v2_scatter(df, outdir, args.v1_name, args.v2_name)

    print(f"Saved plots to {outdir}")


if __name__ == '__main__':
    main()
