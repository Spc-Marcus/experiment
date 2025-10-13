import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter


def load_results(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    expected = {
        'Error-Rate': 'error_rate',
        'Haplotype': 'haplotype',
        'Final-Clusters-V1': 'final_clusters_v1',
        'Final-Clusters-V2': 'final_clusters_v2',
        'Time-ILP-V1': 'time_v1',
        'Time-ILP-V2': 'time_v2',
    }
    missing = [c for c in expected.keys() if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in CSV: {missing}")
    df = df.rename(columns=expected)
    df['error_rate'] = df['error_rate'].astype(float)
    df['haplotype'] = df['haplotype'].astype(int)
    # Ensure cluster columns are numeric
    df['final_clusters_v1'] = pd.to_numeric(df['final_clusters_v1'], errors='coerce')
    df['final_clusters_v2'] = pd.to_numeric(df['final_clusters_v2'], errors='coerce')
    return df


def remove_outliers(df: pd.DataFrame, columns: list = None, percentile: float = 1.0) -> pd.DataFrame:
    if columns is None:
        columns = ['final_clusters_v1', 'final_clusters_v2']
    df_clean = df.copy()
    initial = len(df_clean)
    for col in columns:
        if col in df_clean.columns:
            thresh = df_clean[col].quantile(1 - percentile / 100)
            df_clean = df_clean[df_clean[col] <= thresh]
    removed = initial - len(df_clean)
    if removed:
        print(f"Removed {removed} outliers ({removed/initial*100:.2f}% of rows)")
    return df_clean


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    # Differences: clusters found - haplotype (reference)
    df = df.copy()
    df['cluster_diff_v1'] = df['final_clusters_v1'] - df['haplotype']
    df['cluster_diff_v2'] = df['final_clusters_v2'] - df['haplotype']
    # Difference between V1 and V2 in absolute clusters found
    df['clusters_v1_minus_v2'] = df['final_clusters_v1'] - df['final_clusters_v2']
    return df


def plot_diff_by_error(df: pd.DataFrame, version_col: str, version_label: str, outpath: str):
    # Group by error rate and create violin/box plots of the differences
    groups = df.groupby('error_rate')
    error_rates = sorted(df['error_rate'].unique())
    data = [groups.get_group(er)[version_col].dropna().values if er in groups.groups else np.array([]) for er in error_rates]

    if not any(len(d) for d in data):
        print(f"No data to plot for {version_label}")
        return

    plt.figure(figsize=(12, 6))
    # positions like in plot_one.py (1-based)
    positions = [i + 1 for i in range(len(error_rates))]
    counts = [len(d) for d in data]
    labels = [f"{er:.3f}" for er in error_rates]

    # choose color similar to plot_one
    color_map = {
        'cluster_diff_v1': 'lightblue',
        'cluster_diff_v2': 'lightgreen',
        'clusters_v1_minus_v2': 'lightgray'
    }
    color = color_map.get(version_col, '#88c999')

    # Always draw violin plots (no boxplot option)
    parts = plt.violinplot(data, positions=positions, showmeans=False, showmedians=True)
    for pc in parts['bodies']:
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
    # style median lines if present
    if 'cmedians' in parts:
        parts['cmedians'].set_edgecolor('firebrick')
        parts['cmedians'].set_linewidth(1.2)

    plt.xticks(positions, labels, rotation=45)
    plt.xlabel('Error rate')
    total_n = sum(counts)
    plt.ylabel('Clusters found - Haplotype')
    plt.title(f'Cluster difference by error rate ({version_label}) — total n={total_n}')
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color='red', linestyle='--', linewidth=1, label='Difference = 0')
    # remove default legend to avoid showing unwanted markers
    # plt.legend()

    # Annotate counts per difference value under each x-group
    ax = plt.gca()
    y_min, y_max = ax.get_ylim()
    rng = y_max - y_min if (y_max - y_min) != 0 else 1.0
    # extend lower ylim to leave space for annotations
    new_ymin = y_min - 0.18 * rng
    ax.set_ylim(new_ymin, y_max)
    ann_y = y_min - 0.06 * rng
    for i, d in enumerate(data):
        if len(d) == 0:
            continue
        ctr = Counter(d)
        ann = ', '.join([f'{int(k)}:{v}' for k, v in sorted(ctr.items())])
        # wrap long annotations if needed
        if len(ann) > 60:
            # split into chunks
            parts = []
            cur = ''
            for token in ann.split(', '):
                if len(cur) + len(token) + 2 > 60:
                    parts.append(cur)
                    cur = token
                else:
                    cur = (cur + ', ' + token) if cur else token
            if cur:
                parts.append(cur)
            for j, part in enumerate(parts):
                plt.text(positions[i], ann_y - j * 0.04 * rng, part, ha='center', va='top', fontsize=7)
        else:
            plt.text(positions[i], ann_y, ann, ha='center', va='top', fontsize=8)

    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()
    print(f"Saved plot: {outpath}")


def write_summary(df: pd.DataFrame, outdir: str, v1_name: str, v2_name: str):
    summary_path = os.path.join(outdir, 'summary_cluster_diffs.txt')
    cols = [
        ('cluster_diff_v1', v1_name),
        ('cluster_diff_v2', v2_name),
        ('clusters_v1_minus_v2', f'{v1_name} - {v2_name}')
    ]
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write('Summary of cluster differences grouped by error rate\n')
        f.write('Generated by plot_cluster_differences.py\n\n')
        for col, label in cols:
            f.write(f'=== {label} ({col}) ===\n')
            grp = df.groupby('error_rate')[col]
            f.write('error_rate\tn\tmean\tmedian\tstd\tmin\t25%\t75%\tmax\n')
            for er, series in grp:
                ser = series.dropna()
                if ser.empty:
                    continue
                n = ser.size
                mean = ser.mean()
                median = ser.median()
                std = ser.std()
                mn = ser.min()
                q25 = ser.quantile(0.25)
                q75 = ser.quantile(0.75)
                mx = ser.max()
                f.write(f'{er:.6f}\t{n}\t{mean:.4f}\t{median:.4f}\t{std:.4f}\t{mn}\t{q25}\t{q75}\t{mx}\n')
            f.write('\n')
        # global stats
        f.write('=== Global statistics ===\n')
        for col, label in cols:
            ser = df[col].dropna()
            if ser.empty:
                f.write(f'{label}: no data\n')
                continue
            f.write(f'{label}: n={ser.size}, mean={ser.mean():.4f}, median={ser.median():.4f}, std={ser.std():.4f}, min={ser.min()}, max={ser.max()}\n')
        f.write('\n')
        # True-positive summary (Final-Clusters == Haplotype) for V1 and V2
        f.write('=== True-positive summary (Final-Clusters == Haplotype) ===\n')
        for name, col in [(v1_name, 'final_clusters_v1'), (v2_name, 'final_clusters_v2')]:
            total_n = len(df)
            total_tp = int((df[col] == df['haplotype']).sum())
            total_pct = total_tp / total_n * 100 if total_n > 0 else 0.0
            f.write(f'{name}: exact_count={total_tp}\ttotal={total_n}\texact_pct={total_pct:.2f}%\n')
            f.write('Per error_rate:\n')
            f.write('error_rate\tn\ttrue_positive_count\ttrue_positive_pct\n')
            for er, g in df.groupby('error_rate'):
                n = len(g)
                tp = int((g[col] == g['haplotype']).sum())
                pct = tp / n * 100 if n > 0 else 0.0
                f.write(f'{er:.6f}\t{n}\t{tp}\t{pct:.2f}\n')
            f.write('\n')
    print(f'Wrote summary: {summary_path}')


def main():
    ap = argparse.ArgumentParser(description='Plot cluster differences (FinalClusters - Haplotype) per error rate for V1 and V2')
    ap.add_argument('--csv', required=True, help='Path to results CSV')
    ap.add_argument('--outdir', default=None, help='Output directory for plots')
    ap.add_argument('--v1-name', default='V1', help='Label for version 1')
    ap.add_argument('--v2-name', default='V2', help='Label for version 2')
    ap.add_argument('--remove-outliers', action='store_true', help='Remove top outliers by percentile')
    ap.add_argument('--outlier-percentile', type=float, default=1.0, help='Percentile to remove from the top (default 1.0)')
    args = ap.parse_args()

    df = load_results(args.csv)
    df = prepare_data(df)

    # If requested, remove outliers based on the difference columns (not raw cluster counts)
    if args.remove_outliers:
        diff_cols = ['cluster_diff_v1', 'cluster_diff_v2', 'clusters_v1_minus_v2']
        df = remove_outliers(df, columns=diff_cols, percentile=args.outlier_percentile)

    if args.outdir is None:
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        outdir = os.path.join(os.path.dirname(args.csv) or '.', 'figures', ts)
    else:
        outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Plot for V1
    out_v1 = os.path.join(outdir, f'cluster_diff_{args.v1_name}.png')
    plot_diff_by_error(df, 'cluster_diff_v1', args.v1_name, out_v1)

    # Plot for V2
    out_v2 = os.path.join(outdir, f'cluster_diff_{args.v2_name}.png')
    plot_diff_by_error(df, 'cluster_diff_v2', args.v2_name, out_v2)

    # Plot difference V1 - V2 in absolute clusters
    out_v1v2 = os.path.join(outdir, f'clusters_{args.v1_name}_minus_{args.v2_name}.png')
    plot_diff_by_error(df, 'clusters_v1_minus_v2', f"{args.v1_name} - {args.v2_name}", out_v1v2)

    # Write textual summary
    write_summary(df, outdir, args.v1_name, args.v2_name)

    print(f"Plots and summary written to {outdir}")


if __name__ == '__main__':
    main()
