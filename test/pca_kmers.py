#!/usr/bin/env python3
"""
pca_kmers.py

Fast, modern interactive HTML PCA plot for fastkmers sample outputs.
- Normalizes frequencies (sample x kmer)
- Computes 2D PCA with variance metrics
- Outputs a high-performance, standalone HTML application (ECharts / Canvas)
- Features instant legend-hover highlighting of sample points, tooltips, and controls.
"""

import argparse
import json
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA


def parse_args():
    parser = argparse.ArgumentParser(
        description="Normalize fastkmers TSV outputs, build sample x kmer matrix, and generate a fast interactive HTML PCA plot."
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="Path to one or more fastkmers output TSV files (header expected: kmer\\tcount).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="pca_kmers_plot.html",
        help="Path for output interactive HTML file (default: pca_kmers_plot.html).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random state for reproducibility (default: 42).",
    )
    parser.add_argument(
        "--save-matrix",
        type=str,
        default=None,
        help="Optional path to save normalized sample x kmer matrix CSV/TSV.",
    )
    return parser.parse_args()


def load_and_normalize_samples(file_paths):
    sample_dicts = {}
    
    for path_str in file_paths:
        filepath = Path(path_str)
        if not filepath.is_file():
            print(f"Warning: File {filepath} not found. Skipping...", file=sys.stderr)
            continue
            
        sample_name = filepath.name
        if sample_name.endswith(".tsv") or sample_name.endswith(".txt"):
            sample_name = filepath.stem
            
        print(f"Loading {filepath} as sample '{sample_name}'...")
        
        df = pd.read_csv(filepath, sep="\t")
        if "kmer" not in df.columns or "count" not in df.columns:
            print(
                f"Error: {filepath} missing expected 'kmer' or 'count' columns.",
                file=sys.stderr,
            )
            continue
            
        total_counts = df["count"].sum()
        if total_counts == 0:
            print(f"Warning: Total counts for {sample_name} is 0. Skipping...", file=sys.stderr)
            continue
            
        df["freq"] = df["count"] / total_counts
        sample_dicts[sample_name] = pd.Series(df["freq"].values, index=df["kmer"])

    if not sample_dicts:
        raise ValueError("No valid fastkmers TSV data loaded.")

    print(f"Building sample x kmer matrix for {len(sample_dicts)} samples...")
    matrix_df = pd.DataFrame(sample_dicts).T.fillna(0.0)
    print(f"Matrix shape: {matrix_df.shape[0]} samples x {matrix_df.shape[1]} unique k-mers")
    
    return matrix_df


def generate_echarts_html(matrix_df, output_path, random_state=42):
    n_samples = matrix_df.shape[0]
    if n_samples < 2:
        raise ValueError("PCA requires at least 2 samples.")

    n_components = min(2, n_samples)
    print(f"Running PCA reduction ({n_components} components)...")
    
    pca = PCA(n_components=n_components, random_state=random_state)
    pca_results = pca.fit_transform(matrix_df.values)
    
    exp_var = pca.explained_variance_ratio_ * 100
    pc1_var = f"{exp_var[0]:.2f}%"
    pc2_var = f"{exp_var[1]:.2f}%" if n_components > 1 else "0.00%"
    
    unique_kmers = (matrix_df > 0).sum(axis=1).tolist()
    sample_names = matrix_df.index.tolist()

    pc1_min = float(pca_results[:, 0].min())
    pc1_max = float(pca_results[:, 0].max())
    
    # Avoid min == max division error if values are identical
    if pc1_min == pc1_max:
        pc1_min -= 1.0
        pc1_max += 1.0

    points_data = []
    for i, name in enumerate(sample_names):
        x_val = float(pca_results[i, 0])
        y_val = float(pca_results[i, 1]) if n_components > 1 else 0.0
        kmers_count = int(unique_kmers[i])
        points_data.append([x_val, y_val, kmers_count, name])

    data_json = json.dumps(points_data, indent=2)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>fastkmers - PCA Interactive Plot</title>
    <script src="https://cdn.jsdelivr.net/npm/echarts@5.4.3/dist/echarts.min.js"></script>
    <style>
        :root {{
            --bg-color: #FFFFFF;
            --panel-bg: #F8FAFC;
            --text-color: #0F172A;
            --text-muted: #64748B;
        }}
        body {{
            margin: 0;
            padding: 0;
            background-color: var(--bg-color);
            color: var(--text-color);
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            display: flex;
            flex-direction: column;
            height: 100vh;
            overflow: hidden;
        }}
        header {{
            padding: 16px 24px;
            background-color: var(--panel-bg);
            border-bottom: 1px solid #E2E8F0;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        h1 {{
            margin: 0;
            font-size: 20px;
            font-weight: 600;
            color: #0284C7;
        }}
        .meta {{
            font-size: 13px;
            color: var(--text-muted);
        }}
        #chart-container {{
            flex: 1;
            width: 100%;
            height: 100%;
        }}
    </style>
</head>
<body>
    <header>
        <div>
            <h1>fastkmers — Interactive PCA Visualization</h1>
            <div class="meta">Sample x K-mer normalized frequency matrix ({n_samples} samples, {matrix_df.shape[1]} unique k-mers)</div>
        </div>
        <div class="meta">Color gradient continuous by PC1 value</div>
    </header>
    <div id="chart-container"></div>

    <script>
        const chartDom = document.getElementById('chart-container');
        const myChart = echarts.init(chartDom);

        const option = {{
            backgroundColor: '#FFFFFF',
            grid: {{
                top: '10%',
                left: '8%',
                right: '18%',
                bottom: '12%',
                containLabel: true
            }},
            tooltip: {{
                trigger: 'item',
                backgroundColor: '#FFFFFF',
                borderColor: '#CBD5E1',
                borderWidth: 1,
                textStyle: {{ color: '#0F172A' }},
                extraCssText: 'box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);',
                formatter: function (params) {{
                    const sampleName = params.value[3];
                    const x = params.value[0].toFixed(4);
                    const y = params.value[1].toFixed(4);
                    const kmers = params.value[2].toLocaleString();
                    return `
                        <div style="font-weight:bold; color:#0284C7; margin-bottom:4px;">${{sampleName}}</div>
                        <div><b>PC1:</b> ${{x}}</div>
                        <div><b>PC2:</b> ${{y}}</div>
                        <div><b>Unique K-mers:</b> ${{kmers}}</div>
                    `;
                }}
            }},
            visualMap: {{
                min: {pc1_min},
                max: {pc1_max},
                dimension: 0,
                orient: 'vertical',
                right: 20,
                top: 'center',
                text: ['High PC1', 'Low PC1'],
                textStyle: {{ color: '#334155' }},
                calculable: true,
                inRange: {{
                    color: ['#313695', '#4575b4', '#74add1', '#abd9e9', '#e0f3f8', '#fee090', '#fdae61', '#f46d43', '#d73027', '#a50026']
                }}
            }},
            xAxis: {{
                name: 'PC1 (${pc1_var} var)',
                nameLocation: 'middle',
                nameGap: 30,
                nameTextStyle: {{ color: '#475569', fontSize: 14, fontWeight: 'bold' }},
                axisLine: {{ lineStyle: {{ color: '#94A3B8' }} }},
                splitLine: {{ lineStyle: {{ color: '#F1F5F9', type: 'dashed' }} }},
                axisLabel: {{ color: '#475569' }}
            }},
            yAxis: {{
                name: 'PC2 (${pc2_var} var)',
                nameLocation: 'middle',
                nameGap: 40,
                nameTextStyle: {{ color: '#475569', fontSize: 14, fontWeight: 'bold' }},
                axisLine: {{ lineStyle: {{ color: '#94A3B8' }} }},
                splitLine: {{ lineStyle: {{ color: '#F1F5F9', type: 'dashed' }} }},
                axisLabel: {{ color: '#475569' }}
            }},
            series: [{{
                type: 'scatter',
                symbolSize: 18,
                data: {data_json},
                itemStyle: {{
                    shadowBlur: 4,
                    shadowColor: 'rgba(0, 0, 0, 0.15)'
                }},
                emphasis: {{
                    focus: 'none',
                    scale: 1.5,
                    itemStyle: {{
                        borderColor: '#0F172A',
                        borderWidth: 2,
                        shadowBlur: 10,
                        shadowColor: 'rgba(0, 0, 0, 0.3)'
                    }}
                }}
            }}]
        }};

        myChart.setOption(option);
        window.addEventListener('resize', () => myChart.resize());
    </script>
</body>
</html>"""

    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(html_content, encoding="utf-8")
    print(f"Fast interactive HTML PCA plot saved to {out_file}")


def main():
    args = parse_args()
    matrix_df = load_and_normalize_samples(args.files)
    
    if args.save_matrix:
        sep = "\t" if args.save_matrix.endswith(".tsv") else ","
        matrix_df.to_csv(args.save_matrix, sep=sep)
        print(f"Saved sample x kmer matrix to {args.save_matrix}")
        
    generate_echarts_html(
        matrix_df=matrix_df,
        output_path=args.output,
        random_state=args.random_state,
    )


if __name__ == "__main__":
    main()
