#!/usr/bin/env python3
import argparse
import base64
import csv
import html
import io
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


BARCODE_RANK_MAX_POINTS = 20000
BEAD_CHUNK_ROWS = 1_000_000
PARAMETER_DISPLAY_ORDER = [
    "sample_id",
    "fastq",
    "full_length_fastq",
    "tso_seq",
    "rtp_seq",
    "ref_dir",
    "out_dir",
    "threads",
    "cluster_threads",
    "exp_cells",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-html", required=True)
    parser.add_argument("--report-metrics-tsv", required=True)
    parser.add_argument("--rna-qc-metrics-tsv", required=True)
    parser.add_argument("--saturation-tsv", required=True)
    parser.add_argument("--read-qc-json", required=True)
    parser.add_argument("--parameters-tsv", required=True)
    parser.add_argument("--per-cell-qc-tsv", required=True)
    parser.add_argument("--barcode-counts-3p-tsv", default=None)
    parser.add_argument("--barcode-counts-5p-tsv", default=None)
    parser.add_argument("--whitelist-3p", default=None)
    parser.add_argument("--whitelist-5p", default=None)
    parser.add_argument("--read-assigned-cell", default=None)
    parser.add_argument("--glycine-stats", default=None)
    parser.add_argument("--skip-glycine", action="store_true")
    parser.add_argument("--knee-plot-3p", default=None)
    parser.add_argument("--knee-plot-5p", default=None)
    parser.add_argument("--saturation-png", default=None)
    parser.add_argument("--rna-violin-png", default=None)
    return parser.parse_args()


def read_tsv(path, required=True):
    if not path or not Path(path).exists():
        if required:
            raise FileNotFoundError(path)
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t").fillna("")


def value_to_float(value):
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text == "":
        return None
    if text.endswith("%"):
        text = text[:-1]
    try:
        return float(text)
    except ValueError:
        return None


def metric(report_df, qc_df, name, section=None):
    if not report_df.empty and {"Metric", "Value"}.issubset(report_df.columns):
        sub = report_df[report_df["Metric"].astype(str) == name]
        if section is not None and "Section" in sub.columns:
            sub = sub[sub["Section"].astype(str) == section]
        if not sub.empty:
            return value_to_float(sub.iloc[0]["Value"])
    if not qc_df.empty and {"Metric", "Value"}.issubset(qc_df.columns):
        sub = qc_df[qc_df["Metric"].astype(str) == name]
        if not sub.empty:
            return value_to_float(sub.iloc[0]["Value"])
    return None


def metric_text(report_df, qc_df, name, section=None):
    if not report_df.empty and {"Metric", "Formatted_value"}.issubset(report_df.columns):
        sub = report_df[report_df["Metric"].astype(str) == name]
        if section is not None and "Section" in sub.columns:
            sub = sub[sub["Section"].astype(str) == section]
        if not sub.empty and str(sub.iloc[0]["Formatted_value"]).strip():
            return str(sub.iloc[0]["Formatted_value"])
    value = metric(report_df, qc_df, name, section)
    return format_number(value)


def format_number(value, digits=2):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:,.{digits}f}"


def format_integer(value):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{int(round(value)):,}"


def format_percent(value, digits=2):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    return f"{value * 100:.{digits}f}%"


def compact_number(value):
    value = value_to_float(value)
    if value is None or not math.isfinite(value):
        return "NA"
    abs_value = abs(value)
    if abs_value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if abs_value >= 1_000:
        return f"{value / 1_000:.1f}K"
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.2f}"


def read_parameters(path):
    df = read_tsv(path, required=False)
    if df.empty or not {"Parameter", "Value"}.issubset(df.columns):
        return {}
    return {
        str(row["Parameter"]): str(row["Value"])
        for _, row in df.iterrows()
        if str(row["Parameter"]).strip()
    }


def html_table(rows, headers=("Metric", "Value"), class_name="report-table"):
    header_html = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = []
    for row in rows:
        body.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(h, '')))}</td>" for h in headers)
            + "</tr>"
        )
    return f"<table class=\"{class_name}\"><thead><tr>{header_html}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def load_report_template():
    template_path = Path(__file__).resolve().parent / "report_template.html"
    text = template_path.read_text(encoding="utf-8", errors="ignore")

    # The PBMC-derived template includes one malformed style/script splice and a
    # large "Relocated from body" CSS blob that breaks HTML parsing in some
    # browsers, causing a blank white report page.
    text = re.sub(
        r"""
        \s*/\*\s*Relocated\ from\ body\s*\*/.*?
        </script><script>"use\ strict";\s*
        </style>
        """,
        "\n    </style>\n",
        text,
        count=1,
        flags=re.S | re.X,
    )

    # The template shell should have a real <body>; some captured PBMC exports
    # omit it and rely on browser recovery, which is fragile once we inject our
    # own report body and scripts.
    head_end = text.find("</head>")
    if head_end != -1:
        body_after_head = re.search(r"<body\b", text[head_end:], flags=re.I)
        if body_after_head is None:
            text = text.replace("</head>", "</head>\n<body>", 1)

    return text


def pbmc_metric_rows(rows):
    body = "".join(
        f"""
        <tr>
          <td>{html.escape(str(row.get("Metric", "")))}</td>
          <td class="metric-value">{html.escape(str(row.get("Value", "")))}</td>
        </tr>
        """
        for row in rows
    )
    return f'<div class="stats-table-container"><table class="stats-table"><tbody>{body}</tbody></table></div>'


def pbmc_metric_rows_3col(rows):
    body = "".join(
        f"""
        <tr>
          <td>{html.escape(str(row.get("Metric", "")))}</td>
          <td class="metric-value read-count-value">{html.escape(str(row.get("Read count", "")))}</td>
          <td class="metric-value percent-value">{html.escape(str(row.get("Percent", "")))}</td>
        </tr>
        """
        for row in rows
    )
    return (
        '<div class="stats-table-container"><table class="stats-table sequencing-table">'
        '<colgroup><col class="metric-column"><col class="read-count-column">'
        '<col class="percent-column"></colgroup>'
        '<thead><tr><th>Metric</th><th>Read count</th><th>Percent</th></tr></thead>'
        f"<tbody>{body}</tbody></table></div>"
    )


def pbmc_summary_cards(card_rows):
    return "".join(
        f"""
        <div class="summary-card">
          <div class="summary-card-content">
            <div class="value">{html.escape(str(row.get("value", "NA")))}</div>
            <div class="label">{html.escape(str(row.get("label", "")))}</div>
          </div>
        </div>
        """
        for row in card_rows
    )


def summary_violin_cards(cards):
    blocks = []
    for card in cards:
        image_html = (
            f'<img src="{html.escape(card.get("src", ""))}" alt="{html.escape(card.get("title", ""))}" '
            'class="violin-image" />'
            if card.get("src")
            else f'<div class="empty-plot">{html.escape(card.get("title", "NA"))}</div>'
        )
        blocks.append(
            f"""
            <div class="summary-plot-col">
              <div class="violin-wrap">
                {image_html}
              </div>
            </div>
            """
        )
    return "".join(blocks)


def pbmc_title_block(title, help_id, help_html, sample_id=None):
    title_text = f"{title} : {sample_id}" if sample_id else title
    return f"""
        <div class="section-heading">
            <h2>{html.escape(title_text)}</h2>
            <button class="help-button" type="button" onclick="show('{html.escape(help_id)}')" aria-label="Show help">?</button>
            <div id="{html.escape(help_id)}" class="help-panel">
                {help_html}
            </div>
        </div>
    """


def pbmc_section_bar(title, help_id, help_html, width_px=240):
    return f"""
        <div class="section-heading">
            <h2>{html.escape(title)}</h2>
            <button class="help-button" type="button" onclick="show('{html.escape(help_id)}')" aria-label="Show help">?</button>
            <div id="{html.escape(help_id)}" class="help-panel">
                {help_html}
            </div>
        </div>
    """


def help_dl(items):
    parts = ['<dl class="help-list">']
    for title, desc in items:
        parts.append(f"<dt>{html.escape(title)}</dt>")
        parts.append(f"<dd>{html.escape(desc)}</dd>")
    parts.append("</dl>")
    return "".join(parts)


def build_sample_rows(params):
    rows = []
    for key in PARAMETER_DISPLAY_ORDER:
        if key in params and params[key] != "":
            rows.append({"Metric": key, "Value": params[key]})
    return rows


def parse_glycine_read_counts(path):
    if not path or not Path(path).exists():
        return None, None
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    values = {}
    for line in text.splitlines():
        parts = re.split(r"[\t: ]+", line.strip())
        if len(parts) < 2:
            continue
        key = parts[0].strip()
        for token in parts[1:]:
            if re.fullmatch(r"\d+(?:\.\d+)?", token):
                values.setdefault(key, float(token))
                break
    read_count = values.get("Read_count")
    length_filtered = values.get("Length-filtered", 0.0)
    qc_filtered = values.get("QC-filtered", 0.0)
    if read_count is None:
        return None, None
    return int(read_count), max(0, int(read_count - length_filtered - qc_filtered))


def parse_glycine_clean_reads(path):
    return parse_glycine_read_counts(path)[1]


def build_read_summary(report_df, qc_df, skip_glycine, glycine_stats):
    ordered = [
        ("Full length", "Full length"),
        ("Barcode-valid", "Barcode-valid"),
        ("Cell-assigned", "Cell-assigned"),
        ("Gene assigned", "Gene assigned"),
        ("Transcript assigned", "Transcript assigned"),
    ]
    full_length = metric(report_df, qc_df, "Full length", "Read assignment summary")
    if full_length is None:
        full_length = metric(report_df, qc_df, "Full length reads")
    glycine_raw_reads, glycine_clean_reads = parse_glycine_read_counts(glycine_stats)
    clean_reads = int(full_length) if skip_glycine and full_length is not None else glycine_clean_reads
    if clean_reads is None and full_length is not None:
        clean_reads = int(full_length)

    raw_reads = metric(report_df, qc_df, "Input reads")
    if raw_reads is None:
        raw_reads = full_length if skip_glycine else glycine_raw_reads
    if raw_reads is None:
        raw_reads = clean_reads or full_length

    rows = [{"Metric": "Raw reads", "Read count": format_number(raw_reads), "Percent": "100.00%"}]
    if clean_reads is not None:
        clean_ratio = clean_reads / raw_reads if raw_reads else None
        rows.append({
            "Metric": "Clean reads",
            "Read count": format_number(clean_reads),
            "Percent": format_percent(clean_ratio),
        })
    denominator = raw_reads
    for label, metric_name in ordered:
        count = metric(report_df, qc_df, metric_name, "Read assignment summary")
        percent = count / denominator if count is not None and denominator else None
        rows.append(
            {
                "Metric": label,
                "Read count": format_number(count),
                "Percent": format_percent(percent),
            }
        )
    return rows, denominator


def load_whitelist(path):
    if not path or not Path(path).exists():
        return set()
    return {
        line.strip().split(",")[0]
        for line in Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.strip()
    }


def count_data_rows(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        total = sum(1 for _ in handle)
    return max(0, total - 1)


def barcode_rank_payload(path, whitelist_path):
    if not path or not Path(path).exists():
        return None
    total = count_data_rows(path)
    if total == 0:
        return {"rank": [], "count": [], "threshold": None, "displayed_points": 0, "original_points": 0}
    if total <= BARCODE_RANK_MAX_POINTS:
        keep_ranks = set(range(1, total + 1))
    else:
        log_start = 0.0
        log_end = math.log(total)
        keep_ranks = {
            max(1, min(total, int(round(math.exp(log_start + (log_end - log_start) * i / (BARCODE_RANK_MAX_POINTS - 1))))))
            for i in range(BARCODE_RANK_MAX_POINTS)
        }
        keep_ranks.add(1)
        keep_ranks.add(total)

    whitelist = load_whitelist(whitelist_path)
    ranks = []
    counts = []
    threshold = None
    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rank, row in enumerate(reader, start=1):
            barcode = str(row.get("barcode", "")).strip()
            count = value_to_float(row.get("count", 0)) or 0.0
            if whitelist and barcode in whitelist:
                threshold = count if threshold is None else min(threshold, count)
            if rank in keep_ranks:
                ranks.append(rank)
                counts.append(count)
    return {
        "rank": ranks,
        "count": counts,
        "threshold": threshold,
        "displayed_points": len(ranks),
        "original_points": total,
    }


def beads_per_droplet_payload(path):
    if not path or not Path(path).exists():
        return {"x": [], "y": [], "n_cells": 0}
    cell_barcodes = defaultdict(set)
    usecols = None
    header = pd.read_csv(path, nrows=0)
    candidates = [c for c in ["cell_id", "BC5n", "BC3n"] if c in header.columns]
    if {"cell_id", "BC5n", "BC3n"}.issubset(candidates):
        usecols = candidates
    for chunk in pd.read_csv(path, chunksize=BEAD_CHUNK_ROWS, usecols=usecols):
        if not {"cell_id", "BC5n", "BC3n"}.issubset(chunk.columns):
            return {"x": [], "y": [], "n_cells": 0}
        for cell_id, bc5, bc3 in zip(chunk["cell_id"], chunk["BC5n"], chunk["BC3n"]):
            cell = str(cell_id).strip()
            if not cell:
                continue
            for barcode in (bc5, bc3):
                bc = str(barcode).strip()
                if bc and bc.lower() != "nan":
                    cell_barcodes[cell].add(bc)
    hist = Counter(len(v) for v in cell_barcodes.values())
    xs = sorted(hist)
    return {"x": xs, "y": [hist[x] for x in xs], "n_cells": len(cell_barcodes)}


def dataframe_payload(df):
    return {col: df[col].tolist() for col in df.columns}


def safe_read_json(path):
    if not path or not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def summary_rows(report_df, qc_df):
    items = [
        ("Estimated cells", "Estimated cells"),
        ("Input reads", "Input reads"),
        ("Reads per cell (mean)", "Reads per cell (mean)"),
        ("UMIs per cell (median)", "UMIs per cell (median)"),
        ("Genes per cell (median)", "Genes per cell (median)"),
        ("Unique genes", "Unique genes"),
        ("Unique isoforms", "Unique isoforms"),
    ]
    return [
        {"Metric": label, "Value": format_integer(metric(report_df, qc_df, name))}
        for label, name in items
    ]


def mapping_rows(report_df, qc_df, denominator):
    aligned = metric(report_df, qc_df, "Aligned BAM reads")
    unmapped = metric(report_df, qc_df, "Unmapped")
    if denominator is None:
        denominator = metric(report_df, qc_df, "Input reads")
    aligned_ratio = aligned / denominator if aligned is not None and denominator else None
    unmapped_ratio = unmapped / denominator if unmapped is not None and denominator else None
    return [
        {"Metric": "Aligned BAM reads / total reads", "Value": format_percent(aligned_ratio)},
        {"Metric": "Unmapped / total reads", "Value": format_percent(unmapped_ratio)},
        {"Metric": "Unique genes", "Value": metric_text(report_df, qc_df, "Unique genes")},
        {"Metric": "Unique isoforms", "Value": metric_text(report_df, qc_df, "Unique isoforms")},
    ]


def per_cell_payload(path):
    df = read_tsv(path, required=False)
    if df.empty:
        return {"reads": [], "umis": [], "genes": [], "mito_percent": []}
    keep = min(len(df), 20000)
    if len(df) > keep:
        df = df.sample(n=keep, random_state=1)
    payload = {}
    for col in ["reads", "umis", "genes", "mito_percent"]:
        payload[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce").fillna(0).tolist()
    return payload


def violin_data_uri(values, color, title, ylabel):
    series = pd.to_numeric(pd.Series(values, dtype=float), errors="coerce")
    series = series[series.notna()]
    if series.empty:
        return ""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return ""

    fig, ax = plt.subplots(figsize=(3.5, 5.6), dpi=150)
    violin = ax.violinplot(series.tolist(), showmeans=True, showmedians=True, showextrema=False, widths=0.75)
    for body in violin["bodies"]:
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.70)
    if "cmeans" in violin:
        violin["cmeans"].set_color("#333333")
        violin["cmeans"].set_linewidth(1.2)
    if "cmedians" in violin:
        violin["cmedians"].set_color("#111111")
        violin["cmedians"].set_linewidth(1.3)

    ax.set_title(title, fontsize=11, pad=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_xticks([])
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8)
    ax.set_facecolor("white")
    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#888888")
    fig.patch.set_facecolor("white")
    fig.tight_layout()

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def new_report_markup(sections):
    return f"""
<style>
  .report-section {{ margin-bottom: 1.5rem; }}
  .section-heading {{
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: 0.75rem;
    margin-bottom: 1.25rem;
  }}
  .section-heading h2 {{
    margin: 0;
    font-size: 1.35rem;
    font-weight: 650;
    color: var(--text-primary);
  }}
  .help-button {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 28px;
    height: 28px;
    padding: 0;
    border: 0;
    border-radius: 50%;
    color: #fff;
    background: var(--gradient-bg);
    box-shadow: 0 2px 8px rgba(32, 85, 138, 0.25);
    font-weight: 700;
  }}
  .help-panel {{
    display: none;
    flex-basis: 100%;
    padding: 1rem 1.1rem;
    border: 1px solid var(--border-color);
    border-radius: 10px;
    background: var(--bg-surface-hover);
    color: var(--text-secondary);
  }}
  .help-list {{ margin: 0; }}
  .help-list dt {{
    margin-top: 0.7rem;
    color: var(--text-primary);
    font-weight: 650;
  }}
  .help-list dt:first-child {{ margin-top: 0; }}
  .help-list dd {{ margin: 0.15rem 0 0; }}
  .summary-cards {{ margin-bottom: 1.5rem; }}
  .report-grid {{
    display: grid;
    gap: 1rem;
    align-items: stretch;
  }}
  .report-grid.two {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
  .report-grid.three {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
  .summary-detail-grid {{
    display: grid;
    grid-template-columns: minmax(260px, 0.8fr) minmax(0, 1.7fr);
    gap: 1rem;
    align-items: stretch;
  }}
  .summary-violin-grid {{
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: 0.75rem;
  }}
  .plot-panel,
  .table-panel {{
    box-sizing: border-box;
    min-width: 0;
    max-width: 100%;
    padding: 1rem;
    border: 1px solid var(--border-color);
    border-radius: 10px;
    background: var(--bg-surface);
    box-shadow: var(--shadow-sm);
  }}
  .plot-panel {{ overflow: hidden; }}
  .stats-table-container {{ width: 100%; overflow-x: auto; }}
  .stats-table {{
    width: 100%;
    border-collapse: collapse;
    table-layout: fixed;
  }}
  .stats-table th,
  .stats-table td {{
    padding: 0.75rem 0.65rem;
    text-align: left;
    vertical-align: middle;
  }}
  .stats-table thead tr,
  .stats-table tbody tr {{ border-bottom: 1px solid var(--border-color); }}
  .stats-table th {{
    color: var(--text-secondary);
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}
  .stats-table tbody tr:last-child {{ border-bottom: 0; }}
  .stats-table .metric-value {{
    text-align: right;
    color: var(--primary-dark);
    font-weight: 650;
    white-space: nowrap;
  }}
  .stats-table td.metric-value {{ display: table-cell; }}
  .sequencing-table .metric-column {{ width: 50%; }}
  .sequencing-table .read-count-column {{ width: 30%; }}
  .sequencing-table .percent-column {{ width: 20%; }}
  .sequencing-table th:nth-child(2),
  .sequencing-table th:nth-child(3) {{ text-align: right; white-space: nowrap; }}
  .sequencing-table .read-count-value,
  .sequencing-table .percent-value {{ white-space: nowrap; }}
  .dynamic-note {{
    margin-top: 0.5rem;
    color: var(--text-muted);
    font-size: 0.8rem;
  }}
  .dynamic-plot,
  .dynamic-plot-wide,
  .dynamic-plot-lg {{
    box-sizing: border-box;
    width: 100%;
    max-width: 100%;
    min-width: 0;
    overflow: hidden;
  }}
  .dynamic-plot .plot-container,
  .dynamic-plot-wide .plot-container,
  .dynamic-plot-lg .plot-container,
  .dynamic-plot .svg-container,
  .dynamic-plot-wide .svg-container,
  .dynamic-plot-lg .svg-container {{ max-width: 100%; }}
  .dynamic-plot {{ height: 360px; }}
  .dynamic-plot-wide,
  .dynamic-plot-lg {{ height: 400px; }}
  .violin-wrap {{
    display: flex;
    align-items: center;
    justify-content: center;
    min-height: 390px;
    max-width: 100%;
    overflow: hidden;
  }}
  .violin-image {{
    display: block;
    width: 100%;
    max-width: 100%;
    height: 390px;
    object-fit: contain;
  }}
  .empty-plot {{ color: var(--text-muted); }}
  @media (max-width: 900px) {{
    .report-grid.two,
    .report-grid.three,
    .summary-detail-grid,
    .summary-violin-grid {{ grid-template-columns: 1fr; }}
  }}
</style>
<div class="tab-content">
  <div class="tab-pane active" id="summary-tab">
    <section class="summary-section report-section" data-library-content="gene-expression">
      {sections["sample_title"]}
      <div class="summary-cards">{sections["summary_cards"]}</div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["read_qc_bar"]}
      <div class="report-grid three">
        <div class="plot-panel"><div id="read-qc-quality" class="dynamic-plot"></div></div>
        <div class="plot-panel"><div id="read-qc-length" class="dynamic-plot"></div></div>
        <div class="plot-panel"><div id="read-qc-yield" class="dynamic-plot"></div></div>
      </div>
    </section>
  </div>

  <div class="tab-pane" id="cells-tab">
    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["summary_bar"]}
      <div class="summary-detail-grid">
        <div class="table-panel">{sections["summary_rows"]}</div>
        <div class="summary-violin-grid">{sections["summary_violin_plots"]}</div>
      </div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["beads_bar"]}
      <div class="report-grid two">
        <div class="plot-panel">
          <div id="barcode-rank" class="dynamic-plot"></div>
          <div class="dynamic-note" id="barcode-rank-note"></div>
        </div>
        <div class="plot-panel">
          <div id="beads-per-droplet" class="dynamic-plot"></div>
          <div class="dynamic-note" id="beads-note"></div>
        </div>
      </div>
    </section>
  </div>

  <div class="tab-pane" id="library-tab">
    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["sequencing_bar"]}
      <div class="report-grid two">
        <div class="table-panel" id="ReadSummaryTable">{sections["read_summary_rows"]}</div>
        <div class="plot-panel"><div id="read-assignment-plot" class="dynamic-plot-wide"></div></div>
      </div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["mapping_bar"]}
      <div class="table-panel" id="MappingTable">{sections["mapping_rows"]}</div>
    </section>

    <section class="section-card report-section" data-library-content="gene-expression">
      {sections["saturation_bar"]}
      <div class="report-grid three">
        <div class="plot-panel"><div id="saturation-genes" class="dynamic-plot-lg"></div></div>
        <div class="plot-panel"><div id="saturation-umis" class="dynamic-plot-lg"></div></div>
        <div class="plot-panel"><div id="saturation-rate" class="dynamic-plot-lg"></div></div>
      </div>
    </section>
  </div>
</div>
"""


def build_html(args, payload, sections):
    report_body = new_report_markup(sections) + f"""
<script>
  const payload = {json.dumps(payload, ensure_ascii=False)};
  const plotConfig = {{
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["autoScale2d", "hoverClosestCartesian", "hoverCompareCartesian", "lasso2d", "sendDataToCloud", "toggleSpikelines", "zoomIn2d", "zoomOut2d"]
  }};
  const baseTemplate = {{
    layout: {{
      colorway: ["#636efa", "#EF553B", "#00cc96", "#ab63fa", "#FFA15A", "#19d3f3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"],
      font: {{color: "#2a3f5f"}},
      paper_bgcolor: "white",
      plot_bgcolor: "#E5ECF6",
      hoverlabel: {{align: "left"}},
      title: {{x: 0.05}}
    }}
  }};
  const baseLayout = {{
    autosize: true,
    margin: {{t: 60, l: 60, r: 20, b: 55}},
    template: baseTemplate,
    plot_bgcolor: "rgba(0,0,0,0)",
    paper_bgcolor: "white",
    font: {{color: "#2a3f5f"}},
    xaxis: {{
      automargin: true,
      gridcolor: "lightgray",
      zeroline: true,
      zerolinecolor: "gray",
      zerolinewidth: 1,
      title: {{font: {{size: 12}}, standoff: 15}}
    }},
    yaxis: {{
      automargin: true,
      gridcolor: "lightgray",
      zeroline: true,
      zerolinecolor: "gray",
      zerolinewidth: 1,
      title: {{font: {{size: 12}}, standoff: 15}}
    }},
    legend: {{tracegroupgap: 0}}
  }};
  function show(id) {{
    const el = document.getElementById(id);
    if (!el) return;
    el.style.display = el.style.display === "none" || !el.style.display ? "block" : "none";
  }}
  document.querySelectorAll('.nav-tabs a').forEach(link => {{
    link.addEventListener('click', ev => {{
      ev.preventDefault();
      document.querySelectorAll('.nav-tabs li').forEach(li => li.classList.remove('active'));
      link.parentElement.classList.add('active');
      document.querySelectorAll('#myTabContent .tab-pane').forEach(p => p.classList.remove('active', 'in'));
      const target = document.querySelector(link.getAttribute('href'));
      if (target) target.classList.add('active', 'in');
      setTimeout(() => window.dispatchEvent(new Event('resize')), 50);
    }});
  }});
  function plotIf(id, data, layout) {{
    const el = document.getElementById(id);
    if (!el) return;
    Plotly.newPlot(id, data, Object.assign({{}}, baseLayout, layout || {{}}), plotConfig);
    requestAnimationFrame(() => {{
      if (el.offsetParent !== null) Plotly.Plots.resize(el);
    }});
  }}
  function resizeReportPlots(container) {{
    (container || document).querySelectorAll(".js-plotly-plot").forEach(plot => {{
      if (plot.offsetParent !== null) Plotly.Plots.resize(plot);
    }});
  }}
  window.addEventListener("resize", () => resizeReportPlots(document));
  if (window.ResizeObserver) {{
    const reportPlotObserver = new ResizeObserver(entries => {{
      entries.forEach(entry => resizeReportPlots(entry.target));
    }});
    document.querySelectorAll(".plot-panel").forEach(panel => reportPlotObserver.observe(panel));
  }}
  function pbmcLineLayout(title, xTitle, yTitle, extra) {{
    return Object.assign({{
      title: {{text: title, font: {{size: 14}}}},
      xaxis: {{title: {{text: xTitle, font: {{size: 12}}}}}},
      yaxis: {{title: {{text: yTitle}}}}
    }}, extra || {{}});
  }}
  plotIf("read-assignment-plot", [{{
    x: payload.readSummary.labels,
    y: payload.readSummary.counts,
    type: "bar",
    marker: {{color: "#337ab7"}},
    hovertemplate: "%{{x}}<br>%{{y:,.0f}} reads<extra></extra>"
  }}], {{
    height: 420,
    title: {{text: "Read assignment summary", font: {{size: 14}}}},
    xaxis: {{title: {{text: ""}}}},
    yaxis: {{title: {{text: "Reads"}}}}
  }});

  const rq = payload.readQc || {{}};
  if (rq.quality) {{
    plotIf("read-qc-quality", [{{x: rq.quality.bins || [], y: rq.quality.counts || [], type: "bar", marker: {{color: "#337ab7"}}}}], pbmcLineLayout("Read quality", "Mean Q score", "Reads", {{
      height: 360
    }}));
  }}
  if (rq.length) {{
    plotIf("read-qc-length", [{{x: rq.length.bins_kb || [], y: rq.length.counts || [], type: "bar", marker: {{color: "#169c9c"}}}}], pbmcLineLayout("Read length", "Read length (kb)", "Reads", {{
      height: 360
    }}));
  }}
  if (rq.yield_above_length) {{
    plotIf("read-qc-yield", [{{x: rq.yield_above_length.x_kb || [], y: rq.yield_above_length.y_gb || [], type: "scatter", mode: "lines", line: {{color: "#e18435", width: 3}}}}], pbmcLineLayout("Base yield above read length", "Read length cutoff (kb)", "Yield above cutoff (Gb)", {{
      height: 360
    }}));
  }}

  const rankData = [];
  if (payload.barcodeRank3p) {{
    rankData.push({{
      x: payload.barcodeRank3p.rank,
      y: payload.barcodeRank3p.count,
      type: "scatter",
      mode: "lines",
      name: "3' barcode",
      line: {{color: "#1358A2", width: 3}},
      hovertemplate: "3' barcode<br>Rank: %{{x:,.0f}}<br>Reads: %{{y:,.0f}}<extra></extra>"
    }});
    if (payload.barcodeRank3p.threshold) {{
      rankData.push({{
        x: payload.barcodeRank3p.rank,
        y: payload.barcodeRank3p.rank.map(() => payload.barcodeRank3p.threshold),
        type: "scatter",
        mode: "lines",
        name: "3' threshold",
        line: {{color: "#1358A2", width: 1.5, dash: "dot"}},
        hovertemplate: "3' threshold<br>Reads: %{{y:,.0f}}<extra></extra>"
      }});
    }}
  }}
  if (payload.barcodeRank5p) {{
    rankData.push({{
      x: payload.barcodeRank5p.rank,
      y: payload.barcodeRank5p.count,
      type: "scatter",
      mode: "lines",
      name: "5' barcode",
      line: {{color: "#EF7C23", width: 3}},
      hovertemplate: "5' barcode<br>Rank: %{{x:,.0f}}<br>Reads: %{{y:,.0f}}<extra></extra>"
    }});
    if (payload.barcodeRank5p.threshold) {{
      rankData.push({{
        x: payload.barcodeRank5p.rank,
        y: payload.barcodeRank5p.rank.map(() => payload.barcodeRank5p.threshold),
        type: "scatter",
        mode: "lines",
        name: "5' threshold",
        line: {{color: "#EF7C23", width: 1.5, dash: "dot"}},
        hovertemplate: "5' threshold<br>Reads: %{{y:,.0f}}<extra></extra>"
      }});
    }}
  }}
  plotIf("barcode-rank", rankData, {{
    height: 360,
    title: {{text: "Barcode rank plot", font: {{size: 14}}}},
    xaxis: {{
      title: {{text: "Barcode in Rank-descending Order", font: {{size: 12}}}},
      type: "log",
      fixedrange: true,
      color: "black",
      linecolor: "black",
      linewidth: 1,
      showline: true
    }},
    yaxis: {{
      title: {{text: "Read counts"}},
      type: "log",
      fixedrange: true,
      color: "black",
      linecolor: "black",
      linewidth: 1,
      showline: true
    }},
    legend: {{orientation: "h", x: 0.02, y: 1.12}}
  }});
  const rankNote = document.getElementById("barcode-rank-note");
  if (rankNote) {{
    const n3 = payload.barcodeRank3p ? payload.barcodeRank3p.original_points : 0;
    const n5 = payload.barcodeRank5p ? payload.barcodeRank5p.original_points : 0;
    rankNote.textContent = `Displayed with log-rank sampling. 3' points: ${{n3.toLocaleString()}}; 5' points: ${{n5.toLocaleString()}}.`;
  }}

  const beadPalette = ["rgb(102,194,165)", "rgb(252,141,98)", "rgb(141,160,203)", "rgb(231,138,195)", "rgb(166,216,84)", "rgb(255,217,47)", "rgb(229,196,148)", "rgb(179,179,179)"];
  const beadData = (payload.beads.x || []).map((x, i) => ({{
    x: [String(x)],
    y: [payload.beads.y[i] || 0],
    type: "bar",
    width: 0.9,
    name: `${{x}}  ${{payload.beads.y[i] || 0}}`,
    marker: {{color: beadPalette[i % beadPalette.length]}},
    hovertemplate: "Merge Beads Num: %{{x}}<br>Count: %{{y:,.0f}}<extra></extra>"
  }}));
  plotIf("beads-per-droplet", beadData, {{
    height: 360,
    barmode: "relative",
    title: {{text: `Total cell number ${{(payload.beads.n_cells || 0).toLocaleString()}}`, font: {{color: "black", family: "Arial", size: 18}}, x: 0.05, y: 0.95}},
    margin: {{t: 50, l: 40, r: 10, b: 50}},
    xaxis: {{title: {{text: "Number of beads per droplet", font: {{size: 12}}, standoff: 10}}, showgrid: false}},
    yaxis: {{title: {{text: "Count", font: {{size: 13}}, standoff: 10}}, showgrid: false}},
    legend: {{font: {{family: "Arial", size: 10}}, x: 0.8, y: 1}}
  }});
  const beadsNote = document.getElementById("beads-note");
  if (beadsNote) beadsNote.textContent = `${{(payload.beads.n_cells || 0).toLocaleString()}} cells summarized from read_assigned_cell.csv.`;

  const sat = payload.saturation || {{}};
  plotIf("saturation-genes", [{{x: sat.reads_per_cell || [], y: sat.genes_per_cell || [], type: "scatter", mode: "lines", line: {{color: "#337ab7", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Median Genes per Cell", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Median Genes per Cell", standoff: 0}}}}
  }});
  plotIf("saturation-umis", [{{x: sat.reads_per_cell || [], y: sat.umis_per_cell || [], type: "scatter", mode: "lines", line: {{color: "#337ab7", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Median UMI counts per cell", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Median UMI Counts per Cell"}}}}
  }});
  plotIf("saturation-rate", [{{x: sat.reads_per_cell || [], y: (sat.saturation || []).map(v => v * 100), type: "scatter", mode: "lines", line: {{color: "#337ab7", width: 3}}, showlegend: false, hovertemplate: "x=%{{x}}<br>y=%{{y}}<extra></extra>"}}], {{
    height: 400,
    title: {{text: "Sequencing saturation", font: {{size: 14}}}},
    xaxis: {{title: {{text: "Mean Reads per Cell", font: {{size: 12}}}}, tickformat: "~s", hoverformat: ",.0f"}},
    yaxis: {{title: {{text: "Sequencing Saturation"}}, range: [0, 100]}}
  }});
</script>
"""
    template = load_report_template()
    return (
        template
        .replace("__SAMPLE__", html.escape(args.sample_id))
        .replace("__REPORT_BODY__", report_body)
    )


def main():
    args = parse_args()
    report_df = read_tsv(args.report_metrics_tsv)
    qc_df = read_tsv(args.rna_qc_metrics_tsv)
    saturation_df = read_tsv(args.saturation_tsv)
    params = read_parameters(args.parameters_tsv)
    read_qc = safe_read_json(args.read_qc_json)

    read_summary_rows, total_reads = build_read_summary(report_df, qc_df, args.skip_glycine, args.glycine_stats)
    per_cell = per_cell_payload(args.per_cell_qc_tsv)
    summary_violin_rows = summary_violin_cards([
        {
            "title": "Reads per cell",
            "src": violin_data_uri(per_cell.get("reads", []), "#337ab7", "Reads per cell", "Reads"),
        },
        {
            "title": "UMIs per cell",
            "src": violin_data_uri(per_cell.get("umis", []), "#169c9c", "UMIs per cell", "UMIs"),
        },
        {
            "title": "Genes per cell",
            "src": violin_data_uri(per_cell.get("genes", []), "#4f9d69", "Genes per cell", "Genes"),
        },
    ])

    summary = summary_rows(report_df, qc_df)
    summary_card_rows = []
    for source_label, display_label in [
        ("Estimated cells", "Estimated number of cells"),
        ("UMIs per cell (median)", "Median UMI counts per cell"),
        ("Genes per cell (median)", "Median genes per cell"),
        ("Reads per cell (mean)", "Mean reads per cell"),
    ]:
        value = next((row["Value"] for row in summary if row["Metric"] == source_label), "NA")
        summary_card_rows.append({"label": display_label, "value": value})

    payload = {
        "readSummary": {
            "labels": [row["Metric"] for row in read_summary_rows if row["Metric"] not in {"Raw reads", "Clean reads"}],
            "counts": [value_to_float(row["Read count"]) or 0 for row in read_summary_rows if row["Metric"] not in {"Raw reads", "Clean reads"}],
        },
        "readQc": read_qc,
        "saturation": dataframe_payload(saturation_df) if not saturation_df.empty else {},
        "barcodeRank3p": barcode_rank_payload(args.barcode_counts_3p_tsv, args.whitelist_3p),
        "barcodeRank5p": barcode_rank_payload(args.barcode_counts_5p_tsv, args.whitelist_5p),
        "beads": beads_per_droplet_payload(args.read_assigned_cell),
    }

    sections = {
        "sample_title": pbmc_title_block(
            "Summary",
            "sample-information",
            help_dl([
                ("Estimated number of cells", "The number of barcodes identified as real cells in the sequencing data after barcode merging and cell calling."),
                ("Median UMI counts per cell", "The median number of unique molecular identifiers detected per cell among all identified cells."),
                ("Median genes per cell", "The median number of unique genes with detectable expression in each cell."),
                ("Mean reads per cell", "The average number of sequencing reads per cell, calculated by dividing the total reads assigned to cells by the number of identified cells."),
            ]),
        ),
        "summary_cards": pbmc_summary_cards(summary_card_rows),
        "summary_bar": pbmc_section_bar(
            "Summary",
            "sumary-detail",
            help_dl([
                ("Left", "Core run-level summary values calculated from the final cell set."),
                ("Right", "Per-cell read, UMI, and gene distributions drawn from per-cell QC outputs."),
            ]),
            width_px=190,
        ),
        "summary_rows": pbmc_metric_rows(summary),
        "summary_violin_plots": summary_violin_rows,
        "read_qc_bar": pbmc_section_bar(
            "Read QC",
            "read-qc-detail",
            help_dl([
                ("Read quality", "Distribution of per-read mean base quality scores."),
                ("Read length", "Distribution of full-length read lengths."),
                ("Base yield above read length", "Total retained bases after applying each minimum read-length cutoff."),
            ]),
            width_px=180,
        ),
        "beads_bar": pbmc_section_bar(
            "Beads to cells",
            "bead-detail",
            help_dl([
                ("Left", "Barcode rank plot showing abundance-ranked corrected barcodes. The 3' and 5' barcode curves are overlaid in different colors."),
                ("Right", "Distribution of the number of unique corrected barcodes associated with each final cell."),
            ]),
            width_px=230,
        ),
        "sequencing_bar": pbmc_section_bar(
            "Sequencing",
            "sequence-detail",
            help_dl([
                ("Raw reads", "Total reads in --fastq when Glycine runs, or in --full-length-fastq when Glycine is skipped. This is the denominator for every percentage in the table."),
                ("Clean reads", "Input reads retained after glycine filtering when enabled; otherwise equal to full-length reads."),
                ("Full length", "Reads retained as full-length cDNA reads before downstream cell assignment."),
                ("Barcode-valid / Cell-assigned / Gene assigned / Transcript assigned", "Read counts carried through barcode correction, cell assignment, gene assignment, and transcript assignment."),
            ]),
            width_px=180,
        ),
        "read_summary_rows": pbmc_metric_rows_3col(read_summary_rows),
        "mapping_bar": pbmc_section_bar(
            "Mapping & Annotation",
            "mapping-detail",
            help_dl([
                ("Aligned BAM reads / total reads", "Fraction of total reads that are present in the aligned BAM output."),
                ("Unmapped / total reads", "Fraction of total reads that remain unmapped."),
                ("Unique genes / Unique isoforms", "Numbers of unique genes and isoforms detected in the final outputs."),
            ]),
            width_px=320,
        ),
        "mapping_rows": pbmc_metric_rows(mapping_rows(report_df, qc_df, total_reads)),
        "saturation_bar": pbmc_section_bar(
            "Saturation",
            "saturation-detail",
            help_dl([
                ("Median genes per cell", "Downsampling preview of median genes per cell using the same known-gene definition as the summary metric."),
                ("Median UMI counts per cell", "Downsampling preview of median UMI counts per cell."),
                ("Sequencing saturation", "Estimated saturation as sequencing depth increases."),
            ]),
            width_px=190,
        ),
    }
    html_text = build_html(args, payload, sections)
    Path(args.output_html).write_text(html_text, encoding="utf-8")
    print(f"Wrote report: {args.output_html}")


if __name__ == "__main__":
    main()
