"""
Data Processor for NGS differential expression results.
Takes validated clean records and computes derived fields
used by the dashboard and API.

Derived fields:
  - significant: padj < 0.05 AND abs(log2fc) > 1.0
  - regulation: upregulated / downregulated / not_significant
  - log2fc_bin: binned fold-change for distribution plots
  - neg_log10_pvalue: for volcano plot y-axis
  - abs_log2fc: absolute fold change for ranking
"""
import os
import math
import json
import pandas as pd
import structlog

from src.config import settings

log = structlog.get_logger()

# DE significance thresholds — industry standard
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0  # = 2-fold change


def compute_derived_fields(record: dict) -> dict:
    """Add computed fields to a single clean record."""
    log2fc = float(record["log2fc"])
    pvalue = float(record["pvalue"])
    padj = float(record["padj"])

    # Significance: padj < 0.05 AND absolute fold change > 1 (2x)
    significant = (padj < PADJ_THRESHOLD) and (abs(log2fc) > LOG2FC_THRESHOLD)

    # Regulation direction
    if significant and log2fc > 0:
        regulation = "upregulated"
    elif significant and log2fc < 0:
        regulation = "downregulated"
    else:
        regulation = "not_significant"

    # Volcano plot y-axis: -log10(pvalue)
    # Guard against pvalue = 0 (use a floor)
    pvalue_floor = max(pvalue, 1e-300)
    neg_log10_pvalue = round(-math.log10(pvalue_floor), 4)

    # Absolute fold change for ranking
    abs_log2fc = round(abs(log2fc), 4)

    # Log2FC bins for distribution histogram
    if log2fc <= -4:
        log2fc_bin = "≤ -4"
    elif log2fc <= -2:
        log2fc_bin = "-4 to -2"
    elif log2fc <= -1:
        log2fc_bin = "-2 to -1"
    elif log2fc < 1:
        log2fc_bin = "-1 to 1"
    elif log2fc < 2:
        log2fc_bin = "1 to 2"
    elif log2fc < 4:
        log2fc_bin = "2 to 4"
    else:
        log2fc_bin = "≥ 4"

    return {
        **record,
        "significant": significant,
        "regulation": regulation,
        "neg_log10_pvalue": neg_log10_pvalue,
        "abs_log2fc": abs_log2fc,
        "log2fc_bin": log2fc_bin,
    }


def process_clean_records(clean_records: list[dict]) -> pd.DataFrame:
    """
    Process all clean records into an analysis-ready DataFrame.
    Saves to Parquet in the processed data directory.
    """
    log.info("processing_records", count=len(clean_records))

    enriched = [compute_derived_fields(r) for r in clean_records]
    df = pd.DataFrame(enriched)

    # Sort by significance then by absolute fold change descending
    df = df.sort_values(
        ["significant", "abs_log2fc"],
        ascending=[False, False],
    ).reset_index(drop=True)

    # Save to Parquet
    os.makedirs(settings.processed_data_path, exist_ok=True)
    parquet_path = os.path.join(
        settings.processed_data_path, "de_results.parquet"
    )
    df.to_parquet(parquet_path, index=False)
    log.info("parquet_saved", path=parquet_path, rows=len(df))

    # Also save summary stats to JSON for the API
    summary = _compute_summary(df)
    summary_path = os.path.join(
        settings.processed_data_path, "summary.json"
    )
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    log.info("summary_saved", path=summary_path)

    return df


def _compute_summary(df: pd.DataFrame) -> dict:
    """Compute dataset-level summary statistics."""
    total = len(df)
    significant = df["significant"].sum()
    upregulated = (df["regulation"] == "upregulated").sum()
    downregulated = (df["regulation"] == "downregulated").sum()
    not_significant = (df["regulation"] == "not_significant").sum()

    top_up = (
        df[df["regulation"] == "upregulated"]
        .nlargest(5, "abs_log2fc")[["gene_name", "log2fc", "padj"]]
        .to_dict(orient="records")
    )
    top_down = (
        df[df["regulation"] == "downregulated"]
        .nlargest(5, "abs_log2fc")[["gene_name", "log2fc", "padj"]]
        .to_dict(orient="records")
    )

    return {
        "accession": df["accession"].iloc[0] if total > 0 else "",
        "condition": df["condition"].iloc[0] if total > 0 else "",
        "total_genes": int(total),
        "significant_genes": int(significant),
        "upregulated": int(upregulated),
        "downregulated": int(downregulated),
        "not_significant": int(not_significant),
        "pct_significant": round(significant / total * 100, 2) if total else 0,
        "avg_completeness": round(
            df["completeness_score"].mean(), 4
        ) if "completeness_score" in df.columns else 0,
        "top_upregulated": top_up,
        "top_downregulated": top_down,
        "padj_threshold": PADJ_THRESHOLD,
        "log2fc_threshold": LOG2FC_THRESHOLD,
    }


def load_processed_data() -> pd.DataFrame:
    """Load the processed Parquet file. Used by API and dashboard."""
    parquet_path = os.path.join(
        settings.processed_data_path, "de_results.parquet"
    )
    if not os.path.exists(parquet_path):
        raise FileNotFoundError(
            "Processed data not found at {parquet_path}. "
            "Run the pipeline first: python3 -m src.processor"
        )
    return pd.read_parquet(parquet_path)


def load_summary() -> dict:
    """Load the summary JSON. Used by API and dashboard."""
    summary_path = os.path.join(
        settings.processed_data_path, "summary.json"
    )
    if not os.path.exists(summary_path):
        raise FileNotFoundError(
            "Summary not found at {summary_path}. "
            "Run the pipeline first: python3 -m src.processor"
        )
    with open(summary_path) as f:
        return json.load(f)


if __name__ == "__main__":
    import structlog
    structlog.configure()

    from src.fetcher import fetch_de_results
    from src.validator import validate_batch

    records = fetch_de_results("GSE183947")
    results = validate_batch(records)
    df = process_clean_records(results["clean"])

    summary = load_summary()

    print("\n=== Processing Summary ===")
    print(f"Total genes       : {summary['total_genes']}")
    print(f"Significant       : {summary['significant_genes']}")
    print(f"Upregulated       : {summary['upregulated']}")
    print(f"Downregulated     : {summary['downregulated']}")
    print(f"Not significant   : {summary['not_significant']}")
    print(f"% Significant     : {summary['pct_significant']}%")
    print("\nTop upregulated genes:")
    for g in summary["top_upregulated"]:
        print(f"  {g['gene_name']:12} log2FC={g['log2fc']:.2f}  padj={g['padj']:.2e}")
    print("\nTop downregulated genes:")
    for g in summary["top_downregulated"]:
        print(f"  {g['gene_name']:12} log2FC={g['log2fc']:.2f}  padj={g['padj']:.2e}")
    print("\nProcessed Parquet: data/processed/de_results.parquet")
