"""
FastAPI REST API for NGS Results Explorer.
Serves differential expression results from processed Parquet data.

Endpoints:
  GET /health          — service health + dataset stats
  GET /genes           — paginated gene list with filters
  GET /genes/{gene_id} — single gene detail
  GET /summary         — dataset-level summary statistics
  GET /top-expressed   — top N DE genes by fold change
  GET /volcano-data    — all data points for volcano plot
  GET /distribution    — log2FC bin distribution
  GET /quarantine      — quarantined records with rejection reasons
"""
from contextlib import asynccontextmanager
from typing import Optional
import json
import os

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import pandas as pd
import structlog

from src.processor import load_processed_data, load_summary
from src.config import settings

log = structlog.get_logger()

# ── Global data store (loaded once at startup) ──────────────────
_df: Optional[pd.DataFrame] = None
_summary: Optional[dict] = None


def get_df() -> pd.DataFrame:
    global _df
    if _df is None:
        _df = load_processed_data()
    return _df


def get_summary() -> dict:
    global _summary
    if _summary is None:
        _summary = load_summary()
    return _summary


# ── Lifespan — load data on startup ─────────────────────────────
@asynccontextmanager
async def lifespan(app: FastAPI):
    log.info("api_startup", version=settings.app_version)
    get_df()
    get_summary()
    log.info("data_loaded", rows=len(get_df()))
    yield
    log.info("api_shutdown")


# ── App ──────────────────────────────────────────────────────────
app = FastAPI(
    title="NGS Results Explorer API",
    description=(
        "REST API for querying RNA-Seq differential expression results. "
        "Dataset: GSE183947 — Human breast cancer tumour vs normal. "
        "Built to demonstrate production-grade bioinformatics data engineering."
    ),
    version=settings.app_version,
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


# ── Pydantic response schemas ────────────────────────────────────
class HealthResponse(BaseModel):
    status: str
    version: str
    total_genes: int
    significant_genes: int
    accession: str


class GeneRecord(BaseModel):
    gene_id: str
    gene_name: str
    log2fc: float
    pvalue: float
    padj: float
    base_mean: float
    regulation: str
    significant: bool
    neg_log10_pvalue: float
    abs_log2fc: float
    log2fc_bin: str
    completeness_score: float


class SummaryResponse(BaseModel):
    accession: str
    condition: str
    total_genes: int
    significant_genes: int
    upregulated: int
    downregulated: int
    not_significant: int
    pct_significant: float
    padj_threshold: float
    log2fc_threshold: float
    top_upregulated: list[dict]
    top_downregulated: list[dict]


class VolcanoPoint(BaseModel):
    gene_id: str
    gene_name: str
    log2fc: float
    neg_log10_pvalue: float
    regulation: str
    significant: bool
    padj: float


class DistributionBin(BaseModel):
    bin: str
    count: int


# ── Endpoints ────────────────────────────────────────────────────

@app.get("/", include_in_schema=False)
def root():
    return {
        "service": "NGS Results Explorer API",
        "version": settings.app_version,
        "docs": "/docs",
        "dataset": "GSE183947 — Breast cancer RNA-Seq",
    }


@app.get("/health", response_model=HealthResponse, tags=["System"])
def health():
    """Service health check with dataset statistics."""
    try:
        df = get_df()
        summary = get_summary()
        return HealthResponse(
            status="healthy",
            version=settings.app_version,
            total_genes=len(df),
            significant_genes=int(df["significant"].sum()),
            accession=summary["accession"],
        )
    except Exception as e:
        log.error("health_check_failed", error=str(e))
        raise HTTPException(status_code=503, detail="Data unavailable")


@app.get("/summary", response_model=SummaryResponse, tags=["Analysis"])
def summary():
    """Dataset-level differential expression summary statistics."""
    return get_summary()


@app.get("/genes", response_model=list[GeneRecord], tags=["Genes"])
def list_genes(
    regulation: Optional[str] = Query(
        None,
        description="Filter by regulation: upregulated, downregulated, not_significant"
    ),
    significant_only: bool = Query(
        False,
        description="Return only statistically significant genes"
    ),
    limit: int = Query(50, ge=1, le=500, description="Max records to return"),
    offset: int = Query(0, ge=0, description="Pagination offset"),
):
    """
    Paginated list of genes with optional filters.
    Sorted by significance then absolute fold change descending.
    """
    df = get_df()

    if regulation:
        valid = {"upregulated", "downregulated", "not_significant"}
        if regulation not in valid:
            raise HTTPException(
                status_code=400,
                detail=f"regulation must be one of {valid}"
            )
        df = df[df["regulation"] == regulation]

    if significant_only:
        df = df[df["significant"]]

    page = df.iloc[offset: offset + limit]
    return page.to_dict(orient="records")


@app.get("/genes/{gene_id}", response_model=GeneRecord, tags=["Genes"])
def get_gene(gene_id: str):
    """Full detail for a single gene by Ensembl ID."""
    df = get_df()
    match = df[df["gene_id"] == gene_id]
    if match.empty:
        raise HTTPException(
            status_code=404,
            detail=f"Gene '{gene_id}' not found"
        )
    return match.iloc[0].to_dict()


@app.get(
    "/top-expressed",
    response_model=list[GeneRecord],
    tags=["Analysis"]
)
def top_expressed(
    n: int = Query(10, ge=1, le=50, description="Number of top genes"),
    regulation: Optional[str] = Query(
        None,
        description="upregulated or downregulated"
    ),
):
    """Top N differentially expressed genes by absolute fold change."""
    df = get_df()
    df = df[df["significant"]]

    if regulation:
        df = df[df["regulation"] == regulation]

    top = df.nlargest(n, "abs_log2fc")
    return top.to_dict(orient="records")


@app.get(
    "/volcano-data",
    response_model=list[VolcanoPoint],
    tags=["Visualisation"]
)
def volcano_data():
    """
    All gene data points formatted for volcano plot rendering.
    x = log2FoldChange, y = -log10(pvalue)
    """
    df = get_df()
    cols = [
        "gene_id", "gene_name", "log2fc",
        "neg_log10_pvalue", "regulation", "significant", "padj"
    ]
    return df[cols].to_dict(orient="records")


@app.get(
    "/distribution",
    response_model=list[DistributionBin],
    tags=["Visualisation"]
)
def distribution():
    """Log2 fold-change bin distribution for histogram rendering."""
    df = get_df()
    bin_order = [
        "≤ -4", "-4 to -2", "-2 to -1",
        "-1 to 1", "1 to 2", "2 to 4", "≥ 4"
    ]
    counts = df["log2fc_bin"].value_counts()
    return [
        {"bin": b, "count": int(counts.get(b, 0))}
        for b in bin_order
    ]


@app.get("/quarantine", tags=["Quality"])
def quarantine_records():
    """
    Records that failed validation — preserved with rejection reasons.
    Quarantine-not-delete: no data is ever silently dropped.
    """
    path = os.path.join(settings.quarantine_path, "quarantine.csv")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path)
    return df.to_dict(orient="records")
