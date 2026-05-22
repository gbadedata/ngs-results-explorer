<div align="center">

# NGS Results Explorer

**A production-grade RNA-Seq differential expression pipeline with validation engine, REST API, and interactive volcano plot dashboard.**

Built to demonstrate real-world bioinformatics data engineering - ingesting, validating, processing, and visualising NGS analysis results exactly as commercial platforms deliver to researchers.

[![CI](https://github.com/gbadedata/ngs-results-explorer/actions/workflows/ci.yml/badge.svg)](https://github.com/gbadedata/ngs-results-explorer/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-009688?logo=fastapi&logoColor=white)
![Plotly](https://img.shields.io/badge/Plotly_Dash-3F4F75?logo=plotly&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-150458?logo=pandas&logoColor=white)
![pytest](https://img.shields.io/badge/36_tests-passing-38A169)

**Dataset:** GSE183947 · Human breast cancer tumour vs normal · NCBI GEO

</div>

---

## The Problem

Every RNA-Seq pipeline produces the same problem downstream. After alignment and quantification, researchers receive raw differential expression results - gene identifiers, fold changes, adjusted p-values - typically as flat files with no validation, no quality enforcement, and no way to explore the data without writing code. Bad records enter analyses silently. Biologically implausible values go undetected. Bench scientists cannot query results without help from a bioinformatician.

NGS Results Explorer solves this by building a rigorous data engineering layer directly on top of DE results. It fetches real RNA-Seq data from NCBI GEO, applies nine biological validation rules to every gene record, routes failed records to a quarantine zone with precise rejection reasons, structures clean records into an analysis-ready Parquet dataset with computed fields, and exposes everything through a self-documenting REST API and an interactive browser-based dashboard complete with a volcano plot, fold-change distribution, and sortable gene table.

---

## Architecture

```text
NCBI GEO  (GSE183947 · Human breast cancer RNA-Seq)
     │
     ▼
GEO Fetcher  (requests, tenacity, exponential back-off, local cache)
     │
     ▼
Raw gene records  (gene_id, gene_name, log2fc, pvalue, padj, base_mean)
     │
     ▼
Validation Engine  -  9 rules across 4 categories
     ├── Passed (51 genes)  →  Parquet  (data/processed/)
     └── Failed  (2 genes)  →  Quarantine CSV + rejection reason
          │
          ▼
     Data Processor
     Computes: significant flag, regulation direction,
     neg_log10_pvalue, log2FC bins, abs_log2fc
          │
          ├──────────────────────────────┐
          ▼                              ▼
     FastAPI REST API             Plotly Dash Dashboard
     8 endpoints                  Volcano plot
     Swagger /docs                Top DE genes table
     Pydantic schemas             log2FC distribution
     Pagination + filters         QC summary panel
          │
          ▼
     GitHub Actions CI  (flake8 + pytest · 36 tests passing)
```

---

## Validation Engine

Nine rules applied to every gene record. Rules are pure Python functions - independently testable, composable, and extensible without changing the core validation logic.

| Category | Rule | What It Catches |
|---|---|---|
| Identity | gene_id not null | Missing primary identifier |
| Identity | Ensembl ID format (ENSG + 11 digits) | Raw internal NCBI UIDs leaking as gene IDs |
| Identity | gene_name not null | Missing gene symbol |
| Statistics | pvalue in [0, 1] | Impossible p-values |
| Statistics | padj in [0, 1] | Impossible adjusted p-values |
| Statistics | log2FC in [−50, 50] | Biologically implausible fold changes |
| Metrics | base_mean ≥ 0 | Negative read counts - impossible by definition |
| Completeness | condition present | Missing experimental context |
| Completeness | accession present | Missing study provenance |

Failed records are **quarantined with rejection reasons - never silently dropped.** In this run, 2 of 53 records were quarantined: one raw NCBI internal UID leaking as a gene ID, and one record with a missing gene symbol.

---

## Results

| Metric | Value | Notes |
|---|---|---|
| Records ingested | 53 | NCBI GEO accession GSE183947 |
| Passed validation | 51 (96.23%) | Written to Parquet |
| Quarantined | 2 (3.77%) | Preserved with rejection reasons |
| Significant DE genes | 20 (39.22%) | padj < 0.05 and |log2FC| > 1.0 |
| Upregulated | 12 | ERBB2, MYC, ESR1, CCNE1, CDC20 ... |
| Downregulated | 8 | BRCA1, PTEN, NF1, BRCA2, APC ... |
| Avg completeness score | 100% | All clean records fully populated |
| Validation rules | 9 | Across 4 categories |
| Unit tests | 36 passing | 84% coverage on the validation engine |

The biological results validate the pipeline. Every upregulated gene is a well-established oncogene in breast cancer - ERBB2 defines the HER2-amplified subtype, MYC and CCNE1 are canonical cell-cycle drivers, ESR1 encodes the oestrogen receptor central to luminal breast cancer biology. Every downregulated gene is a known tumour suppressor - BRCA1, BRCA2, and PTEN are among the most studied in the field.

---

## API Endpoints

| Method | Endpoint | Description |
|---|---|---|
| GET | /health | Service health check with dataset statistics |
| GET | /summary | Dataset-level differential expression summary |
| GET | /genes | Paginated gene list with regulation and significance filters |
| GET | /genes/{gene_id} | Full record for a single gene by Ensembl ID |
| GET | /top-expressed | Top N DE genes ranked by absolute fold change |
| GET | /volcano-data | All data points formatted for volcano plot rendering |
| GET | /distribution | log2FC bin distribution for histogram rendering |
| GET | /quarantine | Quarantined records with rejection reasons |

Interactive Swagger documentation available at `http://localhost:8000/docs`.

---

## Running Locally

```bash
# Clone and install
git clone https://github.com/gbadedata/ngs-results-explorer.git
cd ngs-results-explorer
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Run the full pipeline
python3 -m src.fetcher      # Fetch GEO data and cache locally
python3 -m src.validator    # Validate all records, write quarantine
python3 -m src.processor    # Compute derived fields, write Parquet

# Start the API
uvicorn src.api:app --reload --port 8000
# Open http://localhost:8000/docs

# Start the dashboard
python3 -m src.dashboard
# Open http://localhost:8050

# Run tests
pytest -v
```

---

## Technology Stack

| Category | Technology | Role |
|---|---|---|
| Language | Python 3.12 | Primary language throughout |
| Data processing | pandas 2.2, pyarrow 18 | DataFrame operations and Parquet serialisation |
| API | FastAPI 0.109, Pydantic v2 | REST endpoints, schema validation, OpenAPI docs |
| Dashboard | Plotly Dash 2.14, DBC 1.5 | Volcano plot, distribution chart, gene table |
| HTTP client | requests 2.31, tenacity 8.2 | GEO fetching with exponential back-off retry |
| Testing | pytest 7.4, pytest-cov 4.1 | 36 unit tests, 84% validator coverage |
| Logging | structlog 24.1 | Structured JSON logs with correlation fields |
| Configuration | pydantic-settings 2.1 | Environment-based config, .env support |
| CI | GitHub Actions | Lint and test on every push to main |

---

## Project Structure

```text
ngs-results-explorer/
├── src/
│   ├── config.py        Pydantic settings - environment-based configuration
│   ├── fetcher.py       NCBI GEO data fetcher with retry logic
│   ├── validator.py     9-rule validation engine with quarantine writer
│   ├── processor.py     DE results processor, Parquet output, summary JSON
│   ├── api.py           FastAPI REST API - 8 endpoints
│   └── dashboard.py     Plotly Dash interactive dashboard
├── tests/
│   ├── test_validator.py    27 validator unit tests
│   └── test_processor.py    9 processor unit tests
├── data/
│   ├── raw/             Cached GEO JSON - immutable once written
│   ├── processed/       Parquet dataset and summary JSON
│   └── quarantine/      Failed records with rejection reasons
├── .github/workflows/
│   └── ci.yml           GitHub Actions: flake8 + pytest on every push
├── requirements.txt
├── pyproject.toml
└── README.md
```

---

## Why This Matters for NGS Platforms

Bioinformatics Platforms solve the hard problem of running the analysis. What sits downstream - validating the output, structuring it for query, making it accessible to scientists who cannot write Python - is where data engineering earns its place.

Every design decision in this project maps directly to that space.

The nine validation rules exist because NCBI GEO submissions are contributed by thousands of independent research groups with inconsistent conventions. Raw internal UIDs leak through as gene IDs. Gene symbols go missing. P-values occasionally exceed 1.0 due to floating-point artefacts in DESeq2 edge cases. A platform that ingests these results without validation is importing noise alongside signal.

The quarantine-not-delete pattern exists because in a regulated biotech environment - any lab working under GxP guidelines - you cannot discard data without a documented reason. Every rejected record in this pipeline is preserved with a human-readable rejection reason and a UTC timestamp. An auditor can inspect the quarantine zone and trace exactly why each record was excluded.

The volcano plot exists because it is the first thing every RNA-Seq researcher opens when they get results back. It is the standard. Building it from scratch on real DE data demonstrates that the domain knowledge is there, not just the engineering.

The Parquet output exists because columnar storage is what makes downstream analytical queries fast at scale. A platform processing thousands of samples per day cannot afford to scan CSV files. Parquet with Snappy compression is the format every serious analytical pipeline uses.

The FastAPI layer exists because bioinformaticians work in Python notebooks. They call REST APIs. They do not log into dashboards to run queries. A clean, self-documenting API with pagination and filter parameters is what makes a data platform actually usable in a research environment.

---

*Data: NCBI GEO GSE183947 · Python · FastAPI · Plotly Dash · pandas · GitHub Actions*

