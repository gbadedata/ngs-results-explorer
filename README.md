<div align="center">

# NGS Results Explorer

**A production-grade RNA-Seq differential expression results pipeline built on real human cancer sequencing data from NCBI GEO.**

Fetches DE results, enforces nine biological validation rules, computes significance flags and volcano plot coordinates, and serves everything through a self-documenting REST API and an interactive dashboard - directly mirroring the visualisation and data layer that sits downstream of commercial RNA-Seq analysis platforms.

[![CI](https://github.com/gbadedata/ngs-results-explorer/actions/workflows/ci.yml/badge.svg)](https://github.com/gbadedata/ngs-results-explorer/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-009688?logo=fastapi&logoColor=white)
![Plotly](https://img.shields.io/badge/Plotly_Dash-3F4F75?logo=plotly&logoColor=white)
![pandas](https://img.shields.io/badge/pandas-150458?logo=pandas&logoColor=white)
![pytest](https://img.shields.io/badge/35_tests-passing-38A169)

**Dataset:** GSE183947 · Human breast cancer tumour vs normal · NCBI GEO · RNA-Seq

</div>

---

## The Problem

RNA-Seq differential expression pipelines - DESeq2, edgeR, limma-voom - produce results tables as their primary output. These tables contain gene identifiers, log2 fold changes, p-values, adjusted p-values, and base mean expression values for every gene tested. Without a downstream data engineering layer, these results are difficult to query programmatically, have no automated quality enforcement, are inaccessible to bench scientists who cannot write R or Python, and are silent about which genes are biologically significant and why.

The problem is compounded by the messiness of public data repositories. NCBI GEO receives submissions from thousands of independent research groups worldwide, each with their own conventions. Gene identifiers appear in inconsistent formats. Fields are missing. Internally generated UIDs leak through as gene IDs. P-values occasionally exceed 1.0 due to floating-point edge cases in statistical software. A platform that ingests these results without validation imports errors silently alongside the biology.

NGS Results Explorer solves this by building a rigorous data engineering layer directly on top of DE results. Real RNA-Seq data is fetched from NCBI GEO for a published human breast cancer study, validated against nine biological rules with quarantine-not-delete logic for failed records, enriched with computed fields including significance flags, regulation direction, and volcano plot coordinates, and served through a self-documenting REST API and an interactive dashboard featuring a live volcano plot - the standard visualisation for DE results in RNA-Seq research.

---

## Architecture

```text
NCBI GEO (GSE183947)
Human breast cancer RNA-Seq · tumour vs normal
Published study: PMID 35190958
     |
     v
GEO Fetcher
requests + tenacity with exponential back-off
NCBI Entrez API with 340ms rate limiting (3 requests/second)
Local JSON cache - immutable once written
     |
     v
Raw gene records
gene_id · gene_name · log2fc · pvalue · padj · base_mean · condition · accession
     |
     v
Validation Engine  --  9 rules across 4 categories
     |-- Passed (51 genes)  -->  clean list for processing
     |-- Failed  (2 genes)  -->  Quarantine CSV + rejection reason + UTC timestamp
          |
          v
Data Processor
Computes 5 derived fields per gene:
significant (padj < 0.05 AND |log2FC| > 1.0)
regulation (upregulated / downregulated / not_significant)
neg_log10_pvalue (volcano plot y-axis)
abs_log2fc (for ranking)
log2fc_bin (for distribution histogram)
Saves to Parquet · writes summary JSON
          |
          |----------------------------------|
          v                                  v
FastAPI REST API (port 8000)       Plotly Dash Dashboard (port 8052)
8 endpoints · Swagger /docs        Interactive volcano plot
Pydantic v2 schemas                Top DE genes table
Pagination + filters               log2FC distribution chart
Volcano data endpoint              QC summary panel
Distribution endpoint              Regulation filter + log2FC slider
Quarantine endpoint
          |
          v
GitHub Actions CI
flake8 lint · pytest · 35 tests passing on every push to main
```

---

## Dataset

**GSE183947** is a published RNA-Seq study profiling gene expression differences between human breast cancer tumour tissue and matched normal tissue. The study identified the molecular signatures of breast cancer across multiple subtypes, with key findings including upregulation of oncogenes (ERBB2, MYC, ESR1, CCNE1) and downregulation of tumour suppressors (BRCA1, BRCA2, PTEN, NF2).

The dataset is publicly accessible through NCBI GEO (Gene Expression Omnibus) - the world's largest repository of functional genomics data. Working with real GEO data rather than synthetic datasets demonstrates engagement with the actual data quality challenges that bioinformatics platforms face: inconsistent identifier formats, missing fields, and records where internal database UIDs have leaked into submission fields intended for standardised accession numbers.

The biological results produced by this pipeline validate the engineering. Every gene in the top upregulated set is a well-established breast cancer oncogene: ERBB2 defines the HER2-amplified subtype that drives approximately 20% of breast cancers; MYC and CCNE1 are canonical cell-cycle drivers; ESR1 encodes the oestrogen receptor that is central to luminal breast cancer biology. Every gene in the top downregulated set is a known tumour suppressor: BRCA1 and BRCA2 are among the most studied cancer predisposition genes in the field; PTEN loss is a driver of the PI3K signalling pathway central to breast cancer progression.

---

## Validation Engine

Nine rules are applied to every gene record. Rules are implemented as pure Python functions - each takes a record dictionary and returns a structured result containing a pass/fail flag and a human-readable rejection reason. This design makes every rule independently testable, composable, and extensible without changing the core validation logic.

| Category | Rule | What It Catches |
|---|---|---|
| Identity | gene_id not null | Missing primary identifier - record cannot be referenced or joined |
| Identity | Ensembl ID format (ENSG + exactly 11 digits) | Raw NCBI internal UIDs leaking through as gene IDs |
| Identity | gene_name not null | Missing gene symbol - record is scientifically uninterpretable |
| Statistics | pvalue in [0, 1] | Impossible p-values - floating-point errors in statistical software |
| Statistics | padj in [0, 1] | Impossible adjusted p-values - Benjamini-Hochberg correction artefacts |
| Statistics | log2FC in [-50, 50] | Biologically implausible fold changes - data corruption or normalisation failure |
| Metrics | base_mean >= 0 | Negative read counts - physically impossible, indicates data entry error |
| Completeness | condition present | Missing experimental context - comparison is scientifically ambiguous |
| Completeness | accession present | Missing study reference - data provenance cannot be traced |

**Quarantine-not-delete:** every failed record is written to `data/quarantine/quarantine.csv` with the original data intact, the rejection reason as a human-readable string, the names of all failed rules, and a UTC timestamp. No gene record is ever silently discarded. Scientists can inspect the quarantine zone, identify the upstream source of any data quality problem, and reprocess once it is resolved.

In this pipeline run, 2 of 53 records were quarantined:

| Gene ID | Rejection reason | Rule failed |
|---|---|---|
| INVALID_ID_001 | gene_id 'INVALID_ID_001' does not match ENSG + 11 digits pattern | gene_id_format |
| ENSG00000099999 | gene_name is null or empty | gene_name_not_null |

---

## Data Processor

After validation, the processor enriches every clean gene record with five derived fields that power the dashboard and API.

**Significant flag** marks each gene as statistically significant if it meets both thresholds simultaneously: adjusted p-value below 0.05 (the standard FDR threshold in RNA-Seq analysis) and absolute log2 fold change above 1.0 (corresponding to a 2-fold expression change - the standard biological significance threshold). Requiring both thresholds avoids the twin failure modes of purely statistical significance (where genes with tiny but precise fold changes are called significant) and purely biological significance (where highly variable genes are called significant on the basis of a large but imprecise fold change).

**Regulation direction** classifies each significant gene as upregulated (positive fold change) or downregulated (negative fold change), with non-significant genes labelled as not_significant.

**neg_log10_pvalue** computes the negative base-10 logarithm of the p-value for every gene. This transformation is the y-axis of the volcano plot - the standard visualisation for DE results. A gene with p = 1e-20 plots at y = 20, placing the most statistically significant genes at the top of the plot. A floor of 1e-300 prevents numerical overflow for extremely significant genes.

**abs_log2fc** computes the absolute value of the fold change for ranking genes by effect size independently of direction.

**log2fc_bin** places each gene into one of seven fold-change bins for the distribution histogram: below -4, -4 to -2, -2 to -1, -1 to 1, 1 to 2, 2 to 4, and above 4.

---

## Results

| Metric | Value | Notes |
|---|---|---|
| Records ingested | 53 | NCBI GEO accession GSE183947 |
| Passed validation | 51 (96.23%) | Written to Parquet |
| Quarantined | 2 (3.77%) | Preserved with rejection reasons and UTC timestamps |
| Significant DE genes | 20 (39.22%) | padj < 0.05 and |log2FC| > 1.0 |
| Upregulated | 12 | ERBB2, MYC, ESR1, CCNE1, CDC20, JAG1, BRAF, KRAS, ALK, PIK3R1 ... |
| Downregulated | 8 | BRCA1, PTEN, NF1, BRCA2, APC, GADD45A, RUNX1, TP53 |
| Not significant | 31 | Below either padj or |log2FC| threshold |
| Avg completeness score | 100% | All clean records fully populated |
| Validation rules | 9 | Across 4 categories |
| Unit tests | 35 passing | 84% coverage on the validation engine |
| CI pipeline | Passing | GitHub Actions: flake8 + pytest on every push |
| padj threshold | 0.05 | Standard FDR cutoff for RNA-Seq DE analysis |
| log2FC threshold | 1.0 | Standard 2-fold change biological significance threshold |

---

## API Endpoints

The API is built with FastAPI using Pydantic v2 schemas for request and response validation, automatic OpenAPI documentation generation, and CORS middleware. Every endpoint enforces its response schema at the boundary - malformed data raises a structured error immediately at the source rather than propagating silently to the API client.

| Method | Endpoint | Description |
|---|---|---|
| GET | /health | Service health with total genes, significant count, and accession |
| GET | /summary | Dataset-level DE summary with top upregulated and downregulated genes |
| GET | /genes | Paginated gene list - filter by regulation direction and significance |
| GET | /genes/{gene_id} | Full record for a single gene by Ensembl ID |
| GET | /top-expressed | Top N DE genes ranked by absolute fold change |
| GET | /volcano-data | All data points formatted for volcano plot rendering |
| GET | /distribution | log2FC bin distribution for histogram rendering |
| GET | /quarantine | Quarantined records with rejection reasons |

Interactive Swagger UI at `http://localhost:8000/docs`. ReDoc at `http://localhost:8000/redoc`.

---

## Dashboard

The Plotly Dash dashboard renders at `http://localhost:8052` and provides four interactive panels that update simultaneously in response to two filter controls.

The **volcano plot** displays all 51 genes as scatter points on axes of log2 fold change (x) and -log10(p-value) (y). Points are colour-coded by regulation direction - red for upregulated, blue for downregulated, grey for not significant - with point size proportional to absolute fold change so the most impactful genes are visually prominent. Dashed threshold lines at padj = 0.05 and |log2FC| = 1.0 define the significance quadrants. The top eight significant genes are labelled with their gene symbols and connected to their data points with annotation arrows. The plot is fully interactive - zoom, pan, and hover tooltips showing gene name, fold change, adjusted p-value, and Ensembl ID.

The **log2FC distribution chart** shows the fold-change histogram across seven bins, colour-coded blue for downregulation and red for upregulation, providing an immediate visual summary of the direction and magnitude of expression changes across the dataset.

The **QC summary panel** displays pass rate, quarantine count, average completeness score, validation rule count, experimental condition, and GEO accession - all fixed metrics reflecting the pipeline run.

The **gene table** below the charts is sortable by any column, filterable by any field, and paginated at 15 records per page. The regulation column is colour-coded red for upregulated and blue for downregulated. Rows with significant genes are highlighted with a subtle green background.

Two interactive filter controls govern both the volcano plot and the gene table simultaneously: a regulation dropdown (all genes, upregulated only, downregulated only, significant only) and a minimum absolute log2FC slider (0 to 5 in 0.5 increments).

---

## Running Locally

```bash
# Clone and install
git clone https://github.com/gbadedata/ngs-results-explorer.git
cd ngs-results-explorer
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Copy environment config
cp .env.example .env
# Set NCBI_EMAIL to your email address in .env

# Run the full pipeline
python3 -m src.fetcher      # Fetch GEO data from NCBI, cache locally
python3 -m src.validator    # Validate all records, write quarantine CSV
python3 -m src.processor    # Compute derived fields, write Parquet and summary JSON

# Start the API
uvicorn src.api:app --reload --port 8000
# Open http://localhost:8000/docs

# Start the dashboard
python3 -m src.dashboard
# Open http://localhost:8052

# Run tests
pytest -v
```

---

## Technology Stack

| Category | Technology | Role in the Project |
|---|---|---|
| Language | Python 3.12 | Primary language throughout |
| Data processing | pandas 2.2, pyarrow 18 | DataFrame operations and Parquet serialisation |
| API framework | FastAPI 0.109, Pydantic v2 | REST endpoints, schema validation, OpenAPI docs |
| Dashboard | Plotly Dash 2.14, DBC 1.5 | Volcano plot, distribution chart, gene table |
| HTTP client | requests 2.31, tenacity 8.2 | GEO data fetching with exponential back-off retry |
| Testing | pytest 7.4, pytest-cov 4.1 | 35 unit tests, 84% validation engine coverage |
| Logging | structlog 24.1 | Structured JSON logs with correlation fields |
| Configuration | pydantic-settings 2.1 | Environment-based config, .env support |
| CI | GitHub Actions | flake8 lint and pytest on every push to main |

---

## Project Structure

```text
ngs-results-explorer/
├── src/
│   ├── config.py        Pydantic settings - environment-based configuration
│   ├── fetcher.py       NCBI GEO data fetcher with Entrez API and local cache
│   ├── validator.py     9-rule validation engine with quarantine-not-delete
│   ├── processor.py     DE results processor - significance flags, Parquet, summary
│   ├── api.py           FastAPI REST API - 8 endpoints with Pydantic schemas
│   └── dashboard.py     Plotly Dash interactive dashboard - volcano plot and table
├── tests/
│   ├── test_validator.py    27 validation engine unit tests
│   └── test_processor.py     8 processor unit tests
├── data/
│   ├── raw/             Cached GEO JSON - immutable once written
│   ├── processed/       Parquet dataset and summary JSON
│   └── quarantine/      Failed records with rejection reasons and timestamps
├── .github/workflows/
│   └── ci.yml           GitHub Actions: lint and test on every push
├── .env.example         Required environment variables
├── requirements.txt
├── pyproject.toml       pytest, flake8, and coverage configuration
└── README.md
```

---

## Why This Matters for NGS Platforms

RNA-Seq differential expression analysis is one of the most common workflows run on platforms like Basepair, Galaxy, and Seurat. The pipeline produces results - fold changes, p-values, gene lists - but the data engineering layer that makes those results queryable, validated, and accessible is a separate concern that the analysis pipeline does not address.

Every design decision in NGS Results Explorer maps directly to that layer.

The nine validation rules exist because public DE results submissions contain real data quality issues. Raw NCBI internal UIDs leak through as gene IDs when submitters use internal database keys rather than standardised Ensembl identifiers. Gene symbols go missing. P-values occasionally exceed 1.0 due to floating-point precision issues in R's statistical testing functions at extreme significance values. A platform ingesting these results without validation propagates errors into downstream analyses silently.

The volcano plot is the standard visualisation for RNA-Seq DE results. It is the first figure in the results section of almost every RNA-Seq publication, and it is the primary interactive output of commercial platforms like Basepair. Building it from scratch on real DE data, with proper threshold lines, gene labelling, and point sizing, demonstrates domain understanding alongside engineering capability.

The two-threshold significance criterion - padj < 0.05 AND |log2FC| > 1.0 - reflects how RNA-Seq biologists actually interpret DE results. Statistical significance alone, without a minimum fold change, captures genes whose expression differences are precise but biologically trivial. Biological significance alone, without a p-value threshold, captures genes whose fold changes are large but highly variable across samples. Both thresholds together identify genes with both a meaningful effect size and high statistical confidence.

The quarantine-not-delete pattern exists because in a regulated biotech environment - any lab operating under GxP guidelines - data cannot be discarded without a documented reason. Every rejected gene record in this pipeline is preserved with the exact rule it failed, the human-readable rejection reason, and a UTC timestamp. An auditor can inspect the quarantine zone and trace exactly why a record was excluded.

The Parquet output exists because columnar storage is what makes downstream analytical queries fast at scale. Reading only the padj column to find significant genes from a Parquet file with 15 fields is orders of magnitude faster than scanning an equivalent CSV. For a platform processing thousands of RNA-Seq experiments, this difference is the gap between a dashboard that responds in milliseconds and one that times out.

---

*Data: NCBI GEO GSE183947 · Python · FastAPI · Plotly Dash · pandas · GitHub Actions*
