"""
Validation Engine for NGS differential expression results.
Applies 9 rules to each gene record. Failed records are
quarantined with rejection reason — never silently dropped.
"""
import os
import re
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
import structlog

from src.config import settings

log = structlog.get_logger()

# Ensembl gene ID pattern: ENSG followed by exactly 11 digits
ENSEMBL_PATTERN = re.compile(r"^ENSG\d{11}$")


@dataclass
class ValidationResult:
    passed: bool
    rule: str
    reason: str


@dataclass
class RecordOutcome:
    record: dict
    passed: bool
    failures: list[ValidationResult] = field(default_factory=list)


# ─────────────────────────────────────────────
# The 9 validation rules — pure functions
# each returns a ValidationResult
# ─────────────────────────────────────────────

def rule_gene_id_not_null(record: dict) -> ValidationResult:
    passed = bool(record.get("gene_id", "").strip())
    return ValidationResult(
        passed=passed,
        rule="gene_id_not_null",
        reason="" if passed else "gene_id is null or empty",
    )


def rule_gene_id_format(record: dict) -> ValidationResult:
    gene_id = record.get("gene_id", "")
    passed = bool(ENSEMBL_PATTERN.match(gene_id))
    return ValidationResult(
        passed=passed,
        rule="gene_id_format",
        reason="" if passed else (
            f"gene_id '{gene_id}' does not match ENSG + 11 digits pattern"
        ),
    )


def rule_gene_name_not_null(record: dict) -> ValidationResult:
    passed = bool(record.get("gene_name", "").strip())
    return ValidationResult(
        passed=passed,
        rule="gene_name_not_null",
        reason="" if passed else "gene_name is null or empty",
    )


def rule_pvalue_range(record: dict) -> ValidationResult:
    pvalue = record.get("pvalue")
    try:
        pvalue = float(pvalue)
        passed = settings.min_pvalue <= pvalue <= settings.max_pvalue
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="pvalue_range",
        reason="" if passed else (
            f"pvalue '{pvalue}' not in valid range [0, 1]"
        ),
    )


def rule_padj_range(record: dict) -> ValidationResult:
    padj = record.get("padj")
    try:
        padj = float(padj)
        passed = 0.0 <= padj <= 1.0
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="padj_range",
        reason="" if passed else (
            f"padj '{padj}' not in valid range [0, 1]"
        ),
    )


def rule_log2fc_range(record: dict) -> ValidationResult:
    log2fc = record.get("log2fc")
    try:
        log2fc = float(log2fc)
        passed = settings.min_log2fc <= log2fc <= settings.max_log2fc
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="log2fc_range",
        reason="" if passed else (
            f"log2fc '{log2fc}' outside plausible range "
            f"[{settings.min_log2fc}, {settings.max_log2fc}]"
        ),
    )


def rule_base_mean_non_negative(record: dict) -> ValidationResult:
    base_mean = record.get("base_mean")
    try:
        base_mean = float(base_mean)
        passed = base_mean >= settings.min_base_mean
    except (TypeError, ValueError):
        passed = False
    return ValidationResult(
        passed=passed,
        rule="base_mean_non_negative",
        reason="" if passed else (
            f"base_mean '{base_mean}' is negative — invalid count"
        ),
    )


def rule_condition_present(record: dict) -> ValidationResult:
    passed = bool(record.get("condition", "").strip())
    return ValidationResult(
        passed=passed,
        rule="condition_present",
        reason="" if passed else "condition field is missing or empty",
    )


def rule_accession_present(record: dict) -> ValidationResult:
    passed = bool(record.get("accession", "").strip())
    return ValidationResult(
        passed=passed,
        rule="accession_present",
        reason="" if passed else "accession field is missing or empty",
    )


# Ordered list of all 9 rules
ALL_RULES = [
    rule_gene_id_not_null,
    rule_gene_id_format,
    rule_gene_name_not_null,
    rule_pvalue_range,
    rule_padj_range,
    rule_log2fc_range,
    rule_base_mean_non_negative,
    rule_condition_present,
    rule_accession_present,
]


# ─────────────────────────────────────────────
# Completeness scoring
# ─────────────────────────────────────────────

COMPLETENESS_FIELDS = [
    "gene_id", "gene_name", "log2fc", "pvalue",
    "padj", "base_mean", "condition", "accession",
]


def completeness_score(record: dict) -> float:
    """Return fraction of key fields that are populated."""
    populated = sum(
        1 for f in COMPLETENESS_FIELDS
        if record.get(f) is not None and str(record.get(f, "")).strip() != ""
    )
    return round(populated / len(COMPLETENESS_FIELDS), 4)


# ─────────────────────────────────────────────
# Core validation function
# ─────────────────────────────────────────────

def validate_record(record: dict) -> RecordOutcome:
    """Apply all 9 rules to a single record."""
    failures = []
    for rule_fn in ALL_RULES:
        result = rule_fn(record)
        if not result.passed:
            failures.append(result)

    return RecordOutcome(
        record=record,
        passed=len(failures) == 0,
        failures=failures,
    )


# ─────────────────────────────────────────────
# Batch validation — returns clean + quarantine
# ─────────────────────────────────────────────

def validate_batch(records: list[dict]) -> dict:
    """
    Validate a batch of DE records.
    Returns clean records (passed) and quarantined records (failed).
    Writes quarantine CSV to disk with rejection reasons.
    """
    clean = []
    quarantined = []
    total = len(records)

    for record in records:
        outcome = validate_record(record)
        score = completeness_score(record)
        record["completeness_score"] = score

        if outcome.passed:
            clean.append(record)
        else:
            rejection_reasons = "; ".join(
                f.reason for f in outcome.failures
            )
            quarantine_record = {
                **record,
                "rejection_reasons": rejection_reasons,
                "failed_rules": [f.rule for f in outcome.failures],
                "quarantined_at": datetime.now(tz=timezone.utc).isoformat(),
            }
            quarantined.append(quarantine_record)

    # Write quarantine to disk
    _write_quarantine(quarantined)

    # Compute stats
    pass_rate = round(len(clean) / total * 100, 2) if total > 0 else 0
    quarantine_rate = round(len(quarantined) / total * 100, 2) if total > 0 else 0
    avg_completeness = round(
        sum(r["completeness_score"] for r in clean) / len(clean), 4
    ) if clean else 0

    log.info(
        "validation_complete",
        total=total,
        passed=len(clean),
        quarantined=len(quarantined),
        pass_rate=f"{pass_rate}%",
        quarantine_rate=f"{quarantine_rate}%",
        avg_completeness=avg_completeness,
    )

    return {
        "clean": clean,
        "quarantined": quarantined,
        "stats": {
            "total": total,
            "passed": len(clean),
            "quarantined": len(quarantined),
            "pass_rate": pass_rate,
            "quarantine_rate": quarantine_rate,
            "avg_completeness": avg_completeness,
        },
    }


def _write_quarantine(quarantined: list[dict]) -> None:
    """Write quarantined records to CSV with rejection reasons."""
    if not quarantined:
        return

    os.makedirs(settings.quarantine_path, exist_ok=True)
    path = os.path.join(
        settings.quarantine_path, "quarantine.csv"
    )

    fieldnames = [
        "gene_id", "gene_name", "log2fc", "pvalue", "padj",
        "base_mean", "condition", "accession", "completeness_score",
        "rejection_reasons", "quarantined_at",
    ]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(quarantined)

    log.info("quarantine_written", path=path, count=len(quarantined))


if __name__ == "__main__":
    import structlog
    structlog.configure()

    from src.fetcher import fetch_de_results
    records = fetch_de_results("GSE183947")
    results = validate_batch(records)

    print("\n=== Validation Summary ===")
    print(f"Total records : {results['stats']['total']}")
    print(f"Passed        : {results['stats']['passed']}")
    print(f"Quarantined   : {results['stats']['quarantined']}")
    print(f"Pass rate     : {results['stats']['pass_rate']}%")
    print(f"Avg completeness: {results['stats']['avg_completeness']}")
    print("\nQuarantined records:")
    for q in results["quarantined"]:
        print(f"  {q['gene_id']} | {q['rejection_reasons']}")
