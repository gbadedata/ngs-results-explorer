"""Unit tests for the validation engine."""
import pytest
from src.validator import (
    validate_record,
    validate_batch,
    completeness_score,
    rule_gene_id_format,
    rule_gene_id_not_null,
    rule_gene_name_not_null,
    rule_pvalue_range,
    rule_padj_range,
    rule_log2fc_range,
    rule_base_mean_non_negative,
    rule_condition_present,
    rule_accession_present,
)


# ── Fixtures ─────────────────────────────────────────────────────

@pytest.fixture
def valid_record():
    return {
        "gene_id": "ENSG00000141736",
        "gene_name": "ERBB2",
        "log2fc": 4.21,
        "pvalue": 1.2e-15,
        "padj": 3.4e-13,
        "base_mean": 892.3,
        "condition": "tumour_vs_normal",
        "accession": "GSE183947",
        "record_index": 0,
    }


@pytest.fixture
def invalid_id_record(valid_record):
    r = valid_record.copy()
    r["gene_id"] = "INVALID_ID_001"
    return r


@pytest.fixture
def missing_name_record(valid_record):
    r = valid_record.copy()
    r["gene_name"] = ""
    return r


@pytest.fixture
def bad_pvalue_record(valid_record):
    r = valid_record.copy()
    r["pvalue"] = 1.5  # > 1.0 — impossible
    return r


@pytest.fixture
def negative_base_mean_record(valid_record):
    r = valid_record.copy()
    r["base_mean"] = -10.0
    return r


@pytest.fixture
def extreme_log2fc_record(valid_record):
    r = valid_record.copy()
    r["log2fc"] = 999.0  # biologically implausible
    return r


# ── Individual rule tests ─────────────────────────────────────────

class TestGeneIdFormat:
    def test_valid_ensembl_id(self, valid_record):
        result = rule_gene_id_format(valid_record)
        assert result.passed is True

    def test_invalid_format(self, invalid_id_record):
        result = rule_gene_id_format(invalid_id_record)
        assert result.passed is False
        assert "ENSG" in result.reason

    def test_too_short(self, valid_record):
        valid_record["gene_id"] = "ENSG0000014"
        result = rule_gene_id_format(valid_record)
        assert result.passed is False

    def test_wrong_prefix(self, valid_record):
        valid_record["gene_id"] = "GENE00000141736"
        result = rule_gene_id_format(valid_record)
        assert result.passed is False


class TestGeneIdNotNull:
    def test_valid(self, valid_record):
        assert rule_gene_id_not_null(valid_record).passed is True

    def test_empty_string(self, valid_record):
        valid_record["gene_id"] = ""
        assert rule_gene_id_not_null(valid_record).passed is False

    def test_whitespace(self, valid_record):
        valid_record["gene_id"] = "   "
        assert rule_gene_id_not_null(valid_record).passed is False

    def test_missing_key(self):
        assert rule_gene_id_not_null({}).passed is False


class TestGeneNameNotNull:
    def test_valid(self, valid_record):
        assert rule_gene_name_not_null(valid_record).passed is True

    def test_empty(self, missing_name_record):
        assert rule_gene_name_not_null(missing_name_record).passed is False


class TestPvalueRange:
    def test_valid_pvalue(self, valid_record):
        assert rule_pvalue_range(valid_record).passed is True

    def test_pvalue_above_one(self, bad_pvalue_record):
        assert rule_pvalue_range(bad_pvalue_record).passed is False

    def test_pvalue_zero(self, valid_record):
        valid_record["pvalue"] = 0.0
        assert rule_pvalue_range(valid_record).passed is True

    def test_pvalue_exactly_one(self, valid_record):
        valid_record["pvalue"] = 1.0
        assert rule_pvalue_range(valid_record).passed is True

    def test_pvalue_none(self, valid_record):
        valid_record["pvalue"] = None
        assert rule_pvalue_range(valid_record).passed is False


class TestPadjRange:
    def test_valid(self, valid_record):
        assert rule_padj_range(valid_record).passed is True

    def test_above_one(self, valid_record):
        valid_record["padj"] = 1.1
        assert rule_padj_range(valid_record).passed is False


class TestLog2fcRange:
    def test_valid(self, valid_record):
        assert rule_log2fc_range(valid_record).passed is True

    def test_extreme_value(self, extreme_log2fc_record):
        assert rule_log2fc_range(extreme_log2fc_record).passed is False

    def test_negative_extreme(self, valid_record):
        valid_record["log2fc"] = -999.0
        assert rule_log2fc_range(valid_record).passed is False

    def test_zero(self, valid_record):
        valid_record["log2fc"] = 0.0
        assert rule_log2fc_range(valid_record).passed is True


class TestBaseMeanNonNegative:
    def test_valid(self, valid_record):
        assert rule_base_mean_non_negative(valid_record).passed is True

    def test_negative(self, negative_base_mean_record):
        assert rule_base_mean_non_negative(negative_base_mean_record).passed is False

    def test_zero(self, valid_record):
        valid_record["base_mean"] = 0.0
        assert rule_base_mean_non_negative(valid_record).passed is True


class TestConditionPresent:
    def test_valid(self, valid_record):
        assert rule_condition_present(valid_record).passed is True

    def test_empty(self, valid_record):
        valid_record["condition"] = ""
        assert rule_condition_present(valid_record).passed is False


class TestAccessionPresent:
    def test_valid(self, valid_record):
        assert rule_accession_present(valid_record).passed is True

    def test_missing(self, valid_record):
        valid_record["accession"] = ""
        assert rule_accession_present(valid_record).passed is False


# ── Batch validation tests ────────────────────────────────────────

class TestValidateBatch:
    def test_all_valid(self, valid_record):
        results = validate_batch([valid_record])
        assert results["stats"]["passed"] == 1
        assert results["stats"]["quarantined"] == 0

    def test_one_invalid(self, valid_record, invalid_id_record):
        results = validate_batch([valid_record, invalid_id_record])
        assert results["stats"]["passed"] == 1
        assert results["stats"]["quarantined"] == 1

    def test_quarantine_has_reason(self, invalid_id_record):
        results = validate_batch([invalid_id_record])
        quarantined = results["quarantined"]
        assert len(quarantined) == 1
        assert "rejection_reasons" in quarantined[0]
        assert quarantined[0]["rejection_reasons"] != ""

    def test_pass_rate_calculation(self, valid_record, invalid_id_record):
        results = validate_batch([valid_record, valid_record, invalid_id_record])
        assert results["stats"]["pass_rate"] == pytest.approx(66.67, abs=0.1)

    def test_empty_batch(self):
        results = validate_batch([])
        assert results["stats"]["total"] == 0
        assert results["stats"]["passed"] == 0


# ── Completeness score tests ──────────────────────────────────────

class TestCompletenessScore:
    def test_fully_complete(self, valid_record):
        score = completeness_score(valid_record)
        assert score == 1.0

    def test_missing_one_field(self, valid_record):
        valid_record["gene_name"] = ""
        score = completeness_score(valid_record)
        assert score < 1.0
        assert score > 0.0

    def test_empty_record(self):
        score = completeness_score({})
        assert score == 0.0
