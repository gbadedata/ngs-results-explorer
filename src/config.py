"""Configuration management using pydantic-settings."""
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    ncbi_email: str = "gbadedata@gmail.com"
    ncbi_api_key: str = ""
    app_version: str = "1.0.0"
    log_level: str = "INFO"

    # GEO dataset to fetch
    geo_accession: str = "GSE183947"

    # Data paths
    raw_data_path: str = "data/raw"
    processed_data_path: str = "data/processed"
    quarantine_path: str = "data/quarantine"

    # Validation thresholds
    max_pvalue: float = 1.0
    min_pvalue: float = 0.0
    max_log2fc: float = 50.0
    min_log2fc: float = -50.0
    min_base_mean: float = 0.0

    class Config:
        env_file = ".env"


settings = Settings()
