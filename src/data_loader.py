"""
src/data_loader.py
------------------
Loads the curated ClinVar + GWAS Catalog dataset and exposes
helper functions used by both the Streamlit app and the Gradio demo.

Data sources
------------
- ClinVar  : clinvar.ncbi.nlm.nih.gov  (rsID, clinical significance)
- GWAS Cat : ebi.ac.uk/gwas            (p-value, odds ratio, allele freq)
- PubMed   : linked via pubmed_id field
"""

from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

# ── colour palette shared across all plots ──────────────────────────────────
CONDITION_COLORS = {
    "Type 2 Diabetes":     "#3B82F6",   # blue
    "Coronary Artery Disease": "#EF4444",  # red
}

SIG_COLORS = {
    "Pathogenic":          "#DC2626",
    "Likely Pathogenic":   "#F97316",
    "Risk Factor":         "#3B82F6",
    "Uncertain Significance": "#9CA3AF",
    "Benign":              "#22C55E",
}

SIG_ORDER = ["Pathogenic", "Likely Pathogenic", "Risk Factor",
             "Uncertain Significance", "Benign"]


# ── loader ───────────────────────────────────────────────────────────────────

_DATA_PATH = Path(__file__).parent.parent / "data" / "clinvar_gwas_variants.csv"


def load_data() -> pd.DataFrame:
    """Load and clean the curated variant dataset."""
    df = pd.read_csv(_DATA_PATH)

    # coerce numerics
    df["allele_freq"] = pd.to_numeric(df["allele_freq"], errors="coerce")
    df["odds_ratio"]  = pd.to_numeric(df["odds_ratio"],  errors="coerce")
    df["p_value"]     = pd.to_numeric(df["p_value"],     errors="coerce")
    df["chromosome"]  = df["chromosome"].astype(str)

    # derived columns
    df["neg_log10_p"]    = -np.log10(df["p_value"].clip(lower=1e-60))
    df["risk_direction"] = df["odds_ratio"].apply(
        lambda x: "Protective (OR < 1)" if x < 1 else "Risk-increasing (OR ≥ 1)"
    )
    df["pubmed_url"] = df["pubmed_id"].apply(
        lambda i: f"https://pubmed.ncbi.nlm.nih.gov/{int(i)}/" if pd.notna(i) else ""
    )

    return df


def filter_data(
    df: pd.DataFrame,
    conditions: list[str] | None = None,
    sig_levels: list[str] | None = None,
    min_af: float = 0.0,
    max_af: float = 1.0,
) -> pd.DataFrame:
    """Apply sidebar filters and return filtered copy."""
    out = df.copy()
    if conditions:
        out = out[out["condition"].isin(conditions)]
    if sig_levels:
        out = out[out["clinical_significance"].isin(sig_levels)]
    out = out[(out["allele_freq"] >= min_af) & (out["allele_freq"] <= max_af)]
    return out


def summary_stats(df: pd.DataFrame) -> dict:
    return {
        "total_variants": len(df),
        "unique_genes":   df["gene"].nunique(),
        "conditions":     df["condition"].nunique(),
        "pathogenic":     len(df[df["clinical_significance"] == "Pathogenic"]),
        "likely_path":    len(df[df["clinical_significance"] == "Likely Pathogenic"]),
        "risk_factor":    len(df[df["clinical_significance"] == "Risk Factor"]),
    }
