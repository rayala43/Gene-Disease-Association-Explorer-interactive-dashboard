"""
app.py  —  Gene–Disease Association Explorer (Streamlit)
=========================================================
Interactive dashboard to explore gene–rsID–disease relationships
for Type 2 Diabetes and Coronary Artery Disease.

Data: ClinVar + GWAS Catalog (30 curated, literature-validated variants)

Run:
    streamlit run app.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import streamlit as st
import pandas as pd
from src.data_loader import load_data, filter_data, summary_stats, SIG_ORDER, CONDITION_COLORS
from src.charts import (
    build_network,
    build_bubble_scatter,
    build_heatmap,
    build_variants_per_gene,
    build_sig_distribution,
    build_chromosome_scatter,
)

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Gene–Disease Association Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Custom CSS ───────────────────────────────────────────────────────────────
st.markdown("""
<style>
    .main { background-color: #f8fafc; }
    .block-container { padding-top: 1.5rem; padding-bottom: 2rem; }
    .metric-card {
        background: white;
        border-radius: 12px;
        padding: 18px 20px;
        box-shadow: 0 1px 4px rgba(0,0,0,0.08);
        text-align: center;
    }
    .metric-card .num  { font-size: 32px; font-weight: 700; }
    .metric-card .lbl  { font-size: 12px; color: #64748b; text-transform: uppercase; letter-spacing: .5px; }
    .section-title { font-size: 18px; font-weight: 600; color: #1e293b; margin: 8px 0 4px; }
    .data-source { font-size: 11px; color: #94a3b8; margin-top: 4px; }
    div[data-testid="stMetric"] { background: white; border-radius: 10px;
                                   padding: 14px; box-shadow: 0 1px 3px rgba(0,0,0,0.07); }
</style>
""", unsafe_allow_html=True)


# ── Load data ─────────────────────────────────────────────────────────────────
@st.cache_data
def get_data():
    return load_data()

df_full = get_data()

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🧬 Explorer Filters")
    st.markdown("---")

    selected_conditions = st.multiselect(
        "Condition",
        options=df_full["condition"].unique().tolist(),
        default=df_full["condition"].unique().tolist(),
    )

    selected_sig = st.multiselect(
        "Clinical Significance",
        options=SIG_ORDER,
        default=SIG_ORDER,
    )

    af_range = st.slider(
        "Allele Frequency Range",
        min_value=0.0, max_value=1.0,
        value=(0.0, 1.0), step=0.01,
        format="%.2f",
    )

    st.markdown("---")
    st.markdown("### 🔍 Gene Search")
    gene_search = st.text_input("Search gene (e.g. TCF7L2)", placeholder="Gene symbol…")

    st.markdown("---")
    st.markdown(
        "**Data Sources**  \n"
        "🔬 [ClinVar](https://clinvar.ncbi.nlm.nih.gov)  \n"
        "📊 [GWAS Catalog](https://ebi.ac.uk/gwas)  \n"
        "📖 [PubMed](https://pubmed.ncbi.nlm.nih.gov)  \n\n"
        "*30 literature-validated variants*  \n"
        "*15 T2D · 15 CAD genes*"
    )

    st.markdown("---")
    st.markdown(
        "Built by **Rayala Madhu Bhanu Varma**  \n"
        "Lead Data Scientist · GenepoweRx"
    )

# ── Apply filters ─────────────────────────────────────────────────────────────
df = filter_data(
    df_full,
    conditions=selected_conditions if selected_conditions else None,
    sig_levels=selected_sig if selected_sig else None,
    min_af=af_range[0],
    max_af=af_range[1],
)

if gene_search:
    df = df[df["gene"].str.upper().str.contains(gene_search.upper())]

# ── Header ────────────────────────────────────────────────────────────────────
st.markdown("# 🧬 Gene–Disease Association Explorer")
st.markdown(
    "Interactive visualization of gene–rsID relationships for **Type 2 Diabetes** and "
    "**Coronary Artery Disease** · Data: ClinVar + GWAS Catalog · "
    f"Showing **{len(df)}** of {len(df_full)} variants"
)
st.markdown("---")

# ── KPI cards ────────────────────────────────────────────────────────────────
stats = summary_stats(df)
c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Total Variants",  stats["total_variants"])
c2.metric("Unique Genes",    stats["unique_genes"])
c3.metric("Conditions",      stats["conditions"])
c4.metric("Pathogenic",      stats["pathogenic"],   delta="High priority")
c5.metric("Likely Path.",    stats["likely_path"],  delta="Moderate")
c6.metric("Risk Factor",     stats["risk_factor"],  delta="GWAS hits")

st.markdown("---")

# ── Tab layout ───────────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "🕸️ Network Graph",
    "🫧 Bubble Chart",
    "🗺️ Heatmap",
    "📊 Bar Charts",
    "🧬 Chromosome Map",
])

with tab1:
    st.markdown('<p class="section-title">Gene–Disease Association Network</p>', unsafe_allow_html=True)
    st.caption("Nodes = genes (blue circles) and diseases (colored diamonds). Edge thickness = number of shared variants.")
    if df.empty:
        st.warning("No data matches the current filters.")
    else:
        st.plotly_chart(build_network(df), use_container_width=True)

with tab2:
    st.markdown('<p class="section-title">Odds Ratio vs Allele Frequency</p>', unsafe_allow_html=True)
    st.caption("Each bubble is a variant. Bubble size = -log₁₀(p-value). Hover for full details including mechanism.")
    if df.empty:
        st.warning("No data matches the current filters.")
    else:
        st.plotly_chart(build_bubble_scatter(df), use_container_width=True)

with tab3:
    st.markdown('<p class="section-title">Gene × Condition Risk Heatmap</p>', unsafe_allow_html=True)
    st.caption("Mean Odds Ratio per gene–condition pair. Deeper red = stronger disease association.")
    if df.empty:
        st.warning("No data matches the current filters.")
    else:
        st.plotly_chart(build_heatmap(df), use_container_width=True)

with tab4:
    st.markdown('<p class="section-title">Variant Counts & Significance Distribution</p>', unsafe_allow_html=True)
    col_a, col_b = st.columns(2)
    with col_a:
        st.caption("Variants per gene, coloured by condition")
        if not df.empty:
            st.plotly_chart(build_variants_per_gene(df), use_container_width=True)
    with col_b:
        st.caption("Clinical significance breakdown by condition")
        if not df.empty:
            st.plotly_chart(build_sig_distribution(df), use_container_width=True)

with tab5:
    st.markdown('<p class="section-title">Chromosomal Distribution</p>', unsafe_allow_html=True)
    st.caption("Manhattan-style plot. Dashed red = genome-wide significance (5×10⁻⁸). Hover for rsID and gene.")
    if df.empty:
        st.warning("No data matches the current filters.")
    else:
        st.plotly_chart(build_chromosome_scatter(df), use_container_width=True)

# ── Variant Table ─────────────────────────────────────────────────────────────
st.markdown("---")
st.markdown("### 📋 Variant Detail Table")

display_cols = ["rsid", "gene", "condition", "clinical_significance",
                "allele_freq", "odds_ratio", "p_value", "mechanism", "phenotype_detail"]

display_df = df[display_cols].copy()
display_df["allele_freq"]  = display_df["allele_freq"].map(lambda x: f"{x:.3f}")
display_df["odds_ratio"]   = display_df["odds_ratio"].map(lambda x: f"{x:.2f}")
display_df["p_value"]      = display_df["p_value"].map(lambda x: f"{x:.2e}")
display_df.columns         = ["rsID", "Gene", "Condition", "Clinical Sig.",
                               "Allele Freq.", "OR", "p-value", "Mechanism", "Notes"]

col_search, col_export = st.columns([4, 1])
with col_export:
    csv_bytes = df[display_cols].to_csv(index=False).encode()
    st.download_button(
        "⬇️ Download CSV",
        data=csv_bytes,
        file_name="gene_disease_variants.csv",
        mime="text/csv",
    )

st.dataframe(
    display_df,
    use_container_width=True,
    height=320,
    column_config={
        "rsID": st.column_config.TextColumn("rsID", width=100),
        "Gene": st.column_config.TextColumn("Gene", width=80),
        "OR":   st.column_config.TextColumn("OR", width=60),
    }
)

st.caption(
    "Sources: ClinVar (clinvar.ncbi.nlm.nih.gov) · GWAS Catalog (ebi.ac.uk/gwas) · "
    "All rsIDs and effect sizes from published peer-reviewed studies."
)
