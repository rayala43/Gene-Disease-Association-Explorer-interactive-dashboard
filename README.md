# 🧬 Gene–Disease Association Explorer

An interactive multi-view dashboard for exploring gene–rsID relationships in **Type 2 Diabetes** and **Coronary Artery Disease**, built with Streamlit (full app) and Gradio (HuggingFace-deployable demo).

**Data**: 30 curated, literature-validated variants from ClinVar + GWAS Catalog across 26 genes.

---

## Live Demo

| Interface | Command | Best for |
|---|---|---|
| Streamlit | `streamlit run app.py` | Full interactive dashboard |
| Gradio | `python gradio_demo.py` | Quick demo / HuggingFace Spaces |

---

## Visualizations

| Chart | What it shows |
|---|---|
| **🕸️ Network Graph** | Bipartite gene ↔ disease graph; edge thickness = # shared variants |
| **🫧 Bubble Chart** | Odds Ratio vs Allele Frequency; bubble size = −log₁₀(p-value) |
| **🗺️ Heatmap** | Mean OR per gene × condition; deeper red = stronger association |
| **📊 Bar Charts** | Variant counts per gene + clinical significance distribution |
| **🧬 Chromosome Map** | Manhattan-style plot with genome-wide significance thresholds |

All charts are interactive — hover for full variant details (rsID, gene, mechanism, p-value, OR).

---

## Dataset

`data/clinvar_gwas_variants.csv` — 30 variants, all from published literature:

| Field | Description |
|---|---|
| `rsid` | dbSNP rsID |
| `gene` | HGNC gene symbol |
| `chromosome / position` | GRCh38 coordinates |
| `clinical_significance` | ClinVar classification |
| `allele_freq` | Population allele frequency (GWAS Catalog) |
| `odds_ratio` | Effect size |
| `p_value` | GWAS association p-value |
| `mechanism` | Biological mechanism |
| `pubmed_id` | Source publication |

**T2D genes** (15): TCF7L2, PPARG, KCNJ11, GCK, HNF1A, CDKAL1, IGF2BP2, SLC30A8, HHEX, CDKN2B, ADCY5, GCKR, MTNR1B (×2), GCK

**CAD genes** (15): APOE, PCSK9, LPA, APOB (×2), LDLR, CETP, NPC1L1, HMGCR, SH2B3, CDKN2B, CENPW, CELSR2, TRIB1, FADS1

---

## Quick Start

```bash
git clone https://github.com/YOUR_USERNAME/gene-disease-explorer
cd gene-disease-explorer
pip install -r requirements.txt

# Streamlit (full dashboard)
streamlit run app.py

# Gradio (lightweight demo)
python gradio_demo.py
```

---

## Project Structure

```
gene_explorer/
├── app.py                         # Streamlit main app (5 tabs + sidebar)
├── gradio_demo.py                 # Gradio demo (HuggingFace-ready)
├── src/
│   ├── data_loader.py             # CSV loader, filter helper, summary stats
│   └── charts.py                  # All 6 Plotly chart builders
├── data/
│   └── clinvar_gwas_variants.csv  # 30 curated ClinVar + GWAS variants
├── requirements.txt
└── README.md
```

---

## Dashboard Features

- **Sidebar filters**: condition, clinical significance, allele frequency range, gene search
- **5-tab layout**: network, bubble, heatmap, bar charts, chromosome map
- **KPI cards**: live-updating variant counts and risk-tier breakdown
- **Variant detail table**: sortable, with CSV download
- **Streamlit caching**: fast re-renders on filter changes

---

## Deploy to HuggingFace Spaces (Gradio)

1. Create a new Space → Gradio SDK
2. Upload `gradio_demo.py`, `src/`, `data/`, `requirements.txt`
3. Space auto-launches — public URL in seconds

---

## Tech Stack

| Layer | Library |
|---|---|
| Main UI | Streamlit |
| Demo UI | Gradio |
| Charts | Plotly |
| Network | NetworkX |
| Data | pandas, numpy |

---

## Background

Built from real-world clinical genomics experience processing patient VCF files, gene–disease associations, and generating clinical visualizations at GenepoweRx. All rsIDs, odds ratios, and p-values are sourced from peer-reviewed GWAS studies and ClinVar records.

---

## License

MIT
