"""
gradio_demo.py  —  Gene–Disease Association Explorer (Gradio)
=============================================================
Lightweight Gradio demo — great for HuggingFace Spaces deployment.
Exposes the same charts via dropdowns and sliders.

Run:
    python gradio_demo.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import gradio as gr
import pandas as pd
from src.data_loader import load_data, filter_data, SIG_ORDER, CONDITION_COLORS
from src.charts import (
    build_network,
    build_bubble_scatter,
    build_heatmap,
    build_variants_per_gene,
    build_sig_distribution,
    build_chromosome_scatter,
)

df_full = load_data()
ALL_CONDITIONS = df_full["condition"].unique().tolist()
ALL_GENES      = sorted(df_full["gene"].unique().tolist())


# ── Core update function ──────────────────────────────────────────────────────

def update_all(conditions, sig_levels, af_min, af_max, gene_filter):
    df = filter_data(
        df_full,
        conditions=conditions or None,
        sig_levels=sig_levels or None,
        min_af=af_min,
        max_af=af_max,
    )
    if gene_filter and gene_filter != "All":
        df = df[df["gene"] == gene_filter]

    if df.empty:
        empty = gr.update(value=None)
        return empty, empty, empty, empty, empty, empty, "⚠️ No variants match the current filters."

    # summary table
    summary = pd.DataFrame([{
        "Metric": "Total Variants",    "Value": len(df),
    }, {
        "Metric": "Unique Genes",      "Value": df["gene"].nunique(),
    }, {
        "Metric": "Pathogenic",        "Value": len(df[df["clinical_significance"] == "Pathogenic"]),
    }, {
        "Metric": "Likely Pathogenic", "Value": len(df[df["clinical_significance"] == "Likely Pathogenic"]),
    }, {
        "Metric": "Risk Factor",       "Value": len(df[df["clinical_significance"] == "Risk Factor"]),
    }])

    summary_md = "\n".join(
        f"| **{r['Metric']}** | {r['Value']} |"
        for _, r in summary.iterrows()
    )
    summary_md = "| Metric | Value |\n|---|---|\n" + summary_md

    return (
        build_network(df),
        build_bubble_scatter(df),
        build_heatmap(df),
        build_variants_per_gene(df),
        build_chromosome_scatter(df),
        build_sig_distribution(df),
        summary_md,
    )


# ── Gradio UI ─────────────────────────────────────────────────────────────────

with gr.Blocks(
    title="Gene–Disease Association Explorer",
    theme=gr.themes.Soft(primary_hue="blue"),
    css="""
        .gradio-container { max-width: 1200px !important; }
        h1 { color: #1e3a5f; }
        .gr-button-primary { background: #2563eb !important; }
    """
) as demo:

    gr.Markdown("""
    # 🧬 Gene–Disease Association Explorer
    Explore gene–rsID relationships for **Type 2 Diabetes** and **Coronary Artery Disease**
    using curated ClinVar + GWAS Catalog data (30 literature-validated variants).

    > Built by **Rayala Madhu Bhanu Varma** · Lead Data Scientist · GenepoweRx
    """)

    # ── Controls ──────────────────────────────────────────────────────────────
    with gr.Row():
        with gr.Column(scale=1):
            cond_cb = gr.CheckboxGroup(
                label="🏥 Conditions",
                choices=ALL_CONDITIONS,
                value=ALL_CONDITIONS,
            )
            sig_cb = gr.CheckboxGroup(
                label="🔬 Clinical Significance",
                choices=SIG_ORDER,
                value=SIG_ORDER,
            )
        with gr.Column(scale=1):
            af_slider = gr.Slider(
                label="Allele Frequency Range",
                minimum=0.0, maximum=1.0,
                value=1.0, step=0.01,
            )
            af_min_slider = gr.Slider(
                label="Allele Frequency Min",
                minimum=0.0, maximum=1.0,
                value=0.0, step=0.01,
            )
            gene_dd = gr.Dropdown(
                label="🧬 Filter by Gene",
                choices=["All"] + ALL_GENES,
                value="All",
            )
            run_btn = gr.Button("🔄 Update Charts", variant="primary")

        with gr.Column(scale=1):
            summary_md = gr.Markdown(label="Summary", value="Click **Update Charts** to begin.")

    gr.Markdown("---")

    # ── Tabs ──────────────────────────────────────────────────────────────────
    with gr.Tabs():
        with gr.Tab("🕸️ Network Graph"):
            gr.Markdown("Gene nodes (circles) connect to disease nodes (diamonds). Edge thickness = shared variants.")
            net_plot = gr.Plot(label="Gene–Disease Network")

        with gr.Tab("🫧 Bubble Chart"):
            gr.Markdown("X = Allele Frequency · Y = Odds Ratio · Bubble size = significance. Hover for details.")
            bubble_plot = gr.Plot(label="OR vs AF Bubble Chart")

        with gr.Tab("🗺️ Heatmap"):
            gr.Markdown("Mean Odds Ratio per gene × condition. Darker red = stronger association.")
            heat_plot = gr.Plot(label="Gene × Disease Heatmap")

        with gr.Tab("📊 Bar Charts"):
            with gr.Row():
                bar_plot  = gr.Plot(label="Variants per Gene")
                sig_plot  = gr.Plot(label="Clinical Significance Distribution")

        with gr.Tab("🧬 Chromosome Map"):
            gr.Markdown("Manhattan-style. Dashed lines = significance thresholds.")
            chrom_plot = gr.Plot(label="Chromosomal Distribution")

    # ── Wire up ───────────────────────────────────────────────────────────────
    outputs = [net_plot, bubble_plot, heat_plot, bar_plot, chrom_plot, sig_plot, summary_md]
    inputs  = [cond_cb, sig_cb, af_min_slider, af_slider, gene_dd]

    run_btn.click(fn=update_all, inputs=inputs, outputs=outputs)

    # auto-load on start
    demo.load(
        fn=lambda: update_all(ALL_CONDITIONS, SIG_ORDER, 0.0, 1.0, "All"),
        outputs=outputs,
    )


if __name__ == "__main__":
    demo.launch(share=False, server_port=7860)
