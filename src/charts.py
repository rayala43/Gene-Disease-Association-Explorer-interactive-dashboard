"""
src/charts.py
-------------
All Plotly chart builders. Each function returns a go.Figure.
Charts are consumed by both app.py (Streamlit) and gradio_demo.py.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import networkx as nx

from src.data_loader import CONDITION_COLORS, SIG_COLORS, SIG_ORDER


# ── 1. Network graph ─────────────────────────────────────────────────────────

def build_network(df: pd.DataFrame) -> go.Figure:
    """
    Bipartite network: gene nodes ←→ condition nodes, edges weighted by
    number of variants. Node size = degree; edge width = variant count.
    """
    G = nx.Graph()

    # add gene nodes
    for gene in df["gene"].unique():
        G.add_node(gene, node_type="gene")

    # add condition nodes
    for cond in df["condition"].unique():
        G.add_node(cond, node_type="condition")

    # add edges
    edge_weights: dict[tuple, int] = {}
    for _, row in df.iterrows():
        key = (row["gene"], row["condition"])
        edge_weights[key] = edge_weights.get(key, 0) + 1

    for (gene, cond), w in edge_weights.items():
        G.add_edge(gene, cond, weight=w)

    # layout: genes on left half, conditions on right
    pos = {}
    genes      = [n for n, d in G.nodes(data=True) if d["node_type"] == "gene"]
    conditions = [n for n, d in G.nodes(data=True) if d["node_type"] == "condition"]

    for i, g in enumerate(sorted(genes)):
        pos[g] = (0.05, (i + 1) / (len(genes) + 1))
    for i, c in enumerate(sorted(conditions)):
        pos[c] = (0.95, (i + 1) / (len(conditions) + 1))

    # edge traces
    edge_traces = []
    for (u, v, data) in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        w = data["weight"]
        edge_traces.append(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None],
            mode="lines",
            line=dict(width=max(1, w * 2), color="rgba(150,150,180,0.45)"),
            hoverinfo="skip",
            showlegend=False,
        ))

    # node traces — genes
    gene_degrees = {g: G.degree(g) for g in genes}
    gx = [pos[g][0] for g in genes]
    gy = [pos[g][1] for g in genes]
    gsz = [20 + gene_degrees[g] * 8 for g in genes]

    gene_trace = go.Scatter(
        x=gx, y=gy,
        mode="markers+text",
        text=genes,
        textposition="middle left",
        textfont=dict(size=12, color="#1e293b"),
        marker=dict(size=gsz, color="#3B82F6", line=dict(width=2, color="white"),
                    symbol="circle"),
        hovertemplate="<b>%{text}</b><br>Connections: %{marker.size}<extra></extra>",
        name="Gene",
    )

    # node traces — conditions
    cx = [pos[c][0] for c in conditions]
    cy = [pos[c][1] for c in conditions]
    ccolors = [CONDITION_COLORS.get(c, "#8B5CF6") for c in conditions]
    cond_trace = go.Scatter(
        x=cx, y=cy,
        mode="markers+text",
        text=conditions,
        textposition="middle right",
        textfont=dict(size=13, color="#1e293b", family="Arial Black"),
        marker=dict(size=36, color=ccolors, line=dict(width=2, color="white"),
                    symbol="diamond"),
        hovertemplate="<b>%{text}</b><extra></extra>",
        name="Condition",
    )

    fig = go.Figure(data=edge_traces + [gene_trace, cond_trace])
    fig.update_layout(
        title=dict(text="Gene–Disease Association Network", font=dict(size=16), x=0.5),
        showlegend=False,
        xaxis=dict(visible=False, range=[-0.1, 1.1]),
        yaxis=dict(visible=False, range=[0, 1]),
        margin=dict(l=10, r=10, t=50, b=10),
        height=520,
        plot_bgcolor="white",
        paper_bgcolor="white",
        hoverlabel=dict(bgcolor="white", font_size=12),
    )
    return fig


# ── 2. Manhattan-style scatter: OR vs -log10(p) ──────────────────────────────

def build_bubble_scatter(df: pd.DataFrame) -> go.Figure:
    """
    Bubble scatter: x = allele frequency, y = odds ratio,
    size = -log10(p-value), color = clinical significance.
    """
    fig = go.Figure()

    for sig in SIG_ORDER:
        sub = df[df["clinical_significance"] == sig]
        if sub.empty:
            continue
        fig.add_trace(go.Scatter(
            x=sub["allele_freq"],
            y=sub["odds_ratio"],
            mode="markers",
            name=sig,
            marker=dict(
                size=sub["neg_log10_p"].clip(5, 40),
                color=SIG_COLORS.get(sig, "#94A3B8"),
                opacity=0.8,
                line=dict(width=1, color="white"),
                sizemode="diameter",
            ),
            customdata=np.stack([
                sub["rsid"], sub["gene"], sub["condition"],
                sub["mechanism"], sub["neg_log10_p"].round(1),
                sub["p_value"].apply(lambda x: f"{x:.2e}")
            ], axis=1),
            hovertemplate=(
                "<b>%{customdata[0]}</b> — %{customdata[1]}<br>"
                "Condition: %{customdata[2]}<br>"
                "OR: %{y:.2f} &nbsp;|&nbsp; AF: %{x:.3f}<br>"
                "p-value: %{customdata[5]} &nbsp;(-log₁₀p = %{customdata[4]})<br>"
                "Mechanism: %{customdata[3]}"
                "<extra></extra>"
            ),
        ))

    fig.add_hline(y=1.0, line_dash="dash", line_color="gray",
                  annotation_text="OR = 1 (no effect)", annotation_position="top right")

    fig.update_layout(
        title=dict(text="Odds Ratio vs Allele Frequency<br>"
                        "<sup>Bubble size = statistical significance (-log₁₀ p)</sup>",
                   font=dict(size=15), x=0.5),
        xaxis=dict(title="Allele Frequency (population)", tickformat=".0%", gridcolor="#f1f5f9"),
        yaxis=dict(title="Odds Ratio (effect size)", gridcolor="#f1f5f9"),
        legend=dict(title="Clinical Significance", orientation="v"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=480,
        margin=dict(l=60, r=20, t=80, b=60),
        hoverlabel=dict(bgcolor="white", font_size=12),
    )
    return fig


# ── 3. Risk heatmap: gene × condition ────────────────────────────────────────

def build_heatmap(df: pd.DataFrame) -> go.Figure:
    """
    Heatmap of mean odds ratio per gene × condition cell.
    Genes on Y-axis, conditions on X-axis.
    """
    pivot = df.pivot_table(
        index="gene", columns="condition",
        values="odds_ratio", aggfunc="mean"
    ).fillna(0)

    # sort genes by max OR across conditions
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=True).index]

    conditions = pivot.columns.tolist()
    genes      = pivot.index.tolist()
    z          = pivot.values

    # custom text for each cell
    text_matrix = []
    for i, gene in enumerate(genes):
        row = []
        for j, cond in enumerate(conditions):
            val = z[i][j]
            row.append(f"{val:.2f}" if val > 0 else "—")
        text_matrix.append(row)

    fig = go.Figure(go.Heatmap(
        z=z,
        x=conditions,
        y=genes,
        text=text_matrix,
        texttemplate="%{text}",
        textfont=dict(size=10),
        colorscale=[
            [0.0,  "#f0f9ff"],
            [0.3,  "#bae6fd"],
            [0.6,  "#f97316"],
            [1.0,  "#dc2626"],
        ],
        colorbar=dict(
            title="Mean OR",
            tickformat=".2f",
            len=0.8,
        ),
        hovertemplate="<b>%{y}</b> × %{x}<br>Mean OR: %{z:.3f}<extra></extra>",
        zmin=0.8, zmax=df["odds_ratio"].max() + 0.1,
    ))

    fig.update_layout(
        title=dict(text="Gene × Condition Risk Heatmap<br>"
                        "<sup>Mean Odds Ratio — deeper red = stronger association</sup>",
                   font=dict(size=15), x=0.5),
        xaxis=dict(title="", tickfont=dict(size=12)),
        yaxis=dict(title="", tickfont=dict(size=11), autorange="reversed"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=max(380, len(genes) * 26 + 120),
        margin=dict(l=90, r=20, t=80, b=40),
    )
    return fig


# ── 4. Bar charts ─────────────────────────────────────────────────────────────

def build_variants_per_gene(df: pd.DataFrame) -> go.Figure:
    counts = (
        df.groupby(["gene", "condition"])
          .size()
          .reset_index(name="count")
          .sort_values("count", ascending=True)
    )
    fig = px.bar(
        counts, x="count", y="gene", color="condition",
        orientation="h",
        color_discrete_map=CONDITION_COLORS,
        text="count",
        labels={"count": "Number of Variants", "gene": "Gene", "condition": "Condition"},
    )
    fig.update_traces(textposition="outside")
    fig.update_layout(
        title=dict(text="Variants per Gene by Condition", font=dict(size=15), x=0.5),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend_title="Condition",
        height=max(360, len(counts["gene"].unique()) * 26 + 120),
        margin=dict(l=80, r=40, t=60, b=40),
    )
    return fig


def build_sig_distribution(df: pd.DataFrame) -> go.Figure:
    counts = (
        df.groupby(["clinical_significance", "condition"])
          .size()
          .reset_index(name="count")
    )
    fig = px.bar(
        counts, x="clinical_significance", y="count",
        color="condition", barmode="group",
        color_discrete_map=CONDITION_COLORS,
        category_orders={"clinical_significance": SIG_ORDER},
        labels={"clinical_significance": "Clinical Significance",
                "count": "Number of Variants", "condition": "Condition"},
    )
    fig.update_layout(
        title=dict(text="Clinical Significance Distribution", font=dict(size=15), x=0.5),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend_title="Condition",
        height=380,
        margin=dict(l=60, r=20, t=60, b=80),
        xaxis_tickangle=-25,
    )
    return fig


def build_chromosome_scatter(df: pd.DataFrame) -> go.Figure:
    """Manhattan-style chromosome position plot."""
    chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chrom_map   = {c: i for i, c in enumerate(chrom_order)}
    df2 = df.copy()
    df2["chrom_num"] = df2["chromosome"].map(chrom_map).fillna(len(chrom_order))

    fig = go.Figure()
    for cond, color in CONDITION_COLORS.items():
        sub = df2[df2["condition"] == cond]
        fig.add_trace(go.Scatter(
            x=sub["chrom_num"] + np.random.uniform(-0.25, 0.25, len(sub)),
            y=sub["neg_log10_p"],
            mode="markers",
            name=cond,
            marker=dict(size=10, color=color, opacity=0.8,
                        line=dict(width=1, color="white")),
            customdata=np.stack([sub["rsid"], sub["gene"],
                                 sub["p_value"].apply(lambda x: f"{x:.2e}")], axis=1),
            hovertemplate=(
                "<b>%{customdata[0]}</b> — %{customdata[1]}<br>"
                "Chr %{x:.0f} | -log₁₀p = %{y:.1f}<br>"
                "p = %{customdata[2]}<extra></extra>"
            ),
        ))

    # genome-wide significance line
    fig.add_hline(y=-np.log10(5e-8), line_dash="dash", line_color="#dc2626",
                  annotation_text="GW significance (5×10⁻⁸)",
                  annotation_position="top right")
    fig.add_hline(y=-np.log10(1e-5), line_dash="dot", line_color="#f97316",
                  annotation_text="Suggestive (10⁻⁵)",
                  annotation_position="top right")

    fig.update_layout(
        title=dict(text="Chromosomal Distribution of Variants<br>"
                        "<sup>Y-axis = statistical significance (-log₁₀ p-value)</sup>",
                   font=dict(size=15), x=0.5),
        xaxis=dict(
            title="Chromosome",
            tickvals=list(range(len(chrom_order))),
            ticktext=chrom_order,
            tickfont=dict(size=10),
            gridcolor="#f1f5f9",
        ),
        yaxis=dict(title="-log₁₀(p-value)", gridcolor="#f1f5f9"),
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend_title="Condition",
        height=420,
        margin=dict(l=60, r=20, t=80, b=60),
        hoverlabel=dict(bgcolor="white", font_size=12),
    )
    return fig
