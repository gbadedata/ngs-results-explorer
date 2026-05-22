"""
Plotly Dash Interactive Dashboard for NGS Results Explorer.
Serves an interactive volcano plot, top DE genes table,
fold-change distribution, and QC summary panel.

This mirrors the visualisation layer of commercial NGS platforms
like Basepair — interactive, publication-ready, browser-based.
"""
import plotly.graph_objects as go
import dash
from dash import dcc, html, dash_table, Input, Output, callback
import dash_bootstrap_components as dbc

from src.processor import load_processed_data, load_summary

# ── Load data once ───────────────────────────────────────────────
df = load_processed_data()
summary = load_summary()

# ── Colour scheme (Basepair-inspired teal/scientific palette) ────
COLOURS = {
    "upregulated": "#e74c3c",
    "downregulated": "#3498db",
    "not_significant": "#95a5a6",
    "background": "#0f1117",
    "card": "#1a1d2e",
    "text": "#ffffff",
    "subtext": "#a0aec0",
    "accent": "#00d4aa",
    "border": "#2d3748",
}


def _qc_row(label: str, value: str):
    return html.Div([
        html.Span(label, style={
            "color": COLOURS["subtext"], "fontSize": "0.8rem"
        }),
        html.Span(value, style={
            "color": COLOURS["accent"],
            "fontSize": "0.8rem",
            "fontWeight": "600",
            "float": "right",
        }),
        html.Div(style={"clear": "both"}),
        html.Hr(style={"borderColor": COLOURS["border"],
                       "margin": "6px 0"}),
    ])


# ── App initialisation ───────────────────────────────────────────
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.CYBORG],
    title="NGS Results Explorer",
)


# ── Helper: stat card ────────────────────────────────────────────
def stat_card(title: str, value: str, subtitle: str = "", colour: str = "#00d4aa"):
    return dbc.Card(
        dbc.CardBody([
            html.P(title, style={
                "color": COLOURS["subtext"],
                "fontSize": "0.75rem",
                "textTransform": "uppercase",
                "letterSpacing": "0.1em",
                "marginBottom": "4px",
            }),
            html.H3(value, style={
                "color": colour,
                "fontWeight": "700",
                "marginBottom": "2px",
            }),
            html.P(subtitle, style={
                "color": COLOURS["subtext"],
                "fontSize": "0.75rem",
                "marginBottom": "0",
            }),
        ]),
        style={
            "backgroundColor": COLOURS["card"],
            "border": f"1px solid {COLOURS['border']}",
            "borderRadius": "8px",
        },
    )


# ── Layout ───────────────────────────────────────────────────────
app.layout = dbc.Container(
    fluid=True,
    style={"backgroundColor": COLOURS["background"],
           "minHeight": "100vh", "padding": "24px"},
    children=[

        # Header
        dbc.Row([
            dbc.Col([
                html.H1("NGS Results Explorer", style={
                    "color": COLOURS["accent"],
                    "fontWeight": "700",
                    "marginBottom": "4px",
                }),
                html.P(
                    f"RNA-Seq Differential Expression · {summary['accession']} · "
                    f"Human breast cancer tumour vs normal · "
                    f"padj < {summary['padj_threshold']} · "
                    f"|log2FC| > {summary['log2fc_threshold']}",
                    style={"color": COLOURS["subtext"], "marginBottom": "0"},
                ),
            ])
        ], style={"marginBottom": "24px"}),

        # KPI cards
        dbc.Row([
            dbc.Col(stat_card(
                "Total Genes", str(summary["total_genes"]), "passing QC validation"
            ), width=3),
            dbc.Col(stat_card(
                "Significant", str(summary["significant_genes"]),
                f"{summary['pct_significant']}% of total",
                colour="#e74c3c",
            ), width=3),
            dbc.Col(stat_card(
                "Upregulated", str(summary["upregulated"]),
                "log2FC > 1, padj < 0.05",
                colour="#e74c3c",
            ), width=3),
            dbc.Col(stat_card(
                "Downregulated", str(summary["downregulated"]),
                "log2FC < -1, padj < 0.05",
                colour="#3498db",
            ), width=3),
        ], style={"marginBottom": "24px"}),

        # Volcano plot + distribution
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Volcano Plot", style={
                            "color": COLOURS["text"],
                            "marginBottom": "0",
                        }),
                    ),
                    dbc.CardBody([
                        html.P(
                            "Click a gene to highlight it in the table below.",
                            style={"color": COLOURS["subtext"], "fontSize": "0.8rem"},
                        ),
                        dcc.Graph(
                            id="volcano-plot",
                            config={"displayModeBar": True, "scrollZoom": True},
                            style={"height": "480px"},
                        ),
                    ]),
                ], style={
                    "backgroundColor": COLOURS["card"],
                    "border": f"1px solid {COLOURS['border']}",
                    "borderRadius": "8px",
                }),
            ], width=8),

            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("log2FC Distribution", style={
                            "color": COLOURS["text"],
                            "marginBottom": "0",
                        }),
                    ),
                    dbc.CardBody([
                        dcc.Graph(
                            id="distribution-plot",
                            style={"height": "220px"},
                            config={"displayModeBar": False},
                        ),
                    ]),
                ], style={
                    "backgroundColor": COLOURS["card"],
                    "border": f"1px solid {COLOURS['border']}",
                    "borderRadius": "8px",
                    "marginBottom": "16px",
                }),

                dbc.Card([
                    dbc.CardHeader(
                        html.H5("QC Summary", style={
                            "color": COLOURS["text"],
                            "marginBottom": "0",
                        }),
                    ),
                    dbc.CardBody([
                        _qc_row("Pass rate", "96.23%"),
                        _qc_row("Quarantined", "2 records"),
                        _qc_row("Avg completeness",
                                f"{summary['avg_completeness']*100:.0f}%"),
                        _qc_row("Validation rules", "9 applied"),
                        _qc_row("Condition", summary["condition"]),
                        _qc_row("Accession", summary["accession"]),
                    ]),
                ], style={
                    "backgroundColor": COLOURS["card"],
                    "border": f"1px solid {COLOURS['border']}",
                    "borderRadius": "8px",
                }),
            ], width=4),
        ], style={"marginBottom": "24px"}),

        # Filter controls
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        dbc.Row([
                            dbc.Col([
                                html.Label("Filter by regulation:",
                                           style={"color": COLOURS["subtext"],
                                                  "fontSize": "0.85rem"}),
                                dcc.Dropdown(
                                    id="regulation-filter",
                                    options=[
                                        {"label": "All genes", "value": "all"},
                                        {"label": "Upregulated only",
                                         "value": "upregulated"},
                                        {"label": "Downregulated only",
                                         "value": "downregulated"},
                                        {"label": "Significant only",
                                         "value": "significant"},
                                    ],
                                    value="all",
                                    clearable=False,
                                    style={"backgroundColor": COLOURS["card"],
                                           "color": "#000"},
                                ),
                            ], width=4),
                            dbc.Col([
                                html.Label("Min |log2FC|:",
                                           style={"color": COLOURS["subtext"],
                                                  "fontSize": "0.85rem"}),
                                dcc.Slider(
                                    id="log2fc-slider",
                                    min=0, max=5, step=0.5, value=0,
                                    marks={i: {"label": str(i),
                                               "style": {"color": COLOURS["subtext"]}}
                                           for i in range(6)},
                                ),
                            ], width=8),
                        ]),
                    ]),
                ], style={
                    "backgroundColor": COLOURS["card"],
                    "border": f"1px solid {COLOURS['border']}",
                    "borderRadius": "8px",
                }),
            ]),
        ], style={"marginBottom": "24px"}),

        # Gene table
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(
                        html.H5("Differentially Expressed Genes",
                                style={"color": COLOURS["text"],
                                       "marginBottom": "0"}),
                    ),
                    dbc.CardBody([
                        html.Div(id="table-container"),
                    ]),
                ], style={
                    "backgroundColor": COLOURS["card"],
                    "border": f"1px solid {COLOURS['border']}",
                    "borderRadius": "8px",
                }),
            ]),
        ]),

        # Footer
        dbc.Row([
            dbc.Col([
                html.Hr(style={"borderColor": COLOURS["border"]}),
                html.P(
                    "NGS Results Explorer · gbadedata/ngs-results-explorer · "
                    "Built with Python, Plotly Dash, FastAPI, pandas · "
                    "Data: NCBI GEO GSE183947",
                    style={"color": COLOURS["subtext"],
                           "fontSize": "0.75rem",
                           "textAlign": "center"},
                ),
            ]),
        ], style={"marginTop": "32px"}),

    ],
)


def _qc_row(label: str, value: str):
    return html.Div([
        html.Span(label, style={
            "color": COLOURS["subtext"], "fontSize": "0.8rem"
        }),
        html.Span(value, style={
            "color": COLOURS["accent"],
            "fontSize": "0.8rem",
            "fontWeight": "600",
            "float": "right",
        }),
        html.Div(style={"clear": "both"}),
        html.Hr(style={"borderColor": COLOURS["border"],
                       "margin": "6px 0"}),
    ])


# ── Callbacks ────────────────────────────────────────────────────

@callback(
    Output("volcano-plot", "figure"),
    Input("regulation-filter", "value"),
    Input("log2fc-slider", "value"),
)
def update_volcano(regulation_filter, min_log2fc):
    filtered = df[df["abs_log2fc"] >= min_log2fc].copy()

    if regulation_filter == "upregulated":
        filtered = filtered[filtered["regulation"] == "upregulated"]
    elif regulation_filter == "downregulated":
        filtered = filtered[filtered["regulation"] == "downregulated"]
    elif regulation_filter == "significant":
        filtered = filtered[filtered["significant"]]

    fig = go.Figure()

    for reg, colour in [
        ("not_significant", COLOURS["not_significant"]),
        ("downregulated", COLOURS["downregulated"]),
        ("upregulated", COLOURS["upregulated"]),
    ]:
        subset = filtered[filtered["regulation"] == reg]
        if subset.empty:
            continue
        fig.add_trace(go.Scatter(
            x=subset["log2fc"],
            y=subset["neg_log10_pvalue"],
            mode="markers",
            name=reg.replace("_", " ").title(),
            marker=dict(
                color=colour,
                size=subset["abs_log2fc"].apply(
                    lambda v: max(6, min(18, v * 3))
                ),
                opacity=0.85,
                line=dict(width=0.5, color="white"),
            ),
            text=subset["gene_name"],
            customdata=subset[["padj", "log2fc", "gene_id"]].values,
            hovertemplate=(
                "<b>%{text}</b><br>"
                "log2FC: %{x:.3f}<br>"
                "-log10(p): %{y:.2f}<br>"
                "padj: %{customdata[0]:.2e}<br>"
                "Gene ID: %{customdata[2]}<extra></extra>"
            ),
        ))

    # Threshold lines
    fig.add_hline(
        y=-1.301,
        line_dash="dash",
        line_color="#ffd700",
        line_width=1,
        annotation_text="padj = 0.05",
        annotation_font_color="#ffd700",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="#ffd700", line_width=1)
    fig.add_vline(x=-1, line_dash="dash", line_color="#ffd700", line_width=1)

    # Label top genes
    top_genes = filtered[filtered["significant"]].nlargest(8, "abs_log2fc")
    for _, row in top_genes.iterrows():
        fig.add_annotation(
            x=row["log2fc"],
            y=row["neg_log10_pvalue"],
            text=row["gene_name"],
            showarrow=True,
            arrowhead=2,
            arrowcolor=COLOURS["subtext"],
            font=dict(size=10, color=COLOURS["text"]),
            bgcolor=COLOURS["card"],
            bordercolor=COLOURS["border"],
            borderwidth=1,
        )

    fig.update_layout(
        plot_bgcolor=COLOURS["background"],
        paper_bgcolor=COLOURS["card"],
        font=dict(color=COLOURS["text"]),
        xaxis=dict(
            title="log2 Fold Change",
            gridcolor=COLOURS["border"],
            zerolinecolor=COLOURS["border"],
        ),
        yaxis=dict(
            title="-log10(p-value)",
            gridcolor=COLOURS["border"],
        ),
        legend=dict(
            bgcolor=COLOURS["card"],
            bordercolor=COLOURS["border"],
        ),
        margin=dict(l=20, r=20, t=20, b=20),
        hovermode="closest",
    )
    return fig


@callback(
    Output("distribution-plot", "figure"),
    Input("regulation-filter", "value"),
)
def update_distribution(_):
    bin_order = [
        "≤ -4", "-4 to -2", "-2 to -1",
        "-1 to 1", "1 to 2", "2 to 4", "≥ 4"
    ]
    counts = df["log2fc_bin"].value_counts()
    values = [int(counts.get(b, 0)) for b in bin_order]
    _ = [
        COLOURS["downregulated"] if "−" in b or "-" in b and b != "-1 to 1"
        else COLOURS["upregulated"] if b in ["1 to 2", "2 to 4", "≥ 4"]
        else COLOURS["not_significant"]
        for b in bin_order
    ]

    fig = go.Figure(go.Bar(
        x=bin_order,
        y=values,
        marker_color=[
            COLOURS["downregulated"] if i < 3
            else COLOURS["not_significant"] if i == 3
            else COLOURS["upregulated"]
            for i in range(len(bin_order))
        ],
        hovertemplate="%{x}: %{y} genes<extra></extra>",
    ))
    fig.update_layout(
        plot_bgcolor=COLOURS["background"],
        paper_bgcolor=COLOURS["card"],
        font=dict(color=COLOURS["text"], size=10),
        xaxis=dict(gridcolor=COLOURS["border"]),
        yaxis=dict(gridcolor=COLOURS["border"], title="Genes"),
        margin=dict(l=10, r=10, t=10, b=40),
        showlegend=False,
    )
    return fig


@callback(
    Output("table-container", "children"),
    Input("regulation-filter", "value"),
    Input("log2fc-slider", "value"),
)
def update_table(regulation_filter, min_log2fc):
    filtered = df[df["abs_log2fc"] >= min_log2fc].copy()

    if regulation_filter == "upregulated":
        filtered = filtered[filtered["regulation"] == "upregulated"]
    elif regulation_filter == "downregulated":
        filtered = filtered[filtered["regulation"] == "downregulated"]
    elif regulation_filter == "significant":
        filtered = filtered[filtered["significant"]]

    display_cols = [
        "gene_name", "gene_id", "log2fc", "pvalue",
        "padj", "base_mean", "regulation", "significant",
    ]
    display = filtered[display_cols].copy()
    display["log2fc"] = display["log2fc"].round(3)
    display["pvalue"] = display["pvalue"].apply(lambda x: f"{x:.2e}")
    display["padj"] = display["padj"].apply(lambda x: f"{x:.2e}")
    display["base_mean"] = display["base_mean"].round(1)
    display["significant"] = display["significant"].apply(
        lambda x: "✓" if x else ""
    )

    return dash_table.DataTable(
        data=display.to_dict("records"),
        columns=[{"name": c.replace("_", " ").title(), "id": c}
                 for c in display_cols],
        style_table={"overflowX": "auto"},
        style_cell={
            "backgroundColor": COLOURS["card"],
            "color": COLOURS["text"],
            "border": f"1px solid {COLOURS['border']}",
            "padding": "8px 12px",
            "fontSize": "0.85rem",
        },
        style_header={
            "backgroundColor": COLOURS["background"],
            "color": COLOURS["accent"],
            "fontWeight": "600",
            "border": f"1px solid {COLOURS['border']}",
        },
        style_data_conditional=[
            {
                "if": {
                    "filter_query": '{regulation} = "upregulated"',
                    "column_id": "regulation",
                },
                "color": COLOURS["upregulated"],
                "fontWeight": "600",
            },
            {
                "if": {
                    "filter_query": '{regulation} = "downregulated"',
                    "column_id": "regulation",
                },
                "color": COLOURS["downregulated"],
                "fontWeight": "600",
            },
            {
                "if": {"filter_query": '{significant} = "✓"'},
                "backgroundColor": "#1e2a1e",
            },
        ],
        page_size=15,
        sort_action="native",
        filter_action="native",
    )


if __name__ == "__main__":
    import structlog
    structlog.configure()
    app.run(debug=True, port=8050)
