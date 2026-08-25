"""Interactive HTML report generation for HELIOS runs."""

from html import escape
from pathlib import Path

from organelle_pipeline.models import IsomerQuantResult, VariantCall


def write_html_report(
    path: Path,
    sample_name: str,
    isomer_result: IsomerQuantResult,
    heteroplasmy_calls: list[VariantCall],
) -> None:
    go = __import__("plotly.graph_objects", fromlist=["Figure"])
    plot = __import__("plotly.offline", fromlist=["plot"]).plot

    path.parent.mkdir(parents=True, exist_ok=True)

    # ---- isomer proportion chart ------------------------------------------
    isomer_fig = go.Figure(
        data=[
            go.Bar(
                x=[isomer_result.isomer_a_name, isomer_result.isomer_b_name],
                y=[
                    isomer_result.isomer_a_fraction * 100.0,
                    isomer_result.isomer_b_fraction * 100.0,
                ],
                marker_color=["#1f77b4", "#ff7f0e"],
            )
        ]
    )
    isomer_fig.update_layout(
        title="Isomer Proportions",
        xaxis_title="Isomer",
        yaxis_title="Percent of informative reads",
        yaxis_range=[0, 100],
    )
    isomer_div = plot(isomer_fig, include_plotlyjs="cdn", output_type="div")

    ssc_text = (
        f"{isomer_result.ssc_start}-{isomer_result.ssc_end} ({isomer_result.ssc_source})"
        if isomer_result.ssc_start is not None
        else "not determined"
    )
    agreement_text = (
        f"{isomer_result.remap_agreement_delta:.4f}"
        if isomer_result.remap_agreement_delta is not None
        else "NA"
    )
    metrics = [
        (
            "Informative Isomer Reads",
            str(isomer_result.assigned_a + isomer_result.assigned_b),
        ),
        ("Total Reads Seen", str(isomer_result.total_reads_seen)),
        ("Estimation Method", isomer_result.method),
        ("Dual-Ref Agreement Δ", agreement_text),
        ("SSC Region", ssc_text),
        ("Heteroplasmy Calls", str(len(heteroplasmy_calls))),
    ]

    # ---- heteroplasmy scatter ---------------------------------------------
    variant_fig = go.Figure()
    if heteroplasmy_calls:
        variant_fig.add_trace(
            go.Scatter(
                x=[call.position for call in heteroplasmy_calls],
                y=[call.alt_fraction * 100.0 for call in heteroplasmy_calls],
                mode="markers",
                text=[f"{call.ref}>{call.alt}" for call in heteroplasmy_calls],
                marker={"color": "#2ca02c", "size": 8},
            )
        )
    variant_fig.update_layout(
        title="Heteroplasmy Allele Fractions",
        xaxis_title="Genome position",
        yaxis_title="Alt allele fraction (%)",
    )
    variant_div = plot(variant_fig, include_plotlyjs=True, output_type="div")

    metrics_html = "".join(
        '<div class="tile"><div>{}</div><div class="value">{}</div></div>'.format(
            escape(label), escape(value)
        )
        for label, value in metrics
    )
    assumptions_html = "".join(f"<li>{escape(item)}</li>" for item in isomer_result.assumptions)
    variant_table = _variant_table(heteroplasmy_calls)

    html = f"""
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>HELIOS Report - {escape(sample_name)}</title>
    <style>
      body {{
        font-family: "IBM Plex Sans", "Segoe UI", sans-serif;
        background: linear-gradient(180deg, #f4f8fb 0%, #e7eff7 100%);
        margin: 0;
        color: #0f172a;
      }}
      main {{ max-width: 1100px; margin: 0 auto; padding: 24px 16px 48px; }}
      h1 {{ margin-bottom: 8px; }}
      section {{
        background: rgba(255, 255, 255, 0.85);
        border-radius: 14px;
        padding: 16px;
        margin-top: 16px;
        box-shadow: 0 8px 24px rgba(15, 23, 42, 0.08);
      }}
      table {{ border-collapse: collapse; width: 100%; font-size: 0.85rem; }}
      th, td {{ border: 1px solid #cbd5e1; padding: 6px 8px; text-align: left; }}
      th {{ background: #f8fafc; }}
      .kpi {{
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
        gap: 10px;
      }}
      .tile {{
        background: #f8fafc;
        border: 1px solid #dbeafe;
        border-radius: 12px;
        padding: 10px;
      }}
      .tile .value {{
        font-size: 1.25rem;
        font-weight: 700;
        word-break: break-all;
      }}
      .flag-homopolymer {{ color: #b45309; font-weight: 600; }}
      .flag-strand {{ color: #b91c1c; font-weight: 600; }}
    </style>
  </head>
  <body>
    <main>
      <h1>HELIOS Heteroplasmy &amp; Isomer Quantification</h1>
      <p><strong>Sample:</strong> {escape(sample_name)}</p>

      <section>
        <h2>Run Metrics</h2>
        <div class="kpi">{metrics_html}</div>
      </section>

      <section>
        <h2>Isomer Proportion Plot</h2>
        {isomer_div}
      </section>

      <section>
        <h2>Heteroplasmy Plot</h2>
        {variant_div}
      </section>

      <section>
        <h2>Isomer Assumptions</h2>
        <ul>{assumptions_html}</ul>
      </section>

      <section>
        <h2>Heteroplasmy Calls</h2>
        {variant_table}
      </section>
    </main>
  </body>
</html>
"""

    path.write_text(html, encoding="utf-8")


def _variant_table(calls: list[VariantCall]) -> str:
    if not calls:
        return "<p>No heteroplasmy variants passed thresholds.</p>"
    rows = []
    for call in calls:
        strand_flag = ""
        if (
            call.strand_bias_pvalue is not None
            and call.strand_bias_pvalue < 1e-3
            and (call.strand_alt_fwd == 0 or call.strand_alt_rev == 0)
        ):
            strand_flag = ' class="flag-strand"'
        homopolymer_flag = ' class="flag-homopolymer"' if call.homopolymer_context else ""
        rows.append(
            "<tr>"
            f"<td>{escape(call.contig)}</td>"
            f"<td>{call.position}</td>"
            f"<td>{escape(call.ref)}</td>"
            f"<td>{escape(call.alt)}</td>"
            f"<td>{call.depth}</td>"
            f"<td>{call.alt_count}</td>"
            f"<td>{round(call.alt_fraction * 100.0, 4)}</td>"
            f"<td>{escape(call.variant_type)}</td>"
            f"<td>{call.strand_alt_fwd}/{call.strand_alt_rev}</td>"
            f"<td{strand_flag}>"
            + (f"{call.strand_bias_pvalue:.3g}" if call.strand_bias_pvalue is not None else "NA")
            + "</td>"
            f"<td{homopolymer_flag}>{'yes' if call.homopolymer_context else 'no'}</td>"
            "<td>"
            + (
                f"{call.mean_alt_read_position:.3f}"
                if call.mean_alt_read_position is not None
                else "NA"
            )
            + "</td>"
            "</tr>"
        )
    body = "".join(rows)
    return (
        "<table><thead><tr>"
        "<th>Contig</th><th>Position</th><th>Ref</th><th>Alt</th><th>Depth</th>"
        "<th>Alt Count</th><th>Alt %</th><th>Type</th><th>Alt Strand F/R</th>"
        "<th>Strand-Bias p</th><th>Homopolymer</th><th>Mean Alt Read Pos</th>"
        "</tr></thead><tbody>"
        f"{body}</tbody></table>"
    )
