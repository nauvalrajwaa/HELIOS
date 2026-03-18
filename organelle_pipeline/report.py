from html import escape
from pathlib import Path

from organelle_pipeline.models import IsomerQuantResult
from organelle_pipeline.models import VariantCall


def write_html_report(
    path: Path,
    sample_name: str,
    isomer_result: IsomerQuantResult,
    heteroplasmy_calls: list[VariantCall],
) -> None:
    go = __import__("plotly.graph_objects", fromlist=["Figure"])
    plot = __import__("plotly.offline", fromlist=["plot"]).plot

    path.parent.mkdir(parents=True, exist_ok=True)

    assembly_only = isomer_result.method == "assembly_only_candidates"
    isomer_section_heading = "Isomer Candidate Summary" if assembly_only else "Isomer Proportion Plot"
    report_heading = "HELIOS Plastome Assembly Evaluation" if assembly_only else "HELIOS Heteroplasmy & Isomer Quantification"
    mode_note = (
        "<p><strong>Mode:</strong> Assembly-only evaluation. No FASTQ reads were supplied, so HELIOS reports candidate plastome isomers without read-backed proportions or heteroplasmy calls.</p>"
        if assembly_only
        else "<p><strong>Mode:</strong> Read-backed analysis with isomer quantification and optional heteroplasmy calling.</p>"
    )

    if assembly_only:
        isomer_div = (
            "<div class=\"empty\">"
            f"<p><strong>Candidate A:</strong> {escape(isomer_result.isomer_a_name)}</p>"
            f"<p><strong>Candidate B:</strong> {escape(isomer_result.isomer_b_name)}</p>"
            "<p>Read-backed isomer proportions are unavailable in assembly-only mode.</p>"
            "</div>"
        )
        metrics = [
            ("Evaluation Mode", "Assembly-only"),
            ("Candidate Isomers", "2"),
            ("Reads Processed", "N/A"),
            ("Heteroplasmy Calls", "N/A"),
        ]
    else:
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
        metrics = [
            ("Informative Isomer Reads", str(isomer_result.assigned_a + isomer_result.assigned_b)),
            ("Ambiguous Reads", str(isomer_result.ambiguous)),
            ("Unassigned Reads", str(isomer_result.unassigned)),
            ("Heteroplasmy Calls", str(len(heteroplasmy_calls))),
        ]

    variant_fig = go.Figure()
    if heteroplasmy_calls:
        x_values = [call.position for call in heteroplasmy_calls]
        y_values = [call.alt_fraction * 100.0 for call in heteroplasmy_calls]
        labels = [f"{call.ref}>{call.alt}" for call in heteroplasmy_calls]
        variant_fig.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="markers",
                text=labels,
                marker={"color": "#2ca02c", "size": 8},
            )
        )
    variant_fig.update_layout(
        title="Heteroplasmy Allele Fractions",
        xaxis_title="Genome position",
        yaxis_title="Alt allele fraction (%)",
    )

    variant_div = plot(variant_fig, include_plotlyjs=not assembly_only, output_type="div")
    metrics_html = "".join(
        f'<div class="tile"><div>{escape(label)}</div><div class="value">{escape(value)}</div></div>'
        for label, value in metrics
    )
    assumptions_html = "".join(f"<li>{escape(item)}</li>" for item in isomer_result.assumptions)
    variant_table = _variant_table(heteroplasmy_calls, assembly_only=assembly_only)

    html = f"""
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <title>HELIOS Report - {escape(sample_name)}</title>
    <style>
      body {{ font-family: "IBM Plex Sans", "Segoe UI", sans-serif; background: linear-gradient(180deg, #f4f8fb 0%, #e7eff7 100%); margin: 0; color: #0f172a; }}
      main {{ max-width: 1100px; margin: 0 auto; padding: 24px 16px 48px; }}
      h1 {{ margin-bottom: 8px; }}
      section {{ background: rgba(255, 255, 255, 0.85); border-radius: 14px; padding: 16px; margin-top: 16px; box-shadow: 0 8px 24px rgba(15, 23, 42, 0.08); }}
      table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
      th, td {{ border: 1px solid #cbd5e1; padding: 6px 8px; text-align: left; }}
      th {{ background: #f8fafc; }}
      .kpi {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 10px; }}
      .tile {{ background: #f8fafc; border: 1px solid #dbeafe; border-radius: 12px; padding: 10px; }}
      .tile .value {{ font-size: 1.25rem; font-weight: 700; }}
      .empty {{ border: 1px dashed #94a3b8; border-radius: 12px; padding: 16px; background: #f8fafc; }}
    </style>
  </head>
  <body>
    <main>
      <h1>{report_heading}</h1>
      <p><strong>Sample:</strong> {escape(sample_name)}</p>
      {mode_note}

      <section>
        <h2>Run Metrics</h2>
        <div class="kpi">{metrics_html}</div>
      </section>

      <section>
        <h2>{isomer_section_heading}</h2>
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


def _variant_table(calls: list[VariantCall], assembly_only: bool = False) -> str:
    if assembly_only:
        return "<p>Heteroplasmy calling is disabled in assembly-only mode.</p>"
    if not calls:
        return "<p>No heteroplasmy variants passed thresholds.</p>"
    rows = []
    for call in calls:
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
            "</tr>"
        )
    body = "".join(rows)
    return (
        "<table><thead><tr>"
        "<th>Contig</th><th>Position</th><th>Ref</th><th>Alt</th><th>Depth</th><th>Alt Count</th><th>Alt %</th><th>Type</th>"
        "</tr></thead><tbody>"
        f"{body}</tbody></table>"
    )
