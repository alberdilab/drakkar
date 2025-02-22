
import pandas as pd
import plotly.graph_objects as go
import argparse
import os

def generate_html_chunk(df, bar_height=28):
    """Generates an HTML chunk containing two interactive stacked barplots."""

    has_host = "bases_host" in df.columns  # Check if 'bases_host' exists

    # Sort samples in REVERSE alphabetical order
    df = df.sort_values(by="sample", ascending=False)

    # Extract values for each category
    sample_names = df["sample"]
    bases_discarded = df["bases_discarded"]
    bases_metagenomic = df["bases_metagenomic"]
    bases_host = df["bases_host"] if has_host else [0] * len(df)

    # Convert raw bases to gigabases (Gb)
    bases_discarded_gb = bases_discarded / 1e9
    bases_metagenomic_gb = bases_metagenomic / 1e9
    bases_host_gb = (bases_host / 1e9) if has_host else None

    # Compute total bases per sample
    total_bases = bases_discarded + bases_metagenomic + bases_host

    # Compute relative percentages
    pct_discarded = (bases_discarded / total_bases) * 100
    pct_metagenomic = (bases_metagenomic / total_bases) * 100
    pct_host = (bases_host / total_bases) * 100 if has_host else None

    # Dynamic chart height (~30% reduced)
    chart_height = max(len(df) * bar_height, 250)  # Minimum height 250px

    # Create figure for relative abundances
    fig_relative = go.Figure()
    fig_relative.add_trace(go.Bar(
        y=sample_names, x=pct_discarded, name="Discarded",
        orientation='h', marker=dict(color="red"),
        hoverinfo='x+name'
    ))

    if has_host:
        fig_relative.add_trace(go.Bar(
            y=sample_names, x=pct_host, name="Host",
            orientation='h', marker=dict(color="orange"),
            hoverinfo='x+name'
        ))

    fig_relative.add_trace(go.Bar(
        y=sample_names, x=pct_metagenomic, name="Metagenomic",
        orientation='h', marker=dict(color="green"),
        hoverinfo='x+name'
    ))

    fig_relative.update_layout(
        barmode='stack',
        showlegend=True,
        xaxis=dict(title="Percentage of Bases"),
        yaxis=dict(title="Samples"),
        legend=dict(orientation="h", yanchor="top", y=1.02, xanchor="center", x=0.5),
        template="plotly_white",
        width=600,
        height=chart_height
    )

    # Create figure for raw base counts in Gb
    fig_raw = go.Figure()
    fig_raw.add_trace(go.Bar(
        y=sample_names, x=bases_discarded_gb, name="Discarded",
        orientation='h', marker=dict(color="red"),
        hoverinfo='x+name'
    ))

    if has_host:
        fig_raw.add_trace(go.Bar(
            y=sample_names, x=bases_host_gb, name="Host",
            orientation='h', marker=dict(color="orange"),
            hoverinfo='x+name'
        ))

    fig_raw.add_trace(go.Bar(
        y=sample_names, x=bases_metagenomic_gb, name="Metagenomic",
        orientation='h', marker=dict(color="green"),
        hoverinfo='x+name'
    ))

    fig_raw.update_layout(
        barmode='stack',
        showlegend=True,
        xaxis=dict(title="Total Bases (Gb)"),
        yaxis=dict(title="Samples"),
        legend=dict(orientation="h", yanchor="top", y=1.02, xanchor="center", x=0.5),
        template="plotly_white",
        width=600,
        height=chart_height
    )

    # Convert plots to HTML divs
    html_chunk = f"""
    <h2>Preprocessing Statistics</h2>
    <div>
        <p>Relative and absolute fractions of discarded, host, and metagenomic DNA bases. Note that the metagenomic fraction might contain host DNA if no reference genome was used for removing host DNA or the employed genome was of a distantly related organism. The metagenomic fraction might also include non bacterial and archaeal DNA that will be assembled but won't contribute to the reconstruction of metagenome-assembled genomes (MAGs). In consequence, the below results should be interpreted cautiously.</p>

        <div class="row">
            <div class="column">
                <div class="plot-container">{fig_relative.to_html(full_html=False, include_plotlyjs="cdn")}</div>
            </div>
            <div class="column">
                <div class="plot-container">{fig_raw.to_html(full_html=False, include_plotlyjs="cdn")}</div>
            </div>
        </div>

        <style>
            .row {{
                display: flex;
                justify-content: space-between;
                align-items: center;
            }}
            .column {{
                width: 48%;
            }}
            .plot-container {{
                margin-bottom: 20px;
            }}
        </style>
    </div>
    """
    return html_chunk

def update_html_report(input_html, input_data, output_html):
    """Updates an existing HTML report by appending interactive stacked barplots."""

    # Read the existing HTML file
    with open(input_html, "r") as f:
        html_content = f.read()

    # Split into lines and remove the last three lines (</div>, </body>, </html>)
    html_lines = html_content.strip().split("\n")

    # Ensure the last line is </html>, otherwise raise an error
    if len(html_lines) < 3 or not html_lines[-1].strip().endswith("</html>"):
        raise ValueError("Invalid HTML file: missing required closing tags.")

    # Remove the last three lines: </div>, </body>, </html>
    html_content = "\n".join(html_lines[:-3])

    # Load the sequencing data
    file_ext = os.path.splitext(input_data)[-1].lower()
    df = pd.read_csv(input_data, sep="\t" if file_ext in [".tsv", ".txt"] else ",")

    # Ensure required columns exist
    required_cols = {"sample", "bases_discarded", "bases_metagenomic"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_cols}")

    # Generate the new HTML chunk
    html_chunk = generate_html_chunk(df)

    # Append the new section and re-add the removed lines
    updated_html = f"{html_content}\n<section>\n{html_chunk}\n</section>\n</div>\n</body>\n</html>"

    # Write the updated HTML file
    with open(output_html, "w") as f:
        f.write(updated_html)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Append an interactive base composition section to an HTML report.")
    parser.add_argument("-r", "--input_html", required=True, help="Path to the existing HTML report.")
    parser.add_argument("-i", "--input_data", required=True, help="Path to the input CSV/TSV file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the updated output HTML report.")

    args = parser.parse_args()

    update_html_report(args.input_html, args.input_data, args.output)
