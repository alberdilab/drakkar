
import argparse

def create_initial_html(output_html, project_name):
    """Generates the initial HTML report structure."""
    html_content = f"""<html>
    <head>
        <title>{project_name} - Sequencing Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; padding: 20px; background-color: #f8f9fa; }}
            h1, h2 {{ text-align: center; color: #333; }}
            .container {{ max-width: 1200px; margin: auto; background: white; padding: 20px; border-radius: 8px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1); padding-bottom: 20px; }}
            .content {{ padding: 20px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>DRAKKAR Report: {project_name}</h1>
            <div class="content">
                <p>This report contains a summary of the results of your DRAKKAR run.</p>
            </div>
        </div>
    </body>
    </html>"""

    with open(output_html, "w") as f:
        f.write(html_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an initial HTML report structure.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output HTML file.")
    parser.add_argument("-p", "--project", required=True, help="Project name to be used as the HTML title.")

    args = parser.parse_args()

    create_initial_html(args.output, args.project)
