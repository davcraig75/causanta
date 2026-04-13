#!/usr/bin/env python3
"""
Render CAUSANTA documentation from Markdown to HTML.

Requires: pip install markdown

Usage:
    python scripts/render_docs.py           # Render all docs
    python scripts/render_docs.py tutorial  # Render specific doc
"""

import sys
from pathlib import Path

# Check for markdown library
try:
    import markdown
    from markdown.extensions.tables import TableExtension
    from markdown.extensions.fenced_code import FencedCodeExtension
    from markdown.extensions.toc import TocExtension
except ImportError:
    print("Error: markdown library not installed.")
    print("Install with: pip install markdown")
    sys.exit(1)


def get_html_template(title: str, css_path: str = "assets/style.css") -> str:
    """Return HTML template with linked CSS."""
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title} - CAUSANTA</title>
    <link rel="stylesheet" href="{css_path}">
</head>
<body>
    <div class="container">
        {{content}}
    </div>
</body>
</html>
'''


def render_markdown_file(md_path: Path, html_path: Path, css_path: str = "assets/style.css"):
    """Render a markdown file to HTML."""
    # Read markdown
    md_content = md_path.read_text(encoding="utf-8")

    # Extract title from first h1
    title = "CAUSANTA"
    for line in md_content.split("\n"):
        if line.startswith("# "):
            title = line[2:].strip()
            break

    # Configure markdown extensions
    md = markdown.Markdown(extensions=[
        TableExtension(),
        FencedCodeExtension(),
        TocExtension(permalink=False),
        'md_in_html',
    ])

    # Convert to HTML
    html_content = md.convert(md_content)

    # Wrap in template
    template = get_html_template(title, css_path)
    full_html = template.replace("{content}", html_content)

    # Write output
    html_path.write_text(full_html, encoding="utf-8")
    print(f"Rendered: {md_path.name} -> {html_path.name}")


def main():
    docs_dir = Path(__file__).parent.parent / "docs"

    # Documentation files to render
    doc_files = {
        "tutorial": ("tutorial.md", "tutorial.html"),
        "analysis_guide": ("analysis_guide.md", "analysis_guide.html"),
        "immune_system": ("immune_system.md", "immune_system.html"),
    }

    # Check which docs to render
    if len(sys.argv) > 1:
        targets = [arg.lower() for arg in sys.argv[1:]]
    else:
        targets = list(doc_files.keys())

    # Render each target
    for target in targets:
        if target not in doc_files:
            print(f"Unknown doc: {target}")
            print(f"Available: {', '.join(doc_files.keys())}")
            continue

        md_name, html_name = doc_files[target]
        md_path = docs_dir / md_name
        html_path = docs_dir / html_name

        if not md_path.exists():
            print(f"Warning: {md_path} not found, skipping")
            continue

        render_markdown_file(md_path, html_path)

    print("\nDone. HTML files use docs/assets/style.css")


if __name__ == "__main__":
    main()
