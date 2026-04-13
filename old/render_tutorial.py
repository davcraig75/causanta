#!/usr/bin/env python3
"""Render tutorial.md to tutorial.html with GitHub-style formatting."""

import markdown
from pathlib import Path


def render_html(md_path: Path, html_path: Path) -> None:
    """Render markdown file to styled HTML using markdown library."""

    md_content = md_path.read_text()

    # Convert markdown to HTML with extensions
    md = markdown.Markdown(extensions=[
        'tables',
        'fenced_code',
        'toc',
        'nl2br',
    ])
    body_content = md.convert(md_content)

    # Wrap ASCII art blocks in proper formatting
    # The markdown library doesn't handle our ASCII boxes well, so we post-process
    import re

    # Find fenced code blocks that contain box-drawing characters
    def style_ascii_boxes(match):
        code = match.group(1)
        if '┌' in code or '├' in code or '│' in code:
            return f'<div class="diagram-box"><pre>{code}</pre></div>'
        return f'<pre><code>{code}</code></pre>'

    body_content = re.sub(
        r'<pre><code[^>]*>(.*?)</code></pre>',
        style_ascii_boxes,
        body_content,
        flags=re.DOTALL
    )

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CAUSANTA Tutorial</title>
    <style>
        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            line-height: 1.6;
            color: #24292f;
            background: #ffffff;
            margin: 0;
            padding: 0;
        }}

        .container {{
            max-width: 980px;
            margin: 0 auto;
            padding: 45px;
        }}

        h1 {{
            font-size: 2em;
            font-weight: 600;
            padding-bottom: 0.3em;
            border-bottom: 1px solid #d0d7de;
            margin-top: 24px;
            margin-bottom: 16px;
        }}

        h1:first-child {{
            margin-top: 0;
        }}

        h2 {{
            font-size: 1.5em;
            font-weight: 600;
            padding-bottom: 0.3em;
            border-bottom: 1px solid #d0d7de;
            margin-top: 24px;
            margin-bottom: 16px;
        }}

        h3 {{
            font-size: 1.25em;
            font-weight: 600;
            margin-top: 24px;
            margin-bottom: 16px;
        }}

        h4 {{
            font-size: 1em;
            font-weight: 600;
            margin-top: 24px;
            margin-bottom: 16px;
        }}

        p {{
            margin-bottom: 16px;
        }}

        a {{
            color: #0969da;
            text-decoration: none;
        }}

        a:hover {{
            text-decoration: underline;
        }}

        code {{
            font-family: ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, monospace;
            font-size: 85%;
            background: rgba(175, 184, 193, 0.2);
            padding: 0.2em 0.4em;
            border-radius: 6px;
        }}

        pre {{
            font-family: ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, monospace;
            font-size: 85%;
            background: #f6f8fa;
            padding: 16px;
            overflow: auto;
            border-radius: 6px;
            line-height: 1.45;
            margin: 16px 0;
            border: 1px solid #d0d7de;
        }}

        pre code {{
            background: transparent;
            padding: 0;
            font-size: 100%;
            border-radius: 0;
        }}

        .diagram-box {{
            margin: 24px 0;
            overflow-x: auto;
        }}

        .diagram-box pre {{
            background: #f6f8fa;
            border: 1px solid #d0d7de;
            font-size: 12px;
            line-height: 1.35;
            white-space: pre;
            overflow-x: auto;
            padding: 20px;
        }}

        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 16px 0;
            display: block;
            overflow-x: auto;
        }}

        th, td {{
            padding: 8px 16px;
            border: 1px solid #d0d7de;
            text-align: left;
        }}

        th {{
            background: #f6f8fa;
            font-weight: 600;
        }}

        tr:hover td {{
            background: #f6f8fa;
        }}

        ul, ol {{
            padding-left: 2em;
            margin: 16px 0;
        }}

        li {{
            margin: 4px 0;
        }}

        li > ul, li > ol {{
            margin: 4px 0;
        }}

        hr {{
            border: 0;
            height: 4px;
            background: #d0d7de;
            margin: 32px 0;
        }}

        blockquote {{
            padding: 0 1em;
            color: #57606a;
            border-left: 4px solid #d0d7de;
            margin: 16px 0;
        }}

        strong {{
            font-weight: 600;
        }}

        em {{
            font-style: italic;
        }}

        .footer {{
            margin-top: 48px;
            padding-top: 24px;
            border-top: 1px solid #d0d7de;
            text-align: center;
            color: #57606a;
            font-size: 14px;
        }}

        /* Navigation/TOC styling */
        .toc {{
            background: #f6f8fa;
            border: 1px solid #d0d7de;
            border-radius: 6px;
            padding: 16px 24px;
            margin: 24px 0;
        }}

        .toc ul {{
            list-style: none;
            padding-left: 0;
        }}

        .toc li {{
            padding: 4px 0;
        }}

        .toc a {{
            color: #0969da;
        }}

        /* Responsive */
        @media (max-width: 767px) {{
            .container {{
                padding: 15px;
            }}

            .diagram-box pre {{
                font-size: 10px;
                padding: 12px;
            }}

            table {{
                font-size: 14px;
            }}

            th, td {{
                padding: 6px 10px;
            }}
        }}

        /* Print styles */
        @media print {{
            body {{
                font-size: 11pt;
                color: black;
            }}

            .container {{
                max-width: none;
                padding: 0;
            }}

            pre, .diagram-box pre {{
                border: 1px solid #999;
                page-break-inside: avoid;
                font-size: 9pt;
            }}

            h1, h2, h3, h4 {{
                page-break-after: avoid;
            }}

            a {{
                color: black;
                text-decoration: underline;
            }}
        }}

        /* Highlight for table of contents */
        #table-of-contents + ol,
        #table-of-contents + ul {{
            background: #f6f8fa;
            border: 1px solid #d0d7de;
            border-radius: 6px;
            padding: 16px 16px 16px 40px;
        }}
    </style>
</head>
<body>
    <div class="container">
        {body_content}

        <div class="footer">
            <p><strong>CAUSANTA Tutorial</strong></p>
            <p>Causal Analysis Using Somatic And Neighborhood Tissue Architecture</p>
        </div>
    </div>
</body>
</html>'''

    html_path.write_text(html)
    print(f"Rendered: {html_path}")
    print(f"  Size: {len(html):,} bytes")


if __name__ == '__main__':
    docs_dir = Path(__file__).parent
    md_path = docs_dir / 'tutorial.md'
    html_path = docs_dir / 'tutorial.html'

    if not md_path.exists():
        print(f"Error: {md_path} not found")
        exit(1)

    render_html(md_path, html_path)
