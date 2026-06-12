#!/usr/bin/env python3
from __future__ import annotations

import html
import os
import re
import shutil
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import markdown


ROOT = Path(__file__).resolve().parents[1]
WIKI_ROOT = Path(
    os.environ.get("CABSFLEX_WIKI_ROOT", ROOT.parent / "CABSflex_standalone.wiki")
).resolve()
OUTPUT_ROOT = ROOT / "docs"
ASSETS_ROOT = OUTPUT_ROOT / "assets"
SITE_CSS = "assets/site.css"
SITE_JS = "assets/site.js"
VIDEO_EXTENSIONS = {".mp4", ".webm", ".ogg"}

PAGE_ORDER = [
    "Home",
    "Installation",
    "Modeling-Workflow",
    "CABS-Model",
    "Restraints",
    "Sampling-Temperature",
    "Flexibility-Modes",
    "All-Atom-Reconstruction",
    "Protein-Flexibility",
    "Peptide-Modeling",
    "Peptide-Protein-Docking",
    "Examples",
    "Gallery",
    "Options",
    "Visualization-Guide",
    "Generating-Report",
    "Published-Applications",
    "Advanced-Data",
    "Contact-and-Updates",
    "References",
    "Project-Links",
]

PAGE_LABELS = {
    "Home": "Overview",
    "Installation": "Installation",
    "Modeling-Workflow": "Modeling Workflow",
    "CABS-Model": "CABS Model",
    "Restraints": "Restraints",
    "Sampling-Temperature": "Sampling and Temperature",
    "Flexibility-Modes": "Flexibility Modes",
    "All-Atom-Reconstruction": "All-Atom Reconstruction",
    "Protein-Flexibility": "Protein Flexibility",
    "Peptide-Modeling": "Peptide Modeling",
    "Peptide-Protein-Docking": "Peptide–Protein Docking",
    "Examples": "Examples and Tutorials",
    "Gallery": "Gallery",
    "Published-Applications": "Published Applications",
    "Options": "Command-line Options",
    "Visualization-Guide": "Visualization Guide",
    "Generating-Report": "Generating Report",
    "Advanced-Data": "Internal Data Structures",
    "Contact-and-Updates": "Contact and Updates",
    "References": "References",
    "Project-Links": "Project Links",
}

SECTION_LABELS = {
    "Start": [
        "Home",
        "Installation",
    ],
    "Core Concepts": [
        "Modeling-Workflow",
        "CABS-Model",
        "Restraints",
        "Sampling-Temperature",
        "Flexibility-Modes",
        "All-Atom-Reconstruction",
    ],
    "User Guide": [
        "Protein-Flexibility",
        "Peptide-Modeling",
        "Peptide-Protein-Docking",
    ],
    "Resources": [
        "Examples",
        "Options",
        "Visualization-Guide",
        "Generating-Report",
        "Published-Applications",
        "Advanced-Data",
        "Contact-and-Updates",
        "References",
    ],
}

WIKI_LINK_RE = re.compile(r"\]\(([A-Za-z0-9_-]+)(#[^)]+)?\)")
ALERT_RE = re.compile(r"^(\s*)> \[!(NOTE|TIP|WARNING|IMPORTANT|CAUTION)\]\s*$", re.IGNORECASE)
YOUTUBE_THUMB_RE = re.compile(
    r'<a href="(?P<href>https?://(?:www\.)?(?:youtube\.com/watch\?v=[^"]+|youtu\.be/[^"?]+)[^"]*)">'
    r'\s*<img alt="(?P<alt>[^"]*)" src="(?P<src>[^"]*)"[^>]*>\s*</a>',
    re.IGNORECASE | re.DOTALL,
)
LOCAL_VIDEO_THUMB_RE = re.compile(
    r'<a href="(?P<href>(?:videos|media|assets/videos)/[^"]+\.(?:mp4|webm|ogg))">'
    r'\s*<img alt="(?P<alt>[^"]*)" src="(?P<src>[^"]*)"[^>]*>\s*</a>',
    re.IGNORECASE | re.DOTALL,
)
LOCAL_VIDEO_LINK_RE = re.compile(
    r'<p><a href="(?P<href>(?:videos|media|assets/videos)/[^"]+\.(?:mp4|webm|ogg))">'
    r'(?P<label>.*?)</a></p>',
    re.IGNORECASE | re.DOTALL,
)
BITBUCKET_IMAGE_RE = re.compile(r"https://bitbucket\.org/[^)\s]+")
UNRESOLVED_ASSETS: set[str] = set()


def page_output_name(page_name: str) -> str:
    return "index.html" if page_name == "Home" else f"{page_name.lower()}.html"


def youtube_embed_url(url: str) -> str | None:
    parsed = urlparse(url)
    host = parsed.netloc.lower()
    if "youtu.be" in host:
        video_id = parsed.path.strip("/")
    elif "youtube.com" in host:
        video_id = parse_qs(parsed.query).get("v", [None])[0]
    else:
        return None
    if not video_id:
        return None
    return (
        f"https://www.youtube.com/embed/{video_id}"
        f"?rel=0&modestbranding=1&playsinline=1&loop=1&playlist={video_id}"
    )


def build_sidebar(current_page: str, toc: str) -> str:
    parts = ['<nav class="sidebar-nav" aria-label="Documentation">']
    for section, pages in SECTION_LABELS.items():
        parts.append(f"<section><h2>{html.escape(section)}</h2><ul>")
        for page in pages:
            href = page_output_name(page)
            cls = ' class="is-current"' if page == current_page else ""
            label = PAGE_LABELS.get(page, page)
            parts.append(
                f'<li><a{cls} href="{href}">{html.escape(label)}</a>'
            )
            if page == current_page and toc and "toc-empty" not in toc:
                parts.append(f'<div class="sidebar-toc">{toc}</div>')
            parts.append('</li>')
        parts.append("</ul></section>")
    parts.append("</nav>")
    return "".join(parts)


def rewrite_internal_links(text: str) -> str:
    def repl(match: re.Match[str]) -> str:
        target = match.group(1)
        anchor = match.group(2) or ""
        if target in PAGE_LABELS:
            return f"]({page_output_name(target)}{anchor})"
        return match.group(0)

    return WIKI_LINK_RE.sub(repl, text)


def resolve_legacy_image(url: str) -> str:
    file_name = Path(urlparse(url).path).name
    candidates = [file_name, re.sub(r"^\d+-", "", file_name)]
    for candidate in candidates:
        if (WIKI_ROOT / "images" / candidate).exists():
            return f"images/{candidate}"
    UNRESOLVED_ASSETS.add(url)
    return url


def rewrite_bitbucket_images(text: str) -> str:
    return BITBUCKET_IMAGE_RE.sub(lambda match: resolve_legacy_image(match.group(0)), text)


def rewrite_github_alerts(text: str) -> str:
    lines = text.splitlines()
    rewritten: list[str] = []
    idx = 0
    while idx < len(lines):
        alert_match = ALERT_RE.match(lines[idx])
        if not alert_match:
            rewritten.append(lines[idx])
            idx += 1
            continue

        indent = alert_match.group(1)
        kind = alert_match.group(2).lower()
        idx += 1
        # Convert to admonition extension syntax: !!! kind
        # Add a blank line before the admonition if it's not at the start
        if rewritten and rewritten[-1].strip():
             rewritten.append("")
        rewritten.append(f"{indent}!!! {kind}")
        while idx < len(lines) and lines[idx].strip().startswith(">"):
            # Strip the "> " and preserve the rest of the line, then indent it for the admonition block
            line_content = re.sub(r"^\s*>\s?", "", lines[idx])
            rewritten.append(f"{indent}    {line_content}")
            idx += 1
        # Add an empty line to end the block
        rewritten.append("")
    return "\n".join(rewritten)


def preprocess_markdown(text: str) -> str:
    # Remove any existing manually defined header navigation line
    lines = text.splitlines()
    cleaned_lines = []
    for line in lines:
        if "hero-video-card" not in text and "Preprint" in line and ("GitHub" in line or "GitLab mirror" in line or "Code" in line):
            continue
        cleaned_lines.append(line)
    text = "\n".join(cleaned_lines)

    # Prepend the standardized header navigation at the top of the body
    header_nav = (
        '<span style="font-size: 1.15em;">'
        '[**<img src="https://img.icons8.com/material-outlined/24/000000/github.png" '
        'height="24" style="vertical-align: middle;"> GitHub**](Project-Links#GitHub) | '
        '[**<img src="https://img.icons8.com/color/24/gitlab.png" height="24" '
        'style="vertical-align: middle;"> GitLab mirror**](Project-Links#GitLab) | '
        '[**Preprint**](References#Preprint)</span>\n\n'
    )
    
    # If the page starts with the banner image, put the header nav right under the banner
    banner_match = re.match(r'^\s*(!\[.*?\]\(images/[^)]+\))(\s*\n+)?', text)
    if "hero-video-card" in text:
        # Custom hero block has its own layout and links, do not prepend standardized links
        pass
    elif banner_match:
        insert_pos = banner_match.end()
        text = text[:insert_pos] + header_nav + text[insert_pos:]
    else:
        text = header_nav + text

    text = rewrite_internal_links(text)
    text = rewrite_bitbucket_images(text)
    text = rewrite_github_alerts(text)
    
    # Support common emoji shortcodes (like GitLab/GitHub wiki)
    emoji_map = {
        ":arrow_up:": "⬆",
        ":arrow_down:": "⬇",
        ":arrow_left:": "⬅",
        ":arrow_right:": "➡",
    }
    for shortcode, glyph in emoji_map.items():
        text = text.replace(shortcode, glyph)
        
    return text


def convert_video_thumbnails(html_text: str) -> str:
    def repl(match: re.Match[str]) -> str:
        href = match.group("href")
        embed = youtube_embed_url(href)
        if not embed:
            return match.group(0)
        alt = html.escape(match.group("alt") or "Embedded video")
        return (
            '<figure class="video-embed">'
            f'<iframe src="{embed}" title="{alt}" loading="lazy" '
            'allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" '
            'allowfullscreen></iframe>'
            f"<figcaption>{alt}</figcaption>"
            "</figure>"
        )

    return YOUTUBE_THUMB_RE.sub(repl, html_text)


def convert_local_video_links(html_text: str) -> str:
    def thumb_repl(match: re.Match[str]) -> str:
        href = html.escape(match.group("href"))
        poster = html.escape(match.group("src"))
        # Adjust for new thumbnails subdirectory
        if poster.startswith("videos/") and not poster.startswith("videos/thumbnails/"):
            poster = poster.replace("videos/", "videos/thumbnails/", 1)
        alt = html.escape(match.group("alt") or "Embedded video")
        ext = Path(match.group("href")).suffix.lower()
        mime = f"video/{'ogg' if ext == '.ogg' else ext.lstrip('.')}"
        return (
            '<figure class="video-embed video-local">'
            f'<video autoplay muted loop preload="metadata" playsinline poster="{poster}">'
            f'<source src="{href}" type="{mime}">'
            "Your browser does not support the video tag."
            "</video>"
            f"<figcaption>{alt}</figcaption>"
            "</figure>"
        )

    def link_repl(match: re.Match[str]) -> str:
        href = html.escape(match.group("href"))
        label = re.sub(r"<.*?>", "", match.group("label")).strip() or "Embedded video"
        ext = Path(match.group("href")).suffix.lower()
        mime = f"video/{'ogg' if ext == '.ogg' else ext.lstrip('.')}"
        safe_label = html.escape(label)
        return (
            '<figure class="video-embed video-local">'
            f'<video autoplay muted loop preload="metadata" playsinline>'
            f'<source src="{href}" type="{mime}">'
            "Your browser does not support the video tag."
            "</video>"
            f"<figcaption>{safe_label}</figcaption>"
            "</figure>"
        )

    html_text = LOCAL_VIDEO_THUMB_RE.sub(thumb_repl, html_text)
    return LOCAL_VIDEO_LINK_RE.sub(link_repl, html_text)


def normalize_media_blocks(html_text: str) -> str:
    html_text = re.sub(
        r"<p>(<figure class=\"video-embed.*?</figure>)</p>",
        r"\1",
        html_text,
        flags=re.DOTALL,
    )
    html_text = re.sub(
        r'^<p><img alt="CABS-flex logo" src="([^"]+)"></p>',
        r'<figure class="page-banner"><a href="index.html"><img alt="CABS-flex logo" src="\1"></a></figure>',
        html_text,
        count=1,
        flags=re.MULTILINE,
    )
    return html_text


def add_external_link_targets(html_text: str) -> str:
    def repl(match: re.Match[str]) -> str:
        attrs = match.group(1)
        if re.search(r'href="https?://', attrs) and 'target=' not in attrs:
            return f'<a target="_blank" rel="noopener noreferrer" {attrs}>'
        return match.group(0)
    return re.sub(r'<a ([^>]+)>', repl, html_text)


def link_pdb_ids(html_text: str) -> str:
    parts = re.split(r'(<[^>]+>)', html_text)
    # Heuristic: starts with 1-9, 4 chars long, must contain at least one letter to avoid years.
    pdb_pattern = re.compile(r'\b([1-9][0-9A-Za-z]{3})(?:_([A-Z0-9]))?\b')
    for i in range(len(parts)):
        if i % 2 == 0:  # text node
            def pdb_repl(match: re.Match[str]) -> str:
                pdb_id = match.group(1)
                if not any(c.isalpha() for c in pdb_id):
                    return match.group(0)
                # Avoid common false positives
                if pdb_id.upper() in ('12GB', '8GB', 'RAM', 'MC-C', 'MC-S'):
                    return match.group(0)
                url = f"https://www.rcsb.org/structure/{pdb_id.upper()}"
                return f'<a href="{url}" target="_blank" rel="noopener noreferrer">{match.group(0)}</a>'
            parts[i] = pdb_pattern.sub(pdb_repl, parts[i])
    return "".join(parts)


def link_uniprot_ids(html_text: str) -> str:
    parts = re.split(r'(<[^>]+>)', html_text)
    # Pattern: UniProt: followed by ID
    uniprot_pattern = re.compile(r'UniProt:\s*([A-Z0-9]{6,10})')
    for i in range(len(parts)):
        if i % 2 == 0:  # text node
            def uniprot_repl(match: re.Match[str]) -> str:
                uniprot_id = match.group(1)
                url = f"https://www.uniprot.org/uniprotkb/{uniprot_id}/entry"
                return f'UniProt: <a href="{url}" target="_blank" rel="noopener noreferrer">{uniprot_id}</a>'
            parts[i] = uniprot_pattern.sub(uniprot_repl, parts[i])
    return "".join(parts)


def markdown_to_html(text: str) -> tuple[str, str, str]:
    md = markdown.Markdown(
        extensions=[
            "extra",
            "admonition",
            "attr_list",
            "toc",
        ],
        extension_configs={"toc": {"permalink": True}},
        output_format="html5",
    )
    body = md.convert(text)
    body = convert_video_thumbnails(body)
    body = convert_local_video_links(body)
    body = normalize_media_blocks(body)
    body = add_external_link_targets(body)
    body = link_pdb_ids(body)
    body = link_uniprot_ids(body)
    toc = md.toc or "<p class=\"toc-empty\">No table of contents for this page.</p>"
    
    # Strip the redundant top-level H1 link (page title) from the table of contents
    toc_match = re.search(
        r'<div class="toc">\s*<ul>\s*<li><a href="[^"]+">[^<]+</a>\s*<ul>\s*(.*?)\s*</ul>\s*</li>\s*</ul>\s*</div>',
        toc,
        re.DOTALL
    )
    if toc_match:
        toc = f'<div class="toc"><ul>{toc_match.group(1)}</ul></div>'

    title = "CABS-flex Documentation"
    match = re.search(r"<h1[^>]*>(.*?)</h1>", body, re.DOTALL)
    if match:
        title = html.unescape(re.sub(r"<.*?>", "", match.group(1))).replace("¶", "").strip()
    return title, toc, body


def page_template(page: str, title: str, toc: str, body: str) -> str:
    page_title = f"{title} | CABS-flex Docs"
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{html.escape(page_title)}</title>
    <link rel="stylesheet" href="{SITE_CSS}?v=2.4">
    <script src="{SITE_JS}?v=2.4"></script>
  </head>
  <body id="site-top">
    <div class="site-shell">
      <aside class="sidebar" id="sidebar">
        <a class="brand" href="index.html">CABS-flex Standalone 3 Docs</a>
        {build_sidebar(page, toc)}
      </aside>
      <main class="main-content">
        <header class="topbar">
          <a class="mobile-brand" href="index.html">CABS-flex Docs</a>
          <button class="nav-toggle" type="button" aria-expanded="false" aria-controls="sidebar">Menu</button>
        </header>
        <div class="content-grid">
          <article class="doc-content">
            {body}
          </article>
        </div>
      </main>
    </div>
  </body>
</html>
"""


def build_site() -> None:
    if not WIKI_ROOT.exists():
        raise SystemExit(f"Wiki repository not found at {WIKI_ROOT}")

    OUTPUT_ROOT.mkdir(exist_ok=True)
    (OUTPUT_ROOT / "images").mkdir(exist_ok=True)
    (OUTPUT_ROOT / "uploads").mkdir(exist_ok=True)
    ASSETS_ROOT.mkdir(exist_ok=True)

    for folder in ("images", "uploads", "videos"):
        src = WIKI_ROOT / folder
        dst = OUTPUT_ROOT / folder
        if dst.exists():
            shutil.rmtree(dst)
        if src.exists():
            shutil.copytree(src, dst)

    for page in PAGE_ORDER:
        src = WIKI_ROOT / f"{page}.md"
        if not src.exists():
            src_upper = WIKI_ROOT / f"{page}.MD"
            if src_upper.exists():
                src = src_upper
            else:
                raise SystemExit(f"Missing wiki page: {src}")
        text = preprocess_markdown(src.read_text(encoding="utf-8"))
        title, toc, body = markdown_to_html(text)
        out = OUTPUT_ROOT / page_output_name(page)
        out.write_text(page_template(page, title, toc, body), encoding="utf-8")

    report = OUTPUT_ROOT / "migration-report.txt"
    if UNRESOLVED_ASSETS:
        lines = [
            "Unresolved legacy Bitbucket resources",
            "=====================================",
            "",
            "These asset URLs were referenced by the legacy wiki but were not found in",
            "the local wiki repository checkout. They remain external in the generated HTML.",
            "",
        ]
        lines.extend(sorted(UNRESOLVED_ASSETS))
        report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    elif report.exists():
        report.unlink()


if __name__ == "__main__":
    build_site()
