# HTML Documentation Site

This directory contains a static HTML version of the legacy wiki.

## Regenerate the site

Run from the repository root:

```bash
python scripts/build_html_wiki.py
```

The generator expects the legacy wiki repository to exist as a sibling checkout at:

```text
../CABSflex_standalone.wiki
```

After generation, check `docs/migration-report.txt`. It lists legacy Bitbucket-linked resources that were referenced by the old wiki but were not present in the local wiki repo.

## Hosting recommendation

Because the new site is committed as plain HTML, CSS, JavaScript, and static assets, **GitHub Pages** or **GitLab Pages** is a better fit than Read the Docs.

Read the Docs is optimized for source-driven documentation builds. It can still host generated HTML indirectly, but that adds unnecessary build indirection for this repo.

## Video support

YouTube thumbnail links from the legacy wiki are converted into embedded `<iframe>` players during generation.

Local repo-hosted videos are also supported. Put them under `../CABSflex_standalone.wiki/videos/` and reference them in wiki pages using either a direct link:

```md
[Simulation movie](videos/example.mp4)
```

or a poster thumbnail link:

```md
[![Simulation movie](images/example-poster.jpg)](videos/example.mp4)
```

The generator copies `videos/` into `docs/videos/` and emits an HTML5 `<video>` player.

For future pages, you can either:

- use the same wiki-style thumbnail link format, or
- add raw HTML embed snippets directly into the source content before generation.
