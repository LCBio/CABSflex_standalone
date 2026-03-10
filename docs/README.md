# HTML Documentation Site

This directory contains the site shell and local build instructions for the generated HTML documentation.

## Regenerate the site

Run from the repository root:

```bash
python scripts/build_html_wiki.py
```

The generator expects the legacy wiki repository to exist as a sibling checkout at:

```text
../CABSflex_standalone.wiki
```

Alternatively, set `CABSFLEX_WIKI_ROOT` to point at a different wiki checkout path before running the generator.

After generation, check `docs/migration-report.txt`. It lists legacy Bitbucket-linked resources that were referenced by the old wiki but were not present in the local wiki repo.

## Hosting recommendation

The site is generated in CI from the separate wiki repository and then published to Pages. Generated HTML and copied media are not kept in version control in the main repository.

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

The generator copies `videos/` into `docs/videos/` at build time and emits an HTML5 `<video>` player.

Generated HTML, images, uploads, videos, and the migration report are intentionally ignored in the main repository.

For future pages, you can either:

- use the same wiki-style thumbnail link format, or
- add raw HTML embed snippets directly into the source content before generation.
