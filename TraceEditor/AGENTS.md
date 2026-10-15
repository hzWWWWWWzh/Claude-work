# Agent guide

Start here, then read the doc that matches your task:

- **Design philosophy & contribution rules** → [`CONTRIBUTING.md`](CONTRIBUTING.md)
  (keep it simple, match the legacy app, UI-via-`.ui`-only, docstring style)
- **Architecture / internals** → [`docs/architecture.md`](docs/architecture.md)
  (UI layout, Z-correction logic, fit/batch/selection flow, shortcuts)
- **File format** → [`docs/file-format.md`](docs/file-format.md)
- **API reference (read these instead of whole source files)** →
  [`docs/api/`](docs/api/), generated from docstrings via
  `uv run pydoc-markdown`
- **User-facing overview & workflow** → [`README.md`](README.md)
- **Tests** → [`tests/README.md`](tests/README.md)

The legacy tkinter app in `oldhdf_versions/` is the canonical reference and must
not be modified — read the corresponding original before implementing a feature.

## Committing

When asked to commit, always use concise commit messages and do **not** add
Claude/Anthropic co-authorship trailers (no `Co-Authored-By: Claude ...`).
