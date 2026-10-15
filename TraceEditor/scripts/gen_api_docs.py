#!/usr/bin/env python
"""Generate agent-readable Markdown API docs from docstrings.

Renders one Markdown file per public module into docs/api/ using
pydoc-markdown, plus an index. Run after changing public docstrings:

    uv run python scripts/gen_api_docs.py

These files are for reading the API surface without opening whole source files.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
API_DIR = ROOT / "docs" / "api"

# Public modules worth documenting (import path -> output filename).
MODULES = {
    "trace_editor.core.h5_data": "h5_data.md",
    "trace_editor.core.legacy_import": "legacy_import.md",
    "trace_editor.core.inspect_h5": "inspect_h5.md",
    "trace_editor.core.force_spectroscopy": "force_spectroscopy.md",
    "trace_editor.core.batch_fit": "batch_fit.md",
    "trace_editor.core.force_calc": "force_calc.md",
    "trace_editor.core.auto_fit_pipeline": "auto_fit_pipeline.md",
    "trace_editor.core.file_browser": "file_browser.md",
    "trace_editor.core.file_query": "file_query.md",
    "trace_editor.gui.main_window": "main_window.md",
}


def render(module: str) -> str:
    """Run pydoc-markdown for one module and return its Markdown."""
    result = subprocess.run(
        ["pydoc-markdown", "-I", "src", "-m", module, "--render-toc"],
        cwd=ROOT, capture_output=True, text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        raise SystemExit(f"pydoc-markdown failed for {module}")
    return result.stdout


def main() -> None:
    API_DIR.mkdir(parents=True, exist_ok=True)
    index_lines = [
        "# API reference",
        "",
        "Generated from docstrings by `scripts/gen_api_docs.py` "
        "(`uv run python scripts/gen_api_docs.py`). Do not edit by hand.",
        "",
    ]
    for module, filename in MODULES.items():
        (API_DIR / filename).write_text(render(module))
        index_lines.append(f"- [`{module}`]({filename})")
        print(f"wrote docs/api/{filename}")

    index_lines.append("")
    (API_DIR / "README.md").write_text("\n".join(index_lines))
    print("wrote docs/api/README.md")


if __name__ == "__main__":
    main()
