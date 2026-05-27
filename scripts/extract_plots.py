"""Extract the rendered figures from notebook.ipynb into ./plots/.

Usage:
    uv run --with nbformat python scripts/extract_plots.py
"""

from __future__ import annotations

import base64
import pathlib
import sys

import nbformat

# Map of rendered-figure index (in order of appearance in the notebook) to the
# semantic filename we save it under. Keep in sync with plots/README.md.
NAMES: list[str] = [
    "01_basic_coverage.png",
    "02_simple_genes.png",
    "03_multi_exon.png",
    "04_overlapping_lanes.png",
    "05_messy_mixed.png",
    "06_interactive_static_fallback.png",
]


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parent.parent
    notebook_path = repo / "notebook.ipynb"
    out_dir = repo / "plots"
    out_dir.mkdir(exist_ok=True)

    nb = nbformat.read(notebook_path, as_version=4)
    figure_idx = 0
    for cell in nb.cells:
        for output in cell.get("outputs", []):
            data = output.get("data", {})
            if "image/png" not in data:
                continue
            if figure_idx >= len(NAMES):
                print(
                    f"warning: more figures ({figure_idx + 1}+) than NAMES entries; "
                    "add a label in scripts/extract_plots.py",
                    file=sys.stderr,
                )
                return 1
            target = out_dir / NAMES[figure_idx]
            target.write_bytes(base64.b64decode(data["image/png"]))
            print(target.relative_to(repo))
            figure_idx += 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
