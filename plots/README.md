# Plot snapshots

Reference renders of `notebook.ipynb` at the current point in the prototype.
We commit these so we can `git log -- plots/` to walk through how the design
evolved and `git show <sha>:plots/<file>` to recover an older look.

Whenever we iterate on visuals, re-render the notebook and overwrite these
files; git history keeps the previous versions.

| File | Source cell | What it shows |
|------|-------------|---------------|
| `01_basic_coverage.png` | cell 1 | Plain coverage track, no genes |
| `02_simple_genes.png`   | cell 2 | 3 non-overlapping single-exon genes, dotted boundary guides |
| `03_multi_exon.png`     | cell 3 | Single multi-exon gene with introns, chevrons, arrow tip |
| `04_overlapping_lanes.png` | cell 4 | Overlapping genes auto-stacked into lanes |
| `05_messy_mixed.png`    | cell 5 | Kitchen sink: overlapping + multi-exon + tiny + both strands |

The interactive viewer lives in `notebook_interactive.ipynb` and renders as a
live widget — no static snapshot is captured here.

## GenomeView snapshots (`notebook_genomeview.ipynb`)

Saved by `genomeview.save(...)` calls inline in that notebook. SVG (vector)
is GenomeView's native output; rasterize with `rsvg-convert` /  `inkscape`
if you need a PNG.

| File | What it shows |
|------|---------------|
| `genomeview_01_basic.svg`        | Coverage + gene track (counterpart to `02_simple_genes.png`) |
| `genomeview_02_stranded.svg`     | Stranded coverage via `GraphTrack.add_series` (forward blue, reverse red) |
| `genomeview_03_multi_region.svg` | Two non-contiguous regions in one document |

To regenerate after a code change:

```bash
uv run --with nbconvert jupyter nbconvert --to notebook --execute notebook.ipynb --output notebook.ipynb
uv run --with nbformat python scripts/extract_plots.py
```
