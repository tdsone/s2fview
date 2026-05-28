# Plot snapshots

Static snapshots of `notebook.ipynb` saved via `plotly.graph_objects.Figure.write_image`
(which uses `kaleido` under the hood). The notebook renders the live
interactive widget; these PNGs are the "what it looks like" reference for
PRs and git history.

| File | What it shows |
|------|---------------|
| `messy_mixed.png`     | Full-range view: coverage + colored sequence strip + gene track. Letters hidden (no pixel room). |
| `zoomed_letters.png`  | Zoomed-in view: same figure with the letter trace forced visible to demonstrate the on-zoom behavior. |

Regenerate by re-running `notebook.ipynb` end-to-end in a kernel with
`kaleido` available, then calling
`fig.write_image('plots/messy_mixed.png', width=1200, height=600, scale=2)`.
