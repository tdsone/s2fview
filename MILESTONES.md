# Design milestones

A running log of commits that represent design states worth coming back to.
Walk the list to find the look you want; pull it back with
`git checkout <hash>` or `git show <hash>:<file>`.

| Commit | Notes |
|--------|-------|
| [`8c46725`](https://github.com/tdsone/s2fview/commit/8c46725) | **"I really like the interactive version"** — interactive viewer with the zoom/pan toolbar, DNA sequence track (lazy letters on zoom), gene annotations with strand colors + arrow tips, all in Inter. |
| [`40c85e5`](https://github.com/tdsone/s2fview/commit/40c85e5) | **Pivot: Plotly becomes the primary interactive path.** Rendering speed and hover lag are much better in the Plotly viewer (`notebook_plotly.ipynb`), so future interactive work targets `s2fview.plotly_track`. The matplotlib + ipympl path stays for static figures and as a reference, but won't get new features by default. |
