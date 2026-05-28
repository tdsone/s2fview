# Design milestones

A running log of commits that represent design states worth coming back to.
Walk the list to find the look you want; pull it back with
`git checkout <hash>` or `git show <hash>:<file>`.

| Commit | Notes |
|--------|-------|
| [`8c46725`](https://github.com/tdsone/s2fview/commit/8c46725) | **"I really like the interactive version"** — historical: the matplotlib + ipympl viewer (now removed) when it had a zoom/pan toolbar, lazy DNA letters, strand-colored gene arrows, and Inter bundled with the package. Recover via `git show 8c46725:src/s2fview/__init__.py` etc. |
| [`40c85e5`](https://github.com/tdsone/s2fview/commit/40c85e5) | **Pivot to Plotly.** Hover lag and rendering speed were significantly better in the Plotly viewer, so future interactive work targeted `s2fview.plotly_track`. The matplotlib + ipympl, GenomeView, and igv-notebook paths kept as reference. |
| _post-consolidation_ | **Repo simplified to Plotly only.** Matplotlib, ipympl, GenomeView, and igv-notebook viewers + their notebooks were removed; `s2fview.plotly_track.plotly_coverage_track` was promoted into `s2fview.coverage_track`. Use `git log` from this point forward for the single-implementation history. |
