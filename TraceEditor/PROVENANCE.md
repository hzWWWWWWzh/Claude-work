# Provenance

Code-only snapshot of **TraceEditor**, the PySide6 rewrite of the magnetic-tweezers
trace analysis app.

| | |
|---|---|
| Upstream | `github.com/JvN2/TraceEditor` (private) |
| Branch | `TE-rewrite` |
| Commit | `0e8bdc4` — "use lambda dna length in pipepline", 2026-08-12 |
| Version | `0.5.1.dev29+g0e8bdc499` |
| Rewrite author | Barton Muller |
| Snapshot taken | 2026-10-15 |

Upstream git history is **not** included here — this is a flat snapshot of the
working tree at that commit. For history, `git log` against the upstream repo.

## What was left out

`data/` (71 MB of measurement `.h5`/`.tdms`), `hairpin/` (68 MB unzipping analysis
package) and `figures/` (generated PDFs) were excluded to keep this repository light.
`tests/data/real_data.h5` was kept, because without it two tests skip instead of
running.

## Running it

See [`INSTALL.md`](INSTALL.md). Needs Python >= 3.13. Verified 2026-08-27 on
macOS 15 / Python 3.13.15: 242/242 tests pass, GUI opens real `.h5` files, and the
full drift -> force -> zWLC pipeline runs clean.

Known gaps in this version: **File > Export**, **Parameters > Plot Histogram**,
**Process > Free DNA** and **Process > Ruptures** are status-bar stubs, and the old
Tkinter app's power-spectrum force calibration (**Force > Calibrate**) was never
ported — `force_spectroscopy.power_spectrum()` is an empty placeholder. Use
`oldhdf_versions/TraceEditor.py` for that.
