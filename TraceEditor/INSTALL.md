# TraceEditor — install & run

Snapshot of branch `TE-rewrite` at commit **`0e8bdc4`** ("use lambda dna length
in pipepline", 2026-08-12) — the newest commit on
`github.com/JvN2/TraceEditor` as of 2026-08-27. Version string
`0.5.1.dev29+g0e8bdc499`.

Source tree only: no `.git` history, no `.venv`. A prebuilt wheel is included
in `dist/` so you do not need either.

## Requirements

- **Python ≥ 3.13** (verified on 3.13.15). 3.12 and older will not work.
- ~1.5 GB free disk for the virtual environment (PySide6 is most of it).
- No compiler needed — every dependency ships prebuilt wheels for
  macOS arm64/x86-64, Windows x86-64 and Linux x86-64.

## Install (pip)

```bash
cd TraceEditor_BM_0e8bdc4_20260827

python3.13 -m venv .venv
source .venv/bin/activate                  # Windows: .venv\Scripts\activate

pip install -r requirements.txt            # dependencies
pip install --no-deps dist/*.whl           # the app itself
                                           # Windows: pip install --no-deps (Get-Item dist\*.whl).FullName

trace_editor                               # launch the GUI
```

`requirements.txt` — the 16 direct dependencies, pinned to the versions that
were installed and verified on 2026-08-27.
`requirements-lock.txt` — all 73 packages of that same environment, if you
want a byte-identical rebuild (`pip install -r requirements-lock.txt`).

## Install (uv — the project's own recommendation, much faster)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh          # macOS / Linux, once
# powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"   # Windows

cd TraceEditor_BM_0e8bdc4_20260827
uv venv --python 3.13
uv pip install -r requirements.txt
uv pip install --no-deps dist/*.whl
uv run --no-project trace_editor
```

> **Why not plain `uv sync` / `pip install -e .`?**
> `pyproject.toml` derives the version from git tags via `hatch-vcs`, and this
> archive has no `.git`, so a from-source build fails with
> `LookupError: Error getting the version from source 'vcs'`. The bundled wheel
> sidesteps that. If you do want an editable dev install from this snapshot,
> set the version by hand first:
> ```bash
> SETUPTOOLS_SCM_PRETEND_VERSION_FOR_TRACE_EDITOR=0.5.1.dev29 pip install -e .
> ```
> Or `git clone https://github.com/JvN2/TraceEditor.git` (private repo — needs
> access) and run `uv sync` there, which is the normal developer workflow.

## Verify the install

```bash
pip install pytest pytest-qt pytest-cov     # already in requirements.txt
pytest
```

Expect **242 passed** in well under a minute. GUI tests render offscreen, so
no window appears and no display is required.

## Note for OneDrive / Dropbox users

Do **not** put `.venv` inside a synced folder. It is ~500 MB, it will be
uploaded on every sync, and on-demand placeholder files silently corrupt the
environment's binaries. Keep the venv outside the synced tree — e.g.
`python3.13 -m venv ~/.venvs/trace_editor`.

## Verified on 2026-08-27

macOS 15 (arm64), Python 3.13.15, PySide6 6.11.1:

- dependency install clean, all imports fine
- 242/242 tests pass (7.4 s)
- GUI opens real `.h5` files (85 ROIs), plots, and cycles ROIs without error
- full `zwlc_pipeline` (ExponentialDrift → calculate_force → zWLC) runs over
  54 ROIs in 0.5 s with 0 errors and writes 30 parameters per ROI
- legacy `.hdf` → `.h5` import works
- batch fit, folder scan, query language, plot presets, palettes,
  copy-to-clipboard all work

### Known gaps in this version

Four menu entries are **stubs** that only print to the status bar:
**File > Export**, **Parameters > Plot Histogram**, **Process > Free DNA**,
**Process > Ruptures**. Free DNA's maths does exist in
`core/force_spectroscopy.py` (`calc_free_DNA`, `get_free_DNA`) — it is just not
wired to the GUI.

There is **no power-spectrum force calibration** in this version — the old
Tkinter app's *Force > Calibrate* was not ported, and
`force_spectroscopy.power_spectrum()` is an empty placeholder. For that, use
`oldhdf_versions/TraceEditor.py`, which needs the extra packages listed at the
bottom of `requirements.txt`.

An `.h5` written in the older lowercase `parameters`/`traces` layout (no root
attributes) is accepted by **File > Open** without an error, but then makes
*Select > All*, *Select > None* and *Process > Remove All Processing* fail
silently with an uncaught `KeyError: 'num_rois'` in the console — for every
open file, not just the bad one. Close that file to recover.
