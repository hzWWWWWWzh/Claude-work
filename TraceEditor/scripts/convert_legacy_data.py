"""Convert every legacy `.hdf` file in `data/` to the modern `.h5` format.

If a `.h5` with the matching name already exists, it is deleted and
regenerated fresh from the `.hdf` (rather than trusting a possibly-stale
prior conversion). The `.hdf` is deleted only after its `.h5` is written
successfully.

Usage: uv run python scripts/convert_legacy_data.py
"""

from pathlib import Path

from trace_editor.core.legacy_import import import_legacy_hdf

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def main():
    hdf_paths = sorted(DATA_DIR.glob("*.hdf"))
    if not hdf_paths:
        print(f"No .hdf files found in {DATA_DIR}")
        return

    converted, failed = 0, []
    for hdf_path in hdf_paths:
        h5_path = hdf_path.with_suffix(".h5")
        if h5_path.exists():
            h5_path.unlink()

        try:
            import_legacy_hdf(hdf_path, h5_path)
        except Exception as e:
            print(f"FAILED  {hdf_path.name}: {e}")
            failed.append(hdf_path)
            continue

        hdf_path.unlink()
        print(f"OK      {hdf_path.name} -> {h5_path.name}")
        converted += 1

    print(f"\n{converted} converted, {len(failed)} failed out of {len(hdf_paths)}")


if __name__ == "__main__":
    main()
