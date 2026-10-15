"""Scan a folder of `.h5` files and extract cheap preview metadata for the
Open File browser dialog, without loading trace data.
"""

import re
from pathlib import Path

from .h5_data import H5Data


_HASHTAG_RE = re.compile(r"#(\w+)")


def extract_tags(notes: str) -> list:
    """Return all `#tag` tokens found in `notes`, lowercased, in order of
    first appearance, with duplicates removed.

    Free-form: any `#word` in Notes becomes a tag — there is no fixed
    vocabulary. E.g. `"Fail — #lowdensity #retry"` yields
    `["lowdensity", "retry"]`.

    :param notes: free-text Notes attribute.
    :return: list of tag strings, without the leading `#`.
    """
    seen = []
    for match in _HASHTAG_RE.finditer(notes):
        tag = match.group(1).lower()
        if tag not in seen:
            seen.append(tag)
    return seen


def scan_h5_metadata(path) -> dict:
    """Read root attributes and selection counts for one `.h5` file.

    Cheap: only reads root attrs and the small `/Parameters` array, never
    trace data.

    :param path: path to a `.h5` file.
    :return: dict with keys `path`, `filename`, `axis`, `tags` (list),
        `selected`, `num_rois`, `notes`, `date`, `user`, `lut_filename`, and
        `error` (empty on success).
    """
    path = Path(path)
    row = {
        "path": str(path),
        "filename": path.name,
        "axis": "",
        "tags": [],
        "selected": 0,
        "num_rois": 0,
        "notes": "",
        "date": "",
        "user": "",
        "lut_filename": "",
        "error": "",
    }
    try:
        data = H5Data(path)
        notes = data.get_notes()
        row["axis"] = data.get_axis()
        row["notes"] = notes
        row["tags"] = extract_tags(notes)
        row["date"] = data.get_date()
        row["user"] = data.get_user()
        row["num_rois"] = data.num_rois()
        row["selected"] = data.selected_count()
        row["lut_filename"] = data.get_lut_filename()
    except Exception as e:
        row["error"] = str(e)
    return row


def scan_folder(folder) -> list:
    """Scan all `.h5` files directly inside `folder` (non-recursive).

    :param folder: directory to scan.
    :return: list of metadata dicts from `scan_h5_metadata`, sorted by filename.
    """
    folder = Path(folder)
    paths = sorted(folder.glob("*.h5"))
    return [scan_h5_metadata(p) for p in paths]
