"""
Open an `.h5` file in the H5Web VS Code extension, for inspecting the raw
HDF5 structure (all groups/datasets/attributes, not just what TraceEditor's
UI exposes).

Qt-free so it's testable without a QApplication; `gui/main_window.py` wraps
this with the dialog shown on failure.
"""

import shutil
import subprocess

VSCODE_H5WEB_EXTENSION_ID = "h5web.vscode-h5web"
MYHDF5_FALLBACK_URL = "https://myhdf5.hdfgroup.org"


class InspectH5Error(Exception):
    """Raised when the file can't be opened in VS Code's H5Web viewer.

    :ivar reason: `"no_code_cli"` or `"no_extension"`, for callers that want
        to tailor their message/fallback UI.
    """

    def __init__(self, message: str, reason: str):
        super().__init__(message)
        self.reason = reason


def open_in_vscode_h5web(filepath: str):
    """Open `filepath` in VS Code's H5Web viewer via the `code` CLI.

    :param filepath: path to the `.h5` file to inspect.
    :raises InspectH5Error: if the `code` CLI isn't on PATH, or the H5Web
        extension isn't installed.
    """
    code_cli = shutil.which("code")
    if code_cli is None:
        raise InspectH5Error(
            "The 'code' command was not found on PATH.", reason="no_code_cli"
        )

    result = subprocess.run(
        [code_cli, "--list-extensions"], capture_output=True, text=True, timeout=10
    )
    installed = result.stdout.lower().splitlines()
    if VSCODE_H5WEB_EXTENSION_ID not in installed:
        raise InspectH5Error(
            "The H5Web VS Code extension is not installed.", reason="no_extension"
        )

    subprocess.run([code_cli, "-r", filepath], timeout=10)
