"""Tests for core/inspect_h5.py (Qt-free `code` CLI dispatch)."""

from unittest.mock import patch, MagicMock

import pytest

from trace_editor.core.inspect_h5 import open_in_vscode_h5web, InspectH5Error

pytestmark = pytest.mark.io


class TestOpenInVscodeH5web:
    def test_no_code_cli_raises(self):
        with patch("shutil.which", return_value=None):
            with pytest.raises(InspectH5Error) as exc_info:
                open_in_vscode_h5web("/tmp/foo.h5")
        assert exc_info.value.reason == "no_code_cli"

    def test_extension_missing_raises(self):
        with patch("shutil.which", return_value="/usr/bin/code"), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="some.other.extension\n")
            with pytest.raises(InspectH5Error) as exc_info:
                open_in_vscode_h5web("/tmp/foo.h5")
        assert exc_info.value.reason == "no_extension"

    def test_opens_file_when_extension_present(self):
        with patch("shutil.which", return_value="/usr/bin/code"), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="h5web.vscode-h5web\n")
            open_in_vscode_h5web("/tmp/foo.h5")
        # First call lists extensions, second opens the file.
        assert mock_run.call_count == 2
        open_call = mock_run.call_args_list[1]
        assert open_call.args[0] == ["/usr/bin/code", "-r", "/tmp/foo.h5"]

    def test_extension_check_is_case_insensitive(self):
        with patch("shutil.which", return_value="/usr/bin/code"), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="H5Web.VSCode-H5Web\n")
            open_in_vscode_h5web("/tmp/foo.h5")
        assert mock_run.call_count == 2
