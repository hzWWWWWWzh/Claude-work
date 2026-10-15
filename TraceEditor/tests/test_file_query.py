"""Tests for the Open File browser's query language (core/file_query.py)."""

import pytest

from trace_editor.core.file_query import parse_query, evaluate_query, QueryError

pytestmark = pytest.mark.io


ROW_SHIFT = {
    "axis": "Shift_mm", "tags": ["fail", "lowdensity"], "num_rois": 46,
    "selected": 15, "notes": "Fail — low density", "filename": "data_071.h5",
    "lut_filename": "data_013", "user": "wang", "date": "2026-06-05",
}
ROW_FOCUS = {
    "axis": "Focus_mm", "tags": [], "num_rois": 12, "selected": 12,
    "notes": "shift = 10.8 recorded before focus scan", "filename": "data_072.h5",
    "lut_filename": "", "user": "wang", "date": "2026-06-05",
}


class TestQueryEvaluation:

    @pytest.mark.parametrize("query,row,expected", [
        ("axis = shift", ROW_SHIFT, True),
        ("axis = shift", ROW_FOCUS, False),          # Notes containing "shift" must not match
        ("axis = focus", ROW_FOCUS, True),
        ("AXIS = SHIFT", ROW_SHIFT, True),            # case-insensitive keywords/values
        ("not axis = focus", ROW_SHIFT, True),
        ("not axis = focus", ROW_FOCUS, False),
        ("tag = fail", ROW_SHIFT, True),
        ("tag = fail", ROW_FOCUS, False),
        ("tag != fail", ROW_FOCUS, True),
        ("selected >= 10", ROW_SHIFT, True),
        ("selected >= 20", ROW_SHIFT, False),
        ("total >= 40", ROW_SHIFT, True),
        ("total >= 40", ROW_FOCUS, False),
        ("num_rois >= 40", ROW_SHIFT, True),          # alias
        ("lut = data_013", ROW_SHIFT, True),
        ("lut = data_013", ROW_FOCUS, False),
        ("axis = shift and not tag = fail", ROW_SHIFT, False),
        ("axis = shift and tag = fail", ROW_SHIFT, True),
        ("(axis = shift or axis = focus) and selected > 10", ROW_FOCUS, True),
        ("", ROW_SHIFT, True),                        # empty query matches everything
    ])
    def test_matches(self, query, row, expected):
        node = parse_query(query)
        assert evaluate_query(node, row) is expected

    def test_unknown_field_raises(self):
        node = parse_query("bogus_field = x")
        with pytest.raises(QueryError):
            evaluate_query(node, ROW_SHIFT)

    @pytest.mark.parametrize("bad_query", [
        "axis =",
        "(axis = shift",
        "axis shift",       # missing operator
        "and axis = shift",  # leading connective
    ])
    def test_malformed_query_raises_on_parse(self, bad_query):
        with pytest.raises(QueryError):
            parse_query(bad_query)

    def test_empty_query_returns_none_node(self):
        assert parse_query("") is None
        assert parse_query("   ") is None
