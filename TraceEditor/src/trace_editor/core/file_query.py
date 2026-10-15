"""Small boolean query language for filtering the Open File browser table.

Grammar (case-insensitive keywords and field names):

    expr  := term (('and' | 'or') term)*
    term  := ['not'] atom
    atom  := field OP value | '(' expr ')'
    OP    := '=' | '==' | '!=' | '>=' | '<=' | '>' | '<'

Fields: axis, tag, notes, filename, selected, total, lut, user, date.
`tag` matches if the value is a substring of any of the row's tags.

For string fields, `=`/`==` is a case-insensitive substring match (so
``axis = shift`` matches ``Shift_mm``) and `!=` is its negation. `>`/`<`/`>=`/`<=`
on string fields compare alphabetically. For numeric fields (`selected`,
`total`), all operators compare numerically.

No `eval`/`ast` of arbitrary Python — this is a small hand-written
tokenizer/recursive-descent parser so query text can never execute code.
"""

import re


# Maps user-typed field names to keys in a `scan_h5_metadata` row dict.
# `tag` is handled specially (see `evaluate_query`) since a row's "tags" is a
# list, not a single value.
_FIELD_ALIASES = {
    "tags": "tag",
    "total": "num_rois",
    "lut": "lut_filename",
}

_TOKEN_RE = re.compile(r"""
    \s*(?:
        (?P<op>==|!=|>=|<=|=|>|<)
      | (?P<lparen>\()
      | (?P<rparen>\))
      | (?P<word>[^\s()=<>!]+)
    )
""", re.VERBOSE)


class QueryError(ValueError):
    """Raised for unparseable query text."""


def _tokenize(text: str) -> list:
    tokens = []
    pos = 0
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if not m or m.end() == pos:
            if text[pos:].strip() == "":
                break
            raise QueryError(f"Unexpected text: {text[pos:]!r}")
        pos = m.end()
        kind = m.lastgroup
        tokens.append((kind, m.group(kind)))
    return tokens


class _Parser:
    def __init__(self, tokens):
        self.tokens = tokens
        self.pos = 0

    def _peek(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else (None, None)

    def _next(self):
        tok = self._peek()
        self.pos += 1
        return tok

    def parse(self):
        node = self._expr()
        if self._peek() != (None, None):
            raise QueryError(f"Unexpected trailing text near {self._peek()[1]!r}")
        return node

    def _expr(self):
        node = self._term()
        while True:
            kind, word = self._peek()
            if kind == "word" and word.lower() in ("and", "or"):
                self._next()
                rhs = self._term()
                node = (word.lower(), node, rhs)
            else:
                break
        return node

    def _term(self):
        kind, word = self._peek()
        if kind == "word" and word.lower() == "not":
            self._next()
            return ("not", self._term())
        return self._atom()

    def _atom(self):
        kind, word = self._peek()
        if kind == "lparen":
            self._next()
            node = self._expr()
            kind2, _ = self._next()
            if kind2 != "rparen":
                raise QueryError("Missing closing ')'")
            return node
        if kind != "word":
            raise QueryError(f"Expected a field name, got {word!r}")
        field = word.lower()
        field = _FIELD_ALIASES.get(field, field)
        self._next()
        op_kind, op = self._peek()
        if op_kind != "op":
            raise QueryError(f"Expected an operator after {field!r}, got {op!r}")
        self._next()
        val_kind, value = self._next()
        if val_kind != "word":
            raise QueryError(f"Expected a value after {op!r}")
        return ("cmp", field, op, value)


def parse_query(text: str):
    """Parse query text into an evaluatable expression tree.

    :param text: query text, e.g. ``axis = shift and not tag = fail``.
    :return: an opaque expression tree for `evaluate_query`.
    :raises QueryError: on unparseable text.
    """
    tokens = _tokenize(text)
    if not tokens:
        return None
    return _Parser(tokens).parse()


def _compare(op: str, actual, value: str) -> bool:
    if isinstance(actual, (int, float)):
        try:
            value_num = float(value)
        except ValueError:
            return False
        if op in ("=", "=="):
            return actual == value_num
        if op == "!=":
            return actual != value_num
        if op == ">=":
            return actual >= value_num
        if op == "<=":
            return actual <= value_num
        if op == ">":
            return actual > value_num
        if op == "<":
            return actual < value_num
        return False

    actual_str = str(actual).lower()
    value_str = value.lower()
    if op in ("=", "=="):
        return value_str in actual_str
    if op == "!=":
        return value_str not in actual_str
    if op == ">=":
        return actual_str >= value_str
    if op == "<=":
        return actual_str <= value_str
    if op == ">":
        return actual_str > value_str
    if op == "<":
        return actual_str < value_str
    return False


def evaluate_query(node, row: dict) -> bool:
    """Evaluate a parsed query tree against one `scan_h5_metadata` row dict.

    :param node: result of `parse_query`, or `None` (matches everything).
    :param row: a row dict as returned by `core.file_browser.scan_h5_metadata`.
    :return: whether `row` matches the query.
    :raises QueryError: if the query references an unknown field.
    """
    if node is None:
        return True
    kind = node[0]
    if kind == "and":
        return evaluate_query(node[1], row) and evaluate_query(node[2], row)
    if kind == "or":
        return evaluate_query(node[1], row) or evaluate_query(node[2], row)
    if kind == "not":
        return not evaluate_query(node[1], row)
    if kind == "cmp":
        _, field, op, value = node
        if field == "tag":
            tags = row.get("tags", [])
            matches_any = any(_compare("=", tag, value) for tag in tags)
            return not matches_any if op == "!=" else matches_any
        if field not in row:
            raise QueryError(f"Unknown field: {field!r}")
        return _compare(op, row[field], value)
    raise QueryError(f"Unknown node: {node!r}")
