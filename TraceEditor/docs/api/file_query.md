# Table of Contents

* [trace\_editor.core.file\_query](#trace_editor.core.file_query)
  * [QueryError](#trace_editor.core.file_query.QueryError)
  * [parse\_query](#trace_editor.core.file_query.parse_query)
  * [evaluate\_query](#trace_editor.core.file_query.evaluate_query)

<a id="trace_editor.core.file_query"></a>

# trace\_editor.core.file\_query

Small boolean query language for filtering the Open File browser table.

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

<a id="trace_editor.core.file_query.QueryError"></a>

## QueryError Objects

```python
class QueryError(ValueError)
```

Raised for unparseable query text.

<a id="trace_editor.core.file_query.parse_query"></a>

#### parse\_query

```python
def parse_query(text: str)
```

Parse query text into an evaluatable expression tree.

**Arguments**:

- `text`: query text, e.g. ``axis = shift and not tag = fail``.

**Raises**:

- `QueryError`: on unparseable text.

**Returns**:

an opaque expression tree for `evaluate_query`.

<a id="trace_editor.core.file_query.evaluate_query"></a>

#### evaluate\_query

```python
def evaluate_query(node, row: dict) -> bool
```

Evaluate a parsed query tree against one `scan_h5_metadata` row dict.

**Arguments**:

- `node`: result of `parse_query`, or `None` (matches everything).
- `row`: a row dict as returned by `core.file_browser.scan_h5_metadata`.

**Raises**:

- `QueryError`: if the query references an unknown field.

**Returns**:

whether `row` matches the query.

