#!/bin/bash
cd "$(dirname "$(readlink -f "$0")")"

GIT_TERMINAL_PROMPT=0 git fetch --quiet 2>/dev/null

BEHIND=$(git rev-list HEAD..@{u} --count 2>/dev/null)

if [ "$BEHIND" -gt 0 ] 2>/dev/null; then
    zenity --warning --text="TraceEditor is $BEHIND commit(s) behind the remote.\nConsider running: git pull" --title="Update available"
fi

/data/software/uv run trace_editor