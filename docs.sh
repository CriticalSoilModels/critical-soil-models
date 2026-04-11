#!/usr/bin/env bash
set -euo pipefail

# Build all documentation:
#   - FORD API reference  → build/doc/
#   - Typst theory PDFs   → docs/theory/*.pdf
#
# Requires: ford (activate the fpm conda env first), typst

# ford exits non-zero when warn=true finds undocumented items; treat as warnings only
ford fpm.toml || echo "ford: documentation warnings above (non-fatal)"

find docs -name '*.typ' | while read -r typ; do
    pdf="${typ%.typ}.pdf"
    echo "typst: $typ → $pdf"
    typst compile "$typ" "$pdf"
done
