#!/bin/bash
# Build a self-contained arXiv submission bundle: paper/arxiv/arxiv_submission.tar.gz
# (arXiv does not run BibTeX, so the .bbl is included; figures are copied local.)
set -euo pipefail
cd "$(dirname "$0")"

# make sure the .bbl is current
latexmk -pdf -interaction=nonstopmode main >/dev/null

BUILD=arxiv
rm -rf "$BUILD" && mkdir -p "$BUILD"

# tex with a local graphicspath
sed 's|\\graphicspath{{../figures/}}|\\graphicspath{{./}}|' main.tex > "$BUILD/main.tex"
cp main.bbl "$BUILD/"

# copy exactly the figures the manuscript includes
grep -o 'includegraphics\[[^]]*\]{[^}]*}' main.tex | sed 's/.*{\(.*\)}/\1/' | sort -u | while read -r f; do
    cp "../figures/$f" "$BUILD/"
done

# verify the bundle compiles standalone
( cd "$BUILD" && pdflatex -interaction=nonstopmode main >/dev/null && pdflatex -interaction=nonstopmode main >/dev/null )
PAGES=$(python3 -c "import re;print(re.findall(r'Output written on main.pdf \((\d+) pages',open('$BUILD/main.log',encoding='latin-1').read())[0])")

( cd "$BUILD" && tar czf arxiv_submission.tar.gz main.tex main.bbl ./*.png )
echo "bundle: $BUILD/arxiv_submission.tar.gz  (standalone compile: ${PAGES} pages)"
tar tzf "$BUILD/arxiv_submission.tar.gz"
