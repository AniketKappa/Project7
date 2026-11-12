#!/usr/bin/env bash
# Render report/report.html (made by 13_make_report.py) to the final PDF with
# WeasyPrint, as in the pdf-design-skill. WeasyPrint lives in its own small
# conda environment ("pdfdesign": weasyprint + fonttools from conda-forge) and
# the soft style needs the Space Grotesk / IBM Plex Sans instances built by the
# skill's scripts/setup_fonts.sh (fonts in ~/.fonts).
#   bash scripts/14_render_report.sh
set -euo pipefail
TOP=$(cd "$(dirname "$0")/.." && pwd)
eval "$(conda shell.bash hook 2>/dev/null || conda shell.zsh hook)"
conda activate pdfdesign
cd "$TOP/report"
python -m weasyprint report.html "$TOP/Group4_FinalProjectReport.pdf"
echo "wrote $TOP/Group4_FinalProjectReport.pdf"
