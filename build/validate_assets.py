#!/usr/bin/env python3
"""Check that every asset app.py names is present on disk."""
import re, pathlib, sys
APP = pathlib.Path(__file__).resolve().parent.parent / "app.py"
D   = APP.parent / "data"
src = APP.read_text(encoding="utf8")
bad = []

def check(p, why):
    if not (D / p).exists():
        bad.append(f"{why}: data/{p}")

# figures resolved from the caption-dict keys
for m in re.finditer(r'^    "Figure (\d+) —', src, re.M):
    check(f"figures/main/Figure_{m.group(1)}.pdf", "main figure")
for m in re.finditer(r'^    "Figure S(\d+) —', src, re.M):
    check(f"figures/supplementary/Figure_S{m.group(1)}.pdf", "supp figure")
# tree section
for m in re.finditer(r'\("(Figure_S\d+\.pdf)",\s*(?:"([^"]+)"|None),\s*(?:"([^"]+)"|None),\s*(?:"([^"]+)"|None)\)', src):
    check(f"figures/supplementary/{m.group(1)}", "tree figure")
    if m.group(2): check(f"orthology/{m.group(2)}", "orthology")
    if m.group(3): check(f"trees/{m.group(3)}", "newick")
    if m.group(4): check(f"sequences/{m.group(4)}", "fasta")
# additional files
for m in re.finditer(r'^    \("(Additional_file_[^"]+)",', src, re.M):
    check(f"supplementary/{m.group(1)}", "additional file")
# top-level data files referenced by name
for m in re.finditer(r'"((?:BSF|results_|condition_|chemo_)[A-Za-z0-9_\-.]+\.csv)"', src):
    check(m.group(1), "data table")
# stats reports
for m in re.finditer(r'"stats_reports"\s*/\s*"([^"]+)"', src):
    check(f"stats_reports/{m.group(1)}", "stats report")

n_main = len(re.findall(r'^    "Figure \d+ —', src, re.M))
n_supp = len(re.findall(r'^    "Figure S\d+ —', src, re.M))
n_af   = len(re.findall(r'^    \("Additional_file_', src, re.M))
print(f"  {n_main} main figures, {n_supp} supplementary figures, {n_af} additional files referenced")
if bad:
    print(f"  {len(bad)} MISSING:")
    for b in bad: print("   -", b)
    sys.exit(1)
print("  all referenced assets present")
