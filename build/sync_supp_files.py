#!/usr/bin/env python3
"""Regenerate the _SUPP_FILES_ORDERED block in app.py from data/supplementary.

Large CSVs are stored gzipped so the repository stays small enough for a
reliable Streamlit Cloud cold start; the download button serves the .gz.
"""
import pathlib, re, sys

APP = pathlib.Path(__file__).resolve().parent.parent / "app.py"
DIR = APP.parent / "data" / "supplementary"

DESC = {
 1:"Appendage DE: antenna vs tarsi (DESeq2)",
 2:"Appendage DE: antenna vs maxillary palp (DESeq2)",
 3:"Appendage DE: tarsi vs maxillary palp (DESeq2)",
 4:"Antenna, sex: virgin female vs virgin male (DESeq2)",
 5:"Antenna, mating: mated female vs virgin female (DESeq2)",
 6:"Maxillary palp, sex: virgin female vs virgin male (DESeq2)",
 7:"Maxillary palp, mating: mated female vs virgin female (DESeq2)",
 8:"Tarsi, sex: virgin female vs virgin male (DESeq2)",
 9:"Tarsi, mating: mated female vs virgin female (DESeq2)",
10:"Normalised count matrix, all 27 libraries",
11:"GO annotations for all 24,828 genes",
12:"Amino-acid identity matrix — CSP",
13:"Amino-acid identity matrix — GR",
14:"Amino-acid identity matrix — IR",
15:"Amino-acid identity matrix — OBP",
16:"Amino-acid identity matrix — OR",
17:"Amino-acid identity matrix — PPK",
18:"Amino-acid identity matrix — TRP",
19:"Chemosensory family sizes across Diptera",
20:"Orthology assignments — OR",
21:"Orthology assignments — GR",
22:"Orthology assignments — IR",
23:"Orthology assignments — OBP",
24:"Orthology assignments — PPK",
25:"Orthology assignments — TRP",
26:"Orthology assignments — CSP",
27:"BLASTp hits for unknown-domain genes",
28:"Amino-acid identity of 1:1 orthologs",
29:"Antenna, mated sex: mated female vs virgin male (DESeq2)",
30:"Maxillary palp, mated sex: mated female vs virgin male (DESeq2)",
31:"Tarsi, mated sex: mated female vs virgin male (DESeq2)",
32:"Multivariate statistics: PCA, PERMANOVA, variation partitioning",
33:"GO over-representation — significant terms",
34:"GO over-representation — all tests",
35:"Identical-protein clusters",
36:"Newick trees for all seven families and the combined tree",
37:"Transcript model to gene locus map (567 models, 393 loci)",
38:"Statistical tests underlying Figure 3",
39:"FPKM matrix for all 27 libraries (derived from Additional file 10)",
}
MIME = {".csv":"text/csv", ".gz":"application/gzip", ".zip":"application/zip",
        ".txt":"text/plain"}

rows = []
for p in sorted(DIR.iterdir(), key=lambda q: int(re.match(r"Additional_file_(\d+)", q.name).group(1))):
    n = int(re.match(r"Additional_file_(\d+)", p.name).group(1))
    if n not in DESC:
        sys.exit(f"no description for Additional file {n}")
    note = " [gzip]" if p.suffix == ".gz" else (" [zip]" if p.suffix == ".zip" else "")
    mb = p.stat().st_size / 1e6
    size = f"{mb:.1f} MB" if mb >= 0.1 else f"{p.stat().st_size/1e3:.0f} kB"
    rows.append((p.name, f"Additional file {n} — {DESC[n]} ({size}){note}", MIME[p.suffix]))

out = ["_SUPP_FILES_ORDERED = ["]
for fn, desc, mime in rows:
    out.append(f'    ("{fn}",')
    out.append(f'     "{desc}", "{mime}"),')
out.append("]")

src = APP.read_text(encoding="utf8")
pat = re.compile(r"^_SUPP_FILES_ORDERED = \[.*?^\]", re.S | re.M)
if not pat.search(src):
    sys.exit("could not find _SUPP_FILES_ORDERED in app.py")
src = pat.sub("\n".join(out), src, count=1)
# keep the advertised range in step with what is actually bundled
hi = max(int(re.match(r"Additional_file_(\d+)", fn).group(1)) for fn, _, _ in rows)
src = re.sub(r"(Additional Files \(AF1\u2013AF)\d+(\))", lambda m: m.group(1) + str(hi) + m.group(2), src)
APP.write_text(src, encoding="utf8")
print(f"  {len(rows)} additional files listed; total "
      f"{sum(p.stat().st_size for p in DIR.iterdir())/1e6:.1f} MB")
