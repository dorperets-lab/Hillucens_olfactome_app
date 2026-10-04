#!/usr/bin/env python3
"""Regenerate the figure-caption dictionaries in app.py from the manuscript.

The manuscript is the single source of truth.  Run this after any caption edit
in the .docx and the app can never drift from the paper:

    python3 build/sync_captions.py

It rewrites the MAIN_FIGURE_CAPTIONS and SUPP_FIGURE_CAPTIONS blocks in place.
Only surviving text is taken: runs marked struck-through or red (C00000) are
the manual deletion marks, so they are skipped.
"""
import zipfile, re, html, pathlib, sys

DOCX = pathlib.Path("/mnt/c/Users/dorpe/OneDrive/Desktop/PhD_Projects/BSF/"
                    "Updated_GDrive_version/Submission_v5_2026-09/01_Manuscript/"
                    "Perets_et_al_2026_v5_marked.docx")
APP  = pathlib.Path(__file__).resolve().parent.parent / "app.py"

def final_text(p):
    out = []
    for run in re.findall(r"<w:r[ >].*?</w:r>", p, re.S):
        rpr = re.search(r"<w:rPr>.*?</w:rPr>", run, re.S)
        if rpr and ("<w:strike" in rpr.group(0) or "C00000" in rpr.group(0)):
            continue
        out += re.findall(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", run, re.S)
    return html.unescape(re.sub(r"\s+", " ", "".join(out))).strip()

xml   = zipfile.ZipFile(DOCX).read("word/document.xml").decode("utf8")
caps  = {}
for p in re.findall(r"<w:p[ >].*?</w:p>", xml, re.S):
    t = final_text(p)
    m = re.match(r"^Figure\s*(S?\d+)\.\s*(.*)$", t)
    if m and len(t) > 40:
        caps.setdefault(m.group(1), t)

# Short labels for the radio buttons.  The caption body always comes from the
# manuscript; only these menu labels are curated, so they stay readable.
SHORT = {
    "1": "Phylogeny of chemosensory families",
    "2": "Transcriptome overview",
    "3": "Appendage-specific divergence",
    "4": "Sex & mating remodelling",
    "5": "Hierarchical model",
    "S1": "Tree: all families", "S2": "Tree: CSP", "S3": "Tree: GR",
    "S4": "Tree: IR", "S5": "Tree: OBP", "S6": "Tree: OR",
    "S7": "Tree: PPK", "S8": "Tree: TRP",
    "S9":  "GO: antenna-biased", "S10": "GO: palp-biased",
    "S11": "GO: tarsi-biased",   "S12": "GO domain composition",
    "S13": "Heatmap: OR", "S14": "Heatmap: GR", "S15": "Heatmap: IR",
    "S16": "Heatmap: OBP", "S17": "Heatmap: PPK", "S18": "Heatmap: CSP",
    "S19": "Heatmap: TRP",
    "S20": "GO: sex-biased", "S21": "GO: mating-responsive",
    "S22": "GO: sex-biased only", "S23": "GO: mating-responsive only",
    "S24": "GO: sex & mating", "S25": "DE across appendages",
    "S26": "Mated-female-specific",
}

def emit(name, keys):
    lines = [f"{name} = {{"]
    for k in keys:
        body  = caps[k]
        label = f"Figure {k} — {SHORT[k]}"
        md = re.sub(r"^(Figure\s*S?\d+\.)", r"**\1**", body)
        md = md.replace('"', "'").replace("\\", "/")
        lines.append(f'    "{label}": (')
        for chunk in re.findall(r".{1,96}(?:\s|$)", md):
            lines.append(f'        "{chunk.strip()} "')
        lines.append("    ),")
    lines.append("}")
    return "\n".join(lines)

main_keys = [k for k in caps if not k.startswith("S")]
supp_keys = [k for k in caps if k.startswith("S")]
main_keys.sort(key=int)
supp_keys.sort(key=lambda s: int(s[1:]))
missing = [k for k in [str(i) for i in range(1, 6)] + [f"S{i}" for i in range(1, 27)] if k not in caps]
if missing:
    sys.exit(f"captions missing from the manuscript: {missing}")

src = APP.read_text(encoding="utf8")
for name, keys in (("MAIN_FIGURE_CAPTIONS", main_keys), ("SUPP_FIGURE_CAPTIONS", supp_keys)):
    pat = re.compile(rf"^{name} = \{{.*?^\}}", re.S | re.M)
    if not pat.search(src):
        sys.exit(f"could not find the {name} block in app.py")
    src = pat.sub(lambda _m, n=name, k=keys: emit(n, k), src, count=1)
APP.write_text(src, encoding="utf8")
print(f"  synced {len(main_keys)} main + {len(supp_keys)} supplementary captions from")
print(f"  {DOCX.name}")
for k in main_keys + supp_keys:
    print(f"    Figure {k:>3s} — {SHORT[k]}")
