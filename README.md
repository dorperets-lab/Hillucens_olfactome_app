# *Hermetia illucens* olfactome browser

Interactive browser for the chemosensory transcriptomes of the antenna, maxillary
palp and tarsi of the black soldier fly, *Hermetia illucens*.

Companion to **Perets *et al.* 2026**. Every table the app reads is one of the
Additional files published with the paper, so the figures you build here are the
figures in the manuscript.

**Live app:** <https://hillucensolfactomeapp-rk3chuhczmmvbpfq2ctvan.streamlit.app/>

A phone screen is too narrow for the heatmaps and trees. Open **Open on a
computer** in the sidebar to mail yourself the link.

## What is in it

| Section | Contents |
|---|---|
| Phylogenetic Trees | ML trees for all seven families and the combined tree, with orthology assignments and sequences |
| Transcriptome Overview | PCA, correlation structure and variation partitioning across 27 libraries |
| Appendage Comparison | Pairwise differential expression, GO composition, chemosensory classification |
| Sex & Mating Analysis | Volcano plots, Venn diagrams, sex-bias reorganisation, GO term browser |
| Chemosensory Heatmap | Per-family expression across the nine appendage × state groups |
| Data Explorer | Every table behind the paper: lollipops for all 12 contrasts, the mated-female-specific scatter, GO over-representation, gene lookup with BLAST evidence, multivariate statistics, amino-acid identity, comparative counts — all downloadable as CSV |
| Figures & Data | All 5 main and 26 supplementary figures, and all 38 Additional files |

The **Data Explorer** is where the data itself lives. Every view filters, searches
and downloads as CSV, and the significance thresholds are yours to set, so a
reader can check any number in the paper against the table it came from.

## Run it locally

You need Python 3.11. Roughly 2 GB of RAM is enough.

```bash
git clone https://github.com/dorperets-lab/Hillucens_olfactome_app.git
cd Hillucens_olfactome_app

python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install -r requirements.txt
streamlit run app.py
```

It opens at <http://localhost:8501>. Nothing leaves your machine, and there is no
cold start — worth doing if you are working through the data rather than
glancing at it.

To point the sidebar's email helper at your own instance:

```bash
APP_URL="http://localhost:8501" streamlit run app.py
```

## How the numbers are defined

Two conventions matter if you compare the app against the Additional files
directly:

- **Significance** is `adjusted p < 0.001` and `|log₂FC| ≥ 1`. Both are
  adjustable in the sidebar; the defaults are the manuscript's.
- **Sign.** The DESeq2 tables carry `log2FoldChange` with the sign opposite to
  the direction named in the file. The app negates it on load, as the figure
  scripts do, so a positive fold change always means "higher in the first group
  named".
- **Mating-responsive** means significant in the *same direction against both
  virgin groups* — mated female versus virgin female **and** mated female versus
  virgin male. The looser single-contrast definition gives substantially larger
  sets (230 rather than 33 genes in the tarsi).

With the defaults, the app reproduces the manuscript: 535 / 555 / 263 sex-biased
and 146 / 318 / 33 mating-responsive genes in antenna / maxillary palp / tarsi.

## Keeping it in step with the manuscript

Three scripts regenerate the parts of `app.py` that must not drift from the paper:

```bash
python3 build/sync_captions.py      # figure captions, read from the .docx
python3 build/sync_supp_files.py    # the Additional-file download list
python3 build/validate_assets.py    # check every referenced asset exists
```

`sync_captions.py` reads the manuscript directly and keeps only surviving text,
skipping the struck-through and red runs used to mark deletions.

## Citation

> Perets D. *et al.* (2026). Chemosensory transcriptomes of the antenna,
> maxillary palp and tarsi of *Hermetia illucens*.

## Licence

Code released for review alongside the manuscript. The underlying data are the
Additional files of Perets *et al.* 2026; please cite the paper if you use them.
