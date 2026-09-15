# Reproducing Fig1 / Fig3 from scratch (on a new computer)

The bundle `fig1_fig3_repro_pack/` contains **all code, layout files and static assets**
needed to rebuild Fig1 and Fig3 from scratch on another machine (when packed with
`--with-data` it also contains the input data).

---

## 1. Requirements

- **Python 3.10** (no R needed — the Fig1 and Fig3 pipelines never use R)
- Python packages (pip):
  ```
  pip install scanpy anndata h5py infercnvpy matplotlib numpy pandas scipy
  pip install seaborn statsmodels adjustText pypdf pypdfium2 reportlab Pillow openpyxl
  ```
- Data disk: by default at **`D:/bulk-download/`** (the `translate` root in
  `config/paths.py` resolves against it; if you keep the data elsewhere, change the root
  path in `config/paths.py` — every script loads its data through it)

## 2. Directory placement

```
D:/bulk-download/                 <- data (when packed with --with-data it ships as
                                     bulk_download/ inside the bundle; copy it to D:)
<anywhere>/
  └─ repo/                        <- everything in the bundle except bulk_download
                                     (keep the relative structure)
      ├─ tools/  config/  results/assembled/fig1.json fig3.json ...
      ├─ docs/复现指南_fig1_fig3.md (this file)
```

## 3. One-command run

From the **repository root**:

```bash
python tools/rebuild_fig1.py          # Fig1 full (re-render 1a/1b/1c/1f + crop + compose)
python tools/rebuild_fig1.py --fast   # fast: reuse existing 1a/1c outputs, run 1b/1f + compose
python tools/rebuild_fig3.py          # Fig3 full (re-render 3a+3e / 3b+3f / 3g + crop + compose)
```

Outputs:
- `results/assembled/fig1/Fig1_composed.pdf` (174 × 214 mm)
- `results/assembled/fig3/Fig3_composed.pdf` (174 × 337 mm)

## 4. Panel → generator (used when changing content)

| Panel | Generator | Input data |
|---|---|---|
| Fig1a | `tools/10_00_fig1a_only.py` (truncated script) | `GSE182434/adata_processed.h5ad`, `annotation/`, `cnv/`, `Macro_mono_cellmarker.xlsx` |
| Fig1b | `tools/fig1b_dotplot.py` (hand-written dot plot) | `adata_processed.h5ad`, `comparison/cluster_markers.csv` |
| Fig1c | `tools/10_00_fig1a_only.py` (same run) | same as Fig1a (inferCNVpy) |
| Fig1d/1e | **static assets** (svglib-converted; not re-rendered) | — |
| Fig1f | `tools/fig1f_volcano.py` (volcano extracted from `02_02` and re-rendered) | `GSE182434/cellcomm/volcano_data_….csv` |
| Fig3a+3e | `tools/04_01_sankey_only.py` (forced copy) | `GSE182434/scenic/*.csv`, `trajectory/*.csv` |
| Fig3b+3f | `tools/04_01_lineplot_only.py` (forced copy) | same as above |
| Fig3c/3d | **static assets** | — |
| Fig3g | `tools/suoxiao_dot_only.py` | `GSE182434/scenic/MAFB_gseapy_*.csv` |

## 5. Changing the layout

Panel coordinates, row spacing and panel-letter positions are all fixed in
`results/assembled/fig1.json` / `fig3.json` (explicit-coordinate `panels` mode;
`x_mm`/`y_mm` are measured from the page top / page left).
**Changing the layout = editing the JSON; changing the content = editing the panel script.**
Composer: `python tools/fig_compose_pdf.py results/assembled/fig1.json`.
Note that fig3's panel `g` already has a −3 mm left shift baked into its `x_mm`, and the
f–g row spacing is 2.5 mm.

## 6. Known caveats

- If the `PY310` variable in the rebuild scripts does not exist on your machine, it falls
  back automatically to the interpreter running the script; make sure that interpreter has
  all the packages listed above.
- The `pad` / `th` arguments of `vector_crop` are the finalised values (1a uses
  `pad=10`/`th=252` to prevent top clipping; light-grey text needs `th ≥ 250`).
- In design mode the composer `fig_compose_pdf.py` skips `tight_layout` for scripts that
  already call `subplots_adjust` explicitly (the protective patch is built into
  `config/retrofit.py`).
- If the target PDF is locked by a viewer, the composer automatically writes `*_new.pdf`
  instead.
