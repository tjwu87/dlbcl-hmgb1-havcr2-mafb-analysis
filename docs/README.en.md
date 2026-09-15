# docs/ — index of audit and mapping documents

| File | Language | Content |
|---|---|---|
| `panel_mapping_final.en.md` | **English** | Authoritative main-figure ↔ panel ↔ code mapping (Fig1–Fig7), with confidence ratings and known discrepancies |
| `figure_panel_map.en.md` | **English** | Panel composition and provenance of each main figure, with size diagnostics and layout advice, plus a re-export checklist |
| `复现指南_fig1_fig3.en.md` | **English** | Step-by-step guide to rebuilding Fig1 and Fig3 from scratch on a new computer |
| `reproduction_status_report_20260913.en.md` | **English** | Reproduction status report: per-figure entry point, fidelity measurements, environment, fixes and remaining to-dos |
| `figure_code_map_final.md` | Chinese | Per-panel mapping for main figures **and** supplementary figures, with per-panel verification status |
| `code_version_audit.md` | Chinese | Audit of code versions: which script version produces which figure version |
| `复现状态报告_20260913.md` | Chinese | Chinese original of the reproduction status report (the authoritative record for line numbers and measurements) |
| `panel_mapping_final.md` | Chinese | Chinese original of the panel ↔ code mapping |
| `figure_panel_map.md` | Chinese | Chinese original of the panel composition / provenance document |
| `复现指南_fig1_fig3.md` | Chinese | Chinese original of the Fig1/Fig3 reproduction guide |

## A note on language

These documents are **internal audit records** produced while tidying this code base up
for release. The Chinese originals are kept because they are the primary evidence trail:
they carry the exact line numbers, byte-level comparisons and measurements behind every
claim made in this repository.

Where an English version exists it is a faithful translation of the same content; where
only a Chinese file is listed, that document has not yet been translated.

## What these documents are for

- **Traceability** — if you want to know why a particular threshold, cohort list or panel
  layout looks the way it does, the answer is recorded here.
- **Boundary-setting** — the known discrepancies sections state plainly where a
  reproduced figure differs from the published version and why (for example, a panel
  source that comes from an earlier version of the authors' asset library).

They are **not** required to run the pipeline: `run_all.py` and the per-section READMEs
are self-contained.
