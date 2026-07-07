# scRNA-seq Analysis: nAChR Expression in Dopaminergic Neurons (Schizophrenia)

Single-cell RNA-sequencing analysis quantifying nicotinic acetylcholine receptor (nAChR) subunit expression in ventral tegmental area (VTA) dopaminergic neurons, performed as part of a minor project on cholinergic dysfunction in schizophrenia (IIT Delhi, Kusuma School of Biological Sciences, Nov 2025). This dataset feeds directly into a companion [NEURON biophysical model](https://github.com/Dhawaljha/pfc-microcircuit-neuron), where the receptor densities derived here were used to constrain simulated α7 and α4β2 nAChR conductances.

## Dataset

Source: [GSE235149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235149) (GEO), profiling ventral tegmental area (VTA) neuron transcriptomes.

## Pipeline

Analysis was performed in RStudio using the **Seurat** package:

1. Log-normalization, variable feature selection, PCA, and nearest-neighbour graph construction on the full dataset.
2. Dopaminergic neuron identification via UMAP projection, filtering for cells co-expressing canonical markers **Th**, **Ddc**, and **Slc6a3** (~58,000 DA cells).
3. Extraction of ligand-gated ion channel (LGIC) subunit expression — percentage of expressing cells (PctExpr) and average log-normalized expression (AvgLogNorm) — ranked to the top ~100 ion-channel-related genes across nicotinic, glutamatergic, GABAergic, and purinergic families.
4. Statistical filtering: a background comparison score (tolerance 0.05) to exclude low-level/stochastic transcription, followed by a binomial detection test (baseline detection probability p₀ = 0.014%) to confirm consistent expression across the DA population.
5. nAChR-specific evaluation of *Chrna*/*Chrnb* gene families, identifying **Chrna7**, **Chrna4**, and **Chrnb2** as significantly expressed — corresponding to the α7 and α4β2 receptor subtypes.

## Files

| File | Description |
|---|---|
| `AChR_DA_expression.csv` | Per-gene AvgLogNormExpr / PctExpressing / NumCellsExpressing for nAChR subunits in DA neurons |
| `Chrna_background_comparison.csv` | Background transcription comparison scores used for statistical filtering |
| `Chrna_binomial_test.csv` | Binomial detection test results confirming significant nAChR subunit expression |
| `DA_summary_LGIC_EI_balance.csv` | Summary excitatory/inhibitory ligand-gated channel balance in DA neurons |
| `LGIC_channel_estimated_channels_per_positive_cell.csv` | Estimated channel counts per expressing cell |
| `LGIC_channel_expression_metrics.csv` / `_generated.csv` | Channel-level expression metrics (assembled subunit combinations) |
| `LGIC_channel_level_summary.csv` | Channel-level summary: PctExpr, NumCells, MeanUMI for nAChR assemblies (α7, α4β2, α6β2β3) and AMPA/NMDA |
| `LGIC_gene_level_UMI_summary.csv` | Gene-level UMI count summary across the top ligand-gated ion channel genes |
| `LGIC_top100_DA.csv` | Top ~100 ranked ion-channel-related genes in the DA population |
| `LGIC_top20_DA_by_avgexpr.csv` | Top 20 LGIC subunits in DA neurons ranked by average log-normalized expression |
| `LGIC_top20_DA_by_pctexpr.csv` | Top 20 LGIC subunits in DA neurons ranked by percent of cells expressing |

These tables correspond to the expression data summarized in Figure 3 of the accompanying project report.

## Code

`VTA_DA_nAChR_analysis.R` is the full Seurat pipeline used to generate these tables, recovered from the original R session history and verified to reproduce the same expression values found in the CSVs above. It covers: downloading and unpacking GSE235149 from GEO, building per-sample Seurat objects, merging and normalizing, filtering to dopaminergic neurons (Th/Slc6a3/Ddc/Slc18a2+), scoring the LGIC gene panel (including the Chrna/Chrnb nAChR subunits), and ranking subunits by expression.

## Key finding

**Chrna7**, **Chrna4**, and **Chrnb2** are significantly expressed in VTA dopaminergic neurons, consistent with published transcriptomic and postmortem evidence for α7/α4β2 nAChR involvement in schizophrenia. These empirically-derived expression levels were used to set biologically grounded receptor densities in the companion NEURON simulation.
