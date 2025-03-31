# EAE Multiomics 2025

<h3>Distinct transcriptomic and epigenomic responses of mature oligodendrocytes during disease progression in a mouse model of multiple sclerosis</h3>

_Chao Zheng<sup>1*</sup>, Bastien Hervé<sup>1*</sup>, Mandy Meijer<sup>1,3</sup>, Leslie Ann Rubio Rodríguez-Kirby<sup>1</sup>, André Ortlieb Guerreiro Cacais<sup>2</sup>, Petra Kukanja<sup>1</sup>, Mukund Kabbe<sup>1</sup>, Tony Jimenez-Beristain<sup>1</sup>, Tomas Olsson<sup>2</sup>, Eneritz Agirre<sup>1</sup>, Gonçalo Castelo-Branco<sup>1</sup>_

_1 Laboratory of Molecular Neurobiology, Department of Medical Biochemistry and Biophysics, Karolinska Institutet, 17177 Stockholm, Sweden_ <br />
_2 Neuroimmunology Unit, Department of Clinical Neuroscience, Center for Molecular Medicine, Karolinska Institute, 171 76 Stockholm, Sweden_ <br />
_3 Max Delbrück Center for Molecular Medicine in the Helmholtz Association (MDC Berlin), Berlin, Germany_ <br />
_* Equal contribution_ <br />
_Correspondence: goncalo.castelo-branco@ki.se_

## Abstract

Multiple sclerosis (MS) is a chronic demyelinating autoimmune disease that targets mature oligodendrocytes (MOLs) and their myelin. MOLs are transcriptionally heterogeneous and can transition to immune-like states in the context of MS. However, the intricacies of their dynamics throughout disease progression remain poorly understood. Here, we employed simultaneous single-cell multiome ATAC and RNA sequencing targeting oligodendroglia (OLGs) from the experimental autoimmune encephalomyelitis (EAE) MS mouse model at different stages of the disease course. We found that the transition to immune OLG states appear already at the early stages of EAE and persist to the late stages of the disease, consistent with epigenetic memory of previous neuroinflammation. Interestingly, transcription factor activity suggested immunosuppression in MOLs at early stages of EAE and we also observed a transitory regenerative cholesterol-associated program in MOLs at this stage. Importantly, different MOLs exhibit a differential responsiveness to EAE, with MOL2 exhibiting a stronger transcriptional immune response than MOL5/6. Moreover, we observed divergent responses at the epigenetic level of MOL2 and MOL5/6 during disease evolution. Thus, our single-cell multiomic resource highlights dynamic and distinct responses of OLG subpopulations to the evolving environment in EAE, which might modulate their response to regenerative therapeutic interventions in MS. 

## Raw data processing

A total of 7 batches were collected from the sequencing facility. Fastq files from 28 samples were processed throughout the 10X genomics standard pipeline. Gene expression and chromatin accessibility libraries were inputted into cellranger-arc2 ‘count’ v2.0.2 with default settings to align the biological readouts on the associated mm10 reference genome v2020-A-2.0.0. Sample aggregation of both transcriptomic and genomics metrics was done using the ‘aggr’ of the same Cellranger executable file, without normalization.

## Documents

**EAE_multiomics_2025** : Jupyter notebooks containing R and/or Python code to recreate figures containing all cell types

**EAE_multiomics_2025_OL** : Jupyter notebooks containing R and/or Python code to recreate figures containing only oligodendrocytes lineage

**MOL_GRN_CZ_BHpeak_selection.R** : R code to generate gene regulatory network and transcription factors plots

**atac_ifn_nf.sh** and **rna_ifn_nf.sh** : Bash code to run Nextflow pipeline to align primary OPC cell culture with and without interferon-γ treatment in bulk atac and rna datasets

**bulk_atac_rna.Rmd** : R markdown to find differentially expressed genes and accessible genes from primary OPC cell culture with and without interferon-γ treatment

## Visualization

All plots from the manuscript are displayed in the notebooks.

To have a quick look at each jupyter notebook in .html format, convert each file as follow :

```
jupyter nbconvert --to html EAE_multiomics_2025.ipynb
```
