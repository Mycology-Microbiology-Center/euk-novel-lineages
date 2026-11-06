# 25 novel lineages of Eukaryotes

<img src="assets/MycologyMicrobiologyCenter_logo.png" alt="Mycology and Microbiology Center" width="140" align="right">

Data and code accompanying the manuscript:

> Tedersoo L. et al. 2026 "Discovering 25 novel phyla that fill gaps in the eukaryotic tree of life"


## About the study

We generated long-read (PacBio HiFi) metabarcoding data covering much of the eukaryotic rRNA operon (SSU-ITS-LSU; primers Euk575F and 21Rngs) from eDNA of >1,100 soil, sediment, water and other environmental samples. Combined with unidentified long reads from the [EUKARYOME](https://eukaryome.org) database, these data revealed 25 deeply branching protist lineages of phylum-level divergence.


## Repository layout

```
alignments/  # per-lineage SSU-ITS-LSU alignments (ClipKIT-trimmed) and ITS alignments
trees/       # IQ-TREE tree files (Newick) for the global and per-lineage analyses
code/        # code used for the analysis      
```

