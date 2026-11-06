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

## Main software

| Step | Software |
|---|---|
| Read processing | [NextITS](https://github.com/vmikk/NextITS) v1.1.0, [LIMA](https://lima.how/) v2.12.0, [cutadapt](https://github.com/marcelm/cutadapt) v5.1, [ITSx](https://microbiology.se/software/itsx/) v1.1.3 ([HMM profiles](https://github.com/Mycology-Microbiology-Center/ITSx_HMMs)), [VSEARCH](https://github.com/torognes/vsearch), [UNCROSS2](https://doi.org/10.1101/400762) |
| Taxonomic identification | BLAST+ v2.16.0 |
| Alignment and trimming | [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.526, [AliView](https://github.com/AliView/AliView) v1.26, [ClipKIT](https://github.com/JLSteenwyk/ClipKIT) v1.4.0 |
| Phylogenetics | [IQ-TREE](https://iqtree.github.io/) v2.2.5 |
| Species delimitation | [VSEARCH](https://github.com/torognes/vsearch), [ASAPy](https://github.com/iTaxoTools/ASAPy) v0.1.2 |
| Metagenome scaffolds | [Infernal](http://eddylab.org/infernal/) v1.1.5 (Rfam [RF01960](https://rfam.org/family/RF01960), [RF02543](https://rfam.org/family/RF02543)) |



