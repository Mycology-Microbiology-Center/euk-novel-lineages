# Discovering 25 novel protist lineages that fill gaps in the eukaryotic tree of life

[![DOI](https://zenodo.org/badge/1448284500.svg)](https://doi.org/10.5281/zenodo.24340941)

<img src="assets/MycologyMicrobiologyCenter_logo.png" alt="Mycology and Microbiology Center" width="140" align="right">

Data and code accompanying the manuscript:

> Tedersoo L. et al. 2026 "Discovering 25 novel protist lineages that fill gaps in the eukaryotic tree of life" (under review)


## About the study

We generated long-read (PacBio HiFi) metabarcoding data covering much of the eukaryotic rRNA operon (SSU-ITS-LSU; primers Euk575F and 21Rngs) from eDNA of >1,100 soil, sediment, water and other environmental samples. Combined with unidentified long reads from the [EUKARYOME](https://eukaryome.org) database, these data revealed 25 deeply branching protist lineages of phylum-level divergence.


## Repository layout

```
alignments/  # per-lineage SSU-ITS-LSU alignments (ClipKIT-trimmed) and ITS alignments
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


## Code

| Script | Step |
|---|---|
| `code/blastn.sh` | BLASTn search against EUKARYOME |
| `code/phylogeny/build_tree.sh` | MAFFT alignment, ClipKIT trimming, IQ-TREE inference |
| `code/species_delimitation/01_pairwise_distances.sh` | Dereplication and all-vs-all ITS identities (VSEARCH) |
| `code/species_delimitation/02_distances_to_mega.R` | Distance matrices (MEGA CSV for ASAPy; R `dist` objects) |
| `code/species_delimitation/03_run_asap.sh` | ASAP species delimitation (ASAPy) |
| `code/species_delimitation/04_parse_asap.R` | Summary of ASAP partitions |
| `code/species_delimitation/05_cluster_counts.R` | Number of putative species by hierarchical clustering at lineage-specific thresholds (`thresholds.tsv`) |
| `code/jgi/01_scaffold_cart.R` | Scaffold lists for the JGI IMG/MER Scaffold Cart |
| `code/jgi/02_extract_rrna.sh` | Combining exported scaffolds; SSU and LSU search with Infernal (Rfam RF01960, RF02543) |
| `code/jgi/03_rrna_operons.R` | Extraction of SSU, ITS and LSU from complete operons |
| `code/jgi/04_blast_ends.sh` | BLASTn of SSU, LSU, ITS and the first/last 500 bp of SSU and LSU |

Each script has a header with usage, inputs and outputs.

### BLASTn

BLASTn settings used for all searches:
`-task blastn -strand both -max_target_seqs 10 -max_hsps 1`.
All other parameters at BLAST+ defaults. To run the search:

```bash
./code/blastn.sh query.fasta EUKARYOME.fasta results_prefix 8
```

### Alignment trimming and tree building

Sites with more than 85% gaps were removed with [ClipKIT](https://github.com/JLSteenwyk/ClipKIT) v1.4.0. Trees were built using [IQ-TREE](https://iqtree.github.io/) v2.2.5:

```bash
# Trim all sites that are above a 0.85 threshold of gappyness
clipkit -m gappy -g 0.85 Euk6i20c.fasta

# Build tree
iqtree -s Euk6i20c.fasta.clipkit -m TEST -bb 1000 -alrt 1000 -T 3
```

`code/phylogeny/build_tree.sh` runs the same steps, starting from unaligned sequences. Manual alignment curation in AliView is not scripted.

### Species delimitation and richness

```bash
cd code/species_delimitation
./01_pairwise_distances.sh ../../alignments_its/ dist_out 8
Rscript 02_distances_to_mega.R dist_out/dist mega
find mega -name "*.csv" | parallel -j1 \
  "ASAP_GENERAL_SEQUENCE_LEN=600 ASAP_GENERAL_MEGA=true ./03_run_asap.sh {} asap/{/.}"
Rscript 04_parse_asap.R asap ASAP
Rscript 05_cluster_counts.R mega/distances.rds thresholds.tsv cluster_counts.tsv
```

### JGI IMG/MER metagenome scaffolds

```bash
cd code/jgi
Rscript 01_scaffold_cart.R selected_scaffolds.tsv scaffold_ids
# upload scaffold_ids/*.txt to the IMG/MER Scaffold Cart, export scaffold data, download *.tar.gz
./02_extract_rrna.sh scaffold_downloads rrna 8
Rscript 03_rrna_operons.R rrna rrna
./04_blast_ends.sh rrna EUKARYOME_SSU.fasta EUKARYOME_LSU.fasta EUKARYOME_ITS.fasta blast 8
```


## Data availability

- Raw PacBio HiFi reads: ENA BioProject PRJEB129568
- Legitype and representative sequences: ENA accessions OZ474986-OZ475166; UNITE; EUKARYOME (v2.0 onwards)
- Holotype material and DNA samples: University of Tartu natural history collections (TUE)
- Archived release: Zenodo [DOI:10.5281/zenodo.24340941](https://doi.org/10.5281/zenodo.24340941)

## Citation

Tedersoo L. et al. 2026 "Discovering 25 novel protist lineages that fill gaps in the eukaryotic tree of life" (under review)


