# CPAP

**CPAP - Clustering and Phylogenetic Analyzer of Proteins/Nucleotides**

ver 0.3.0 (2025-11-18)

## Description

CPAP generates heatmap-based visualization according to ...

- protein similarity/%-identity using blastp/diamond/mmseqs
- nucleotide similarity/%-identity using blastn

CPAP first computes sequence identity based on BLASTp results, then converted into an identity matrix and visualizes it by a heatmap and a dendrogram.

## Installation

```bash
# 1. Clone the repository
git clone <repo_url> && cd CPAP

# 2. Create conda environment (micromamba, mamba, or conda)
micromamba env create -f environment.yaml

# 3. Install Ruby gems
micromamba run -n CPAP_v0.3.0 bundle install
```

Or, run `bash install.sh` to perform steps 2-3 in one go.

> **Note:** You can use `mamba` or `conda` in place of `micromamba`.

## Usage

```
$ CPAP [options] <input fasta> <output dir>
```

- `<input fasta>` should be in FASTA format.
- `<output dir>` should not exist.

## Dependencies

- BLAST+
- diamond
- mmseqs2
- ruby (ver >=3.0)
- R (ver >=3.0)
    - R package `gplots`     -- is used for heatmap drawing. (for install, try `install.packages('gplots')` in R terminal)
    - R package `phylogram`  -- is used for generation of a newick formatted dendrogram.
    - R package `dendextend` -- [optional] This is requried when a phylogenetic tree is NOT given (without `--phy` or `--phy-as-chronogram`).
    - R package `ape`        -- [optional] This is requried when a phylogenetic tree is given (with `--phy` or `--phy-as-chronogram`).
    - R package `phangorn`   -- [optional] This is requried when a phylogenetic tree is given (with `--phy` or `--phy-as-chronogram`).

## Options

### general

```
-h, --help
-v, --version
--overwrite   (default: off) -- overwrite output directory
--no-heatmap  (default: off) -- skip heatmap/dendrogram steps and stop after identity/similarity matrix generation
--mode   [blastp|diamond|mmseqs|blastn] (default: blastp) -- choose blastp/diamond/mmseqs (for protein) or blastn (for nucleotide) for similarity search, currently only blastp/blastn is implemented.
```

### similarity measure

```
--measure     [identity|sim-score] (default: identity) -- choice of similarity measure for heatmap generation (identity: %identity of blastp, sim-score: normalized score like Sg used in ViPTree).
```

### computing

```
--ncpus       [int] (default: 4)        -- the number of cpus to use
```

### blastp

```
--min-aln-len [int] (default: 20)       -- minimum amino acid length of acceptable HSP
--min-idt     [int] (default: 20)       -- %-identity of acceptable HSP
--dbsize      [int] (default: 100000000)
--matrix      [str] (default: BLOSUM62)
--evalue      [num] (default: 0.01)
```

### blastn

```
--min-aln-len [int] (default: 60)       -- minimum nucleotide length of acceptable HSP
--min-idt     [int] (default: 60)       -- %-identity of acceptable HSP
--dbsize      [int] (default: 100000000)
--evalue      [num] (default: 0.01)
```

### provide a phylogenetic tree used as a dendrogram

```
--fphy               [str] (default: disabled) -- to provide a phylogenetic tree file (newick or nexus format). The tree is directly converted to a dendrogram of the heatmap.
--fphy-as-chronogram [str] (default: disabled) -- to provide a phylogenetic tree file (newick or nexus format). The tree is converted to chronogram and used as a dendrogram of the heatmap.
```

### clustering

```
--clust-method [ALL|average|ward.D|ward.D2|single|complete|mcquitty|median|centroid]  (default: ALL)   -- method of hclust function in R. 'ALL' performs all methods.
```

## Output files

### blastp identity mode

```
result/idt.tsv                        -- blastp %identity matrix
result/idt-heat.<clust-method>.pdf    -- heatmap pdf file
result/idt.ordered.<clust-method>.tsv -- blastp %identity matrix, sorted in the same order of the dendrogram.
```
