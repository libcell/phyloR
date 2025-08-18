
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

[![R-CMD-check](https://github.com/libcell/phyloR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/libcell/phyloR/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Phylogenetic Analysis Toolkit for R

**phyloR** provides an integrated workflow for molecular phylogenetics, from sequence retrieval to tree comparison. Key features include:

- 🧬 Ortholog identification from major databases
- 🔍 Multiple sequence alignment and trimming
- 🌳 Phylogenetic tree construction (ML/NJ/MP)
- ↔️ Tree comparison and visualization
- 📊 Integrated tidyverse compatibility

## Installation

### From GitHub (Latest Version)
``` r
# Install via remotes (recommended), or devtools
if (!require("remotes")) install.packages("remotes")
remotes::install_github("libcell/phyloR", dependencies = TRUE)

# For stable releases, specify version tag:
# remotes::install_github("libcell/phyloR@v0.1.0")
```

### Alternatively, using **install_ciblab()** to install phyloR from our Labsite 

``` r
source("https://ciblab.net/pub/install_ciblab.R")
install_ciblab("phyloR")
```

## Quick Start

``` r
library(phyloR)

# Get orthologs
ortho <- get_orthologs(gene = "COX1", 
                       species = c("human", "mouse", "zebrafish"))

# Align sequences
aln <- align_sequences(ortho, method = "mafft")

# Build tree
tree <- build_tree(aln, method = "RAxML")

# Visualize
plot_tree(tree) + 
  ggtree::theme_tree2()
```



## Contact

If you have any question, please email to Feifei Li (<libcell@cqnu.edu.cn>) or raise an issue for that.

## Citation

Zou Y, Zhang Z, Zeng Y, *et al*. Common methods for phylogenetic tree construction and their implementation in R. ***Bioengineering***, 2024, 11(5): 480.
