
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

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

### Platform Requirements

- ***for Windows users:*** Please ensure that Rtools is installed and properly configured before installing from source.

- ***for MacOS users:*** Some dependencies (e.g., rgl) require XQuartz to be installed.

- ***for Linux users:*** Make sure to install system libraries required by Bioconductor packages (e.g., libxml2, libcurl, libssl).

### From GitHub (Latest Version)
``` r
# Install via remotes (recommended), or devtools
if (!require("remotes")) install.packages("remotes")
remotes::install_github("libcell/phyloR", dependencies = TRUE)
```

### Alternatively, using **install_ciblab()** to install phyloR from our Labsite 

``` r
source("https://ciblab.net/pub/install_ciblab.R")
install_ciblab("phyloR")
```

## Quick Start

### Examples
``` r
library(phyloR)

# Get orthologs
taxinfo1 <- get_orthologs(gene_id = "K00161",
                          id.type = "KO_id",
                          species.list = "Homo sapiens",
                          species.type = "scientificname")

# Build gene tree
dna <- c("K01939", "K03644", "K00797", "K00927", "K00088", "K02257", "K00164",
         "K00820", "K06158", "K00008")
data_dir <- system.file("extdata", "sequences", package = "phyloR")
tree1 <- coalescent_tree(seq.file = dna,
                         seq.type = "protein",
                         data_dir = data_dir,
                         tree_method = "NJ",
                         show_tree = TRUE)

# Obtain species tree
# Example 1:
species1 <- c("Homo sapiens", "Pan troglodytes", "Mus musculus",
              "Rattus norvegicus","Canis lupus familiaris", "Felis catus")
tree1 <- species_tree(species = species1, species.type = "scientificname")

# Example 2:
species2 <- c("9606", "9598", "10090", "9615", "9685", "10116")
tree2 <- species_tree(species = species2, species.type = "taxonomic_id")

# Example 3:
species3 = c("ath", "gmx", "zma", "osa",
             "dme", "cel", "mmu", "rno",
             "hsa", "mcc", "ssc", "bta",
             "gga", "xla", "sce", "ece")
tree3 <- species_tree(species = species3, species.type = "abbspname")
```

### Detailed Guides
| Topic                      | Command/Resource                  |
|----------------------------|-----------------------------------|
| 🧬 Ortholog Retrieval    | `vignette("orthologs")`           |
| 🔍 Sequence Alignment  | `vignette("alignment-methods")`   |
| 🌳 Tree Construction   | `vignette("tree-building")`       |
| ↔️ Tree Comparison        | `vignette("tree-comparison")`     |

### Online Resources
- [ [WEB] Website Docs ](https://libcell.github.io/phyloR/)
- [ [PDF] User Manual ](docs/phyloR_manual.pdf)
- [ [VID] Tutorial Videos ](https://youtube.com/playlist?list=XXX)

## Contact

If you have any question, please email to Feifei Li (<libcell@cqnu.edu.cn>) or raise an issue for that.

## Citation
If you use phyloR in your research, please cite the following publication:

Zou Y, Zhang Z, Zeng Y, *et al*. Common methods for phylogenetic tree construction and their implementation in R. ***Bioengineering***, 2024, 11(5): 480. 