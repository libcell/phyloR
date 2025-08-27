
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-red.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![codecov](https://codecov.io/gh/libcell/phyloR/branch/main/graph/badge.svg)](https://codecov.io/gh/libcell/phyloR)
[![codecov](https://codecov.io/gh/libcell/phyloR/branch/main/graph/badge.svg)](https://codecov.io/gh/libcell/phyloR)
[![Codecov](https://codecov.io/github/libcell/phyloR/coverage.svg?branch=main)](https://codecov.io/github/libcell/phyloR?branch=main)
[![Codecov](https://codecov.io/github/ropensci/drake/coverage.svg?branch=master)](https://codecov.io/github/ropensci/drake?branch=master)

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

### Alternatively, using **install_ciblab()** to install phyloR from CibLab repository 

``` r
source("https://ciblab.net/pub/install_ciblab.R")
install_ciblab("phyloR")
```

## Quick Start

### Examples
``` r
# Define specific gene and species 
gene_id <- "K00826"
species.list <- c("hsa", "sce", "dme", "cel",
                  "xla", "gga", "ssc", "rno",
                  "mmu", "mcc", "gmx", "bta",
                  "ece", "zma", "osa", "ath")

# Create a new directory for sequences files
temp_dir <- tempdir()
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir, recursive = TRUE)
}

# Retrieve orthologous gene information for the provided species
species_info <- get_orthologs(gene_id = gene_id,
                              id.type = "ko_id",
                              species.list = species.list,
                              species.type = "abbspname")

# Process species names and get the KEGG IDs
species <- tolower(species_info[, 3])
gene_ids <- paste(species, species_info[, 1], sep = ":")

# Retrieve sequences for KEGG IDs base on the sequence type
seqset <- get_kegg_sequences(gene_ids = gene_ids,
                          id.type = "kegg_id",
                          seq.type = "DNA")

# Write the sequences to a FASTA file
names(seqset) <- spnames
output_path <- file.path(temp_dir, "sequences.fasta")
seqinr::write.fasta(seqset, names = names(seq), file.out = output_path)

# Prepare output file path for processed sequences
data_file <- paste(gene_id, "fasta", sep = ".")
output_file <- file.path(temp_dir, data_file)


# Align and trim the sequences using the selected alignment method
processed_seq <- align_trim(seq.file = output_path,
                            seq.type = "DNA",
                            method = "ClustalW",
                            output_file = output_file)

# Construct a phylogenetic tree based on the processed sequences using the selected method
tree <- gene_tree(seq.file = output_file,
                  seq.type = "DNA",
                  tree_method = "NJ",
                  show_tree = TRUE)

# Construct the gene tree using DNA sequences of K00826

tree <- gene_tree_fetch(gene_id = "K00826",
                        species.list = c("hsa", "sce", "dme", "cel",
                                         "xla", "gga", "ssc", "rno",
                                         "mmu", "mcc", "gmx", "bta",
                                         "ece", "zma", "osa", "ath"),
                        species.type = "abbspname",
                        seq.type = "DNA",
                        tree_method = "NJ",
                        show_tree = TRUE)
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