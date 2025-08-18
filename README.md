
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

[![R-CMD-check](https://github.com/libcell/phyloR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/libcell/phyloR/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Phylogenetic Analysis Toolkit for R

**phyloR** provides an integrated workflow for molecular phylogenetics, from sequence retrieval to tree comparison. Key features include:

- [DNA] Ortholog identification from major databases
- [Search] Multiple sequence alignment and trimming
- [Tree] Phylogenetic tree construction (ML/NJ/MP)
- <-> Tree comparison and visualization
- [Plot] Integrated tidyverse compatibility

## Installation

You can install the development version of phyloR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("libcell/phyloR")
```

Alternatively, you can use **install_ciblab()** to install it under the MacOS or Linux environment: 

``` r
source("https://ciblab.net/pub/install_ciblab.R")
install_ciblab("phyloR")
```

## Usage

The documentation of phyloR is availabel at our
[Simsite](http://www.ciblab.net/software/phyloR/). In addition, you can also learn the usage of every simulation methods we have collected.

Users can also download our phyloR [Docker
Image](https://hub.docker.com/repository/docker/duohongrui/simpipe/general)
and use
[**simpipe2docker**](https://github.com/duohongrui/simpipe2docker)
package for linking R environmrnt and the Docker container. For more details, please refer to
[Simsite](http://www.ciblab.net/software/Simsite/).

## Contact

If you have any question, please email to Feifei Li (<libcell@cqnu.edu.cn>) or raise an issue for that.

## Citation

Zou Y, Zhang Z, Zeng Y, *et al*. Common methods for phylogenetic tree construction and their implementation in R. ***Bioengineering***, 2024, 11(5): 480.
