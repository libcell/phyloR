
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

The goal of simpipe is to establish a standard pipeline to estimating
parameters from real datasets, simulating single-cell RNA-seq datasets
and evaluating the simulated data from general and functional
perspectives.

## Installation

You can install the development version of simpipe from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("libcell/phyloR")
```

Or, you also can use **install.packages()** to install our package with following R script under the Windows Operation System: 

``` r
install.packages(
  "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.zip",
  repos = NULL,       # 不去CRAN，而是直接用给定URL
  type = "win.binary" # 指定为Windows二进制包
)
```

Alternatively, you can use **install.packages()** to install it under the MacOS or Linux environment: 

``` r
install.packages(
  "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Usage

The documentation of simpipe is availabel at our
[Simsite](http://www.ciblab.net/software/Simsite/). In addition, you can
also learn the usage of every simulation methods we have collected.

Users can also download our simpipe [Docker
Image](https://hub.docker.com/repository/docker/duohongrui/simpipe/general)
and use
[**simpipe2docker**](https://github.com/duohongrui/simpipe2docker)
package for linking R environmrnt and the Docker container. For more
details, please refer to
[Simsite](http://www.ciblab.net/software/Simsite/).

## Contact

If you have any question, please email to Hongrui Duo
(<libcell@cqnu.edu.cn>) or raise an issue for that.

## Citation

Duo H, Li Y, Lan Y, *et al*. Systematic evaluation with practical
guidelines for single-cell and spatially resolved transcriptomics data
simulation under multiple scenarios. ***Genome Biology***, 2024, 25(1):
145.
