
<img src="man/figures/phyloR_logo.png" align="right" width = "158px" height="183px"/>

# phyloR: An Automated R Package for Sequence Retrievd and Phylogenetic Analysis in Ecological and Environment Genomics.

PhyloR is an R package used to obtain orthologs and sequences, including sequence alignment and trimming, to construct phylogenetic or species trees, and even to compare different trees.

## Installation

You can install the development version of phyloR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("libcell/phyloR")
```

Or, you can also install it using the following R script: 

``` r
# -------------------------------
# Auto installer for phyloR
# - Windows   : install binary (.zip)
# - macOS/Linux: install source  (.tar.gz)
# Dependencies are resolved via remotes (CRAN + Bioconductor repos)
# -------------------------------

install_phyloR <- function(
  win_url = "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.zip",      # Windows binary
  src_url = "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.tar.gz",   # Source for macOS/Linux
  upgrade = "never"   # do not force upgrade of already-installed packages
) {
  # 0) Detect OS
  # .Platform$OS.type returns "windows" or "unix" (macOS/Linux fall into "unix")
  is_windows <- identical(.Platform$OS.type, "windows")

  # 1) Ensure CRAN + Bioconductor repos are available
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  options(repos = BiocManager::repositories())  # make CRAN + Bioc available

  # 2) Ensure remotes is available
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # 3) Choose the proper URL based on OS
  pkg_url <- if (is_windows) win_url else src_url

  message(sprintf("Detected OS: %s", if (is_windows) "Windows" else "non-Windows (macOS/Linux)"))
  message(sprintf("Installing phyloR from: %s", pkg_url))

  # 4) Install from URL with dependency resolution
  #    - dependencies = TRUE installs Depends/Imports (and Suggests unless you restrict)
  #    - upgrade controls whether to upgrade existing packages
  remotes::install_url(
    url = pkg_url,
    dependencies = TRUE,
    upgrade = upgrade
  )

  message("Installation finished. Try: library(phyloR)")
}

# Run the installer
install_phyloR()
```

Alternatively, you can use **install_ciblab()** to install it under the MacOS or Linux environment: 

``` r
source("https://ciblab.net/pub/install_ciblab.R")
install_ciblab()
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
