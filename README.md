
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
install.packages(
  "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.tar.gz",
  repos = NULL,
  type = "source"
)
```

``` r
install_phyloR <- function(
  url = "https://www.ciblab.net/pub/pkg/phyloR_0.1.0.zip",
  upgrade = "never"   # 不强制升级用户现有包
){
  # 1) 配好仓库（含 Bioc）
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  options(repos = BiocManager::repositories())  # CRAN + Bioc 全部可用

  # 2) 安装 remotes
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # 3) 从 URL 安装，并自动装依赖（Imports/Depends/Suggests[可控]）
  remotes::install_url(
    url,
    dependencies = TRUE,     # 自动装 Imports/Depends；Suggests 也会装（可设 c("Depends","Imports")）
    upgrade = upgrade
  )

  message("Done. Try: library(phyloR)")
}

install_phyloR()
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
