
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
# ----------------------------------------------------------
# install_phyloR(): install packages from local dir or ciblab
# Features:
#   - Auto-detect highest local version in working directory
#   - OS-aware remote format fallback (.zip on Windows, .tar.gz otherwise)
#   - Resolves dependencies via remotes with CRAN + Bioconductor repos
# ----------------------------------------------------------

install_ciblab <- function(
  pkg        = "phyloR",                          # package name
  version    = "0.1.0",                               # version; if NULL, auto-detect highest local version
  base_url   = "https://www.ciblab.net/pub/pkg",   # remote base URL
  upgrade    = "never",                            # remotes::install_* upgrade policy
  deps       = c("Depends", "Imports"),            # which dependency classes to install
  verbose    = TRUE                                # print progress
) {
  # ---- 0) Small helpers ----
  is_windows <- identical(.Platform$OS.type, "windows")

  # Pattern to match local files: pkg_1.2.3(.9000 etc).(tar.gz|tgz|zip|tar.bz2)
  build_local_pattern <- function(pkg) {
    # version segment: digits and dots + optional suffixes; keep it simple, rely on package_version() later
    sprintf("^%s_([0-9][0-9\\.\\-]*[0-9A-Za-z\\.]*)\\.(tar\\.gz|tgz|zip|tar\\.bz2)$", pkg)
  }

  # Find highest local version file in current working directory
  find_latest_local <- function(pkg) {
    pat <- build_local_pattern(pkg)
    lf  <- list.files(".", pattern = pat, ignore.case = FALSE, full.names = FALSE)
    if (!length(lf)) return(NULL)

    # Extract version substrings by regex capture group
    m <- regexec(pat, lf)
    cap <- regmatches(lf, m)
    vers <- vapply(cap, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    keep <- !is.na(vers)
    lf   <- lf[keep]
    vers <- vers[keep]

    # Use package_version to compare semantically
    ord <- order(package_version(vers), decreasing = TRUE)
    list(file = lf[ord][1], version = vers[ord][1])
  }

  # ---- 1) Ensure repos (CRAN + Bioc) and remotes are available ----
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  options(repos = BiocManager::repositories())

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # ---- 2) Try local highest version if 'version' not given ----
  local_choice <- NULL
  if (is.null(version)) {
    local_choice <- find_latest_local(pkg)
    if (!is.null(local_choice)) {
      version <- local_choice$version
      if (verbose) {
        message(sprintf("Found local %s, highest version: %s", pkg, version))
        message("Installing from local file: ", local_choice$file)
      }
      # Install local file directly
      ok <- try(
        remotes::install_local(
          path = local_choice$file,
          dependencies = deps,
          upgrade = upgrade
        ),
        silent = TRUE
      )
      if (!inherits(ok, "try-error")) {
        if (verbose) message("Installation finished. Try: library(", pkg, ")")
        return(invisible(TRUE))
      } else {
        if (verbose) message("Local install failed, will try remote fallback…")
      }
    } else {
      if (verbose) message("No local tarball/zip found for ", pkg, " in current directory.")
    }
  }

  # ---- 3) Remote install (OS-aware format) ----
  # If version is still NULL (no local file), you can set a default or stop with message.
  if (is.null(version)) {
    stop("No local package found and 'version' not specified.\n",
         "Please provide 'version' (e.g., '0.1.0') or place a file like ",
         sprintf("'%s_0.1.0.tar.gz' in the working directory.", pkg))
  }

  # Compose remote URLs for both formats and try preferred first
  file_zip <- sprintf("%s/%s_%s.zip",    base_url, pkg, version)     # Windows binary
  file_src <- sprintf("%s/%s_%s.tar.gz", base_url, pkg, version)     # Source

  urls <- if (is_windows) c(file_zip, file_src) else c(file_src, file_zip)

  if (verbose) {
    message(sprintf("Detected OS: %s", if (is_windows) "Windows" else "non-Windows (macOS/Linux)"))
    message(sprintf("Installing %s %s from ciblab…", pkg, version))
  }

  last_err <- NULL
  for (u in urls) {
    if (verbose) message("Trying: ", u)
    ok <- try(
      remotes::install_url(
        url = u,
        dependencies = deps,
        upgrade = upgrade
      ),
      silent = TRUE
    )
    if (!inherits(ok, "try-error")) {
      if (verbose) message("Installation finished. \nTry: library(", pkg, ")")
      return(invisible(TRUE))
    } else {
      last_err <- ok
      if (verbose) message("Failed with this URL, will try the fallback format…")
    }
  }

  stop(sprintf("Failed to install %s %s from ciblab.\nLast error:\n%s",
               pkg, version, as.character(last_err)))
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
