#' Construct Phylogenetic Tree from Processed Sequence Files with Coalescent Method
#'
#' It constructs a phylogenetic tree from multiple sequence files with coalescent method,
#' either DNA or protein sequences, by applying various tree-building methods
#' such as Maximum Likelihood (ML), Neighbor Joining (NJ), UPGMA, Maximum Parsimony (MP),
#' and Bayesian Inference (BI). It processes sequence files, constructs individual trees
#' based on the specified method, and combines them into a supertree.
#' The tree can be displayed if specified.
#'
#' @param seq.files A character vector of sequence file names (without extension)
#'        to be used for constructing the phylogenetic tree.
#' @param seq.type The type of sequence to be processed, either "DNA" or "protein".
#'        Default is "DNA".
#' @param data_dir The directory containing the sequence files. This is a required parameter.
#' @param show_tree Logical, if TRUE the constructed supertree is displayed. Default is TRUE.
#' @param tree_method The method for constructing the phylogenetic tree. Options include:
#'        "ML" (Maximum Likelihood), "NJ" (Neighbor Joining), "UPGMA", "MP" (Maximum Parsimony),
#'        and "BI" (Bayesian Inference). Default is "ML".
#'
#' @return A phylogenetic tree object (of class `phylo`), representing the supertree constructed from
#'         the individual trees built from the sequence files.
#'
#' @details It performs the following steps:
#'   1. Retrieves sequence data from the provided files and processes them based on the specified sequence type.
#'   2. Builds an individual phylogenetic tree for each sequence file using the chosen tree-building method.
#'   3. Combines the individual trees into a supertree.
#'   4. Optionally, displays the supertree if the `show_tree` parameter is TRUE.
#'   5. Returns the supertree object.
#'
#' @return If `show_tree` is TRUE, the function plots the coalescence tree.
#'         If FALSE, it returns the constructed phylogenetic tree as a "multiPhylo" object.
#'
#' @importFrom adegenet fasta2DNAbin
#' @importFrom ape dist.dna dist.aa nj di2multi prop.clades
#' @importFrom phangorn as.phyDat modelTest pml_bb upgma pratchet acctran superTree
#' @importFrom babette bbt_run_from_model
#' @importFrom beautier create_inference_model
#' @importFrom beastierinstall install_beast2
#' @importFrom beastier create_beast2_options
#'
#' @examples
#' dna <- c("K01939", "K03644", "K00797", "K00927", "K00088", "K02257", "K00164",
#'          "K00820", "K06158", "K00008")
#' data_dir <- system.file("extdata", "sequences", package = "phyloR")
#' tree1 <- coalescent_tree(seq.file = dna,
#'                          seq.type = "protein",
#'                          data_dir = data_dir,
#'                          tree_method = "NJ",
#'                          show_tree = TRUE)
#' tree2 <- coalescent_tree(seq.file = dna,
#'                          seq.type = "DNA",
#'                          data_dir = data_dir,
#'                          tree_method = "UPGMA")
#'
#' @export
coalescent_tree <- function(seq.files,
                            data_dir = NULL,
                            seq.type = "DNA",
                            tree_method = "ML",
                            show_tree = TRUE) {

  # Check if the data directory is provided
  if(is.null(data_dir)){
    stop("Please input a valid 'data_dir'.")
  }

  mtree <- list()

  # Loop through sequence files
  for (i in seq_along(seq.files)) {
    seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")

    # DNA sequence processing
    if (toupper(seq.type) == "DNA") {
      dna <- adegenet::fasta2DNAbin(seq.file)

      # Tree construction based on the selected method for DNA sequences
      if (toupper(tree_method) == "ML") {
        phydata <- phangorn::as.phyDat(dna)
        modeltest <- phangorn::modelTest(phydata)
        pml <- phangorn::pml_bb(modeltest)
        pml <- pml$tree
      } else if (toupper(tree_method) == "NJ") {
        dist_matrix <- ape::dist.dna(dna)
        pml <- ape::nj(dist_matrix)
      } else if (toupper(tree_method) == "UPGMA") {
        dist_matrix <- ape::dist.dna(dna)
        pml <- phangorn::upgma(dist_matrix)
      } else if (toupper(tree_method) == "MP") {
        # Maximum Parsimony (MP) method using ratchet algorithm
        phydata <- phangorn::as.phyDat(dna)
        treeRatchet <- phangorn::pratchet(phydata, trace = 0)
        treeRatchet <- phangorn::acctran(treeRatchet, phydata)
        treeRatchet <- ape::di2multi(treeRatchet)

        if (inherits(treeRatchet, "multiPhylo")) {
          treeRatchet <- unique(treeRatchet)
        }
        pml <- treeRatchet
      } else if (toupper(tree_method) == "BI") {
        # Bayesian Inference (BI) method using BEAST2
        remotes::install_github("richelbilderbeek/beastierinstall")
        requireNamespace("beastierinstall")

        install_beast2_with_check <- function() {
          is_installed <- Sys.which("beast")
          if (nzchar(is_installed)) {
            message("BEAST2 is already installed!")
          } else {
            tryCatch({
              beastierinstall::install_beast2()
              message("BEAST2 installation was successful!")
            }, error = function(e) {
              message("Installation failed: ", e$message)
            })
          }
        }

        install_beast2_with_check()

        # Run BEAST2 and obtain the tree
        outputs <- babette::bbt_run_from_model(seq.file,
                                               inference_model = beautier::create_inference_model(),
                                               beast2_options = beastier::create_beast2_options())
        valid_extensions <- c("fas", "fasta")
        file_extension <- tools::file_ext(seq.file)
        if (!(file_extension %in% valid_extensions)) {
          stop(paste("Invalid file extension for", seq.file, ". Only 'fas' and 'fasta' files are allowed."))
        }
        if(file_extension == "fas"){
          file_name <- sub("\\.fas$", "", basename(seq.file))
        }
        if(file_extension == "fasta"){
          file_name <- sub("\\.fasta$", "", basename(seq.file))
        }
        tmp <- paste(file_name, "trees", sep = "_")
        treeBI <- outputs[[tmp]][[which.max(outputs$estimates$posterior)]]
        con.trees <- outputs[[tmp]][3:length(outputs[[tmp]])]
        posterior <- round(ape::prop.clades(treeBI, con.trees)/9999*100)
        pml <- treeBI
      } else {
        stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
      }

      mtree[[i]] <- pml

    } else if (tolower(seq.type) == "protein") {

      # Protein sequence processing
      protein <- adegenet::fasta2DNAbin(seq.file)

      # Tree construction based on the selected method for protein sequences
      if (toupper(tree_method) == "ML") {
        phydata <- phangorn::as.phyDat(protein)
        modeltest <- phangorn::modelTest(phydata)
        pml <- phangorn::pml_bb(modeltest)
        pml <- pml$tree
      } else if (toupper(tree_method) == "NJ") {
        dist_matrix <- ape::dist.aa(protein)
        pml <- ape::nj(dist_matrix)
      } else if (toupper(tree_method) == "UPGMA") {
        dist_matrix <- ape::dist.aa(protein)
        pml <- phangorn::upgma(dist_matrix)
      } else if (toupper(tree_method) == "MP") {
        # Maximum Parsimony (MP) method using ratchet algorithm
        phydata <- phangorn::as.phyDat(protein)
        treeRatchet <- phangorn::pratchet(phydata, trace = 0)
        treeRatchet <- phangorn::acctran(treeRatchet, phydata)
        treeRatchet <- ape::di2multi(treeRatchet)

        if (inherits(treeRatchet, "multiPhylo")) {
          treeRatchet <- unique(treeRatchet)
        }
        pml <- treeRatchet
      } else if (toupper(tree_method) == "BI") {
        # Bayesian Inference (BI) method using BEAST2
        remotes::install_github("richelbilderbeek/beastierinstall")
        requireNamespace("beastierinstall")

        install_beast2_with_check <- function() {
          is_installed <- Sys.which("beast")
          if (nzchar(is_installed)) {
            message("BEAST2 is already installed!")
          } else {
            tryCatch({
              beastierinstall::install_beast2()
              message("BEAST2 installation was successful!")
            }, error = function(e) {
              message("Installation failed: ", e$message)
            })
          }
        }

        install_beast2_with_check()

        # Run BEAST2 and obtain the tree
        outputs <- babette::bbt_run_from_model(seq.file,
                                               inference_model = beautier::create_inference_model(),
                                               beast2_options = beastier::create_beast2_options())
        valid_extensions <- c("fas", "fasta")
        file_extension <- tools::file_ext(seq.file)
        if (!(file_extension %in% valid_extensions)) {
          stop(paste("Invalid file extension for", seq.file, ". Only 'fas' and 'fasta' files are allowed."))
        }
        if(file_extension == "fas"){
          file_name <- sub("\\.fas$", "", basename(seq.file))
        }
        if(file_extension == "fasta"){
          file_name <- sub("\\.fasta$", "", basename(seq.file))
        }
        tmp <- paste(file_name, "trees", sep = "_")
        treeBI <- outputs[[tmp]][[which.max(outputs$estimates$posterior)]]
        con.trees <- outputs[[tmp]][3:length(outputs[[tmp]])]
        posterior <- round(ape::prop.clades(treeBI, con.trees)/9999*100)
        pml <- treeBI
      } else {
        stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
      }

      mtree[[i]] <- pml
    } else {
      stop("Unknown sequence type, please choose 'DNA' or 'protein'.")
    }
  }

  # Combine individual trees into a supertree
  class(mtree) <- "multiPhylo"
  supertree <- phangorn::superTree(mtree)

  # Plot or return the supertree based on the `show_tree` parameter
  if(show_tree == TRUE){
    plot(supertree)
    return(supertree)
  } else {
    return(supertree)
  }
}
