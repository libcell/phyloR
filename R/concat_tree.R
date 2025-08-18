#' Construct Phylogenetic Tree from Processed Sequence Files with Concatenate Method
#'
#' It reads multiple sequence files, performs alignment, and constructs a
#' phylogenetic tree based on the concatenated sequences. It supports both DNA and
#' protein sequences and offers various tree construction methods including Maximum
#' Likelihood (ML), Neighbor Joining (NJ), UPGMA, Maximum Parsimony (MP), and
#' Bayesian Inference (BI).
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
#' @return A phylogenetic tree object of class `phylo` representing the tree constructed
#'         from the sequences. If `show_tree` is TRUE, the tree is plotted.
#'
#' @details This function performs the following steps:
#'   1. Reads sequence files from the specified directory and converts them to DNAbin
#'      or AAbin format depending on the sequence type.
#'   2. Aligns and concatenates the sequences.
#'   3. Constructs a phylogenetic tree based on the selected tree method.
#'   4. Optionally displays the tree if `show_tree` is set to TRUE.
#'   5. Returns the phylogenetic tree object.
#'
#' @importFrom adegenet fasta2DNAbin
#' @importFrom ape dist.dna dist.aa nj di2multi prop.clades as.DNAbin as.AAbin
#' @importFrom phangorn as.phyDat modelTest pml_bb upgma pratchet acctran
#' @importFrom methods new
#' @importFrom apex concatenate
#' @importFrom babette bbt_run_from_model
#' @importFrom beautier create_inference_model
#' @importFrom beastier create_beast2_options
#' @importFrom beastierinstall install_beast2
#'
#' @examples
#' dna <- c("K01939", "K03644", "K00797", "K00927", "K00088", "K02257", "K00164",
#'          "K00820", "K06158", "K00008")
#' data_dir <- system.file("extdata", "sequences", package = "phyloR")
#' tree1 <- concat_tree(seq.file = dna,
#'                      seq.type = "DNA",
#'                      data_dir = data_dir,
#'                      tree_method = "NJ")
#' tree2 <- concat_tree(seq.file = dna,
#'                      seq.type = "DNA",
#'                      data_dir = data_dir,
#'                      tree_method = "UPGMA")
#'
#' @export
concat_tree <- function(seq.files,
                        data_dir = NULL,
                        seq.type = "DNA",
                        tree_method = "ML",
                        show_tree = TRUE) {

  # Check if the data directory is provided
  if(is.null(data_dir)){
    stop("Please input a valid 'data_dir'.")
  }

  # Handle DNA sequences case
  if (toupper(seq.type) == "DNA") {
    DNAbin_list <- list()

    # Loop through each sequences file and convert it to DNAbin format
    for (i in seq_along(seq.files)) {
      seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")
      dna_data <- adegenet::fasta2DNAbin(seq.file)
      DNAbin_list[[i]] <- dna_data
    }

    # DNA sequences processing
    multidna <- methods::new("multidna", DNAbin_list)
    phydata <- phangorn::as.phyDat(apex::concatenate(multidna))

    # Tree construction based on the selected method for DNA sequences
    if (toupper(tree_method) == "ML") {
      modeltest <- phangorn::modelTest(phydata)
      pml <- phangorn::pml_bb(modeltest)
    } else if (toupper(tree_method) == "NJ") {
      dist_matrix <- ape::dist.dna(ape::as.DNAbin(phydata))
      pml <- ape::nj(dist_matrix)
    } else if (toupper(tree_method) == "UPGMA") {
      dist_matrix <- ape::dist.dna(ape::as.DNAbin(phydata))
      pml <- phangorn::upgma(dist_matrix)
    } else if (toupper(tree_method) == "MP") {
      # Maximum Parsimony (MP) method using ratchet algorithm
      treeRatchet <- phangorn::pratchet(phydata, trace = 0)
      treeRatchet <- phangorn::acctran(treeRatchet, phydata)
      treeRatchet <- ape::di2multi(treeRatchet)

      if(inherits(treeRatchet, "multiPhylo")){
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

    # Plot or return the tree based on the `show_tree` parameter
    if(show_tree == TRUE) {
      plot(pml)
      return(pml)
    }

    return(pml)

  } else if (tolower(seq.type) == "protein") {
    AAbin_list <- list()

    # Loop through each sequences file and convert it to AAbin format
    for (i in seq_along(seq.files)) {
      seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")
      aa_data <- adegenet::fasta2DNAbin(seq.file)
      AAbin_list[[i]] <- aa_data
    }

    # Protein sequences processing
    multidna <- methods::new("multidna", AAbin_list)
    phydata <- phangorn::as.phyDat(apex::concatenate(multidna))

    # Tree construction based on the selected method for protein sequences
    if (toupper(tree_method) == "ML") {
      modeltest <- phangorn::modelTest(phydata)
      pml <- phangorn::pml_bb(modeltest)
    } else if (toupper(tree_method) == "NJ") {
      dist_matrix <- ape::dist.aa(ape::as.AAbin(phydata))
      pml <- ape::nj(dist_matrix)
    } else if (toupper(tree_method) == "UPGMA") {
      dist_matrix <- ape::dist.aa(ape::as.AAbin(phydata))
      pml <- phangorn::upgma(dist_matrix)
    } else if (toupper(tree_method) == "MP") {
      # Maximum Parsimony (MP) method using ratchet algorithm
      treeRatchet <- phangorn::pratchet(phydata, trace = 0)
      treeRatchet <- phangorn::acctran(treeRatchet, phydata)
      treeRatchet <- ape::di2multi(treeRatchet)

      if(inherits(treeRatchet, "multiPhylo")){
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

    # Plot or return the tree based on the `show_tree` parameter
    if(show_tree == TRUE) {
      plot(pml)
      return(pml)
    }
    return(pml)
  } else {
    stop("Unknown sequence type, please choose 'DNA' or 'protein'.")
  }
}

