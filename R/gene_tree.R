#' Construct Phylogenetic Tree from Processed Single Gene Sequences File
#'
#' This function generates a phylogenetic tree from sequence data using various construction methods.
#' **Important**: The function returns an UNROOTED tree. Users MUST root the tree manually using
#' `ape::root()` (outgroup rooting, recommended) or `phangorn::midpoint()` (midpoint rooting)
#' before biological interpretation.
#'
#' @param seq.file A string specifying the path to the sequence file (FASTA format).
#' @param seq.type A string specifying the type of sequences. Options: "DNA" or "protein". Default "DNA".
#' @param tree_method Tree construction method. Options: "NJ" (Neighbor-Joining),
#'        "MP" (Maximum-Parsimony), "ML" (Maximum-Likelihood), "UPGMA", "BI" (Bayesian-Inference).
#'        Default "NJ".
#' @param model Character string specifying the evolutionary model. If NULL, defaults:
#'        \itemize{
#'          \item For NJ/UPGMA with DNA: "K80"
#'          \item For ML with DNA: automatically tests all available models and selects the best
#'          \item For ML with protein: automatically tests all available models and selects the best
#'          \item For BI with DNA: "JC69"
#'          \item For NJ/UPGMA with protein: no model selection available
#'          \item For MP: no model selection available
#'        }
#' @param show_tree Logical, whether to plot the tree. Default TRUE.
#'
#' @return A phylogenetic tree object of class `phylo` (UNROOTED). The tree should be rooted
#'         before further analysis using `ape::root()` (for outgroup rooting) or
#'         `phangorn::midpoint()` (for midpoint rooting).
#'
#' @details Supported methods:
#'   \itemize{
#'     \item \strong{NJ}: Neighbor-Joining, distance-based, fast. Suitable for large datasets.
#'     \item \strong{MP}: Maximum-Parsimony, uses parsimony ratchet algorithm. Finds the tree with
#'           the minimum number of evolutionary changes.
#'     \item \strong{ML}: Maximum-Likelihood, includes automatic model selection. Generally the most
#'           accurate but computationally intensive.
#'     \item \strong{UPGMA}: Hierarchical clustering based on distances. Assumes molecular clock.
#'     \item \strong{BI}: Bayesian Inference using BEAST2 (DNA only). Provides posterior probabilities
#'           for clade support.
#'   }
#'
#' @importFrom ape dist.dna njs dist.aa di2multi prop.clades root is.rooted keep.tip as.AAbin
#' @importFrom adegenet fasta2DNAbin
#' @importFrom phangorn as.phyDat pratchet acctran modelTest pml_bb upgma midpoint read.phyDat
#' @importFrom babette bbt_run_from_model
#' @importFrom beautier create_inference_model create_jc69_site_model create_hky_site_model
#'             create_tn93_site_model create_gtr_site_model
#' @importFrom beastier create_beast2_options
#' @importFrom tools file_ext
#'
#' @examples
#' # Example of constructing trees using different methods from a DNA sequence file
#' DNA_seq <- system.file("extdata", "DNA_seq.fas", package = "phyloPipeR")
#' tree1 <- gene_tree(seq.file = DNA_seq,
#'                    seq.type = "DNA",
#'                    tree_method = "UPGMA")
#'
#' # Example of constructing trees using different methods from a protein sequence file
#' protein_seq <- system.file("extdata", "protein_seq.fas", package = "phyloPipeR")
#' tree2 <- gene_tree(seq.file = protein_seq,
#'                    seq.type = "protein",
#'                    tree_method = "UPGMA")
#'
#' # Example with NJ method and custom model
#' tree3 <- gene_tree(seq.file = DNA_seq,
#'                    seq.type = "DNA",
#'                    tree_method = "NJ",
#'                    model = "TN93")
#'
#'
#' @export
gene_tree <- function(seq.file,
                      seq.type = "DNA",
                      tree_method = "NJ",
                      model = NULL,
                      show_tree = TRUE) {

  # Check if the sequence file exists
  if (!file.exists(seq.file)) {
    stop("File not found: ", seq.file)
  }

  # Process input type
  seq.type <- toupper(seq.type)

  tree <- NULL

  # Read the sequence data based the sequence type
  if (!seq.type %in% c("DNA", "PROTEIN")) {
    stop("Invalid sequence type. Choose either 'DNA' or 'protein'.")
  }

  if(seq.type == "DNA"){
    seq <- adegenet::fasta2DNAbin(seq.file)
  }
  if(seq.type == "PROTEIN"){
    seq <- phangorn::read.phyDat(seq.file, format = "fasta", type = "AA")
    seq <- ape::as.AAbin(seq)
  }

  # Check for duplicate sequences
  seq_chars <- as.character(seq)
  seq_strings <- apply(seq_chars, 1, paste, collapse = "")
  seq_counts <- table(seq_strings)
  duplicate_groups <- seq_counts[seq_counts > 1]

  # If duplicates exist, report and stop execution
  if(length(duplicate_groups) > 0) {
    err_msg <- paste("Found", length(duplicate_groups), "group(s) of duplicate sequences:\n")
    for(i in 1:length(duplicate_groups)) {
      err_msg <- paste0(err_msg, "\nGroup ", i, ": ", duplicate_groups[i], " sequence(s)\n")
      group_seqs <- names(seq_strings)[seq_strings == names(duplicate_groups)[i]]
      err_msg <- paste0(err_msg, "  IDs: ", paste(group_seqs, collapse = ", "), "\n")
    }
    stop(err_msg)
  }

  # --- Neighbor-Joining (NJ) method ---
  if(toupper(tree_method) == "NJ"){
    if(seq.type == "DNA") {
      valid_models <- c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81",
                        "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin",
                        "indel", "indelblock")

      if(is.null(model)) {
        model <- "K80"
      }
      if(!model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
      dist <- ape::dist.dna(seq, model = model)
      tree <- ape::njs(dist)
    }

    if(seq.type == "PROTEIN"){
      if(!is.null(model)){
        stop("For protein sequences, the NJ method currently has no selectable tree-building models.")
      }
      dist <- ape::dist.aa(seq)
      tree <- ape::njs(dist)
    }
  }

  # --- Maximum-Parsimony (MP) method ---
  else if(toupper(tree_method) == "MP") {
    if(!is.null(model)){
      stop("No tree-building models are currently available for the Maximum Parsimony (MP) method.")
    }
    phydata <- phangorn::as.phyDat(seq)
    treeRatchet <- phangorn::pratchet(phydata, trace = 0)
    treeRatchet <- phangorn::acctran(treeRatchet, phydata)
    treeRatchet <- ape::di2multi(treeRatchet)
    if(inherits(treeRatchet, "multiPhylo")){
      treeRatchet <- unique(treeRatchet)
    }
    tree <- treeRatchet
  }

  # --- Maximum-Likelihood (ML) method ---
  else if(toupper(tree_method) == "ML") {
    phydata <- phangorn::as.phyDat(seq)

    if(seq.type == "DNA"){
      valid_models <- c("all", "JC", "F81", "K80", "HKY", "SYM", "GTR")

      if(!is.null(model) && !model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
      modeltest <- phangorn::modelTest(phydata, model = model)
    } else if(seq.type == "PROTEIN"){
      valid_models <- c("all", "JTT", "WAG", "LG")

      if(!is.null(model) && !model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
      modeltest <- phangorn::modelTest(phydata, model = model)
    }
    pml <- phangorn::pml_bb(modeltest)
    tree <- pml$tree
  }

  # --- Unweighted Pair Group Method with Arithmetic Mean (UPGMA) ---
  else if(toupper(tree_method) == "UPGMA"){
    if(seq.type == "DNA"){
      valid_models <- c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81",
                        "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin",
                        "indel", "indelblock")

      if(is.null(model)) {
        model <- "K80"
      }
      if(!model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
      dist_matrix <- ape::dist.dna(seq, model = model, pairwise.deletion = TRUE)

      # Handle missing values
      dist_matrix[is.nan(dist_matrix)] <- 0
      dist_matrix[is.infinite(dist_matrix)] <- 0
      dist_matrix[is.na(dist_matrix)] <- 0

      tree <- phangorn::upgma(dist_matrix)
    } else if(seq.type == "PROTEIN"){
      if(!is.null(model)){
        stop("For protein sequences, the UPGMA method currently has no selectable tree-building models.")
      }
      dist_matrix <- ape::dist.aa(seq)

      # Handle missing values
      dist_matrix[is.na(dist_matrix)] <- 0
      dist_matrix[is.nan(dist_matrix)] <- 0
      dist_matrix[is.infinite(dist_matrix)] <- 0

      tree <- phangorn::upgma(dist_matrix)
    }
  }

  # --- Bayesian-Inference (BI) method ---
  else if(toupper(tree_method) == "BI") {
    if(seq.type == "PROTEIN"){
      stop("The Bayesian Inference (BI) method is not applicable for constructing phylogenetic trees from protein sequences.")
    }

    if(seq.type == "DNA"){
      # Check if BEAST2 is installed
      is_installed <- Sys.which("beast")
      if (!nzchar(is_installed)) {
        message("BEAST2 not found. Attempting to install...")
        if(requireNamespace("beastierinstall", quietly = TRUE)) {
          tryCatch({
            beastierinstall::install_beast2()
            message("BEAST2 installation was successful!")
          }, error = function(e) {
            stop("Failed to install BEAST2: ", e$message,
                 "\nPlease install BEAST2 manually from https://www.beast2.org/")
          })
        } else {
          stop("Please install BEAST2 manually from https://www.beast2.org/")
        }
      }

      valid_models <- c("JC69", "HKY", "TN93", "GTR")

      # Set default model if not specified
      if(is.null(model)) {
        model <- "JC69"
      }

      # Validate model
      if(!model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }

      # Create inference model based on selected model
      if(model == "JC69") {
        inference_model <- beautier::create_inference_model(
          site_model = beautier::create_jc69_site_model()
        )
      } else if(model == "HKY") {
        inference_model <- beautier::create_inference_model(
          site_model = beautier::create_hky_site_model()
        )
      } else if(model == "TN93") {
        inference_model <- beautier::create_inference_model(
          site_model = beautier::create_tn93_site_model()
        )
      } else if(model == "GTR") {
        inference_model <- beautier::create_inference_model(
          site_model = beautier::create_gtr_site_model()
        )
      }

      # Run BEAST2 and obtain the tree
      outputs <- babette::bbt_run_from_model(seq.file,
                                             inference_model = inference_model,
                                             beast2_options = beastier::create_beast2_options())

      # Extract tree from outputs
      file_name <- sub("\\.(fas|fasta)$", "", basename(seq.file))
      tree_pattern <- paste(file_name, "trees", sep = "_")

      if(!tree_pattern %in% names(outputs)) {
        stop("Failed to extract tree from BEAST2 output. Expected output name: ", tree_pattern)
      }

      treeBI <- outputs[[tree_pattern]]
      if(is.list(treeBI) && length(treeBI) > 0) {
        # Find the tree with maximum posterior probability
        posterior_vals <- sapply(treeBI, function(x) attr(x, "posterior"))
        tree <- treeBI[[which.max(posterior_vals)]]
      } else {
        tree <- treeBI
      }
    }
  } else {
    stop("Unknown tree method. Please choose 'NJ', 'MP', 'ML', 'UPGMA', or 'BI'.")
  }

  # Check if tree was successfully created
  if(is.null(tree)) {
    stop("Failed to construct tree. Please check your input parameters.")
  }

  # Delete BEAST2 generated files if they exist (only for BI method)
  if(toupper(tree_method) == "BI") {
    unlink(list.files(pattern = "\\.(log|trees)$"), recursive = TRUE)
  }

  # Plot the tree if requested
  if(show_tree == TRUE){
    plot(tree)
  }

  return(tree)
}
