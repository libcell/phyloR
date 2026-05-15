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
#' @param tree_method The method for constructing the phylogenetic tree. Options include:
#'        "ML" (Maximum Likelihood), "NJ" (Neighbor Joining), "UPGMA", "MP" (Maximum Parsimony),
#'        and "BI" (Bayesian Inference). Default is "ML".
#' @param model The evolutionary model to use for tree construction. For DNA sequences with ML method,
#'        options include "all", "JC", "F81", "K80", "HKY", "SYM", "GTR". For DNA sequences with NJ/UPGMA,
#'        options include "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92",
#'        "TN93", "GG95", "logdet", "paralin", "indel", "indelblock". For DNA sequences with BI method,
#'        options include "JC69", "HKY", "TN93", "GTR". For protein sequences with ML method,
#'        options include "all", "JTT", "WAG", "LG". Default is NULL (auto-selected).
#' @param show_tree Logical, if TRUE the constructed supertree is displayed. Default is TRUE.
#'
#' @return A phylogenetic tree object (of class `phylo`), representing the supertree constructed from
#'         the individual trees built from the sequence files. If `show_tree` is TRUE, the tree is plotted.
#'
#' @details This function performs the following steps:
#'   1. Retrieves sequence data from the provided files and processes them based on the specified sequence type.
#'   2. Builds an individual phylogenetic tree for each sequence file using the chosen tree-building method.
#'   3. Combines the individual trees into a supertree using the superTree function.
#'   4. Optionally, displays the supertree if the `show_tree` parameter is TRUE.
#'   5. Returns the supertree object.
#'
#' @importFrom adegenet fasta2DNAbin
#' @importFrom ape dist.dna dist.aa nj di2multi prop.clades as.AAbin
#' @importFrom phangorn as.phyDat modelTest pml_bb upgma pratchet acctran superTree read.phyDat
#' @importFrom babette bbt_run_from_model
#' @importFrom beautier create_inference_model create_jc69_site_model create_hky_site_model
#'             create_tn93_site_model create_gtr_site_model
#' @importFrom beastierinstall install_beast2
#' @importFrom beastier create_beast2_options
#'
#' @examples
#' \dontrun{
#' dna <- c("K01939", "K03644", "K00797", "K00927", "K00088", "K02257", "K00164",
#'          "K00820", "K06158", "K00008")
#' data_dir <- system.file("extdata", "sequences", package = "phyloPipeR")
#' tree1 <- coalescent_tree(seq.files = dna,
#'                          seq.type = "DNA",
#'                          data_dir = data_dir,
#'                          tree_method = "NJ",
#'                          show_tree = TRUE)
#' tree2 <- coalescent_tree(seq.files = dna,
#'                          seq.type = "DNA",
#'                          data_dir = data_dir,
#'                          tree_method = "UPGMA")
#' }
#'
#' @export
coalescent_tree <- function(seq.files,
                            data_dir = NULL,
                            seq.type = "DNA",
                            tree_method = "ML",
                            model = NULL,
                            show_tree = TRUE) {

  # Check if the data directory is provided
  if(is.null(data_dir)){
    stop("Please input a valid 'data_dir'.")
  }

  # Process input type
  seq.type <- toupper(seq.type)

  mtree <- list()

  # Loop through sequence files
  for (i in seq_along(seq.files)) {
    seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")

    # Check if the sequence file exists
    if (!file.exists(seq.file)) {
      warning("File not found: ", seq.files[i])
      next
    }

    # DNA sequence processing
    if (seq.type == "DNA") {
      dna <- adegenet::fasta2DNAbin(seq.file)

      seq_chars <- as.character(dna)
      seq_strings <- apply(seq_chars, 1, paste, collapse = "")
      seq_counts <- table(seq_strings)
      duplicate_groups <- seq_counts[seq_counts > 1]

      # If duplicates exist, report and stop execution
      if(length(duplicate_groups) > 0) {
        err_msg <- paste("In file", seq.files[i], "found",
                         length(duplicate_groups), "group(s) of duplicate sequences:\n")
        for(j in 1:length(duplicate_groups)) {
          err_msg <- paste0(err_msg, "\nGroup ", j, ": ", duplicate_groups[j], " sequence(s)\n")
          group_seqs <- names(seq_strings)[seq_strings == names(duplicate_groups)[j]]
          err_msg <- paste0(err_msg, "  IDs: ", paste(group_seqs, collapse = ", "), "\n")
        }
        stop(err_msg)
      }

      # Tree construction based on the selected method for DNA sequences
      if (toupper(tree_method) == "ML") {
        valid_models <- c("all", "JC", "F81", "K80", "HKY", "SYM", "GTR")

        if(!is.null(model) && !model %in% valid_models) {
          stop("Invalid model. Please choose one of: ",
               paste(valid_models, collapse = ", "))
        }
        phydata <- phangorn::as.phyDat(dna)
        modeltest <- phangorn::modelTest(phydata, model = model)
        pml <- phangorn::pml_bb(modeltest)
        tree <- pml$tree
      } else if (toupper(tree_method) == "NJ") {
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
        dist_matrix <- ape::dist.dna(dna, model = model)
        tree <- ape::nj(dist_matrix)

      } else if (toupper(tree_method) == "UPGMA") {
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
        dist_matrix <- ape::dist.dna(dna, model = model)
        tree <- phangorn::upgma(dist_matrix)
      } else if (toupper(tree_method) == "MP") {
        if(!is.null(model)){
          stop("No tree-building models are currently available for the Maximum Parsimony (MP) method.")
        }
        # Maximum Parsimony (MP) method using ratchet algorithm
        phydata <- phangorn::as.phyDat(dna)
        treeRatchet <- phangorn::pratchet(phydata, trace = 0)
        treeRatchet <- phangorn::acctran(treeRatchet, phydata)
        treeRatchet <- ape::di2multi(treeRatchet)

        if (inherits(treeRatchet, "multiPhylo")) {
          treeRatchet <- unique(treeRatchet)
        }
        tree <- treeRatchet
      } else if (toupper(tree_method) == "BI") {
        # Bayesian Inference (BI) method using BEAST2
        valid_models <- c("JC69", "HKY", "TN93", "GTR")

        if(is.null(model)) {
          model <- "JC69"
        }
        if(!model %in% valid_models) {
          stop("Invalid model. Please choose one of: ",
               paste(valid_models, collapse = ", "))
        }

        # Check if BEAST2 is installed
        is_installed <- Sys.which("beast")
        if (!nzchar(is_installed)) {
          message("BEAST2 not found. Attempting to install...")
          tryCatch({
            beastierinstall::install_beast2()
            message("BEAST2 installation was successful!")
          }, error = function(e) {
            stop("Failed to install BEAST2: ", e$message)
          })
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
        tree_names <- names(outputs)
        tree_pattern <- paste0(seq.files[i], "_trees")
        tree_index <- grep(tree_pattern, tree_names)

        if(length(tree_index) == 0) {
          stop("Failed to extract tree from BEAST2 output.")
        }

        treeBI <- outputs[[tree_index[1]]]
        if(is.list(treeBI) && length(treeBI) > 0) {
          tree <- treeBI[[which.max(sapply(treeBI, function(x) attr(x, "posterior")))]]
        } else {
          tree <- treeBI
        }

        # Delete all .log and .trees files
        unlink(list.files(pattern = "\\.(log|trees)$"), recursive = TRUE)

      } else {
        stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
      }

      mtree[[i]] <- tree

    } else if (seq.type == "PROTEIN") {

      # Protein sequence processing
      protein <- phangorn::read.phyDat(seq.file, format = "fasta", type = "AA")
      protein_abin <- ape::as.AAbin(protein)

      seq_chars <- as.character(protein_abin)
      seq_strings <- apply(seq_chars, 1, paste, collapse = "")
      seq_counts <- table(seq_strings)
      duplicate_groups <- seq_counts[seq_counts > 1]

      # If duplicates exist, report and stop execution
      if(length(duplicate_groups) > 0) {
        err_msg <- paste("In file", seq.files[i], "found",
                         length(duplicate_groups), "group(s) of duplicate sequences:\n")
        for(j in 1:length(duplicate_groups)) {
          err_msg <- paste0(err_msg, "\nGroup ", j, ": ", duplicate_groups[j], " sequence(s)\n")
          group_seqs <- names(seq_strings)[seq_strings == names(duplicate_groups)[j]]
          err_msg <- paste0(err_msg, "  IDs: ", paste(group_seqs, collapse = ", "), "\n")
        }
        stop(err_msg)
      }

      # Tree construction based on the selected method for protein sequences
      if (toupper(tree_method) == "ML") {
        valid_models <- c("all", "JTT", "WAG", "LG")

        if(!is.null(model) && !model %in% valid_models) {
          stop("Invalid model. Please choose one of: ",
               paste(valid_models, collapse = ", "))
        }
        phydata <- phangorn::as.phyDat(protein_abin)
        modeltest <- phangorn::modelTest(phydata, model = model)
        pml <- phangorn::pml_bb(modeltest)
        tree <- pml$tree
      } else if (toupper(tree_method) == "NJ") {
        if(!is.null(model)){
          stop("For protein sequences, the NJ method currently has no selectable tree-building models.")
        }
        dist_matrix <- ape::dist.aa(protein_abin)
        tree <- ape::nj(dist_matrix)
      } else if (toupper(tree_method) == "UPGMA") {
        if(!is.null(model)){
          stop("For protein sequences, the UPGMA method currently has no selectable tree-building models.")
        }
        dist_matrix <- ape::dist.aa(protein_abin)
        tree <- phangorn::upgma(dist_matrix)
      } else if (toupper(tree_method) == "MP") {
        # Maximum Parsimony (MP) method using ratchet algorithm
        if(!is.null(model)){
          stop("No tree-building models are currently available for the Maximum Parsimony (MP) method.")
        }
        phydata <- phangorn::as.phyDat(protein_abin)
        treeRatchet <- phangorn::pratchet(phydata, trace = 0)
        treeRatchet <- phangorn::acctran(treeRatchet, phydata)
        treeRatchet <- ape::di2multi(treeRatchet)

        if (inherits(treeRatchet, "multiPhylo")) {
          treeRatchet <- unique(treeRatchet)
        }
        tree <- treeRatchet
      } else if (toupper(tree_method) == "BI") {
        # Bayesian Inference (BI) method using BEAST2
        stop("The Bayesian Inference (BI) method is not applicable for constructing phylogenetic trees from protein sequences.")
      } else {
        stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
      }

      mtree[[i]] <- tree
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
