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
#' @param tree_method The method for constructing the phylogenetic tree. Options include:
#'        "ML" (Maximum Likelihood), "NJ" (Neighbor Joining), "UPGMA", "MP" (Maximum Parsimony),
#'        and "BI" (Bayesian Inference). Default is "ML".
#' @param model The evolutionary model to use for tree construction. For DNA sequences with ML method,
#'        options include "all", "JC", "F81", "K80", "HKY", "SYM", "GTR". For DNA sequences with NJ/UPGMA,
#'        options include "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92",
#'        "TN93", "GG95", "logdet", "paralin", "indel", "indelblock". For protein sequences with ML method,
#'        options include "all", "JTT", "WAG", "LG". Default is NULL (auto-selected).
#' @param show_tree Logical, if TRUE the constructed supertree is displayed. Default is TRUE.
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
#' tree1 <- concat_tree(seq.files = dna,
#'                      seq.type = "DNA",
#'                      data_dir = data_dir,
#'                      tree_method = "NJ")
#' tree2 <- concat_tree(seq.files = dna,
#'                      seq.type = "DNA",
#'                      data_dir = data_dir,
#'                      tree_method = "UPGMA")
#' }
#'
#' @export
concat_tree <- function(seq.files,
                        data_dir = NULL,
                        seq.type = "DNA",
                        tree_method = "ML",
                        model = NULL,
                        show_tree = TRUE) {
  # Process input type
  seq.type <- toupper(seq.type)

  # Check if the data directory is provided
  if(is.null(data_dir)){
    stop("Please input a valid 'data_dir'.")
  }

  # Handle DNA sequences case
  if (seq.type == "DNA") {
    DNAbin_list <- list()

    # Loop through each sequences file and convert it to DNAbin format
    for (i in seq_along(seq.files)) {
      seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")

      # Check if the sequence file exists
      if (!file.exists(seq.file)) {
        warning("File not found: ", seq.files[i])
        next
      }

      dna_data <- adegenet::fasta2DNAbin(seq.file)

      seq_chars <- as.character(dna_data)
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
      DNAbin_list[[i]] <- dna_data
    }

    # DNA sequences processing
    multidna <- methods::new("multidna", DNAbin_list)
    phydata <- phangorn::as.phyDat(apex::concatenate(multidna))

    # Tree construction based on the selected method for DNA sequences
    if (toupper(tree_method) == "ML") {
      valid_models <- c("all", "JC", "F81", "K80", "HKY", "SYM", "GTR")

      if(!is.null(model) && !model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
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
      dist_matrix <- ape::dist.dna(ape::as.DNAbin(phydata), model = model)
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
      dist_matrix <- ape::dist.dna(ape::as.DNAbin(phydata), model = model)
      tree <- phangorn::upgma(dist_matrix)
    } else if (toupper(tree_method) == "MP") {
      # Maximum Parsimony (MP) method using ratchet algorithm
      treeRatchet <- phangorn::pratchet(phydata, trace = 0)
      treeRatchet <- phangorn::acctran(treeRatchet, phydata)
      treeRatchet <- ape::di2multi(treeRatchet)

      if(inherits(treeRatchet, "multiPhylo")){
        treeRatchet <- unique(treeRatchet)
      }
      tree <- treeRatchet

    } else if (toupper(tree_method) == "BI") {
      # Bayesian Inference (BI) method using BEAST2
      if(length(seq.files) > 1) {
        stop("Bayesian Inference (BI) method currently only supports a single sequence file.")
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

      seq.file <- paste0(data_dir, "/", seq.files[1], ".fas")
      if (!file.exists(seq.file)) {
        stop("Sequence file not found: ", seq.file)
      }

      # Run BEAST2 and obtain the tree
      outputs <- babette::bbt_run_from_model(seq.file,
                                             inference_model = beautier::create_inference_model(),
                                             beast2_options = beastier::create_beast2_options())

      # Extract tree from outputs
      tree_names <- names(outputs)
      tree_pattern <- paste0(seq.files[1], "_trees")
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

    # Plot or return the tree based on the `show_tree` parameter
    if(show_tree == TRUE) {
      plot(tree)
      return(tree)
    }

    return(tree)

  } else if (seq.type == "PROTEIN") {
    AAbin_list <- list()

    # Loop through each sequences file and convert it to AAbin format
    for (i in seq_along(seq.files)) {
      seq.file <- paste0(data_dir, "/", seq.files[i], ".fas")

      # Check if the sequence file exists
      if (!file.exists(seq.file)) {
        warning("File not found: ", seq.files[i])
        next
      }

      phyDat_data <- phangorn::read.phyDat(seq.file, format = "fasta", type = "AA")
      aa_data <- ape::as.AAbin(phyDat_data)

      seq_chars <- as.character(aa_data)
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
      AAbin_list[[i]] <- aa_data
    }

    # Protein sequences processing
    multidna <- methods::new("multidna", AAbin_list)
    phydata <- phangorn::as.phyDat(apex::concatenate(multidna))

    # Tree construction based on the selected method for protein sequences
    if (toupper(tree_method) == "ML") {
      valid_models <- c("all", "JTT", "WAG", "LG")

      if(!is.null(model) && !model %in% valid_models) {
        stop("Invalid model. Please choose one of: ",
             paste(valid_models, collapse = ", "))
      }
      modeltest <- phangorn::modelTest(phydata, model = model)
      tree <- phangorn::pml_bb(modeltest)
    } else if (toupper(tree_method) == "NJ") {
      if(!is.null(model)){
        stop("For protein sequences, the NJ method currently has no selectable tree-building models.")
      }
      dist_matrix <- ape::dist.aa(ape::as.AAbin(phydata))
      tree <- ape::nj(dist_matrix)

    } else if (toupper(tree_method) == "UPGMA") {
      if(!is.null(model)){
        stop("For protein sequences, the UPGMA method currently has no selectable tree-building models.")
      }
      dist_matrix <- ape::dist.aa(ape::as.AAbin(phydata))
      tree <- phangorn::upgma(dist_matrix)
    } else if (toupper(tree_method) == "MP") {
      if(!is.null(model)){
        stop("No tree-building models are currently available for the Maximum Parsimony (MP) method.")
      }
      # Maximum Parsimony (MP) method using ratchet algorithm
      treeRatchet <- phangorn::pratchet(phydata, trace = 0)
      treeRatchet <- phangorn::acctran(treeRatchet, phydata)
      treeRatchet <- ape::di2multi(treeRatchet)

      if(inherits(treeRatchet, "multiPhylo")){
        treeRatchet <- unique(treeRatchet)
      }
      tree <- treeRatchet
    } else if (toupper(tree_method) == "BI") {
      stop("The Bayesian Inference (BI) method is not applicable for constructing phylogenetic trees from protein sequences.")
    } else {
      stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
    }

    # Plot or return the tree based on the `show_tree` parameter
    if(show_tree == TRUE) {
      plot(tree)
      return(tree)
    }
    return(tree)
  } else {
    stop("Unknown sequence type, please choose 'DNA' or 'protein'.")
  }
}
