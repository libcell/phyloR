#' Construct Phylogenetic Tree from Processed Single Gene Sequences File
#'
#' It generates a phylogenetic tree from sequence data using various tree construction methods,
#' including Neighbor-Joining (NJ), Maximum-Parsimony (MP), Maximum-Likelihood (ML), Unweighted Pair Group
#' Method with Arithmetic Mean (UPGMA), and Bayesian-Inference (BI). It supports both DNA and protein
#' sequences and allows for bootstrapping and tree visualization.
#'
#' @param seq.file A string specifying the path to the sequence file (FASTA format).
#' @param seq.type A string specifying the type of sequences in the file. Options are "DNA" or "protein". Default is "DNA".
#' @param tree_method A string specifying the tree construction method. Options include:
#'        "NJ" (Neighbor-Joining), "MP" (Maximum-Parsimony), "ML" (Maximum-Likelihood), "UPGMA", "BI" (Bayesian-Inference).
#'        Default is "NJ".
#' @param show_tree Logical, if TRUE the constructed phylogenetic tree is displayed. Default is TRUE.
#' @param root Character string or NULL, specify outgroup for rooting. If NULL, automatic outgroup selection is performed
#'
#' @return A phylogenetic tree object of class `phylo` representing the tree constructed using the selected method.
#'         If `show_tree` is TRUE, the tree is plotted with bootstrap support values.
#'
#' @details The function supports the following tree construction methods:
#'   - **Neighbor-Joining (NJ)**: Constructs a tree using the NJ algorithm, with options for DNA or protein sequences.
#'   - **Maximum-Parsimony (MP)**: Uses the MP algorithm for tree construction with parsimony optimization.
#'   - **Maximum-Likelihood (ML)**: Implements ML-based tree construction, including model selection and bootstrapping.
#'   - **Unweighted Pair Group Method with Arithmetic Mean (UPGMA)**: Generates a UPGMA tree based on pairwise distances.
#'   - **Bayesian-Inference (BI)**: Constructs a Bayesian phylogenetic tree using BEAST2 and supports posterior probabilities.
#'
#' @importFrom ape dist.dna njs dist.aa di2multi prop.clades root is.rooted keep.tip
#' @importFrom adegenet fasta2DNAbin
#' @importFrom phangorn as.phyDat pratchet acctran modelTest pml_bb upgma midpoint
#' @importFrom babette bbt_run_from_model
#' @importFrom beastierinstall install_beast2
#' @importFrom remotes install_github
#' @importFrom beautier create_inference_model
#' @importFrom beastier create_beast2_options
#'
#' @examples
#' # Example of constructing trees using different methods from a DNA sequence file
#' DNA_seq <- system.file("extdata", "DNA_seq.fas", package = "phyloR")
#' tree1 <- gene_tree(seq.file = DNA_seq,
#'                    seq.type = "DNA",
#'                    tree_method = "UPGMA")
#'
#' # Example of constructing trees using different methods from a protein sequence file
#' protein_seq <- system.file("extdata", "protein_seq.fas", package = "phyloR")
#' tree2 <- gene_tree(seq.file = protein_seq,
#'                    seq.type = "protein",
#'                    tree_method = "UPGMA")
#'
#' @export
gene_tree <- function(seq.file,
                      seq.type = "DNA",
                      tree_method = "NJ",
                      root = NULL,
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
    seq <- adegenet::fasta2DNAbin(seq.file)
  }

  # --- Neighbor-Joining (NJ) method ---
  if(toupper(tree_method) == "NJ"){
    if(seq.type == "DNA"){
      dist <- ape::dist.dna(seq)
      treenj <- ape::njs(dist)
      tree <- treenj
    }
    if(seq.type == "PROTEIN"){
      dist <- ape::dist.aa(seq)
      treenj <- ape::njs(dist)
      tree <- treenj
    }
  }

  # --- Maximum-Parsimony (MP) method ---
  else if(toupper(tree_method) == "MP") {
    phydata <- phangorn::as.phyDat(seq)
    treeRatchet <- phangorn::pratchet(phydata, trace = 0)
    treeRatchet <- phangorn::acctran(treeRatchet, phydata)
    treeRatchet  <- ape::di2multi(treeRatchet)
    if(inherits(treeRatchet, "multiPhylo")){
      treeRatchet <- unique(treeRatchet)}
    tree <- treeRatchet
  }

  # --- Maximum-Likelihood (ML) method ---
  else if(toupper(tree_method) == "ML") {
    phydata <- phangorn::as.phyDat(seq)
    if(seq.type == "DNA"){
      modeltest <- phangorn::modelTest(phydata)
    }
    if(seq.type == "PROTEIN"){
      modeltest <- phangorn::modelTest(phydata, model = "WAG")
    }
    pml <- phangorn::pml_bb(modeltest)
    tree <- pml$tree
  }

  # --- Unweighted Pair Group Method with Arithmetic Mean (UPGMA) ---
  else if(toupper(tree_method) == "UPGMA"){
    if(seq.type == "DNA"){
      dist_matrix <- ape::dist.dna(seq, model = "JC69")

      dist_matrix[is.na(dist_matrix)] <- 0
      dist_matrix[is.nan(dist_matrix)] <- 0
      dist_matrix[is.infinite(dist_matrix)] <- 0

      upgma_tree <- phangorn::upgma(dist_matrix)
    }
    if(seq.type == "PROTEIN"){
      dist_matrix <- ape::dist.aa(seq)

      dist_matrix[is.na(dist_matrix)] <- 0
      dist_matrix[is.nan(dist_matrix)] <- 0
      dist_matrix[is.infinite(dist_matrix)] <- 0

      upgma_tree <- phangorn::upgma(dist_matrix)

    }
    tree <- upgma_tree
  }

  # ---Bayesian-Inference (BI) method ---
  else if(toupper(tree_method) == "BI") {
    remotes::install_github("richelbilderbeek/beastierinstall")
    requireNamespace("beastierinstall")

    #Check for BEAST2 if not already installed
    install_beast2_with_check <- function() {
      is_installed <- Sys.which("beast")
      if (nzchar(is_installed)) {
        message("BEAST2 is already installed!")
      } else {
        tryCatch({
          beastierinstall::install_beast2()
          message("BEAST2 installation was successful!")
        }, error = function(e) {
          message("No installation required: ", e$message)
        })
      }
    }
    install_beast2_with_check()

    outputs <- babette::bbt_run_from_model(seq.file,
                                           inference_model = beautier::create_inference_model(),
                                           beast2_options = beastier::create_beast2_options())

    # Check the file format
    file_extension <- tools::file_ext(seq.file)
    if (!(file_extension %in% c("fas", "fasta"))) {
      stop("Invalid file extension.")
    }
    file_name <- sub("\\.(fas|fasta)$", "", basename(seq.file))
    tmp <- paste(file_name, "trees", sep = "_")

    treeBI <- outputs[[tmp]][[which.max(outputs$estimates$posterior)]]
    con.trees <- outputs[[tmp]][3:length(outputs[[tmp]])]
    posterior <- round(ape::prop.clades(treeBI, con.trees)/9999*100)
    tree <- treeBI

    # Delete all .log and .trees files
    unlink(list.files(pattern = "\\.(log|trees)$"), recursive = TRUE)

  }
  else {
    stop("Unknown tree method. Please choose 'ML', 'NJ', 'UPGMA', 'MP', or 'BI'.")
  }

  auto_select_outgroup <- function(tree, alignment) {

    # Convert alignment to matrix for distance calculation
    if (inherits(alignment, "DNAbin")) {
      alignment_matrix <- as.matrix(alignment)
    } else {
      alignment_matrix <- alignment
    }

    # Find common taxa between tree and alignment
    tree_tips <- tree$tip.label
    alignment_taxa <- rownames(alignment_matrix)

    if (!all(tree_tips %in% alignment_taxa)) {
      # Keep only common taxa
      common_taxa <- intersect(tree_tips, alignment_taxa)

      if (length(common_taxa) < 3) {
        # Not enough taxa for meaningful analysis - use midpoint
        warning("Insufficient overlap between tree and alignment. Using midpoint rooting.")
        return(phangorn::midpoint(tree))
      }

      # Subset tree and alignment to common taxa
      tree <- ape::keep.tip(tree, common_taxa)
      alignment_matrix <- alignment_matrix[common_taxa, , drop = FALSE]
    }

    # Calculate genetic distance matrix
    # Using K80 model (Kimura 2-parameter) as default
    tryCatch({
      dist_mat <- as.matrix(ape::dist.dna(alignment_matrix,
                                          model = "K80",
                                          pairwise.deletion = TRUE))
    }, error = function(e) {
      warning("Error calculating distance matrix: ", e$message, ". Using midpoint rooting.")
      return(phangorn::midpoint(tree))
    })

    # Handle NA values in distance matrix
    if (any(is.na(dist_mat))) {
      # Replace NA with 0 for calculation
      dist_mat[is.na(dist_mat)] <- 0
    }

    # Calculate average distance for each taxon (ignore self-distance)
    diag(dist_mat) <- NA
    avg_dist <- apply(dist_mat, 1, mean, na.rm = TRUE)

    # Check if all distances are NA
    if (all(is.na(avg_dist)) || length(avg_dist) == 0) {
      warning("Unable to calculate meaningful distances. Using midpoint rooting.")
      return(phangorn::midpoint(tree))
    }

    # ===== FIX: Handle multiple taxa with same maximum distance =====
    # Find the maximum distance value
    max_distance <- max(avg_dist, na.rm = TRUE)

    # Find ALL taxa with this maximum distance
    candidates <- names(avg_dist)[avg_dist == max_distance]

    # If multiple candidates, select the first one
    if (length(candidates) > 1) {
      warning("Multiple taxa have the same maximum average distance (",
              round(max_distance, 4), "). Selecting: ", candidates[1])
      outgroup_candidate <- candidates[1]
    } else {
      outgroup_candidate <- candidates
    }

    # Calculate distances between candidate and ingroup
    ingroup_taxa <- setdiff(names(avg_dist), outgroup_candidate)

    if (length(ingroup_taxa) < 2) {
      # Not enough taxa in ingroup
      warning("Too few taxa for outgroup analysis. Using midpoint rooting.")
      return(phangorn::midpoint(tree))
    }

    # Average distance from candidate to ingroup
    candidate_to_ingroup <- mean(dist_mat[ingroup_taxa, outgroup_candidate], na.rm = TRUE)

    # Average internal distance within ingroup
    ingroup_mat <- dist_mat[ingroup_taxa, ingroup_taxa]
    diag(ingroup_mat) <- NA
    ingroup_internal <- mean(ingroup_mat, na.rm = TRUE)

    # Check if candidate is sufficiently distant
    # Condition: outgroup distance > 1.5 × internal ingroup distance
    if (is.finite(candidate_to_ingroup) &&
        is.finite(ingroup_internal) &&
        candidate_to_ingroup > ingroup_internal * 1.5) {
      # Valid outgroup found
      return(outgroup_candidate)
    } else {
      # Candidate not sufficiently distant - use midpoint rooting
      return(phangorn::midpoint(tree))
    }
  }
  if(!is.null(root)){
    tree <- ape::root(tree, root, resolve.root = TRUE)
  }
  if(is.null(root) && !ape::is.rooted(tree)){
    outgroup_result <- auto_select_outgroup(tree, seq)

    # Safety check: ensure tree has tip labels
    if(is.null(tree$tip.label) || length(tree$tip.label) == 0) {
      stop("Tree has no tip labels. Cannot root tree.")
    }

    if (is.character(outgroup_result) && length(outgroup_result) > 0) {
      # Extract and clean outgroup name
      outgroup <- outgroup_result[1]
      outgroup <- as.character(trimws(outgroup))

      # Validate outgroup name
      if(is.na(outgroup) || outgroup == "" || nchar(outgroup) == 0) {
        message("Invalid outgroup name. Using midpoint rooting.")
        tree <- phangorn::midpoint(tree)
      } else if (!outgroup %in% tree$tip.label) {
        message("Outgroup '", outgroup, "' not found in tree tips. Using midpoint rooting.")
        tree <- phangorn::midpoint(tree)
      } else {
        message("Using automatically selected outgroup: ", outgroup)

        # Safely call ape::root with error handling
        tree_try <- tryCatch({
          rooted_tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)

          # Validate rooted tree
          if(is.null(rooted_tree$tip.label) || length(rooted_tree$tip.label) == 0) {
            stop("Rooted tree has no tips")
          }
          rooted_tree
        }, error = function(e) {
          message("Error rooting tree with outgroup '", outgroup, "': ", conditionMessage(e))
          return(NULL)
        })

        if(!is.null(tree_try)) {
          tree <- tree_try
        } else {
          message("Using midpoint rooting instead.")
          tree <- phangorn::midpoint(tree)
        }
      }
    } else if (inherits(outgroup_result, "phylo")) {
      # auto_select_outgroup returned a midpoint-rooted tree
      message("Using midpoint rooted tree from auto_select_outgroup")
      tree <- outgroup_result
    } else {
      # No suitable outgroup found
      message("No suitable outgroup found. Using midpoint rooting.")
      tree <- phangorn::midpoint(tree)
    }
  }
  if(show_tree == TRUE){
    plot(tree)
  }
  return(tree)
}
