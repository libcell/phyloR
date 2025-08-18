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
#' @importFrom ape dist.dna njs boot.phylo nj dist.aa drawSupportOnEdges di2multi prop.clades
#' @importFrom adegenet fasta2DNAbin
#' @importFrom phangorn as.phyDat pratchet acctran plotBS modelTest pml_bb bootstrap.pml upgma
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
#'                    tree_method = "NJ")
#'
#' @export
gene_tree <- function(seq.file,
                      seq.type = "DNA",
                      tree_method = "NJ",
                      show_tree = TRUE) {


  # Check if the sequence file exists
  if (!file.exists(seq.file)) {
    stop("File not found: ", seq.file)
  }

  # Read the sequence data based the sequence type
  if (!seq.type %in% c("DNA", "protein", "dna", "Protein")) {
    stop("Invalid sequence type. Choose either 'DNA' or 'protein'.")
  }
  if(toupper(seq.type) == "DNA"){
    seq <- adegenet::fasta2DNAbin(seq.file)
  }
  if(tolower(seq.type) == "protein"){
    seq <- adegenet::fasta2DNAbin(seq.file)
  }

  # --- Neighbor-Joining (NJ) method ---
  if(toupper(tree_method) == "NJ"){
    if(toupper(seq.type) == "DNA"){
      dist <- ape::dist.dna(seq)
      treenj <- ape::njs(dist)
      bstrees <- ape::boot.phylo(treenj, seq,
                                 FUN = function(x) ape::nj(ape::dist.dna(x)),
                                 trees = FALSE)
    }
    if(tolower(seq.type) == "protein"){
      dist <- ape::dist.aa(seq)
      treenj <- ape::njs(dist)
      bstrees <- ape::boot.phylo(treenj, seq,
                                 FUN = function(x) ape::nj(ape::dist.aa(x)),
                                 trees = FALSE)
    }

    if(show_tree == TRUE){
      plot(treenj)
      ape::drawSupportOnEdges(bstrees)
    }
    return(treenj)
  }

  # --- Maximum-Parsimony (MP) method ---
  if(toupper(tree_method) == "MP") {
    phydata <- phangorn::as.phyDat(seq)
    treeRatchet <- phangorn::pratchet(phydata, trace = 0)
    treeRatchet <- phangorn::acctran(treeRatchet, phydata)
    treeRatchet  <- ape::di2multi(treeRatchet)
    if(inherits(treeRatchet, "multiPhylo")){
      treeRatchet <- unique(treeRatchet)}
    if(show_tree == TRUE) {
      phangorn::plotBS(treeRatchet)
    }
    return(treeRatchet)
  }

  # --- Maximum-Likelihood (ML) method ---
  if(toupper(tree_method) == "ML") {
    phydata <- phangorn::as.phyDat(seq)
    if(toupper(seq.type) == "DNA"){
      modeltest <- phangorn::modelTest(phydata)
    }
    if(tolower(seq.type) == "protein"){
      modeltest <- phangorn::modelTest(phydata, model = "WAG")
    }
    pml <- phangorn::pml_bb(modeltest)
    bstrees <- phangorn::bootstrap.pml(pml)
    if(show_tree == TRUE) {
      phangorn::plotBS(pml$tree, bstrees)
    }
    return(pml$tree)
  }

  # --- Unweighted Pair Group Method with Arithmetic Mean (UPGMA) ---
  if(toupper(tree_method) == "UPGMA"){
    if(toupper(seq.type) == "DNA"){
      dist_matrix <- ape::dist.dna(seq, model = "JC69")
      upgma_tree <- phangorn::upgma(dist_matrix)
      bootstrap <- ape::boot.phylo(upgma_tree,
                                   seq,
                                   FUN = function(x) phangorn::upgma(ape::dist.dna(x)),
                                   trees = FALSE)
    }
    if(tolower(seq.type) == "protein"){
      dist_matrix <- ape::dist.aa(seq)
      upgma_tree <- phangorn::upgma(dist_matrix)
      bootstrap <- ape::boot.phylo(upgma_tree,
                                   seq,
                                   FUN = function(x) phangorn::upgma(ape::dist.aa(x)),
                                   trees = FALSE)
    }
    if(show_tree == TRUE){
      plot(upgma_tree)
      ape::drawSupportOnEdges(bootstrap)
    }
    return(upgma_tree)
  }

  # ---Bayesian-Inference (BI) method ---
  if(toupper(tree_method) == "BI") {
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
          message("Installation failed: ", e$message)
        })
      }
    }
    install_beast2_with_check()

    outputs <- babette::bbt_run_from_model(seq.file,
                                           inference_model = beautier::create_inference_model(),
                                           beast2_options = beastier::create_beast2_options())

    # Check the file format
    valid_extensions <- c("fas", "fasta")
    file_extension <- tools::file_ext(seq.file)
    if (!(file_extension %in% valid_extensions)) {
      stop(paste("Invalid file extension for", seq.file, " Only 'fas' and 'fasta' files are allowed."))
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
    if(show_tree == TRUE) {
      plot(treeBI)
      ape::drawSupportOnEdges(posterior)
    }
    return(treeBI)
  }
}
