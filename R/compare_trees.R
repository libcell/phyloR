#' Compare Two Phylogenetic Trees
#'
#' This function compares two phylogenetic trees by calculating their
#' similarity based on three different metrics: Robinson-Foulds (RF) similarity,
#' edge similarity, and subtree similarity. The comparisons are returned as a data frame.
#'
#' @param tree1 A phylogenetic tree of class 'phylo'.
#' @param tree2 A phylogenetic tree of class 'phylo'.
#'
#' @return A data frame containing the following similarity metrics:
#' \describe{
#'   \item{RF_similarity}{Robinson-Foulds similarity, a measure of tree topology similarity (range: 0 to 1).}
#'   \item{edge_similarity}{Edge similarity, representing the proportion of common edges between the trees (range: 0 to 1).}
#'   \item{subtree_similarity}{Subtree similarity, based on the proportion of identical subtrees between the trees (range: 0 to 1).}
#' }
#'
#' @details
#' - **Robinson-Foulds Similarity** is calculated by measuring the RF distance and converting it into a similarity score. The RF distance compares the trees based on their edge partitioning.
#' - **Edge Similarity** is determined by counting the common edges between the two trees and calculating the proportion of matching edges.
#' - **Subtree Similarity** is determined by comparing the subtrees rooted at internal nodes, counting how many of them are identical between the two trees.
#'
#'
#' @examples
#' # Compare two phylogenetic trees
#' tree1 <- ape::rtree(100)
#' tree2 <- phangorn::rSPR(tree1, 3)
#' x <- compare_trees(tree1, tree2)
#'
#' @importFrom phangorn RF.dist rSPR
#' @importFrom ape Ntip Nnode extract.clade rtree
#'
#' @export
compare_trees <- function(tree1, tree2) {
  # Ensure both trees are of class 'phylo'
  if (!inherits(tree1, "phylo") || !inherits(tree2, "phylo")) {
    stop("Both input trees must be of class 'phylo'.")
  }

  # Initialize a data frame to store comparison results
  relation <- data.frame(RF_similarity = NA,
                         edge_similarity = NA,
                         subtree_similarity = NA,
                         stringsAsFactors = FALSE)

  # Calculate Robinson-Foulds similarity
  # Robinson-Foulds distance using phangorn package, then convert to similarity (range 0, 1]
  max_RF <- max(length(tree1$edge), length(tree2$edge)) - 1  # Max possible RF distance
  RF_distance <- phangorn::RF.dist(tree1, tree2)
  relation$RF_similarity <- max(0, 1 - RF_distance / max_RF)  # Convert RF to similarity range [0, 1]

  # Function to calculate edge similarity based on common edges
  calculate_edge_similarity <- function(tree1, tree2) {
    # Convert edges to a unique identifier (string)
    tree1_edges <- apply(tree1$edge, 1, function(x) paste(x, collapse = "-"))
    tree2_edges <- apply(tree2$edge, 1, function(x) paste(x, collapse = "-"))

    # Find common edges between the two trees
    common_edges <- length(intersect(tree1_edges, tree2_edges))

    # The edge similarity is the proportion of common edges, mapped to range (0, 1]
    edge_similarity <- common_edges / max(length(tree1_edges), length(tree2_edges))
    return(max(0, edge_similarity))  # Map to [0, 1]
  }

  # Calculate edge similarity
  relation$edge_similarity <- calculate_edge_similarity(tree1, tree2)

  # Function to calculate subtree similarity based on identical subtrees
  calculate_subtree_similarity <- function(tree1, tree2) {
    subtree_common_count <- 0
    num_tips <- ape::Ntip(tree1)  # Number of tips in the tree
    num_nodes <- ape::Nnode(tree1)  # Number of internal nodes in the tree

    # Loop over the internal nodes and extract subtrees to compare
    for (i in (num_tips + 1):(num_tips + num_nodes)) {
      # Ensure the node is valid before extracting the subtree
      if (i <= length(tree1$edge)) {
        subtree1 <- ape::extract.clade(tree1, i)
        subtree2 <- ape::extract.clade(tree2, i)

        # Compare the subtrees
        if (identical(subtree1, subtree2)) {
          subtree_common_count <- subtree_common_count + 1
        }
      }
    }

    # Calculate subtree similarity as the proportion of common subtrees, mapped to range (0, 1]
    total_subtrees <- num_nodes  # Total internal nodes/subtrees to compare
    subtree_similarity <- subtree_common_count / total_subtrees
    return(max(0, subtree_similarity))  # Map to [0, 1]
  }

  # Calculate subtree similarity
  relation$subtree_similarity <- calculate_subtree_similarity(tree1, tree2)

  # Return the results as a data frame
  return(relation)
}

















