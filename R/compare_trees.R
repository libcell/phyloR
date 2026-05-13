#' Compare Two Phylogenetic Trees
#'
#' It compares two phylogenetic trees by calculating their
#' similarity based on four different metrics: Robinson-Foulds (RF) similarity,
#' edge similarity, and entanglement value. The comparisons
#' are returned as a data frame containing the computed similarity values.
#'
#' @param tree1 A phylogenetic tree of class 'phylo', representing the first tree to compare.
#' @param tree2 A phylogenetic tree of class 'phylo', representing the second tree to compare.
#'
#' @return A data frame containing the following similarity metrics:
#' \describe{
#'   \item{RF_similarity}{Robinson-Foulds similarity, a measure of tree topology similarity (range: 0 to 1). A higher value indicates greater similarity in tree topology.}
#'   \item{edge_similarity}{Edge similarity, representing the proportion of common edges between the two trees (range: 0 to 1). A value of 1 indicates identical edge structures.}
#'   \item{entanglement_value}{Entanglement value, a metric quantifying the structural similarity of the trees, with higher values indicating greater similarity.}
#' }
#'
#' @details
#' - **Robinson-Foulds Similarity**: This metric is based on the Robinson-Foulds (RF) distance, which compares the trees by evaluating their edge partitions. The RF similarity score is computed by normalizing the RF distance with respect to the maximum possible RF distance between the two trees.
#' - **Edge Similarity**: This metric measures the proportion of edges that are common between the two trees. It is calculated by comparing the edges of both trees and counting how many are identical.
#' - **Entanglement Value**: This metric calculates the structural similarity of the trees based on their dendrogram representations. A higher entanglement value indicates a greater similarity in tree structure.
#'
#' @examples
#' # Compare two phylogenetic trees
#' DNA_seq <- system.file("extdata", "DNA_seq.fas", package = "phyloPipeR")
#' tree1 <- gene_tree(seq.file = DNA_seq,
#'                    seq.type = "DNA",
#'                    tree_method = "UPGMA")
#' species = c("ath", "gmx", "zma", "osa",
#'             "dme", "cel", "mmu", "rno",
#'             "hsa", "mcc", "ssc", "bta",
#'             "gga", "xla", "sce", "ece")
#' tree2 <- species_tree(species = species, species.type = "abbspname")
#' relation <- compare_trees(tree1, tree2)
#'
#' @importFrom phangorn RF.dist
#' @importFrom ape Ntip Nnode extract.clade node.depth.edgelength is.ultrametric
#' @importFrom stats as.dendrogram
#' @importFrom dendextend dendlist entanglement
#'
#' @export
compare_trees <- function(tree1, tree2) {

  # Ensure both trees are of class 'phylo'
  if (!inherits(tree1, "phylo") || !inherits(tree2, "phylo")) {
    stop("Both input trees must be of class 'phylo'.")
  }

  # Initialize result dataframe
  relation <- data.frame(RF_similarity = NA,
                         edge_similarity = NA,
                         entanglement_value = NA,
                         stringsAsFactors = FALSE)

  # Calculate Robinson-Foulds similarity
  max_RF <- max(length(tree1$edge), length(tree2$edge)) - 1
  RF_distance <- phangorn::RF.dist(tree1, tree2)
  relation$RF_similarity <- round(max(0, 1 - RF_distance / max_RF), 4)

  # Function to calculate edge similarity
  calculate_edge_similarity <- function(tree1, tree2) {
    # Convert edges to unique strings
    tree1_edges <- apply(tree1$edge, 1, function(x) paste(x, collapse = "-"))
    tree2_edges <- apply(tree2$edge, 1, function(x) paste(x, collapse = "-"))

    # Find common edges
    common_edges <- length(intersect(tree1_edges, tree2_edges))

    # Edge similarity as the proportion of common edges
    edge_similarity <- common_edges / max(length(tree1_edges), length(tree2_edges))
    return(max(0, edge_similarity))
  }

  # Calculate edge similarity
  relation$edge_similarity <- round(calculate_edge_similarity(tree1, tree2), 4)

  # Align tip nodes if trees are not ultrametric (equal tip depths)
  if(!ape::is.ultrametric(tree1)){
    align_tip_nodes <- function(phy) {
      aligned_tree <- phy
      node_depths <- ape::node.depth.edgelength(phy)
      tip_depths <- node_depths[1:length(phy$tip.label)]
      max_depth <- max(tip_depths)
      shortage <- max_depth - tip_depths

      # Store the node labels before adjusting tree
      original_node_labels <- phy$node.label

      # Adjust edge lengths to make tree ultrametric
      for (i in 1:length(phy$tip.label)) {
        edge_index <- which(phy$edge[, 2] == i)
        if (length(edge_index) > 0) {
          aligned_tree$edge.length[edge_index] <- phy$edge.length[edge_index] + shortage[i]
        }
      }

      # Restore the node labels after adjustment
      aligned_tree$node.label <- original_node_labels

      return(aligned_tree)
    }
    tree1 <- align_tip_nodes(tree1)
  }

  if(!ape::is.ultrametric(tree2)){
    align_tip_nodes <- function(phy) {
      aligned_tree <- phy
      node_depths <- ape::node.depth.edgelength(phy)
      tip_depths <- node_depths[1:length(phy$tip.label)]
      max_depth <- max(tip_depths)
      shortage <- max_depth - tip_depths

      # Store the node labels before adjusting tree
      original_node_labels <- phy$node.label

      # Adjust edge lengths to make tree ultrametric
      for (i in 1:length(phy$tip.label)) {
        edge_index <- which(phy$edge[, 2] == i)
        if (length(edge_index) > 0) {
          aligned_tree$edge.length[edge_index] <- phy$edge.length[edge_index] + shortage[i]
        }
      }

      # Restore the node labels after adjustment
      aligned_tree$node.label <- original_node_labels

      return(aligned_tree)
    }
    tree2 <- align_tip_nodes(tree2)
  }


  # Convert trees to dendrogram objects
  dend1 <- stats::as.dendrogram(tree1)
  dend2 <- stats::as.dendrogram(tree2)

  if(length(labels(dend1)) != length(tree1$tip.label)){
    stop("The tip.label of tree1 are lost!")
  }

  if(length(labels(dend2)) != length(tree2$tip.label)){
    stop("The tip.label of tree2 are lost!")
  }

  # Create dendlist and calculate entanglement value
  dend_list <- dendextend::dendlist(dend1, dend2)
  entanglement_value <- dendextend::entanglement(dend_list)

  # Store entanglement value in the result data frame
  relation$entanglement_value <- round(entanglement_value, 4)

  # Return the comparison results
  return(relation)
}

