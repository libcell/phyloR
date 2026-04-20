#' Construct Phylogenetic Tree from Species Data
#'
#' It generates a phylogenetic tree based on a list of species names or taxonomic IDs.
#' The tree is constructed using the taxonomy information retrieved from a specified database (e.g., NCBI).
#' It supports three input types for species: scientific names, taxonomic IDs, and abbreviated species names.
#'
#' @param species A character vector containing species names or taxonomic IDs. The input can be:
#'        - scientific names (e.g., "Homo sapiens"),
#'        - taxonomic IDs (numeric IDs),
#'        - abbreviated species names (e.g., "hsa" for Homo sapiens).
#' @param species.type A string indicating the type of input species. Options are:
#'        "scientificname" (default), "taxonomic_id", or "abbspname" (abbreviated species name).
#' @param db A string specifying the database to use for retrieving taxonomy information. Default is "ncbi".
#' @param show_tree Logical, if TRUE the constructed phylogenetic tree is plotted. Default is TRUE.
#'
#' @return A phylogenetic tree object of class `phylo` representing the tree constructed from the species data.
#'         If `show_tree` is TRUE, the tree is displayed.
#'
#' @details It retrieves species classification data from a database (e.g., NCBI) based on the species names or IDs
#' provided. It then constructs a phylogenetic tree using the taxonomy classification. The supported species input types
#' are scientific names, taxonomic IDs, and abbreviated species names.
#' @importFrom taxize classification class2tree
#' @examples
#'
#' # Example 1:
#' species1 <- c("Homo sapiens", "Pan troglodytes", "Mus musculus",
#'               "Rattus norvegicus","Canis lupus familiaris", "Felis catus")
#' tree1 <- species_tree(species = species1, species.type = "scientificname")
#'
#' # Example 2:
#' species2 <- c("9606", "9598", "10090", "9615", "9685", "10116")
#' tree2 <- species_tree(species = species2, species.type = "taxonomic_id")
#'
#' # Example 3:
#' species3 = c("ath", "gmx", "zma", "osa",
#'              "dme", "cel", "mmu", "rno",
#'              "hsa", "mcc", "ssc", "bta",
#'              "gga", "xla", "sce", "ece")
#' tree3 <- species_tree(species = species3, species.type = "abbspname")
#'
#'
#' @export
species_tree <- function(species,
                         species.type = "scientificname",
                         db = "ncbi",
                         show_tree = TRUE) {

  # Process based on species type
  if(species.type == "scientificname") {
    find.id <- function(a) {which(species_tbl[, 4] == a)}
    ids <- NULL
    for(i in 1:length(species)){
      position <- unlist(sapply(species[i], find.id))
      if(length(position) == 0){
        warning(paste(species[i], "No valid species found.", sep = ":"))
        position <- position
      }
      if(length(position) >= 1){
        position <- position[1]
      }
      id <- species_tbl[position, 2]
      ids <- c(ids, id)
    }
    classification <- taxize::classification(ids, db = db)
    tr <- taxize::class2tree(classification)
  }

  if(species.type == "taxonomic_id") {
    find.id <- function(a) {which(species_tbl[, 2] == a)}
    ids <- NULL
    for(i in 1:length(species)){
      position <- unlist(sapply(species[i], find.id))
      if(length(position) == 0){
        warning(paste(species[i], "No valid species found.", sep = ":"))
        position <- position
      }
      if(length(position) >= 1){
        position <- position[1]
      }
      id <- species_tbl[position, 2]
      ids <- c(ids, id)
    }
    classification <- taxize::classification(ids, db = db)
    tr <- taxize::class2tree(classification)
  }

  if(species.type == "abbspname"){
  find.id <- function(a) {which(species_tbl[, 3] == a)}
  ids <- NULL
  for(i in 1:length(species)){
    position <- unlist(sapply(species[i], find.id))
    if(length(position) == 0){
      warning(paste(species[i], "No valid species found.", sep = ":"))
      position <- position
    }
    if(length(position) >= 1){
      position <- position[1]
    }
    id <- species_tbl[position, 2]
    ids <- c(ids, id)
  }
    classification <- taxize::classification(ids, db = db)
    tr <- taxize::class2tree(classification)
  }

  # Error if species.type is invalid
  if (!(species.type %in% c("scientificname", "taxonomic_id", "abbspname"))) {
    warning("Please ensure that the species.type is 'scientificname', 'taxonomic_id' or 'abbspname'!")
  }
  phylo_tree <- tr$phylo
  if(show_tree == TRUE) {
    plot(phylo_tree)
    return(phylo_tree)}
  else{
    return(phylo_tree)
  }
}
