#' Retrieve Gene Sequences from KEGG Database
#'
#' It retrieves DNA or protein sequences for a given list of gene identifiers
#' from the KEGG database based on the provided gene IDs, species, and sequence type.
#' It can handle different types of gene identifiers and allows filtering by species.
#'
#' @param gene_ids A character vector of gene identifiers (e.g., KO IDs, NCBI gene IDs, or KEGG gene IDs).
#' @param id.type A character string indicating the type of the provided gene IDs.
#'                Options are "KO_id", "ncbi_id", or "kegg_id" (default is "KO_id").
#' @param seq.type A character string specifying the type of sequence to retrieve.
#'                 Options are "DNA" or "protein" (default is "protein").
#' @param species.list A character vector of species names to filter results by species.Normally not NULL only if gene_ids is "KO_id".
#'                     If NULL, no filtering is applied (default is NULL).
#' @param species.type A character string specifying the type of species identifier.
#'                     Options are "scientificname", "taxonomic_id", or "abbspname" (default is "scientificname").
#'
#' @return A list of gene sequences (either DNA or protein) corresponding to the input gene IDs.
#'         If `seq.type` is "DNA", it returns nucleotide sequences, otherwise protein sequences.
#'
#' @importFrom KEGGREST keggGet
#' @importFrom rentrez entrez_summary
#' @importFrom pbapply pblapply
#'
#' @examples
#' # Example 1: Retrieve protein sequences for a list of KO IDs
#' seq_set1 <- get_kegg_sequences(gene_ids = "K00161",
#'                                id.type = "KO_id",
#'                                seq.type = "DNA",
#'                                species.list = c("hsa", "ptr"),
#'                                species.type = "abbspname")
#'
#' # Example 2: Retrieve DNA sequences for a list of NCBI gene IDs
#' seq_set2 <- get_kegg_sequences(gene_ids = c("5160", "465525"),
#'                                id.type = "ncbi_id")
#'
#' # Example 3: Retrieve protein sequences for a list of KEGG IDS
#' seq_set3 <- get_kegg_sequences(gene_ids = c("hsa:5160", "ptr:465525"),
#'                                id.type = "kegg_id",
#'                                seq.type = "DNA")
#'
#'
#' @export
get_kegg_sequences <- function(gene_ids,
                               id.type = "KO_id",
                               seq.type = "protein",
                               species.list = NULL,
                               species.type = "scientificname") {

  # Helper function to select species from a species table based on species name
  select.species <- function(ortholog.list, species, species.type = "scientificname") {

    # Search for species based on the species type
    if(species.type == "scientificname"){
      find.sciname <- function(a) {which(species_tbl[, 4] == a)}
      position <- unlist(sapply(species, find.sciname))
      spname <- species_tbl[position, 3]
      abbspname <- toupper(spname)
      find.species <- function(a) {which(ortholog.list[, 3] == a)}
      res <- NULL
      for (i in 1:length(abbspname)) {
        number <- unlist(sapply(abbspname[i], find.species))
        if (length(number) == 0) {
          if(length(grep(abbspname[i], toupper(species_tbl[, 3]))) == 0){
            warning(paste(species[i], "Please ensure that the species name is correct!", sep = ":"))
          }
          else{
            warning(paste(species[i], "No valid species found in the current gene.", sep = ":"))
          }
          return(NULL)
        }
        if(length(number >= 1)){
          number <- number[1]
        }
        orglist <- ortholog.list[number, ]
        res <- rbind(res, orglist)
      }
      return(res)
    }
    if(species.type == "taxonomic_id"){
      find.id <- function(a) {which(species_tbl[, 2] == a)}
      position <- as.vector(sapply(species, find.id))
      spname <- species_tbl[position, 3]
      abbspname <- toupper(spname)
      find.species <- function(a) {which(ortholog.list[, 3] == a)}
      res <- NULL
      for (i in 1:length(abbspname)) {
        number <- unlist(sapply(abbspname[i], find.species))
        if (length(number) == 0) {
          if(length(grep(abbspname[i], toupper(species_tbl[, 3]))) == 0){
            warning(paste(species[i], "Please ensure that the species name is correct!", sep = ":"))
          }
          else{
            warning(paste(species[i], "No valid species found in the current gene.", sep = ":"))
          }
          return(NULL)
        }
        if(length(number >= 1)){
          number <- number[1]
        }
        orglist <- ortholog.list[number, ]
        res <- rbind(res, orglist)
      }
      return(res)
    }
    if(species.type == "abbspname"){
      abbspname <- toupper(species)
      find.species <- function(a) {which(ortholog.list[, 3] == a)}
      res <- NULL
      for (i in 1:length(abbspname)) {
        number <- unlist(sapply(abbspname[i], find.species))
        if (length(number) == 0) {
          if(length(grep(abbspname[i], toupper(species_tbl[, 3]))) == 0){
            warning(paste(species[i], "Please ensure that the species name is correct!", sep = ":"))
          }
          else{
            warning(paste(species[i], "No valid species found in the current gene.", sep = ":"))
          }
          return(NULL)
        }
        if(length(number >= 1)){
          number <- number[1]
        }
        orglist <- ortholog.list[number, ]
        res <- rbind(res, orglist)      }
      return(res)
    }
    if (!(species.type %in% c('scientificname', 'taxonomic_id', 'abbspname'))) {
      warning("Please ensure that the species.type is 'scientificname', 'taxonomic_id' or 'abbspname'!")
    }
  }

  # KO_id handling
  if(id.type == "KO_id"){
    gene_ids <- toupper(gene_ids)

    # Check if gene_ids are valid in KO_list
    if(length(grep(gene_ids, ko_ids)) != 0){
      gene_ids <- gene_ids
    } else{
      stop("Please ensure that the gene_ids is correct!")
    }

    # Retrieve gene data from KEGG
    res <- KEGGREST::keggGet(gene_ids)
    genes <- res[[1]]$GENES
    x <- NULL

    # Use pblapply to add progress bar for processing genes
    x <- pbapply::pblapply(1:length(genes), function(s) {
      species <- strsplit(genes[s], ":")[[1]][1]
      entrez_ids <- strsplit(genes[s], ":")[[1]][2]
      id_seq <- strsplit(entrez_ids, " ")[[1]]
      id_seq <- id_seq[-1]
      tmp_list <- lapply(1:length(id_seq), function(i) {
        id <- id_seq[i]
        tmp <- strsplit(id, "\\(")[[1]]
        if (length(tmp) > 2) {
          new_tmp <- NULL
          new_tmp[1] <- tmp[1]
          new_tmp[2] <- paste(tmp[-1], collapse = "(")
          tmp <- new_tmp
        }
        tmp <- gsub(")$", "", tmp)
        if (length(tmp) == 1) tmp <- c(tmp, NA)
        tmp <- c(tmp, species)
        return(tmp)
      })
      do.call(rbind, tmp_list)
    })
    x <- do.call(rbind, x)
    colnames(x) <- c("Entrez_ID", "Gene_Symbol", "Species")
    rownames(x) <- NULL

    # If species list is provided, filter results by species
    if(!is.null(species.list)){
      select.list <- select.species(x, species = species.list, species.type = species.type)
    }
    if(is.null(species.list)){
      select.list <- x
    }

    # Build KEGG ID from species and gene
    kegg_id <- paste(tolower(select.list[, 3]), select.list[, 1], sep = ":")

    # Retrieve the sequence (DNA or protein)
    if(toupper(seq.type) == "DNA"){
      ntseq <- NULL
      for (i in 1:length(kegg_id)) {
        seq <- KEGGREST::keggGet(kegg_id[i], option = "ntseq")
        ntseq <- c(ntseq, seq)
        Sys.sleep(0.34)
      }
      return(ntseq)
    }
    if(tolower(seq.type) == "protein"){
      aaseq <- NULL
      for (i in 1:length(kegg_id)) {
        seq <- KEGGREST::keggGet(kegg_id[i], option = "aaseq")
        aaseq <- c(aaseq, seq)
        Sys.sleep(0.34)
      }
      return(aaseq)
    }
    if(!(seq.type %in% c("DNA", "protein", "dna", "Protein"))){
      stop("Unknown sequence type, please choose 'DNA' or 'protein'")
    }
  }

  # NCBI gene ID handling
  if(id.type == "ncbi_id"){
    ntseq <- NULL
    aaseq <- NULL
    seq <- NULL
    seq <- pbapply::pblapply(1:length(gene_ids), function(i) {
      summary <- rentrez::entrez_summary(db = "gene", id = gene_ids[i])
      otheraliases <- summary$otheraliases
      spname <- summary$organism$scientificname
      find.sciname <- function(a) {which(species_tbl[, 4] == a)}
      position <- as.vector(sapply(spname, find.sciname))
      abbspname <- species_tbl[position, 3]
      kegg_id <- NULL
      if (!is.null(otheraliases) && length(otheraliases) > 0) {
        locus_tag <- strsplit(otheraliases, ",")[[1]][1]
        kegg_id1 <- paste(abbspname, gene_ids[i], sep = ":")
        kegg_id2 <- paste(abbspname, locus_tag, sep = ":")

        result <- tryCatch({
          res <- KEGGREST::keggGet(kegg_id1)
        }, error = function(e) {
          message("Error with kegg_id1: ", e$message)
          return(NULL)
        })
        if (is.null(result)) {
          result <- tryCatch({
            res <- KEGGREST::keggGet(kegg_id2)
          }, error = function(e) {
            message("Error with kegg_id2: ", e$message)
            return(NULL)
          })
        }
        if (!is.null(result)) {
          kegg_id <- kegg_id1
        } else {
          kegg_id <- kegg_id2
        }
      }
      if(toupper(seq.type) == "DNA"){
        ntseq <- KEGGREST::keggGet(kegg_id, option = "ntseq")
        return(ntseq)
      }
      if(tolower(seq.type) == "protein"){
        aaseq <- KEGGREST::keggGet(kegg_id, option = "aaseq")
        return(aaseq)
      }
      if(!(seq.type %in% c("DNA", "protein", "dna", "Protein"))){
        stop("Unknown sequence type, please choose 'DNA' or 'protein'")
      }
    }, cl = 2)
    seq <- unlist(seq)
    return(seq)
  }

  # KEGG ID handling
  if(id.type == "kegg_id"){
    ntseq <- NULL
    aaseq <- NULL
    seq <- NULL
    seq <- pbapply::pblapply(1:length(gene_ids), function(i) {
      if(toupper(seq.type) == "DNA"){
        ntseq <- KEGGREST::keggGet(gene_ids[i], option = "ntseq")
        return(ntseq)
      }
      if(tolower(seq.type) == "protein"){
        aaseq <- KEGGREST::keggGet(gene_ids[i], option = "aaseq")
        return(aaseq)
      }
      if(!(seq.type %in% c("DNA", "protein", "dna", "Protein"))){
        stop("Unknown sequence type, please choose 'DNA' or 'protein'")
      }
    }, cl = 2)
    seq <- unlist(seq)
  }
  if(!(id.type) %in% c("KO_id", "ncbi_id", "kegg_id")) {
    warning("Please ensure that id.type is 'KO_id', 'ncbi_id' or 'kegg_id'!")
  }
  return(seq)
}
