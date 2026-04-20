#' Get Orthologs for a Gene
#'
#' It retrieves orthologous genes for a given gene ID. The function supports different gene ID types (ko_id, ncbi_id, kegg_id) and allows filtering by species.
#'
#' @param gene_id A character string specifying the gene identifier (e.g., ko_id, ncbi_id, or kegg_id).
#' @param id.type A character string specifying the type of the gene identifier. Options are "ko_id", "ncbi_id", and "kegg_id". Default is "ko_id".
#' @param species.list A character vector or string specifying the species or species identifiers to filter the results. Default is NULL.
#' @param species.type A character string specifying the type of species identifier. Options are "scientificname", "taxonomic_id", and "abbspname". Default is "scientificname".
#'
#' @return A data frame containing the orthologous genes with columns for Entrez ID, Gene Symbol, and Species. If species is specified, only the orthologs for that species will be returned.
#'
#' @importFrom KEGGREST keggGet
#' @importFrom rentrez entrez_summary
#' @importFrom pbapply pblapply
#'
#' @examples
#' taxinfo1 <- get_orthologs(gene_id = "K00161",
#'                           id.type = "ko_id",
#'                           species.list = "Homo sapiens",
#'                           species.type = "scientificname")
#
#' taxinfo2 <- get_orthologs(gene_id = "5160",
#'                           id.type = "ncbi_id",
#'                           species.list = "hsa",
#'                           species.type = "abbspname")
#'
#' taxinfo3 <- get_orthologs(gene_id = "hsa:5160",
#'                           id.type = "kegg_id",
#'                           species.list = "9606",
#'                           species.type = "taxonomic_id")
#'
#' @export
get_orthologs <- function(gene_id,
                          id.type = "ko_id",
                          species.list = NULL,
                          species.type = "scientificname") {

  # Function to select species from the ortholog list based on given species names or identifiers
  select.species <- function(ortholog.list, species, species.type = "scientificname"){
    species.type <- tolower(species.type)
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

  # Process input type
  id.type <- toupper(id.type)

  # If the gene_id is ko_id
  if(id.type == "KO_ID") {
    gene_id <- as.character(toupper(gene_id))

    # Check if gene_id are valid in KO_list
    if(length(grep(gene_id, ko_ids)) != 0) {
      gene_id <- gene_id
    } else {
      stop("Please ensure that the gene_id is correct!")
    }

    # Retrieve gene data from KEGG
    res <- KEGGREST::keggGet(gene_id)
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
    if(!is.null(species.list)) {
      select.list <- select.species(x, species.list, species.type = species.type)
      return(select.list)
    }
    if(is.null(species.list)) {
      return(x)
    }
  }

  # If the gene_id is ncbi_id
  if (id.type == "NCBI_ID") {
    res <- tryCatch({
      rentrez::entrez_summary(db = "gene", id = gene_id)
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (res == FALSE) {
      warning(paste0("No valid gene_id found: ", gene_id))
    }
    summary <- rentrez::entrez_summary(db = "gene", id = gene_id)
    otheraliases <- summary$otheraliases
    spname <- summary$organism$scientificname
    find.sciname <- function(a) {which(species_tbl[, 4] == a)}
    position <- as.vector(sapply(spname, find.sciname))
    abbspname <- species_tbl[position, 3]
    kegg_id <- NULL
    if (!is.null(otheraliases) && length(otheraliases) > 0) {
      locus_tag <- strsplit(otheraliases, ",")[[1]][1]
      kegg_id1 <- paste(abbspname, gene_id, sep = ":")
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
      gene_infor <- KEGGREST::keggGet(kegg_id)
      ko <- gene_infor[[1]]$ORTHOLOGY
      ko_id <- names(ko)
      res <- KEGGREST::keggGet(ko_id)
      genes <- res[[1]]$GENES
      x <- NULL
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
      if(!is.null(species.list)) {
        select.list <- select.species(x, species.list, species.type = species.type)
        return(select.list)
      }
      if(is.null(species.list)) {
        return(x)
      }
    }
  }

  # If the gene_id is kegg_id
  if(id.type == "KEGG_ID") {
    res <- tryCatch({
      KEGGREST::keggGet(gene_id)
      TRUE
      }, error = function(e) {
      FALSE
    })
    if (res == FALSE) {
      warning(paste0("No valid gene_id found: ", gene_id))
    }
    gene_infor <- KEGGREST::keggGet(gene_id)
    ko <- gene_infor[[1]]$ORTHOLOGY
    ko_id <- names(ko)
    res <- KEGGREST::keggGet(ko_id)
    genes <- res[[1]]$GENES
    x <- NULL
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
    if(!is.null(species.list)) {
      select.list <- select.species(x, species.list, species.type = species.type)
      return(select.list)
    }
    if(is.null(species.list)) {
      return(x)
    }
  }

  # Return an error if the id.type is not valid
  if(!(id.type) %in% c("KO_ID", "NCBI_ID", "KEGG_ID")) {
    warning("Please ensure that id.type is 'ko_id', 'ncbi_id' or 'kegg_id'!")
  }
}

