#' Construct Phylogenetic Tree from Single Gene
#'
#' It retrieves orthologous gene for a list of species, performs sequence
#' alignment, and constructs a phylogenetic tree based on the concatenated sequences.
#' It supports different sequence types and alignment methods, with multiple tree construction options.
#'
#' @param gene_id A string specifying the gene identifier (e.g., ko_id).
#' @param id.type The type of gene ID used(e.g., ko_id, ncbi_id, or kegg_id), default is "ko_id" (KEGG Orthology ID).
#' @param species.list A character vector specifying the species to be included in the analysis.
#'        If NULL, all species available for each gene ID will be used.
#' @param species.type The type of species identifier to be used, Options are "scientificname", "taxonomic_id", and "abbspname".
#'        default is "scientificname".
#' @param seq.type The type of sequence to be retrieved, default is "DNA".
#'        Other options may include "protein".
#'               Options are "ClustalW" (default), "Muscle", or "ClustalOmega".
#' @param gapOpening The penalty for opening a gap in the alignment.
#'        Default is "default", which uses algorithm-specific default values.
#' @param gapExtension The penalty for extending an existing gap.
#'        Default is "default", which uses algorithm-specific default values.
#' @param maxiters The maximum number of refinement iterations.
#'        Default is "default", which uses algorithm-specific default values.
#' @param gap.end Fraction of gaps tolerated at the ends of the alignment (0-1). Default is 0.5.
#' @param gap.mid Fraction of gaps tolerated inside the alignment (0-1). Default is 0.9.
#' @param align_method The alignment method to be used, default is "ClustalW".
#'        Options are "ClustalW", "Muscle" or "ClustalOmega". Default is "ClustalW".
#' @param tree_method The method for constructing the phylogenetic tree. Options include:
#'        "ML" (Maximum Likelihood), "NJ" (Neighbor Joining), "UPGMA", "MP" (Maximum Parsimony),
#'        and "BI" (Bayesian Inference). Default is "ML".
#' @param model Character string specifying the evolutionary model. If NULL, defaults:
#'        \itemize{
#'          \item For NJ/UPGMA with DNA: "K80"
#'          \item For ML with DNA: automatically tests all available models and selects the best
#'          \item For ML with protein: automatically tests all available models and selects the best
#'          \item For BI with DNA: "JC69"
#'          \item For NJ/UPGMA with protein: no model selection available
#'          \item For MP: no model selection available
#'        }
#' @param show_tree Logical, if TRUE the constructed tree is displayed using a suitable tree plotting method.
#'        Default is TRUE.
#'
#' @return A phylogenetic tree object (of class `phylo`), representing the coalescent tree constructed from
#'         the aligned and trimmed sequences.
#'
#' @details It performs the following steps:
#'   1. Retrieves orthologous gene sequences for the provided gene IDs and species.
#'   2. Aligns and trims the sequences based on the specified alignment method.
#'   3. Constructs a phylogenetic tree using the processed sequence data and the selected tree-building method.
#'   4. Optionally, displays the tree if the `show_tree` parameter is set to TRUE.
#'   5. Cleans up temporary files and directories after tree construction.
#'
#' @importFrom KEGGREST keggGet
#' @importFrom rentrez entrez_summary
#' @importFrom pbapply pblapply
#' @importFrom KEGGREST keggGet
#' @importFrom rentrez entrez_summary
#' @importFrom msa msa msaConvert
#' @importFrom microseq readFasta msaTrim writeFasta
#' @importFrom bios2mds export.fasta
#' @importFrom Biostrings readDNAStringSet readAAStringSet
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
#' \dontrun{
#' # Example usage:
#' tree <- gene_tree_fetch(gene_id = "K00826",
#'                         species.list = c("hsa", "sce", "dme", "cel",
#'                                          "xla", "gga", "ssc", "rno",
#'                                          "mmu", "mcc", "gmx", "bta",
#'                                          "ece", "zma", "osa", "ath"),
#'                         species.type = "abbspname",
#'                         seq.type = "DNA",
#'                         tree_method = "NJ",
#'                         model = NULL,
#'                         show_tree = TRUE)
#' }
#'
#' @export
gene_tree_fetch <- function(gene_id,
                            id.type = "ko_id",
                            species.list,
                            species.type = "scientificname",
                            seq.type = "DNA",
                            align_method = "ClustalW",
                            gapOpening = "default",
                            gapExtension = "default",
                            maxiters = "default",
                            gap.end = 0.5,
                            gap.mid = 0.9,
                            tree_method = "NJ",
                            model = NULL,
                            show_tree = TRUE){

  # Create a new directory for sequences files
  temp_dir <- tempdir()
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }

  # Retrieve orthologous gene information for the provided species list
  species_info <- get_orthologs(gene_id = gene_id,
                                id.type = id.type,
                                species.list = species.list,
                                species.type = species.type)

  # Process species names and get the KEGG IDs
  species <- tolower(species_info[, 3])
  gene_ids <- paste(species, species_info[, 1], sep = ":")
  find.species <- function(a) {which(species_tbl[, 3] == a)}
  spnames <- NULL
  for(s in 1:length(species)){
    position <- unlist(sapply(species[s], find.species))
    if(length(position) == 0){
      warning(paste(species[s], "No valid species found.", sep = ":"))
      next
    }
    position <- position[1]
    spe <- species_tbl[position, 4]
    spnames <- c(spnames, spe)
  }

  # Retrieve sequences for KEGG IDs base on the sequence type
  seq <- get_kegg_sequences(gene_ids = gene_ids,
                            id.type = "kegg_id",
                            seq.type = seq.type)

  # Ensure spnames and seq lengths are equal
  if (length(spnames) != length(seq)) {
    stop("Error: Length of spnames does not match the length of seq.")
  }

  names(seq) <- spnames

  # Write the sequences to a FASTA file
  output_path <- file.path(temp_dir, "sequences.fasta")
  seqinr::write.fasta(seq, names = names(seq), file.out = output_path, nbchar = 60)

  # Prepare output file path for processed sequences
  if(toupper(seq.type) == "DNA"){
    file.name <- paste("nt", gene_id, sep = "_")
    data_file <- paste(file.name, "fasta", sep = ".")
  }

  if(toupper(seq.type) == "PROTEIN"){
    file.name <- paste("aa", gene_id, sep = "_")
    data_file <- paste(file.name, "fasta", sep = ".")
  }

  output_file <- file.path(temp_dir, data_file)

  # Align and trim the sequences using the selected alignment method
  processed_seq <- align_trim(seq.file = output_path,
                              seq.type = seq.type,
                              method = align_method,
                              gapOpening = gapOpening,
                              gapExtension = gapExtension,
                              maxiters = maxiters,
                              gap.end = gap.end,
                              gap.mid = gap.mid,
                              output_file = output_file)

  # Remove the temporary sequence file
  file.remove(output_path)

  # Construct a phylogenetic tree based on the processed sequences using the selected method
  tree <- gene_tree(seq.file = output_file,
                    seq.type = seq.type,
                    tree_method = tree_method,
                    model = model,
                    show_tree = show_tree)
  return(tree)

  # Clean up the directory
  unlink(temp_dir, recursive = TRUE)
}
