#' Align and Trim DNA or Protein Sequences
#'
#' Performs multiple sequence alignment using specified algorithms (ClustalW, MUSCLE, or ClustalOmega),
#' then trims the alignment by removing gap-rich regions from both ends and internal positions.
#' The function supports both DNA and protein sequences and can either return the trimmed alignment
#' as an R object or save it directly to a FASTA file.
#'
#' @param seq.file A character string specifying the path to the input sequence file.
#'        The file must be in FASTA format (extensions: .fas, .fasta, or .fa).
#' @param seq.type A character string specifying the type of sequences. Options are "DNA" (default) or "PROTEIN".
#'        Case-insensitive (e.g., "dna", "DNA", "protein" are all accepted).
#' @param method A character string specifying the alignment method to use.
#'        Options are "ClustalW" (default), "Muscle", or "ClustalOmega".
#' @param gapOpening The penalty for opening a gap in the alignment.
#'        Default is "default", which uses algorithm-specific default values.
#' @param gapExtension The penalty for extending an existing gap.
#'        Default is "default", which uses algorithm-specific default values.
#' @param maxiters The maximum number of refinement iterations.
#'        Default is "default", which uses algorithm-specific default values.
#' @param gap.end Fraction of gaps tolerated at the ends of the alignment (0-1). Default is 0.5.
#' @param gap.mid Fraction of gaps tolerated inside the alignment (0-1). Default is 0.9.
#' @param output_file A character string specifying the path to the output file where the
#'        trimmed alignment will be saved in FASTA format. If NULL (default), the function
#'        returns the trimmed alignment as an R object. If the file already exists, the user
#'        will be prompted to confirm overwriting.
#'
#' @return If `output_file` is NULL, returns a data frame (from `microseq::readFasta`) containing
#'        the trimmed alignment with two columns: `Header` (sequence names) and `Sequence`
#'        (trimmed sequence strings). If `output_file` is specified, the function writes the
#'        trimmed alignment to the file and returns nothing (invisible NULL).
#'
#' @details
#' The function internally uses `msa::msa()` for multiple sequence alignment. The alignment
#' result is then converted to a format compatible with `bios2mds::export.fasta()` via
#' `msa::msaConvert()`, temporarily written to disk as "alignment.fas", and read back by
#' `microseq::readFasta()` for trimming. The temporary file is automatically removed after use.
#'
#' Trimming is performed by `microseq::msaTrim()`, which removes:
#' \itemize{
#'   \item Columns from the ends where gap proportion exceeds `gap.end`
#'   \item Internal columns where gap proportion exceeds `gap.mid`
#' }
#'
#' @importFrom msa msa msaConvert
#' @importFrom microseq readFasta msaTrim writeFasta
#' @importFrom bios2mds export.fasta
#' @importFrom Biostrings readDNAStringSet readAAStringSet
#'
#' @examples
#' # Align and trim DNA sequences
#' DNA_seq <- system.file("extdata", "DNA_sequences.fas", package = "phyloPipeR")
#' align_trim(DNA_seq,
#'            seq.type = "DNA",
#'            method = "ClustalW",
#'            gapOpening = 10,
#'            gapExtension = 5)
#'
#' # Align and trim protein sequences
#' protein_seq <- system.file("extdata", "protein_sequences.fas", package = "phyloPipeR")
#' align_trim(protein_seq,
#'            seq.type = "PROTEIN",
#'            method = "Muscle",
#'            gap.end = 0.3,
#'            gap.mid = 0.95)
#'
#' @export
align_trim <- function(seq.file,
                       seq.type = "DNA",
                       method = "ClustalW",
                       gapOpening="default",
                       gapExtension="default",
                       maxiters="default",
                       gap.end = 0.5,
                       gap.mid = 0.9,
                       output_file = NULL){

  # Check if the sequence file exists
  if (!file.exists(seq.file)) {
    stop("File not found: ", seq.file)
  }

  # Additional check for file format (FASTA)
  file_ext <- tools::file_ext(seq.file)
  if (!(file_ext %in% c("fas", "fasta", "fa"))) {
    stop("The file must be in FASTA format. Invalid file extension: ", file_ext)
  }

  # Check if the sequence type is valid
  seq.type <- toupper(seq.type)
  if (!seq.type %in% c("DNA", "PROTEIN")) {
    stop("Invalid sequence type. Please choose 'DNA' or 'protein'.")
  }

  # Read the sequence file based on the sequence type
  if(seq.type == "DNA"){
    mySeqs <- Biostrings::readDNAStringSet(seq.file)
  } else if(seq.type == "PROTEIN"){
    mySeqs <- Biostrings::readAAStringSet(seq.file)
  }

  # Perform multiple sequence alignment
  alignment <- msa::msa(mySeqs,
                        method = method,
                        type = tolower(seq.type),
                        gapOpening = gapOpening,
                        gapExtension = gapExtension,
                        maxiters = maxiters)

  # Convert alignment to bios2mds format and export it to a temporary FASTA file
  tryCatch({
    alignment_set <- msa::msaConvert(alignment, type = "bios2mds::align")
    bios2mds::export.fasta(alignment_set,
                           outfile = "alignment.fas",
                           ncol = 60,
                           open = "w")
  }, error = function(e) {
    stop("Failed to convert the alignment to bios2mds format. Please check the alignment method.")
  })

  # Read the temporary file and trim the alignment
  aligned <- microseq::readFasta("alignment.fas")
  aln_trimmed <- microseq::msaTrim(aligned,
                                   gap.end = gap.end,
                                   gap.mid = gap.mid)

  # Remove the temporary FASTA file
  tryCatch({
    file.remove("alignment.fas")
  }, error = function(e) {
    warning("Failed to remove temporary file: alignment.fas. Please delete it manually.")
  })

  # Return processed sequences or save it to a file
  if(is.null(output_file)){
    return(aln_trimmed)
  }

  if(!is.null(output_file)){
    # Check if file exists and prompt for overwrite
    if (file.exists(output_file)) {
      overwrite <- readline(prompt = paste("File", output_file, "already exists. Do you want to overwrite? (y/n): "))
      if (tolower(overwrite) != "y") {
        stop("Operation cancelled. File not overwritten.")
      }
    }
    microseq::writeFasta(aln_trimmed, out.file = output_file)
  }
}


