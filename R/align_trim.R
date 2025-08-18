#' Align and Trim DNA or Protein Sequences
#'
#' It reads DNA or protein sequences from a file, aligns them using a selected alignment method, trims the alignment,
#' and optionally saves the trimmed alignment to a file.
#'
#' @param seq.file A character string specifying the path to the input sequence file (FASTA format).
#' @param seq.type A character string specifying the type of sequences. Options are "DNA" or "protein".
#'        Default is "DNA".
#' @param method A character string specifying the alignment method to use. Options are "ClustalW", "Muscle", "Mafft" or "T-coffee".
#'        Default is "ClustalW".
#' @param output_file A character string specifying the path to the output file where the trimmed alignment will be saved.
#'        If NULL (default), the function returns the trimmed alignment as a sequence object.
#'
#' @return If `output_file` is NULL, returns the trimmed alignment as an object of class `DNAMultipleAlignment` or `AAMultipleAlignment`,
#'        depending on the sequence type. Otherwise, the trimmed alignment is written to the specified output file.
#'
#' @importFrom msa msa msaConvert
#' @importFrom microseq readFasta msaTrim writeFasta
#' @importFrom bios2mds export.fasta
#' @importFrom Biostrings readDNAStringSet readAAStringSet
#'
#' @examples
#' # Align and trim DNA sequences
#' DNA_seq <- system.file("extdata", "DNA_sequences.fas", package = "phyloR")
#' align_trim(DNA_seq, seq.type = "DNA", method = "ClustalW")
#'
#' # Align and trim protein sequences
#' protein_seq <- system.file("extdata", "protein_sequences.fas", package = "phyloR")
#' align_trim(protein_seq, seq.type = "protein", method = "ClustalW")
#'
#' @export
align_trim <- function(seq.file,
                       seq.type = "DNA",
                       method = "ClustalW",
                       output_file = NULL){


  # Read the sequence file based the sequence type
  if(toupper(seq.type) == "DNA"){
    mySeqs <- Biostrings::readDNAStringSet(seq.file)
  }
  if(tolower(seq.type) == "protein"){
    mySeqs <- Biostrings::readAAStringSet(seq.file)
  }

  # Perform multiple sequence alignment
  alignment <- msa::msa(mySeqs, method = method)

  # Convert alignment to bios2mds format and export it to a temporary FASTA file
  alignment_set <- msa::msaConvert(alignment, type = "bios2mds::align")
  bios2mds::export.fasta(alignment_set,
                         outfile = "alignment.fas",
                         ncol = 60,
                         open = "w")

  # Read the temporary file and trim the alignment
  aligned <- microseq::readFasta("alignment.fas")
  aln_trimmed <- microseq::msaTrim(aligned)

  # Remove the temporary FASTA file
  file.remove("alignment.fas")

  # Return processed sequences or save it to a file
  if(is.null(output_file)){
    return(aln_trimmed)
  }
  if(!is.null(output_file)){
    microseq::writeFasta(aln_trimmed, out.file = output_file)
  }
}

