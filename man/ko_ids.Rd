\name{ko_ids}
\alias{ko_ids}
\title{KO IDs List}
\description{
  A dataset containing a list of all KO (KEGG Orthology) IDs. This dataset is used for biological pathway analysis, gene function annotation, and other related bioinformatics tasks such as pathway enrichment analysis and functional genomics studies.
  KEGG Orthology (KO) IDs are used to link genes and their associated pathways in the KEGG database.
}
\usage{
  ko_ids
}
\format{
  A character vector containing the KO IDs.
  Each element in the vector corresponds to a unique KO identifier,
  which typically consists of a 'K' prefix followed by a numeric code (e.g., K00001).
  The data does not contain any duplicate entries, ensuring each KO ID is unique.
}
\examples{
  # Load the ko_ids dataset
  data(ko_ids)

  # View the first few KO IDs
  head(ko_ids)

  # Check the length of the dataset
  length(ko_ids)
}
\source{
  Data derived from the KEGG database (Kyoto Encyclopedia of Genes and Genomes).
}


