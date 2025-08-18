\name{species_tbl}
\alias{species_tbl}
\title{Species Information Table}
\description{
  A data frame containing detailed information about various species, including their KEGG organization ID (KEGG_org_ID), NCBI organization ID (NCBI_org_ID), abbreviation, scientific name, and classification.
  This dataset is useful for species comparison, genomics studies, and bioinformatics workflows, where identification and classification of species are required.
}
\usage{
  species_tbl
}
\format{
  A data frame with the following columns:
  \describe{
    \item{KEGG_org_ID}{A character vector containing the KEGG organization IDs for each species (e.g., 'T01001' for Homo sapiens).}
    \item{NCBI_org_ID}{A character vector containing the NCBI organization IDs for each species (e.g., '9606' for Homo sapiens).}
    \item{abbreviation}{A character vector with the abbreviations used for each species (e.g., 'hsa' for Homo sapiens).}
    \item{scientific_name}{A character vector with the full scientific names of each species (e.g., 'Homo sapiens').}
  }
}
\examples{
  # Load the species_tbl dataset
  data(species_tbl)

  # View the first few rows of the species table
  head(species_tbl)

  # Check the structure of the dataset
  str(species_tbl)

  # Access a specific column, such as the scientific names
  species_tbl$scientific_name

  # Filter species by KEGG ID (e.g., Homo sapiens)
  species_tbl[species_tbl$KEGG_org_ID == "hsa", ]
}
\source{
  Data derived from KEGG and NCBI databases.
}


