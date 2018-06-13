#' Pairs of kegg network.
#'
#' A dataset containing kegg IDs in the KEGG database with all networks.
#'
#' @format A data frame with 57070 rows and 2 variables:
#' \describe{
#'   \item{r1}{KEGG IDs}
#'   \item{r2}{KEGG IDs, which have a connection with KEGG ID in the first column}
#'   ...
#' }
"kegg_network"


#' Inchikey database.
#'
#' A dataset containing PubChem CIDs, InchiKey in the PubChem database.
#'
#' @format A data frame with 101494 rows and 2 variables:
#' \describe{
#'   \item{CID}{PubChem CIDs}
#'   \item{InchiKey}{Inchikeys}
#'   ...
#' }
"InchiKey"


#' Example of input dataset, in which colnames does not meet requirement.
#'
#' A dataset which can be used as input dataset and its row names do not match the default row names.
#'
#' @format A data frame with 20 rows and 6 variables:
#' \describe{
#'   \item{Query.Mass}{Mass of compounds.}
#'   \item{Name}{Names of putative IDs.}
#'   \item{Formula}{Formulas of putative IDs.}
#'   \item{Exact.Mass}{Exact mass of putative IDs.}
#'   \item{PubChem.CID}{PubChem IDs of putative IDs.}
#'   \item{KEGG.ID}{KEGG IDs of putative IDs.}
#'   ...
#' }
"demo1"


#' Example of input dataset, in which colnames does not meet requirement.
#'
#' A dataset which can be used as input dataset and its row names do not match the default row names.
#'
#'
#' @format A data frame with 3592 rows and 6 variables:
#' \describe{
#'   \item{Query.Mass}{Mass of compounds.}
#'   \item{Name}{Names of putative IDs.}
#'   \item{Formula}{Formulas of putative IDs.}
#'   \item{Exact.Mass}{Exact mass of putative IDs.}
#'   \item{PubChem.CID}{PubChem IDs of putative IDs.}
#'   \item{KEGG.ID}{KEGG IDs of putative IDs.}
#'   ...
#' }
"demo2"



#' MetID: A package for Network-based prioritization of putative metabolite IDs.
#'
#' The foo package provides one important functions:
#' get_scores_for_LC_MS
#'
#' @section Foo functions:
#' get_scores_for_LC_MS: Get scores for metabolite putative IDs by LC-MS.
#'
#' @docType package
#' @name MetID
NULL
