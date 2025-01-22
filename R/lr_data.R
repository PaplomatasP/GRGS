#' Ligand-Receptor Pairs Dataset
#'
#' A dataset containing ligand-receptor pairs for human genes.
#' This dataset is used in the genetic algorithm to prioritize features related to cell-cell communication.
#'
#' @format A data frame with 2558 rows and 2 variables:
#' \describe{
#'   \item{ligand_gene_symbol}{Symbol of the ligand gene (e.g., "TGFB1")}
#'   \item{receptor_gene_symbol}{Symbol of the receptor gene (e.g., "TGFBR1")}
#' }
#' @references
#' Αναφορά στη δημοσίευση που περιγράφει τα δεδομένα (π.χ., CellChat, CellPhoneDB, κλπ).
#' @examples
#' data(lr_data)
#' head(lr_data)
"lr_data"
