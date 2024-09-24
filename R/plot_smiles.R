#' @title plot_smiles
#' @description
#' Input a smiles vector return a molecule pic.
#'
#' @param smiles A named vector.
#'
#' @return pic.
#' @export
#'
#' @examples
#' plot_smiles(smiles = "C-C-C-C")
plot_smiles <- function(smiles){
  if(is.null(names(smiles))) names(smiles) <- smiles
  sdf <- ChemmineR::smiles2sdf(smiles)
  ChemmineR::plot(sdf)
}
