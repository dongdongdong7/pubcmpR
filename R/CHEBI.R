#' @title load_chebiCmpTb
#' @description
#' CHEBI
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' chebiCmpTb <- load_chebiCmpTb()
load_chebiCmpTb <- function(){
  chebiCmpTb_path <- system.file("extdata", "chebiCmpTb.rds", package = "pubcmpR")
  message("Load chebiCmpTb...")
  if(file.exists(chebiCmpTb_path)) chebiCmpTb <- tibble::tibble(readRDS(chebiCmpTb_path))
  else stop("Can not find chebiCmpTb, please redownload!")
  return(chebiCmpTb)
}
