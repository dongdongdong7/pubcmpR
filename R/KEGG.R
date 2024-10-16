#' @title load_keggCmpTb
#' @description
#' Load the KEGG Compounds Database.
#'
#' @return keggCmpTb.
#' @author Barry Song.
#' @export
#'
#' @examples
#' keggCmpTb <- load_keggCmpTb()
load_keggCmpTb <- function(){
  keggCmpTb_path <- system.file("extdata", "keggCmpTb.rds", package = "pubcmpR")
  message("Load keggCmpTb...")
  if(file.exists(keggCmpTb_path)) keggCmpTb <- tibble::tibble(readRDS(keggCmpTb_path))
  else stop("Can not find keggCmpTb, please redownload!")
  return(keggCmpTb)
}
