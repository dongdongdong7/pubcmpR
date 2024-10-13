#' @title load_lipidmapsCmpTb
#'
#' @return A lipidmapsCmpTb.
#' @export
#'
#' @examples
#' lipidmapsCmpTb <- load_lipidmapsCmpTb()
load_lipidmapsCmpTb <- function(){
  lipidmapsCmpTb_path <- system.file("extdata", "lipidmapsCmpTb.rds", package = "pubcmpR")
  message("Load lipidmapsCmpTb...")
  if(file.exists(lipidmapsCmpTb_path)) lipidmapsCmpTb <- tibble::tibble(readRDS(lipidmapsCmpTb_path))
  else stop("Can not find lipidmapsCmpTb, please redownload!")
  return(lipidmapsCmpTb)
}
