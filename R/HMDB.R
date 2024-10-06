# Enter a HMDB id and return a .xml url.
.hmxPath <- function(id = "HMDB0000001")
{
  prefix = "http://www.hmdb.ca/metabolites/"
  sub("__PRE__", prefix, sub("%%ID%%", id, "__PRE__%%ID%%.xml"))
}
# Enter a HMDB id and return a page list.
.hmxToList <- function(id = "HMDB0000001"){
  stopifnot(is.atomic(id),
            length(id) == 1)
  txt = readLines(.hmxPath(id = id))
  prs = XML::xmlTreeParse(txt, asText = TRUE)
  XML::xmlToList(prs)
}
#' @title load_hmdbCmpTb
#'
#' @return hmdbCmpTb
#' @export
#'
#' @examples
#' hmdbCmpTb <- load_hmdbCmpTb()
load_hmdbCmpTb <- function(){
  hmdbCmpTb_path <- system.file("extdata", "hmdbCmpTb.rds", package = "pubcmpR")
  message("Load hmdbCmpTb...")
  if(file.exists(hmdbCmpTb_path)) hmdbCmpTb <- tibble::tibble(readRDS(hmdbCmpTb_path))
  else stop("Can not find hmdbCmpTb, please redownload!")
  return(hmdbCmpTb)
}
