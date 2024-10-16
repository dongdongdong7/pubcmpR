.load_idMappingTb <- function(){
  idMappingTb_path <- system.file("extdata", "idMappingTb.rds", package = "pubcmpR")
  message("Load idMappingTb...")
  if(file.exists(idMappingTb_path)) idMappingTb <- tibble::tibble(readRDS(idMappingTb_path)) %>% dplyr::select(CAS, DTXCID, CID, KEGG, ChEBI, HMDB, Drugbank, Name) %>% dplyr::distinct(CAS, DTXCID, CID, KEGG, ChEBI, HMDB, Drugbank, Name, .keep_all = TRUE)
  return(idMappingTb)
}
#' @title idMapping
#' @description
#' This function rely on metaboliteIDmapping data.
#' DTXCID: Comptox Chemical Dashboard
#' CID: Pubchem
#' CAS: CAS Registry numbers
#' HMDB: Human Metabolome Database
#' ChEBI: Chemical Entities of Biological Interest
#' KEGG: KEGG Compounds
#' Drugbank: Drugbank
#'
#' @param id id vector with database name.
#' @param to database name.
#' @param unique If TRUE, An id corresponds to only one conversion ID.
#'
#' @return A id vector.
#' @author Barry Song.
#' @export
#'
#' @examples
#' idMapping(id = "C18707")
#' id <- c(idMappingTb$CAS[1:5], idMappingTb$DTXSID[6:10], idMappingTb$DTXCID[4:11], idMappingTb$CID[12:15], idMappingTb$KEGG[12:15], idMappingTb$HMDB[16:19], idMappingTb$Drugbank[20:25])
#' idTb <- idMapping(id = id)
#' id <- c(idMappingTb$CAS[1:1000])
#' idTb <- idMapping(id = id, to = "HMDB", unique = FALSE)
idMapping <- function(id, to = c("HDMB", "CID", "CAS", "ChEBI", "KEGG", "Drugbank", "DTXCID")[1], unique = TRUE){
  idMappingTb_path <- system.file("extdata", "idMappingTb.rds", package = "pubcmpR")
  message("Load idMappingTb...")
  if(file.exists(idMappingTb_path) & !exists("idMappingTb")) idMappingTb <<- tibble::tibble(readRDS(idMappingTb_path)) %>% dplyr::select(CAS, DTXCID, CID, KEGG, ChEBI, HMDB, Drugbank, Name) %>% dplyr::distinct(CAS, DTXCID, CID, KEGG, ChEBI, HMDB, Drugbank, Name, .keep_all = TRUE)

  message("Mapping...")
  if(is.null(names(id))){
    # names(id) %in% c("DTXCID", "CID", "CAS", "HMDB", "ChEBI", "KEGG", "Drugbank")
    names(id) <- NA
    names(id)[stringr::str_detect(id, "-")] <- "CAS"
    names(id)[stringr::str_detect(id, "DTXCID")] <- "DTXCID"
    names(id)[stringr::str_detect(id, "^\\d+$")] <- "CID"
    names(id)[stringr::str_detect(id, "CHEBI:")] <- "ChEBI"
    names(id)[stringr::str_detect(id, "^C\\d+")] <- "KEGG"
    names(id)[stringr::str_detect(id, "HMDB")] <- "HMDB"
    names(id)[stringr::str_detect(id, "^(DB)\\d+")] <- "Drugbank"
  }

  idTb <- tibble::tibble(id = id, from = names(id))
  idTbList <- lapply(unique(idTb$from), function(x) {
    if(is.na(x)) tmp <- idTb %>% dplyr::filter(is.na(from))
    else tmp <- idTb %>% dplyr::filter(from == x)
    return(tmp)
  })

  idTbList <- lapply(1:length(idTbList), function(i) {
    idTbTmp <- idTbList[[i]]
    if(any(is.na(idTbTmp$from))){
      idTbTmp <- tibble::tibble(id = idTbTmp$id, to = NA)
      names(idTbTmp) <- c("id", to)
      return(idTbTmp)
    }
    else if(idTbTmp$from[1] == to){
      idTbTmp <- tibble::tibble(id = idTbTmp$id, to = idTbTmp$id)
      names(idTbTmp) <- c("id", to)
      return(idTbTmp)
    }
    else {
      idTbTmp <- dplyr::left_join(idTbTmp, idMappingTb, by = c("id" = idTbTmp$from[1])) %>%
        dplyr::select(id, dplyr::all_of(to)) %>%
        dplyr::distinct_all(.keep_all = TRUE)
      dup_index <- which(duplicated(purrr::pluck(idTbTmp[, "id"], 1), fromLast = FALSE) | duplicated(purrr::pluck(idTbTmp[, "id"], 1), fromLast = TRUE))
      delete_index <- dup_index[which(is.na(purrr::pluck(idTbTmp[dup_index, 2], 1)))]
      if(length(delete_index)!= 0) idTbTmp <- idTbTmp[-delete_index, ]
      if(unique){
        idTbTmp <- idTbTmp %>% dplyr::distinct(id, .keep_all = TRUE)
      }
      return(idTbTmp)
    }
  })
  idTb_new <- purrr::list_rbind(idTbList)
  return(idTb_new)
}
