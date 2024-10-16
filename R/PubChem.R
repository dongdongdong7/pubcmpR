# .cid2sdf_online(cids = c(333, 3040203, 40030312), existTest = FALSE)
.cid2sdf_online <- function(cids, existTest = FALSE){
  if(existTest){
    pubchemServerURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    urls <- sapply(cids, function(x) {
      paste(pubchemServerURL, "compound", "cid", x, "SDF", sep = "/")
    })
    existsRes <- RCurl::url.exists(urls)
    if(!all(existsRes)){
      stop(paste(names(existsRes[!existsRes]), "not found!"))
    }
  }
  ChemmineR::pubchemCidToSDF(cids = cids)
}

.sdf2smiles <- function(SDFset){
  SMIset <- ChemmineR::sdf2smiles(SDFset)
  smiles <- SMIset@smilist
  nameVec <- sub("\r", "", names(smiles))
  smiles <- as.character(smiles)
  names(smiles) <- nameVec
  return(smiles)
}

# test <- .search_pubchem("CCCCC", "smiles")
.search_pubchem <- function(input, type = c("cid",
                                            "name",
                                            "smiles",
                                            "inchi",
                                            "inchikey")[1], synonyms = FALSE){
  pubchem_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
  domain <- "compound"
  url <- paste(pubchem_url, domain, sep = "/")
  if(!is.atomic(input)) stop("Input length should be one")

  if(type == "cid"){
    url <- paste(url, "cid", sep = "/")
    sdf_url <- paste(url, input, "SDF", sep = "/")
    synonyms_url <- paste(url, input, "synonyms/TXT", sep = "/")
  }else if(type == "name"){
    url <- paste(url, "name", sep = "/")
    sdf_url <- paste(url, input, "SDF", sep = "/")
    synonyms_url <- paste(url, input, "synonyms/TXT", sep = "/")
  }else if(type == "smiles"){
    cid <- readLines(paste(url, "smiles", input, "cids/TXT", sep = "/"))
    sdf_url <- paste(url, "cid", cid, "SDF", sep = "/")
    synonyms_url <- paste(url, "cid", cid, "synonyms/TXT", sep = "/")
  }else if(type == "inchi"){
    inchis <- input
    inchi_escaped <- stringi::stri_replace_all_fixed(str = inchis,
                                                     pattern = c("/", "[", "]", "&", ",", "(", ")", "=",
                                                                 "+"), replacement = c("%2F", "%5B", "%5D", "%26",
                                                                                       "%2C", "%28", "%29", "%3D", "%2B"), vectorize_all = FALSE)
    inchi_request_base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/TXT?inchi="
    req_url <- paste0(inchi_request_base_url, inchi_escaped)
    cid <- readLines(req_url)
    sdf_url <- paste(url, "cid", cid, "SDF", sep = "/")
    synonyms_url <- paste(url, "cid", cid, "synonyms/TXT", sep = "/")
  }else if(type == "inchikey"){
    cid <- readLines(paste(url, "inchikey", input, "cids/TXT", sep = "/"))
    sdf_url <- paste(url, "cid", cid, "SDF", sep = "/")
    synonyms_url <- paste(url, "cid", cid, "synonyms/TXT", sep = "/")
  }

  sdf <- ChemmineR::read.SDFset(readLines(sdf_url))[[1]]
  if(!synonyms){
    df_colnames <- c("cid", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey")
  }else{
    synonyms_vec <- readLines(synonyms_url)
    df_colnames <- c("cid", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey", "synonyms")
  }

  df <- as.data.frame(matrix(NA,ncol = length(df_colnames)))
  colnames(df) <- df_colnames

  df$cid <- unname(sdf@datablock["PUBCHEM_COMPOUND_CID"])
  df$name <- unname(sdf@datablock["PUBCHEM_IUPAC_TRADITIONAL_NAME"])
  df$formula <- unname(sdf@datablock["PUBCHEM_MOLECULAR_FORMULA"])
  df$exact_mass <- as.numeric(unname(sdf@datablock["PUBCHEM_EXACT_MASS"]))
  df$smiles <- unname(sdf@datablock["PUBCHEM_OPENEYE_CAN_SMILES"])
  df$inchi <- unname(sdf@datablock["PUBCHEM_IUPAC_INCHI"])
  df$inchikey <- unname(sdf@datablock["PUBCHEM_IUPAC_INCHIKEY"])

  if(synonyms) df$synonyms <- list(synonyms_vec)

  return(df)
}

#' @title massSearch_CID
#' @description
#' Searching compounds with exact mass.
#'
#' @param mass A numeric vector of length 1 or 2.
#' @param tol_mass Use tol_mass to compute a mass range when the mass vector length is 1.
#'
#' @return A character vector with cids.
#' @author Barry Song
#' @export
#'
#' @examples
#' massSearch_CID(mass = 433.323)
massSearch_CID <- function(mass, tol_mass = 0.001){
  url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/exact_mass/range/"
  if(is.atomic(mass) & is.numeric(mass)){
    mass_range <- c(mass - tol_mass, mass + tol_mass)
    cids_url <- paste(url, mass_range[1], mass_range[2], "cids", "TXT", sep = "/")
  }else if(!is.atomic(mass) & is.numeric(mass) & length(mass) == 2){
    mass_range <- mass
    cids_url <- paste(url, mass_range[1], mass_range[2], "cids", "TXT", sep = "/")
  }else stop("Input is wrong!")
  cids <- readLines(cids_url)
  return(cids)
}

#' @title search_pubchem_online
#' @description
#' Search PubChem database online.
#'
#' @param inputs A character vector.
#' @param type Supports five search types, cid, name, smiles, inchi and inchikey.
#' @param synonyms Whether to return synonyms in the result.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' search_pubchem_online(inputs = c("333", "4221"))
#' search_pubchem_online(inputs = c("Water"), type = "name")
#' search_pubchem_online(inputs = massSearch_CID(mass = 433.323, tol_mass = 0.0001), type = "cid")
search_pubchem_online <- function(inputs, type = c("cid", "name", "smiles", "inchi", "inchikey")[1], synonyms = FALSE){
  num <- length(inputs)
  message(sprintf("You're searching for %d small molecules in PubChem.", num))
  pb <- txtProgressBar(max = num, style = 3)
  dfList <- lapply(1:num, function(i) {
    if(i %% 5 == 0) Sys.sleep(1)
    input <- inputs[i]
    setTxtProgressBar(pb, i)
    tryCatch(.search_pubchem(input = input, type = type, synonyms = synonyms),
             error = function(e) {
               message(e)
               NULL
             })
  })
  df <- purrr::list_rbind(dfList) %>% tibble::tibble()
  message("\n")
  return(df)
}
