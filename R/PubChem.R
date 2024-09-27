# cid2sdf_online(cids = c(333, 3040203, 40030312), existTest = FALSE)
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
