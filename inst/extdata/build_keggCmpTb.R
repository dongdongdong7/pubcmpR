#BiocManager::install("KEGGREST")
#BiocManager::install("ChemmineOB")
#BiocManager::install("ChemmineR")
# KEGG COMPOUND Database
# 241014

KEGGREST::listDatabases()
cmp_id <- KEGGREST::keggList("compound") # 19406

#mols <- KEGGREST::keggGet(names(cmp_id)[13058], option = "mol")

dfList <- lapply(1:length(cmp_id), function(i) {
  print(paste0(i, " / ", length(cmp_id)))
  Sys.sleep(runif(n = 1, min = 0, max = 1))
  id <- names(cmp_id[i])
  name <- unname(cmp_id[i])
  molstr <- tryCatch(KEGGREST::keggGet(id, option = "mol"),
                     error = function(e) {
                       print(e)
                       NA
                     })
  if(is.na(molstr)){
    df <- as.data.frame(matrix(c(id, name, NA, NA, NA, NA, NA), ncol = 7))
    colnames(df) <- c("id", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey")
    return(df)
  }
  else{
    mol <- c("", "", "", stringr::str_split(molstr, "\n")[[1]]) # Title line (can be blank but line must exist) (3 lines)
    sdf <- tryCatch(ChemmineR::read.SDFset(ChemmineR::read.SDFstr(mol)),
                    warning = function(w) {
                      print(w)
                      sdf <- ChemmineR::read.SDFset(ChemmineR::read.SDFstr(mol))
                      valid <- ChemmineR::validSDF(sdf)
                      sdf[valid]
                    },
                    error = function(e) {
                      print(e)
                      NA
                    })
    if(is.na(sdf)) {
      df <- as.data.frame(matrix(c(id, name, NA, NA, NA, NA, NA), ncol = 7))
      colnames(df) <- c("id", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey")
      return(df)
    }
    if(length(sdf) == 0){
      df <- as.data.frame(matrix(c(id, name, NA, NA, NA, NA, NA), ncol = 7))
      colnames(df) <- c("id", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey")
      return(df)
    }
    formula <- tryCatch(unname(ChemmineR::MF(sdf)),
                        error = function(e) {
                          print(e)
                          NA
                        })
    exact_mass <- tryCatch(unname(ChemmineR::MW(sdf)),
                           error = function(e) {
                             print(e)
                             NA
                           })
    smiles <- tryCatch(as.character(ChemmineR::sdf2smiles(sdf)[[1]]),
                       error = function(e) {
                         print(e)
                         ""
                       })
    if(smiles == "") smiles <- NA
    inchi <- tryCatch(rinchi::get.inchi(smiles),
                      error = function(e) {
                        print(e)
                        NA
                      })
    inchikey <- tryCatch(rinchi::get.inchi.key(smiles),
                         error = function(e) {
                           print(e)
                           NA
                         })
    df <- as.data.frame(matrix(c(id, name, formula, exact_mass, smiles, inchi, inchikey), ncol = 7))
    colnames(df) <- c("id", "name", "formula", "exact_mass", "smiles", "inchi", "inchikey")
    return(df)
  }
})
df <- purrr::list_rbind(dfList)
name_vec <- unname(sapply(df$name, function(x) {
  stringr::str_split(x, ";")[[1]][1]
}))
synonym_list <- lapply(df$name, function(x) {
  synonym_tmp <- stringr::str_split(x, ";")[[1]][-1]
  if(length(synonym_tmp) == 0) return(NA)
  else return(synonym_tmp)
})
keggCmpTb <- tibble::tibble(df)
keggCmpTb$name <- name_vec
keggCmpTb$synonym <- synonym_list
keggCmpTb <- keggCmpTb[which(!(is.na(keggCmpTb$formula) & is.na(keggCmpTb$exact_mass) & is.na(keggCmpTb$smiles))), ]
saveRDS(keggCmpTb, file = "./inst/extdata/keggCmpTb.rds")
