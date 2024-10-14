# https://www.ebi.ac.uk/chebi/downloadsForward.do
# ChEBI
# 241013
library(magrittr)
sdf_file <- "D:/fudan/Projects/2024/pubcmpR/Progress/raw_database/ChEBI/ChEBI_complete_3star.sdf"
sdf_set <- ChemmineR::read.SDFset(sdf_file, skipErrors = TRUE)
df_colnames <- unique(purrr::list_c(lapply(1:length(sdf_set), function(i) {
  names(sdf_set@SDF[[i]]@datablock)
})))
df_colnames <- c("ChEBI ID", "ChEBI Name", "Secondary ChEBI ID",
                 "InChI", "InChIKey", "SMILES", "Formulae", "Monoisotopic Mass",
                 "IUPAC Names", "Synonyms",
                 "CAS Registry Numbers", "KEGG COMPOUND Database Links", "LIPID MAPS instance Database Links",
                 "PubChem Database Links", "Chemspider Database Links",
                 "FooDB Database Links", "HMDB Database Links", "Wikipedia Database Links",
                 "SwissLipids Database Links", "DrugBank Database Links", "YMDB Database Links")
dfList <- lapply(1:length(sdf_set), function(i) {
  df <- as.data.frame(matrix(sdf_set@SDF[[i]]@datablock[df_colnames], ncol = length(df_colnames)))
  colnames(df) <- df_colnames
  return(df)
})
df <- purrr::list_rbind(dfList)
tb <- purrr::list_rbind(dfList) %>%
  tibble::tibble()

for(i in c(3, 9:length(df_colnames))){
  print(i)
  colname <- df_colnames[i]
  if(any(stringr::str_detect(df[, colname][!is.na(df[, colname])], " __ "))){
    tb[[i]] <- lapply(1:nrow(df), function(j) {
      stringr::str_split(df[j, i], " __ ")[[1]]
    })
  }
}
tb$`PubChem Database Links` <- lapply(1:nrow(tb), function(i) {
  tmp <- tb$`PubChem Database Links`[[i]][stringr::str_detect(tb$`PubChem Database Links`[[i]], "CID")]
  if(length(tmp) == 0) return(NA)
  else return(stringr::str_replace(tmp, "CID:\\s+", ""))
})

chebiCmpTb <- tb %>%
  dplyr::filter(!stringr::str_detect(Formulae, "\\(*+\\)[n|x|y]")) %>%
  dplyr::filter(!stringr::str_detect(Formulae, "[R|X]")) %>%
  dplyr::filter(!stringr::str_detect(Formulae, "\\[*+\\][n|x|y]")) %>%
  dplyr::filter(!stringr::str_detect(Formulae, "\\(*+\\)\\d")) %>%
  dplyr::filter(!stringr::str_detect(Formulae, "\\[*+\\]\\d")) %>%
  dplyr::filter(!stringr::str_detect(Formulae, "\\(*+\\)[ran|.]"))

mass_vec <- as.vector(MetaboCoreUtils::mz2mass(MetaboCoreUtils::formula2mz(chebiCmpTb[which(is.na(chebiCmpTb$`Monoisotopic Mass`)), ]$Formulae)))
chebiCmpTb$`Monoisotopic Mass`[which(is.na(chebiCmpTb$`Monoisotopic Mass`))] <- mass_vec
chebiCmpTb <- chebiCmpTb %>%
  dplyr::filter(!is.na(`Monoisotopic Mass`))
chebiCmpTb$`Monoisotopic Mass` <- as.numeric(chebiCmpTb$`Monoisotopic Mass`)

saveRDS(chebiCmpTb, file = "./inst/extdata/chebiCmpTb.rds")
