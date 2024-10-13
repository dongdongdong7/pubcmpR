# https://lipidmaps.org/databases/lmsd/download
# Last updated: 2024-10-13
library(magrittr)
sdf_file <- "D:/fudan/Projects/2024/pubcmpR/Progress/raw_database/LipidMaps/LMSD.sdf"
sdf_set <- ChemmineR::read.SDFset(sdf_file)
df_colnames <- unique(purrr::list_c(lapply(1:length(sdf_set), function(i) {
  names(sdf_set@SDF[[i]]@datablock)
})))

dfList <- lapply(1:length(sdf_set), function(i) {
  df <- as.data.frame(matrix(NA, ncol = 21))
  colnames(df) <- df_colnames
  df$LM_ID <- sdf_set@SDF[[i]]@datablock["LM_ID"]
  df$NAME <- sdf_set@SDF[[i]]@datablock["NAME"]
  df$SYSTEMATIC_NAME <- sdf_set@SDF[[i]]@datablock["SYSTEMATIC_NAME"]
  df$CATEGORY <- sdf_set@SDF[[i]]@datablock["CATEGORY"]
  df$MAIN_CLASS <- sdf_set@SDF[[i]]@datablock["MAIN_CLASS"]
  df$EXACT_MASS <- sdf_set@SDF[[i]]@datablock["EXACT_MASS"]
  df$FORMULA <- sdf_set@SDF[[i]]@datablock["FORMULA"]
  df$INCHI_KEY <- sdf_set@SDF[[i]]@datablock["INCHI_KEY"]
  df$INCHI <- sdf_set@SDF[[i]]@datablock["INCHI"]
  df$SMILES <- sdf_set@SDF[[i]]@datablock["SMILES"]
  df$ABBREVIATION <- sdf_set@SDF[[i]]@datablock["ABBREVIATION"]
  df$SYNONYMS <- sdf_set@SDF[[i]]@datablock["SYNONYMS"]
  df$PUBCHEM_CID <- sdf_set@SDF[[i]]@datablock["PUBCHEM_CID"]
  df$CHEBI_ID <- sdf_set@SDF[[i]]@datablock["CHEBI_ID"]
  df$KEGG_ID <- sdf_set@SDF[[i]]@datablock["KEGG_ID"]
  df$HMDB_ID <- sdf_set@SDF[[i]]@datablock["HMDB_ID"]
  df$SWISSLIPIDS_ID <- sdf_set@SDF[[i]]@datablock["SWISSLIPIDS_ID"]
  df$SUB_CLASS <- sdf_set@SDF[[i]]@datablock["SUB_CLASS"]
  df$LIPIDBANK_ID <- sdf_set@SDF[[i]]@datablock["LIPIDBANK_ID"]
  df$PLANTFA_ID <- sdf_set@SDF[[i]]@datablock["PLANTFA_ID"]
  df$CLASS_LEVEL4 <- sdf_set@SDF[[i]]@datablock["CLASS_LEVEL4"]
  return(df)
})
df <- purrr::list_rbind(dfList)

lipidmapsCmpTb <- df %>%
  dplyr::select(LM_ID, NAME, SYSTEMATIC_NAME, SYNONYMS, ABBREVIATION, FORMULA, EXACT_MASS,
                SMILES, INCHI, INCHI_KEY,
                CATEGORY, MAIN_CLASS, SUB_CLASS, CLASS_LEVEL4,
                PUBCHEM_CID, CHEBI_ID, KEGG_ID, HMDB_ID, SWISSLIPIDS_ID, LIPIDBANK_ID, PLANTFA_ID) %>%
  tibble::tibble()
lipidmapsCmpTb$EXACT_MASS <- as.numeric(lipidmapsCmpTb$EXACT_MASS)
saveRDS(lipidmapsCmpTb, file = "./inst/extdata/lipidmapsCmpTb.rds")
