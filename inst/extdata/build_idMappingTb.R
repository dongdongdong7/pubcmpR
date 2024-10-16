#BiocManager::install("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
datasets <- query( ah, "metaboliteIDmapping")
data1 <- ah[["AH79817"]]
data2 <- ah[["AH83115"]]
data3 <- ah[["AH91792"]]

data2$ChEBI <- unname(sapply(data2$ChEBI, function(x) {
  if(is.na(x)) return(NA)
  else paste0("CHEBI:", x)
}))

saveRDS(data2, file = "./inst/extdata/idMappingTb.rds")

