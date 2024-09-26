#BiocManager::install("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub()
datasets <- query( ah, "metaboliteIDmapping")
data1 <- ah[["AH79817"]]
data2 <- ah[["AH83115"]]
data3 <- ah[["AH91792"]]
saveRDS(data2, file = "./idMappingTb.rds")
