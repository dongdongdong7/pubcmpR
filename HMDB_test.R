BiocManager::install("hmdbQuery")
library(hmdbQuery)
lk1 = HmdbEntry(prefix = "http://www.hmdb.ca/metabolites/",
                id = "HMDB0000001")
diseases(lk1)[1,]

prefix = "http://www.hmdb.ca/metabolites/";id = "HMDB0000001";keepFull = TRUE
imp = .hmxToList(prefix = prefix, id = id)
nimp = names(imp)

# metaboliteIDmapping
# biodbHmdb
# hmdbQuery
