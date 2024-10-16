# pubcmpR

The purpose of this R package is to create an R interface for easy access to small molecule information in public databases.

## Install

```R
devtools::install_github("dongdongdong7/pubcmpR")
```

## HMDB

### Load HMDB offline

```R
hmdbCmpTb <- load_hmdbCmpTb()
hmdbCmpTb
```

### Search HMDB online

```R
hmdbCmpTb_test <- search_hmdb_online(id = c("HMDB0000001", "HMDB0000002"))
hmdbCmpTb_test
```

## ChEBI

### Load ChEBI offline

```R
chebiCmpTb <- load_chebiCmpTb()
chebiCmpTb
```

## KEGG

### Load KEGG offline

```R
keggCmpTb <- load_keggCmpTb()
keggCmpTb
```

## Lipid Maps

### Load LMSD offline

```R
lipidmapsCmpTb <- load_lipidmapsCmpTb()
lipidmapsCmpTb
```

## PubChem

### Search PubChem online

```R
massSearch_CID(mass = 433.323)
search_pubchem_online(inputs = c("333", "4221"))
search_pubchem_online(inputs = c("Water"), type = "name")
search_pubchem_online(inputs = massSearch_CID(mass = 433.323, tol_mass = 0.0001), type = "cid")
```

## ID Mapping Function

```R
idMapping(id = "C18707", to = "CAS")
```

