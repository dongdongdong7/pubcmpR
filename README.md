# pubcmpR

The purpose of this R package is to create an R interface for easy access to small molecule information in public databases.

## HMDB

### Loading HMDB offline

```R
hmdbCmpTb <- load_hmdbCmpTb()
hmdbCmpTb
```

```R
# A tibble: 217,920 × 35
   accession   secondary_accession name                  synonyms  chemical_formula average_molecular_we…¹
   <chr>       <list>              <chr>                 <list>    <chr>                             <dbl>
 1 HMDB0000001 <chr [7]>           1-Methylhistidine     <chr>     C7H11N3O2                         169. 
 2 HMDB0000002 <chr [3]>           1,3-Diaminopropane    <chr>     C3H10N2                            74.1
 3 HMDB0000005 <chr [3]>           2-Ketobutyric acid    <chr>     C4H6O3                            102. 
 4 HMDB0000008 <chr [1]>           2-Hydroxybutyric acid <chr>     C4H8O3                            104. 
 5 HMDB0000010 <chr [5]>           2-Methoxyestrone      <chr [7]> C19H24O3                          300. 
 6 HMDB0000011 <chr [3]>           3-Hydroxybutyric acid <chr>     C4H8O3                            104. 
 7 HMDB0000012 <chr [1]>           Deoxyuridine          <chr>     C9H12N2O5                         228. 
 8 HMDB0000014 <chr [1]>           Deoxycytidine         <chr>     C9H13N3O4                         227. 
 9 HMDB0000015 <chr [1]>           Cortexolone           <chr>     C21H30O4                          346. 
10 HMDB0000016 <chr [1]>           Deoxycorticosterone   <chr>     C21H30O3                          330. 
# ℹ 217,910 more rows
# ℹ abbreviated name: ¹​average_molecular_weight
# ℹ 29 more variables: monisotopic_molecular_weight <dbl>, iupac_name <chr>, traditional_iupac <chr>,
#   cas_registry_number <chr>, smiles <chr>, inchi <chr>, inchikey <chr>, taxonomy <list>,
#   ontology <list>, cellular_locations <list>, biospecimen_locations <list>, tissue_locations <list>,
#   pathway <list>, diseases <list>, foodb_id <chr>, kegg_id <chr>, chemspider_id <chr>, chebi_id <chr>,
#   pubchem_compound_id <chr>, pdb_id <chr>, biocyc_id <chr>, drugbank_id <chr>, …
# ℹ Use `print(n = ...)` to see more rows
```

### Search HMDB online

```R
hmdbCmpTb_test <- search_hmdb_online(id = c("HMDB0000001", "HMDB0000002"))
hmdbCmpTb_test
```

```R
# A tibble: 2 × 35
  accession   secondary_accession name               synonyms   chemical_formula average_molecular_weight
  <chr>       <list>              <chr>              <list>     <chr>                               <dbl>
1 HMDB0000001 <chr [7]>           1-Methylhistidine  <chr [16]> C7H11N3O2                           169. 
2 HMDB0000002 <chr [3]>           1,3-Diaminopropane <chr [13]> C3H10N2                              74.1
# ℹ 29 more variables: monisotopic_molecular_weight <dbl>, iupac_name <chr>, traditional_iupac <chr>,
#   cas_registry_number <chr>, smiles <chr>, inchi <chr>, inchikey <chr>, taxonomy <list>,
#   ontology <list>, cellular_locations <list>, biospecimen_locations <list>, tissue_locations <list>,
#   pathway <list>, diseases <list>, foodb_id <chr>, kegg_id <chr>, chemspider_id <chr>, chebi_id <chr>,
#   pubchem_compound_id <chr>, pdb_id <lgl>, biocyc_id <chr>, drugbank_id <chr>,
#   phenol_explorer_compound_id <lgl>, wikipedia_id <chr>, knapsack_id <chr>, bigg_id <lgl>,
#   metlin_id <chr>, vmh_id <lgl>, protein_associations <list>
```

