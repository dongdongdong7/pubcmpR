#' @title load_hmdbCmpTb
#'
#' @return hmdbCmpTb
#' @export
#'
#' @examples
#' hmdbCmpTb <- load_hmdbCmpTb()
load_hmdbCmpTb <- function(){
  hmdbCmpTb_path <- system.file("extdata", "hmdbCmpTb.rds", package = "pubcmpR")
  message("Load hmdbCmpTb...")
  if(file.exists(hmdbCmpTb_path)) hmdbCmpTb <- tibble::tibble(readRDS(hmdbCmpTb_path))
  else stop("Can not find hmdbCmpTb, please redownload!")
  return(hmdbCmpTb)
}

# Enter a HMDB id and return a .xml url.
.hmxPath <- function(id = "HMDB0000001")
{
  prefix = "http://www.hmdb.ca/metabolites/"
  sub("__PRE__", prefix, sub("%%ID%%", id, "__PRE__%%ID%%.xml"))
}
# Enter a HMDB id and return a page list.
.hmxToList <- function(id = "HMDB0000001"){
  stopifnot(is.atomic(id),
            length(id) == 1)
  doc <- xml2::read_xml(.hmxPath(id = id))
  xml2::as_list(doc)$metabolite
}
# Ontology to list
.process_ontology_root <- function(root){
  extract_elements <- function(descendant) {
    result <- list(term = NULL, level = NULL)  # 初始化结果列表

    # 检查是否为列表
    if (is.list(descendant)) {
      # 提取 a, b, c
      if ("term" %in% names(descendant)) {
        result$term <- c(result$term, descendant$term)
      }
      if ("level" %in% names(descendant)) {
        result$level <- c(result$level, descendant$level)
      }

      # 递归提取 descants
      if ("descendants" %in% names(descendant)) {
        descendatnsNum <- length(descendant$descendants)
        for (i in 1:descendatnsNum) {
          descant_item <- descendant$descendants[[i]]
          nested_result <- extract_elements(descant_item)
          result$term <- c(result$term, nested_result$term)
          result$level <- c(result$level, nested_result$level)
        }
      }
    }

    return(result)
  }
  if(!is.list(root)) return(NULL)
  tmp <- extract_elements(root)

  df <- tibble::tibble(term = as.character(tmp$term), level = as.integer(tmp$level))

  return(df)
}

#' @title search_hmdb_online
#' @description
#' Search HMDB database online.
#'
#' @param id A id vector.
#'
#' @return A hmdbCmpTb tibble.
#' @export
#'
#' @examples
#' hmdbCmpTb_test <- search_hmdb_online(id = c("HMDB0000001", "HMDB0000002"))
search_hmdb_online <- function(id){
  metabolitesList <- lapply(id, function(x) {
    Sys.sleep(runif(n = 1, min = 0, max = 2))
    message(x)
    .hmxToList(x)
  })
  num <- length(metabolitesList)
  accession_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    metabolite$accession[[1]]
  })
  secondary_accession_list <- lapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    unname(sapply(metabolite$secondary_accessions, function(x) {
      if(stringr::str_detect(x[[1]], "\n")) return(NA)
      else return(x[[1]])
    }))
  })
  name_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    metabolite$name[[1]]
  })
  synonyms_list <- lapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    unname(sapply(metabolite$synonyms, function(x) {
      if(stringr::str_detect(x[[1]], "\n")) return(NA)
      else return(x[[1]])
    }))
  })
  chemical_formula_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$chemical_formula) == 1, metabolite$chemical_formula[[1]], NA)
  })
  average_molecular_weight_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    as.numeric(ifelse(length(metabolite$average_molecular_weight) == 1, metabolite$average_molecular_weight[[1]], NA))
  })
  monisotopic_molecular_weight_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    as.numeric(ifelse(length(metabolite$monisotopic_molecular_weight) == 1, metabolite$monisotopic_molecular_weight[[1]], NA))
  })
  iupac_name_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$iupac_name) == 1, metabolite$iupac_name[[1]], NA)
  })
  traditional_iupac_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$traditional_iupac) == 1, metabolite$traditional_iupac[[1]], NA)
  })
  cas_registry_number_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$cas_registry_number) == 1, metabolite$cas_registry_number[[1]], NA)
  })
  smiles_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$smiles) == 1, metabolite$smiles[[1]], NA)
  })
  inchi_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$inchi) == 1, metabolite$inchi[[1]], NA)
  })
  inchikey_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$inchikey) == 1, metabolite$inchikey[[1]], NA)
  })
  taxonomy_list <- lapply(1:num, function(i) {
    #print(paste0(i, " / ", num))
    metabolite <- metabolitesList[[i]]
    if(is.null(metabolite$taxonomy$alternative_parents)){
      alternative_parents <- NA
    }else{
      alternative_parents <- unname(sapply(metabolite$taxonomy$alternative_parents, function(x) {
        if(stringr::str_detect(x[[1]], "\n")) return(NA)
        else return(x[[1]])
      }))
    }
    if(is.null(metabolite$taxonomy$substituents)){
      substituents <- NA
    }else{
      substituents <- unname(sapply(metabolite$taxonomy$substituents, function(x) {
        if(stringr::str_detect(x[[1]], "\n")) return(NA)
        else return(x[[1]])
      }))
    }
    if(is.null(metabolite$taxonomy$external_descriptors)){
      external_descriptors <- NA
    }else{
      external_descriptors <- unname(sapply(metabolite$taxonomy$external_descriptors, function(x) {
        if(stringr::str_detect(x[[1]], "\n")) return(NA)
        else return(x[[1]])
      }))
    }
    list(kingdom = ifelse(length(metabolite$taxonomy$kingdom) == 1, metabolite$taxonomy$kingdom[[1]], NA),
         super_class = ifelse(length(metabolite$taxonomy$super_class) == 1, metabolite$taxonomy$super_class[[1]], NA),
         class = ifelse(length(metabolite$taxonomy$class) == 1, metabolite$taxonomy$class[[1]], NA),
         sub_class = ifelse(length(metabolite$taxonomy$sub_class) == 1, metabolite$taxonomy$sub_class[[1]], NA),
         direct_parent = ifelse(length(metabolite$taxonomy$direct_parent) == 1, metabolite$taxonomy$direct_parent[[1]], NA),
         alternative_parents = alternative_parents,
         substituents = substituents,
         molecular_framework = ifelse(length(metabolite$taxonomy$molecular_framework) == 1, metabolite$taxonomy$molecular_framework[[1]], NA),
         external_descriptors = external_descriptors)
  })
  ontology_list <- lapply(1:num, function(i) {
    #print(paste0(i, " / ", num))
    metabolite <- metabolitesList[[i]]
    ontologyList <- lapply(metabolite$ontology, .process_ontology_root)
    ontology <- list(`Physiological effect` = NULL, Disposition = NULL, Process = NULL, Role = NULL)
    for(j in 1:length(ontologyList)){
      #print(j)
      if(is.null(ontologyList[[j]])) next
      ontology[[which(names(ontology) == ontologyList[[j]][1, ]$term)]] <- ontologyList[[j]]
    }
    return(ontology)
  })
  cellular_locations_list <- lapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    unname(sapply(metabolite$biological_properties$cellular_locations, function(x) {
      if(stringr::str_detect(x[[1]], "\n")) return(NA)
      else return(x[[1]])
    }))
  })
  biospecimen_locations_list <- lapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    unname(sapply(metabolite$biological_properties$biospecimen_locations, function(x) {
      if(stringr::str_detect(x[[1]], "\n")) return(NA)
      else return(x[[1]])
    }))
  })
  tissue_locations_list <- lapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    unname(sapply(metabolite$biological_properties$tissue_locations, function(x) {
      if(stringr::str_detect(x[[1]], "\n")) return(NA)
      else return(x[[1]])
    }))
  })
  pathway_list <- lapply(1:num, function(i) {
    #print(paste0(i, " / ", num))
    metabolite <- metabolitesList[[i]]
    if(length(metabolite$biological_properties$pathways) > 200) return(NULL) # current metabolite
    pathway <- purrr::list_rbind(lapply(1:length(metabolite$biological_properties$pathways), function(j) {
      tmpRow <- metabolite$biological_properties$pathways[[j]]
      if(!is.list(tmpRow)) return(NULL)
      tibble::tibble(name = ifelse(length(tmpRow$name) == 1, tmpRow$name[[1]], ""),
                     smpdb_id = ifelse(length(tmpRow$smpdb_id) == 1, tmpRow$smpdb_id[[1]], ""),
                     kegg_map_id = ifelse(length(tmpRow$kegg_map_id) == 1, tmpRow$kegg_map_id[[1]], ""))
    }))
    if(nrow(pathway) == 0) return(NULL)
    return(pathway)
  })
  diseases_list <- lapply(1:num, function(i) {
    #print(paste0(i, " / ", num))
    metabolite <- metabolitesList[[i]]
    diseases <- purrr::list_rbind(lapply(1:length(metabolite$diseases), function(j) {
      diseaseTmp <- metabolite$diseases[[j]]
      if(!is.list(diseaseTmp)) return(NULL)
      tibble::tibble(name = ifelse(length(diseaseTmp$name) == 1, diseaseTmp$name[[1]], ""),
                     omim_id = ifelse(length(diseaseTmp$omim_id) == 1, diseaseTmp$omim_id[[1]], ""),
                     pubmed_id = list(sapply(1:length(diseaseTmp$references), function(k) {
                       ifelse(length(diseaseTmp$references[[k]]$pubmed_id) == 1, diseaseTmp$references[[k]]$pubmed_id[[1]], "")
                     })))
    }))
    if(nrow(diseases) == 0) return(NULL)
    return(diseases)
  })
  foodb_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$foodb_id) == 1, metabolite$foodb_id[[1]], NA)
  })
  kegg_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$kegg_id) == 1, metabolite$kegg_id[[1]], NA)
  })
  chemspider_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$chemspider_id) == 1, metabolite$chemspider_id[[1]], NA)
  })
  chebi_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$chebi_id) == 1, metabolite$chebi_id[[1]], NA)
  })
  pubchem_compound_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$pubchem_compound_id) == 1, metabolite$pubchem_compound_id[[1]], NA)
  })
  pdb_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$pdb_id) == 1, metabolite$pdb_id[[1]], NA)
  })
  biocyc_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$biocyc_id) == 1, metabolite$biocyc_id[[1]], NA)
  })
  drugbank_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$drugbank_id) == 1, metabolite$drugbank_id[[1]], NA)
  })
  phenol_explorer_compound_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$phenol_explorer_compound_id) == 1, metabolite$phenol_explorer_compound_id[[1]], NA)
  })
  wikipedia_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$wikipedia_id) == 1, metabolite$wikipedia_id[[1]], NA)
  })
  knapsack_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$knapsack_id) == 1, metabolite$knapsack_id[[1]], NA)
  })
  bigg_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$bigg_id) == 1, metabolite$bigg_id[[1]], NA)
  })
  metlin_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$metlin_id) == 1, metabolite$metlin_id[[1]], NA)
  })
  vmh_id_vec <- sapply(1:num, function(i) {
    metabolite <- metabolitesList[[i]]
    ifelse(length(metabolite$vmh_id) == 1, metabolite$vmh_id[[1]], NA)
  })
  protein_associations_list <- lapply(1:num, function(i) {
    #print(paste0(i, " / ", num))
    metabolite <- metabolitesList[[i]]
    protein_associations <- purrr::list_rbind(lapply(1:length(metabolite$protein_associations), function(j) {
      tmp <- metabolite$protein_associations[[j]]
      if(!is.list(tmp)) return(NULL)
      tibble::tibble(protein_accession = ifelse(length(tmp$protein_accession) == 1, tmp$protein_accession[[1]], ""),
                     name = ifelse(length(tmp$name), tmp$name[[1]], ""),
                     uniprot_id = ifelse(length(tmp$uniprot_id) == 1, tmp$uniprot_id[[1]], ""),
                     gene_name = ifelse(length(tmp$gene_name) == 1, tmp$gene_name[[1]], ""),
                     protein_type = ifelse(length(tmp$protein_type) == 1, tmp$protein_type[[1]], ""))
    }))
    if(nrow(protein_associations) == 0) return(NULL)
    return(protein_associations)
  })

  tibble::tibble(accession = accession_vec, secondary_accession = secondary_accession_list, name = name_vec,
                 synonyms = synonyms_list, chemical_formula = chemical_formula_vec,
                 average_molecular_weight = average_molecular_weight_vec, monisotopic_molecular_weight = monisotopic_molecular_weight_vec,
                 iupac_name = iupac_name_vec, traditional_iupac = traditional_iupac_vec,
                 cas_registry_number = cas_registry_number_vec,
                 smiles = smiles_vec, inchi = inchi_vec, inchikey = inchikey_vec,
                 taxonomy = taxonomy_list, ontology = ontology_list,
                 cellular_locations = cellular_locations_list, biospecimen_locations = biospecimen_locations_list, tissue_locations = tissue_locations_list,
                 pathway = pathway_list, diseases = diseases_list,
                 foodb_id = foodb_id_vec, kegg_id = kegg_id_vec, chemspider_id = chemspider_id_vec,
                 chebi_id = chebi_id_vec, pubchem_compound_id = pubchem_compound_id_vec, pdb_id = pdb_id_vec,
                 biocyc_id = biocyc_id_vec, drugbank_id = drugbank_id_vec, phenol_explorer_compound_id = phenol_explorer_compound_id_vec,
                 wikipedia_id = wikipedia_id_vec, knapsack_id = knapsack_id_vec, bigg_id = bigg_id_vec,
                 metlin_id = metlin_id_vec, vmh_id = vmh_id_vec,
                 protein_associations = protein_associations_list)
}
