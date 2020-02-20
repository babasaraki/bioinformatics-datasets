# query built gene sets 

library(dplyr)
library(magrittr)

## Read built gene sets
disease.hs <- read.table("geo-hs-disease-gene-sets.tsv", sep = "\t", header = T, stringsAsFactors = F)
ko.mm <- read.table("geo-mm-ko-gene-sets.tsv", sep = "\t", header = T, stringsAsFactors = F)

## Functions

### List Datasets
# Lists compiled datasets by disease and dataset ID with unique gene counts.
#
# sig.only    Only count genes with FDR < 0.05. Default is FALSE.
# pos.only    Only count genes with log2FC > 0. Default is FALSE.
# by.dataset  Separate counts by dataset. Default is TRUE.
#
# returns a tibble
#
# examples: listDatasets()
#           listDatasets(T)
#
listDatasets <- function(sig.only = FALSE, pos.only = FALSE, by.dataset = TRUE){
  
  disease.hs %>%
    dplyr::filter(if (pos.only) log2fc > 0 else is.numeric(log2fc)) %>%
    dplyr::filter(if (sig.only) fdr < 0.05 else is.numeric(fdr)) %>%
    group_by(dataset_annot) %>%
    {if (by.dataset) count(., dataset_id, name = "unique_genes")
      else summarise(., unqiue_genes = n_distinct(symbol))
    }
}

### Genes By Dataset
# Retrieve genes and associated data for a given dataset.
#
# dataset.id  ID of dataset for genes.
# sig.only    Only count genes with FDR < 0.05. Default is FALSE.
# pos.only    Only count genes with log2FC > 0. Default is FALSE.
#
# returns a datafame of symbol, log2fc, pv and fdr
#
# examples: genesByDataset("GSE6357")
#           genesByDataset("GSE6357", T,T)
#
genesByDataset <- function(dataset.id = NULL, sig.only = FALSE, pos.only = FALSE){
  
  if(is.null(dataset.id)){
    stop("Must provide a dataset ID. See listDatasets().")
  }
  
  disease.hs %>%
    dplyr::filter(dataset_id == dataset.id) %>%
    dplyr::filter(if (pos.only) log2fc > 0 else is.numeric(log2fc)) %>%
    dplyr::filter(if (sig.only) fdr < 0.05 else is.numeric(fdr)) %>%
    dplyr::select(symbol, log2fc, pv, fdr)
}