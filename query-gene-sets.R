# query built gene sets 

library(dplyr)
library(magrittr)

## Read built gene sets
geo.disease.hs <- read.table("geo-hs-disease-gene-sets.tsv", sep = "\t", header = T, stringsAsFactors = F)
geo.ko.mm <- read.table("geo-mm-ko-gene-sets.tsv", sep = "\t", header = T, stringsAsFactors = F)
harmonizome.disease.hs <- read.table("harmonizome-hs-disease-gene-sets.tsv", sep = "\t", header = T, stringsAsFactors = F)

## Functions

### List GEO Datasets
# Lists compiled datasets by disease and dataset ID with unique gene counts.
#
# sig.only    Only count genes with FDR < 0.05. Default is FALSE.
# pos.only    Only count genes with log2FC > 0. Default is FALSE.
# by.dataset  Separate counts by dataset. Default is TRUE.
#
# returns a tibble
#
# examples: listGeoSets()
#           listGeoSets(T)
#
listGeoSets <- function(sig.only = FALSE, pos.only = FALSE, by.dataset = TRUE){
  
  geo.disease.hs %>%
    dplyr::filter(if (pos.only) log2fc > 0 else is.numeric(log2fc)) %>%
    dplyr::filter(if (sig.only) fdr < 0.05 else is.numeric(fdr)) %>%
    group_by(dataset_annot) %>%
    {if (by.dataset) count(., dataset_id, name = "unique_genes")
      else summarise(., unqiue_genes = n_distinct(symbol))
    }
}

### Genes By GEO Dataset
# Retrieve genes and associated data for a given dataset.
#
# dataset.id  ID of dataset for genes.
# sig.only    Only count genes with FDR < 0.05. Default is FALSE.
# pos.only    Only count genes with log2FC > 0. Default is FALSE.
#
# returns a datafame of symbol, log2fc, pv, fdr and sv
#
# Details: sv are "standardizedValues" calculated as -log10(pv) * sign(log2fc).
#   Theses values can be used for GSEA ranking or defining positive/negative significantly
#   differentially expressed gene sets.
#
# examples: genesByGeoSet("GSE3467")
#           genesByGeoSet("GSE3467", T,T)
#
genesByGeoSet <- function(dataset.id = NULL, sig.only = FALSE, pos.only = FALSE){
  
  if(is.null(dataset.id)){
    stop("Must provide a dataset ID. See listGeoSets().")
  }
  
  geo.disease.hs %>%
    dplyr::filter(dataset_id == dataset.id) %>%
    dplyr::filter(if (pos.only) log2fc > 0 else is.numeric(log2fc)) %>%
    dplyr::filter(if (sig.only) fdr < 0.05 else is.numeric(fdr)) %>%
    dplyr::select(symbol, log2fc, pv, fdr, sv)
}

### List Harmonizome Datasets
# Lists compiled datasets by disease and dataset ID with unique gene counts.
# Details: "standardizedValues" provided by Harmonizome are -log10(p-value) * sign(log2FC).
#   Theses values can be used for GSEA ranking or defining positive/negative significantly
#   differentially expressed gene sets.
#
# pos.only    Only count genes with standardizedValue > 0. Default is FALSE.
# by.dataset  Separate counts by dataset. Default is TRUE.
#
# returns a tibble
#
# examples: listHarmonizomeSets()
#           listHarmonizomeSets(T)
#
listHarmonizomeSets <- function(pos.only = FALSE, by.dataset = TRUE){
  
  harmonizome.disease.hs %>%
    dplyr::filter(if (pos.only) sv > 0 else is.numeric(sv)) %>%
    group_by(dataset_annot) %>%
    {if (by.dataset) count(., dataset_id, name = "unique_genes")
      else summarise(., unqiue_genes = n_distinct(symbol))
    }
}

### Genes By Harmonizome Dataset
# Retrieve genes and associated data for a given dataset.
# Details: "standardizedValues" provided by Harmonizome are -log10(p-value) * sign(log2FC).
#   Theses values can be used for GSEA ranking or defining positive/negative significantly
#   differentially expressed gene sets.
#
# dataset.id  ID of dataset for genes.
# pos.only    Only count genes with standardizedValue > 0. Default is FALSE.
#
# returns a datafame of symbol and standardizedValue
#
# examples: genesByHarmonizomeSet("GSE3467")
#           genesByHarmonizomeSet("GSE3467", T)
#
genesByHarmonizomeSet <- function(dataset.id = NULL, pos.only = FALSE){
  
  if(is.null(dataset.id)){
    stop("Must provide a dataset ID. See listGeoSets().")
  }
  
  harmonizome.disease.hs %>%
    dplyr::filter(dataset_id == dataset.id) %>%
    dplyr::filter(if (pos.only) sv > 0 else is.numeric(sv)) %>%
    dplyr::select(symbol, sv)
}