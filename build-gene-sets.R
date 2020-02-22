# build gene sets from public data

library(dplyr)
library(magrittr)

## List of dirs to process
buildGeoSet <- function(){
  dir.list <- c("geo-hs-disease-gene-sets", "geo-mm-ko-gene-sets")
  
  lapply(dir.list, function(d){
    file.list <- list.files(d, "\\.tsv")
    
    df.data <- data.frame(symbol = character(),
                          log2fc = numeric(),
                          pv = numeric(),
                          fdr = numeric(),
                          dataset_id = character(),
                          dataset_annot = character())
    
    lapply(file.list, function(f){
      fn.a <- strsplit(f,".tsv")[[1]]
      fn.sp <- strsplit(fn.a, "__")[[1]]
      id <- fn.sp[1]
      annot <- fn.sp[2]
      data <- read.table(file.path(d, f) , header = T, sep = "\t", stringsAsFactors = F)
      
      #Filter for significant p-values and non-blank symbols
      sig.data <- data %>%
        dplyr::filter(!Gene.symbol == "" & P.Value < 0.05) %>%
        dplyr::select(Gene.symbol, logFC, P.Value, adj.P.Val) %>%
        set_colnames(c("symbol", "log2fc", "pv", "fdr")) 
      #Select median probe to handle duplicates
      med.data <- sig.data %>%
        dplyr::group_by(symbol) %>%
        dplyr::summarise(log2fc = median(log2fc)) %>%
        mutate(med.value = TRUE)
      this.data <- merge(sig.data, med.data, by=c("symbol", "log2fc")) %>%
        arrange(desc(log2fc)) %>%
        mutate(sv = sign(log2fc)*-log(pv,10)) %>% ## Add standardValue
        dplyr::select(-med.value)
      
      this.data$dataset_id = id
      this.data$dataset_annot = annot
      df.data <<- rbind(df.data, this.data)
    })
    
    write.table(df.data, paste0(d, ".tsv"), sep = "\t", row.names = F)
  })
}

## Harmonizome
buildHarmonizomeSets <- function() {
  gams <- read.table("gene_attribute_matrix_standardized.txt", sep="\t", stringsAsFactors = F, header = T)
  attr.list <- read.table("attribute_list_entries.txt", sep="\t", stringsAsFactors = F, header = T, quote="\"")
  
  df.data <- data.frame(symbol = character(),
                        sv = numeric(),
                        dataset_id = character(),
                        dataset_annot = character())
  
  lapply(names(gams)[-c(1:3)], function(d){
    annot.gse <- substring(d,2)
    annot <- attr.list[which(attr.list$GSE == annot.gse),"Disease"]
    id <- paste0("GSE",annot.gse)
    
    #Filter for significant p-values
    this.data <- gams %>%
      dplyr::filter(!!as.name(d) > -log(0.05,10) |
                      !!as.name(d) < log(0.05,10)) %>%
      dplyr::arrange(desc(!!as.name(d))) %>%
      dplyr::select("GeneSym", d) %>%
      set_colnames(c("symbol","sv"))

    this.data$dataset_id = id
    this.data$dataset_annot = annot
    df.data <<- rbind(df.data, this.data)
  })
  
  write.table(df.data, "harmonizome-hs-disease-gene-sets.tsv", sep = "\t", row.names = F)

}