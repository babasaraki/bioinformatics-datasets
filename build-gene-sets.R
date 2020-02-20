# build gene sets from public data

library(dplyr)
library(magrittr)

## List of dirs to process
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
    fn <- f
    fn.a <- strsplit(f,".tsv")[[1]]
    fn.sp <- strsplit(fn.a, "__")[[1]]
    id <- fn.sp[1]
    annot <- fn.sp[2]
    data <- read.table(file.path(d, fn) , header = T, sep = "\t", stringsAsFactors = F)
    this.data <- data %>%
      arrange(desc(logFC)) %>%
      dplyr::filter(!Gene.symbol == "" & P.Value < 0.05) %>%
      dplyr::select(Gene.symbol, logFC, P.Value, adj.P.Val) %>%
      set_colnames(c("symbol", "log2fc", "pv", "fdr"))
    this.data$dataset_id = id
    this.data$dataset_annot = annot
    df.data <<- rbind(df.data, this.data)
  })
  
  write.table(df.data, paste0(d, ".tsv"), sep = "\t", row.names = F)
})
