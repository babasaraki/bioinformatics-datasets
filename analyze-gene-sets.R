# analyze gene sets 

## CONFIG VARS
n_genes_threshold = 400
dataset_source = "har" #geo or har
pos_only = TRUE
string = TRUE

## Libs
load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "tidyverse",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}


## Functions
source("query-gene-sets.R")

## Analysis

# Prep WP
#wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt", destpath = "_results")
wp.hs.gmt <-"_results/wikipathways-20200210-gmt-Homo_sapiens.gmt"
wp2gene <- clusterProfiler::read.gmt(wp.hs.gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME
wpid2name<-unique(wpid2name)
wpid2name <- wpid2name %>%
  mutate(name = substr(name, 1,60))

# Prep datasets
dataset.df <- {if(dataset_source == "geo") listGeoSets(T,pos_only) 
  else if(dataset_source == "har") listHarmonizomeSets(pos_only) } %>%
  dplyr::filter(unique_genes >= n_genes_threshold)
dataset.names <- paste(dataset.df$dataset_id, dataset.df$dataset_annot, sep = "__")
  
# Initialize
res.cnts <- data.frame(set = character(0),
                      genes = integer(0),
                      wp_pathways = integer(0)
) 

# Clear prior results
res.wp.files <- list.files(file.path("_results",dataset_source), pattern=".*wp\\.csv")
lapply(res.wp.files, function(fn){
  file.remove(file.path("_results",dataset_source,fn))
})

# Loop
lapply(dataset.names, function(d){
  d.id <- strsplit(d, "__")[[1]][1]
  
  ### get gene set
  d.genes <- {if(dataset_source == "geo") genesByGeoSet(d.id, T, T) 
    else if(dataset_source == "har") genesByHarmonizomeSet(d.id, T) } 
  d.gene.list <- na.omit(unique(d.genes$symbol))
  d.gene.list.entrez <- bitr(d.gene.list, fromType = "SYMBOL",
                             toType = c("ENTREZID"),
                             OrgDb = org.Hs.eg.db)

  ### STRING network
  if(string){
    library(RCy3)
    commandsRun(paste0('string protein query query="',paste(d.gene.list.entrez[,1],collapse = ","),
                       '" limit=0 cutoff=0.4'))
    createSubnetwork(edges = "all")
    net.gene.list <- getTableColumns(columns="query term")[,1]
    d.gene.list.entrez <- dplyr::filter(d.gene.list.entrez, SYMBOL %in% net.gene.list)
  }
  
  ## Filter and subset
  ewp.genes <- na.omit(d.gene.list.entrez$ENTREZID[1:n_genes_threshold])
    
  
  ### WikiPathways Analysis
  ewp <- clusterProfiler::enricher(
    ewp.genes,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05, #p.adjust cutoff
    minGSSize = 5,
    maxGSSize = 300,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  ewp <- DOSE::setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  head(ewp, 20)
  
  ### plots
  # dotplot(ewp, showCategory = 10, orderBy="GeneRatio")
  # dotplot(ewp, showCategory = 10, orderBy="p.adjust", x = ~-p.adjust)
  # emapplot(ewp.dec, showCategory = 40)
  # png(paste(d,"wikipathways_dot.png", sep = "_"), width = 2200, height = 1300, res=300)
  
  # print(dotplot(ewp, showCategory = 10, orderBy="GeneRatio") )
  # dev.off()
  
  # my_x_sort = "p.adjust"
  # my_result = ewp
  # df <- fortify(my_result, showCategory = 15)
  # df <- dplyr::mutate(df, x = -log10(eval(parse(text=my_x_sort))))
  # idx <- order(df[["x"]], decreasing = TRUE)
  # df$Description <- factor(df$Description, levels=rev(unique(df$Description[idx])))
  # print(ggplot(df, aes_string(x = "x", y="Description", size="Count", color="p.adjust")) +
  #         geom_point() +
  #         xlab("-log(p.adjust)") +
  #         scale_color_continuous(low="red", high="blue", name = my_x_sort, guide=guide_colorbar(reverse=TRUE)) +
  #         ## scale_color_gradientn(name = color, colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
  #         ylab(NULL) + ggtitle("") + theme_dose(12) + scale_size(range=c(3, 8)))
  # dev.off()
  
  ### saves
  write.csv(ewp, file = file.path("_results", dataset_source,paste(d,"wp.csv",sep="__")), row.names = F)
  wp.cnt <- nrow(ewp)
  
  ## Summary objects
  this.cnt <- data.frame(set = d,
                         genes = length(d.gene.list.entrez[,2]),
                         wp_pathways = wp.cnt
  )
  res.cnts <<- rbind(res.cnts, this.cnt)
  write.csv(res.cnts, file.path("_results", dataset_source,"summary_counts.csv"), row.names = F)
})



## Summarize
# collect pathway results in useful objects for plotting
res.wp.files <- list.files(file.path("_results",dataset_source), pattern=".*wp\\.csv")

res.mat <- data.frame(ID = character(), Title = character()) #for heatmap

lapply(res.wp.files, function(r.file) { 
  res.name <- unlist(strsplit(r.file, "__wp"))[1]
  
  res <- read.csv(file.path("_results",dataset_source, r.file), stringsAsFactors = F)
  this.res <- res %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, Description, p.adjust) %>%
    set_colnames(c("ID","Title", res.name))
  
  if(nrow(this.res) > 0){
    res.mat <<- res.mat %>%
      full_join(this.res, by=c("ID","Title"))
  }
})


##############
## Heatmap of pathway results
#########
library(pheatmap)  

# replace pvalues with negative log values.
# replace NA with 0 to avoid NA error; should remove from legend
# (optional) use filter by sum to reduce complexity of result display
res.mat.log <- res.mat[,3:ncol(res.mat)] %>%
  mutate_all(log10) %>%
  replace(is.na(.), 20) %>% ## workaround for NA error
  mutate_if(is.numeric, funs(-.)) %>%
  mutate(sum = rowSums(.)) %>%
  magrittr::set_rownames(res.mat[,2]) %>%
  rownames_to_column('temp') %>%  #temporarily stash rownames while filtering
  #filter(sum > 4) %>%  #NOTE: CHANGE TO ADJUST NUMBER OF PATHWAYS
  column_to_rownames('temp') %>% #restore rownames
  dplyr::select(-sum)

wp.heat <- as.matrix(res.mat.log)

# proportionally set font sizes
fontsize_row = 14 - nrow(wp.heat) / 25
fontsize_col = 14 - ncol(wp.heat) / 15

# plot
pheatmap(wp.heat, 
         scale = "none",
         col=colorRampPalette(brewer.pal(n = 9, name ="OrRd"))(6),
         fontsize_row=fontsize_row, 
         fontsize_col=fontsize_col,
         border_color=NA)

