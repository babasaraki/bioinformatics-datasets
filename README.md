# bioinformatics-datasets
public datasets for testing, benchmarking and comparison

Note: The PDF Supplemental tables provide an annotated list of GEO datasets for disease-specific and KO experiments. These can be processed by the protocol below to generate DE gene sets (see TSV files in subfolders).

## GEO Protocol
1. Go to https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE22873. Replace GEO ID.
1. Expand Samples
1. Define two groups, e.g., cancer/normal or ko/wt
1. Assign samples to groups
1. Click on Top 250
1. Confirm directionality, i.e., gene with increased expression in cancer vs normal should have a positive log2FC.
1. Click on Save all results
1. Save as TSV into appropriate subfolder

## Compiling datasets
1. Source `build-gene-sets.R` 
1. Run `buildGeoSets()` to generate a top-level TSV file compiled from each subfolder of collected datasets.

## Get Sets of DE Genes
Source 'query-gene-sets.R` to access functions to extract gene sets and associated data values per dataset.

For example, `listDatasets(T)` will summarize the unique DE genes with FDR < 0.5 per dataset per disease. And `genesByDataset("GSE6357", T, T)` will return a dataframe of DE genes with FDR < 0.05 and log2FC > 0.

## Harmonizome Protocol
1. Download [gene-attribute-matrix-standardized](https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/geodisease/gene_attribute_matrix_standardized.txt.gz) and [attribute-list](https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/geodisease/attribute_list_entries.txt.gz) txt.gz files.
1. Unzip files
1. Source `build-gene-sets.R` 
1. Run `buildHarminizomeSets()` to generate a top-level TSV file compiled from downloaded files.
1. Delete unzipped files
