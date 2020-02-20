# bioinformatics-datasets
public datasets for testing, benchmarking and comparison

Note: The PDF Supplemental tables provide an annotated list of GEO datasets for disease-specific and KO experiments. These can be processed by the protocol below to generate DE gene sets (see TSV files in subfolders).

## GEO Protocol
1. Go to https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE22873. Replace GEO ID.
2. Expand Samples
3. Define two groups, e.g., cancer/normal or ko/wt
4. Assign samples to groups
5. Click on Top 250
6. Click on Save all results
7. Save as TSV into appropriate subfolder

## Compiling datasets
Run `build-gene-sets.R` to generate a top-level TSV file compiled from each subfolder of collected datasets.

## Get Sets of DE Genes
Source 'query-gene-sets.R` to access functions to extract gene sets and associated data values per dataset.

For example, `listDatasets(T)` will summarize the unique DE genes with FDR < 0.5 per dataset per disease. And `genesByDataset("GSE6357", T, T)` will return a dataframe of DE genes with FDR < 0.05 and log2FC > 0.
