# consensusClustR
R Package for bootstrapped consensus clustering of scRNA-seq Data

Use the consensusClust function to cluster scRNA-seq data over a number of resolutions and bootstraps. 

Briefly this function:
1) Forms a number of bootstraps (samples with replacement) of a scRNA-seq data matrix
2) Performs PCA and computes shared nearest neighbour distances for these boostraps, then clusters them over a number of resolutions, using either leiden of louvain clustering. 
3) Produce a distance matrix between cells based on how often each pair of cells are clustered together, in all the bootstraps. By default, only the assignments from the resolution which produces the most distinct clusters (highest silhouette score) per each bootstrap are used to determine co-clustering distances. Results from all resolutions can be used to produce this distance matrix if desired, by changing 'mode' parameter to 'granular'. 
4) The cell-cell co-clustering distance matrix is clustered over the same range of resolutions. 
5) The user is returned a list of the assignments from the resolution which produced the most distinct clusters (highest silhouette score) as a character vector in the same order as colnames of the input data matrix: 'assignments', the value of the resolution which produced the most distinct clusters: 'res', and a dendrogram showing the relatedness of output cluster assignments, based on the co-clustering distance matrix: 'clusterDendrogram'.

Input can be a gene by cell data matrix, a Seurat object, or a SingleCellExperiment object. For a Seurat or SingleCellExperiment object, change the 'assay' parameter to reflect the name of the assay you want to pull data from.

Installation:
```
if(!require("remotes")){
  install.packages("remotes")
}
remotes::install_github("AndyCGraham/consensusClustR")
```

Example use on a gene by cell data matrix called data, using 5 PCs and 12 cpus:

```
results <- consensusClust(data, pcNum = 5, threads = 12)
```
