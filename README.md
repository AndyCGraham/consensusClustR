# consensusClustR
R Package for bootstrapped consensus clustering of scRNA-seq Data

Use the consensusClust function to cluster scRNA-seq data over a number of resolutions and bootstraps. 

Briefly this function:
1) Forms a number of bootstraps (samples with replacement) of a principal component matrix, to replicate resampling cells from the same samples.
2) Computes shared nearest neighbour graphs for these boostraps.
3) Clusters these graphs over a number of resolutions, using either leiden or louvain community detection algorithms (leiden by default). 
4) Produces a distance matrix between cells based on how often each pair of cells are clustered together, in all the bootstraps. By default, only the assignments from the resolution which produces the most distinct clusters (highest silhouette score) for each bootstrap, is used to determine co-clustering distances. Results from all resolutions can be used to produce this distance matrix if desired, by changing the 'mode' parameter to 'granular'. 
4) Clusters the cell-cell co-clustering distance matrix over the same range of resolutions. 
5) Returns the user a list containing:
'assignments': the assignments from the resolution which produced the most distinct clusters (highest silhouette score) as a character vector in the same order as colnames of the input data matrix.
'res': the value of the resolution which produced the most distinct clusters.
'clusterDendrogram': a dendrogram showing the relatedness of output cluster assignments, based on the co-clustering distance matrix.

Input can be a cell by principal component cell embeddings matrix from principal component analysis (e.g. SeuratObj@reductions$pca@cell.embeddings), a Seurat object, or a SingleCellExperiment object, from which the function will extract the pca cell embeddings if it's already present.

Installation:
```
if(!require("remotes")){
  install.packages("remotes")
}
remotes::install_github("AndyCGraham/consensusClustR")
```
  
Example use on a cell by PC principal component matrix called pca, using 5 PCs and 12 cpus:
```
library(consensusClustR)
results <- consensusClust(pca, pcNum = 5, threads = 12)
```
  
Example use on a Seurat object called data, using 5 PCs and 12 cpus:
```
library(consensusClustR)
library(Seurat)

#Run PCA as normal with Seurat
data = RunPCA(data)

#Use elbow plot todetermine optimal PC number as normal
ElbowPlot(data, ndims = 50)

#Run consensus clustering
results <- consensusClust(data, pcNum = 5, threads = 12)

#Assign to Seurat Obj metadata
data$consensusClusts = results$assignments

#Make consensusClusts the cell identities
Idents(data) = data$consensusClusts

#Show cluster dendrogram plot, to visualise cluster relatedness
plot(results$clusterDendogram)
```
