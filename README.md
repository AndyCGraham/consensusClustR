# consensusClustR
R Package for iterative bootstrapped consensus clustering of scRNA-seq Data

Use the consensusClust function to cluster scRNA-seq data over a number of resolutions, using bootstrap sampling and statistical testing to reject clusterings which may have been produced by sampling noise (i.e. by chance). By setting the iterate parameter to TRUE, consensusClust will subcluster your matrix (using principal components newly computed for each cluster) until no statistically distinct clusters are found.

Briefly this function:
1) Either processes a provided count matrix (normalising, regressing out unwanted variables), or extracts this data from a Seurat or SIngleCellExperiment object.
2) Computes principal components (pcs) of this matrix.
3) Forms a number of bootstraps (samples with replacement) of the pc matrix, to replicate resampling cells from the same samples, and computes shared nearest neighbour graphs for these boostraps.
4) Clusters these bootstrap nearest neighbout graphs over a number of resolutions, using either leiden or louvain community detection algorithms (leiden by default). 
5) Produces a distance matrix between cells based on how often each pair of cells are clustered together, in all the bootstraps. By default, only the assignments from the resolution which produces the most distinct clusters (highest silhouette score) for each bootstrap, is used to determine co-clustering distances. Results from all resolutions can be used to produce this distance matrix if desired, by changing the 'mode' parameter to 'granular'. 
6) Clusters the cell-cell co-clustering distance matrix over the same range of resolutions, choosing the assignments giving the most distinct clusters (highest silhouette score), and merging clusters whose cells are often co-clustered (bootstrap stability < 0.4).
7) If iterate=TRUE, repeats this process on all subclusters, until no statistically distinct clusters are detected. 
8) Returns the user a list containing:
'assignments': the cluster (and subcluster if applicable) assignments for each cell as a character vector in the same order as colnames of the input data matrix.
'clusterDendrogram': a dendrogram showing the relatedness of output cluster assignments, based on the co-clustering distance matrix.
'clustree': A visualisation of subcluster-parent cluster relationships, made using the clustree package.

Input can be a count matrix, or a Seurat/SingleCellExperiment object from which the function will extract the counts, normalised counts, variable features, and pca cell embeddings, if it's already present.

Installation:
```
if(!require("remotes")){
  install.packages("remotes")
}
remotes::install_github("AndyCGraham/consensusClustR")
```
  
Example use on a gene by cell scRNA-seq count matrix called counts
```
library(consensusClustR)
#Using 5 PCs:
results <- consensusClust(counts, pcNum = 5, threads = 12)

#With iterative subclustering
results <- consensusClust(counts, iterate=TRUE, pcNum = 5)

#With biocparallel parallel processing
results <- consensusClust(counts, iterate=TRUE, pcNum = 5, BPPARAM = MulticoreParam(7,RNGseed = 123))

#Without known pcNum
results <- consensusClust(counts, iterate=TRUE, BPPARAM = MulticoreParam(7,RNGseed = 123))
```
  
Example use on a Seurat object called data:
```
library(consensusClustR)

#Run consensus clustering
results <- consensusClust(data, iterate=TRUE, BPPARAM = MulticoreParam(7,RNGseed = 123))

#Assign to Seurat Obj metadata
data$consensusClusts = results$assignments

#Make consensusClusts the cell identities
Idents(data) = data$consensusClusts

#Show cluster dendrogram plot, to visualise cluster relatedness
plot(results$clusterDendogram)

#Show the clustree representation of subcluster-parent cluster relationships
restults$clustree
```
