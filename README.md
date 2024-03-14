# consensusClustR
R Package for iterative bootstrapped consensus clustering of scRNA-seq Data

Use the consensusClust function to cluster scRNA-seq data over a number of resolutions, using bootstrap sampling and statistical testing to reject clusterings which may have been produced by sampling noise (i.e. by chance). By setting the iterate parameter to TRUE, consensusClust will subcluster your matrix (using principal components newly computed for each cluster) until no statistically distinct clusters are found.

Briefly this function:
1) Either processes a provided count matrix (normalising, regressing out unwanted variables), or extracts this data from a Seurat or SingleCellExperiment object.
2) Computes principal components (pcs) of this matrix.
3) Forms a number of bootstraps (samples with replacement) of the pc matrix, to replicate resampling cells from the same samples, and computes shared nearest neighbour graphs for these boostraps.
4) Clusters these bootstrap nearest neighbour graphs over a number of resolutions, using either leiden or louvain community detection algorithms (leiden by default). 
5) Produces a distance matrix between cells based on how often each pair of cells are clustered together in the bootstraps. By default, only the assignments from the resolution which produces the most distinct clusters (highest silhouette score) for each bootstrap, is used to determine co-clustering distances. Results from all resolutions can be used to produce this distance matrix if desired, by changing the 'mode' argument to 'granular'. 
6) Clusters the cell-cell co-clustering distance matrix over the same range of resolutions, choosing the assignments giving the most distinct clusters (highest silhouette score), and merging clusters whose cells are often co-clustered (bootstrap stability < 0.4).
7) Tests if the distinctness if output clustering, would be unlikely to occur if all cells really came from one subpopulation, and the clusters have solely been formed due to the noise caused by random incomplete sampling of RNAs from the cells and other technical factors (the null hypothesis). This is done by creating a normally distributed null distribution of silhouette scores from clustering of new count matrices, simulated based on the real count matrix under the assumption that all cells come from one cluster. The probability of finding the silhouette score of the real data is then tested under this null distribution. If pval < alpha (0.05 by default), clustering is retained. If pval >= alpha, clustering is rejected, and all cells are assigned to the same cluster.
8) If iterate=TRUE, repeats this process on all identified subclusters, until no statistically distinct clusters are detected. 
9) Returns the user a list containing:
'assignments': the cluster (and subcluster if applicable) assignments for each cell as a character vector in the same order as colnames of the input data matrix.
'clusterDendrogram': a dendrogram showing the relatedness of output cluster assignments, based on the co-clustering distance matrix.
'clustree': A visualisation of subcluster-parent cluster relationships, made using the clustree package.

Input can be a count matrix, or a Seurat/SingleCellExperiment object from which the function will extract the counts, normalised counts, variable features, and pca cell embeddings, if they are present in the object.

Installation:
```
#Install devtools for installs from github 
if(!require("devtools")){
  install.packages("devtools")
}
#Install consensusClustR
devtools::install_github("AndyCGraham/consensusClustR")
```
  
Example use on a gene by cell scRNA-seq count matrix called counts
```
library(consensusClustR)
#Using 5 PCs:
results <- consensusClust(counts, pcNum = 5, threads = 12)

#With iterative subclustering
results <- consensusClust(counts, iterate=TRUE, pcNum = 5)

#With biocparallel parallel processing - for Linux or Mac
results <- consensusClust(counts, iterate=TRUE, pcNum = 5, BPPARAM = MulticoreParam(7,RNGseed = 123))

#With biocparallel parallel processing - for windows (cannot use multicore forked processing)
results <- consensusClust(counts, iterate=TRUE, pcNum = 5, BPPARAM = SnowParam(7,RNGseed = 123))

#Without known pcNum - will be estimated by the function
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
results$clustree
```
