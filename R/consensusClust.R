#' Perform consensus clustering on a scRNA-seq matrix
#'
#' @param data matrix of scaled counts, or a Seurat or SingleCellExperiment object from which these can be extracted
#' @param pcNum number of principal components to consider when finding nearest neighbours
#' @param nboots number of boostraps to perform
#' @param clusterFun Which of the igraph community detection algorithms to use for clustering: "leiden" or "louvain"
#' @param bootSize fraction of total cell number to use per bootstrap
#' @param resRange Numeric vector of resolutions to cluster over
#' @param pcaMethod Method for computing pca (specify irbla for faster, less accurate PCA)
#' @param kNum number of nearest neighbours for knn graph construction
#' @param mode How to deal with cluster assignments from different resolutions, either pick the one with highest silhouette score ('robust')
#' or use all when deciding consensus clusters ('granular')
#' @param threads How many cpus to use in parallel
#' @param assay Name of seurat or SingleCellExperiment assay, if providing an object rather than scaled counts
#' 
#' @return vector: character vector of consensus clustering assignments
#' @export
#' 
#' @importFrom bluster neighborsToSNNGraph
#' @importFrom igraph cluster_leiden E
#' @importFrom parallel mclapply
#' @importFrom dbscan kNN
#' @importFrom cluster  silhouette
#' @importFrom RcppXPtrUtils  cppXPtr
#' @importFrom parallelDist parDist
#' @importFrom SummarizedExperiment assay
#' 
consensusClust <- function(data, pcNum=15, nboots=200, clusterFun="leiden", bootSize=0.8, resRange = seq.int(0.05, 1, by = 0.05),  
                            pcaMethod = "irbla", kNum=30, mode = "robust", threads=1, assay="RNA", ...) {
  
  #Check input is correct
  stopifnot("`data` must be a matrix, sparse matrix of type dgCMatrix, seurat object, or single-cell experiment object." = 
              any((class(data)=="Seurat") | (class(data)=="SingleCellExperiment") | (class(data)[1]=="matrix") |
               (class(data)=="dgCMatrix") ) )
  stopifnot("`pcNum` must be a positive integer, smaller than cell number." = 
              all((length(pcNum)==1) & (pcNum%%1==0) & (pcNum > 0) & (pcNum < ncol(data))) )
  stopifnot("`nboots` must be a positive integer." = 
              all((nboots%%1==0) & (nboots > 0) & (length(nboots)==1)) )
  stopifnot("`clusterFun` must be a either leiden or louvain." = 
              clusterFun %in% c("leiden", "louvain") )
  stopifnot("`bootSize` must be a single number." = 
              all((length(bootSize)==1) & (is.numeric(bootSize))) )
  stopifnot("`resRange` must be a numeric vector." = 
              is.numeric(resRange) )
  stopifnot("`pcaMethod` must be either 'irbla' or 'svd'." = 
              pcaMethod %in% c("irbla", "svd") )
  stopifnot("`kNum` must be a positive integer." = 
              all((kNum%%1==0) & (kNum > 0)))
  stopifnot("`mode` must be either 'robust' or 'granular'." = 
              mode %in% c("robust", "granular") )
  
  if(class(data)=="Seurat"){
    data = data[[assay]]$scale.data
  } else if(class(data)=="SingleCellExperiment"){
    data = assay(sce, assay)
  }
  
  #Get cluster assignments for bootstrapped selections of cells
  clustAssignments = mclapply(1:nboots, \(boot){
    getClustAssignments( data[,sample(colnames(data), bootSize*length(colnames(data)), replace = TRUE)],
    pcNum = pcNum, resRange = resRange,kNum=kNum, clusterFun = clusterFun, pcaMethod = pcaMethod, cellOrder = colnames(data), mode = mode)
    }, mc.cores = threads )
  
  ##Calculate cell-cell distance matrix based on these co-clusterings
  #Get clustering assignments in the right format (dataframe)
  clustAssignments = do.call(cbind, clustAssignments)
  rownames(clustAssignments) = colnames(data)
  
  #Alter NAs to -1, which will be ignored by the distance function
  clustAssignments[is.na(clustAssignments)] = -1
  
  #Function to compute jaccard similarity between non-binary sets (ignoring -1s), implemented in C++
  jaccardDist <- cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) { 
                              float overlap, U;
                              double jaccard;
                              overlap = arma::accu((A == B) && (A != -1));
                              U = arma::accu((A != -1) && (B != -1));
                              jaccard = overlap / U;
                              return jaccard; }",
                                        depends = c("RcppArmadillo"))
  
  #Calculate distance matrix
  jaccardDist = 1-parDist(clustAssignments, method="custom", func = jaccardDist, threads = threads) 
  
  ##Adjacency graph from similarity matrix
  knn = kNN(jaccardDist, k = kNum, search = "kdtree")$id
  snnGraph = neighborsToSNNGraph(knn, type = "rank")
  
  ##Cluster adjacency graph at different resolutions
  if(clusterFun=="leiden"){
    clustAssignments = mclapply(resRange, \(res){
      cluster_leiden(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                     beta = 0.01, n_iterations = 2)$membership 
    }, mc.cores = threads)
  } else if(clusterFun=="louvain"){
    clustAssignments = mclapply(resRange, \(res){
      cluster_louvain(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                     beta = 0.01, n_iterations = 2)$membership 
    }, mc.cores = threads)
  }
  
  #Assess silhouette score per cluster for every resolution
  clustScores = rank(unlist(mclapply(clustAssignments, \(res){
    if(length(unique(res))>1){
      mean(aggregate(sil_width ~ cluster, silhouette(x = res, dist = jaccardDist), mean)[,2], na.rm = T)
    } else {
      0.15
    }
  }, mc.cores = threads)), ties.method ="first")
  
  #Take assignments from resolution with highest score
  finalAssignments = clustAssignments[[which(clustScores == max(clustScores))]]
  
  #Compute cluster dendrogram
  if(length(unique(finalAssignments)) > 1){
    dendrogram = determineHierachy(as.matrix(jaccardDist), finalAssignments)
  } else {
    dendrogram = NULL
  }
  
  
  return(list(assignments = finalAssignments, res = resRange[which(clustScores == max(clustScores))], clusterDendrogram = dendrogram))
  
}




#' Conducts PCA, and clusters cells based on resulting PCs, returning assignments as a dataframe
#' @param scaledCounts  Matrix of scaled counts
#' @param pcNum number of principal components to consider when finding nearest neighbours
#' @param clusterFun Which of the igraph community detection algorithms to use for clustering: "leiden" or "louvain"
#' @param resRange Numeric vector of resolutions to cluster over
#' @param pcaMethod Method for computing pca (specify irbla for faster, less accurate PCA)
#' @param kNum number of nearest neighbours for knn graph construction
#' @param mode How to deal with cluster assignments from different resolutions, either pick the one with highest silhouette score ('robust')
#' or use all when deciding consensus clusters ('granular')
#' @param cellOrder Vector of non-boostrapped cellnames in the order of the non-boostrapped data matrix
#' 
#' @export
#' 
#' @importFrom bluster clusterRows  NNGraphParam  approxSilhouette
#' @importFrom irlba  irlba
#' @importFrom cluster  silhouette
#' 
getClustAssignments <- function(scaledCounts, pcNum, clusterFun="leiden", resRange, pcaMethod, kNum, mode = "robust", cellOrder, ...) {
  
  if(pcaMethod == "irbla"){
    pcs = irlba::irlba(scaledCounts, pcNum)$v
  } else {
    pcs = svd(scaledCounts, 
              nv = pcNum)$v
    }
  
  #Cluster adjacency matrix, returning assignments named by cell name
  clustAssignments = purrr::map(resRange, function(res){
    assignments = setNames(
      clusterRows(pcs, BLUSPARAM=NNGraphParam(k=kNum, cluster.fun=clusterFun, cluster.args=list(resolution=res))),
      colnames(scaledCounts) 
    )
    
    if(mode == "robust"){
      if(length(unique(assignments))>1){
        clustScore = mean(aggregate(width ~ cluster, approxSilhouette(pcs, assignments), mean)[,2], na.rm = T)
      } else {
        clustScore = 0.15
      }
      
      }
    
    assignments = assignments[match(cellOrder, names(assignments))]
    
    if(mode == "robust"){
      return(list(assignments, clustScore))
    } else{
    return(assignments)
      }
  })
  
  if(mode == "robust"){
    clustScores = rank(unlist(purrr::map(clustAssignments, \(res) res[[2]] )), ties.method ="first")
    clustAssignments = clustAssignments[[which(clustScores == max(clustScores))]][[1]]
  } else {
    clustAssignments = do.call(cbind, clustAssignments)
  }
  
  return(clustAssignments)
}

#' Conducts PCA, and clusters cells based on resulting PCs, returning assignments as a dataframe
#' @param distanceMatrix  Matrix of cell-cell distances
#' @param assignments clustering assignments
#' 
#' @export
determineHierachy <- function(distanceMatrix, assignments) {
  # Create a distance matrix for the clusters
  clusterDistanceMatrix <- matrix(0, nrow = length(unique(assignments)), 
                                    ncol = length(unique(assignments)))
  
  # Fill the cluster distance matrix with the distances between cluster centroids
  for (clust1 in unique(assignments)) {
    for (clust2 in unique(assignments)[!unique(assignments) == clust1]) {
      if(clusterDistanceMatrix[clust1, clust2] == 0){
        clust1Samples = which(assignments == clust1)
        clust2Samples = which(assignments == clust2)
        dist = mean(distanceMatrix[clust1Samples, clust2Samples], na.rm = T)
        clusterDistanceMatrix[clust1, clust2] = dist
        clusterDistanceMatrix[clust2, clust1] = dist
      }
    }
  }
  
  clusterDistanceMatrix = as.dist(clusterDistanceMatrix)
  # Perform hierarchical clustering on the cluster distance matrix
  hclustRes <- hclust(clusterDistanceMatrix, method = "complete")
  
  return(hclustRes)
}

