#' Performs iterative consensus clustering of a scRNA-seq matrix, to find robust, statistically distinct cell clusters.
#'
#' @param counts gene by cell matrix of scRNA-seq counts, or a Seurat or SingleCellExperiment object from which 
#' these can be extracted.
#' @param iterate Boolean specifying whether to continuously subcluster until no statistically distinct cluster are found. Default: FALSE.
#' @param alpha Significance threshold for difference between silhouette scores generated from clustering real vs null data.
#' @param pca A cell by pc pricipal component matrix. If left NULL, will be computed internally.
#' @param pcNum number of principal components to consider when finding nearest neighbours. If left as "find", willl be determined internally.
#' @param pcaMethod Method to compute PCA. Either irlba (default) - fast and memory efficient, 'prcomp' - base R prcomp, or 'svd' - base R svd.
#' @param scale Boolean specifying whether normalised counts be scaled prior to pca computation. Default: TRUE.
#' @param center Boolean specifying whether normalised counts be centered prior to pca computation. Default: TRUE.
#' @param sizeFactors Either a vector of precomputed sizefactors to normalize by, or a character vector to specify how to compute these (see
#' shifted_log_transform function from the transformGamPoi package for details). Default: 'deconvolution' - pooled size factor estimation using 
#' the scran package (requires that scran is installed on your machine).  
#' @param variableFeatures Boolean array specifying whether counts are variable features or name of single cell experiment row data columns 
#' with this array. If not provided, will be determined using scry devianceFeatureSelection. If Seurat object is provided which contains veriable
#' features, these will be used for the first clustering, unless variableFeatures is set to 'Find', in which case scry devianceFeatureSelection will be 
#' used. Setting to 'None' will use all features.
#' @param nVarFeatures Number of variable features to conduct pca on. Default is 2000.
#' @param nboots number of boostraps to perform.
#' @param varsToRegress Variables to regress from counts, prior to pca computation. Either a dataframe of variables, with each column being a seperate
#' variable to regress and each row being a cell in counts, or a character vector of names of columns present in a Seurat or SingleCellExperiment 
#' object, which is provided to counts.  
#' @param regressMethod method to regress variables from counts, either "lm" - linear regression of normalised counts, 'glmgampoi' - alpha gamma
#' regression of raw counts using the glmgampoi packaged (requires that glmgampoi is installed on your machine), or 'poisson' - poisson 
#' regression of raw counts.
#' @param skipFirstRegression Either a character vector specifying varsToRegress variables not to regress on the first clustering (e.g. don't regress 
#' the effect of percent of mitochondrial genes or read depth when multiple cell types are likely present, as these factors are likely confounded with 
#' cell type), or a boolean specifying whether to skip regression for all varsToRegress variables during the first clustering. Default: FALSE.
#' @param bootsize Size of bootstrap samples as a proportion of original cell number. Default: 0.9.
#' @param nboots Number of bootstrap samplings to perform, to identify robust clusters. Increasing will increase runtime, decreasing may identify less 
#' robust clusters, or reduce ability to identify true, less distinct clusters. Default: 100.
#' @param clusterFun Which of the igraph community detection algorithms to use for clustering: "leiden" or "louvain".
#' @param resRange Numeric vector of resolutions to cluster over. Default: c(0.01, 0.03, 0.05, 0.07, 0.10, seq.int(0.11, 1.5, length.out=10)).
#' @param kNum number of nearest neighbours for knn graph construction.
#' @param silhouetteThresh Threshold silhouette score, below which to test whether such a score could arise from sampling (poisson) noise.
#' Deafult: 0.4, as we have found sampling noise rarely if ever produces clusters this distinct. Raise to 1 to test all clusterings, at the
#' cost of increased runtime.
#' @param minSize Minimum size of a cluster to try to subcluster. Default: 50 cells. Reduce to 0 to try subcluster everything, at the cost of
#' increased runtime.
#' @param mode How to deal with cluster assignments from different resolutions, either pick the one with highest silhouette score ('robust')
#' or use all when deciding consensus clusters ('granular').
#' @param assay Name of Seurat assay to extract counts from if proving a Seurat object.
#' @param BPPARAM Biocparallel paramters to run in parallel.
#' 
#' @return list containing:
#' 'assignments': character vector of consensus clustering assignments;
#' 'res': the clustering resolution which produced these optimal assignments; and,
#' 'dendrogram': dendrogram showing the relatedness of output cluster assignments, based on the co-clustering distance matrix
#' @export
#' 
#' @importFrom bluster neighborsToSNNGraph pairwiseRand approxSilhouette
#' @importFrom igraph cluster_leiden E
#' @importFrom irlba prcomp_irlba
#' @importFrom BiocParallel bplapply SerialParam
#' @importFrom dbscan kNN
#' @importFrom RcppXPtrUtils  cppXPtr
#' @importFrom parallelDist parDist
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment counts SingleCellExperiment
#' @importFrom MASS fitdistr
#' @importFrom Seurat VariableFeatures
#' @importFrom scry devianceFeatureSelection
#' @importFrom scran getDenoisedPCs modelGeneVarByPoisson
#' @importFrom transformGamPoi shifted_log_transform
#' @importFrom scDesign3 construct_data fit_marginal fit_copula extract_para
#' @importFrom sparseMatrixStats rowMeans2
#' @importFrom stringr str_split
#' @importFrom clustree clustree
#' 
#' @examples
#' library(consensusClustR)
#' 
#' #Example fake counts matrix
#' ncells = 500
#' ngenes = 20000
#' counts = matrix(rpois(ncells*ngenes, 5), ncol=ncells)
#' colnames(counts) = c(1:ncells)
#' rownames(counts) = c(1:ngenes)
#' 
#' #Cluster:
#' results <- consensusClust(counts)
#' 
#' #Cluster with prefined number of PCs (for the first clustering if iterating):
#' results <- consensusClust(counts, pcNum = 7)
#' 
#' #Iteratively cluster until no statistically significant subclusters are found
#' results <- consensusClust(counts, iterate=TRUE)
#' 
#' #Iteratively cluster until no statistically significant subclusters are found, and 10 cores, ensuring reproducible random seed
#' library(BiocParallel)
#' results <- consensusClust(counts, iterate=TRUE, BPPARAM = MulticoreParam(RNGseed = 123, workers=10))
#' 
#' #' #Using 5 PCs, and provinding a SingleCellExperiment experiment object 'data' with raw counts in the 'counts' assay,
#' and normalised features in the "logcounts" assay, iterating until completion:
#' #Requires Seurat package
#' library(Seurat)
#' data = log2(counts + 1) # Scale and logcounts
#' colnames(data) = c(1:ncells)
#' rownames(data) = c(1:ngenes)
#' seurat = CreateSeuratObject(counts = counts, data = data) 
#' results = consensusClust(seurat, pcNum = 5, iterate = T)
#' 
#' #Using 5 PCs, and provinding a SingleCellExperiment experiment object 'data' with raw counts in the 'counts' assay,
#' and normalised features in the "logcounts" assay, iterating until completion:
#' #Requires SingleCellExperiment package from bioconductor
#' library(SingleCellExperiment)
#' colnames(data) = c(1:ncells)
#' rownames(data) = c(1:ngenes)
#' sce = SingleCellExperiment(assays=list(counts=counts, logcounts=data)) 
#' results = consensusClust(sce, pcNum = 5, iterate = T)
#' 
  consensusClust <- function(counts, iterate=FALSE, alpha=0.05, pca=NULL, pcNum="find", pcaMethod = "irlba", scale=T, center=T, 
                                    sizeFactors="deconvolution", variableFeatures=NULL, nVarFeatures=2000, varsToRegress=NULL, 
                                    regressMethod = "lm", skipFirstRegression=F, nboots=100, bootSize=0.9, 
                                  clusterFun="leiden", resRange = c(0.01, 0.03, 0.05, 0.07, 0.10, seq.int(0.11, 1.5, length.out=10)),
                                  kNum=30, silhouetteThresh = 0.4, minSize = 50, assay="RNA", mode = "robust", 
                                  BPPARAM=SerialParam(RNGseed = seed), seed=123, depth=1, ...) {
  
  #Check inputs are correct
  stopifnot("`counts` must be a matrix, sparse matrix of type dgCMatrix, seurat object, or single-cell experiment object." = 
              class(counts)[1] %in% c("Seurat", "SingleCellExperiment", "matrix", "dgCMatrix") )
  if(!is.null(pca)){
    stopifnot("`pca` must be a matrix or sparse matrix of type dgCMatrix." = 
                class(pca)[1] %in% c("matrix", "dgCMatrix") )
  }
  if(all(!is.null(pcNum), !pcNum == "find")){
    stopifnot("`pcNum` must be a positive integer, smaller than cell number." = 
                all((length(pcNum)==1) & (pcNum%%1==0) & (pcNum > 0) & (pcNum < ncol(data))) )
  }
  stopifnot("`pcaMethod` must be one of 'irlba', 'svd', or 'prcomp.'" = 
                  pcaMethod %in% c("irlba", "svd", "prcomp") )
  stopifnot("`nboots` must be a positive integer." = 
              all((nboots%%1==0) & (nboots > 0) & (length(nboots)==1)) )
  stopifnot("`clusterFun` must be a either leiden or louvain." = 
              clusterFun %in% c("leiden", "louvain") )
  stopifnot("`bootSize` must be a single number." = 
              all((length(bootSize)==1) & (is.numeric(bootSize))) )
  stopifnot("`resRange` must be a numeric vector." = 
              is.numeric(resRange) )
  stopifnot("`kNum` must be a positive integer." = 
              all((kNum%%1==0) & (kNum > 0)))
  stopifnot("`mode` must be either 'robust' or 'granular'." = 
              mode %in% c("robust", "granular") )
  
  #Set random seed for sampling
  set.seed(123)
  
  #Get Counts and PCA from Seurat or SingleCellExperiment object
  #First Seurat
  if(class(counts)[1]=="Seurat"){
    if(all(is.null(variableFeatures), !is.null(VariableFeatures(counts)))){
      variableFeatures = VariableFeatures(counts)
    } 
    if(!is.null(variableFeatures)){
      if (variableFeatures[1] =="find"){
        variableFeatures = NULL
      }
    }
    
    #Get vars to regress if metadata column
    if(any(varsToRegress %in% colnames(counts@meta.data))){
      varNames = varsToRegress[varsToRegress %in% colnames(counts@meta.data)]
      varsToRegress = as.data.frame(counts@meta.data[,colnames(counts@meta.data) %in% varsToRegress])
      colnames(varsToRegress) = varNames
      rm(varNames)
    }
    
    #Get PCA if in object
    if(!is.null(counts@reductions$pca)){
      pcaList = counts@reductions$pca@feature.loadings
      pca = counts@reductions$pca@cell.embeddings
    }
    
    #Get scaled or normalised counts if in object
    if(dim(counts[[assay]]$scale.data)[1] > 0){
      normCounts = counts[[assay]]$scale.data
      scaleData=TRUE
    } else if (dim(counts[[assay]]$data)[1] > 0){
      normCounts = counts[[assay]]$data
    }
    
    #Get counts
    counts = counts[[assay]]$counts
    
    #Get variable features if in object
    if(!is.null(variableFeatures)){
    variableFeaturesCounts = rownames(counts) %in% variableFeatures
    if(!is.null(normCounts)){
      variableFeatures = rownames(normCounts) %in% variableFeatures
    }
    }
    #If not SCE
  } else if(class(counts)[1]=="SingleCellExperiment"){
    if(all(is.null(variableFeatures), exists("rowData(counts)[,variableFeatures]"))){
      variableFeatures = rowData(counts)[,variableFeatures]
    } 
    if(!is.null(variableFeatures)){
      if (variableFeatures[1] =="find"){
        variableFeatures = NULL
      }
    }
    
    #Get vars to regress if metadata column
    if(any(varsToRegress %in% colnames(colData(counts)))){
      varNames = varsToRegress[varsToRegress %in% colData(counts)]
      varsToRegress = colData(counts)[,colData(counts) %in% varsToRegress]
      colnames(varsToRegress) = varNames
      rm(varNames)
    }
    
    #Get PCA if in object
    if(!is.null(reducedDims(counts)$PCA)){
      pca = reducedDims(counts)$PCA
    }
    
    #Get normalised counts if in object
    if(!is.null(sce@assays@data$logcounts)){
      normCounts = sce@assays@data$logcounts
    }
    
    #Get counts
    counts = counts(counts)
  }
  
  #Normalise counts if not provided
  if(!exists("normCounts")){
    # if(sizeFactors[1]=="deconvolution"){
    #   sizeFactors = calculateSumFactors(counts, BPPARAM = BPPARAM)
    #   # stabilize size factors to have geometric mean of 1
    #   zeroSFs = any(is.nan(sizeFactors) | sizeFactors <= 0)
    #   sizeFactors[zeroSFs] <- NA
    #   if(zeroSFs){
    #     sizeFactors <- sizeFactors/exp(mean(log(sizeFactors), na.rm=TRUE))
    #     sizeFactors[zeroSFs] <- 0.001
    #   }else{
    #     sizeFactors <- sizeFactors/exp(mean(log(sizeFactors)))
    #   }
    # }
    normCounts = shifted_log_transform(counts, size_factors = sizeFactors)  
  }
  
  #Find variable features if required
  if(is.null(variableFeatures)){
    if(all(exists("variableFeaturesCounts"), !exists("variableFeatures"))){
      variableFeatures = variableFeaturesCounts
    } else if(is.null(variableFeatures)){
      deviance = devianceFeatureSelection(counts)
      variableFeatures = deviance >= -sort(-deviance, partial=nVarFeatures)[nVarFeatures]
      variableFeaturesCounts = variableFeatures
    }
    #Subset normCounts to variable features
    normCounts = normCounts[variableFeatures,]
  }
  
  #Regress out unwanted variables and find pcNum if desired
  if(!is.null(varsToRegress)){
    if(all(isTRUE(skipFirstRegression))){
      skipFirstRegression = colnames(varsToRegress)
    }
  }
  if(!all(colnames(varsToRegress) %in% skipFirstRegression)){
    if(pcNum == "find"){
      #Estmate pcNum with getDenoisedPCs for large clusters
      if(ncol(counts) > 400){
        #Model gene variance before regression
        var.stats <- modelGeneVarByPoisson(normCounts, design=varsToRegress)
        pcNum = ncol(getDenoisedPCs(normCounts, var.stats, subset.row=NULL)$components)
      } else {
        #Otherise use variance explained
        pca = irlba(t(normCounts), 30, scale=if(center){rowSds(normCounts)}else{NULL}, center=if(center){rowMeans2(normCounts)}else{NULL})
        pcNum = max(which(sapply(1:30, \(pcNum) sum(pca[["d"]][1:pcNum])/ sum(pca[["d"]]) ) > 0.25)[1], 5)
        pca = pca$u %*% diag(sqrt(pca$d))
        rownames(pca) = colnames(normCounts)
      }
    }
    normCounts = regressFeatures(normCounts, varsToRegress, regressMethod = regressMethod, BPPARAM = BPPARAM, seed=seed)
  } else if(pcNum == "find"){ #Else just find pcNum if desired
    #Model gene variance
    var.stats <- modelGeneVarByPoisson(normCounts)
    pcNum = ncol(getDenoisedPCs(normCounts, var.stats, subset.row=NULL)$components)
  }
  
  #Compute PCA if not extracted from an input Seurat or SCE object
  if(is.null(pca)){
    pcaList = prcomp_irlba(t(normCounts), pcNum, scale=if(center){rowSds(normCounts)}else{NULL}, center=if(center){rowMeans2(normCounts)}else{NULL})
    pca = pcaList$x
    rownames(pca) = colnames(normCounts)
  } 
  
  #Restrict pca to chosen features if given
  pca = pca[,1:pcNum]
  
  #Get cluster assignments for bootstrapped selections of cells
  clustAssignments = bplapply(1:nboots, \(boot){
    getClustAssignments( pca[sample(rownames(pca), bootSize*length(rownames(pca)), replace = TRUE),],
                         resRange = resRange, kNum=kNum, clusterFun = clusterFun, cellOrder = rownames(pca), mode = mode, seed=seed)
  }, BPPARAM = BPPARAM )
  
  ##Calculate cell-cell distance matrix based on these co-clusterings
  #Get clustering assignments in the right format (dataframe)
  clustAssignments = do.call(cbind, clustAssignments)
  rownames(clustAssignments) = colnames(normCounts)
  
  #Remove normCounts and deviance if present
  suppressWarnings( rm(normCounts) )
  gc(verbose = F)
  
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
  jaccardDist = 1-parDist(clustAssignments, method="custom", func = jaccardDist, threads = BPPARAM$workers) 
  
  ##Adjacency graph from similarity matrix
  knn = kNN(jaccardDist, k = kNum, search = "kdtree")$id
  snnGraph = neighborsToSNNGraph(knn, type = "rank")
  
  ##Cluster adjacency graph at different resolutions
  if(clusterFun=="leiden"){
    finalAssignments = lapply(resRange, \(res){
      cluster_leiden(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                     beta = 0.01, n_iterations = 2)$membership 
    })
  } else if(clusterFun=="louvain"){
    finalAssignments = lapply(resRange, \(res){
      cluster_louvain(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                     beta = 0.01, n_iterations = 2)$membership 
    })
  }
  
  #Assess silhouette score per cluster for every resolution
  clustScores = rank(unlist(lapply(finalAssignments, \(res){
    if(all(length(unique(res))>1, length(unique(res))<length(res)/10)){
      silhouette = mean(approxSilhouette(pca, res)[,3], na.rm = T)
    } else if(length(unique(res))==length(res)){
      -1
      } else {
      0.15
    }
  })), ties.method ="last")
  
  #Take assignments from resolution with highest score
  finalAssignments = finalAssignments[[which(clustScores == max(clustScores))]]
  
  #If more than one output cluster, lets test the robustness and statistical seperation of these clusters
  if(length(unique(finalAssignments)) > 1){
    
    #Compute cluster stability score in the bootstraps
    compare <- function(...) pairwiseRand(..., mode="ratio", adjusted=TRUE)
    collated = lapply(1:ncol(clustAssignments), \(boot) 
                      compare(finalAssignments[clustAssignments[,boot] != -1], clustAssignments[,boot][clustAssignments[,boot] != -1]))
    stabilityMat = apply(simplify2array(collated), 2, rowMeans2, na.rm = TRUE) #Reduce("+", stabilityMat) / length(stabilityMat)
    diag(stabilityMat) = 1
    dimnames(stabilityMat) = list(unique(finalAssignments), unique(finalAssignments))
    stabilityMat[is.na(stabilityMat)] = 1
    
    #Merge clusters with low stablity in the bootstraps
    while(min(stabilityMat) < 0.4){
      finalAssignments[finalAssignments == max(which(stabilityMat == min(stabilityMat), arr.ind = TRUE))] = 
        unique(finalAssignments[finalAssignments == min(which(stabilityMat == min(stabilityMat), arr.ind = TRUE))])      
      stabilityMat[which(stabilityMat == min(stabilityMat), arr.ind = TRUE)] = 1
    }
    
    #Merge clusters so small it's hard to assess their seperation to the nearest cluster
    while(min(table(finalAssignments)) < kNum){
      clustDist = determineHierachy(as.matrix(jaccardDist), finalAssignments, return = "distance")
      diag(clustDist) = 1
      finalAssignments[finalAssignments == names(which.min(table(finalAssignments)))] = 
        colnames(clustDist)[which.min(clustDist[names(which.min(table(finalAssignments))),])]                         
    }
    
    #If still more than one cluster lets test whether it's likely such strong clustering would not appear due to sampling noise
    if(length(unique(finalAssignments)) > 1){
    
    ##If more than one cluster, and a low silhouette score or small clusters test if these assignments are better than noise
    silhouette = mean(approxSilhouette(pca, finalAssignments)[,3], na.rm = T)
    
    if(any(silhouette <= silhouetteThresh, min(table(finalAssignments)<50))){
      
      #Make SCE object of variable features, and model paramters of these features, assuming they come from one distribution, with scDesign3
      sce = SingleCellExperiment(assays=list(counts = counts[variableFeaturesCounts, ]))
      colData(sce) = cbind( colData(sce), finalAssignments, varsToRegress)
      colData(sce)$single = rep(1)
      
      data <- construct_data(sce = sce, assay_use = "counts", celltype = "single", pseudotime = NULL, spatial = NULL, 
                             other_covariates = "logUMI", corr_by = "1")
      marginals <- fit_marginal(data = data, mu_formula = "1", sigma_formula = "1", family_use = "poisson", usebam = T, n_cores = 1)
      copula <- fit_copula(sce = sce, assay_use = "counts", marginal_list = marginals, family_use = "poisson", copula = "gaussian",
          n_cores = 1, input_data = data$dat)
      params <- extract_para(sce = sce, marginal_list = marginals, family_use = "poisson", new_covariate = data$newCovariate,
            data = data$dat, n_cores = 1)
      
      # Generate null distribution for Silhouette score based on simulating a new dataset based on these single population paramters, 
      # clustering, and calculating the silhouette score.
      nullDist = unlist(bplapply(1:20, function(i) {
        generateNullStatistic(sce, params, data, copula, pcNum=pcNum, varsToRegress = varsToRegress,regressMethod = regressMethod, 
                              scale=scale, center = center, kNum=kNum, clusterFun = clusterFun, seed=seed)
      }, BPPARAM = BPPARAM))
      
      # Test the statistical signifcance of the difference in silhouette scores between the NULL and real clusterings - are your clusters
      # significantly better connected than those geneerated if we assume the data is truly from a single population?
      fit <- fitdistr(nullDist,'normal')
      pval <- 1-pnorm(silhouette,mean=fit$estimate[1],sd=fit$estimate[2])
  
    } else {
      #If silhouette score is below the threshold, we can save time by assuming these clusters are real - adjust threshold to avoid this 
      #behaviour
      pval = 0
    }
    #If subclusters are much better than those returned by chance keep subclustering until assignments are no better than chance
    if(all(pval < alpha, iterate)){
      clustersToSubcluster = unique(finalAssignments)[sapply(unique(finalAssignments), \(cluster) sum(finalAssignments == cluster)) > minSize]
      clustersToSubcluster = setNames(clustersToSubcluster, clustersToSubcluster)
      subassignments = bplapply(clustersToSubcluster, function(cluster){
        #Subset vars to regress
        newVarsToRegress = as.data.frame(varsToRegress[finalAssignments == cluster, ])
        colnames(newVarsToRegress) = colnames(varsToRegress)
        consensusClust(counts[,finalAssignments == cluster], pcaMethod=pcaMethod, nboots=nboots, clusterFun=clusterFun,
                              bootSize=bootSize, resRange = resRange, kNum=kNum, mode = mode, variableFeatures=NULL,
                              scale=scale, varsToRegress=newVarsToRegress, regressMethod=regressMethod, depth=depth+1,iterate=T,
                              #sizeFactors = if(length(sizeFactors)>1){sizeFactors[finalAssignments == cluster]}else{sizeFactors},
                              sizeFactors = "deconvolution", BPPARAM=SerialParam(RNGseed = seed), ...)$assignments
      }, BPPARAM=BPPARAM)
      
      #Add subcluster annotations
      for (cluster in clustersToSubcluster[sapply(subassignments, \(subclusters) length(unique(subclusters))) > 1]){
        finalAssignments[finalAssignments == cluster] = paste0(finalAssignments[finalAssignments == cluster], "_", subassignments[[as.character(cluster)]])
      }
      
      #At top level assess cluster relationships
      if(all(depth==1)){
        #Compute cluster dendrogram
        dendrogram = determineHierachy(as.matrix(jaccardDist), finalAssignments)
        
        #Use clustree to visualise parent-child relationships
        clustree = str_split(finalAssignments, "_")
        maxlen <- max(lengths(clustree))
        clustree = lapply(clustree, \(cell) c(if(length(cell) > 1){ c(cell[1], sapply(2:length(cell), \(res) 
                                                                                      paste(cell[1:res], collapse ="_") ) )
                                                                      } else { cell } , rep(NA, maxlen - length(cell))))  
        clustree = as.data.frame(do.call(rbind, clustree))
        #Fill with previous depth if not subclustered
        for(depth in 2:ncol(clustree)){
          clustree[,depth] = coalesce2(clustree[,depth], clustree[,depth-1])
        }
        colnames(clustree) = gsub("V", "Cluster", colnames(clustree))
        if(ncol(clustree) > 1){
          clustree = clustree(clustree, prefix = "Cluster")
        } else {
          clustree = NULL
        }
      } else{
        dendrogram = NULL
        clustree = NULL
      }
      
    } else if (pval >= alpha) {
      #Else if cluster assignments are no better than chance, then don't return assignments as likely overclustered
      finalAssignments = rep(1, length(finalAssignments))
      dendrogram = NULL
    } else {
      #If not iterating just compute dendrogram and return assignments
      dendrogram = determineHierachy(as.matrix(jaccardDist), finalAssignments)
      clustree = NULL
    }
    
    }
  }
  
  #If there were no real clusters found, don't make a dendrogram or clustree
    if(length(unique(finalAssignments)) == 1) {
      finalAssignments = rep(1, length(finalAssignments))
      dendrogram = NULL
      clustree = NULL
    }
  
  return(list(assignments = finalAssignments, clusterDendrogram = dendrogram, clustree=clustree))
  
}

#' Conducts PCA, and clusters cells based on resulting PCs, returning assignments as a dataframe
#' @param scaledCounts  Matrix of scaled counts
#' @param pcNum number of principal components to consider when finding nearest neighbours
#' @param clusterFun Which of the igraph community detection algorithms to use for clustering: "leiden" or "louvain"
#' @param resRange Numeric vector of resolutions to cluster over
#' @param kNum number of nearest neighbours for knn graph construction
#' @param mode How to deal with cluster assignments from different resolutions, either pick the one with highest silhouette score ('robust')
#' or use all when deciding consensus clusters ('granular')
#' @param cellOrder Vector of non-boostrapped cellnames in the order of the non-boostrapped data matrix
#' 
#' @export
#' 
#' @importFrom bluster clusterRows NNGraphParam approxSilhouette
#' 
getClustAssignments <- function(pca, pcNum, clusterFun="leiden", resRange, kNum, mode = "robust", cellOrder, seed, ...) {
  
  #Cluster adjacency matrix, returning assignments named by cell name
  clustAssignments = lapply(resRange, function(res){
    assignments = suppressWarnings( setNames(
      clusterRows(pca, BLUSPARAM=NNGraphParam(k=kNum, cluster.fun=clusterFun, cluster.args=list(resolution=res))),
      rownames(pca) 
    ) )
    
    if(mode == "robust"){
      if(length(unique(assignments))>1){
        clustScore = mean(approxSilhouette(pca, assignments)[,3], na.rm = T)
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
    clustScores = rank(unlist(lapply(clustAssignments, \(res) res[[2]] )), ties.method ="first")
    clustAssignments = clustAssignments[[which(clustScores == max(clustScores))]][[1]]
  } else {
    clustAssignments = do.call(cbind, clustAssignments)
  }
  
  return(clustAssignments)
}

#' Computes dendrogram between clusters
#' @param distanceMatrix  Matrix of cell-cell distances
#' @param assignments clustering assignments
#' 
#' @export
determineHierachy <- function(distanceMatrix, assignments, return = "hclust") {
  
  # Create a distance matrix for the clusters
  clusterDistanceMatrix <- matrix(0, nrow = length(unique(assignments)), 
                                    ncol = length(unique(assignments)),
                                  dimnames = list(unique(assignments), unique(assignments)))
  
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
  
  if(return == "distance"){
    return(clusterDistanceMatrix)
  }
  
  clusterDistanceMatrix = as.dist(clusterDistanceMatrix)
  
  # Perform hierarchical clustering on the cluster distance matrix
  hclustRes <- as.dendrogram(hclust(clusterDistanceMatrix, method = "complete"))
  
  return(hclustRes)
}


#' Generate a null count matrix from a scDesign3 paramters, process and cluster, and return the silhouette score of the clustering
#' @importFrom bluster approxSilhouette
#' @importFrom irlba prcomp_irlba
#' @importFrom scDesign3 simu_new
#' @importFrom transformGamPoi shifted_log_transform
#' @importFrom BiocParallel SerialParam
#' @noRd
#'
generateNullStatistic <- function(sce, my_para, my_data, my_copula,pcNum, scale, center,varsToRegress,regressMethod, seed, ...) {
  
  # # #Simulate one group dataset based on count matrix
  null = simu_new(
      sce = sce,
      mean_mat = my_para$mean_mat,
      sigma_mat = my_para$sigma_mat,
      zero_mat = my_para$zero_mat,
      quantile_mat = NULL,
      copula_list = my_copula$copula_list,
      n_cores = 1,
      fastmvn=F,
      nonzerovar=T,
      family_use = "poisson",
      input_data = my_data$dat,
      new_covariate = my_data$newCovariate,
      important_feature = my_copula$important_feature,
      filtered_gene = my_data$filtered_gene
    )
  null = shifted_log_transform(null, size_factors = "deconvolution")
  null = regressFeatures(null, varsToRegress, regressMethod = regressMethod, BPPARAM = SerialParam(RNGseed = seed), seed=seed)
  
  pcaNull = prcomp_irlba(t(null), pcNum, scale=if(center){rowSds(null)}else{NULL}, center=if(center){rowMeans2(null)}else{NULL})$x
  
  if(all(is.na(pcaNull))){
    return(0)
  }
  
  rownames(pcaNull) = colnames(null)
  
  assignments = getClustAssignments( pcaNull, cellOrder = rownames(pcaNull), resRange = seq.int(0.01, 3, length.out = 30), seed=seed, ...)
  
  sil = if(all(length(unique(assignments)) > 1, min(table(assignments)) > 15, length(unique(assignments)) < length(assignments)/5 )){
    mean(approxSilhouette(pcaNull, assignments)[,3], na.rm = T)
  } else {
    0
  }
    
  return(sil)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Regress effects from a scRNA-seq matrix
#' @importFrom glmGamPoi glm_gp
#' @noRd
#'
regressFeatures = function(normCounts, variablesToRegress,regressMethod, BPPARAM = SerialParam(RNGseed = seed), seed){
  
  #Regress out effect of variablestoRegress
  normCounts = if(regressMethod == "lm"){
    #Make formula to regress genes vs variablesToRegress
    fmla <- as.formula(paste('GENE ~', paste(colnames(variablesToRegress), collapse = '+')))
    
    data = cbind(variablesToRegress, normCounts[1,])
    colnames(data) = c(colnames(variablesToRegress), "GENE")
    qr = lm(fmla, data, qr = TRUE)$qr
    rm(data)
    geneOrder = rownames(normCounts)
    geneChunks = split(rownames(normCounts), ceiling(seq_along(rownames(normCounts))/(nrow(normCounts)/BPPARAM$workers)))
    geneChunks = bplapply(geneChunks, function(genes){
      res = lapply(genes, \(gene) qr.resid(qr = qr, y = normCounts[gene,]) )
      normData = do.call(rbind,  res)
      rownames(normData) = genes
      return(normData)
    } , BPPARAM=BPPARAM)
    normCounts = do.call(rbind, geneChunks)
    normCounts = normCounts[match(geneOrder, rownames(normCounts)),]
    
  } else if (regressMethod == "glmGamPoi"){
    #Make formula to regress genes vs variablesToRegress
    fmla <- as.formula(paste(' ~', paste(colnames(variablesToRegress), collapse = '+')))
    
    residuals(glm_gp(
      as.matrix(normCounts),
      design = fmla,
      col_data = variablesToRegress,
      reference_level = NULL,
      offset = 0,
      size_factors = 1
    ), type = 'pearson')
    
  } else if (regressMethod == "poisson") {
    #Make formula to regress genes vs variablesToRegress
    fmla <- as.formula(paste(' ~', paste(colnames(variablesToRegress), collapse = '+')))
    
    genes = rownames(normCounts)
    #Prepare data for regression
    normCounts <- cbind(variablesToRegress, t(normCounts))
    
    bplapply(genes, \(gene)
             residuals(object = glm(
              formula = fmla,
              family = 'poisson',
              data = normCounts[,c(gene, colnames(variablesToRegress))]),
              type = 'pearson'
              ), BPPARAM=BPPARAM)
             
    normCounts = do.call(cbind, residuals) 
    rownames(normCounts) = genes
  }
  
  return(normCounts)
}

#' Given two vectors of equal length, replaces missing values in the first with values at the same indices in the second
#' @noRd
#' 
coalesce2 <- function(...) {
  Reduce(function(x, y) {
    i <- which(is.na(x))
    x[i] <- y[i]
    x},
    list(...))
}

