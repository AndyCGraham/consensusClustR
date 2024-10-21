#' Performs iterative consensus clustering of a scRNA-seq matrix, to find robust, statistically distinct cell clusters.
#'
#' @param counts gene by cell matrix of scRNA-seq counts, or a Seurat or SingleCellExperiment object from which 
#' these can be extracted.
#' @param normCounts optional gene by cell matrix of normalised scRNA-seq counts, if inputting a matrix to counts.
#' @param iterate Boolean specifying whether to continuously subcluster until no statistically distinct cluster are found. Default: FALSE.
#' @param interactive Boolean specifying whether you want to choose number of PCs to consider at each iteration from an elbow plot 
#' (doesn't work in Rmarkdown). Default: FALSE.
#' @param pcVar If interactive=FALSE (default behaviour), then then the number of PCs to consider for clustering is chosen as the number
#' which contains pcVar proportion of variance of the top 50 PCs, Default: 0.2.
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
#' @param varsToRegress Variables to regress from counts, prior to pca computation. Either a dataframe of variables, with each column being a seperate
#' variable to regress and each row being a cell in counts, or a character vector of names of columns present in a Seurat or SingleCellExperiment 
#' object, which is provided to counts.  
#' @param regressMethod method to regress variables from counts, either "lm" - linear regression of normalised counts, 'glmgampoi' - alpha gamma
#' regression of raw counts using the glmgampoi packaged (requires that glmgampoi is installed on your machine), or 'poisson' - poisson 
#' regression of raw counts.
#' @param skipFirstRegression Either a character vector specifying varsToRegress variables not to regress on the first clustering (e.g. don't regress 
#' the effect of percent of mitochondrial genes or read depth when multiple cell types are likely present, as these factors are likely confounded with 
#' cell type), or a boolean specifying whether to skip regression for all varsToRegress variables during the first clustering. Default: FALSE.
#' @param nboots Number of bootstrap samplings to perform, to identify robust clusters. Increasing will increase runtime, decreasing may identify less 
#' robust clusters, or reduce ability to identify true, less distinct clusters. Default: 100.
#' @param bootsize Size of bootstrap samples as a proportion of original cell number. Default: 0.9.
#' @param minStability Minimum pairwise rand index between two clusters in boostrap clustering assignments (bootstrap stability). 
#' Clusters whose pairwise rand falls below minStability  are merged. Higher minStability increases cluster merging. Calculkated
#' via the pairwiseRand function of tttthe bluster package.
#' @param test_splits_seperately Boolean specifying whether to test each split in each set of clusterings seperately. Default: FALSE. 
#' @param clusterFun Which of the igraph community detection algorithms to use for clustering: "leiden" or "louvain".
#' @param resRange Numeric vector of resolutions to cluster over. Default: c(0.01, 0.03, 0.05, 0.07, 0.10, seq.int(0.11, 1.5, length.out=10)).
#' @param kNum number of nearest neighbours for knn graph construction.
#' @param silhouetteThresh Threshold silhouette score, below which to test whether such a score could arise from sampling (poisson) noise.
#' Deafult: 0.4, as we have found sampling noise rarely if ever produces clusters this distinct. Raise to 1 to test all clusterings, at the
#' cost of increased runtime.
#' @param minSize Minimum size of a cluster to try to subcluster. Default: 50 cells. Reduce to 0 to try subcluster everything, at the cost of
#' increased runtime.
#' @param assay Name of Seurat assay to extract counts from if proving a Seurat object.
#' @param mode How to deal with cluster assignments from different resolutions, either pick the one with highest silhouette score ('robust')
#' or use all when deciding consensus clusters ('granular').
#' @param BPPARAM Biocparallel paramters to run in parallel.
#' @param seed Random seed to ensure reproducibility. Default: 123.
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
#' @importFrom BiocParallel bplapply SerialParam bptry
#' @importFrom dbscan kNN
#' @importFrom RcppXPtrUtils  cppXPtr
#' @importFrom parallelDist parDist
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment counts SingleCellExperiment colData
#' @importFrom MASS fitdistr
#' @importFrom scry devianceFeatureSelection
#' @importFrom scran getDenoisedPCs modelGeneVarByPoisson calculateSumFactors
#' @importFrom transformGamPoi shifted_log_transform
#' @importFrom scDesign3 construct_data fit_marginal fit_copula extract_para
#' @importFrom sparseMatrixStats rowMeans2 rowSds
#' @importFrom stringr str_split
#' @importFrom clustree clustree
#' @importFrom ggplot2 element_blank  ggplot  geom_point  xlab
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
  consensusClust <- function(counts, normCounts=NULL, iterate=FALSE, interactive=FALSE, pcVar = 0.2, alpha=0.05, pca=NULL, pcNum="find", 
                             pcaMethod = "irlba", scale=T, center=T, sizeFactors="deconvolution", variableFeatures=NULL, nVarFeatures=2000, 
                             varsToRegress=NULL, regressMethod = "lm", skipFirstRegression=F, nboots=100, bootSize=0.9, 
                             minStability=0.175, test_splits_seperately=F, clusterFun="leiden", 
                             resRange = c(seq.int(0.01, 0.3, length.out=10), seq.int(0.25, 1.5, length.out=10)),
                             kNum=c(10,15,20), silhouetteThresh = 0.45, minSize = 50, assay="RNA", mode = "robust", 
                             BPPARAM=SerialParam(RNGseed = seed), seed=123, depth=1, ...) {
  
  #Check inputs are correct
  stopifnot("`counts` must be a matrix, sparse matrix of type dgCMatrix, seurat object, or single-cell experiment object." = 
              class(counts)[1] %in% c("Seurat", "SingleCellExperiment", "matrix", "dgCMatrix") )
  stopifnot("`normCounts` must be a matrix, sparse matrix of type dgCMatrix, or NULL." = 
                class(normCounts)[1] %in% c("NULL", "matrix", "dgCMatrix") )
  stopifnot("`iterate` must be TRUE or FALSE." = 
              class(iterate)[1] %in% c("logical") )
  stopifnot("`interactive` must be TRUE or FALSE." = 
              class(interactive)[1] %in% c("logical") )
  stopifnot("`nboots` must be a positive integer." = 
              all((nboots%%1==0) & (nboots > 0) & (length(nboots)==1)) )
  stopifnot("`pcVar` must be a float > 0 < 1." = 
              all((pcVar > 0) & (pcVar < 1) & (length(pcVar)==1)) )
  stopifnot("`alpha` must be a float > 0 < 1." = 
              all((alpha > 0) & (alpha < 1) & (length(alpha)==1)) )
  stopifnot("`pca` must be a matrix, sparse matrix of type dgCMatrix, or NULL." = 
                class(pca)[1] %in% c("NULL", "matrix", "dgCMatrix") )
  if(all(!is.null(pcNum), !pcNum == "find")){
    stopifnot("`pcNum` must be a positive integer, smaller than cell number, NULL, or 'find'." = 
                all((length(pcNum)==1) & (pcNum%%1==0) & (pcNum > 0) & (pcNum < ncol(data))) )
  }
  stopifnot("`pcaMethod` must be one of 'irlba', 'svd', or 'prcomp.'" = 
                  pcaMethod %in% c("irlba", "svd", "prcomp") )
  stopifnot("`scale` must be TRUE or FALSE." = 
              class(scale)[1] %in% c("logical") )
  stopifnot("`center` must be TRUE or FALSE." = 
              class(center)[1] %in% c("logical") )
  stopifnot("`sizeFactors` must be a numeric vector of size factors or 'deconvolution'." = 
              any(class(sizeFactors)[1] %in% c("character") & sizeFactors[1] == "deconvolution", class(sizeFactors) == "numeric")  )
  stopifnot("`variableFeatures` must be a character vector of names of variable features or a boolean of length nrow(counts) specifying whether 
            each feature is a variable feature ." = 
              any(class(variableFeatures) == "character", class(variableFeatures) == "logical")  )
  stopifnot("`nVarFeatures` must be a positive integer." = 
              all((nVarFeatures%%1==0) & (nVarFeatures > 0) & (length(nVarFeatures)==1)) )
  stopifnot("`regressMethod` must be 'lm' or 'glmGamPoi'." = 
              regressMethod %in% c("glmGamPoi", "lm") )
  stopifnot("`skipFirstRegression` must be TRUE or FALSE." = 
              class(center)[1] %in% c("logical") )
  stopifnot("`nboots` must be a positive integer." = 
              all((nboots%%1==0) & (nboots > 0) & (length(nboots)==1)) )
  stopifnot("`bootSize` must be a float > 0 < 1." = 
              all((bootSize > 0) & (bootSize < 1) & (length(bootSize)==1)) )
  stopifnot("`minStability` must be a float > 0 < 1." = 
              all((minStability > 0) & (minStability < 1) & (length(minStability)==1)) )
  stopifnot("`test_splits_seperately` must be a float > 0 < 1." = 
              all((test_splits_seperately > 0) & (test_splits_seperately < 1) & (length(test_splits_seperately)==1)) )
  stopifnot("`clusterFun` must be a either leiden or louvain." = 
              clusterFun %in% c("leiden", "louvain") )
  stopifnot("`resRange` must be a numeric vector." = 
              is.numeric(resRange) )
  stopifnot("`kNum` must be a positive integer or vector of integers." = 
              any(is.numeric(kNum), all((kNum%%1==0) & (kNum > 0))))
  stopifnot("`silhouetteThresh` must be a float > 0 < 1." = 
              all((silhouetteThresh > 0) & (silhouetteThresh < 1) & (length(silhouetteThresh)==1)) )
  stopifnot("`minSize` must be a positive integer." = 
              all((minSize%%1==0) & (minSize > 0) & (length(minSize)==1)) )
  stopifnot("`assay` must be a single character string." = 
              all(class(assay) == "character", length(assay) == 1)  )
  stopifnot("`mode` must be either 'robust' or 'granular'." = 
              mode %in% c("robust", "granular") )
  stopifnot("`seed` must be a positive integer." = 
              all((seed%%1==0) & (seed > 0) & (length(seed)==1)) )
  
  #Set random seed for sampling
  set.seed(123)
  
  #Get Counts and PCA from Seurat or SingleCellExperiment object
  #First Seurat
  if(class(counts)[1]=="Seurat"){
    if(all(is.null(variableFeatures), !is.null(Seurat::VariableFeatures(counts)))){
      variableFeatures = Seurat::VariableFeatures(counts)
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
      varNames = varsToRegress[varsToRegress %in% colnames(colData(counts))]
      varsToRegress = as.data.frame(colData(counts)[,colnames(colData(counts)) %in% varsToRegress])
      colnames(varsToRegress) = varNames
      rm(varNames)
    }
    
    #Get PCA if in object
    if(!is.null(reducedDims(counts)$PCA)){
      pca = reducedDims(counts)$PCA
    }
    
    #Get normalised counts if in object
    if(!is.null(counts@assays@data$logcounts)){
      normCounts = counts@assays@data$logcounts
    }
    
    #Get counts
    counts = counts(counts)
  }
  
  #Normalise counts if not provided
  if(sizeFactors[1]=="deconvolution"){
    sizeFactors = calculateSumFactors(counts, BPPARAM = BPPARAM)
      # stabilize size factors to have geometric mean of 1
      zeroSFs = any(is.nan(sizeFactors) | sizeFactors <= 0)
      sizeFactors[zeroSFs] <- NA
      if(zeroSFs){
        sizeFactors <- sizeFactors/exp(mean(log(sizeFactors), na.rm=TRUE))
        sizeFactors[zeroSFs] <- 0.001
        }else{
        sizeFactors <- sizeFactors/exp(mean(log(sizeFactors)))
        }
      }
  if(is.null(normCounts)){
    normCounts = shifted_log_transform(counts, size_factors = sizeFactors, pseudo_count = 1)  
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
  }
    
  if(!exists("scaleData")){
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
    
    if(!exists("scaleData")){
      #Regress out unwanted effects
      normCounts = regressFeatures(normCounts, varsToRegress, regressMethod = regressMethod, BPPARAM = BPPARAM, seed=seed)
    } else {
      rm(scaleData)
    }
    
    if(pcNum == "getDenoisedPCs"){
      #Estmate pcNum with getDenoisedPCs for large clusters
      if(ncol(counts) > 400){
        #Model gene variance before regression
        var.stats <- modelGeneVarByPoisson(counts, design=varsToRegress, size.factors=sizeFactors, subset.row=variableFeaturesCounts)
        pcNum = ncol(getDenoisedPCs(normCounts, var.stats, subset.row=NULL)$components)
      } 
    }
  } else if(pcNum == "getDenoisedPCs"){ #Else just find pcNum if desired
    #Model gene variance
    if (ncol(counts) > 400) {
      var.stats <- modelGeneVarByPoisson(counts, size.factors = sizeFactors, subset.row=variableFeaturesCounts)
      pcNum = ncol(getDenoisedPCs(normCounts, var.stats, subset.row=NULL)$components)
    }
  }
  
  #If finding pcNum doesn't work well or cluster is very small, find the elbow of variance explained
  if (any(pcNum == "find", pcNum > 30)) {
    pca = prcomp_irlba(t(as.matrix(normCounts)), 50, scale=if(center){rowSds(normCounts)}else{NULL}, center=if(center){rowMeans2(normCounts)}else{NULL})
    if(!all(is.na(pca))){
      if(interactive){
        elbowPlot = data.frame(stdev=pca[["sdev"]], pc=c(1:length(pca[["sdev"]])))
        elbowPlot = ggplot(elbowPlot, aes(x=pc, y=stdev)) + geom_point() + 
          theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank(), panel.background = element_blank()) +
          xlab("Principal Component")
        print(elbowPlot)
        pcNum = readline("Choose PC Number From Elbow Plot: ")
        pcNum= as.integer(pcNum)
      } else {
        #angles = sapply(4:(length(pca$sdev)-2), \(pc) atan( (mean(pca[["sdev"]][(pc-1):pc]) - mean(pca[["sdev"]][(pc+1):(pc+2)])) )  )
        #tangents = sapply(2:length(angles), \(pc) angles[pc] / angles[pc+1]  )
        #pcNum = which.max(tangents) + 4
        #pcNum = max(which(sapply(1:length(angles), \(pc) angles[pc] / mean(angles[1:3])  ) 
        #                  > 0.1)) + 2
        #pcNum = max(pcNum, 5)
        pcNum = max(which(sapply(1:50, \(pcNum) sum(pca[["sdev"]][1:pcNum])/ sum(pca[["sdev"]]) ) > pcVar)[1], 5)
        
        #pcNum = which.max(sapply(5:45, \(pcNum) ( mean(sapply((pcNum-2):pcNum, \(pc) pca[["sdev"]][pc-1] - pca[["sdev"]][pc])) - 
        #                                            mean(sapply(pcNum:(pcNum+2), \(pc) pca[["sdev"]][pc] - pca[["sdev"]][pc-1])) ) /
        #                           pca[["sdev"]][pcNum])) + 4
      }
      pca = pca$x
      rownames(pca) = colnames(normCounts)
    }
  } 
  
  if(is.null(pca)){
    pca = tryCatch(
      { prcomp_irlba(t(as.matrix(normCounts)), pcNum, scale=if(center){rowSds(normCounts)}else{NULL}, 
                     center=if(center){rowMeans2(normCounts)}else{NULL})$x
      }, error = function(e) {
        NA
      }, warning=function(e){
        NA
      } )
  }
  if(all(is.na(pca))){
    return(list(assignments = rep(1, ncol(counts))))
  }
  #Restrict pca to chosen features if given
  rownames(pca) = colnames(normCounts)
  pca = pca[,1:pcNum]
  
  #Remove normCounts if present
  suppressWarnings( rm(normCounts) )
  
  #Perform bootstrapped clustering if nboots>1
  if(nboots>1){
    
    #Get cluster assignments for bootstrapped selections of cells
    clustAssignments = bplapply(1:nboots, \(boot){
      tryCatch(
        {
          getClustAssignments( pca[sample(rownames(pca), bootSize*nrow(pca), replace = TRUE),],
                               resRange = resRange, kNum=kNum, clusterFun = clusterFun, cellOrder = rownames(pca), mode = mode, seed=seed)
        },
        error = function(e) {
          rep(1, length(rownames(pca)))
        } )
    }, BPPARAM = BPPARAM )
    
    ##Calculate cell-cell distance matrix based on these co-clusterings
    #Get clustering assignments in the right format (dataframe)
    clustAssignments = do.call(cbind, clustAssignments)
    rownames(clustAssignments) = colnames(counts)
    
    #Alter NAs to -1, which will be ignored by the distance function
    clustAssignments[is.na(clustAssignments)] = -1
    
    #Function to compute jaccard similarity between non-binary sets (ignoring -1s), implemented in C++
    jaccardDist = cppXPtr("double customDist(const arma::mat &A, const arma::mat &B) { 
                                float overlap, U;
                                double jaccard;
                                overlap = arma::accu((A == B) && (A != -1));
                                U = arma::accu((A != -1) && (B != -1));
                                jaccard = overlap / U;
                                return jaccard; }",
                                          depends = c("RcppArmadillo"))
    
    #Calculate distance matrix
    jaccardDist = 1-parDist(clustAssignments, method="custom", func = jaccardDist, threads = BPPARAM$workers) 
    
    finalAssignments = unlist(lapply(kNum, function(k){
      ##Adjacency graph from similarity matrix
      knn = kNN(jaccardDist, k = k, search = "kdtree")$id
      snnGraph = neighborsToSNNGraph(knn, type = "rank")
      
      ##Cluster adjacency graph at different resolutions
      if(clusterFun=="leiden"){
        finalAssignments = bplapply(resRange, \(res){
          cluster_leiden(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                         beta = 0.01, n_iterations = 2)$membership 
        }, BPPARAM = BPPARAM )
      } else if(clusterFun=="louvain"){
        finalAssignments = bplapply(resRange, \(res){
          cluster_louvain(snnGraph, objective_function = "modularity", resolution_parameter = res, 
                          beta = 0.01, n_iterations = 2)$membership 
        }, BPPARAM = BPPARAM )
      }
      return(finalAssignments)
    }), recursive = F)
    
    
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
    
    #If more than one output cluster, lets test the robustness of these clusters in the bootstraps
    if(length(unique(finalAssignments)) > 1){
      
      #Merge clusters so small it's hard to assess their seperation to the nearest cluster
      while(min(table(finalAssignments)) < max(kNum[1], 20)){
        clustDist = determineHierachy(as.matrix(jaccardDist), finalAssignments, return = "distance")
        diag(clustDist) = 1
        finalAssignments[finalAssignments == names(which.min(table(finalAssignments)))] = 
          colnames(clustDist)[which.min(clustDist[names(which.min(table(finalAssignments))),])]                         
      }
      
      #Compute cluster stability score in the bootstraps
      compare <- function(...) pairwiseRand(..., mode="ratio", adjusted=TRUE)
      
      #Merge clusters with low stablity in the bootstraps
      collated = lapply(1:ncol(clustAssignments), \(boot) 
                        compare(finalAssignments[clustAssignments[,boot] != -1], clustAssignments[,boot][clustAssignments[,boot] != -1]))
      stabilityMat = tryCatch(
        {
          apply(simplify2array(collated), 2, rowMeans2, na.rm = TRUE) 
        },
        error = function(e) {
          NULL
        } ) 
      if(is.null(stabilityMat)){
        finalAssignments = rep(1, length(finalAssignments))
      } else {
        diag(stabilityMat) = 1
        dimnames(stabilityMat) = list(unique(finalAssignments), unique(finalAssignments))
        stabilityMat[is.na(stabilityMat)] = 1
        
        while(min(stabilityMat) < minStability){
          clustersToMerge = as.numeric(which(stabilityMat == min(stabilityMat), arr.ind = TRUE))
          finalAssignments[finalAssignments == clustersToMerge[2]] = clustersToMerge[1]  
          clustAssignments[clustAssignments == clustersToMerge[2]] = clustersToMerge[1] 
          stabilityMat[clustersToMerge[1], clustersToMerge[2]] = 1
          stabilityMat[clustersToMerge[2], clustersToMerge[1]] = 1
        }
      }
    }
  } else {
      #If not bootstrapping just iterate over kNums and resRange and return highest silhouette score
      finalAssignments = as.character(getClustAssignments( pca, resRange = resRange, kNum=kNum, 
                                              clusterFun = clusterFun, cellOrder = rownames(pca), 
                                              mode = "robust", seed=seed ))
      
      #Merge clusters so small it's hard to assess their seperation to the nearest cluster
      while(min(table(finalAssignments)) < max(kNum[1], 30)){
        clustDist = determineHierachy(as.matrix(dist(pca)), finalAssignments, return = "distance")
        diag(clustDist) = 1
        finalAssignments[finalAssignments == names(which.min(table(finalAssignments)))] = 
          colnames(clustDist)[which.min(clustDist[names(which.min(table(finalAssignments))),])]                         
      }
    }
    
  
  #If still more than one cluster lets test whether it's likely such strong clustering would not appear due to sampling noise
  if(all(length(unique(finalAssignments)) > 1)){
  
    #Clustering quality score
    silhouette = mean(approxSilhouette(pca, finalAssignments)[,3], na.rm = T)
    
    #If below the threshold or small clusters test if this is better than chance
    if(any(silhouette <= silhouetteThresh, min(table(finalAssignments)<50))){
      
      dend = determineHierachy(as.matrix(dist(pca)), finalAssignments)
      
      #Make SCE object of variable features, and model paramters of these features, assuming they come from one distribution, with scDesign3
      sce = SingleCellExperiment(assays=list(counts = counts[variableFeaturesCounts, ]))
      if(!is.null(varsToRegress)){
        colData(sce) = cbind( colData(sce), varsToRegress)
      }
      
      colData(sce)$single = rep(1)
      
      #If more than one split, test each split in order
      #Iterate/Walk the dendrogram
      finalAssignments = testSplits(sce, pca = pca, dend, kNum, alpha, finalAssignments, varsToRegress, BPPARAM=BPPARAM,
                                regressMethod = regressMethod, scale=scale, center = center, test_sep=test_splits_seperately,
                                clusterFun = clusterFun, seed=seed, pcNum = pcNum, silhouette=silhouette, silhouetteThresh=silhouetteThresh)
  
    } 
    
    #If subclusters are much better than those returned by chance keep subclustering until assignments are no better than chance
    if(all(length(unique(finalAssignments)) > 1, iterate, any(sapply(unique(finalAssignments), \(cluster) sum(finalAssignments == cluster)) > minSize))){
      clustersToSubcluster = unique(finalAssignments)[sapply(unique(finalAssignments), \(cluster) sum(finalAssignments == cluster)) > minSize]
      clustersToSubcluster = setNames(clustersToSubcluster, clustersToSubcluster)
      
      subassignments = lapply(clustersToSubcluster, function(cluster){
        if(!is.null(varsToRegress)){
          #Subset vars to regress
          newVarsToRegress = as.data.frame(varsToRegress[finalAssignments == cluster, ])
          colnames(newVarsToRegress) = colnames(varsToRegress)
        } else {
          newVarsToRegress = NULL
        }
        #Reduce some parameters for small clusters
        # if(all(sum(finalAssignments == cluster) < 500, ncol(counts) >= 500)){
        #   kNum = kNum[kNum <= 15]
        #   if(length(kNum)<1){
        #     kNum=c(10,15)
        #   }
        #   pcVar = min(0.2, pcVar)
        # }
        consensusClust(counts[,finalAssignments == cluster], pcaMethod=pcaMethod, 
                       nboots=nboots, clusterFun=clusterFun, bootSize=bootSize, resRange = resRange, kNum=kNum, mode = mode, 
                       variableFeatures=NULL, scale=scale, varsToRegress=newVarsToRegress, regressMethod=regressMethod, depth=depth+1,
                       iterate=T, interactive=interactive, BPPARAM=BPPARAM, pcVar=pcVar,
                       minStability=minStability, test_splits_seperately=test_splits_seperately, ...)$assignments
      })
      
      # Replace assignments with subassignments if required
      if(exists("subassignments")){
        # Replace errors (from pca not being able to be run in tiny clusters etc.) with lack of clustering
        subassignments[sapply(subassignments, \(subcluster) length(subcluster) == 2)] = "1"
        
        # Add subcluster annotations
        for (cluster in clustersToSubcluster[sapply(subassignments, \(subclusters) length(unique(subclusters))) > 1]){
          finalAssignments[finalAssignments == cluster] = paste0(finalAssignments[finalAssignments == cluster], "_", subassignments[[as.character(cluster)]])
        }
      }
      
      #At top level assess cluster relationships
      if(all(depth==1)){
        
        if(nboots>1){
          #Compute cluster dendrogram
          dendrogram = determineHierachy(as.matrix(jaccardDist), finalAssignments)
        } else {
          dendrogram = determineHierachy(as.matrix(dist(pca)), finalAssignments)
        }
        
        #Use clustree to visualise parent-child relationships
        clustree = str_split(finalAssignments, "_")
        maxlen <- max(lengths(clustree))
        clustree = lapply(clustree, \(cell) c(if(length(cell) > 1){ c(cell[1], sapply(2:length(cell), \(res) 
                                                                                      paste(cell[1:res], collapse ="_") ) )
                                                                      } else { cell } , rep(NA, maxlen - length(cell))))  
        clustree = as.data.frame(do.call(rbind, clustree))
        if(ncol(clustree) > 1){
          #Fill with previous depth if not subclustered
          for(depth in 2:ncol(clustree)){
            clustree[,depth] = coalesce2(clustree[,depth], clustree[,depth-1])
          }
          colnames(clustree) = gsub("V", "Cluster", colnames(clustree))
          clustree = clustree(clustree, prefix = "Cluster")
        } else {
          clustree = NULL
        }
      } else{
        dendrogram = NULL
        clustree = NULL
      }
      
    } else if (length(unique(finalAssignments)) == 1) {
      message("Failed Test")
      #Else if cluster assignments are no better than chance, then don't return assignments as likely overclustered
      dendrogram = NULL
      clustree = NULL
    } else {
      #If not iterating just compute dendrogram and return assignments
      if(nboots>1){
        #Compute cluster dendrogram
        dendrogram = determineHierachy(as.matrix(jaccardDist), finalAssignments)
      } else {
        dendrogram = determineHierachy(as.matrix(dist(pca)), finalAssignments)
      }
      clustree = NULL
    }
  
  } else {#{If no clusters found, return this and don't make a dendrogram or clustree
    return(list(assignments = finalAssignments))
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
#' @importFrom bluster clusterRows SNNGraphParam approxSilhouette
#' 
getClustAssignments <- function(pca, clusterFun="leiden", resRange, kNum, mode = "robust", cellOrder, seed, minSize=0, ...) {
  
  #Cluster adjacency matrix, returning assignments named by cell name
  clustAssignments = unlist(lapply(kNum, function(k){
    lapply(resRange, function(res){
      assignments = suppressWarnings( setNames(
        clusterRows(pca, BLUSPARAM=SNNGraphParam(k=k, cluster.fun=clusterFun, 
                                                 type="number",
                                                cluster.args=list(resolution=res))),
        rownames(pca) 
      ) )
      
      if(mode == "robust"){
        if(all(length(unique(assignments))>1, min(table(assignments)) > minSize)){
          clustScore = mean(approxSilhouette(pca, assignments)[,3], na.rm = T)
        } else if(min(table(assignments)) > minSize) {
          clustScore = 0
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
  } ), recursive=FALSE)
    
  
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
determineHierachy <- function(distanceMatrix, assignments, return = "dendrogram") {
  
  # Create a distance matrix for the clusters
  clusterDistanceMatrix <- matrix(0, nrow = length(unique(assignments)), 
                                    ncol = length(unique(assignments)),
                                  dimnames = list(unique(assignments), unique(assignments)))
  
  # Fill the cluster distance matrix with the distances between cluster centroids
  for (clust1 in unique(assignments)) {
    for (clust2 in unique(assignments)[!unique(assignments) == clust1]) {
      if(clusterDistanceMatrix[as.character(clust1), as.character(clust2)] == 0){
        clust1Samples = which(assignments == clust1)
        clust2Samples = which(assignments == clust2)
        dist = mean(distanceMatrix[clust1Samples, clust2Samples], na.rm = T)
        clusterDistanceMatrix[as.character(clust1), as.character(clust2)] = dist
        clusterDistanceMatrix[as.character(clust2), as.character(clust1)] = dist
      }
    }
  }
  
  if(return == "distance"){
    return(clusterDistanceMatrix)
  }
  
  clusterDistanceMatrix = as.dist(clusterDistanceMatrix)
  
  # Perform hierarchical clustering on the cluster distance matrix
  hclustRes <- hclust(clusterDistanceMatrix, method = "complete")
  
  if(return == "hclust"){
    return(hclustRes)
  }
  
  dend = as.dendrogram(hclustRes)
  
  return(dend)
}


#' Generate a null count matrix from a scDesign3 paramters, process and cluster, and return the silhouette score of the clustering
#' @param sce  Single Cell experiment object
#' @param my_para scDesign paramter object
#' @param my_para scDesign paramter object
#' @param my_data scDesign data object
#' @param my_copula scDesign copula object
#' @param pcNum number of PC's
#' @param scale Logical - whether to scale null data prior to PCA
#' @param kNum Clustering kNum
#' @param center Logical - whether to center null data prior to PCA
#' @param varsToRegress Dataframe containing variables to regress in columns
#' @param seed Random seed number. Default: 123.
#' 
#' @importFrom bluster approxSilhouette
#' @importFrom irlba prcomp_irlba
#' @importFrom scDesign3 simu_new
#' @importFrom transformGamPoi shifted_log_transform
#' @importFrom BiocParallel SerialParam
#' @importFrom SingleCellExperiment colData
#' @noRd
#'
generateNullStatistic <- function(sce, my_para, my_data, my_copula, pcNum, scale, 
                                  kNum, center,varsToRegress, seed=123, ...) {
  
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
      family_use = "nb",
      input_data = my_data$dat,
      new_covariate = my_data$newCovariate,
      important_feature = my_copula$important_feature,
      filtered_gene = my_data$filtered_gene
    )
  null = shifted_log_transform(null, size_factors = "deconvolution", pseudo_count = 1)
  if(!is.null(varsToRegress)){
    varsToRegress[1:nrow(varsToRegress), 1:ncol(varsToRegress)] = as.data.frame(colData(sce)[,colnames(colData(sce)) %in% colnames(varsToRegress)])[1:nrow(varsToRegress), 1:ncol(varsToRegress)]
    null = regressFeatures(null, regressMethod="lm", BPPARAM = SerialParam(RNGseed = seed),
                           seed=seed,variablesToRegress = varsToRegress)
  } else {
    null = as.matrix(null)
  }
  
  pcaNull = tryCatch(
    {
      prcomp_irlba(t(null), pcNum, scale=if(center){rowSds(null)}else{NULL}, center=if(center){rowMeans2(null)}else{NULL})$x
    },
    error = function(e) {
      NA
    } )
  
  if(all(is.na(pcaNull))){
    return(0)
  }
  
  rownames(pcaNull) = colnames(null)
  
  assignments = as.character(getClustAssignments( pcaNull, cellOrder = rownames(pcaNull),kNum=kNum,
                                     seed=seed, resRange = c(seq(0.01, 0.3, 0.03), seq(0.3, 2, 0.2)),
                                     minSize = 5, ...))
  
  #Return a score of 0 if only 1 cluster
  if(length(unique(assignments))< 2){
    return(0)
  }
  
  sil = mean(approxSilhouette(pcaNull, assignments)[,3], na.rm = T)
    
  return(sil)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Regress effects from a scRNA-seq matrix
#' @noRd
#'
regressFeatures = function(normCounts, variablesToRegress,regressMethod, BPPARAM, seed){
  
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
      res = bplapply(genes, \(gene) qr.resid(qr = qr, y = normCounts[gene,]) )
      normData = do.call(rbind,  res)
      rownames(normData) = genes
      return(normData)
    } , BPPARAM=BPPARAM)
    normCounts = do.call(rbind, geneChunks)
    normCounts = normCounts[match(geneOrder, rownames(normCounts)),]
    
  } else if (regressMethod == "glmGamPoi"){
    #Make formula to regress genes vs variablesToRegress
    fmla <- as.formula(paste(' ~', paste(colnames(variablesToRegress), collapse = '+')))
    
    residuals(glmGamPoi::glm_gp(
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


#' Generate a null count matrix from a scDesign3 paramters, process and cluster, and return the silhouette score of the clustering
#' @importFrom bluster approxSilhouette
#' @importFrom scDesign3 construct_data fit_marginal  fit_copula  extract_para
#' @importFrom dplyr case_match
#' @importFrom dendextend cutree
#' @export
#'
#' 
testSplits <- function(sce, pca, dend, kNum, alpha, finalAssignments, varsToRegress, 
                       BPPARAM,silhouette, silhouetteThresh,test_sep=F, resRange=seq(0.1, 3.4, 0.15), ...) {
  
  if(test_sep){
    cccm <- cophenetic(dend)
    sps <- sort(unique(cccm), decreasing = T)
    membership <- cutree(dend, h=floor(sps[max(1, which(sps > max(sps)*(8.5/10)))])) #Cut the tree at first splitting point
    
    #Rename assignments to reflect only the current split
    assignments = membership[as.character(finalAssignments)] 
    
    silhouette = mean(approxSilhouette(pca, assignments)[,3], na.rm = T)
  } else {
    assignments = finalAssignments
  }
    
  if(silhouette <= silhouetteThresh){
    
    data <- construct_data(sce = sce, assay_use = "counts", celltype = "single", pseudotime = NULL, 
                           spatial = NULL, other_covariates = colnames(varsToRegress), 
                           corr_by = "1", parallelization="bpmapply", 
                           BPPARAM = BPPARAM)
    marginals <- fit_marginal(data = data, mu_formula = "1", sigma_formula = "1", family_use = "nb", 
                              usebam = T, parallelization="bpmapply", 
                              BPPARAM = BPPARAM, n_cores=BPPARAM$workers)
    copula <- fit_copula(sce = sce, assay_use = "counts", marginal_list = marginals, family_use = "nb", copula = "gaussian",
                         input_data = data$dat, parallelization="bpmapply", 
                         BPPARAM = BPPARAM, n_cores=BPPARAM$workers)
    params <- extract_para(sce = sce, marginal_list = marginals, family_use = "nb", 
                           new_covariate = data$newCovariate, data = data$dat, parallelization="bpmapply", 
                           BPPARAM = BPPARAM, n_cores=BPPARAM$workers)
    
    # Generate null distribution for Silhouette score based on simulating a new dataset based on these single population paramters, 
    # clustering, and calculating the silhouette score.
    # nullDist = tryCatch( {
    #                       unlist(bplapply(1:20, function(i) {
    #                         generateNullStatistic(sce=sce, params, data, copula, kNum=kNum,
    #                                               varsToRegress = varsToRegress, ...)
    #                       }, BPPARAM = BPPARAM)) 
    # }, error = function(e) {
    #   c(1,1,1,1,1)
    # } )
    nullDist =  unlist(bplapply(1:20, function(i) {
                              generateNullStatistic(sce=sce, params, data, copula, kNum=kNum,
                                                    varsToRegress = varsToRegress, ...)
                            }, BPPARAM = BPPARAM))
    # Test the statistical signifcance of the difference in silhouette scores between the NULL and real clusterings - are your clusters
    # significantly better connected than those geneerated if we assume the data is truly from a single population?
    fit <- fitdistr(nullDist,'normal')
    pval <- 1-pnorm(silhouette,mean=fit$estimate[1],sd=fit$estimate[2])
    
    #Test more null samples if significance not clear
    if(all(pval < 0.1, pval >= 0.05)){
      BPPARAM$RNGseed = BPPARAM$RNGseed + 1 #Change RNGseed to get different results
      nullDist2 =  unlist(bplapply(1:20, function(i) {
        generateNullStatistic(sce=sce, params, data, copula, kNum=kNum,
                              varsToRegress = varsToRegress, ...)
      }, BPPARAM = BPPARAM))
      nullDist = c(nullDist, nullDist2)
      fit <- fitdistr(nullDist,'normal')
      pval <- 1-pnorm(silhouette,mean=fit$estimate[1],sd=fit$estimate[2])
    }
    
    #And Again
    if(all(pval < 0.075, pval >= 0.05)){
      BPPARAM$RNGseed = BPPARAM$RNGseed + 1 #Change RNGseed to get different results
      nullDist2 =  unlist(bplapply(1:20, function(i) {
        generateNullStatistic(sce=sce, params, data, copula, kNum=kNum,
                              varsToRegress = varsToRegress, ...)
      }, BPPARAM = BPPARAM))
      nullDist = c(nullDist, nullDist2)
      fit <- fitdistr(nullDist,'normal')
      pval <- 1-pnorm(silhouette,mean=fit$estimate[1],sd=fit$estimate[2])
    }
    
    #If failed test then merge split cluster(s) to closest cluster and test next split if there's one
    if(pval >= alpha){
      if(!test_sep){ #If we're testing multiple clusterings and it fails, reject all
        return(rep(1, length(finalAssignments)))
      }
      while(all(pval >= alpha, length(unique(finalAssignments))>1)){
        names(assignments) = finalAssignments
        clusts_to_merge = sapply(unique(assignments), \(cluster) 
                                  names(membership[membership==cluster])[which.max(table(finalAssignments[finalAssignments %in% names(membership[membership==cluster])]))]
                                 )
        #Join clusters to merge into one cluster if there is multiple
        finalAssignments[finalAssignments %in% clusts_to_merge] = clusts_to_merge[1]
        
        #Keep testing lower splits if present
        if(length(unique(finalAssignments))>1){
          #Make new dend withput merged clusters
          dend = determineHierachy(as.matrix(dist(pca)), finalAssignments)
          cccm <- cophenetic(dend)
          sps <- sort(unique(cccm), decreasing = T)
          membership <- cutree(dend, h=floor(sps[max(1, which(sps > max(sps)*(8.5/10)))])) #Cut the tree at first splitting point
          
          #Rename assignments to reflect only the current split
          assignments = membership[as.character(finalAssignments)] 
          
          silhouette = mean(approxSilhouette(pca, assignments)[,3], na.rm = T)
          
          pval = 1-pnorm(silhouette,mean=fit$estimate[1],sd=fit$estimate[2])
        } 
      } 
      
      #If this merges all clusters, then return failed test
      if(length(unique(finalAssignments))==1){
        return(finalAssignments)
      }
    } 
  }
  
  if(test_sep){
    #Otherwise test the next split(s) (ignoring splits containing merged clusters)
    dend <- cut(dend, h=floor(sps[max(1, which(sps > max(sps)*(8.5/10)))]))$lower
    if(any(sapply(dend, \(split) is.list(split)))){
      lowerAssignments = lapply(dend, \(lowerSplit){
        if(all(is.list(lowerSplit), labels(lowerSplit) %in% finalAssignments)){
          
          if(!is.null(varsToRegress)){
            newVarsToRegress = as.data.frame(varsToRegress[finalAssignments %in% labels(lowerSplit),])
            colnames(newVarsToRegress) = colnames(varsToRegress)
          } else {
            newVarsToRegress = NULL
          }
        
          testSplits(sce[,finalAssignments %in% labels(lowerSplit)], pca[finalAssignments %in% labels(lowerSplit),],
                     lowerSplit, kNum, alpha, finalAssignments[finalAssignments %in% labels(lowerSplit)],
                     varsToRegress = newVarsToRegress, BPPARAM = BPPARAM, test_sep=test_sep,
                     silhouette=silhouette, silhouetteThresh=silhouetteThresh, ...)
        } else {
          NA
        }
      })
      
      #Replace clusters in finalAssignments with results from testing any lower splits
      for(split in 1:length(lowerAssignments)){
        if(!is.na(lowerAssignments[[split]][[1]])){
            finalAssignments[finalAssignments %in% labels(dend[[split]])] = lowerAssignments[[split]]
          }
      }
    }
    
  } 
  #Return whichever clusters survived testing
  return(finalAssignments)
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



