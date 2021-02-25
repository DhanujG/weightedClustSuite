
#library(reticulate)
#library(stats)
#library(dbscan)
#library(FNN)
#library(rtsne)
#library(rcpp)
#library(graphics)
#Library(cluster)
#Library(factoextra)
#Library(dbscan)


#Rcpp::sourceCpp('~/src/distance_optim_calc_func.cpp')

#reticulate::source_python("DBCV.py")
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- ClustObj(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#'
#' irisClust <- findDensityPeakClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#'
#' @useDynLib weightedClustSuite
#' @importFrom Rcpp sourceCpp
#' 
'_PACKAGE'

#' @export
localDensity <- function(weights, distance, dc, gaussian = FALSE) {
  # These implementations are faster by virtue of being written in C++
  # They also avoid the need to convert `distance` to a matrix. 
  if (gaussian) {
    res <- weightedGaussianLocalDensity(weights, distance, attr(distance, "Size"), dc)
  } else {
    res <- weightedNonGaussianLocalDensity(weights, attr(distance, "Size") * sum(weights), distance, attr(distance, "Size"), dc)
  }
  if (is.null(attr(distance, 'Labels'))) {
    names(res) <- NULL
  } else {
    names(res) <- attr(distance, 'Labels')
  }
  res
}

#' @export
distanceToPeak <- function(distance, rho) {
  # This implementation is faster by virtue of being written in C++.
  # It also avoids the need to convert `distance` to a matrix.
  res <- DensityPeak_distancesToPeakC(distance, rho);
  names(res) <- names(rho)
  res
}


## turn 1 distance matrix into i,j coordinates
#' @export
get_ij <- function (k, dist_obj) {
  if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
  n <- attr(dist_obj, "Size")
  valid <- (k >= 1) & (k <= n * (n - 1) / 2)
  k_valid <- k[valid]
  j <- rep.int(NA_real_, length(k))
  j[valid] <- floor(((2 * n + 1) - sqrt((2 * n - 1) ^ 2 - 8 * (k_valid - 1))) / 2)
  i <- j + k - (2 * n - j) * (j - 1) / 2
  cbind(i, j)
}




#' @examples
#' irisDist <- dist(iris[,1:4])
#' DensityPeakEstimateDc(irisDist)


#' @export
DensityPeakEstimateDc <- function(weights, distance, neighborRateLow = 0.01, neighborRateHigh = 0.02) {
  # This implementation uses binary search instead of linear search.
  
  size <- attr(distance, 'Size')
  # If size is greater than 448, there will be >100000 elements in the distance
  # object. Subsampling to 100000 elements will speed performance for very
  # large dist objects while retaining good accuracy in estimating the cutoff
  if (size > 448) {
    distance <- distance[sample.int(length(distance), 100128)]
    size <- 448
  }
  
  low <- min(distance)
  high <- max(distance)
  dc <- 0
  newsize <- sum(weights)
  uniquesize <- size
  size <- newsize
  while (TRUE) {
    dc <- (low + high) / 2
    # neighborRate = average of number of elements of comb per row that are
    # less than dc minus 1 divided by size.
    # This implementation avoids converting `distance` to a matrix. The matrix is
    # symmetrical, so doubling the result from `distance` (half of the matrix) is
    # equivalent. The diagonal of the matrix will always be 0, so as long as dc
    # is greater than 0, we add 1 for every element of the diagonal, which is
    # the same as size
    sum_distance_below_dc <- SumCutOff(weights, distance, attr(distance, "Size"), dc)

    #for (k in 1:uniquesize){

    #  if (distance[k] < dc){
    #    vals <- get_ij(k, distance)
    #    sum_distance_below_dc <- sum_distance_below_dc + (weights[vals[1]]*weights[vals[2]])
    #  }

    #}

    




    neighborRate <- (((sum_distance_below_dc * 2 + (if (0 <= dc) size)) / size - 1)) / size
    if (neighborRate >= neighborRateLow && neighborRate <= neighborRateHigh) break
    
    if (neighborRate < neighborRateLow) {
      low <- dc
    } else {
      high <- dc
    }
  }
  cat('Distance cutoff calculated to', dc, '\n')
  dc
}

#' Calculate clustering attributes based on the DensityPeak Clustering algorithm
#'
#' This function takes a distance matrix and optionally a distance cutoff and
#' calculates the values necessary for clustering based on the algorithm
#' proposed by Alex Rodrigues and Alessandro Laio (see references). The actual
#' assignment to clusters are done in a later step, based on user defined
#' threshold values. If a distance matrix is passed into `distance` the 
#' original algorithm described in the paper is used. If a matrix or data.frame
#' is passed instead it is interpretted as point coordinates and rho will be 
#' estimated based on k-nearest neighbors of each point (rho is estimated as 
#' `exp(-mean(x))` where `x` is the distance to the nearest 
#' neighbors). This can be useful when data is so large that calculating the 
#' full distance matrix can be prohibitive.
#'
#' @details
#' The function calculates rho and delta for the observations in the provided
#' distance matrix. If a distance cutoff is not provided this is first estimated
#' using [DensityPeakEstimateDc()] with default values.
#'
#' The information kept in the Clustering***c object is:
#' \describe{
#'   \item{`rho`}{A vector of local density values}
#'   \item{`delta`}{A vector of minimum distances to observations of higher density}
#'   \item{`distance`}{The initial distance matrix}
#'   \item{`dc`}{The distance cutoff used to calculate rho}
#'   \item{`threshold`}{A named vector specifying the threshold values for rho and delta used for cluster detection}
#'   \item{`peaks`}{A vector of indexes specifying the cluster center for each cluster}
#'   \item{`clusters`}{A vector of cluster affiliations for each observation. The clusters are referenced as indexes in the peaks vector}
#'   \item{`halo`}{A logical vector specifying for each observation if it is considered part of the halo}
#'   \item{`knn_graph`}{kNN graph constructed. It is only applicable to the case where coordinates are used as input. Currently it is set as NA.}
#'   \item{`nearest_higher_density_neighbor`}{index for the nearest sample with higher density. It is only applicable to the case where coordinates are used as input.}
#'   \item{`nn.index`}{indices for each cell's k-nearest neighbors. It is only applicable for the case where coordinates are used as input.}
#'   \item{`nn.dist`}{distance to each cell's k-nearest neighbors. It is only applicable for the case where coordinates are used as input.}
#' }
#' Before running findDensityPeakClusters the threshold, peaks, clusters and halo data is
#' `NA`.
#' 
#' @param distance A distance matrix or a matrix (or data.frame) for the 
#' coordinates of the data. If a matrix or data.frame is used the distances and
#' local density will be estimated using a fast k-nearest neighbor approach.
#' 
#' @param dc A distance cutoff for calculating the local density. If missing it
#' will be estimated with `DensityPeakEstimateDc(distance)`
#'
#' @param gaussian Logical. Should a gaussian kernel be used to estimate the
#' density (defaults to FALSE)
#' 
#' @param verbose Logical. Should the running details be reported  
#'
#' @param ... Additional parameters passed on to [get.knn][FNN::get.knn]
#'
#' @return A Clustering*** object. See details for a description.
#'
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- ClustObj(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#'
#' irisClust <- findDensityPeakClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#'
#' @seealso [DensityPeakEstimateDc()], [findDensityPeakClusters()]
#'
#' @references Rodriguez, A., & Laio, A. (2014). *Clustering by fast search and find of density peaks.* Science, **344**(6191), 1492-1496. doi:10.1126/science.1242072
#'
#' @export
#' 
ClustObj <- function(orig, weights, distance, dc, gaussian=FALSE, verbose = FALSE, ...) {


  #orig = unclass(orig)
  path = paste(getwd(), "/temp.txt", sep = "")
  #write.table(orig, file = path, col.names = F, row.names =F, sep = ",")
  path2 = paste(getwd(), "/temp_weights.txt", sep = "")
  #write.table(weights, file = path2, col.names = F, row.names =F, sep = ",")
  
  new_dist <- (createWeightedDist(distance, attr(distance, 'Size'), weights));

  #norm_weightedDist <- new_dist / (sqrt(sum(new_dist^2)))

  if (class(distance) %in% c('data.frame', 'matrix')) {
    dp_knn_args <- list(mat = distance, verbose = verbose, ...)
    res <- do.call(ClustObj.knn, dp_knn_args)
  } else {
    if (missing(dc)) {
      if (verbose)  message('Calculating the weighted distance cutoff')
      dc <- DensityPeakEstimateDc(weights, distance)
    }
    if (verbose) message('Calculating the weighted local density for each sample based on distance cutoff')
    rho <- localDensity(weights, distance, dc, gaussian = gaussian)
    
    if (verbose) message('Calculating the minimal distance of a weighted sample to another weighted sample with higher density')
    delta <- distanceToPeak(distance, rho)
    
    if (verbose) message('Returning result...')
    res <- list(
      orig = orig,
      size = attr(distance, 'Size'),
      truesize = sum(weights),
      weights = weights,
      fpath = path,
      wpath = path2,
      rho = rho, 
      delta = delta, 
      distance = distance,
      weighted_distance = as.dist(new_dist), 
      dc = dc,
      dbscan_eps = NA,
      dbscan_minpts = NA, 
      threshold = c(rho = NA, delta = NA), 
      peaks = NA,
      centers_spectral = NA,
      mediods = NA, 
      clusters = NA,
      clusters2 = NA,
      clustersDensityPeak = NA,
      clustersSpectral = NA,
      clustersKmediods = NA,
      clustersDBSCAN = NA,
      halo = NA, 
      knn_graph = NA, 
      nearest_higher_density_neighbor = NA, 
      nn.index = NA, 
      nn.dist = NA
    )
    class(res) <- 'ClusteringObject'

  }
  res
}

#' @export
plotDensityPeak <- function(x, ...) {
  UseMethod('plotDensityPeak')
}

#' @export
#' @importFrom graphics plot points
#'
plotDensityPeak.ClusteringObject <- function(x, ...) {
  plot(x$rho, x$delta, main = 'Decision graph', xlab = expression(rho), 
       ylab = expression(delta))
  if (!is.na(x$peaks[1])) {
    points(x$rho[x$peaks], x$delta[x$peaks], col = 2:(1 + length(x$peaks)), 
           pch = 19)
  }
}

#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- DensityPeakCluster(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#'
#' irisClust <- findDensityPeakClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)
#' @export
plotMDS <- function(x, ...) {
  UseMethod('plotMDS')
}

#' @importFrom stats cmdscale
#' @importFrom graphics plot points legend
#' @importFrom stats dist
#' @export
plotMDS.ClusteringObject <- function(x, ...) {
  if (class(x$distance) %in% c('data.frame', 'matrix')) {
    mds <- cmdscale(dist(x$distance))
  } else {
    mds <- cmdscale(x$distance)
  }

  
  if (length(x$peaks) == 1){
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = 'MDS plot of observations')
  } else {
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = 'MDS plot of observations', cex = 0.5, col = 0)
  }
  mds

  #Scale the weights for each point to match their new point size
  if ( max(x$weights)!=  min(x$weights)){
    cex_weights = 2*((x$weights-min(x$weights))/(max(x$weights)-min(x$weights))) + 0.5
  } else {
    cex_weights = (x$weights)/(max(x$weights))*0.5
  }


  if (!is.na(x$peaks[1])) {
    for (i in 1:length(x$peaks)) {
      #print(i)
      ind <- which(x$clusters == i)
      #print(ind)
      #points(mds[ind, 1], mds[ind, 2], col = i + 1, pch = ifelse(x$halo[ind], 1, 19))
      for (index in ind){
        if (index == x$peaks[i]){
          #print("center_found")
          #print(cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (1), pch = 4, cex = cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 2, 17), cex = cex_weights[index])
          
        }
        else {
          #print("other_point")
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 1, 19), cex = cex_weights[index])
        }
        
      }
    }
    legend('topright', legend = c('core', 'halo'), pch = c(19, 1), horiz = TRUE)
  }
  
}
#' @export
plotDensityPeak2D <- function(x, ...) {
  UseMethod('plotDensityPeak2D')
}
#' @export
#' @importFrom stats cmdscale
#' @importFrom graphics plot points legend
#' @importFrom stats dist



plotDensityPeak2D.ClusteringObject <- function(x, ...) {
  
  mds <- x$orig

  
  if (length(x$peaks) == 1){
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = '2D DensityPeak Clustering plot of observations')
  } else {
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = '2D DensityPeak Clustering plot of observations', cex = 0.5, col = 0)
  }
  mds

  #Scale the weights for each point to match their new point size
  if ( max(x$weights)!=  min(x$weights)){
    cex_weights = 2*((x$weights-min(x$weights))/(max(x$weights)-min(x$weights))) + 0.5
  } else {
    cex_weights = (x$weights)/(max(x$weights))*0.5
  }


  if (!is.na(x$peaks[1])) {
    for (i in 1:length(x$peaks)) {
      #print(i)
      ind <- which(x$clusters == i)
      #print(ind)
      #points(mds[ind, 1], mds[ind, 2], col = i + 1, pch = ifelse(x$halo[ind], 1, 19))
      for (index in ind){
        if (index == x$peaks[i]){
          #print("center_found")
          #print(cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (1), pch = 4, cex = cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 2, 17), cex = cex_weights[index])
          
        }
        else {
          #print("other_point")
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 1, 19), cex = cex_weights[index])
        }
        
      }
    }
    legend('topright', legend = c('core', 'halo'), pch = c(19, 1), horiz = TRUE)
  }
  
}
#' Plot observations using t-distributed neighbor embedding and colour by cluster
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- ClustObj(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#' 
#' irisClust <- findDensityPeakClusters(irisClust, rho=2, delta=2)
#' plotTSNE(irisClust)
#' split(iris[,5], irisClust$clusters)
#' @export
plotTSNE <- function(x, ...) {
  UseMethod('plotTSNE')
}



#' @export
#' @importFrom graphics plot points legend
#' @importFrom stats dist
#' @importFrom stats rnorm
#' @importFrom Rtsne Rtsne
plotTSNE.ClusteringObject <- function(x, max_components = 2, ...) {
  if (class(x$distance) %in% c('data.frame', 'matrix')) {
    data <- as.matrix(dist(x$distance))
  } else {
    data <- as.matrix(x$distance)
  } 
  
  # avoid issues related to repetitions
  dup_id <- which(duplicated(data))
  if (length(dup_id) > 0) {
    data[dup_id, ] <- data[dup_id, ] + rnorm(length(dup_id) * ncol(data), sd = 1e-10)
  }
  tsne_res <- Rtsne::Rtsne(as.matrix(data), dims = max_components,
                           pca = T)
  tsne_data <- tsne_res$Y[, 1:max_components]
  
  plot(tsne_data[,1], tsne_data[,2], xlab = '', ylab = '', main = 'tSNE plot of observations')
  if (!is.na(x$peaks[1])) {
    for (i in 1:length(x$peaks)) {
      ind <- which(x$clusters == i)
      points(tsne_data[ind, 1], tsne_data[ind, 2], col = i + 1, pch = ifelse(x$halo[ind], 1, 19))
    }
    legend('topright', legend = c('core', 'halo'), pch = c(19, 1), horiz = TRUE)
  }
}
#' @export
#'


plotSPECTRAL <- function(x, ...) {
  UseMethod('plotSPECTRAL')
}
#' @export
#' @importFrom stats cmdscale
#' @importFrom graphics plot points legend
#' @importFrom stats dist
plotSPECTRAL.ClusteringObject <- function(x, ...) {
  
  mds <- x$orig

  
  if (length(x$centers_spectral) == 1){
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = 'RAW Spectral plot of observations')
  } else {
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = 'RAW Spectral plot of observations', cex = 0.5, col = 0)
  }
  mds

  #Scale the weights for each point to match their new point size
  if ( max(x$weights)!=  min(x$weights)){
    cex_weights = 2*((x$weights-min(x$weights))/(max(x$weights)-min(x$weights))) + 0.5
  } else {
    cex_weights = (x$weights)/(max(x$weights))*0.5
  }


  if (!is.na(x$centers_spectral[1])) {
    for (i in 1:length(x$centers_spectral)) {
      #print(i)
      ind <- which(x$clustersSpectral == i)
      #print(ind)
      #points(mds[ind, 1], mds[ind, 2], col = i + 1, pch = ifelse(x$halo[ind], 1, 19))
      for (index in ind){
        if (index == x$centers_spectral[i]){
          #print("center_found")
          print(cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (1), pch = 4, cex = cex_weights[index])
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 2, 17), cex = cex_weights[index])
          
        }
        else {
          #print("other_point")
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = ifelse(x$halo[index], 1, 19), cex = cex_weights[index])
        }
        
      }
    }
    legend('topright', legend = c('core', 'halo'), pch = c(19, 1), horiz = TRUE)
  }
  
}
#' @export
plot2DCluster <- function(x, ...) {
  UseMethod('plot2DCluster')
}
#' @export
#' @importFrom stats cmdscale
#' @importFrom graphics plot points legend
#' @importFrom stats dist
plot2DCluster.ClusteringObject <- function(x, clusterChoice = x$clusters, type = 'DensityPeak', ...) {
  
  mds <- x$orig

  
  if (length(unique(clusterChoice)) == 1) {
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = paste('2D ', type,  ' Clustering plot of observations'))
  } else {
  plot(mds[,1], mds[,2], xlab = '', ylab = '', main = paste('2D ', type,  ' Clustering plot of observations'), cex = 0.1, col = 0)
  }
  mds

  #Scale the weights for each point to match their new point size
  if ( max(x$weights)!=  min(x$weights)){
    cex_weights = 2*((x$weights-min(x$weights))/(max(x$weights)-min(x$weights))) + 0.5
  } else {
    cex_weights = (x$weights)/(max(x$weights))*0.5
  }


  if (length(unique(clusterChoice)) != 1) {
    for (i in 1:length(unique(clusterChoice))) {
      #print(i)
      ind <- which(clusterChoice == i)
      #print(ind)
      #points(mds[ind, 1], mds[ind, 2], col = i + 1, pch = ifelse(x$halo[ind], 1, 19))
      for (index in ind){
        
        
          #print("other_point")
          points(mds[index, 1], mds[index, 2], col = (i + 1), pch = 19, cex = cex_weights[index])
        
        
      }
    }
    #legend('topright', legend = c('core', 'halo'), pch = c(19, 1), horiz = TRUE)
    
  }
  
}


#' @export
printDensityPeak.ClusteringObject <- function(x, ...) {
  if (is.na(x$peaks[1])) {
    cat('A DensityPeak Cluster object with no clusters defined\n\n')
    cat('Number of observations:', length(x$rho), '\n')
  } else {
    cat('A DensityPeak Cluster object with', length(x$peaks), 'clusters defined\n\n')
    cat('Number of observations:', length(x$rho), '\n')
    cat('Observations in core:  ', sum(!x$halo), '\n\n')
    cat('Parameters:\n')
    cat('dc (distance cutoff)   rho threshold          delta threshold\n')
    cat(formatC(x$dc, width = -22), formatC(x$threshold[1], width = -22), x$threshold[2])
  }
}
#' Detect clusters in a  DensityPeakity Cluster obejct
#' @examples
#' irisDist <- dist(iris[,1:4])
#' irisClust <- ClustObj(irisDist, gaussian=TRUE)
#' plot(irisClust) # Inspect clustering attributes to define thresholds
#'
#' irisClust <- findDensityPeakClusters(irisClust, rho=2, delta=2)
#' plotMDS(irisClust)
#' split(iris[,5], irisClust$clusters)


#' @export
findDensityPeakClusters <- function(x, ...) {
  UseMethod("findDensityPeakClusters")
}

#' @export
findDensityPeakCluster_validationChart <- function(x, ...) {
  UseMethod("findDensityPeakCluster_validationChart")
}

#' @rdname findDensityPeakClusters
#'
#' @param rho The threshold for local density when detecting cluster peaks
#'
#' @param delta The threshold for minimum distance to higher density when detecting cluster peaks
#'
#' @param plot Logical. Should a decision plot be shown after cluster detection
#' 
#' @param peaks A numeric vector indicates the index of density peaks used for clustering. This vector should be retrieved from the decision plot with caution. No checking involved.  
#'
#' @param verbose Logical. Should the running details be reported  
#' @importFrom graphics plot locator
#' @export
findDensityPeakClusters.ClusteringObject <- function(x, rho, delta, plot = FALSE, peaks = NULL, verbose = FALSE, ...) {

  if (class(x$distance) %in% c('data.frame', 'matrix')) {
    peak_ind <- which(x$rho > rho & x$delta > delta)
    x$peaks <- peak_ind
    
    # Assign observations to clusters
    runOrder <- order(x$rho, decreasing = TRUE)
    cluster <- rep(NA, length(x$rho))
    
    #replace certain values in cluster matrix with the cluster centers
    for (i in x$peaks) {
      cluster[i] <- match(i, x$peaks)
    } 

    #for all indexs that arent in the orginal cluster centers
    for (ind in setdiff(runOrder, x$peaks)) { 

      #set target_* to the index where the nearest higher density neighbors of each point are equal to the non cluster center
      target_lower_density_samples <- which(x$nearest_higher_density_neighbor == ind) #all the target cells should have the same cluster id as current higher density cell
      
      cluster[ind] <- cluster[x$nearest_higher_density_neighbor[ind]]
    }
    

    #now the cluster matrix consists of cluster centers [ind] = point and other points of highest near density


    potential_duplicates <- which(is.na(cluster))
    for (ind in potential_duplicates) {
      res <- as.integer(names(which.max(table(cluster[x$nn.index[ind, ]]))))
      
      if (length(res) > 0) {
        cluster[ind] <- res #assign NA samples to the majority of its clusters 
      } else {
        message('try to increase the number of kNN (through argument k) at step of ClustObj.')
        cluster[ind] <- NA
      }
    }
    
    x$clusters <- factor(cluster)
    
    # Calculate core/halo status of observation
    border <- rep(0, length(x$peaks))
    if (verbose) message('Identifying core and halo for each cluster')
    
    for (i in 1:length(x$peaks)) {
      if (verbose) message('the current index of the peak is ', i)
      
      #intersection of 
      connect_samples_ind <- intersect(unique(x$nn.index[cluster == i, ]), which(cluster != i))

      averageRho <- outer(x$rho[cluster == i], x$rho[connect_samples_ind], '+') / 2 

      if (any(connect_samples_ind)) border[i] <- max(averageRho[connect_samples_ind]) 
    }
    x$halo <- x$rho < border[cluster] 
    
    x$threshold['rho'] <- rho
    x$threshold['delta'] <- delta
  } 
  else {
    # Detect cluster peaks
    if (!is.null(peaks)) {
      
      if (verbose) message('peaks are provided, clustering will be performed based on them')
      x$peaks <- peaks
    } else {
      if (missing(rho) || missing(delta)) {
        x$peaks <- NA
        plot(x)
        cat('Click on plot to select thresholds\n')
        threshold <- locator(1)
        if (missing(rho)) rho <- threshold$x
        if (missing(delta)) delta <- threshold$y
        plot = TRUE
      }
      x$peaks <- which(x$rho > rho & x$delta > delta)
      x$threshold['rho'] <- rho
      x$threshold['delta'] <- delta
    }
    if (plot) {
      plot(x)
    }
    
    # Assign observations to clusters
    runOrder <- order(x$rho, decreasing = TRUE)
    cluster <- rep(NA, length(x$rho))
    if (verbose) message('Assigning each sample to a cluster based on its nearest density peak')
    for (i in runOrder) { 
      if ((i %% round(length(runOrder) / 25)) == 0) {
        if (verbose) message(paste('the runOrder index is', i))
      }
      
      if (i %in% x$peaks) {
        cluster[i] <- match(i, x$peaks)
      } else {
        higherDensity <- which(x$rho > x$rho[i])
        cluster[i] <- cluster[higherDensity[which.min(extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), i, higherDensity))]] 
      }
    }
    x$clusters <- cluster
    
    # Calculate core/halo status of observation
    border <- rep(0, length(x$peaks))
    if (verbose) message('Identifying core and halo for each cluster')
    for (i in 1:length(x$peaks)) {
      if (verbose) message('the current index of the peak is ', i)
      
      averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], '+')/2 
      index <- extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), which(cluster == i), which(cluster != i)) <= x$dc 
      if (any(index)) border[i] <- max(averageRho[index]) 
    }
    x$halo <- x$rho < border[cluster] 
  }
  x$halo <- x$rho < border[cluster]
  
  # Sort cluster designations by gamma (= rho * delta)
  gamma <- x$rho * x$delta
  pk.ordr <- order(gamma[x$peaks], decreasing = TRUE)
  x$peaks <- x$peaks[pk.ordr]
  
  x$clusters <- match(x$clusters, pk.ordr)


  if (length(x$peaks) > 1 && (length(x$peaks) < 20) ){

        for (z in 1:x$size){
        x$clusters2[z] = x$peaks[x$clusters[z]]
        }

        #cpath = paste(getwd(), "/temp_cluster.txt", sep = "")
       # write.table(x$clusters2, file = cpath, col.names = F, row.names =F, sep = ",")
        

        
        
        
        #tempDBCV <- DBCV(x$fpath,cpath,x$wpath)
        #print("DBCV is: ")
        #print(tempDBCV)
        


        
      }
  
  x
}


findDensityPeakClusters_dbcv.ClusteringObject <- function(x, rho, delta, plot = FALSE, peaks = NULL, verbose = FALSE, ...) {

  if (class(x$distance) %in% c('data.frame', 'matrix')) {
    peak_ind <- which(x$rho > rho & x$delta > delta)
    x$peaks <- peak_ind
    
    # Assign observations to clusters
    runOrder <- order(x$rho, decreasing = TRUE)
    cluster <- rep(NA, length(x$rho))
    
    #replace certain values in cluster matrix with the cluster centers
    for (i in x$peaks) {
      cluster[i] <- match(i, x$peaks)
    } 

    #for all indexs that arent in the orginal cluster centers
    for (ind in setdiff(runOrder, x$peaks)) { 

      #set target_* to the index where the nearest higher density neighbors of each point are equal to the non cluster center
      target_lower_density_samples <- which(x$nearest_higher_density_neighbor == ind) #all the target cells should have the same cluster id as current higher density cell
      
      cluster[ind] <- cluster[x$nearest_higher_density_neighbor[ind]]
    }
    

    #now the cluster matrix consists of cluster centers [ind] = point and other points of highest near density


    potential_duplicates <- which(is.na(cluster))
    for (ind in potential_duplicates) {
      res <- as.integer(names(which.max(table(cluster[x$nn.index[ind, ]]))))
      
      if (length(res) > 0) {
        cluster[ind] <- res #assign NA samples to the majority of its clusters 
      } else {
        message('try to increase the number of kNN (through argument k) at step of ClustObj.')
        cluster[ind] <- NA
      }
    }
    
    x$clusters <- factor(cluster)
    
    # Calculate core/halo status of observation
    border <- rep(0, length(x$peaks))
    if (verbose) message('Identifying core and halo for each cluster')
    
    for (i in 1:length(x$peaks)) {
      if (verbose) message('the current index of the peak is ', i)
      
      #intersection of 
      connect_samples_ind <- intersect(unique(x$nn.index[cluster == i, ]), which(cluster != i))

      averageRho <- outer(x$rho[cluster == i], x$rho[connect_samples_ind], '+') / 2 

      if (any(connect_samples_ind)) border[i] <- max(averageRho[connect_samples_ind]) 
    }
    x$halo <- x$rho < border[cluster] 
    
    x$threshold['rho'] <- rho
    x$threshold['delta'] <- delta
  } 
  else {
    # Detect cluster peaks
    if (!is.null(peaks)) {
      
      if (verbose) message('peaks are provided, clustering will be performed based on them')
      x$peaks <- peaks
    } else {
      if (missing(rho) || missing(delta)) {
        x$peaks <- NA
        plot(x)
        cat('Click on plot to select thresholds\n')
        threshold <- locator(1)
        if (missing(rho)) rho <- threshold$x
        if (missing(delta)) delta <- threshold$y
        plot = TRUE
      }
      x$peaks <- which(x$rho > rho & x$delta > delta)
      x$threshold['rho'] <- rho
      x$threshold['delta'] <- delta
    }
    if (plot) {
      plot(x)
    }
    
    # Assign observations to clusters
    runOrder <- order(x$rho, decreasing = TRUE)
    cluster <- rep(NA, length(x$rho))
    if (verbose) message('Assigning each sample to a cluster based on its nearest density peak')
    for (i in runOrder) { 
      if ((i %% round(length(runOrder) / 25)) == 0) {
        if (verbose) message(paste('the runOrder index is', i))
      }
      
      if (i %in% x$peaks) {
        cluster[i] <- match(i, x$peaks)
      } else {
        higherDensity <- which(x$rho > x$rho[i])
        cluster[i] <- cluster[higherDensity[which.min(extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), i, higherDensity))]] 
      }
    }
    x$clusters <- cluster
    
    # Calculate core/halo status of observation
    border <- rep(0, length(x$peaks))
    if (verbose) message('Identifying core and halo for each cluster')
    for (i in 1:length(x$peaks)) {
      if (verbose) message('the current index of the peak is ', i)
      
      averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], '+')/2 
      index <- extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), which(cluster == i), which(cluster != i)) <= x$dc 
      if (any(index)) border[i] <- max(averageRho[index]) 
    }
    x$halo <- x$rho < border[cluster] 
  }
  x$halo <- x$rho < border[cluster]
  
  # Sort cluster designations by gamma (= rho * delta)
  gamma <- x$rho * x$delta
  pk.ordr <- order(gamma[x$peaks], decreasing = TRUE)
  x$peaks <- x$peaks[pk.ordr]
  
  x$clusters <- match(x$clusters, pk.ordr)


  if (length(x$peaks) > 1 && (length(x$peaks) < 20) ){

        for (z in 1:x$size){
        x$clusters2[z] = x$peaks[x$clusters[z]]
        }

        cpath = paste(getwd(), "/temp_cluster.txt", sep = "")
        write.table(x$clusters2, file = cpath, col.names = F, row.names =F, sep = ",")
        

        
        
        
        tempDBCV <- DBCV(x$fpath,cpath,x$wpath)
        print("DBCV is: ")
        print(tempDBCV)
        


        
      }
  
  x
}
#' @export
findDensityPeakCluster_validationChart.ClusteringObject <- function(x, rho_step = 0, delta_step = 0, status = FALSE, DBCV = FALSE, plot = FALSE, peaks = NULL, verbose = FALSE, ...) {
  
  #obtain max, min rho
  rho_max = max(x$rho) - 0.01
  rho_min = min(x$rho) + 0.01
  #default rho step size
  if (rho_step == 0){
    rho_step = (rho_max - rho_min)/10
  }
  #obtain max, min delta
  delta_max = max(x$delta) - 0.01
  delta_min = min(x$delta) + 0.01
  #default delta step size
  if (delta_step == 0){
    delta_step = (delta_max - delta_min)/10
  }

  #create data frame
  Rho_Vals <- seq(from = rho_min , to = rho_max , by = rho_step)
  Delta_Vals <- seq(from = delta_min , to = delta_max , by = delta_step)

  if(DBCV == TRUE){
    testClusters <- data.frame(Rho = double(), Delta = double(), Gamma = double(), ClusterCenters = integer(), Unclassified = integer(), NumOutliers = integer(), DBCV = double())
  } else {
    testClusters <- data.frame(Rho = double(), Delta = double(), Gamma = double(), ClusterCenters = integer(), Unclassified = integer(), NumOutliers = integer())

  }
  #implement for loop

  for (rho_temp in Rho_Vals){
    for (delta_temp in Delta_Vals){

      rho = rho_temp
      delta = delta_temp
  
      if (class(x$distance) %in% c('data.frame', 'matrix')) {
        peak_ind <- which(x$rho > rho & x$delta > delta)
        x$peaks <- peak_ind
        
        # Assign observations to clusters
        runOrder <- order(x$rho, decreasing = TRUE)
        cluster <- rep(NA, length(x$rho))
        
        #replace certain values in cluster matrix with the cluster centers
        for (i in x$peaks) {
          cluster[i] <- match(i, x$peaks)
        } 

        #for all indexs that arent in the orginal cluster centers
        for (ind in setdiff(runOrder, x$peaks)) { 

          #set target_* to the index where the nearest higher density neighbors of each point are equal to the non cluster center
          target_lower_density_samples <- which(x$nearest_higher_density_neighbor == ind) #all the target cells should have the same cluster id as current higher density cell
          
          cluster[ind] <- cluster[x$nearest_higher_density_neighbor[ind]]
        }
        

        #now the cluster matrix consists of cluster centers [ind] = point and other points of highest near density


        potential_duplicates <- which(is.na(cluster))
        for (ind in potential_duplicates) {
          res <- as.integer(names(which.max(table(cluster[x$nn.index[ind, ]]))))
          
          if (length(res) > 0) {
            cluster[ind] <- res #assign NA samples to the majority of its clusters 
          } else {
            message('try to increase the number of kNN (through argument k) at step of ClustObj.')
            cluster[ind] <- NA
          }
        }
        
        x$clusters <- factor(cluster)
        
        # Calculate core/halo status of observation
        border <- rep(0, length(x$peaks))
        if (verbose) message('Identifying core and halo for each cluster')
        
        for (i in 1:length(x$peaks)) {
          if (verbose) message('the current index of the peak is ', i)
          
          #intersection of 
          connect_samples_ind <- intersect(unique(x$nn.index[cluster == i, ]), which(cluster != i))

          averageRho <- outer(x$rho[cluster == i], x$rho[connect_samples_ind], '+') / 2 
          
          if (any(connect_samples_ind)) border[i] <- max(averageRho[connect_samples_ind]) 
        }
        x$halo <- x$rho < border[cluster] 
        
        x$threshold['rho'] <- rho
        x$threshold['delta'] <- delta
      } 
      else {
        # Detect cluster peaks
        if (!is.null(peaks)) {
          
          if (verbose) message('peaks are provided, clustering will be performed based on them')
          x$peaks <- peaks
        } else {
          if (missing(rho) || missing(delta)) {
            x$peaks <- NA
            plot(x)
            cat('Click on plot to select thresholds\n')
            threshold <- locator(1)
            if (missing(rho)) rho <- threshold$x
            if (missing(delta)) delta <- threshold$y
            plot = TRUE
          }
          x$peaks <- which(x$rho > rho & x$delta > delta)
          x$threshold['rho'] <- rho
          x$threshold['delta'] <- delta
        }
        if (plot) {
          plot(x)
        }
        
        # Assign observations to clusters
        runOrder <- order(x$rho, decreasing = TRUE)
        cluster <- rep(NA, length(x$rho))
        if (verbose) message('Assigning each sample to a cluster based on its nearest density peak')
        for (i in runOrder) { 
          if ((i %% round(length(runOrder) / 25)) == 0) {
            if (verbose) message(paste('the runOrder index is', i))
          }
          
          if (i %in% x$peaks) {
            cluster[i] <- match(i, x$peaks)
          } else {
            higherDensity <- which(x$rho > x$rho[i])
            cluster[i] <- cluster[higherDensity[which.min(extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), i, higherDensity))]] 
          }
        }
        x$clusters <- cluster
        
        # Calculate core/halo status of observation
        border <- rep(0, length(x$peaks))
        if (verbose) message('Identifying core and halo for each cluster')
        for (i in 1:length(x$peaks)) {
          if (verbose) message('the current index of the peak is ', i)
          
          averageRho <- outer(x$rho[cluster == i], x$rho[cluster != i], '+')/2 
          index <- extractRowColDistfromDistMatrix(x$distance, attr(x$distance, 'Size'), which(cluster == i), which(cluster != i)) <= x$dc 
          if (any(index)) border[i] <- max(averageRho[index]) 
        }
        x$halo <- x$rho < border[cluster] 
      }
      x$halo <- x$rho < border[cluster]
      
      # Sort cluster designations by gamma (= rho * delta)
      gamma <- x$rho * x$delta
      pk.ordr <- order(gamma[x$peaks], decreasing = TRUE)
      x$peaks <- x$peaks[pk.ordr]
      
      x$clusters <- match(x$clusters, pk.ordr)
      
      #assign cluster matrix of raw indices to cluster2
      

      

      if (length(x$peaks) > 1 && (length(x$peaks) < 20) ){

        for (z in 1:x$size){
        x$clusters2[z] = x$peaks[x$clusters[z]]
        }
        if (DBCV == TRUE){


          cpath = paste(getwd(), "/temp_cluster.txt", sep = "")
          write.table(x$clusters2, file = cpath, col.names = F, row.names =F, sep = ",")
        
          #print(x$clusters)
          #message(paste("Running DBCV for cluster number:", length(x$peaks)))
          tempDBCV <- DBCV(x$fpath,cpath,x$wpath)
          #message(paste("DBCV was: ", tempDBCV))


          testClusters[nrow(testClusters) + 1, ] = c(rho, delta, (rho*delta), length(x$peaks), length(x$halo[x$halo == TRUE]), 0, tempDBCV )
        } else {

          cpath = paste(getwd(), "/temp_cluster.txt", sep = "")
          write.table(x$clusters2, file = cpath, col.names = F, row.names =F, sep = ",")
        
          


          testClusters[nrow(testClusters) + 1, ] = c(rho, delta, (rho*delta), length(x$peaks), length(x$halo[x$halo == TRUE]), 0)
        
        }
      }

      if (status == TRUE){
        print(cat("Testing complete for Rho = " , rho , "; Delta = " , delta))
      }
      
    }
  }

  testClusters
}



#' @export
# matrix power operator: computes M^power (M must be diagonalizable)
"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))


#' @export
findSpectralClusters <- function(x, ...) {
  UseMethod("findSpectralClusters")
}
#' @export
findSpectralClusters.ClusteringObject <- function(x, neighbors, ...) {

  #compute an affinity matrix A based on S. A must be made of positive values and be symmetric.

  N  <- nrow(x$orig)
  S <- matrix(rep(NA,N^2), ncol=N)
  og <- x$distance
  #og <- x$weighted_distance
  new <-  attr(og, "Size")

  for(i in 1:N) {
    for(j in 1:N) {
      if(i < j){
        S[i,j] <- 1/ og[(new*(i-1) - i*(i-1)/2 + j-i)]
        S[j,i] <- 1/ og[(new*(i-1) - i*(i-1)/2 + j-i)]
      } else if (i == j){
        S[i,j] <- 0
      }
    }
  }



  N <- length(S[,1])

  if (neighbors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:neighbors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  #A  


  #There is the need of a degree matrix D where each diagonal value is the degree of the respective vertex and all other positions are zero:
  D <- diag(apply(A, 1, sum))

  #compute the unnormalized graph Laplacian (U=D−A) and/or a normalized version (L). The normalized Laplacian
  U <- D - A


  # Assuming we want k clusters, the next step is to find the k smallest eigenvectors (ignoring the trivial constant eigenvector).
  k   <- neighbors
  evL <- eigen(U, symmetric=TRUE)
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]



  #in this transformed space it becomes easy for a standard k-means clustering to find the appropriate clusters:
  km <- stats::kmeans(Z, centers=k, nstart=5)
  plot(x$orig, col=km$cluster)


  
  #print(typeof(km$cluster))
  #print((length(km$cluster)))

  for (z in 1:x$size){
    #print(km$cluster[z])
    x$clustersSpectral[z] = km$cluster[z]
  }
  print(km$centers)
  x$centers_spectral <- km$centers
  


  x
}



#' @export
findKMediodsClusters <- function(x, ...) {
  UseMethod("findKMediodsClusters")
}
#' @export
findKMediodsClusters.ClusteringObject <- function(x, neighbors, ...) {

  dist = x$distance
  kmediods <- cluster::pam(x = dist, k = neighbors)

  x$clustersKmediods <- kmediods$clustering
  x$mediods <- kmediods$mediods
  plot(x$orig, col=x$clustersKmediods)

  x
}




#' @export
findDBSCANClusters <- function(x, eps, ...) {
  UseMethod("findDBSCANClusters")
}
#' @export
findDBSCANClusters.ClusteringObject <- function(x, eps = x$dc, minPts = 5, ...) {


  object <- dbscan::dbscan(x = x$orig, eps = eps , weights = x$weights, minPts = minPts)

  x$dbscan_eps <- object$eps

  x$dbscan_minpts<- object$minPts

  x$clustersDBSCAN <- object$cluster

  plot(x$orig, col=x$clustersDBSCAN)


  x

}

#' Extract DensityPeak cluster membership from a ClusteringObject 
#'
#' This function allows the user to extract the RL Dense cluster membership of all the
#' observations in the given ClusteringObject object. The output can be formatted
#' in two ways as described below. Halo observations can be chosen to be removed
#' from the output.
#'
#' @details
#' Two formats for the output are available. Either a vector of integers
#' denoting for each observation, which cluster the observation belongs to. If
#' halo observations are removed, these are set to NA. The second format is a
#' list with a vector for each group containing the index for the member
#' observations in the group. If halo observations are removed their indexes are
#' omitted. The list format correspond to the following transform of the vector
#' format `split(1:length(clusters), clusters)`, where `clusters` are
#' the cluster information in vector format.
#'
#' @param x The ClusteringObject object. [findDensityPeakClusters()] must have
#' been performed prior to this call to avoid throwing an error.
#'
#' @param ... Currently ignored
#'
#' @return A vector or list with cluster memberships for the observations in the
#' initial distance matrix
#'
#' @export
#'
DensityPeakClusters <- function(x, ...) {
  UseMethod("DensityPeakClusters")
}
#' @rdname DensityPeakClusters
#'
#' @param as.list Should the output be in the list format. Defaults to FALSE
#'
#' @param halo.rm Logical. should halo observations be removed. Defaults to TRUE
#'
#' @export
#'
DensityPeakClusters.ClusteringObject <- function(x, as.list = FALSE, halo.rm = TRUE, ...) {
  if (!DensityPeakClustered(x)) stop('x must be DensityPeakClustered prior to cluster extraction')
  res <- x$clusters
  if (halo.rm) {
    res[x$halo] <- NA
  }
  if (as.list) {
    res <- split(1:length(res), res)
  }
  res
}

#' Check whether a DensityPeake ClusteringObject object have been DensityPeakClustered
#'
#' This function checks whether [findDensityPeakClusters()] has been performed on
#' the given object and returns a boolean depending on the outcome
#'
#' @param x A ClusteringObject object
#'
#' @return `TRUE` if [findDensityPeakClusters()] have been performed, otherwise
#' `FALSE`
#'
#' @export
#'
DensityPeakClustered <- function(x) {
  UseMethod("DensityPeakClustered")
}
#' @rdname DensityPeakClustered
#'
#' @export
#'
DensityPeakClustered.ClusteringObject <- function(x) {
  !any(is.na(x$peaks[1]), is.na(x$clusters[1]), is.na(x$halo[1]))
}

#' Extract labels
#'
#' @noRd
#'
#' @export
#'
labels.ClusteringObject <- function(object, ...) {
  labels(object$distance)
}

#' Fast knn version of ClustObj
#' 
#' This function will be called by ClustObj if a matrix or data.frame is
#' passed in rather than a distance object
#' 
#' @noRd
#' 
#' @importFrom FNN get.knn

ClustObj.knn <- function(mat, k = NULL, verbose = T, ...) {
  if (is.null(k)) {
    k <- round(sqrt(nrow(mat)) / 2) # empirical way to select the number of neighbor points 
    k <- max(10, k) # ensure k is at least 10 
  }  
  
  if (verbose) message('Finding kNN using FNN with ', k, ' neighbors')
  
  dx <- get.knn(mat, k = k, ...)
  
  nn.index <- dx$nn.index
  nn.dist <- dx$nn.dist
  N <- nrow(nn.index)
  
  knn_graph <- NULL
  
  if (verbose)  message('Calculating the local density for each sample based on kNNs ...')
  
  rho <- apply(nn.dist, 1, function(x) {
    exp(-mean(x))
  })
  
  if (verbose) message('Calculating the minimal distance of a sample to another sample with higher density ...')
  
  rho_order <- order(rho)
  
  delta <- vector(mode = 'integer', length = N)
  nearest_higher_density_neighbor <- vector(mode = 'integer', length = N)
  
  
  delta_neighbor_tmp <- DensityPeak_smallest_dist_rho_order_coords(rho[rho_order], as.matrix(mat[rho_order, ]))
  delta[rho_order] <- delta_neighbor_tmp$smallest_dist
  nearest_higher_density_neighbor[rho_order] <- rho_order[delta_neighbor_tmp$nearest_higher_density_sample + 1]
  
  if (verbose) message('Returning result...')
  res <- list(
    rho = rho, 
    delta = delta, 
    distance = mat, 
    dc = NULL, 
    threshold = c(rho = NA, delta = NA), 
    peaks = NA, 
    clusters = NA, 
    halo = NA, 
    knn_graph = knn_graph, 
    nearest_higher_density_neighbor = nearest_higher_density_neighbor, 
    nn.index = nn.index, 
    nn.dist = nn.dist
  )
  class(res) <- 'ClusteringObject'
  res
}


#'@name plotDensityPeakClust
#' @title Plot ClusteringObject results
#' @description Generate a single panel of up to three diagnostic plots for a
#'   \code{ClustObj} object.
#'
#' @param x A ClusteringObject object as produced by \code{\link{ClustObj}}
#' @param type A character vector designating which figures to produce. Valid
#'   options include \code{"dg"} for a decision graph of \eqn{\delta} vs.
#'   \eqn{\rho}, \code{"gg"} for a gamma graph depicting the decrease of
#'   \eqn{\gamma} (= \eqn{\delta} * \eqn{\rho}) across samples, and \code{"mds"},
#'   for a Multi-Dimensional Scaling (MDS) plot of observations. Any combination
#'   of these three can be included in the vector, or to produce all plots,
#'   specify \code{type = "all"}.
#' @param n Number of observations to plot in the gamma graph.
#' @param mds A matrix of scores for observations from a Principal Components
#'   Analysis or MDS. If omitted, and a MDS plot has been requested, one will
#'   be calculated.
#' @param dim.x,dim.y The numbers of the dimensions to plot on the x and y
#'   axes of the MDS plot.
#' @param col Vector of colors for clusters.
#' @param alpha Value in \code{0:1} controlling transparency of points in the
#'   decision graph and MDS plot.
#'
#' @return A panel of the figures specified in \code{type} are produced.
#'   If designated, clusters are color-coded and labelled. If present in
#'   \code{x}, the rho and delta thresholds are designated in the
#'   decision graph by a set of solid black lines.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(iris)
#' data.dist <- dist(iris[, 1:4])
#' pca <- princomp(iris[, 1:4])
#'
#' # Run initial density clustering
#' dens.clust <- ClustObj(data.dist)
#
#' op <- par(ask = TRUE)
#'
#' # Show the decision graph
#' plotDensityPeakClust(dens.clust, type = "dg")
#'
#' # Show the decision graph and the gamma graph
#' plotDensityPeakClust(dens.clust, type = c("dg", "gg"))
#'
#' # Cluster based on rho and delta
#' new.clust <- findDensityPeakClusters(dens.clust, rho = 4, delta = 2)
#'
#' # Show all graphs with clustering
#' plotDensityPeakClust(new.clust, mds = pca$scores)
#'
#' par(op)
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot aes_string geom_text geom_point geom_segment labs
#'   theme_bw theme scale_color_manual geom_line geom_label
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices rainbow
#' @export
#'
plotDensityPeakClust <- function(x, type = "all", n = 20,
                             mds = NULL, dim.x = 1, dim.y = 2,
                             col = NULL, alpha = 0.8) {

  type <- tolower(type)
  if(any(pmatch(type, "all", nomatch = 0))) type <- c("dg", "gg", "mds")

  df <- data.frame(
    rho = x$rho, delta = x$delta, gamma = x$rho * x$delta,
    peaks = FALSE, cluster = factor(x$clusters), halo = x$halo
  )
  df$peaks[x$peaks] <- TRUE

  if(is.null(col)) {
    num.cols <- max(nlevels(df$cluster), 3)
    col <- if(num.cols <= 8) {
      brewer.pal(num.cols, "Set2")
    } else if(num.cols <= 12) {
      brewer.pal(num.cols, "Set3")
    } else rainbow(num.cols + 1)[1:num.cols]
  }

  plots <- list(dg = NULL, gg = NULL, mds = NULL)

  # Plot decision graph (dg)
  if(any(pmatch(type, "dg", nomatch = 0))) {
    plots$dg <- ggplot(df, aes_string(x = "rho", y = "delta"))
    if(!any(is.na(x$threshold))) {
      rho <- x$threshold["rho"]
      delta <- x$threshold["delta"]
      thresh.df <- data.frame(
        x = c(rho, rho),
        y = c(delta, delta),
        xend = c(rho, Inf),
        yend = c(Inf, delta)
      )
      plots$dg <- plots$dg +
        geom_segment(
          aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
          data = thresh.df, inherit.aes = F,
          lineend = "butt"
        )
    }
    if(any(df$peaks)) {
      plots$dg <- plots$dg +
        geom_label(
          aes_string(label = "cluster", color = "cluster"),
          data = df[df$peaks, ],
          fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col)
    }
    plots$dg <- plots$dg  +
      geom_point(
        data = df[!df$peaks, ],
        size = 3, color = "gray50", alpha = alpha
      ) +
      labs(x = expression(rho), y = expression(delta), color = "Cluster") +
      theme(legend.position = "none")
  }

  # Plot gamma graph (gg)
  if(any(pmatch(type, "gg", nomatch = 0))) {
    gg.df <- df[order(df$gamma, decreasing = TRUE), ]
    gg.df <- gg.df[1:n, , drop = FALSE]
    gg.df$Sample <- 1:nrow(gg.df)

    plots$gg <- ggplot(gg.df, aes_string(x = "Sample", y = "gamma")) + geom_line()
    if(any(gg.df$peaks)) {
      plots$gg <- plots$gg +
        geom_label(
          aes_string(label = "cluster", color = "cluster"),
          data = gg.df[gg.df$peaks, , drop = FALSE],
          fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col)
    }
    plots$gg <- plots$gg +
      geom_point(
        data = gg.df[!gg.df$peaks, , drop = FALSE],
        size = 3, color = "gray50"
      ) +
      labs(y = expression(gamma), color = "Cluster") +
      theme(legend.position = "none")
  }

  # Plot MDS (mds)
  if(any(pmatch(type, "mds", nomatch = 0))) {
    if(is.null(mds)) mds <- cmdscale(x$distance, k = max(dim.x, dim.y))
    df$x <- mds[, dim.x]
    df$y <- mds[, dim.y]

    plots$mds <- ggplot()
    plots$mds <- if(all(is.na(df$cluster))) {
      plots$mds +
        geom_point(
          aes_string(x = "x", y = "y"),
          data = df,
          size = 3, color = "gray50", alpha = alpha
        )
    } else {
      plots$mds +
        geom_point(
          aes_string(x = "x", y = "y", color = "cluster"),
          data = df[df$halo, , drop = FALSE],
          shape = 21, size = 3
        ) +
        geom_point(
          aes_string(x = "x", y = "y", color = "cluster"),
          data = df[!df$halo, , drop = FALSE],
          size = 3, alpha = alpha
        ) +
        geom_label_repel(
          aes_string(x = "x", y = "y", label = "cluster", color = "cluster"),
          data = df[df$peaks, , drop = FALSE],
          size = 6, fontface = "bold", alpha = alpha
        ) +
        scale_color_manual(values = col, na.value = "gray50")
    }
    plots$mds <- plots$mds +
      labs(x = paste("Dimension", dim.x), y = paste("Dimension", dim.y)) +
      theme(legend.position = "none")
  }

  has.plot <- !sapply(plots, is.null)
  switch(
    sum(has.plot),
    print(plots[[which(has.plot)]]),
    {
      plots <- plots[has.plot]
      if("mds" %in% names(plots)) plots$nrow <- 2 else plots$ncol <-2
      do.call(grid.arrange, plots)
    },
    {
      plots$layout_matrix <- matrix(c(1, 3, 2, 3), nrow = 2)
      do.call(grid.arrange, plots)
    }
  )
}