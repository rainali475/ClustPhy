#' Compare different PAM clusterings of the same phylogenetic tree using gap
#' statistics.
#'
#' A function that compares a set of clustering schemes by PAM and computes the
#' gap statistics from a given range of k.
#'
#' @param distM A distance matrix of tree leaves.
#' @param k.max An integer that is at least 2, indicating the max number of
#'     clusters desired for the tree.
#' @param clusterings A list of k.max - 1 positive integer vectors. Each
#'     positive integer vector is usually the clustering output of clustPAM, in
#'     which each leaf node is associated with an integer representing the
#'     cluster it is assigned. Each positive integer vector must be the same
#'     length as the number of leaves in the tree. Clustering vectors must
#'     indicate the clustering scheme of the tree from 2 clusters (k = 2) to
#'     k.max clusters (k = k.max). Clustering vectors must be ordered in the
#'     list by increasing number of clusters. Default value of clusterings is
#'     NULL, in which case PAM clustering will be performed automatically on
#'     distM over 2 <= k <= k.max.
#' @param B An integer indicating the number of bootstrap samples. Default
#'     value is 100.
#'
#' @return
#'
#' @examples
#'
#' @references
#'
#' @export
#' @import cluster
comparePAM <- function(distM, k.max, clusterings = NULL, B = 100L) {

  # Check user input
  if (length(dim(distM)) != 2) {
    stop("distM must be a 2-dimensional distance matrix.")
  }

  if (nrow(distM) != ncol(distM)) {
    stop("distM must be a square matrix indicating distances between each leaf.")
  }

  if (any(is.na(distM))) {
    stop("distM should not contain NA or NaN.")
  }

  if (! isSymmetric(distM)) {
    stop("distM should be a symmetric distance matrix.")
  }

  if (! is.numeric(k.max) | k.max < 2 | k.max != as.integer(k.max)) {
    stop("k.max should be an integer that is at least 2.")
  }

  if (! is.null(clusterings)) {
    if (! is.list(clusterings) | length(clusterings) != k.max - 1) {
      stop("clusterings must be a list of length k.max - 1. ")
    }

    for (i in seq_along(clusterings)) {
      if (! is.integer(clusterings[[i]] |
                       as.integer(clusterings[[i]]) != clusterings[[i]])) {
        stop("Each clustering vector must be an integer vector.")
      }
      if (length(clusterings[[i]]) != nrow(distM)) {
        stop("Length of each clustering vector must be equal to the nrow(distM).")
      }
      if (! all(clusterings[[i]] > 0)) {
        stop("Each clustering vector must contain positive integers (at least 1)
           only.")
      }
      if (length(table(clusterings[[i]])) != i + 1) {
        stop(paste0("There should be ", i + 1, " clusters in clusterings[[", i,
                    "]]."))
      }
    }
  }

  if (! is.numeric(B) | B < 1 | B != as.integer(B)) {
    stop("B should be a positive integer that is at least 1.")
  }

  # Define a clustering function used in clusGap()
  clusteringFunc <- function() {}
  if (is.null(clusterings)) {
    clusteringFunc <- function(x, k) {
      return(list(cluster = cluster::pam(x, k, cluster.only=TRUE, diss = TRUE)))
    }
  } else {
    clusteringFunc <- function(x, k) {
      return(list(cluster = clusterings[[k - 1]]))
    }
  }

  return(cluster::clusGap(distM, FUN = clusteringFunc, K.max = k.max, B = B, verbose = FALSE))
}
#' Compare different EM clustering of the same phylogenetic tree using gap
#' statistics.
#'
#' A function that compares a set of clustering schemes by EM and computes the
#' gap statistics from a given range of k.
#'
#' @param distM A distance matrix of tree leaves.
#' @param k.max An integer that is at least 2, indicating the max number of
#'     clusters desired for the tree.
#' @param clusterings A list of k.max - 1 positive integer vectors. Each
#'     positive integer vector is usually the clustering output of clustEM, in
#'     which each leaf node is associated with an integer representing the
#'     cluster it is assigned. Each positive integer vector must be the same
#'     length as the number of leaves in the tree. Clustering vectors must
#'     indicate the clustering scheme of the tree from 2 clusters (k = 2) to
#'     k.max clusters (k = k.max). Clustering vectors must be ordered in the
#'     list by increasing number of clusters. Default value of clusterings is
#'     NULL, in which case EM clustering will be performed automatically on
#'     distM over 2 <= k <= k.max.
#' @param B An integer indicating the number of bootstrap samples. Default
#'     value is 100.
#'
#' @return
#'
#' @examples
#'
#' @references
#'
#' @export
#' @import cluster
#' @import mclust
compareEM.gap <- function(distM, k.max, clusterings = NULL, B = 100L) {

  # Check user input
  if (length(dim(distM)) != 2) {
    stop("distM must be a 2-dimensional distance matrix.")
  }

  if (nrow(distM) != ncol(distM)) {
    stop("distM must be a square matrix indicating distances between each leaf.")
  }

  if (any(is.na(distM))) {
    stop("distM should not contain NA or NaN.")
  }

  if (! isSymmetric(distM)) {
    stop("distM should be a symmetric distance matrix.")
  }

  if (! is.numeric(k.max) | k.max < 2 | k.max != as.integer(k.max)) {
    stop("k.max should be an integer that is at least 2.")
  }

  if (! is.null(clusterings)) {
    if (! is.list(clusterings) | length(clusterings) != k.max - 1) {
      stop("clusterings must be a list of length k.max - 1. ")
    }

    for (i in seq_along(clusterings)) {
      if (! is.numeric(clusterings[[i]] |
                       as.integer(clusterings[[i]]) != clusterings[[i]])) {
        stop("Each clustering vector must be an integer vector.")
      }
      if (length(clusterings[[i]]) != nrow(distM)) {
        stop("Length of each clustering vector must be equal to the nrow(distM).")
      }
      if (! all(clusterings[[i]] > 0)) {
        stop("Each clustering vector must contain positive integers (at least 1)
           only.")
      }
      if (length(table(clusterings[[i]])) != i + 1) {
        stop(paste0("There should be ", i + 1, " clusters in clusterings[[", i,
                    "]]."))
      }
    }
  }

  if (! is.numeric(B) | B < 1 | B != as.integer(B)) {
    stop("B should be a positive integer that is at least 1.")
  }

  # Define a clustering function used in clusGap()
  clusteringFunc <- function() {}
  if (is.null(clusterings)) {
    clusteringFunc <- function(x, k) {
      return(list(cluster = Mclust(x, G = k, verbose = FALSE)$classification))
    }
  } else {
    clusteringFunc <- function(x, k) {
      return(list(cluster = clusterings[[k - 1]]))
    }
  }

  return(cluster::clusGap(distM, FUN = clusteringFunc, K.max = k.max, B = B, verbose = FALSE))
}
