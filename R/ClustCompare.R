#' Compare different clusterings of the same phylogenetic tree using gap
#' statistics.
#'
#' A function that compares a set of clustering schemes by PAM or EM and
#' computes the gap statistics from a given range of k.
#'
#' @param distM A distance matrix of tree leaves.
#' @param k.max An integer that is at least 2, indicating the max number of
#'     clusters desired for the tree.
#' @param clusterings A list of k.max - 1 positive integer vectors. Each
#'     positive integer vector is usually the clustering output of clustPAM or
#'     clustEM, in which each leaf node is associated with an integer
#'     representing the cluster it is assigned. Each positive integer vector
#'     must be the same length as the number of leaves in the tree. Clustering
#'     vectors must indicate the clustering scheme of the tree from 2 clusters
#'     (k = 2) to k.max clusters (k = k.max). Clustering vectors must be
#'     ordered in the list by increasing number of clusters. Default value of
#'     clusterings is NULL, in which case EM or PAM clustering will be performed
#'     on distM over 2 <= k <= k.max depending on the method argument.
#' @param B An integer indicating the number of bootstrap samples. Default
#'     value is 100.
#' @param method A character string that is either "PAM" or "EM", indicating
#'     the clustering method in case clustering is not given. This argument is
#'     ignored if clustering is given. Default value is "PAM".
#'
#' @return An S3 object of "clusGap" with gap statistics.
#'
#' @examples
#'
#' @references
#'
#' @export
#' @import cluster
compareGap <- function(distM,
                       k.max,
                       clusterings = NULL,
                       B = 100L,
                       method = "PAM") {

  # Check user input
  if (! is.numeric(distM) | ! is.matrix(distM)){
    stop("distM must be a numeric matrix.")
  }

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
      if (! is.integer(clusterings[[i]]) &
                       any(as.integer(clusterings[[i]]) != clusterings[[i]])) {
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
  } else if (! all(method %in% c("PAM", "EM")) | length(method) != 1) {
    stop("method should be a character string of one of \"PAM\" or \"EM\".")
  }

  if (! is.numeric(B)) {
    stop("B should be a positive integer that is at least 1.")
  } else if (B < 1 | B != as.integer(B)) {
    stop("B should be a positive integer that is at least 1.")
  }

  # Define a clustering function used in clusGap()
  # clusteringFunc <- function() {}
  if (is.null(clusterings)) {
    clusterings <- vector(mode = "list", length = k.max - 1)
    if (method == "PAM") {
      for (i in 2:k.max) {
        clustering <- cluster::pam(distM, i, cluster.only=TRUE, diss = TRUE)
        clusterings[[i - 1]] <- clustering
      }
      # clusteringFunc <- function(x, k) {
      #   return(list(cluster = cluster::pam(x, k, cluster.only=TRUE, diss = TRUE)))
      # }
    } else { # method is EM
      for (i in 2:k.max) {
        clustering <- Mclust(distM, G = c(i), verbose = FALSE)$classification
        clusterings[[i - 1]] <- clustering
      }
      # clusteringFunc <- function(x, k) {
      #   return(list(cluster = Mclust(x, G = c(k), verbose = FALSE)$classification))
      # }
    }
  }
  clusteringFunc <- function(x, k) {
    if (k == 1) {
      return(list(cluster = numeric(length = nrow(distM))))
    }
    return(list(cluster = clusterings[[k - 1]]))
  }

  return(cluster::clusGap(distM, FUN = clusteringFunc, K.max = k.max, B = B, verbose = FALSE))
}
