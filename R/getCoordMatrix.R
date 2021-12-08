getCoordMatrix <- function(distM, maxDim = NULL) {
  # Get euclidean coordinates from distance matrix
  # Using algorithms outlined in
  # https://math.stackexchange.com/questions/156161/finding-the-coordinates-of-points-from-distance-matrix
  # by Legendre17.
  sq.distM <- distM ^ 2
  # Generate positive semi-finite matrix M
  M <- matrix(data = NA, nrow = nrow(distM), ncol = ncol(distM))
  for (i in 1:nrow(M)) {
    for (j in 1:ncol(M)) {
      M[i, j] <- sq.distM[i, 1] + sq.distM[1, j] - sq.distM[i, j]
    }
  }
  M <- M / 2
  # Use eigendecomposition of M to find coordinates
  eigen.M <- eigen(M, symmetric = TRUE)
  sqrt.eigenvals <- eigen.M$values ^ 0.5
  # Ignore the non-significant dimensions and only keep the valid ones
  valid.dims <- which(! is.na(sqrt.eigenvals))
  coordMatrix <- matrix(data = NA, nrow = nrow(distM), ncol = 0)
  for (i in valid.dims) {
    coordMatrix <- cbind(coordMatrix, eigen.M$vectors[, i] * sqrt.eigenvals[i])
  }
  # Remove the columns that contain only 0s. But make sure the matrix has at
  # least 2 columns
  keep.cols <- ! logical(ncol(coordMatrix))
  for (i in 1:ncol(coordMatrix)) {
    if (all(coordMatrix[, i] == 0)) {
      keep.cols[i] <- FALSE
    }
  }
  if (length(which(keep.cols)) < 2) {
    # keep some columns that contain only 0s so that there is at least 2 columns
    for (i in 1:(2 - length(which(keep.cols)))) {
      keep.cols[which(! keep.cols)[i]] <- TRUE
    }
  }
  coordMatrix <- coordMatrix[, keep.cols]
  # Remove least significant dimensions if there's more dimensions than maxDim
  if (! is.null(maxDim)) {
    if (maxDim < ncol(coordMatrix)) {
      coordMatrix <- coordMatrix[, 1:maxDim]
    }
  }

  return(coordMatrix)
}
