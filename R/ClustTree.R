#' Cluster phylogenetic tree using PAM
#'
#' A function that calculates optimal tree clustering scheme using PAM method
#' given the number of clusters, file directory for a phylogenetic tree or a
#' newick tree character string.
#'
#' @param k A positive integer greater than 1 indicating the number of clusters.
#' @param file A character string indicating the path to the newick tree.
#'     Default value is "".
#' @param text A character string of the tree in newick format. Default value
#'     is NULL. By default, this is ignored. When this argument is assigned a
#'     value, the argument file is ignored.
#'
#' @return Returns an S3 object of class PAMclusts with results.
#' \itemize{
#'     \item distM - A distance matrix of tree leaves.
#'     \item phyloTree - An S3 object of class phylo, containing tree
#'     information.
#'     \item clustering - A positive integer vector indicating the clusters
#'     assigned to each leaf.
#'     \item medoids - A character vector indicating names of leaves that have
#'     been assigned as cluster centers.
#'     \item stats - A matrix containing parameters for each cluster. The
#'     columns are "size", "max_diss", "av_diss", "diameter", and "separation".
#'     Each row represents a computed cluster.
#'     \itemize{
#'         \item size - number of leaves in the cluster.
#'         \item max_diss - maximum dissimilarity between leaves in the cluster
#'         and center of the cluster.
#'         \item av_diss - average dissimilarity between leaves in the cluster
#'         and center of the cluster.
#'         \item diameter - maximum dissimilarity between two leaves in the
#'         cluster.
#'         \item separation - minimum dissimilarity between one leaf of this
#'         cluster and one leaf of another cluster.
#'     }
#' }
#'
#' @examples
#' # Make 6 clusters from a newick tree using PAM
#' pam <- clustPAM(6, text = NwkTree2)
#'
#' @references
#' Kaufman, L., &amp; Rousseeuw, P. J. (2005). Finding groups in data:
#' An introduction to cluster analysis. Wiley.
#'
#' @export
#' @importFrom cluster pam
#' @import ape
clustPAM <- function(k, file = "", text = NULL){

  # Perform checks for user input
  if (is.numeric(k) == FALSE | k <= 1 | k != as.integer(k)) {
    stop("k must be a positive integer greater than 1.")
  }

  phyloTree <- getPhyloTree(file, text)

  # Check if there is more clusters required than the number of leaves in the
  # tree.
  if (k > length(phyloTree$tip.label)) {
    stop("There cannot be more clusters than the number of leaves in the tree.
         Please decrease k.")
  }

  # Get distance matrix from distances between tree leaves
  distMatrix <- ape::cophenetic.phylo(phyloTree)

  # Use PAM to cluster the tree leaves
  clusterResult <- pam(x = distMatrix, k = k, diss = TRUE)

  # Create result S3 object of class PAMclusts
  PAMresults <- list(distM = distMatrix,
                     phyloTree = phyloTree,
                     clustering = clusterResult$clustering,
                     medoids = clusterResult$medoids,
                     stats = clusterResult$clusinfo)
  class(PAMresults) <- "PAMclusts"
  return(PAMresults)
}

#' Cluster phylogenetic tree using EM
#'
#' A function that calculates optimal tree clustering scheme using EM method
#' given the number of clusters, file directory for a phylogenetic tree or a
#' newick tree character string, and an optional set of labeled leaves for
#' semi-supervised clustering.
#'
#' @param k A positive integer greater than 1 indicating the number of clusters.
#' @param file A character string indicating the path to the newick tree.
#'     Default value is "".
#' @param text A character string of the tree in newick format. Default value
#'     is NULL. By default, this is ignored. When this argument is assigned a
#'     value, the argument file is ignored.
#' @param maxDim A positive integer that is at least 2, indicating the max
#'     dimension of the coordinates that represent tree leaves. The dimensions
#'     will be less than or equal to the number of leaves. The default value of
#'     maxDim is NULL, for which the full dimension will be used. maxDim can be
#'     set smaller to decrease runtime of PCA at the cost of discarding the
#'     least significant dimensions beyond maxDim.
#' @param maxPC A positive integer that is at least 2 indicating the maximum
#'     number of dimensions of the reduced coordinates of the tree leaves after
#'     PCA used towards EM clustering. The dimensions of the reduced coordinates
#'     will never exceedthe number of leaves in the tree when maxDim is NULL and
#'     will be at most maxDim if maxDim is given. Default value of maxPC is 5,
#'     as usually most of the variance in the data can be explained by the top
#'     5 principle components. Including too many dimensions can lead to sparse
#'     datapoints and prevent efficient clustering.
#'
#' @return Returns an S3 object of class EMclusts with results.
#' \itemize{
#'     \item distM - A distance matrix of tree leaves.
#'     \item phyloTree - An S3 object of class phylo, containing tree
#'     information.
#'     \item clustering - A positive integer vector indicating the clusters
#'     assigned to each leaf.
#'     \item mean - A numeric matrix with columns corresponding to mean
#'     coordinates of cluster centers.
#'     \item bic - A numeric indicating BIC value of the optimal model.
#'     \item model - A character string describing the optimal model used.
#'     \item dimredResult - An S3 object of class treeDimred with the coordinate
#'     matrix and PCA results.
#' }
#'
#' @examples
#' # Make 6 clusters from a newick tree using EM
#' em <- clustEM(6, text = NwkTree2)
#'
#' @author {Yuzi Li, \email{rainal.li@mail.utoronto.ca}}
#'
#' @references
#' Kaufman, L., &amp; Rousseeuw, P. J. (2005). Finding groups in data:
#' An introduction to cluster analysis. Wiley.
#'
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics
#' and evolutionary analyses in R.” Bioinformatics, 35, 526-528.
#'
#' Scrucca L, Fop M, Murphy TB, Raftery AE (2016). “mclust 5: clustering,
#' classification and density estimation using Gaussian finite mixture models.”
#' The R Journal, 8(1), 289–317. https://doi.org/10.32614/RJ-2016-021.
#'
#' @export
#' @import ape
#' @import mclust
clustEM <- function(k, file = "", text = NULL, maxDim = NULL, maxPC = 5) {

  # Check user input
  if (is.numeric(k) == FALSE | k <= 1 | k != as.integer(k)) {
    stop("k must be a positive integer greater than 1.")
  }

  phyloTree <- getPhyloTree(file, text)

  # Check if there is more clusters required than the number of leaves in the
  # tree.
  if (k > length(phyloTree$tip.label)) {
    stop("There cannot be more clusters than the number of leaves in the tree.
         Please decrease k.")
  }

  # Check maxDim and maxPC
  if (! is.null(maxDim)) {
    if (! is.numeric(maxDim)) {
      stop("maxDim must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
    } else if (maxDim < 2 | as.integer(maxDim) != maxDim) {
      stop("maxDim must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
    }
  }

  if (! is.numeric(maxPC)) {
    stop("maxPC must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
  } else if (maxPC < 2 | as.integer(maxPC) != maxPC) {
    stop("maxPC must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
  }

  # Get distance matrix from distances between tree leaves
  distMatrix <- ape::cophenetic.phylo(phyloTree)

  # Get coordinate matrix from distance matrix
  coordMatrix <- getCoordMatrix(distMatrix, maxDim)

  # Do PCA on coordinate matrix to reduce dimensions
  pca <- prcomp(coordMatrix, retx = TRUE)

  # Use the desired number of dimensions for EM cluster input
  reducedMatrix <- pca$x
  if (ncol(reducedMatrix) > maxPC) {
    reducedMatrix <- reducedMatrix[, 1:maxPC]
  }

  # Use expectation maximization to cluster the coordinates\
  em.mod <- mclust::Mclust(reducedMatrix, G = k, verbose = FALSE)

  # Interpretation of model names
  modNames <- list(E =  "equal variance (univariate)",
                   V = "variable variance (univariate)",
                   EII = "spherical, equal volume",
                   VII = "spherical, unequal volume",
                   EEI = "diagonal, equal volume and shape",
                   VEI = "diagonal, varying volume, equal shape",
                   EVI = "diagonal, equal volume, varying shape",
                   VVI = "diagonal, varying volume and shape",
                   EEE = "ellipsoidal, equal volume, shape, and orientation",
                   EEV = "ellipsoidal, equal volume and equal shape",
                   VEV = "ellipsoidal, equal shape",
                   VVV = "ellipsoidal, varying volume, shape, and orientation")

  # Create an object of class EMclusts with results
  dimredResult <- list(coordM = coordMatrix,
                       PCA = pca)
  class(dimredResult) <- "treeDimRed"
  EMresults <- list(distM = distMatrix,
                    phyloTree = phyloTree,
                    clustering = em.mod$classification,
                    mean = unname(em.mod$parameters$mean),
                    bic = em.mod$bic,
                    model = modNames[[em.mod$modelName]],
                    dimredResult = dimredResult)
  class(EMresults) <- "EMclusts"
  return(EMresults)
}

getPhyloTree <- function(file = "", text = NULL) {

  # Check user inputs
  if (is.null(text)) {
    if (file.exists(file) == FALSE) {
      stop("file should be a valid file path.")
    }
  } else if (is.character(text) == FALSE) {
    stop("text should be a character string of tree in newick format.")
  }

  # Read the tree to an S3 object of class "phylo"
  phyloTree <- ape::read.tree(file = file, text = text)

  # Check if the tree contains branch lengths
  if (! "edge.length" %in% names(phyloTree)) {
    stop("The input newick tree must contain branch lengths.")
  }

  if (any(is.na(phyloTree$edge.length))) {
    stop("There must be a valid branch length for every edge in the tree.
         Also check tree format: brackets cannot be included in node names.")
  }

  # Check if the tree contains duplicate leaf names
  if (any(table(phyloTree$tip.label) > 1)) {
    stop("The input newick tree must not contain duplicate leaf names.")
  }

  return(phyloTree)
}

#[END]
