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
#'
#' @references
#'
#' @export
#' @import cluster
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
  clusterResult <- cluster::pam(x = distMatrix, k = k, diss = TRUE)

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
#' }
#'
#' @examples
#'
#' @references
#'
#' @export
#' @import ape
#' @import mclust
clustEM <- function(k, file = "", text = NULL) {

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

  # Get distance matrix from distances between tree leaves
  distMatrix <- ape::cophenetic.phylo(phyloTree)

  # Use expectation maximization to cluster the coordinates\
  em.mod <- mclust::Mclust(distMatrix, G = k, verbose = FALSE)

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
  EMresults <- list(distM = distMatrix,
                    phyloTree = phyloTree,
                    clustering = em.mod$classification,
                    mean = unname(em.mod$parameters$mean),
                    bic = em.mod$bic,
                    model = modNames[[em.mod$modelName]])
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
