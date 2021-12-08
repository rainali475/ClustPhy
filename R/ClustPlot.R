#' Plot phylogenetic tree clusters as phylogram
#'
#' A function that plots a phylogenetic tree with a given clustering scheme and
#' given options of coloring scheme and whether to show the leaf names.
#'
#' @param tree A phylo tree object.
#' @param clustering An positive integer vector, usually the from the output of
#'     clustPAM, indicating index of specific colors from the palette. The
#'     integers in this vector must be at most the length of palette. Values of
#'     0 indicate that the node is not assigned a cluster and therefore will not
#'     be labeled. If show.centers is not NULL, then all cluster centers
#'     indicated in the show.centers must not have a clustering value of 0.
#'     clustering should be the same length as the number of leaves in the tree.
#' @param col.palette A character vector indicating the colors used for different
#'     clusters. If no value is given, the palette will be generated from
#'     rainbow.
#' @param show.centers A character vector indicating leaf names of the cluster
#'     centers. Default value is NULL.
#' @param center.symbol A character vector indicating the symbol or text
#'     representative of the cluster centers. Default value is " * ". If the length
#'     of center.symbol is less than the number of cluster centers specified in
#'     show.centers, then the symbols are recycled. This argument is ignored
#'     when show.centers is NULL.
#' @param symbol.cex A numeric indicating factor scaling of center symbols.
#'     Default value is 1. This can be a numeric vector indicating scaling of
#'     each symbol. If the vector is shorter than the number of cluster centers,
#'     the numbers will be recycled. This argument is ignored if show.centers is
#'     NULL.
#' @param node.cex A numeric indicating factor scaling of node circles. Default
#'     value is 1. This can be a numeric vector indicating scaling of each
#'     node (that is not a cluster center). If the vector is shorter than the
#'     number of nodes, the numbers will be recycled.
#' @param title A character string indicating the title of the plot. Default
#'     value is "Tree Clusters".
#'
#' @return Returns a phylogram with different colors indicating different
#'     clusters.
#'
#' @examples
#' # Example 1
#' # Make 6 clusters from a newick tree using PAM
#' pam <- clustPAM(6, text = NwkTree2)
#' # Plot the pam clusters
#' plotClustersTree(pam$phyloTree,
#'                  pam$clustering,
#'                  show.centers = pam$medoids,
#'                  center.symbol = pam$medoids)
#'
#' # Example 2
#' # Make 6 clusters from a newick tree using EM
#' em <- clustEM(6, text = NwkTree2)
#' # Plot the em clusters
#' plotClustersTree(em$phyloTree, em$clustering)
#'
#' @export
#' @import ape
plotClustersTree <- function(tree,
                         clustering,
                         col.palette = NULL,
                         show.centers = NULL,
                         center.symbol = " * ",
                         symbol.cex = 1,
                         node.cex = 1,
                         title = "Tree Clusters") {

  # Check user inputs
  col.palette <- getColPalette(tree,
                               clustering,
                               col.palette,
                               show.centers,
                               center.symbol,
                               symbol.cex)

  # Check title
  if (! is.character(title)) {
    stop("title must be a character string indicating the plot title.")
  }

  # Check if the tree contains branch lengths
  if (! "edge.length" %in% names(tree)) {
    stop("The input tree must contain branch lengths.")
  }

  if (any(is.na(tree$edge.length))) {
    stop("There must be a valid branch length for every edge in the tree.
         Also check tree format: brackets cannot be included in node names.")
  }

  # Check if the tree contains duplicate leaf names
  if (any(table(tree$tip.label) > 1)) {
    stop("The input tree must not contain duplicate leaf names.")
  }

  # Plot tree
  ape::plot.phylo(tree, type = "fan", show.tip.label = FALSE, main = title)

  # Get indices of cluster centers
  centerIndices <- which(tree$tip.label %in% show.centers)
  if (! is.null(show.centers)) {
    if (any(clustering[centerIndices] == 0)) {
      warning("All cluster centers must belong to a valid cluster. Ones that
              do not belong to a valid cluster are ignored.",
              call. = FALSE)
    }
  }

  # Get the indices of leaves that are assigned to clusters and are not cluster
  # centers
  clusterMembers <- which(clustering > 0)
  if (! is.null(show.centers)) {
    clusterMembers <- clusterMembers[! clusterMembers %in% centerIndices]
  }

  # Remove leaves without a cluster or are centers from clustering
  leaf.clustering <- clustering[clusterMembers]

  # Add colored tip labels for leaves
  leaf.colors <- col.palette[leaf.clustering]
  ape::tiplabels(text = character(length(leaf.clustering)),
                 tip = clusterMembers,
                 frame = "circle",
                 bg = leaf.colors,
                 cex = node.cex)

  # Add colored tip labels for centers
  if (! is.null(show.centers)) {
    # Organize texts and cex for each leaf label
    center.text <- character(length(show.centers))
    center.cex <- sample(1, length(show.centers), replace = TRUE)
    for (i in seq_along(show.centers)) {
      # Get index of current symbol in center.symbol
      iSymbol <- i %% length(center.symbol)
      if (iSymbol == 0) {
        iSymbol <- length(center.symbol)
      }
      # Get index of current cex in symbol.cex
      iCex <- i %% length(symbol.cex)
      if (iCex == 0) {
        iCex <- length(symbol.cex)
      }
      center.text[i] <- center.symbol[iSymbol]
      center.cex[i] <- symbol.cex[iCex]
    }
    # Get center colors
    center.clustering <- clustering[centerIndices]
    center.colors <- col.palette[center.clustering]
    # Remove centers that are not assigned to a valid cluster (clustering = 0)
    kept.centers <- center.clustering != 0
    center.clustering <- center.clustering[kept.centers]
    center.colors <- center.colors[kept.centers]
    centerIndices <- centerIndices[kept.centers]
    center.text <- center.text[kept.centers]
    center.cex <- center.cex[kept.centers]
    # Add colored tip labels for centers
    ape::tiplabels(text = center.text,
                   tip = centerIndices,
                   frame = "rect",
                   bg = center.colors,
                   cex = center.cex)
  }
}

#' Plot phylogenetic tree clusters as a 2D graph
#'
#' A function that plots clusters on a 2D graph using dimensionality reduction
#' by PCA, with a given clustering and options of coloring scheme and whether
#' to show the leaf names.
#'
#' @param tree A phylo tree object.
#' @param clustering An positive integer vector, usually the from the output of
#'     clustEM, indicating index of specific colors from the palette. The
#'     integers in this vector must be at most the length of palette. Values of
#'     0 indicate that the node is not assigned a cluster and therefore will not
#'     be labeled. If show.centers is not NULL, then all cluster centers
#'     indicated in the show.centers must not have a clustering value of 0.
#'     clustering should be the same length as the number of leaves in the tree.
#' @param col.palette A character vector indicating the colors used for different
#'     clusters. If no value is given, the palette will be generated from
#'     rainbow.
#' @param show.centers A character vector indicating leaf names of the cluster
#'     centers. Default value is NULL.
#' @param center.symbol A character vector indicating the symbol or text
#'     representative of the cluster centers. Default value is " * ". If the length
#'     of center.symbol is less than the number of cluster centers specified in
#'     show.centers, then the symbols are recycled. This argument is ignored
#'     when show.centers is NULL.
#' @param symbol.cex A numeric indicating factor scaling of center symbols.
#'     Default value is 1. This can be a numeric vector indicating scaling of
#'     each symbol. If the vector is shorter than the number of cluster centers,
#'     the numbers will be recycled. This argument is ignored if show.centers is
#'     NULL.
#' @param node.cex A numeric indicating factor scaling of node circles. Default
#'     value is 1. This can be a numeric vector indicating scaling of each
#'     node (that is not a cluster center). If the vector is shorter than the
#'     number of nodes, the numbers will be recycled.
#' @param maxDim A positive integer that is at least 2, indicating the max
#'     dimension of the coordinates that represent tree leaves. The dimensions
#'     will be less than or equal to the number of leaves. The default value of
#'     maxDim is NULL, for which the full dimension will be used. maxDim can be
#'     set smaller to decrease runtime of PCA at the cost of discarding the
#'     least significant dimensions beyond maxDim.
#' @param title A character string indicating the title of the plot. Default
#'     value is "Tree Clusters Biplot".
#'
#' @return Returns a biplot of different clusters with different colors and
#'     an S3 object of class treeDimred with the coordinate matrix and PCA results.
#' \itemize{
#'     \item coordM - The coordinate matrix representation of tree leaves in
#'     n-dimensional space.
#'     \item PCA - A list of 5 elements containing results from principle
#'     component analysis.
#' }
#'
#' @examples
#' pam <- clustPAM(6, text = NwkTree2)
#' dimRedInfoPAM <- plotClusters2D(pam$phyloTree,
#'                                 pam$clustering,
#'                                 show.centers = pam$medoids,
#'                                 center.symbol = pam$medoids)
#'
#' @export
#' @import ape
#' @importFrom graphics points text
#' @importFrom stats prcomp
plotClusters2D <- function(tree,
                           clustering,
                           col.palette = NULL,
                           show.centers = NULL,
                           center.symbol = " * ",
                           symbol.cex = 1,
                           node.cex = 1,
                           maxDim = NULL,
                           title = "Tree Clusters Biplot") {
  # Check user inputs
  col.palette <- getColPalette(tree,
                               clustering,
                               col.palette,
                               show.centers,
                               center.symbol,
                               symbol.cex,
                               node.cex)

  # Check title
  if (! is.character(title)) {
    stop("title must be a character string indicating the plot title.")
  }

  if (! is.null(maxDim)) {
    if (! is.numeric(maxDim)) {
      stop("maxDim must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
    } else if (maxDim < 2 | as.integer(maxDim) != maxDim) {
      stop("maxDim must be a positive integer indicating the maximum
         dimension of the coordinates representation of tree leaves.")
    }
  }

  # Check if the tree contains branch lengths
  if (! "edge.length" %in% names(tree)) {
    stop("The input tree must contain branch lengths.")
  }

  if (any(is.na(tree$edge.length))) {
    stop("There must be a valid branch length for every edge in the tree.
         Also check tree format: brackets cannot be included in node names.")
  }

  # Check if the tree contains duplicate leaf names
  if (any(table(tree$tip.label) > 1)) {
    stop("The input tree must not contain duplicate leaf names.")
  }

  # Check if the tree contains less than 2 leaves
  if (length(tree$tip.label) < 2) {
    stop("There must be at least 2 leaves in the input tree.")
  }

  # Get distance matrix from distances between tree leaves
  distMatrix <- ape::cophenetic.phylo(tree)

  coordMatrix <- getCoordMatrix(distMatrix, maxDim)

  # Perform PCA on the coordinates matrix
  tree.pca <- prcomp(coordMatrix, retx = TRUE)

  # Get values for the first 2 dimensions
  dim1 <- tree.pca$x[, 1]
  dim2 <- tree.pca$x[, 2]

  # Generate an empty plot
  dim1.pvar <- summary(tree.pca)$importance["Proportion of Variance", 1]
  dim1.pvar <- round(dim1.pvar * 100, digits = 1)
  dim2.pvar <- summary(tree.pca)$importance["Proportion of Variance", 2]
  dim2.pvar <- round(dim2.pvar * 100, digits = 1)
  plot(x = c(),
       y = c(),
       xlim = c(min(dim1), max(dim1)),
       ylim = c(min(dim2), max(dim2)),
       xlab = paste0("Dimension 1 (",
                     dim1.pvar,
                     "% explained var.)"),
       ylab = paste0("Dimension 2 (",
                     dim2.pvar,
                     "% explained var.)"),
       main = title)

  # Add points from each cluster to the plot
  numNodes <- length(which(clustering != 0))
  if (! is.null(show.centers)) {
    numNodes <- numNodes - length(show.centers)
  }
  node.cex <- rep(node.cex,length.out = numNodes)
  clusters <- as.integer(names(table(clustering)))
  clusters <- clusters[clusters != 0]
  start.node <- 1
  end.node <- 0
  for (i in clusters) {
    in.clust <- clustering == i
    end.node <- end.node + length(which(in.clust))
    points(dim1[in.clust],
           dim2[in.clust],
           pch = 21,
           bg = col.palette[i],
           cex = node.cex[start.node:end.node])
    start.node <- start.node + length(which(in.clust))
  }

  # Add center symbols
  if (! is.null(show.centers)) {
    # Get indices of cluster centers
    centerIndices <- which(tree$tip.label %in% show.centers)
    if (any(clustering[centerIndices] == 0)) {
      warning("All cluster centers must belong to a valid cluster. Ones that
              do not belong to a valid cluster are ignored.",
              call. = FALSE)
      centerIndices <- centerIndices[clustering[centerIndices] != 0]
    }
    # Label the centers
    for (i in seq_along(centerIndices)) {
      iSymbol <- i %% length(center.symbol)
      if (iSymbol == 0) {
        iSymbol <- length(center.symbol)
      }
      iCex <- i %% length(symbol.cex)
      if (iCex == 0) {
        iCex <- length(symbol.cex)
      }
      curr.text <- center.symbol[iSymbol]
      curr.cex <- symbol.cex[iCex]
      text(x = dim1[centerIndices[i]],
           y = dim2[centerIndices[i]],
           labels = curr.text,
           pos = 3,
           cex = curr.cex)
    }
  }

  # Return the coordinate matrix and PCA results in an S3 object of class
  # "treeDimRed"
  dimredResult <- list(coordM = coordMatrix,
                       PCA = tree.pca)
  class(dimredResult) <- "treeDimRed"
  return(dimredResult)
}

getColPalette <- function(tree,
                          clustering,
                          col.palette = NULL,
                          show.centers = NULL,
                          center.symbol = " * ",
                          symbol.cex = 1,
                          node.cex = 1) {
  # Perform checks for user input
  if (class(tree) != "phylo") {
    stop("tree should be an object of class phylo, as outputed by
         ape::read.tree.")
  }

  if (! is.numeric(clustering) | any(as.integer(clustering) != clustering)) {
    stop("clustering should be an integer vector indicating index of specific
         colors from the col.palette.")
  }

  if (length(clustering) != length(tree$tip.label)) {
    stop("clustering should be the same length as the number of leaves in the
         tree.")
  }

  if (! is.null(col.palette) & ! is.character(col.palette)) {
    stop("col.palette should be a character vector indicating colors for each
         cluster.")
  } else if (is.null(col.palette)) {
    col.palette <- grDevices::rainbow(max(clustering))
  }

  if (max(clustering) > length(col.palette) | min(clustering) < 0) {
    stop("clustering contains invalid index for color in col.palette.")
  }

  if (all(clustering == 0)) {
    stop("clustering cannot be all 0s.")
  }

  if (! is.null(show.centers) & ! is.character(show.centers)) {
    stop("show.centers must be a character vector indicating the leaf names
           of the centers.")
  }

  if (! is.null(show.centers) & ! all(show.centers %in% tree$tip.label)) {
    stop("show.centers contain invalid leaf names.")
  }

  if (! is.character(center.symbol)) {
    stop("center.symbol must be a character vector indicating symbols used to
         represent cluster centers.")
  }

  if (! is.numeric(symbol.cex)) {
    stop("symbol.cex must be a numeric indicating factor scaling of center
         symbols.")
  }

  if (any(symbol.cex <= 0)) {
    stop("symbol.cex contains invalid values.")
  }

  if (! is.numeric(node.cex)) {
    stop("node.cex must be a numeric indicating factor scaling of center
         symbols.")
  }

  if (any(node.cex <= 0)) {
    stop("node.cex contains invalid values.")
  }

  return(col.palette)
}
