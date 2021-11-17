library(ClustPhy)
library(ape)

test_that("plotClustersTree error upon invalid user input", {
  pam <- clustPAM(k = 3, text = NwkTree)
  centers <- pam$medoids
  tree <- pam$phyloTree
  clustering <- pam$clustering
  bad.clustering <- numeric(length(clustering))
  bad.clustering <- bad.clustering[1] <- 3.7

  # tree is not a valid phylo object
  expect_error(plotClustersTree(tree = "not phylo", clustering = clustering))

  # no branch lengths in tree
  expect_error(plotClustersTree(tree = ape::read.tree(text = "(B,(A,C,E),D);"),
                                clustering = clustering))

  # invalid branch lengths
  expect_error(plotClustersTree(tree = ape::read.tree(text = "('B(bad)':6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"),
                                clustering = clustering))

  # redundant leaf names
  expect_error(plotClustersTree(tree = ape::read.tree(text = "(B:6.0,(B:5.0,C:3.0,E:4.0):5.0,D:11.0);"),
                                clustering = clustering))

  # clustering is not a numeric vector
  expect_error(plotClustersTree(tree = tree, clustering = character(length(clustering))))

  # clustering contains non-integer values
  expect_error(plotClustersTree(tree = tree, clustering = bad.clustering))

  # clustering contains negative values
  bad.clustering[1] <- -1
  expect_error(plotClustersTree(tree = tree, clustering = bad.clustering))

  # clustering contains only 0s
  expect_error(plotClustersTree(tree = tree, clustering = numeric(length(clustering))))

  # clustering vector length is not equal to number of leaves
  expect_error(plotClustersTree(tree = tree, clustering = clustering[-1]))

  # col.palette is not a character vector
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                col.palette = numeric(3)))

  # col.palette contains invalid colors
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                col.palette = c("not a color", "", "")))

  # col.palette has fewer colors than the number of clusters
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                col.palette = c("blue", "yellow")))

  # show.center is not a character vector
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                show.center = numeric(3)))

  # show.center contains invalid leaf names
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                show.center = c("bad leafname")))

  # center.symbol is not a character vector
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = 6))

  # symbol.cex is not a number
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = centers,
                                symbol.cex = "not number"))

  # symbol.cex is 0 or less
  expect_error(plotClustersTree(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = centers,
                                symbol.cex = -1))
})

test_that("plotClusters2D error upon invalid user input", {

  pam <- clustPAM(k = 3, text = NwkTree)
  centers <- pam$medoids
  tree <- pam$phyloTree
  clustering <- pam$clustering
  bad.clustering <- numeric(length(clustering))
  bad.clustering <- bad.clustering[1] <- 3.7

  # tree is not a valid phylo object
  expect_error(plotClusters2D(tree = "not phylo", clustering = clustering))

  # no branch lengths in tree
  expect_error(plotClusters2D(tree = ape::read.tree(text = "(B,(A,C,E),D);"),
                                clustering = clustering))

  # invalid branch lengths
  expect_error(plotClusters2D(tree = ape::read.tree(text = "('B(bad)':6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"),
                                clustering = clustering))

  # redundant leaf names
  expect_error(plotClusters2D(tree = ape::read.tree(text = "(B:6.0,(B:5.0,C:3.0,E:4.0):5.0,D:11.0);"),
                                clustering = clustering))

  # less than 2 leaves in tree
  expect_error(plotClusters2D(tree = ape::read.tree(text = "(A);"),
                              clustering = c(1)))

  # clustering is not a numeric vector
  expect_error(plotClusters2D(tree = tree, clustering = character(length(clustering))))

  # clustering contains non-integer values
  expect_error(plotClusters2D(tree = tree, clustering = bad.clustering))

  # clustering contains negative values
  bad.clustering[1] <- -1
  expect_error(plotClusters2D(tree = tree, clustering = bad.clustering))

  # clustering contains only 0s
  expect_error(plotClusters2D(tree = tree, clustering = numeric(length(clustering))))

  # clustering vector length is not equal to number of leaves
  expect_error(plotClusters2D(tree = tree, clustering = clustering[-1]))

  # col.palette is not a character vector
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                col.palette = numeric(3)))

  # col.palette contains invalid colors
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                col.palette = c("not a color", "", "")))

  # col.palette has fewer colors than the number of clusters
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                col.palette = c("blue", "yellow")))

  # show.center is not a character vector
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                show.center = numeric(3)))

  # show.center contains invalid leaf names
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                show.center = c("bad leafname")))

  # center.symbol is not a character vector
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = 6))

  # symbol.cex is not a number
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = centers,
                                symbol.cex = "not number"))

  # symbol.cex is 0 or less
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                show.center = centers,
                                center.symbol = centers,
                                symbol.cex = -1))

  # mexDim is not a number
  expect_error(plotClusters2D(tree = tree,
                                clustering = clustering,
                                maxDim = "not number"))

  # maxDim is less than 2
  expect_error(plotClusters2D(tree = tree,
                              clustering = clustering,
                              maxDim = 1))

  # maxDim is not an integer
  expect_error(plotClusters2D(tree = tree,
                              clustering = clustering,
                              maxDim = 1.5))
})
