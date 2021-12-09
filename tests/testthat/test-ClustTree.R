library(ClustPhy)

test_that("PAM clustering with 3 clusters from a character string of newick tree", {

  pam <- clustPAM(k = 3, text = NwkTree)

  distM <- c(0, 16, 14, 15, 17,
             16, 0, 8, 9, 21,
             14, 8, 0, 7, 19,
             15, 9, 7, 0, 20,
             17, 21, 19, 20, 0)
  distM <- matrix(distM, nrow = 5, byrow = TRUE)
  clustering <- c(1, 2, 2, 2, 3)
  medoids <- c("B", "C", "D")
  stats <- c(1, 0, 0, 0, 14,
             3, 8, 5, 9, 14,
             1, 0, 0, 0, 17)
  stats <- matrix(stats, nrow = 3, byrow = TRUE)

  expect_type(pam, "list")
  expect_s3_class(pam, "PAMclusts")
  expect_length(pam, 5)
  expect_identical(unname(pam$distM), distM)
  expect_s3_class(pam$phyloTree, "phylo")
  expect_equal(unname(pam$clustering), clustering)
  expect_identical(pam$medoids, medoids)
  expect_identical(unname(pam$stats), stats)
})

test_that("EM clustering with 3 clusters from a character string of newick tree", {

  em <- clustEM(k = 3, text = NwkTree)

  distM <- c(0, 16, 14, 15, 17,
             16, 0, 8, 9, 21,
             14, 8, 0, 7, 19,
             15, 9, 7, 0, 20,
             17, 21, 19, 20, 0)
  distM <- matrix(distM, nrow = 5, byrow = TRUE)
  clustering <- c(1, 3, 3, 3, 2)
  trunc.mean <- c(-3, -13, 5,
                  -9, 4, 1,
                  0, 0, 0,
                  0, 0, 0)
  trunc.mean <- matrix(trunc.mean, nrow = 4, byrow = TRUE)

  expect_type(em, "list")
  expect_s3_class(em, "EMclusts")
  expect_length(em, 7)
  expect_identical(unname(em$distM), distM)
  expect_s3_class(em$phyloTree, "phylo")
  expect_equal(unname(em$clustering), clustering)
  expect_identical(trunc(em$mean), trunc.mean)
  expect_identical(trunc(em$bic), -113)
  expect_identical(em$model, "spherical, equal volume")
})

test_that("clustPAM error upon invalid user input", {

  # k provided as a double
  expect_error(em <- clustPAM(k = 2.5, text = NwkTree))

  # k provided as a character
  expect_error(pam <- clustPAM(k = "3", text = NwkTree))

  # k is 1
  expect_error(pam <- clustPAM(k = 1, text = NwkTree))

  # k is greater than the number of leaves in the tree
  expect_error(pam <- clustPAM(k = 8, text = NwkTree))

  # invalid file path
  expect_error(pam <- clustPAM(k = 3, file = "invalid/path"))

  # text tree provided as a numeric
  expect_error(pam <- clustPAM(k = 3, text = 5))

  # multiple trees exist in input text or file
  expect_error(pam <- clustPAM(k = 3,
                               text = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);
                                 (B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"))

  # no branch lengths in input tree
  expect_error(pam <- clustPAM(k = 3, text = "(B,(A,C,E),D);"))

  # invalid tree format
  expect_error(pam <- clustPAM(k = 3, text = "('B(bad)':6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"))

  # tree contains redundant leaf names
  expect_error(pam <- clustPAM(k = 3, text = "(B:6.0,(B:5.0,C:3.0,E:4.0):5.0,D:11.0);"))
})

test_that("clustEM error upon invalid user input", {

  # k provided as a double
  expect_error(em <- clustEM(k = 2.5, text = NwkTree))

  # k provided as a character
  expect_error(em <- clustEM(k = "3", text = NwkTree))

  # k is 1
  expect_error(em <- clustEM(k = 1, text = NwkTree))

  # k is greater than the number of leaves in the tree
  expect_error(em <- clustEM(k = 8, text = NwkTree))

  # invalid file path
  expect_error(em <- clustEM(k = 3, file = "invalid/path"))

  # text tree provided as a numeric
  expect_error(em <- clustEM(k = 3, text = 5))

  # multiple trees exist in input text or file / bad format
  expect_error(em <- clustEM(k = 3,
                             text = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);
                                 (B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"))

  # no branch lengths in input tree
  expect_error(em <- clustEM(k = 3, text = "(B,(A,C,E),D);"))

  # invalid tree format
  expect_error(em <- clustEM(k = 3, text = "('B(bad)':6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"))

  # tree contains redundant leaf names
  expect_error(em <- clustEM(k = 3, text = "(B:6.0,(B:5.0,C:3.0,E:4.0):5.0,D:11.0);"))
})
