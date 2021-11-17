library(ClustPhy)

test_that("compareGap works using pre-defined clusterings", {

  pam2 <- clustPAM(k = 2, text = NwkTree)
  pam3 <- clustPAM(k = 3, text = NwkTree)
  pam4 <- clustPAM(k = 4, text = NwkTree)
  clusterings <- list(pam2$clustering, pam3$clustering, pam4$clustering)
  distM <- pam2$distM

  set.seed(5)
  gapstat <- compareGap(distM = distM,
                        k.max= 4,
                        clusterings = clusterings)

  tab <- c(3.1, 2.7, -0.4, 0.2,
           2.6, 2.4, -0.1, 0.2,
           1.8, 2.0, 0.3, 0.3,
           0.9, 1.3, 0.3, 0.4)
  tab <- matrix(tab, nrow = 4, byrow = TRUE)

  expect_type(gapstat, "list")
  expect_s3_class(gapstat, "clusGap")
  expect_identical(unname(round(gapstat$Tab, digits = 1)), tab)

  set.seed(NULL)
})

test_that("compareGap works for EM method without pre-defined clusterings", {

  em <- clustEM(k = 2, text = NwkTree)
  distM <- em$distM

  set.seed(5)
  gapstat <- compareGap(distM = distM,
                        k.max= 4,
                        method = "EM")

  tab <- c(3.1, 2.7, -0.4, 0.2,
           2.5, 2.5, 0.0, 0.2,
           1.8, 2.0, 0.3, 0.3,
           0.9, 1.3, 0.3, 0.4)
  tab <- matrix(tab, nrow = 4, byrow = TRUE)

  expect_type(gapstat, "list")
  expect_s3_class(gapstat, "clusGap")
  expect_identical(unname(round(gapstat$Tab, digits = 1)), tab)

  set.seed(NULL)
})

test_that("compareGap works for PAM method without pre-defined clusterings", {

  pam <- clustPAM(k = 2, text = NwkTree)
  distM <- pam$distM

  set.seed(5)
  gapstat <- compareGap(distM = distM,
                        k.max= 4,
                        method = "PAM")

  tab <- c(3.1, 2.7, -0.4, 0.2,
           2.6, 2.4, -0.1, 0.2,
           1.8, 2.0, 0.3, 0.3,
           0.9, 1.3, 0.3, 0.4)
  tab <- matrix(tab, nrow = 4, byrow = TRUE)

  expect_type(gapstat, "list")
  expect_s3_class(gapstat, "clusGap")
  expect_identical(unname(round(gapstat$Tab, digits = 1)), tab)

  set.seed(NULL)
})

test_that("compareGap error upon invalid user input", {

  pam2 <- clustPAM(k = 2, text = NwkTree)
  pam3 <- clustPAM(k = 3, text = NwkTree)
  pam4 <- clustPAM(k = 4, text = NwkTree)
  distM <- pam2$distM
  clusterings <- list(pam2$clustering, pam3$clustering, pam4$clustering)

  # distM is not a matrix
  expect_error(compareGap(distM = c(1), k.max = 2))

  # distM is not numeric
  expect_error(compareGap(distM = matrix("not numeric"), k.max = 2))

  # distM is not a square matrix
  expect_error(compareGap(distM = distM[1:2, ], k.max = 2))

  # distM contains NA
  bad.matrix <- distM
  bad.matrix[1, 1] <- NA
  expect_error(compareGap(distM = bad.matrix, k.max = 2))

  # distM is not symmetric
  bad.matrix <- distM
  bad.matrix[1, 3] <- 100
  expect_error(compareGap(distM = bad.matrix, k.max = 2))

  # k.max is less than 2
  expect_error(compareGap(distM = distM, k.max = 1))

  # k.max is more than the number of leaves
  expect_error(compareGap(distM = distM, k.max = 10))

  # clusterings is not a list
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          clusterings = c(pam2$clustering, pam3$clustering)))

  # clusterings length does not equal to k.max - 1
  expect_error(compareGap(distM = distM,
                          k.max = 5,
                          clusterings = clusterings))

  # clusterings contain non-integer values
  bad.clustering <- pam2$clustering
  bad.clustering <- replace(bad.clustering, bad.clustering == 1, 1.5)
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          clusterings = list(bad.clustering, pam3$clustering)))

  # number of clusters in clustering does not match the k range of 2 to k.max
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          clusterings = list(pam2$clustering, pam4$clustering)))

  # clusterings not ordered by increasing number of clusters
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          clusterings = list(pam3$clustering, pam2$clustering)))

  # clustering length does not match number of leaves
  bad.clustering <- pam2$clustering[1:3]
  expect_error(compareGap(distM = distM,
                         k.max = 3,
                         clusterings = list(bad.clustering, pam3$clustering)))

  # method is not one of "PAM" or "EM"
  expect_error(compareGap(distM = distM,
                         k.max = 3,
                         method = "not a valid method"))

  # B is not a numeric
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          B = "not a number"))

  # B is less than 1
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          B = 0))

  # B is not an integer
  expect_error(compareGap(distM = distM,
                          k.max = 3,
                          B = 100.5))
})
