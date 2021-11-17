library(ClustPhy)
library(ape)

test_that("plotGapStat error upon invalid user input", {

  tree <- ape::read.tree(text = NwkTree2)
  distM <- ape::cophenetic.phylo(tree)
  gapStat <- compareGap(distM, k.max = 10)

  # gapStat is not of the class "clusGap"
  expect_error(plotGapStat(gapStat = "not a S3 clusGap object"))

  # method is not one of "globalmax", "firstmax", "Tibs2001SEmax", "firstSEmax"
  expect_error(plotGapStat(gapStat = gapStat,
                           method = "not a valid method"))

  # color is not a character string
  expect_error(plotGapStat(gapStat = gapStat,
                           color = 110000))
})
