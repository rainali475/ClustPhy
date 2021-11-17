#' Plot gap statistics.
#'
#' A function that plots given gap statistics and finds the best k value
#' using a user-specified method.
#'
#' @param gapStat An object of class "clusGap" indicating the gap statistics to
#'     be plotted. This is usually output from comparePAM, compareEM, or
#'     cluster::clusGap.
#' @param method A character string indicating the method used to determine
#'     best k. Default value is "Tibs2001SEmax". Must be one of "globalmax",
#'     "firstmax", "Tibs2001SEmax", "firstSEmax". Tibs2001SEmax is the
#'     method proposed by Tibshirani et al (2001): "the smallest k such that
#'     gap(k) >= gap(k+1) - s_k+1".
#' @param color A character string indicating the color of line graph. Default
#'     value is "steelblue".
#'
#' @return A plot of the gap statistics, with a vertical dashed line indicating
#'     the best k value.
#'
#' @examples
#' pam <- clustPAM(6, text = NwkTree2)
#' set.seed(5)
#' gapStat <- compareGap(distM = pam$distM, k.max = 10, method = "PAM")
#' gapStat$Tab
#' set.seed(NULL)
#' plotGapStat(gapStat)
#'
#' @export
#' @import factoextra
plotGapStat <- function(gapStat, method = "Tibs2001SEmax", color = "steelblue") {

  # Check user input
  if (class(gapStat) != "clusGap") {
    stop("gapStat must be an object of class clusGap.")
  }

  if (! is.character(method) | ! method %in% c("globalmax",
                                              "firstmax",
                                              "Tibs2001SEmax",
                                              "firstSEmax")) {
    stop("method should be a character string and should be one of
         \"globalmax\", \"firstmax\", \"Tibs2001SEmax\", \"firstSEmax\".")
  }

  if (! is.character(color)) {
    stop("color must be a character string indicating the color of line graph.")
  }

  # Generate gap statistics pot using factoextra::fviz_gap_stat
  return(factoextra::fviz_gap_stat(gapStat,
                                   maxSE = list(method = method, SE.factor = 1),
                                   linecolor = color))
}
