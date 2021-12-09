#' Launch Shiny App for ClustPhy
#'
#' A function that launches the Shiny app for ClustPhy.
#' The shiny app allows user to perform phylogenetic clustering via an
#' interactive platform.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#' ClustPhy::runClustPhy()
#' }
#'
#' @author {Yuzi Li, \email{rainal.li@mail.utoronto.ca}}
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runClustPhy <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "ClustPhy")
  shiny::runApp(appDir, display.mode = "normal")
  return()
}

# [END]
