library(shiny)
if (! require(colourpicker)) {
  install.packages("colourpicker")
}
library(colourpicker)

ui <- fluidPage(

  # Add title panel
  titlePanel(tags$h1(tags$b("ClustPhy"), ": Phylogenetic Clustering", tags$hr()),
             windowTitle = "ClustPhy"),

  # Use sidebar layout
  sidebarLayout(

    # Add sidebar panel
    sidebarPanel(

      # Ask for text input of newick tree
      textAreaInput(inputId = "tree",
                    label = "Enter Newick tree text",
                    value = NwkTree2,
                    height = "400px",
                    resize = "both"),

      tags$hr(),

      # Ask for number of desired clusters
      numericInput(inputId = "k",
                   label = "Number of clusters",
                   value = 3,
                   min = 2),

      tags$hr(),

      # Select clustering method. Default: PAM
      selectInput(inputId = "clustMethod",
                  label = "Clustering method",
                  choices = c("k-medoids (PAM)" = "PAM",
                              "Expectation maximization (EM)" = "EM")),

      # EM-specific options
      conditionalPanel(

        condition = "input.clustMethod == 'EM'",

        style = "padding-left: 50px;",

        tags$hr(),

        tags$h4("EM-specific options: "),

        # Adjust maximum number of principle components used for EM
        numericInput(inputId = "maxPC",
                     label = "Enter maxPC",
                     value = 5,
                     min = 2),

        # Show the percentage of total variance explained by the selected
        # number of principle components
        textOutput("pvar"),

        tags$br(),

        tags$p("* maxPC indicates the maximum number of dimensions of",
               tags$br(),
               "the reduced tree leaves coordinates after PCA that is used towards",
               tags$br(),
               "EM clustering. Usually most of the variance in the data can be",
               tags$br(),
               "explained by the top 5 principle components. Including too many",
               tags$br(),
               "dimensions can lead to sparse datapoints and prevent effective",
               tags$br(),
               "clustering.")
      ),

    ),
    mainPanel(

      # Build output tabs
      tabsetPanel(

        # Tab for clusters phylogram
        tabPanel(
          title = "Clusters Phylogram",

          # Wrap splitLayout (float left) inside fluidRow to define a top padding
          fluidRow(

            # Add padding
            style = "padding-top: 50px; padding-left: 50px",

            # Add biplot options
            splitLayout(

              # Define style
              tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                                float: left;
                              }
                              "))),

              # Adjust size of nodes
              sliderInput(inputId = "phyloNodeSize",
                          label = "Adjust node size",
                          value = 1,
                          min = 0.1,
                          max = 10),

              # Customize title
              textInput(inputId = "phyloTitle",
                        label = "Enter title",
                        value = "Tree Clusters Phylogram")
            )
          ),

          # Add PAM-specific options
          conditionalPanel(

            condition = "input.clustMethod == 'PAM'",

            style = "padding-left: 50px;",

            tags$hr(),

            tags$h4("PAM-specific options: "),

            splitLayout(

              # Whether to show medoids on plot
              checkboxInput(inputId = "phyloShowMedoids",
                            label = "Show medoids",
                            value = FALSE),

              # If show medoids, ask to show custom symbol or medoid names
              conditionalPanel(

                condition = "input.phyloShowMedoids",

                verticalLayout(

                  # Show either custom symbol or medoid names
                  selectInput(inputId = "phyloSymbol",
                              label = "Choose labels for medoids",
                              choices = c("Medoid names" = "names",
                                          "Custom symbol" = "custom")),

                  # Control symbol size
                  sliderInput(inputId = "phyloSymbolCex",
                              label = "Adjust label size",
                              value = 1,
                              min = 0.1,
                              max = 10)
                )
              ),

              # If use custom symbol, ask user to enter custom symbol
              conditionalPanel(

                condition = "input.phyloSymbol == 'custom'",

                textInput(inputId = "phyloCustomSymbol",
                          label = "Enter custom label for medoids",
                          value = " * ")
              )
            )
          ),

          plotOutput("phylogram", height = "600px")
        ),

        # Tab for clusters biplot
        tabPanel(
          title = "Clusters Biplot",

          # Wrap splitLayout (float left) inside fluidRow to define a top padding
          fluidRow(

            # Add padding
            style = "padding-top: 50px; padding-left: 50px",

            # Add biplot options
            splitLayout(

              # Define style
              tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                                float: left;
                              }
                              "))),

              # Adjust size of nodes
              sliderInput(inputId = "biplotNodeSize",
                          label = "Adjust node size",
                          value = 1,
                          min = 0.1,
                          max = 10),

              # Customize title
              textInput(inputId = "biplotTitle",
                        label = "Enter title",
                        value = "Tree Clusters Biplot")
            )
          ),

          # Add PAM-specific options
          conditionalPanel(

            condition = "input.clustMethod == 'PAM'",

            style = "padding-left: 50px;",

            tags$hr(),

            tags$h4("PAM-specific options: "),

            splitLayout(

              # Whether to show medoids on plot
              checkboxInput(inputId = "biplotShowMedoids",
                            label = "Show medoids",
                            value = FALSE),

              # If show medoids, ask to show custom symbol or medoid names
              conditionalPanel(

                condition = "input.biplotShowMedoids",

                verticalLayout(

                  # Show either custom symbol or medoid names
                  selectInput(inputId = "biplotSymbol",
                              label = "Choose labels for medoids",
                              choices = c("Medoid names" = "names",
                                          "Custom symbol" = "custom")),

                  # Control symbol size
                  sliderInput(inputId = "biplotSymbolCex",
                              label = "Adjust label size",
                              value = 1,
                              min = 0.1,
                              max = 10)
                )
              ),

              # If use custom symbol, ask user to enter custom symbol
              conditionalPanel(

                condition = "input.biplotSymbol == 'custom'",

                textInput(inputId = "biplotCustomSymbol",
                          label = "Enter custom label for medoids",
                          value = " * ")
              )
            )
          ),

          # Plot biplot
          plotOutput("biplot", height = "600px")
        ),

        # Tab for finding the best k number of clusters via gap statistics
        tabPanel(
          title = "Gap Statistics",

          # Wrap splitLayout (float left) inside fluidRow to define a top padding
          fluidRow(

            # Add top padding
            style = "padding-top: 50px; padding-left: 50px",

            # Add gapstat options
            splitLayout(

              # Define style
              tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                                float: left;
                              }
                              "))),

              # Ask for k.max
              numericInput(inputId = "kmax",
                           label = "Maximum number of clusters",
                           value = 10,
                           min = 2),

              # Pick a color
              colourInput(inputId = "gapstatCol",
                          label = "Select color",
                          value = "steelblue"),

              # Select method for computing optimal k value
              selectInput(inputId = "gapstatMethod",
                          label = "Method for determining optimal k value",
                          choices = c("Global maximum" = "globalmax",
                                      "First maximum" = "firstmax",
                                      "Tibshirani et al (2001) SE maximum" = "Tibs2001SEmax",
                                      "First SE maximum" = "firstSEmax"),
                          selected = "Tibs2001SEmax")
            )
          ),
          # Plot gapstat
          plotOutput("gapstat", height = "600px")
        )
      )
    )
  )
)

server <- function(input, output) {

  # Define clustering based on user selection
  clusts <- reactive({
    if (input$clustMethod == "PAM") {
      clustPAM(k = input$k, text = input$tree)
    } else {
      clustEM(k = input$k, text = input$tree, maxPC = input$maxPC)
    }
  })

  # Add cluster plots to output
  # Output phylogram
  output$phylogram <- renderPlot({
    if (class(clusts()) == "PAMclusts" & input$phyloShowMedoids) {
      # Plot showing medoids
      labels <- clusts()$medoids
      if (input$phyloSymbol == "custom") {
        labels <- input$phyloCustomSymbol
      }
      plotClustersTree(tree = clusts()$phyloTree,
                     clustering = clusts()$clustering,
                     node.cex = input$phyloNodeSize,
                     title = input$phyloTitle,
                     show.centers = clusts()$medoids,
                     center.symbol = labels,
                     symbol.cex = input$phyloSymbolCex)
    } else {
      # Plot without medoids
      plotClustersTree(tree = clusts()$phyloTree,
                     clustering = clusts()$clustering,
                     node.cex = input$phyloNodeSize,
                     title = input$phyloTitle)
    }
  })

  # Output biplot
  output$biplot <- renderPlot({
    if (class(clusts()) == "PAMclusts" & input$biplotShowMedoids) {
      # Plot showing medoids
      labels <- clusts()$medoids
      if (input$biplotSymbol == "custom") {
        labels <- input$biplotCustomSymbol
      }
      plotClusters2D(tree = clusts()$phyloTree,
                     clustering = clusts()$clustering,
                     node.cex = input$biplotNodeSize,
                     title = input$biplotTitle,
                     show.centers = clusts()$medoids,
                     center.symbol = labels,
                     symbol.cex = input$biplotSymbolCex)
    } else {
      # Plot without medoids
      plotClusters2D(tree = clusts()$phyloTree,
                     clustering = clusts()$clustering,
                     node.cex = input$biplotNodeSize,
                     title = input$biplotTitle)
    }
  })

  # Compute percent variance accounted by PCs when EM clustering is used
  output$pvar <- renderText({
    if (class(clusts()) == "EMclusts") {
      pca <- clusts()$dimredResult$PCA
      maxPC <- input$maxPC
      if (input$maxPC > ncol(pca$rotation)) {
        maxPC <- ncol(pca$rotation)
      }
      pvars <- summary(pca)$importance["Proportion of Variance", 1:maxPC]
      sprintf("%#.1f%% of total variance is explained by %d dimensions. ",
              sum(pvars) * 100,
              maxPC)
    } else {
      NULL
    }
  })

  # Add gap statistics plot to output
  output$gapstat <- renderPlot({
    # Calculate gap statistics
    gapstat <- compareGap(distM = clusts()$distM,
                          k.max = input$kmax,
                          method = input$clustMethod)
    # Plot gap statistics and optimal k
    plotGapStat(gapStat = gapstat,
                method = input$gapstatMethod,
                color = input$gapstatCol)
  })
}
shiny::shinyApp(ui = ui, server = server)
# [END]
