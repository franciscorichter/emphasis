library(shiny)
library(paleobuddy)
library(ape)
library(shinyjs)


if (!requireNamespace("TreeSim", quietly = TRUE)) {
  install.packages("TreeSim")
}

library(TreeSim)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Tree Simulation App"),
  sidebarLayout(
    sidebarPanel(
      selectInput("sim_type", "Select type of simulation:", 
                  choices = c("Simulate one tree" = "one_tree",
                              "Simulate multiple trees" = "multiple_trees")
      ),
      selectInput("sim_function", "Select Model:", 
                  choices = c("Time-dependent Birth Rate" = "time_dependent_rate",
                              "Constant Birth Rate" = "constant_rate")),
      textOutput("description"),
      uiOutput("paramInputs"),
      actionButton("simulate", "Simulate")
    ),
    mainPanel(
      tags$h4("Simulation Results"),
      textOutput("simulatedMessage"),
      textOutput("computeTime"),
      textOutput("extantTips"),
      textOutput("fullTreeTips"),
      downloadButton('downloadData', 'Download Phylo Object', style = "display: none;"),
      tags$h4("Visualizations"),
      checkboxInput("showTree", "Show Tree Plot", FALSE),
      plotOutput("treePlot")
    )
  )
)

server <- function(input, output, session) {
  
  simulated <- reactiveVal(FALSE)
  
  observe({
    if (simulated()) {
      shinyjs::show("downloadData")
    } else {
      shinyjs::hide("downloadData")
    }
  })
  
  
  observeEvent(input$simulate, {
    start_time <- Sys.time()
    
    if(input$sim_function == "time_dependent_rate") {
      lambda <- function(t) {
        0.1 + 0.001*t
      }
      mu <- function(t) {
        0.03 * exp(-0.01*t)
      }
      tree <- make.phylo(bd.sim(n0 = input$n0,
                                lambda = input$lambda, 
                                mu = input$mu, 
                                tMax = input$tMax))
      
    } else if(input$sim_function == "constant_rate") {
      tree <- sim.bd.taxa(input$taxa, input$birth_rate, input$death_rate, input$time)[[1]]
    } 
    
    end_time <- Sys.time()
    time_taken <- end_time - start_time
    
    tree_data(tree)  # Save the tree to the reactive value
    
    output$simulatedMessage <- renderText({
      "Simulated Tree"
    })
    
    output$computeTime <- renderText({
      paste("Computation Time:", round(as.numeric(time_taken), 2), "seconds")
    })
    
    output$extantTips <- renderText({
      tree <- tree_data()
      paste("Number of tips in extant tree:", length(tree$tip.label))
    })
    
    output$fullTreeTips <- renderText({
      tree <- tree_data()
      paste("Number of tips in full tree:", length(tree$tip.label) + length(tree$node.label))
    })
    
    simulated(TRUE)
  })
  
  output$simulated <- reactive({
    simulated()
  })
  
  output$paramInputs <- renderUI({
    wellPanel(
      tags$h4("Parameters:"), 
      if (input$sim_function == "time_dependent_rate") {
        list(
          numericInput("n0", "Initial number of species:", 1),
          numericInput("lambda", "lambda:", 50),
          numericInput("mu", "Mu:", 50),
          numericInput("tMax", "Maximum simulation time:", 50)
        )
      } else if (input$sim_function == "constant_rate") {
        list(
          numericInput("taxa", "Number of taxa:", 10),
          numericInput("birth_rate", "Birth rate:", 1),
          numericInput("death_rate", "Death rate:", 0.5),
          numericInput("time", "Simulation time:", 0.1)
        )
      }
    )
  })
  
  output$treePlot <- renderPlot({
    if (input$showTree) {
      tree <- tree_data()
      plot(tree, show.tip.label = FALSE)
      axisPhylo()
    }
  })
  
  output$description <- renderText({
    if (input$sim_function == "time_dependent_rate") {
      return("This simulation model uses a time-dependent birth rate. Parameters like initial number of species, lambda, mu, and maximum simulation time influence the resulting tree.")
    } else if (input$sim_function == "constant_rate") {
      return("This simulation model uses a constant birth rate. Parameters like number of taxa, birth rate, death rate, and simulation time determine the shape and structure of the tree.")
    } else {
      return("")
    }
  })
  
  # Placeholder for the phylo object
  tree_data <- reactiveVal()
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("phylo_object", Sys.Date(), ".RData", sep="")
    },
    content = function(file) {
      temp_tree <- tree_data()
      save(temp_tree, file = file)
    },
    contentType = "application/x-rdata"
  )
  
}

shinyApp(ui = ui, server = server)

