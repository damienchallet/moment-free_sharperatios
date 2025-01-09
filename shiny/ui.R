library(shiny)

shinyUI(fluidPage(
  titlePanel("Moment-free estimation of Sharpe ratios"),
  
 sidebarLayout(
    sidebarPanel(
      helpText("Choose a symbol. 
        Daily data will be downloaded from Yahoo Finance. Computation time is about 10 seconds per year."),
    
      textInput("symb", "Symbol", "TSLA"),
    
      dateRangeInput("dates", 
        "Date range",
        start = "2019-01-01", 
        end = as.character(Sys.Date())),
      numericInput("Tin","calibration length",value = 50),
      actionButton("startButton","Start"),
      br(),
      br(),
      
      helpText("This app provides an interactive way to estimate Sharpe ratios on rolling windows of 252 days using the moment-free methods introduced in "),
               tags$a(href="http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2603682","D. Challet, Sharper Asset Ranking from Total Drawdown Durations, 2015"),

      helpText("The moment-free estimator is more robust, precise (efficient), and does not over-estimate the amplitude of Sharpe ratios in wild market conditions."),
      br(),
      helpText("An R package is available on CRAN"),
      tags$a(href="https://cran.r-project.org/web/packages/sharpeRratio/index.html","sharpeRratio"),
      h5("Author: Damien Challet. Some parts of the code have been taken from Shiny's tutorials"),
      a(href="http://shiny.rstudio.com/tutorial/","http://shiny.rstudio.com/tutorial/"),
      br(),
      h5("Full source code is available at"),
      a(href="https://github.com/damienchallet/moment-free_sharperatios","https://github.com/damienchallet/moment-free_sharperatios")
    ),
    
    mainPanel(plotOutput("plot"), 
              plotOutput("plot2"))
  ))
)