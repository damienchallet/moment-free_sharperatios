library(shiny)


shinyUI(fluidPage(
  titlePanel("Moment-free estimation of Sharpe ratio"),
  
 sidebarLayout(
    sidebarPanel(
      helpText("Select a stock to examine. 
        Information will be collected from Yahoo finance."),
    
      textInput("symb", "Symbol", "SHP.L"),
    
      dateRangeInput("dates", 
        "Date range",
        start = "2014-01-01", 
        end = as.character(Sys.Date())),
      
      br(),
      br(),
      
      helpText("This app provides an interactive way to estimate Sharpe ratios on rolling windows of 252 days using the moment-free methods introduced in "),
               tags$a(href="http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2603682","D. Challet, Moment-free Sharpe ratio estimation with total drawdown durations, 2015"),
      br(),
      br(),
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