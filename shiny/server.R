library(shiny)

# server.R
library(parallel)
library(quantmod)
library(sharpeRratio)

source("helpers.R")

N=252

shinyServer(function(input, output) {
  dataInput <-  eventReactive(input$startButton,{
    if(as.numeric(input$dates[2])-as.numeric(input$dates[1])<N+20){
      input$dates[1]=input$dates[2]-floor(N*365.25/252)-20  #one needs at least as many points as N
      
    }
    print("Downloading data ...")
    myts=getSymbols(input$symb, src = "yahoo", 
                    from = input$dates[1],
                    to = input$dates[2],
                    auto.assign = FALSE,)
    #print(head(myts))
    print("Computing ...")
    price=myts[,6]
    names(price)="price"
    returns=    na.omit(diff(log(price)))
    SNR_DT=estimateSharpeRatio(returns,N=input$Tin)
    SNR_xts=as.xts(SNR_DT)
    m=merge(SNR_xts,price,all=TRUE)
    return(m)
  })
  
  output$plot <- renderPlot({
    myres=dataInput()
    plot(myres$price[paste0(start(myres$SNR),"::",end(myres$SNR))],main=input$symb,log='y',ylab="Adjusted price")
  })
 
  output$plot2 <- renderPlot({
    myres=dataInput()
    plot(myres$SNR0,main="Sharpe ratio estimation",ylab="Annualized Sharpe ratio")
    lines(myres$SNR,col=2,lwd=3)
    legend("topright",legend=c("vanilla estimator","moment-free estimator"),col=1:2,lwd=2)
  })
  
   
})