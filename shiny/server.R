library(shiny)

# server.R
library(parallel)
library(quantmod)
library(ghyp)

N=252
source("helpers.R")



shinyServer(function(input, output) {
  dataInput <- reactive({
    myts=getSymbols(input$symb, src = "yahoo", 
                    from = input$dates[1],
                    to = input$dates[2],
                    auto.assign = FALSE)
    #print(head(myts))
    mySNR=estimateSharpeRatio(diff(log(myts[,6])))
    
    SNR0s=unlist(lapply(mySNR,function(x) x$SNR0))
    SNRs=unlist(lapply(mySNR,function(x) max(-1,min(1,x$SNR))))
    nus=unlist(lapply(mySNR,function(x) x$nu))
    R0bars=unlist(lapply(mySNR,function(x) x$R0bar))
    

    SNRs=xts(SNRs,as.Date(names(SNRs)))*sqrt(252)
    SNR0s=xts(SNR0s,as.Date(names(SNR0s)))*sqrt(252)
    nus=xts(nus,as.Date(names(nus)))
    R0bars=xts(R0bars,as.Date(names(R0bars)))
 
    zeroxts=xts(rep(NA,N),seq(start(SNRs)-N,start(SNRs)-1,length=N))
    SNRs=c(zeroxts,SNRs)
    SNR0s=c(zeroxts,SNR0s)
    nus=c(zeroxts,nus)
    R0bars=c(zeroxts,R0bars)
    return(list(myts=myts[,6],SNRs=SNRs,SNR0s=SNR0s,nus=nus,R0bars=R0bars))
  })
  
  output$plot <- renderPlot({
    myres=dataInput()
    plot(myres$myts[paste0(start(myres$SNRs),"::",end(myres$SNRs))],main=input$symb,log='y',ylab="Adjusted price")
  })
 
  output$plot2 <- renderPlot({
    myres=dataInput()
    plot(myres$SNR0s,main="Sharpe ratio estimation",ylab="Annualized Sharpe ratio")
    lines(myres$SNR0s,lwd=3)
    lines(myres$SNRs,col=2,lwd=3)
    legend("topright",legend=c("vanilla estimator","moment-free estimator"),col=1:2,lwd=2)
  })
  
   
})