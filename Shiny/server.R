library(shiny)
library(shinydashboard)
library(plotly)
library(readr)
library(factoextra)

server <- function(input, output) {
  RNA <- read.csv("../data/GSE74534_RNA-seq_raw_counts.csv")
  colnames(RNA)[1] <- "Gene ID"
  scMT <- readRDS('../data/GSE74534_DNAm.rds')
  A02 <- as.data.frame(read_delim("../data/A02.cov",delim="\t",col_names = F))
  colnames(A02) <- c("Chromosome of gene", "Start position", "End position ",
                     "Methylation percentage", "Count methylated", "count non-methylated")
  id_vec <- unique(scMT$ens_id)
  RNA_sub <- filter(RNA, RNA[,1]%in%id_vec)
  rownames(RNA_sub) <- RNA_sub[,1]
  RNA_sub <- RNA_sub[,-1]
  RNA <- RNA_sub
  m <- median(colSums(RNA))
  tmp <- m/colSums(RNA)
  tmp_mat <- t(matrix(tmp, dim(RNA)[2], dim(RNA)[1]))
  RNA_norm <- floor(RNA*tmp)
  RNA_sd <- apply(RNA_norm, 1, sd)
  RNA_cbind <- cbind(RNA,RNA_sd)
  RNA_cbind_name <- cbind(rownames(RNA_cbind),RNA_cbind)
  colnames(RNA_cbind_name)[1] <- "Gene_id" 
  RNA_500 <- RNA_cbind_name %>% arrange(desc(RNA_sd)) %>% head(500)
  RNA_500 <- RNA_500[,-63]
  RNA_500 <- RNA_500[,-1]
  
  output$scRNA <- DT::renderDataTable({
    dt <- DT::datatable(RNA)
  })
  
  output$summaryR <- renderPrint({
    summary(RNA)
  })
  
  output$DNAm <- DT::renderDataTable({
    dt <- DT::datatable(A02)
  })
  
  output$summaryM <- renderPrint({
    summary(A02)
  })
  
  output$Plotheat <- renderPlotly({
    vals <- unique(scales::rescale(c(volcano)))
    o <- order(vals, decreasing = FALSE)
    cols <- scales::col_numeric("Blues", domain = NULL)(vals)
    colz <- setNames(data.frame(vals[o], cols[o]), NULL)
    p <- plot_ly(z = as.matrix(RNA_500),
                 zmax=10000,zmin=0,
                 colorscale = colz, type = "heatmap")%>%
      layout(
        xaxis = list(title = "Cell")
        , yaxis = list(title = "Gene"))
  })
  
  output$Plotmsd <- renderPlotly({
    mean_sd <- as.data.frame(cbind(colSums(RNA_norm), apply(RNA_norm,2,sd)))
    colnames(mean_sd) <- c("mean","sd")
    p<-plot_ly(data = mean_sd ,y = ~sd, x = ~mean, 
               color = ~sd,mode = 'markers',showlegend=T)%>%
      layout(
        xaxis = list(title = "mean")
        , yaxis = list(title = "sd"))
  })
  
  output$Plotcol <- renderPlotly({
    zero_ratio <- colSums(RNA_norm == 0)/dim(RNA_norm)[1]
    rate_count <- as.data.frame(cbind(zero_ratio, colSums(RNA_norm)))
    p<-plot_ly(data = rate_count ,y = ~rate_count[,1], x = ~rate_count[,2], 
               color = ~zero_ratio,mode = 'markers',showlegend=T)%>%
      layout(xaxis = list(title = "total read counts of each single cell"),
             yaxis = list(title = "zero ratio"))
  })
  
  output$Plotrow <- renderPlotly({
    zero_ratio_r <- rowSums(RNA_norm == 0)/dim(RNA_norm)[2]
    zero_ratio_r <- as.data.frame(cbind(zero_ratio_r, rowSums(RNA_norm)))
    p<-plot_ly(data = zero_ratio_r ,y = ~zero_ratio_r, x = ~V2, 
               color = ~zero_ratio_r,mode = 'markers',showlegend=T)%>%
      layout(xaxis = list(title = "total read counts of each single gene"),
             yaxis = list(title = "zero ratio"))
  })
  
  output$PlotBIC <- renderPlotly({
    BICvalue <- c(379599.7, 392696.1, 386628.4, 393845)
    Kvalue <- c(1,2,3,4)
    BICdata <- as.data.frame(cbind(Kvalue, BICvalue))
    K <- input$K
    BICinp <- BICdata[1:K,]
    p <- plot_ly(data = BICinp, x = ~as.factor(Kvalue), y = ~BICvalue, mode = 'lines+markers',type = 'scatter')%>%
      layout(
        xaxis = list(title = "K"),
        yaxis = list(title = "BIC"))
  })
  
  output$PlotInd <- renderPlotly({
    rdata <- paste("../data/RNA_DNAm_",input$clu,".RData", sep = "")
    load(rdata)
    d1 <- c(sum(Gam_est), 500-sum(Gam_est))
    d2 <- c(1,0)
    Gam <- cbind(d1,d2)
    p1<-plot_ly(data = as.data.frame(Gam), y = ~d1, x= ~as.factor(d2),
                colorscale = "Portland", type = 'bar', name = "Gamma"
    )%>%
      layout(
        xaxis = list(title = "Gamma indicator"),
        yaxis = list(title = ""))
    d1 <- c(sum(Omega_est), 89-sum(Omega_est))
    d2 <- c(1,0)
    Omega <- cbind(d1,d2)
    p2<-plot_ly(data = as.data.frame(Omega), y = ~d1, x= ~as.factor(d2),
                colorscale = "Portland", type = 'bar', name = "Omega"
    )%>%
      layout(
        xaxis = list(title = "Omega indicator"),
        yaxis = list(title = ""))
    p<-subplot(p1,p2)%>%
      layout(
        xaxis = list(title = "Indicator"),
        xaxis2 = list(title = "Indicator"))
  })
  
  output$Plotres <- renderPlotly({
    rdata <- paste("../data/RNA_DNAm_",input$clu,".RData", sep = "")
    load(rdata)
    set.seed(20190330)
    km_result <- kmeans(Z_est, input$clu)
    cell <- as.data.frame(cbind(c(1:61), Z_est[,2],km_result$cluster))
    p<-plot_ly(data = cell, x = ~cell[,1], y = ~cell[,2], symbol =as.factor(cell$V3),
               colorscale = "Portland",
               mode = 'markers',type = 'scatter',showlegend=T
    )%>%
      layout(
        xaxis = list(title = "Cell ID"),
        yaxis = list(title = "Z value"))
  })
}

# # Run the application 
# shinyApp(ui = ui, server = server)

