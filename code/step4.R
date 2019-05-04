setwd("../data")

RNA <- read.csv("GSE74534_RNA-seq_raw_counts.csv")
colnames(RNA)[1] <- "Gene ID"
scMT <- readRDS('GSE74534_DNAm.rds')
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

## Figure 1: Heatmap of ScRNA-seq data
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


## Figure 2: Means and standard deviations of ESCs
mean_sd <- as.data.frame(cbind(colSums(RNA_norm), apply(RNA_norm,2,sd)))
colnames(mean_sd) <- c("mean","sd")
p<-plot_ly(data = mean_sd ,y = ~sd, x = ~mean, 
           color = ~sd,mode = 'markers',showlegend=T)%>%
  layout(
    xaxis = list(title = "mean")
    , yaxis = list(title = "sd"))


## Figure 3: Seurat result
Seu <- CreateSeuratObject(raw.data = as.data.frame(RNA),project = "Seu")
Seu <- ScaleData(object = Seu)
Seu <-RunPCA(object = Seu, pc.genes = rownames(Seu@data))
PCAPlot(Seu)
Seu <- FindClusters(object = Seu, reduction.type = "pca", dims.use = 1:9,  save.SNN = TRUE) #resolution = 0.6,
PrintFindClustersParams(Seu)
PrintCalcParams(Seu,calculation = 'RunPCA', raw = TRUE)
Seu<-RunTSNE(object = Seu, dims.use = 1:9, do.fast = TRUE, perplexity=10)
TSNEPlot(object = Seu)+ 
  xlim(-30,30)+ylim(-20,20)+
  labs(x = "tSNE1", y = "tSNE2")


## Figure 4: Select cluster K
BICvalue <- c(379599.7, 392696.1, 386628.4, 393845)
Kvalue <- c(1,2,3,4)
BICdata <- as.data.frame(cbind(Kvalue, BICvalue))
p<-plot_ly(data = BICdata, x = ~as.factor(Kvalue), y = ~BICvalue, mode = 'lines+markers',type = 'scatter')%>%
  layout(
    # title = "BIC values",
    xaxis = list(title = "K"),
    yaxis = list(title = "BIC"))

## Figure 5: Trace plots
z_line <- as.data.frame(cbind(c(1:2500),collection_z[,3,2]))
p1<-ggplot(z_line,aes(V1,V2))+
  geom_line(colour = "blue")+ ylim(2.6, 2.8) +labs(x = "iteration", y = expression(Z[ik]))

a_line <- as.data.frame(cbind(c(1:2500),collection_alpha[,5,2]))
p2<-ggplot(a_line,aes(V1,V2))+
  geom_line(colour = "blue")+ ylim(-4, 4)+labs(x = "iteration", y = expression(alpha[jk]))


b_line <- as.data.frame(cbind(c(1:2500),collection_beta[,5,2]))
p3<-ggplot(b_line,aes(V1,V2))+
  geom_line(colour = "blue")+ ylim(-5, 5)+labs(x = "iteration", y = expression(beta[gk]))


lam_line <- as.data.frame(cbind(c(1:2500),collection_lambda[,1,2]))
p4<-ggplot(lam_line,aes(V1,V2))+
  geom_line(colour = "blue")+ ylim(-3, 3)+labs(x = "iteration", y = expression(lambda[g1]))

grid.arrange(p1,p2,p3,p4)

## Inference
load("../data/RNA_DNAm_2.RData")

Z_est <- apply(collection_z, c(2,3), mean)

getMode <- function(x){
  return(as.numeric(names(table(x))[table(x) == max(table(x))][1]))
}

Gam_est <- apply(collection_Gam, 2, getMode)
Omega_est <- apply(collection_Omega, 2, getMode)

RNA <- read.csv("RNA_500.csv")
DEgene <- RNA[which(Gam_est == 1),1]

### Figure 6
d1 <- c(sum(Gam_est), 500-sum(Gam_est))
d2 <- c(1,0)
Gam <- cbind(d1,d2)
p1<-plot_ly(data = as.data.frame(Gam), y = ~d1, x= ~as.factor(d2),
            colorscale = "Portland", type = 'bar', name = "Gamma")%>%
  layout(xaxis = list(title = "Gamma indicator"),
         yaxis = list(title = ""))
d1 <- c(sum(Omega_est), 89-sum(Omega_est))
d2 <- c(1,0)
Omega <- cbind(d1,d2)
p2<-plot_ly(data = as.data.frame(Omega), y = ~d1, x= ~as.factor(d2),
            colorscale = "Portland", type = 'bar', name = "Omega")%>%
  layout(xaxis = list(title = "Omega indicator"),
         yaxis = list(title = ""))
p<-subplot(p1,p2)%>%
  layout(xaxis = list(title = "Indicator"),
         xaxis2 = list(title = "Indicator"))

### Figure 7
set.seed(20190330)
km_result <- kmeans(Z_est, 2)
cell <- as.data.frame(cbind(c(1:61), Z_est[,2],km_result$cluster))
p<-plot_ly(data = cell, x = ~cell[,1], y = ~cell[,2], symbol =as.factor(cell$V3),
           colorscale = "Portland", mode = 'markers',type = 'scatter',showlegend=T)%>%
  layout(xaxis = list(title = "Cell ID"),
         yaxis = list(title = "Z value"))
