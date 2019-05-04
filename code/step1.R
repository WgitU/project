setwd("../data")

## scRNA-seq data
RNA <- read.csv("GSE74534_RNA-seq_raw_counts.csv")
scMT <- readRDS('GSE74534_DNAm.rds')
id_vec <- unique(scMT$ens_id)
RNA_sub <- filter(RNA, RNA[,1]%in%id_vec)
scMT_sub <- select(scMT, c(sample, ens_id, weight))
rownames(RNA_sub) <- RNA_sub[,1]
RNA_sub <- RNA_sub[,-1]
RNA <- RNA_sub
head(RNA)
tail(RNA)
dim(RNA)
head(RNA)
summary(RNA)

### zero-inflated effects
colSums(RNA == 0)/dim(RNA)[1]

### Plots for scRNA-seq data
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

#### visualization of 500 most variable genes 
vals <- unique(scales::rescale(c(volcano)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = as.matrix(RNA_500),
        zmax=10000,zmin=0,
        colorscale = colz, type = "heatmap")%>%
  layout(
    title = "Heatmap of scRNA-seq data",
    xaxis = list(title = "Cell")
    , yaxis = list(title = "Gene"))

#### Mean and standard deviations of single cells 
mean_sd <- as.data.frame(cbind(colSums(RNA_norm), apply(RNA_norm,2,sd)))
colnames(mean_sd) <- c("mean","sd")
plot_ly(data = mean_sd ,y = ~sd, x = ~mean, 
        color = ~sd,mode = 'markers',showlegend=T)%>%
  layout(
    title = "Mean and standard deviations of single cells",
    xaxis = list(title = "mean")
    , yaxis = list(title = "sd"))

#### visualization of zero-inflated effects
zero_ratio <- colSums(RNA_norm == 0)/dim(RNA_norm)[1]
rate_count <- as.data.frame(cbind(zero_ratio, colSums(RNA_norm)))
plot_ly(data = rate_count ,y = ~rate_count[,1], x = ~rate_count[,2], 
        color = ~zero_ratio,mode = 'markers',showlegend=T)%>%
  layout(title = "Zero-inflated effects",
         xaxis = list(title = "total read counts of each single cell"),
         yaxis = list(title = "zero ratio"))

zero_ratio_r <- rowSums(RNA_norm == 0)/dim(RNA_norm)[2]
zero_ratio_r <- as.data.frame(cbind(zero_ratio_r, rowSums(RNA_norm)))
plot_ly(data = zero_ratio_r ,y = ~zero_ratio_r, x = ~V2, 
        color = ~zero_ratio_r,mode = 'markers',showlegend=T)%>%
  layout(title = "Zero-inflated effects",
         xaxis = list(title = "total read counts of each single gene"),
         yaxis = list(title = "zero ratio"))

## a example of DNA methylation data
A02 <- as.data.frame(read_delim(A02.cov,delim="\t",col_names = F,col_types =list(col_integer(),col_integer(),col_integer(),col_double(),col_integer(),col_integer())))
A02 <- A02[,c(1,2,5)]
head(A02)
tail(A02)
dim(A02)

## K-means
df<-t(RNA_norm)
fviz_nbclust(df, kmeans, method = "wss",k.max=10) + geom_vline(xintercept = 6, linetype = 2)
set.seed(159753)
km_result <- kmeans(df, 6)
fviz_cluster(km_result, data = df,
             ellipse = F,
             show.clust.cent = F,
             shape = 1,
             repel = F,
             ggtheme = theme_grey()
)

