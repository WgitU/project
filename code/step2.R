setwd("../data")
## Preprocess scRNA-seq data
RNA <- read.csv("GSE74534_RNA-seq_raw_counts.csv")
scMT <- readRDS('GSE74534_DNAm.rds')
RNA_sub <- filter(RNA, RNA[,1]%in%id_vec)
scMT_sub <- select(scMT, c(sample, ens_id, weight))
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
write.table(RNA_500,"RNA_500.csv",row.names=F,col.names=TRUE,sep=",")

## Preprocess DNA methylation data
files<-list.files(pattern = "*.cov")
datalist <- list()
num <- 1
col_name <- NULL
for(i in files) { 
  x <- as.data.frame(read_delim((i),delim="\t",col_names = F,col_types =list(col_integer(),col_integer(),col_integer(),col_double(),col_integer(),col_integer())))
  x <- x[,c(1,2,5)]
  x <- x %>%
    group_by(X1,X2) %>%
    summarise(X5 = sum(X5))
  datalist[[num]] <- assign(i, x)
  num <- num + 1
  na_me<-sub('....$','',i)
  col_name <- c(col_name,na_me)
}
datamat <- Reduce(function(...) merge(..., by=c("X1","X2"), all=F), datalist)
names(datamat)[3:63]<- col_name
DNAm <- datamat
DNAm <-DNAm %>%
  arrange(X1,X2)
DNAm<-na.omit(DNAm)
Prom_cgi <- filter(scMT, scMT$name == "prom_cgi")
Prom_cgi_unique <- Prom_cgi[!duplicated(Prom_cgi$start.x),] %>%
  arrange(chromo.x,start.x)
ind <- NULL
for (i in 1:917){
  sub_i <- Prom_cgi_unique[Prom_cgi_unique$chromo.x == DNAm$X1[i],10:11]
  ind_left <- max(which(DNAm$X2[i] >= sub_i[,1]),na.rm =T)
  ind_right <- min(which(DNAm$X2[i] <= sub_i[,2]),na.rm =T)
  if(ind_left == ind_right) {
    ind <- c(ind, T)
  } else {
    ind <- c(ind, F)
  }
}
DNAm_prom <- DNAm[ind,]
write.table(DNAm_prom, file = "DNAm_pre.txt")