setwd("D:/submission/others/energy/z25/Figure7/relation")


library(psych)
library(reshape2)
library(pheatmap)

df1<- read.csv("factors.csv", header = T, row.names = 1)
df2 <-  read.csv("met.csv", header = T, row.names = 1)


cor <- corr.test(df1, df2, method = "pearson",adjust="none")              #计算相关性矩阵、p值矩阵
cmt <- cor$r
pmt <- cor$p
cmt.out<-cbind(rownames(cmt),cmt)
pmt.out<-cbind(rownames(pmt),pmt)

df <- melt(cmt,value.name="cor")
df$pvalue <- as.vector(pmt)
head(df)


if (!is.null(pmt)) {
  ssmt <- pmt < 0.01
  pmt[ssmt] <- '**'
  smt <- pmt > 0.01& pmt <0.05
  pmt[smt] <- '*'
  pmt[!ssmt&!smt] <- ''
} else {
  pmt <- F
}

colormap1 <- colorRampPalette(colors =c("purple","white", "orange"))(30) #粉-金

pmt
pheatmap(cmt, scale = "none", cluster_row = F, cluster_col = F,                
         display_numbers = pmt, fontsize_number = 12, number_color = "white",
         cellwidth = 20, cellheight = 20,filename="heatmap.pdf", color = colormap1)  

#输出：相关性热图pdf文件
pheatmap(cmt, scale = "none",cluster_row = F, cluster_col = F,
         display_numbers = pmt, fontsize_number = 12, number_color = "white",
         cellwidth = 20, cellheight = 20)

