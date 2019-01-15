##Step 1: download data from GEO
library(GEOquery)
setwd("C:\\Users\\AM\\Desktop\\Diff Analysis")
gset<-getGEO("GSE57820",destdir = ".\\",GSEMatrix = TRUE, AnnotGPL = TRUE)
exprSet<-read.table(gzfile("GSE57820_series_matrix.txt.gz"),comment.char = "!",
                 header = T,row.names=1)#sep=""表示原文件是以空格进行分隔,comment.char = "!"

##Step 2: make proper column names to match toptable
library(illuminaHumanv4.db)
ids=toTable(illuminaHumanv4SYMBOL)#匹配geneid对应的genesymbol
plot(table(sort(table(ids$symbol))))#绘制出匹配状况

table(rownames(exprSet) %in% ids$probe_id)#以表格的形式判断在ids中的数量
dim(exprSet)#显示有多少行和列
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]#通过判断数量判断是否取了
dim(exprSet)

ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]
tmp = by(exprSet,ids$symbol,function(x) rownames(x)[which.max(rowMeans(x))] )#出现多个探针对应一个gensymbol时，取最大平均表达值对应的探针
probes = as.character(tmp)
exprSet=exprSet[rownames(exprSet) %in% probes ,]
dim(exprSet)

rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]#将probe_id转换为genesymbol
exprSet[1:5,1:5]

##Step3: hclust and PCA of samples
exprSet["GAPDH",]
exprSet["ACTB",]
boxplot(exprSet,1)#声明：如果含有批次效应，需使用sv包中的combine函数
group_list=c(rep("transfected",6),rep("untransfected",6))#设定分组信息
### hclust
colnames(exprSet)=paste(group_list,1:12,sep='')# Define nodePar
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))
pdf(file="hclust_samples.pdf")
par(mai=c(1,0.5,1,1))
plot(as.dendrogram(hc), main="hclust of samples",nodePar = nodePar, horiz = TRUE)
dev.off()
### PCA
library(ggfortify)
df=as.data.frame(t(exprSet))
df$group=group_list
pdf(file="PCA_samples.pdf")
autoplot(prcomp( df[,1:(ncol(df)-1)] ), main="PCA of samples",data=df,colour = 'group')
dev.off()

##Step 4: DEGs Analysis
library(limma)
design <- model.matrix(~0+factor(group_list))#序列矩阵
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)#如果相反，应该怎么处理？？
contrast.matrix

fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2) ## default no trend !!!
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG = na.omit(tempOutput)
write.table(DEG,"limma_notrend.results.txt",sep='\t',quote=F,row.names=T,col.names=T)

##Step 5:Volcano plot
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs(logFC)) )#95%以上的观测值
logFC_cutoff<-1###一般要求logFC>1
DEG$change <- as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
diff<-DEG[which((DEG$change == "UP")|(DEG$change == "DOWN")),]
write.table(diff,"diff_genesymbol.foldchange.txt",sep='\t',quote=F,row.names=T)

diff$genesymbol<-rownames(diff)
loc<-match(diff$genesymbol,rownames(exprSet))
DEG_exp<-exprSet[loc,]
genesymbol<-rownames(DEG_exp)
DEG_exp<-cbind(genesymbol,DEG_exp)
write.table(DEG_exp,file="diff_genesymbol.FDR0.05.exprs.txt",sep='\t',quote=F,row.names=F)

pdf(file="volcono.pdf")
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
dev.off()

##Step 6:Heatmap
library(pheatmap)
DEG_exp<-read.table("diff_genesymbol.FDR0.05.exprs.txt",sep='\t',header=T,row.names=1)
group_info <- data.frame(row.names=names(DEG_exp),groups=group_list)
pdf(file="pheatmap.pdf")
pheatmap(DEG_exp,fontsize_row=1,color=colorRampPalette(c("blue", "white", "red"))(256)
         ,fontsize_col=1,scale="row",border_color=NA,cluster_col = FALSE,annotation_col=group_info)
dev.off()
