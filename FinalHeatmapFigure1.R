###Reyka Jayasinghe
###November 19th, 2015
###rjayasin@genome.wustl.edu

###Script to make heatmap for clustering splicing defects across mutations

#Packages to install
#install.packages("gplots")
#install.packages("devtools")

#Load necessary packages
#library("gplots")
#library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#Set a working directory for output files
setwd("/Users/rjayasin/Desktop/")
#loaddata
TEST_large <- read.delim("~/Desktop/OneDrive/Ding_Lab/NMD_Alternative_Splicing/comparison_october2015/Heatmap/TEST_large")
TEST_large <- read.delim("~/Desktop/OneDrive/Ding_Lab/NMD_Alternative_Splicing/comparison_october2015/Heatmap/TEST_large_exp")
#Extract intron retention and cancer information from Heatmap
IRStatus<-as.matrix(TEST_large$Intron.Retention)
CancerTypeStatus<-as.matrix(TEST_large$Cancer.Type)
ExpressionStatus<-as.matrix(TEST_large$Gene.Expression)
Status<-t(cbind(CancerTypeStatus,ExpressionStatus,IRStatus))
#Set Mutation information as rownames and delete associated columns in matrix
rownames(TEST_large)=TEST_large$Sample.Name.and.Gene
TEST_large$Sample.Name.and.Gene<- NULL
TEST_large$Intron.Retention<- NULL
TEST_large$Cancer.Type<- NULL
#TEST_large$Cancer.Type.1<- NULL
TEST_large$Gene.Expression<- NULL
TEST<-as.matrix(TEST_large)

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Set all NA values to 0 for mapping purposes
TEST[is.na(TEST)] <- 0

#Rename column and transpose matrix for use in RowSideColors in heatmap.3 below
colnames(IRStatus)=c("Intron Retention")
rownames(Status)=c("Cancer Type","Gene Expression","Intron Retention")
colnames(TEST)=c("Multi Exon Skipping","Exon Skipping","3' Exon Shrinkage", "5' Exon Shrinkage", "3' Exon Extension", "5' Exon Extension")
rlabme<-Status

#####HEATMAP
#Heatmap using two color gradients for reads 0-upper value and uppervalue to uppervalue2 
uppervalue=30
uppervalue2=60
breaks = seq(0,uppervalue2+1,length.out=1000)
#Setting up gradients for different break cutoffs
gradient1 = colorpanel( sum( breaks[-1]<=1 ), "white", "white") #Everything that has 0 mark it white
gradient2 = colorpanel( sum( breaks[(breaks[-1]>1)]<=uppervalue ),"khaki1","coral1" )#color gradient for lower read counts
gradient3 = colorpanel( sum( breaks[(breaks[-1]>uppervalue)]<=uppervalue2 ),"coral2","brown4" )#color gradient for higher read counts
gradient4 = colorpanel( sum( breaks[-1]>uppervalue2), "black","black" )#everything greater than uppervalue2 is set to black
hm.colors = c(gradient1,gradient2,gradient3,gradient4)


###Final printout of PDF
pdf(file="FinalHeatmapSplicingDefects.pdf")
main_title="Splicing Defects of Highly Conserved Splice Site Mutations"
par(cex.main=0.6)
hh<-heatmap.3(as.matrix(TEST),hclustfun=myclust, distfun=mydist, RowSideColors=rlabme,breaks=breaks,col=hm.colors[1:(length(hm.colors)-2)],margins=c(10,15),cexCol=0.7, labRow=FALSE,main=main_title)
legend("topright",legend=c("High Intron Retention","Low Intron Retention","No Intron Retention","BLCA","BRCA","COAD","GBM","HNSC","KIRC","LAML","LUAD","LUSC","OV","READ","UCEC","STAD","LGG","PRAD","CESC","Expressed","Not Expressed","No expression data"),
       fill=c("blue","cyan","white","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#08B1E6","#FCE452","#B2D49C","#B51C98","green","red","white"), border=TRUE, bty="n", y.intersp = 1, cex=0.6)

dev.off()


#Extract clustering information from heatmap and save to dataframe
sorted <- TEST[match(rev(labels(hh$rowDendrogram)), rownames(TEST)),]
write.csv(sorted,"/Users/rjayasin/Desktop/SortedFile_Dendrogram",quote=FALSE)

