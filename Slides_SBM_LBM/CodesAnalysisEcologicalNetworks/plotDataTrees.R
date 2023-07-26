rm(list=ls())
library(sbm)
library(GGally)
library(network)
library(RColorBrewer)
library(ggpubr)



mythemeTransp <- theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent'), #transparent legend panel
  axis.text=element_text(size=5)
)

#################################################""
data("fungusTreeNetwork")
##################################################  

Mat <- fungusTreeNetwork$fungus_tree
rownames(Mat) <- gsub("\\(.*","",fungusTreeNetwork$fungus_names)
colnames(Mat) <- gsub("\\(.*","",fungusTreeNetwork$tree_names)
U <- order(rownames(Mat))
V <- order(colnames(Mat))
Mat2 <- Mat[U,V]
bipartite::plotweb(Mat2,text.rot=90,method = 'normal',y.width.low = 0.03,y.width.high = 0.03, col.high = "#266D2F",bor.col.high="#266D2F", col.low = "#96691B",bor.col.low="#96691B",low.lablength = 20)

PlotMatrix <- sbm::plotMyMatrix(t(Mat2),dimLabels = list(row = "Tree", col = "Fungis"),plotOptions=list(title='Tree-Fungis',colNames = TRUE, rowNames = TRUE))+ mythemeTransp 
PlotMatrix
#ggsave(PlotMatrix, file="/home/donnet/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_08_IBCChannel_Course/Slides_SBM_LBM/slides_SBM_LBM_intro/plots/Tree_Fungis_matrix.png")


#------------- Degree Fungus
degree <- as.data.frame(c(rowSums(Mat2),colSums(Mat2)))
names(degree) <- 'degree'
degree$Nodes <- as.factor(rep(c('Fungus','Tree'), c(nrow(Mat2),ncol(Mat2)))) 
PlotDeg <- ggplot(degree,aes(x=degree,color=Nodes,fill=Nodes))  + geom_histogram(aes(y=..density..),alpha=0.5) +
  facet_grid( . ~ Nodes) + 
  scale_color_manual(values=c("#96691B", "#266D2F"))+ 
  scale_fill_manual(values=c("#96691B", "#266D2F"))+
  theme(legend.position="none")
ggsave(PlotDeg, file="/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_08_IBCChannel_Course/Slides_SBM_LBM/slides_SBM_LBM_intro/plots/Tree_Fungis_degree.png")


#----------Nestedness
Mat.nested <- nestedrank(Mat2, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=TRUE)$nested.matrix
PlotMatrix <- sbm::plotMyMatrix(t(Mat.nested),dimLabels = list(row = "Tree", col = "Fungis"),plotOptions=list(title='Tree-Fungis',colNames = TRUE, rowNames = TRUE))+ mythemeTransp 
PlotMatrix

Mat <- fungusTreeNetwork$tree_tree

rownames(Mat) <- colnames(Mat) <- gsub("\\(.*","",fungusTreeNetwork$tree_names)
U <- order(rownames(Mat))
Mat <- Mat[U,U]
net = network(Mat, directed = FALSE)
  
PlotNet = ggnet2(net,node.color = "#00BFC4", label = TRUE,label.size = 3,mode = "kamadakawai") +  mythemeTransp
PlotNet
ggsave(PlotNet, file="/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_06_MCMC/plots/Tree_network.png", height = 10 , units = "cm")
ggsave(PlotMatrix, file="/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_06_MCMC/plots/Tree_matrix.png", height = 20 , units = "cm")
  
  
 
