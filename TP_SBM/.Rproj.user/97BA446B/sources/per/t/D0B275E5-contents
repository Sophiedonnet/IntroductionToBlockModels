rm(list=ls())
library(sbm)
library(network)
library(RColorBrewer)

source('functionsPlots.R')


if(Sys.info()[7]=='donnet'){
  whereSave <- "/home/donnet/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_08_IBCChannel_Course/Slides_SBM_LBM/slides_SBM_LBM_intro/plots/"
}else{
  whereSave <- "/home/sophie/Dropbox/WORK_DROPBOX/RECHERCHE_en_cours/CONFERENCE_EXPOSES/EXPOSES/2023/2023_08_IBCChannel_Course/Slides_SBM_LBM/slides_SBM_LBM_intro/plots/"
}



###################################################"
###########" Trial Communities
#####################################################
K = 3; 
nbNodes  <- 50
blockProp <- c(.25, 0.5 ,.25) # group proportions
means <- diag(.4, 3) + 0.05  # connectivity matrix: affiliation network
connectParam <- list(mean = means)



mySampler <- sampleSimpleSBM(nbNodes, blockProp, directed = TRUE, connectParam,dimLabels = 'Species', model = 'bernoulli')

mythemeTransp <- theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)

#--------------------------- raw data
colSet2 <- brewer.pal(n = 8, name = "Set2")
rawNet <- ggnet2(mySampler$networkData,label = TRUE,color=colSet2[1]) + mythemeTransp
ggsave(rawNet, file = paste0(whereSave,"network_raw.png"), height = 10 , units = "cm",width = 10) 
PlotMatrixRaw  <- plotMyMatrix(mySampler$networkData,dimLabels = list(row = "species", col = "species"),plotOptions=list(colNames = TRUE, rowNames = TRUE))+ mythemeTransp 
ggsave(PlotMatrixRaw, file = paste0(whereSave,"matrix_raw.png"), height = 20 , units = "cm")

myPlot <- myPlotForPaper(mySampler,myTitle='')
ggsave(myPlot[[2]], file = paste0(whereSave,"network_SimuCommunities.png"), height = 10 , units = "cm",width = 10)
ggsave(myPlot[[1]], file = paste0(whereSave,"matrix_SimuCommunities.png"), height = 20 , units = "cm")


########################################
###########" Trial  Foodwebs
######################################### 
K = 4; 
nbNodes  <- 50
blockProp <- c(0.2, .25, 0.30 ,.25) # group proportions
alpha <- matrix(0.02,K,K)
alpha[2,1] = 0.5
alpha[3,1] = 0.5
alpha[3:4,2] = 0.4
alpha[3:4,3] = 0.4
diag(alpha) = 0.1


connectParam <- list(mean = alpha)

mySampler <- sampleSimpleSBM(nbNodes, blockProp, directed = TRUE, connectParam,dimLabels = 'Species', model = 'bernoulli')
myPlot <- myPlotForPaper(mySampler,myTitle = NULL)
myPlot
ggsave(myPlot[[2]], file = paste0(whereSave,"network_SimuFoodWeb.png"),  height = 10 , units = "cm")
ggsave(myPlot[[1]], file = paste0(whereSave,"matrix_SimuFoodWeb.png"), height = 20 , units = "cm")





