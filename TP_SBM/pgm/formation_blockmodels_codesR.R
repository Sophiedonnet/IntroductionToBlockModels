## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)


## ----import dataset,  echo = TRUE--------------------------------------------------------------------------------------------------------------
rm(list=ls())
load('fungi_tree_data.Rdata')
ls()


## ----load pack,  echo = TRUE,message=FALSE-----------------------------------------------------------------------------------------------------
library(ggplot2)
library(igraph)
library(blockmodels)
library(alluvial)



## ----source function,  echo = TRUE, message=FALSE----------------------------------------------------------------------------------------------
source('function_for_blockmodels.R')


## ----tree tree bin plot, echo=TRUE-------------------------------------------------------------------------------------------------------------
plotMatrix(Mat = tree_bin,rowFG = 'tree', colFG  = 'tree')


## ----SBM, echo=TRUE, eval = TRUE---------------------------------------------------------------------------------------------------------------
sbm.tree_bin <- BM_bernoulli("SBM_sym",tree_bin)


## ----estimate SBM, echo=TRUE, eval = FALSE-----------------------------------------------------------------------------------------------------

sbm.tree_bin$estimate()


## ----select SBM, echo=TRUE, eval = TRUE--------------------------------------------------------------------------------------------------------
Q = which.max(sbm.tree_bin$ICL)
Q


## ----extract param SBM, echo=TRUE, eval = TRUE-------------------------------------------------------------------------------------------------
paramEstimSBM <- extractParamBM(sbm.tree_bin,Q)
paramEstimSBM$pi

paramEstimSBM$alpha
paramEstimSBM$Z


## ----plot org treeb bin,  echo=TRUE, eval = TRUE-----------------------------------------------------------------------------------------------
plotMatrix(tree_bin,'tree','tree', fileNameSave = NULL, clustering = list(row=paramEstimSBM$Z))


## ----plot BM network treebin,  echo=TRUE -----------------------------------------------------------------------------------------
par(mfrow = c(1,1))
G <- graph_from_adjacency_matrix(paramEstimSBM$alpha, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=paramEstimSBM$pi*100,edge.width=abs(E(G)$weight)*3,vertex.color=1:Q, layout=layout_nicely)


## ----list names blocks,  echo=TRUE, eval = TRUE------------------------------------------------------------------------------------------------
lapply(1:Q,function(q){tree_list[paramEstimSBM$Z == q]})


################################################## SBM poisson
## ----tree tree plot, echo=TRUE-----------------------------------------------------------------------------------------------------------------
plotMatrix(Mat = tree,rowFG = 'tree', colFG  = 'tree')


## ----SBM Poisson, echo=TRUE, eval = TRUE-------------------------------------------------------------------------------------------------------
sbm.tree <- BM_poisson("SBM_sym",tree)


## ----estimate SBM poisson, echo=TRUE, eval = FALSE---------------------------------------------------------------------------------------------
sbm.tree$estimate()


## ----select SBM poisson, echo=TRUE, eval = TRUE------------------------------------------------------------------------------------------------
Q = which.max(sbm.tree$ICL)
Q


## ----extract param SBM poisson, echo=TRUE, eval = TRUE-----------------------------------------------------------------------------------------
paramEstimSBMPoisson <- extractParamBM(sbm.tree,Q)
paramEstimSBMPoisson$pi
paramEstimSBMPoisson$lambda
paramEstimSBMPoisson$Z


## ----plot org tree ,  echo=TRUE, eval = TRUE---------------------------------------------------------------------------------------------------

plotMatrix(tree,'tree','tree', fileNameSave = NULL, clustering = list(row=paramEstimSBMPoisson$Z))


## ----plot BM network tree poisson, echo=TRUE, eval = FALSE-------------------------------------------------------------------------------------

par(mfrow = c(1,1))
G <- graph_from_adjacency_matrix(paramEstimSBMPoisson$lambda, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=paramEstimSBMPoisson$pi*100,edge.width=abs(E(G)$weight)*3,vertex.color=1:Q, layout=layout_nicely)


## ----list names blocks Poisson,  echo=TRUE, eval = TRUE----------------------------------------------------------------------------------------
lapply(1:Q,function(q){tree_list[paramEstimSBMPoisson$Z == q]})


## ----alluvial, echo=TRUE,eval=TRUE-------------------------------------------------------------------------------------------------------------
A <- as.data.frame(table(paramEstimSBM$Z,paramEstimSBMPoisson$Z))
colnames(A)=c('SBM Bern',"SBM Poisson","Freq")
w   <- which(A$Freq != 0)
A <- A[w,]
alluvial(A[,c(1,2)],freq=A$Freq)

################################################## SBM poisson avec covariables
## ----covar SBM,echo=TRUE,eval=FALSE------------------------------------------------------------------------------------------------------------


sbm.cov <- BM_poisson_covariates("SBM_sym",adj = tree, covariates = ListVar)
sbm.cov$estimate()


## ----load SBM covar, echo=FALSE, eval = TRUE---------------------------------------------------------------------------------------------------


## ----select SBM covar, echo=TRUE, eval = TRUE--------------------------------------------------------------------------------------------------
Q = which.max(sbm.cov$ICL)
Q


## ----extract param SBM poisson covar, echo=TRUE, eval = TRUE-----------------------------------------------------------------------------------
paramEstimSBMPoissonCov <- extractParamBM(sbm.cov,Q)
paramEstimSBMPoissonCov$pi
paramEstimSBMPoissonCov$lambda
paramEstimSBMPoissonCov$Z
paramEstimSBMPoissonCov$theta


## ----plot org tree  cov,  echo=TRUE, eval = TRUE-----------------------------------------------------------------------------------------------
plotMatrix(tree,'tree','tree', fileNameSave = NULL, clustering = list(row=paramEstimSBMPoissonCov$Z))


## ----plot BM network tree cov,  echo=TRUE, eval = FALSE----------------------------------------------------------------------------------------
G <- graph_from_adjacency_matrix(paramEstimSBMPoissonCov$lambda, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size = paramEstimSBMPoissonCov$pi*100,edge.width=sqrt(abs(E(G)$weight)),vertex.color = 1:Q, layout = layout_nicely)


## ----list names blocks poisson cov,  echo=TRUE, eval = TRUE------------------------------------------------------------------------------------
lapply(1:Q,function(q){tree_list[paramEstimSBMPoissonCov$Z == q]})


## ----alluvial cov, echo=TRUE,eval=TRUE---------------------------------------------------------------------------------------------------------
B <- as.data.frame(table(paramEstimSBM$Z,paramEstimSBMPoisson$Z,paramEstimSBMPoissonCov$Z))
colnames(B) = c('SBM Bern',"SBM Poisson","SBM Poisson + Cov","Freq")
w   <- which(B$Freq!=0)
B <- B[w,]
alluvial(B[,c(1,2,3)],freq=B$Freq)


################################################## LBM


## ----tree fungis plot, echo=TRUE---------------------------------------------------------------------------------------------------------------
plotMatrix(Mat = fungi_tree,rowFG = 'fungi', colFG  = 'tree')


## ----LBM,echo=TRUE,eval=FALSE------------------------------------------------------------------------------------------------------------------
lbm <- BM_bernoulli("LBM",as.matrix(fungi_tree))
lbm$estimate()



## ----select lBM covar, echo=TRUE, eval = TRUE--------------------------------------------------------------------------------------------------
Q = which.max(lbm$ICL)
Q
paramEstimLBM <- extractParamBM(lbm,Q)
paramEstimLBM$Q



## ----extract param LBM, echo=TRUE, eval = TRUE-------------------------------------------------------------------------------------------------
paramEstimLBM$piRow
paramEstimLBM$piCol
paramEstimLBM$alpha
paramEstimLBM$ZRow
paramEstimLBM$ZCol


## ----plot org  tree fungis,  echo=TRUE, eval = TRUE--------------------------------------------------------------------------------------------
plotMatrix(fungi_tree,'fungi','tree', fileNameSave = NULL, clustering = list(row = paramEstimLBM$ZRow,col = paramEstimLBM$ZCol))

## ----plot BM network tree fungi,  echo=TRUE, eval = TRUE---------------------------------------------------------------------------------------
G <- graph_from_incidence_matrix(paramEstimLBM$alpha, weighted = TRUE)
plot(G,vertex.size=c(paramEstimLBM$piRow*100, paramEstimLBM$piCol*100), vertex.shape=c("circle", "square")[V(G)$type +1], edge.width=abs(E(G)$weight*2),vertex.color=1:Q, layout=layout.bipartite)

