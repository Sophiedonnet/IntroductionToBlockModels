setwd("~/Dropbox/Multiplex/Ecologie/Exposes/2017-06-Inecol-Xalapa/Code_R_plots")

library(ggplot2)
library(igraph)
source('func_sampling.R')


List_type_graph <- c("Affiliation", "Affiliation2","Bipartite","Etoile","Mixte","Nested")

epsilon <- 0.09
rho= .6
n <- 200
M = matrix(epsilon,4,4)
M[4,4]=M[3,3] = 0.7
M[1,1] = M[1,2] = M[2,1] = 0.5

Mmodnest=matrix(epsilon,nrow=6,ncol=6)
Mmodnest[1,1:3]=rho
Mmodnest[2,1:2]=rho
Mmodnest[3,1]=rho
Mmodnest[4,4:6]=rho
Mmodnest[5,4:5]=rho
Mmodnest[6,4]=rho
type_graph="autre"

for (type_graph in List_type_graph){
  print(type_graph)
  alpha <- switch(type_graph,
                  "Affiliation" = c(1/5,2/5,2/5),
                  "Affiliation2" = c(.2,.6,.2),
                  "Nested" = rep(1/4,4),
                  "Bipartite"   = c(1/4,1/4,1/4,1/4),
                  "Etoile"      = c(.15,.55,.15,.15),
                  "Mixte"      = c(.2,.5,.15,.15),
                  "autre"      = rep(1/6,6))
  
  pi <- switch(type_graph,
               "Affiliation" = matrix(c(rho,epsilon,epsilon,epsilon,rho,epsilon,epsilon,epsilon,rho),3,3),
               "Nested" = matrix(c(rho,rho,rho,rho,rho,rho,rho,epsilon,rho,rho,epsilon,epsilon,rho,epsilon,epsilon,epsilon),4,4,byrow=TRUE),
               "Affiliation2" = matrix(c(epsilon,epsilon/10,epsilon/10,epsilon/10,epsilon,epsilon/10,epsilon/10,epsilon/10,epsilon),3,3),
               "Bipartite"   = matrix(c(epsilon,rho,epsilon,epsilon,rho,epsilon,epsilon,epsilon,epsilon,epsilon, epsilon,rho,epsilon,epsilon,rho,epsilon),4,4),
               "Etoile"      = matrix(c(rho,rho,epsilon,epsilon,rho,epsilon,epsilon,epsilon,epsilon,epsilon,rho,rho,epsilon,epsilon,rho,epsilon),4,4),
               "Mixte"      = M,
               "autre" =Mmodnest)
  
  
  Q=length(alpha)
  
  sbm <- rSBM(n,Q,alpha,pi)
  
  
  R = graph_from_adjacency_matrix(sbm$X,mode="undirected")
  layout <- layout.fruchterman.reingold(R)

  plot(R,layout=layout,vertex.size = 7,vertex.label=NA)

  plot(R,layout=layout,vertex.color=sbm$cl,vertex.size=7,vertex.label=NA)

  piS=pi
  piS[pi<.1]=0
  G = graph_from_adjacency_matrix(piS, mode = c("undirected"), weighted = TRUE, diag = TRUE)
  plot.igraph(G,vertex.size=alpha*100,edge.width=E(G)$weight*3,vertex.color=1:Q)  
  
  Adj=sbm$X
  index_plants = rep(1:dim(Adj)[1],each=dim(Adj)[2])
  index_animals =- rep(1:dim(Adj)[2],dim(Adj)[1])
  
  melted_Adj= data.frame(plants=index_plants , animals= index_animals)
  link = rep(-10,dim(Adj)[2]*dim(Adj)[1])
  for (k in 1:(dim(Adj)[2]*dim(Adj)[1])){link[k] = Adj[index_plants[k],-index_animals[k]]}
  melted_Adj$link = link
  melted_Adj$index_plants = index_plants
  melted_Adj$index_animals = index_animals
  
  
  g<- ggplot(data = melted_Adj, aes(x=index_plants, y=index_animals, fill=link)) + geom_tile()
  g
  #reordonnée
  r = order(sbm$cl)
  Z_ord = sort(sbm$cl)
  sep =which(diff(Z_ord)!=0)
  q = r ;
  Adj_ord = Adj[r,q]
  
  melted_Adj_ord = data.frame(plants=index_plants,animals =index_animals)
  link_ord = rep(-10,dim(Adj)[2]*dim(Adj)[1])
  for (k in 1:(dim(Adj)[2]*dim(Adj)[1])){link_ord[k] = Adj_ord[index_plants[k],-index_animals[k]]}
  melted_Adj_ord$link = link_ord
  melted_Adj_ord$index_plants = index_plants
  melted_Adj_ord$index_animals = index_animals
  
  
  
  g<- ggplot(data = melted_Adj_ord, aes(x=index_plants, y=index_animals, fill=link)) + geom_tile() +ggtitle("Reordered adjacency matrix")
  g <- g + theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())
  g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
  g <- g + geom_hline(yintercept=-sep,linetype="dashed",color='orange')+ geom_vline(xintercept=sep,linetype="dashed",color='orange')
  g
  
  
  
  #ex blockmodels
  library(blockmodels)
  SBMex=BM_bernoulli("SBM_sym",Adj)
  SBMex$estimate()
  SBMex  
  k=which.max(SBMex$ICL)
  Zest=apply(SBMex$memberships[[k]]$Z,1,which.max)
  table(Zest,sbm$cl)
library(mclust)
  adjustedRandIndex(Zest,sbm$cl)    

  
  #modularity
  groupmod=cluster_walktrap(R)
  table(groupmod$membership,sbm$cl)  
  modularity(groupmod)
  

  #nestedness
  library(vegan)
  out=nestedtemp(Adj)
  out
  plot(out,kind="incid")
  out
  oecosimu(Adj, nestedchecker, "quasiswap")
  oecosimu(Adj, nestedchecker, "r00", statistic = "C.score")
  
  
  