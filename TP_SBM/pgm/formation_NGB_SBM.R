setwd("/home/ouadah/Bureau/Extension_DROPBOX/vbemapp/gof-network/Data/FungiTree/")
setwd("~/Dropbox/Formation_NGB/")
tree = read.table("/home/donnet/Dropbox/Formation_NGB/data/FungiTree/fungi_tree_interactions.txt")
head(tree)
rowSums(tree[,-52])
dim(tree)
colSums(tree[-155,])
tree[155,]

tree_genet = read.table("/home/donnet/Dropbox/Formation_NGB/data/FungiTree/tree_genetic.txt")
tree_genet = as.matrix(tree_genet)
tree_taxo = read.table("/home/donnet/Dropbox/Formation_NGB/data/FungiTree/tree_taxo.txt")
tree_taxo = as.matrix(tree_taxo)
tree_geo = read.table("/home/donnet/Dropbox/Formation_NGB/data/FungiTree/tree_geo.txt")
tree_geo = as.matrix(tree_geo)

tree_fungi = as.matrix(tree[-155,-52])
tree = t(tree_fungi) %*% tree_fungi
diag(tree) <- 0

tree_bin = tree
tree_bin[tree > 1] = 1
n = dim(tree_bin)[1]


########### PLOT des données
source('function_for_blockmodels.R')
plotMatrix(tree_bin,rowFG = 'trees', colFG = 'trees')
plotMatrix(tree,rowFG = 'trees', colFG = 'trees')
plotMatrix(tree_fungi,rowFG = 'fungis', colFG = 'trees')




#------------------------------------------------
# SBM
#------------------------------------------------

sbm.tree_bin <- BM_bernoulli("SBM_sym",tree_bin)
m = sbm.tree_bin$estimate()

# nombre de groupes estime
Q = which.max(sbm.tree_bin$ICL)
Q

# blabla sur selection groupes via ICL

# proba a posteriori d'appartenance des noeuds aux groupes
membership=sbm.tree_bin$memberships[[Q]]$Z
node.membership=apply(membership, 1, which.max)

# graphique a faire et une fonction ppur
r = order(node.membership)
Z_ord = sort(node.membership)
sep =which(diff(Z_ord)!=0)
q=r
Adj_ord = tree_bin[r,q]

image(t(Adj_ord)[n:1,])
abline(v=sep/51, col="green")
abline(h=sep/51, col="green")

# probabilités de connexion entre groupes
pc=sbm.tree_bin$model_parameters[Q]

# probabilités a priori d'appartenance aux groupes
alpha=colSums(sbm.tree_bin$memberships[[Q]]$Z)/n

# graphe des groupes
library(igraph)
piS=sbm.tree_bin$model_parameters[[Q]]$pi
G = graph_from_adjacency_matrix(piS, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=alpha*100,edge.width=abs(E(G)$weight*2),vertex.color=1:Q)

#------------------------------------------------
# SBM value
#------------------------------------------------
sbm.tree <- BM_poisson("SBM",tree)
m=sbm.tree$estimate()

# nombre de groupes estime
Q=which.max(sbm.tree$ICL)
Q

# blabla sur selection groupes via ICL

# proba a posteriori d'appartenance des noeuds aux groupes
membership=sbm.tree$memberships[[Q]]$Z
node.membership=apply(membership, 1, which.max)

# graphique a faire
r = order(node.membership)
Z_ord = sort(node.membership)
sep =which(diff(Z_ord)!=0)
q=r
Adj_ord = tree_bin[r,q]

image(t(Adj_ord)[n:1,])
abline(v=sep/51, col="green")
abline(h=sep/51, col="green")

# probabilités de connexion entre groupes
pc=sbm.tree$model_parameters[Q]

# probabilités a priori d'appartenance aux groupes
alpha=colSums(sbm.tree$memberships[[Q]]$Z)/n

# graphe des groupes
piS=sbm.tree$model_parameters[[Q]]$lambda
G = graph_from_adjacency_matrix(piS, mode = c("undirected"), weighted = TRUE, diag = TRUE)
plot.igraph(G,vertex.size=alpha*100,edge.width=abs(E(G)$weight),vertex.color=1:Q)

# alluvial flow pour comparer les classifs

for(q in 1:Q)
{
  x=c(tree_genet[which(paramEstimSBM$Z==q),which(paramEstimSBM$Z==q)])
  x=x[-which(x==0)]
  mean(x)
}


#------------------------------------------------
# SBM value avce covariables
#------------------------------------------------

ListVar=list(tree_genet,tree_taxo,tree_geo)
sbm.cov<- BM_poisson_covariates("SBM",tree, ListVar)
# voir fast
m=sbm.cov$estimate()

# exemple qu'avec la genet et alluvial flow

#------------------------------------------------
# LBM
#------------------------------------------------
lbm.tree_fungi<- BM_bernoulli("LBM", tree_fungi)
m=lbm.tree_fungi$estimate()

# nombre de groupes estime : COMPRENDRE ce qui et fait
Q=which.max(lbm.tree_fungi$ICL)
Q
# comment recuperer Q1 et Q2

# proba a posteriori d'appartenance des noeuds aux groupes
# membership
lbm.tree_fungi$memberships[[Q]]$Z1
lbm.tree_fungi$memberships[[Q]]$Z2

# probabilités a priori d'appartenance aux groupes
alphaZ1=colSums(lbm.tree_fungi$memberships[[Q]]$Z1)/n

# probabilités de connexion entre groupes
pc=lbm.tree_fungi$model_parameters[Q][[1]]$pi

# graphe des groupes : trouver comment faire
piS=pc
G = graph_from_incidence_matrix(piS, weighted = TRUE)
plot.igraph(G,edge.width=abs(E(G)$weight),vertex.color=1:Q,  layout=layout.bipartite)

plotNetBM = function(BMobject,Q){
  membership_name <-  BMobject$membership_name
  a <- extractParamBM(BMobject,Q)$alpha

  if (membership_name == 'SBM'){
    G <- graph_from_adjacency_matrix(a, mode = c("undirected"), weighted = TRUE, diag = TRUE)
    g <- plot.igraph(G,vertex.size=alpha*100,edge.width=abs(E(G)$weight)*2,vertex.color=1:Q)
  }

  if (membership_name == 'LBM'){
    G <- graph_from_incidence_matrix(a, weighted = TRUE)
    g <- plot.igraph(G,edge.width=abs(E(G)$weight),vertex.color=1:Q, layout=layout.bipartite)
  }

}
