fungi_tree = read.table("FungiTree/fungi_tree_interactions.txt")
fungi_tree = fungi_tree[-dim(fungi_tree)[1],
                        -dim(fungi_tree)[2]]

tree = t(as.matrix(fungi_tree)) %*% as.matrix(fungi_tree)
diag(tree) <- 0

tree_bin = tree
tree_bin[tree > 1] = 1

tree_genet = read.table("FungiTree/tree_genetic.txt")
tree_genet = as.matrix(tree_genet)
tree_taxo = read.table("FungiTree/tree_taxo.txt")
tree_taxo = as.matrix(tree_taxo)
tree_geo = read.table("FungiTree/tree_geo.txt")
tree_geo = as.matrix(tree_geo)
ListVar = list(tree_genet,tree_taxo,tree_geo)
ListVarscale <- ListVar
ListVarscale[[1]] <- ListVar[[1]]/max(c(ListVar[[1]]))
ListVarscale[[2]] <- ListVar[[2]]/max(c(ListVar[[2]]))
ListVar <- ListVarscale


tree_list= read.table("nouvelles_data/Tree_species_list.csv", sep=";")
tree_list=tree_list[,2]

fungi_list= read.table("nouvelles_data/Fungal_species_list.csv", sep=";")
fungi_list=fungi_list[,2]

fungi_tree_data= list(Name='data_list', fungi_tree=fungi_tree, tree=tree, tree_bin=tree_bin, ListVar=ListVar, tree_list=tree_list, fungi_list=fungi_list)
# liste avec fungi_tree, tree, tree_bin, ListVar, tree_list, fungi_list
save(fungi_tree, tree, tree_bin, ListVar, tree_list, fungi_list,file='fungi_tree_data.Rdata')
