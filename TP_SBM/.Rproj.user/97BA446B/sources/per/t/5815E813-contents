###############################################"

myPlotForPaper = function(mySampler,myTitle = NULL){
  
  if(is.null(myTitle)){myTitle = ''}
  myMat <-mySampler$networkData
  nbNodes <- nrow(myMat)
  row.names(myMat)<- 1:nbNodes
  colnames(myMat)<- 1:nbNodes
  
  ##################################
  U <- order(mySampler$memberships)
  myMat <- myMat[U,U]
  Z <- mySampler$memberships[U]
  #############################
  net = network(myMat, directed = TRUE)
  
  mythemeTransp <- theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'), #transparent legend panel,
    legend.position = "none"
    
  )
  
  K <- max(Z)
  colSet2 <- brewer.pal(n = 8, name = "Set2")
  my_color_nodes_col <- my_color_nodes <- rep(colSet2[1],nbNodes)
  for (k in 2:K){
    xk <- which(Z==k)
    my_color_nodes[nbNodes-xk + 1] <- colSet2[k]
    my_color_nodes_col[xk] <- colSet2[k]
    
  }
  
  
  PlotNet <- ggnet2(net,node.color =  Z, palette = "Set2",label = TRUE,arrow.size = 6, arrow.gap = 0.017,label.size = 3) +  mythemeTransp
  
  
  
  PlotMatrix  <- plotMyMatrix(myMat,dimLabels = list(row = "species", col = "species"),plotOptions=list(title=myTitle,colNames = TRUE, rowNames = TRUE))+ mythemeTransp 
  PlotMatrix <- PlotMatrix + theme(axis.text.x = element_text(color=my_color_nodes_col),axis.text.y = element_text(color=my_color_nodes))
  return(list(PlotMatrix,PlotNet))#figure <- ggarrange(PlotNet , PlotMatrix ,ncol = 2, nrow = 1) 
  #return(figure)
}



