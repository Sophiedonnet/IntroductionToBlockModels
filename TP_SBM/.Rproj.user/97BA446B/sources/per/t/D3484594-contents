rm(list=ls())
library(sbm)
library(GGally)
library(network)
library(RColorBrewer)
library(ggpubr)



pathtodata <- paste0(getwd(),'/data')
my_files <- list.files(pathtodata)



food_webs <- lapply(
  X = seq_along(my_files),
  FUN = function(i) {
    df <- read.csv(file = paste0(pathtodata,"/",my_files[i]), header = TRUE, row.names = 1)
    A <- as.matrix(df)
    return(list(
      net = A,
      nr = nrow(A),
      nc = ncol(A),
      dens = mean(A),
      id = stringr::str_sub(my_files[i], 1, -5))
    )
  }
)

site_names <- c("Martins (M)", "Cooper (NC)", "Herlzier (NC)", "Venlaw (NC)",
                "Berwick (NZ)", "North Col. (NZ)", "Powder (NZ)", "Trib. C (NZ)" )
names(food_webs)<- site_names


PlotMatrix = list()
PlotNet = list()
PlotDeg = list()
mythemeTransp <- theme(
  panel.background = element_rect(fill='transparent'), #transparent panel bg
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  panel.grid.major = element_blank(), #remove major gridlines
  panel.grid.minor = element_blank(), #remove minor gridlines
  legend.background = element_rect(fill='transparent'), #transparent legend bg
  legend.box.background = element_rect(fill='transparent') #transparent legend panel
)



for (i in 1:length(site_names)){
  Li <- food_webs[[i]]
  Mat <- Li$net
  nr <- Li$nr
  U <- order(rownames(Mat))
  Mat <- Mat[U,U]
  net = network(Mat, directed = TRUE)
  
  PlotNet[[i]] = ggnet2(net,node.color = "#00BFC4", label = TRUE,arrow.size = 6, arrow.gap = 0.017,label.size = 3) +  mythemeTransp
  PlotMatrix[[i]] <- plotMyMatrix(Mat,dimLabels = list(row = "species", col = "species"),plotOptions=list(title=site_names[i],colNames = TRUE, rowNames = TRUE))+ mythemeTransp 
  
  degMat <- as.data.frame(rowSums(Mat))
  names(degMat) = c('Degree')
  degMat$type='Obs.'
  
  p = sum(Mat)/(nrow(Mat)^2-nrow(Mat))
  x <- seq(min(degMat$Degree),max(degMat$Degree),by=1)
  df <- with(degMat, data.frame(Degree = x, y = dbinom(x, nrow(Mat)-1,p)))
  
  PlotDeg[[i]]  <- ggplot(degMat)+ geom_histogram(aes(x = Degree, y = ..density..),binwidth = 1) +   mythemeTransp
  PlotDeg[[i]]  <- PlotDeg[[i]]  + geom_line(data = df, aes(x = Degree, y = y), color = "#00BFC4",linewidth=2) + ylab('Density')
  ggsave(PlotDeg[[i]], file=paste0("plots/degree_",gsub(" ", "", site_names[i]),".png"), width = 25,height = 30 , units = "cm" )
  ggsave(PlotNet[[i]], file=paste0("plots/network_",gsub(" ", "", site_names[i]),".png"),  width = 25,height = 30 , units = "cm")
  ggsave(PlotMatrix[[i]], file=paste0("plots/matrix_",gsub(" ", "", site_names[i]),".png"), width = 25 ,height = 30, units = "cm")
  
  
  
}
 


############################################## 
