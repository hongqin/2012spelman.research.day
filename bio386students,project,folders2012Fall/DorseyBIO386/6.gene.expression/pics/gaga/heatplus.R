###### Generate heat maps based a set of  
### microRNA expression data from African-American prostate patients
### and Caucasian ones. 

 #install Heatplus pakcage for heatmap 
 source("http://bioconductor.org/biocLite.R"); 
 biocLite("Heatplus");

 library(Heatplus); #load Heatplus package into R
 library(adimpro)

 g = read.image('gaga.tif')
 m = matrix(as.numeric(g[[1]])[seq(2,1157100,2)], nrow=2436, ncol=2375)
  

 mymethod = 'ward';
 heatmap_2(m, do.dendro = c(TRUE, FALSE), legend = 1, 
  legfrac = 10,
  hclustfun = function(c) hclust( c, method= mymethod),
  col = gray(0:100/100),
  #main = "All experiments"
 )




 data=read.csv("yates.miRNA.csv",row.names=1); #read in microarray data
 data = t(data); #transpose the data
 
 heatmap_2(data); #A simple heat map
 heatmap_2(data, legend = 1 ) #Add a figure legend
 heatmap_2(data, legend = 1, col = RGBColVec(64) ) #customerize the color

############
############ END
############

