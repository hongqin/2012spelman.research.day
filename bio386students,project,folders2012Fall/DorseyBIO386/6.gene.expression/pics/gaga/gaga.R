# source("http://bioconductor.org/biocLite.R"); 
# biocLite("Heatplus");

 library(Heatplus); #load Heatplus package into R
 library(adimpro)
 library(pixmap)

 g = read.pnm("gsmall.pbm")
 m = matrix(g@grey, g@size[1], g@size[2])
 #write.pnm(g, "test.pnm")
  
 m2 = m[sample(1:g@size[1], g@size[1]), ] 
 g2 = g; 
 g2@grey = m2
 write.pnm(g2, "test.pbm")
  
 h = hclust(dist(m2))
 m3 = m2[h$order, ]
 g3 = g;
 g3@grey = m3;
 write.pnm(g3, "hclust.pbm")  
  
  
 mymethod = 'ward';
 heatmap_2(m, do.dendro = c(FALSE, FALSE), legend = 1, 
  #legfrac = 10,
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

