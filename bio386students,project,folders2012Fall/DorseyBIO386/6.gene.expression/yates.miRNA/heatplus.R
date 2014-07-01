###### Generate heat maps based a set of  
### microRNA expression data from African-American prostate patients
### and Caucasian ones. 

 #install Heatplus pakcage for heatmap 
 source("http://bioconductor.org/biocLite.R"); 
 biocLite("Heatplus");

 library(Heatplus); #load Heatplus package into R

 data=read.csv("yates.miRNA.csv",row.names=1); #read in microarray data
 data = t(data); #transpose the data
 
 heatmap_2(data); #A simple heat map
 heatmap_2(data, legend = 1 ) #Add a figure legend
 heatmap_2(data, legend = 1, col = RGBColVec(64) ) #customerize the color

############
############ END
############

