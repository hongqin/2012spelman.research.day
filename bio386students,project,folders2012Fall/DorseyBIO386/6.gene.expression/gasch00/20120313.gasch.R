# use long andn not-long categories

#library(e1071)
#library(pcurve);

gasch = read.csv("gasch00.tab", sep="\t", header=T );
row.names(gasch) = gasch[,1]

sub.g = gasch[, seq(4,176)]
sub.g[ is.na(sub.g) ] = 0; #this means log(mut/wt)=0, i.e. no change of expression
row.names(sub.g) = gasch[,1]

library(Heatplus); #load Heatplus package into R

data= t(sub.g) #transpose, because I want the experiments listed horizontally

#just for a quic demo
#data = t(sub.g[1:100,1:20]) #use a small subset for quick demo

#heatmap_2(data); # A simple heat map
#heatmap_2(data, legend = 1 ) #Add a figure legend
#heatmap_2(data, legend = 1, col = RGBColVec(64) ) #customerize the color
heatmap_2(data, legend = 1, col = rev(RGBColVec(64) )) #customerize the color again


