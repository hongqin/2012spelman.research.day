rm(list=ls());
mh = read.csv("levy-mh.csv")

#generate row labels 
mh$label = paste( mh[,1], mh[,2], mh[,3])
mh$label
mh[10:12, ]

#change -1 to 0.01
for( r in 1:length(mh[,1])) { #row
  for ( c in 4:11) { #columns
    if (mh[r,c] == -1) { mh[r,c] = 0.0001 }
    if (mh[r,c] == 0) { mh[r,c] = 0.001 }
  }
} #row

data = mh[,4:11]
data[1:5,1:4]

row.names(data) = as.character(mh$label)

data = as.matrix(data);
tmp = as.vector(data)
summary(tmp)
#hist(as.numeric(data))

logdata = log10(data)
data.ori =data; 
data = logdata;

#by treatment
dt = dist(data)
ht = hclust(dt)
plot(ht)

#by genes
dg = dist(t(data))
hg = hclust(dg)
plot(hg)

library(Heatplus)

heatmap_2(data, scale="none", do.dendro = c(T, T), legend = 2, 
          legfrac = 10, 
          #hclustfun = function(c) hclust( c, method= 'average'),
          col = rev(RGBColVec(64))
          )


quit("no")

####
#### Qin only worked out the codes before this line
####



#do heatmap with clustering by genes (by rows)
heatmap_2(data,scale="none", legend=4, do.dendro=c(T,F), col=RGBColVec(64), Colv=NA,
          main = "cluster by genes (rows)");

#do heatmap with two-way clustering
heatmap_2(data,scale="none", legend=4, do.dendro=c(T,T), col=RGBColVec(64),
          main="Two way clustering");


