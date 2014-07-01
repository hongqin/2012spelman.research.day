library(Heatplus)

data = read.csv("abcde-expression.csv", row.names=1);
data = as.matrix(data);

#show the data in heatmap without clustering
heatmap_2(data,scale="none", legend=4, do.dendro=c(F,F), col=RGBColVec(64),
Rowv= NA, Colv = NA, main = "microarray:5 genes 10 expts "
);
#ask the students which genes share similar expression patterns?


#do heatmap with clustering by genes (by rows)
heatmap_2(data,scale="none", legend=4, do.dendro=c(T,F), col=RGBColVec(64), Colv=NA,
main = "cluster by genes (rows)");

#do heatmap with two-way clustering
heatmap_2(data,scale="none", legend=4, do.dendro=c(T,T), col=RGBColVec(64),
main="Two way clustering");

