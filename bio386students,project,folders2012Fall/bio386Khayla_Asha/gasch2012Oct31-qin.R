# use long andn not-long categories

gasch = read.csv("gasch00.tab", sep="\t", header=T );
row.names(gasch) = gasch[,1]

sub.g = gasch[, seq(4,176)]
sub.g[ is.na(sub.g) ] = 0; #this means log(mut/wt)=0, i.e. no change of expression
row.names(sub.g) = gasch[,1]

#For research project, we need to calculate CV
gasch.var = apply( sub.g, 1, FUN=var)
gasch$stddev = sqrt(gasch.var)
gasch$mean = apply( sub.g, 1, FUN=mean)
gasch$CV = gasch$stddev / gasch$mean

summary(gasch$CV) #too much negative values due to log2 ratio
summary(gasch$stddev) #let use stddev for the time being. 

#####
###we need merge gasch$stddev with RLS, interaction data, fitness, and K data using match()

### RLS data
lifespan.tb = read.csv("lifespan.csv")

#### PPI
#input protein interaction pairs
pairs = read.csv("pairs.csv");
pairs$ORF1 = as.character( pairs$ORF1 );
pairs$ORF2 = as.character( pairs$ORF2 );

#input genetic interaction pairs
gpairs = read.delim("sgadata_costanzo2009_stringentCutoff_101120.txt", header=F,)
names(gpairs) = c("ORF1", "name1", "ORF2", "name2", NA, NA, NA)
gpairs$ORF1 = as.character( gpairs$ORF1 );
gpairs$ORF2 = as.character( gpairs$ORF2 );

################
#calculate connecting degrees for PIN
ids = c(pairs$ORF1, pairs$ORF2); 
PIN = data.frame(table( ids )); #for protein interaction network, 
PIN$ids = as.character( PIN$ids);
str(PIN);
names(PIN) = c("ids", 'pDegree')
str(PIN)

#Likewise, calculate connecting degrees for GIN
gids = c(gpairs$ORF1, gpairs$ORF2);
GIN = data.frame(table( gids ));
GIN$gids = as.character( GIN$gids );
str(GIN);
names(GIN) = c('ids', 'gDegree')
str(GIN)


#Now, use match to merge gDegree and pDegree
PIN$gDegree = GIN$gDegree[match(PIN$ids, GIN$ids) ]
