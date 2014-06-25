# Dec 14, 2009 
# 4 node butter model: the network is: 1-2, 1-3, 2-3, 2-4, 3-4
# 4 node fully connected model: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
# exponential distribution f(x) = lambda e^(- lambda x) for component

rm(list=ls());

n=4; #four nodes
pois_lambda = n/2; #mean of exponetial distribution rexp( n, lambda)

colors = c("black","red","blue","green","yellow","brown");

#simstep = lambda/ 50; 
popsize = 10^5; #population size
age = rep(NA, popsize);


#simulate and determine the system's age
 for(  p in 1: popsize ) {
  lambda = rpois( 4, pois_lambda );
  x = rexp( 4, lambda);
  age13 = min( x[1], x[3] );
  age24 = min( x[2], x[4] );
  age23 = min( x[2], x[3] );
  age[p] =  max( age13, age24, age23);
 }
 age = age[! is.na(age)]
 
 median.age = median( age );
mean.age = mean( age );

#MAXAGE = 20; 
#mybreak1 = seq(0, lambda*4, by=simstep);
#mybreak2 = seq( lambda*4, MAXAGE, by = 0.2);
#mybreaks = unique(c(mybreak1, mybreak2, max(age))); 
#mybreaks = mybreaks[ ! is.na(mybreaks)]
#collect = data.frame( matrix(nrow= length(mybreaks)-1, ncol=length(ns)+1) );

h = hist(age, br=100, main=paste("n=",n) );

tb = data.frame( cbind( h$mids, h$density, h$counts ) );
names(tb) = c("mids","density","counts");
tb$cumulative.counts = 0; #initialization
for( i in 2:length(tb[,1]) ) {
  tb$cumulative.counts[i] = tb$cumulative.counts[i-1] + tb$counts[i];
}

tb$s = 1 - tb$cumulative.counts / max(tb$cumulative.counts);
tb$m = tb$density / tb$s ;
plot( tb$m ~ tb$mids, xlim = c(0.01, mean.age*4), main=paste( 'simulation of mortality rate', 'n=',n), log='xy' );
plot( tb$m ~ tb$mids, xlim = c(0.01, mean.age*4), main=paste( 'simulation of mortality rate', 'n=',n), log='y' );
 


