# 2012, March 28, 
# Use log-normal distribution to generate exponential rates
# Plan:
# simulated the life span of each elelment
# the maximal element is the life span of each block
# the minimum block lifespan is the life span of the system
# repeat for all systems (individuals in a population)

rm(list=ls());

setwd("C:/Users/Jwill107/Desktop/Dropbox/JChristopherBio320-aging-modeling/1.current/n5000")


mymodel = 'poisson'
#mymodel = 'Fixed#Elements'
nIn = 1400; # nIn is the # of elements

calculate.s.m.v2 = function( lifespan, bin.size=10 ){ #lifespan is simulated data
  my.data = sort( lifespan[!is.na(lifespan)] );
  my.max = max(my.data)
  my.interval.size = bin.size
  deathFreq = c(NA, NA)
  for ( i in 1:round(length(my.data) / my.interval.size))  {
    sub = my.data[ ((i-1)* my.interval.size + 1) : (i* my.interval.size) ]
    deathFreq = rbind( deathFreq, c(mean(sub, na.rm=T), my.interval.size) )
  }
  deathFreq = data.frame(deathFreq)
  deathFreq = deathFreq[-1,] #remove NA
  #deathFreq       = table( my.data )
  deathCumulative = deathFreq
  for( i in 2:length(deathFreq[,1])) {
    deathCumulative[i, 2] = deathCumulative[i-1, 2] + deathCumulative[i, 2]                
  }
  tot = length(my.data)
  s = 1 - deathCumulative[,2]/tot
  currentLive  = tot - deathCumulative[,2]        
  m =  deathFreq[,2] / currentLive; 
  
  #list( s=s, t=unique(my.data));
  #ret = data.frame( cbind(s, m, unique(my.data)));
  ret = data.frame( cbind(s, m, deathFreq[,1]) );
  names(ret) = c("s", "m", "t");
  ret;
}

#test
#lifespan = rnorm(1000, mean=100, sd=10)
#ret = calculate.s.m.v2(lifespan)


m = 15;  # numOfBlocks in a system
Npop = 1E4; # numOfSystems (individuals) 
mymean = 0.1; 
#mysds = c(0.1, 0.3, 0.5, 1, 1.3, 1.5);
mysds = c(0.2);

mysd = 0.2
#for ( mysd in mysds) { #sd loop
#containers
SystemAges = numeric(Npop);#store the lifespan for all individuals 
BlockAges = numeric(m) #buffer for temporary storage

for( i in 1: Npop){  # i-th individual in the population
  
  n = rpois(m, nIn); # numOfElements for each individual
  # this can be changed to power-law distribution
  
  if (mymodel == 'Fixed#Elements' ) {
    n = rep( nIn, m )  #This is the fixed number of elements per units
  }
  
  for( j in 1:m) {#Block loop
    #mu.vec strores the constant failure rates of elements
    #We are not sure whether mu.vec should be inside of Nop loop
    #mu.vec = rlnorm(n[j], mean=0.005, sd=1) #Element constant mortality rates
    
    #In GG01 paper, mu.vec are the same. 
    mu.vec = rep(0.005, n[j])
    
    ElementAges =  rexp(n[j], rate=mu.vec);    
    
    #maximum of elelementAges -> BlackAges
    BlockAges[j] = round( max(ElementAges), 1 );  
  }
  SystemAges[i] = min(BlockAges)  
}

tb = calculate.s.m.v2( SystemAges )
#sub = tb[1:floor(length(tb[,1])/4), ]
sub = tb[ tb$s > 0.25, ]
#with( sub, plot( m ~ t))
#with( sub, plot( s~ t ))

sub$m[sub$m==0] = NA;  #remove zero for log operations
plot (log(sub$m) ~ sub$t, main = paste(mymodel, " sd=", mysd) ) #plot the mortality ~ age
#m = with( sub, lm( log(m) ~ t) ); # linear regression
#m = lm( log(sub$m[sub$t<7]) ~ sub$t[sub$t<7] )
#abline(m, col="red")
#summary(lm( log(sub$m) ~ sub$t))

plot( sub$s ~ sub$t, main=paste(mymodel, " sd=", mysd) )

m = lm( log(sub$m[1:30]) ~ sub$t[1:30] )
m
out = c(m$coef, nIn)
names(out) = c("m0", 'G', 'nIn')
write.csv(t(out), paste('n', nIn, '.csv', sep=''))  

plot( log(sub$m[1:30]) ~ sub$t[1:30])
abline( m , col='red')

#} #sd loop

#x axis is time, y axis is viability
#where log(sub$m) refers to failure rate graph
#where sub$s is viability (capacity for normal growth and development) graph
