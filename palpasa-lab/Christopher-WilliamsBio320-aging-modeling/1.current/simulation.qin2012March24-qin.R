# 2012, March 20
# Use log-normal distribution to generate exponential rates
# Plan:
# simulated the life span of each elelment
# the maximal element is the life span of each block
# the minimum block lifespan is the life span of the system
# repeat for all systems (individuals in a population)

rm(list=ls());

calculate.s.m = function( lifespan ){ #lifespan is simulated data
   my.data = sort( lifespan[!is.na(lifespan)] );
   deathFreq       = table( my.data )
   deathCumulative = deathFreq
   for( i in 2:length(deathFreq)) {
        deathCumulative[i] = deathCumulative[i-1] + deathCumulative[i]                
   }
   tot = length(my.data)
   s = 1 - deathCumulative/tot
   currentLive  = tot - deathCumulative        
   m =  deathFreq / currentLive; 

   #list( s=s, t=unique(my.data));
   ret = data.frame( cbind(s, m, unique(my.data)));
   names(ret) = c("s", "m", "t");
   ret;
}

m = 10;  # numOfBlocks in a system
Npop = 1E4; # numOfSystems (individuals) 
mymean = 0.1; 
#mysds = c(0.1, 0.3, 0.5, 1, 1.3, 1.5);
mysds = c(0.2);

for ( mysd in mysds) { #sd loop

  #containers
  SystemAges = numeric(Npop);#store the lifespan for all individuals 
  BlockAges = numeric(m) #buffer for temporary storage
 
 for( i in 1: Npop){  # i-th individual in the population
   
   n = rpois(m, 50); # numOfElements for each individual
   # this can be changed to power-law distribution

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

tb = calculate.s.m( SystemAges )
#sub = tb[1:floor(length(tb[,1])/4), ]
sub = tb[ tb$s > 0.25, ]
#with( sub, plot( m ~ t))
#with( sub, plot( s~ t ))

sub$m[sub$m==0] = NA;  #remove zero for log operations
plot (log(sub$m) ~ sub$t, main = paste("sd=", mysd) ) #plot the mortality ~ age
#m = with( sub, lm( log(m) ~ t) ); # linear regression
#m = lm( log(sub$m[sub$t<7]) ~ sub$t[sub$t<7] )
#abline(m, col="red")
#summary(lm( log(sub$m) ~ sub$t))

plot( sub$s ~ sub$t, main=paste("sd=", mysd) )

} #sd loop

#x axis is time, y axis is viability
#where log(sub$m) refers to failure rate graph
#where sub$s refers to viability graph
