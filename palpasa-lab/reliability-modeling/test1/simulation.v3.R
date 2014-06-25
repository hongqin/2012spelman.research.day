# Feb 2, 2010
# use log-normal distribution for exponential rates
# Plan:
# simulated the life span of each component
# the maximal component age is the life span of each unit
#the minimum unit lifespan is the life span of the system
# repeat for all systems

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


n = 50; # numOfElements
m = 5; #number of blocks
Npop = 1E3; # numOfSystems (individuals)

mymean = 0.1; 
mysds = c(0.1, 0.3, 0.5, 1, 1.3, 1.5);

#mu.vec strores the constant failure rates of elements

for ( mysd in mysds) { 
#mu.vec = rlnorm(n, mean=2, sd=2) #linear m ~ t
mu.vec = rlnorm(n*m, mean=0.2, sd=1) #sigmoidal s, and m
#mu.vec = rep(0.1, 1000)

mu.vec = rlnorm(n*m, mean=0.2, sd=1) #Element constant mortality rates
#mu.vec = rlnorm(n, mean=0.2, sd=0.2) #linear m ~ t, smaller variance -> more linear
#mu.vec = rlnorm(n, mean=0.1, sd=2)  #larger variance, linear log(m) ~ t!!!YES.
#mu.vec = rlnorm(n, mean=0.1, sd=2)  #larger variance, linear log(m) ~ t!!!YES.
#mu.vec = rlnorm(n, mean=1, sd=25)  #larger variance, linear log(m) ~ t!!!YES.
#mu.vec = rlnorm(n, mean=0.1, sd=5)  #larger variance, linear log(m) ~ t!!!YES.

#containers
BlockAges = numeric(m)
SystemAges = numeric(Npop); 
# large simulation runs (population size) is needed for presentaiton

#simulate every individual
for( i in 1: Npop){
 ElementAges =  rexp(n, rate=mu.vec);
 #summary(ElementAges);
 
 #maximum of elementAges -> BlockAges
 BlockAges[i] = round( max(ElementAges), 1 );   
 #extreme value distribution??
 # for smooth curves, use round( ?, digits=1)
 
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
}

#x axis is time, y axis is viability
#where log(sub$m) refers to failure rate graph
#where sub$s refers to viability graph
