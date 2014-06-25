# Feb 2, 2010
# use normal distribution for exponential rates
# Plan:
# simulated the life span of each component
# the maximal age is the life span of the system
# repeat for all systems

rm(list=ls());

calculate.s.m = function( lifespan ){
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



n = 50; # numOfComponents
Npop = 1E3; # numOfSystems 

#mu.vec = rlnorm(n, mean=2, sd=2) #linear m ~ t
mu.vec = rlnorm(n, mean=0.2, sd=1)
#mu.vec = rlnorm(n, mean=0.2, sd=0.2) #linear m ~ t, smaller variance -> more linear
mu.vec = rlnorm(n, mean=0.2, sd=100)  #larger variance, linear log(m) ~ t!!!YES.

SystemAges = numeric(100);

for( i in 1: Npop){
 componentAges =  rexp(n, rate=mu.vec);
 summary(componentAges);
 SystemAges[i] = floor( max(componentAges) + 0.5 );   #extreme value distribution??
}


tb = calculate.s.m( SystemAges )
#sub = tb[1:floor(length(tb[,1])/4), ]
sub = tb[ tb$s > 0.10, ]
with( sub, plot( m ~ t))
with( sub, plot( s~ t ))

sub$m[sub$m==0] = NA; 
with( sub, plot (log(m) ~ t ))
m = with( sub, lm( log(m) ~ t) );
abline(m, col="red")
summary(lm( log(sub$m) ~ sub$t))
