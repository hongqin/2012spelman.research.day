# Feb 2, 2010
# use possion distribution for exponential rates
# numerically study realibity function of a parallel system 
# Use the equations 
# Ref: Gavrilov & Gavrilova, 2006 Handbook of models of for human aging

rm(list=ls());
single.component.vialbility = function(mu, t) { exp( - mu * t) }
system.viability  = function(n,mu,t) { 1 - (1 - exp(-mu*t)^n ) } 
system.viability.vectorMu.single.t = function(mu.vector, t1) {  
      Etmp = 1 - exp(-mu.vector*t1)
      1 - exp( sum(log(Etmp)) );
      } 

system.viability.vectorMu.vectorT = function( mu.vec, t.vec) {  
      result = 1:length(t.vec);
      for( i in 1:length(t.vec)){
      	   result[i] = system.viability.vectorMu.single.t(mu.vec, t.vec[i]);
      	}
      result; 
      } 

n = 10; #10 node simulation
mu = rpois(n, 10); #mean of exponetial distribution rexp( n, lambda)

maxmu = max(mu);
t = seq(0, 7*maxmu, by= maxmu/50)


s = system.viability.vectorMu.vectorT(mu, t )




