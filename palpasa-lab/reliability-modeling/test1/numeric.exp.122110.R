# Dec 21, 2010
# numerically study realibity function of a parallel system 
# Use the equations 
# Ref: Gavrilov & Gavrilova, 2006 Handbook of models of for human aging

rm(list=ls());
single.component.vialbility = function(mu, t) { exp( - mu * t) }
system.viability = function(n,mu,t) { 1 - (1 - exp(-mu*t)^n ) } 
#for n parallel components
system.m         = function(n,mu,t){
	E = exp(-mu*t)
	A = n*mu*E*(1-E)^(n-1)
	B = 1 - (1 - E)^n
	A/B; 
	}

mu = 1; #mean of exponetial distribution rexp( n, lambda)

ns = c(1, 2, 3, 4 ); #number of parrallel nodes
colors = c("black","red","blue","green","yellow","brown");

t = seq(0,7*mu, by= mu/50)
i = 4; 
 s = system.viability(ns[i], mu, t )
 m = system.m(ns[i],mu,t)
 plot (m ~ t, ty='l', col=colors[i]); 
mylabels = paste("n=",ns, sep='');
#legend( 5, 0.8, mylabels, col=colors, lty=1)
legend( max(t)/2, max(m)*0.8, mylabels, col=colors, lty=1)

for( i in 1:3 ) {
 s = system.viability(ns[i], mu, t )
 m = system.m(ns[i],mu,t)
   lines( m ~ t, col=colors[i]); 
}



