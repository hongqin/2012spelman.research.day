#2012 April 9, log transformation of H2O2
#2012 Feb 29, 

rm(list=ls())
library(nlme)
 
infile = "_h2o2loh.csv";
   
tb = read.csv(infile, colClasses=c("character",NA, NA, "character", rep("numeric",8 ), NA));
names(tb) = c("Strain", "OD600", "Dilution","Date","H2O2stock", "White", "Black", "halfBlack", "quarterBlack", "ThreeQBlack", "QQBlack", "Other", "Notes")
tb$H2O2 = tb$H2O2stock/2
tb$tot = tb$White + tb$Black + tb$halfBlack + tb$quarterBlack + tb$ThreeQBlack + tb$QQBlack + tb$Other
tb.ori = tb; 
tb$Dilution = tb$Dilution / tb$Dilution[1]

###### some manual curations here

######## normalize all data
mycolumns = c("White","Black","halfBlack", "quarterBlack","ThreeQBlack", "QQBlack", "Other","tot"); 
for ( j in mycolumns) {
  tb[,j] = tb[,j] * tb$Dilution
}
tb[,c( "H2O2", "Dilution", "White")]

####### find out means
H2O2 = sort( unique( tb$H2O2) )
#s = H2O2
tbm = data.frame(cbind(H2O2))
for ( i in 1:length(H2O2)) {
  c = H2O2[i]
  tmp = tb[ tb$H2O2==c, ]  
  tbm$tot[i] = mean(tmp$tot, na.rm=T)
  tbm$White[i] = mean(tmp$White, na.rm=T)
  tbm$Black[i] = mean(tmp$Black, na.rm=T) 
  tbm$halfBlack[i] = mean(tmp$halfBlack, na.rm=T)  
  tbm$quarterBlack[i] = mean(tmp$quarterBlack, na.rm=T)
  tbm$ThreeQBlack[i] = mean(tmp$ThreeQBlack, na.rm=T)  
  tbm$QQBlack[i] = mean(tmp$QQBlack, na.rm=T);  
}

###### some manual curations here
tbm$halfBlack[tbm$halfBlack==0 & tbm$H2O2>0] = NA;

###### calculate fractions
tbf = tbm; 
tbf$s = tbf$tot / max(tbf$tot)
for ( j in 3:8) {
  tbf[, j] = tbf[,j] / tbf$tot
}

#### Begin ***********
s = tbf$s
t = tbf$H2O2
#2012 April 9, log transformation
#t = log10(tbf$H2O2); t[1]= -3
#shift = 3
#t = t + shift
plot( s ~ t, type = 'l')

#mystart = list( M=0.001, I = 0.001, G = 0.2 )
#fm <- gnls( s ~ exp( -M*t + (I/G)* ( 1 - exp(G * t) )  ) ,
#             start = list( M=0.001, I = 0.001, G = 0.2 )
#            , verbose=T);
#summary( fm );

#source("/Users/hongqin/lib/R/lifespan.r")
GM.s = function( IGM, t ) {
  I = IGM[1]; G = IGM[2]; M = IGM[3]; 
  log_s = (I/G) *(1 - exp(G* t)) - M*t
  ret  = exp( log_s )
}

myErrors = function( IGM, t, sObs) {
  myEstimate = GM.s(IGM, t)
  # print (IGM); print(myEstimate)
  SumOfErrors = sum( (myEstimate - sObs)^2 ); 
}

#test
 mystart = c(0.5, 80, 3); R = 100; sd=c(0.005, 0.05, 0.005) #for debug
 test =  myErrors( mystart, t, s)
 s2 = GM.s(mystart, t)
 sum( ( s2 - s)^2 )
 test; 
 plot( s ~ t); points( s2 ~ t, col='red')

# c is H2O2 concentration, i.e. t in Gompertz model. 
MetropolisH2O2 <- function(start=IGM, sObs=s, c=t, Runs = 10, sd = c(0.005, 2, 1)) {
  # c = t; start = mystart; Runs = 5; sd = c(0.005, 0.05, 0.005); sObs=s
  parmcount <- length(start)
  sims <- matrix(NA, nrow=Runs, ncol = parmcount)
  p = rep(0, Runs)
  colnames(sims) <- names(start)
  
  oldError <-  myErrors( start, c, sObs) ##
  #sims[1,] <- c( start, oldError); 
  sims[1, ] = start; 
  p[1] = oldError; 
  accepts <- 0
  
  for (i in 2:Runs) {
    newError = NA; y = NA; 
    while( is.na(newError) ) {
      jump <- rnorm(parmcount, mean=0, sd=sd)
      y <- sims[i-1,1:3] + jump
      if (y[1] <=0) { y[1]=10^-15;}
      if (y[2] <=0) { y[2]=0.001; }
      if (y[3] <=0) { y[3]=0; }
      newError <- myErrors( y, c, sObs); 
    }
    
    if (log(runif(1)) <  (oldError - newError)/ 75 ) { # if error is decreasing
      #sims[i,] <- c( y, newlogalpha)
      sims[i, ] = y; 
      oldError <- newError
      accepts <- accepts + 1
      p[i] = newError; 
    } else {
      sims[i,] <- sims[i-1,]
      p[i] = p[i-1]; 
    }
  }
  cat('Accepted ',100*accepts/(Runs-1),'%\n')
  #plot( p)
  data.frame( cbind(sims, p) );
}

mystart = c(1E-8,18, 0); mysd = c(1,2,1); myRuns= 1E3
#myFit = MetropolisH2O2( mystart, sObs=s[1:8], t[1:8], Runs=1E5, sd=c(1, 1, 1))
myFit = MetropolisH2O2( mystart, sObs=s, t, Runs= myRuns, sd=mysd );
myFit
myBestIGM = t(myFit[ myFit$p == min(myFit$p), 1:3 ])

#fm <- gnls( s ~ exp( (I/G)* ( 1 - exp(G * t)) - M*t),start = list( I = 4, G = 23, M=2));
#fm <- gnls( s ~ exp( (I/G)* ( 1 - exp(G * t)) ),start = list( I = 0.001, G = 0.2));
fm <- gnls( s ~ exp( (I/G)* ( 1 - exp(G * t)) - M*t), start = list( I = myBestIGM[1], G = myBestIGM[2], M=myBestIGM[3] ) );

#tFit =  seq(0, max(t), by=0.005)
#sFit = GM.s( t(myBestIGM), tFit )

#for log-transformated
tFit =  seq(min(t), max(t), by=0.005)
#sFit = GM.s( t(myBestIGM), tFit )
sFit = GM.s( t(fm$coef), tFit) 

#sFit
plot( s ~ t, type='p', main=paste( tb$Strain[1], tb$Date[1],'I=', myBestIGM[1],'G=',myBestIGM[2],'M=',myBestIGM[3] ) ) 
points( sFit ~ tFit, col="red"); lines( sFit ~ tFit, col='red')
lines( s ~ t, col="blue")
myBestIGM
fm$coef

#plot(myFit$p)
