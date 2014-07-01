rm(list=ls())

setwd("C:/Users/Jasmine/Dropbox/Bio386_Jasmine")
list.files()
#Answers
#[1] "desktop.ini"                                    
#[2] "growth.fitness.hom.csv"                                                
#[3] "GSE9985_log2CV.csv"                                                         
#[4] "lifespan.csv"                                   
#[5] "pairs.csv"                                      
#[6] "protein half life.csv"                          
#[7] "Sce.Spa.KaKs.csv"                               
#[8] "sgadata_costanzo2009_stringentCutoff_101120.csv"

list.files(, pattern="GSE9985")
#Answer
#[1]"GSE9985_log2CV.csv"  

#gene expression CV table
CV.tb = read.csv("GSE9985_log2CV.csv")  #CV table

CV.tb = read.csv("GSE9985_log2CV.csv")  

CV.tb$ORF = as.character( CV.tb$ORF )
str(CV.tb)
#Answers 
#'data.frame':  5667 obs. of  4 variables:
# $ ORF     : chr  "YKR009C" "YDL157C" "YER059W" "YGL059W" ...
#$ myStddev: num  0.0285 0.0639 0.0583 0.0668 0.0282 ...
#$ myMean  : num  3.63 3.9 3.7 3.66 3.85 ...
#$ myCV    : num  0.00785 0.01638 0.01574 0.01826 0.00734 ...
#RLS table
#RLS.tb = read.csv("lifespan.csv")
#RLS.tb$ORF = as.character(RLS.tb$ORF)
#str(RLS.tb)
#summary(RLS.tb)

#protein network in pairs
PIN.pairs.tb = read.csv("pairs.csv", colClass=c("character","character")) #Protein network in pairs
str(PIN.pairs.tb)
#Answers
#'data.frame':  13667 obs. of  2 variables:
# $ ORF1: chr  "YMR047C" "YOR128C" "YLL028W" "YMR292W" ...
#$ ORF2: chr  "YGR218W" "YCR067C" "YHR188C" "YDR414C" ...

#genetic network in pairs
GIN.pairs.tb = read.csv("sgadata_costanzo2009_stringentCutoff_101120.csv", header=F, 
                        colClass=c("character","character","character","character",NA,NA,NA))
names(GIN.pairs.tb) = c("ORF1", "name1", "ORF2", "name2", "epsilon", "stddev", "p") 
#epsilon = fij - fi * fj 
PositiveGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon>0, ]
NegativeGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon<0, ]
str(PositiveGIN.pairs.tb)
#Answers
#'data.frame':  6956 obs. of  7 variables:
# $ ORF1   : chr  "YAL059W" "YAL056W" "YAL041W_tsq412" "YAL041W_tsq412" ...
#$ name1  : chr  "ECM1" "GPB2" "cdc24-3" "cdc24-3" ...
#$ ORF2   : chr  "YDR359C" "YCR063W" "YBR231C" "YGL029W" ...
#$ name2  : chr  "EAF1" "BUD31" "SWC5" "CGR1" ...
#$ epsilon: num  0.241 0.193 0.258 0.252 0.16 ...
#$ stddev : num  0.1594 0.0166 0.0358 0.1232 0.0905 ...
#$ p      : num  4.63e-02 3.77e-17 2.25e-08 4.36e-02 4.90e-02 ...

str(NegativeGIN.pairs.tb)
#Answers
#'data.frame':  67631 obs. of  7 variables:
# $ ORF1   : chr  "YAL063C" "YAL063C" "YAL063C" "YAL063C" ...
#$ name1  : chr  "FLO9" "FLO9" "FLO9" "FLO9" ...
#$ ORF2   : chr  "YCR028C-A" "YGL024W" "YHR116W" "YLR218C" ...
#$ name2  : chr  "RIM1" "YGL024W" "COX23" "YLR218C" ...
#$ epsilon: num  -0.197 -0.131 -0.131 -0.138 -0.212 ...
#$ stddev : num  0.0273 0.0485 0.0364 0.0435 0.0248 ...
#$ p      : num  1.52e-08 7.20e-03 2.96e-03 8.19e-03 4.55e-10 7.03e-03 9.82e-04 1.80e-03 4.51e-02 6.45e-08 ...

#######################
# First, define a generic function to calculate difference in pairs of proteins, 
diff.Value = function( inpairs, inData ) { #inData must in "ORF', "Value" format
  inpairs$Value1= inData$Value[match(inpairs$ORF1, inData$ORF)];
  inpairs$Value2= inData$Value[match(inpairs$ORF2, inData$ORF)];
  ret = mean( abs( inpairs$Value1 - inpairs$Value2 ), na.rm=T );
} 

#Qin's first attempt for a generic permutation test function
call.diff.Value.in.permuted.pairs = function( orig.pairs, inData, inNsims) {
  Nsims = inNsims; #number of permutations
  permutated.diff.values = numeric( Nsims ); #empty vector to store calculations
  
  #merge two columns to a single columns for 'sample' function
  ids.orig = c( orig.pairs$ORF1, orig.pairs$ORF2 ); 
  len = length( ids.orig );
  
  # Now do N simulations
  for( i in 1:Nsims ) {
    newids = sample( ids.orig ); #permutation is done here. 
    # newids is a single column vector
    
    #now, reformat newids to two columns (random pairs)
    ORF1 = newids[1: (len/2)];   # split the long to two columns, 1st half
    ORF2 = newids[ (len/2+1) : len]; # the second hard
    
    #Convert them back spreadsheet, (the random network in pairs)
    new.pairs = data.frame( cbind(ORF1, ORF2) ); 
    
    #calculate delta.K for one random network
    permutated.diff.values[i] = diff.Value( new.pairs, inData ); 
  }
  return( permutated.diff.values ) 
} #end of permutation function


#################### calculate the observed difference in PIN 
#convert CV data to "ORF" and "Value"
inData = CV.tb[, c(1,4)]; names(inData) = c("ORF", "Value")
diff.CV.obs = diff.Value ( PIN.pairs.tb, inData );

permutated.diff.CV = call.diff.Value.in.permuted.pairs( PIN.pairs.tb, inData, 1000)
hist(permutated.diff.CV)

# calulate p-value
sub = permutated.diff.CV[ permutated.diff.CV < diff.CV.obs ]
p.Value = length(sub) / length(permutated.diff.CV);
if(p.Value == 0) { p.Value = 1 / length(permutated.diff.CV)}

# generate a figure
mylim = c(min(c(permutated.diff.CV, diff.CV.obs))*0.95, max(c(permutated.diff.CV, diff.CV.obs))*1.05 )
hist( permutated.diff.CV, xlim=c(mylim ), br=20, main="PIN" );
arrows( diff.CV.obs, 50, diff.CV.obs, 2, col="red" );
text( diff.CV.obs, 50.5, "obs");

#################### calculate the observed difference in Positive GIN
PositiveGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon>0, c("ORF1", "ORF2") ]
str(PositiveGIN.pairs.tb)

inData = CV.tb[, c(1,4)]; names(inData) = c("ORF", "Value")

diff.CV.obs = diff.Value ( PositiveGIN.pairs.tb, inData );

permutated.diff.CV = call.diff.Value.in.permuted.pairs( PositiveGIN.pairs.tb, inData, 1000)
hist(permutated.diff.CV)

# calulate p-value
sub = permutated.diff.CV[ permutated.diff.CV < diff.CV.obs ]
p.Value = length(sub) / length(permutated.diff.CV);

# generate a figure
mylim = c(min(c(permutated.diff.CV, diff.CV.obs))*0.95, max(c(permutated.diff.CV, diff.CV.obs))*1.05 )
hist( permutated.diff.CV, xlim=c(mylim ), br=20, main="PositiveGIN" );
arrows( diff.CV.obs, 50, diff.CV.obs, 2, col="red" );
text( diff.CV.obs, 50.5, "obs");


#################### calculate the observed difference in Negative GIN
NegativeGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon<0, c("ORF1", "ORF2") ]
str(PositiveGIN.pairs.tb)

inData = CV.tb[, c(1,4)]; names(inData) = c("ORF", "Value")

diff.CV.obs = diff.Value ( NegativeGIN.pairs.tb, inData );

permutated.diff.CV = call.diff.Value.in.permuted.pairs( NegativeGIN.pairs.tb, inData, 1000)
hist(permutated.diff.CV)

# calulate p-value
sub = permutated.diff.CV[ permutated.diff.CV < diff.CV.obs ]
p.Value = length(sub) / length(permutated.diff.CV);

# generate a figure
mylim = c(min(c(permutated.diff.CV, diff.CV.obs))*0.95, max(c(permutated.diff.CV, diff.CV.obs))*1.05 )
hist( permutated.diff.CV, xlim=c(mylim ), br=20, main="PIN" );
arrows( diff.CV.obs, 50, diff.CV.obs, 2, col="red" );
text( diff.CV.obs, 50.5, "obs");

