rm(list=ls())

setwd("~/Dropbox/bio386Khayla_Asha")
list.files()
list.files(, pattern="GSE18334")

#gene expression CV table
CV.tb = read.csv("GSE18334_log2CV.csv")  #CV table
CV.tb$ORF = as.character( CV.tb$ORF )
str(CV.tb)

#RLS table
RLS.tb = read.csv("lifespan.csv")
RLS.tb$ORF = as.character(RLS.tb$ORF)
str(RLS.tb)
summary(RLS.tb)

#fitness table
list.files(, pattern="fit")
fit.tb = read.csv("growth.fitness.hom.2012Nov13.csv")
#fit.tb = read.table("Regression_Tc1_hom.txt", sep='\t')
#fit.tb = read.table("Regression_Tc1_hom_2012Nov13.txt", sep="\t", header=T)
fit.tb$ORF = as.character(fit.tb$orf)
str(fit.tb)
fit.tb[fit.tb$ORF=="YDR001C",] #test

#protein network in pairs
PIN.pairs.tb = read.csv("pairs.csv", colClass=c("character","character")) #Protein network in pairs
str(PIN.pairs.tb)

#genetic network in pairs
GIN.pairs.tb = read.csv("sgadata_costanzo2009_stringentCutoff_101120.csv", header=F, 
                        colClass=c("character","character","character","character",NA,NA,NA))
names(GIN.pairs.tb) = c("ORF1", "name1", "ORF2", "name2", "epsilon", "stddev", "p") 
#epsilon = fij - fi * fj 
PositiveGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon>0, ]
NegativeGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon<0, ]
str(PositiveGIN.pairs.tb)
str(NegativeGIN.pairs.tb)

diff.Value = function( inpairs, inData ) { #inData must in "ORF', "Value" format
inpairs$Value1= inData$Value[match(inpairs$ORF1, inData$ORF)];
inpairs$Value2= inData$Value[match(inpairs$ORF2, inData$ORF)];
ret = mean( abs( inpairs$Value1 - inpairs$Value2 ), na.rm=T );
} 

call.diff.Value.in.permuted.pairs = function( orig.pairs, inData, inNsims) {
Nsims = inNsims; #number of permutations
permutated.diff.values = numeric( Nsims ); #empty vector to store calculations

  #merge two columns to a single columns for 'sample' function
ids.orig = c( orig.pairs$ORF1, orig.pairs$ORF2 ); 
len = length( ids.orig );
  
  # Now do N simulations
for( i in 1:Nsims ) {
newids = sample( ids.orig ); #permut
  
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

mylim = c(min(c(permutated.diff.CV, diff.CV.obs))*0.95, max(c(permutated.diff.CV, diff.CV.obs))*1.05 )
hist( permutated.diff.CV, xlim=c(mylim ), br=20, main="PIN" );
arrows( diff.CV.obs, 50, diff.CV.obs, 2, col="red" );
text( diff.CV.obs, 50.5, "obs");

PositiveGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon>0, c("ORF1", "ORF2") ]
str(PositiveGIN.pairs.tb)

inData = CV.tb[, c(1,4)]; names(inData) = c("ORF", "Value")

diff.CV.obs = diff.Value ( PositiveGIN.pairs.tb, inData );


permutated.diff.CV = call.diff.Value.in.permuted.pairs( PositiveGIN.pairs.tb, inData, 1000)
hist(permutated.diff.CV)


sub = permutated.diff.CV[ permutated.diff.CV < diff.CV.obs ]
p.Value = length(sub) / length(permutated.diff.CV);



  ### RLS, CV
RLS.tb2 = merge(RLS.tb, CV.tb)

summary(lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myStddev))  #p=0.0332
summary(lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myMean)) #p=0.66

m = lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV) #p=0.05618, weak positive correlation
summary(m)
plot( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV  )
abline(m, col='red')
#large CV or Stddev -> more noisy -> less robust -> long RLS (?)

m = lm(RLS.tb2$RLS_Del.pooled  ~ RLS.tb2$myCV) #p=0.15
summary(m)

### CV, fitness
fit.tb2 = merge(CV.tb, fit.tb)

m=lm(fit.tb2$YPD ~ fit.tb2$myCV); summary(m)
plot( YPD ~ myCV, data=fit.tb2)
abline(m, col="red")
#weak posiive, large CV -> more noise, less robust -> hight fitness? 



### PIN
#PIN.pairs.tb = read.csv("pairs.csv", colClass=c("character","character")) #Protein network in pairs
str(PIN.pairs.tb)
longPids = c(PIN.pairs.tb$ORF1, PIN.pairs.tb$ORF2)
degree = table( longPids );
sum(degree); #check the counting result, the length of ids
length(longPids)

PIN = data.frame(degree);
PIN$ORF = as.character( PIN$longPids); #make sure gene names are "characters"
PIN = PIN[, c(3,2)]; #remove the first column of unecessary data
str(PIN)

#genetic network in pairs
#GIN.pairs.tb = read.csv("sgadata_costanzo2009_stringentCutoff_101120.csv", header=F, 
#                        colClass=c("character","character","character","character",NA,NA,NA))
#names(GIN.pairs.tb) = c("ORF1", "name1", "ORF2", "name2", "epislon", "stddev", "p") 
#epsilon = fij - fi * fj 
#PositiveGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon>0, ]
#NegativeGIN.pairs.tb = GIN.pairs.tb[GIN.pairs.tb$epsilon<0, ]
#str(PositiveGIN.pairs.tb)
#str(NegativeGIN.pairs.tb)

positiveGids = c(PositiveGIN.pairs.tb$ORF1, PositiveGIN.pairs.tb$ORF2 )
negativeGids = c(NegativeGIN.pairs.tb$ORF1, NegativeGIN.pairs.tb$ORF2 )

plusGIN = data.frame(table(positiveGids ) );
minusGIN = data.frame(table(negativeGids ) )
names(plusGIN) = c("ORF", "plusGDegree")
names(minusGIN) = c("ORF", "minusGDegree")
plusGIN$ORF = as.character( plusGIN$ORF )
minusGIN$ORF = as.character( minusGIN$ORF )
GIN = merge(plusGIN, minusGIN)

#merge GIN, PIN to RLS
str(GIN)
str(PIN)

RLS.tb2$pDegree =  PIN$Freq[ match(RLS.tb2$ORF, PIN$ORF) ]
RLS.tb2 = merge( RLS.tb2, GIN)
RLS.tb2 = merge(RLS.tb2, fit.tb)

summary(lm(RLS_Del_alpha ~ pDegree, data=RLS.tb2))
summary(lm(RLS_Del_alpha ~ plusGDegree, data=RLS.tb2))
summary(lm(RLS_Del_alpha ~ minusGDegree, data=RLS.tb2))
summary(lm(RLS_Del_alpha ~ myCV, data=RLS.tb2)) # p = 0.006, R2=0.03
summary(lm(RLS_Del_alpha ~ myCV + pDegree + plusGDegree + minusGDegree, data=RLS.tb2)) #p_myCV = 0.018
summary(lm( myCV ~ RLS_Del_alpha + pDegree + plusGDegree + minusGDegree, data=RLS.tb2)) #
summary(lm( RLS_Del_alpha ~ myCV + pDegree + YPE + + plusGDegree + minusGDegree, data=RLS.tb2))   
summary(lm( RLS_Del_alpha ~ YPD, data=RLS.tb2))   
summary(lm( RLS_Del_alpha ~ YPG, data=RLS.tb2))   
summary(lm( RLS_Del_alpha ~ YPE, data=RLS.tb2))   
summary(lm( RLS_Del_alpha ~ myCV  + YPD + YPE + YPG, data=RLS.tb2))   
summary(lm( myCV ~ RLS_Del_alpha  + YPD + YPE + YPG, data=RLS.tb2))   

summary(lm( myCV ~ RLS_Del_alpha  + YPD + YPE + YPG, data=RLS.tb2))   

sub = permutated.diff.CV[ permutated.diff.CV < diff.CV.obs ]
p.K = length(sub) / length(permutated.diff.CV);