rm(list=ls())

setwd("~/Dropbox/bio386Khayla_Asha")
list.files()
list.files(, pattern="GSE")

#gene expression CV table
CV.tb = read.csv("GSE18334_log2CV.csv")  #CV table
#3821: wildtype, Study of the short term (within the first 330 seconds) transcriptional response of S.cerevisiae upon a sudden addition of glucose.
#Keywords: glucose pulse, chemostat culture, glucose catabolite repression

CV.tb = read.csv("GSE18334_log2CV.csv")  #Khayla, small p in +GIN, PIN?
#18332: Profiling oxidative stress in Saccharomyces cerevisiae by oxidation of transcription factor Swi6p
#mutants of Swi6 ?

CV.tb = read.csv("GSE12221_log2CV.csv")  #Jasmine, small p in -GIN, PIN
#12221: Decay profiles of Saccharomyces cerevisiae mRNAs following oxidative stress and DNA damage
#Wildtype strains ?
#oxidative stress (0.3mM hydrogen peroxide) or DNA damage (0.1% methyl methanesulfonate)
#QIN: 12221 and 7645 give opposit patterns. So, removal of DNA stress from 12221 
# and reco comparison will be need to further clarification.


CV.tb = read.csv("GSE33276_log2CV.csv")  #Asha, p<0.001 PIN?, p=0.046 -GIN, 
#33276, wildyptes, Gene expression profiles of S. cerevisiae under heat stress

CV.tb = read.csv("GSE7645_log2CV.csv")  #jessica, small p PIN, p=0.999 +GIN, p<0.001 -GIN
#7645: CHP?, Expression data for Saccharomyces cerevisiae oxidative stress response
#strains


CV.tb$ORF = as.character( CV.tb$ORF )
str(CV.tb)

#RLS table
#RLS.tb = read.csv("lifespan.csv")
#RLS.tb$ORF = as.character(RLS.tb$ORF)
#str(RLS.tb)
#summary(RLS.tb)

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



