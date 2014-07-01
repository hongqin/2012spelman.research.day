# You are given the genetic interaction data in pairwise format, provided by Costanzo 2009 Science.
# Please use the following codes to load this genetic interaction network into R. 

geneticPairsTb = read.delim("sgadata_costanzo2009_stringentCutoff_101120.txt", header=F,)
names(geneticPairsTb) = c("ORF1", "name1", "ORF2", "name2", NA, NA, NA)


# Now, please calculate the number of interactions per gene.
geneticPairsTb$ORF1 = as.character( geneticPairsTb$ORF1 ); 
geneticPairsTb$ORF2 = as.character( geneticPairsTb$ORF2 );

ids = c(geneticPairsTb$ORF1, geneticPairsTb$ORF2);
gNet= data.frame(table( ids) );
head(gNet)
     ids Freq
1 YAL002W  222
2 YAL004W   10
3 YAL005C    8
4 YAL007C    3
5 YAL008W    3
6 YAL010C  176

# Please feel free to carry out more analysis as in our inclass analysis of Fraser02 Science paper. 

# growth fitness data
fitness = read.csv("growth.fitness.hom.csv");
fitness$orf = as.character( fitness$orf ); 

#evolutionary rates
Kdata = read.csv( "Sce.Spa.KaKs.csv")
Kdata$orframe = as.character(Kdata$orframe)

Kdata$gDegree = gNet$freq[ match(Kdata$orframe, gNet$ids)]

head(gNet)
      ids Freq
1 YAL002W  222
2 YAL004W   10
3 YAL005C    8
4 YAL007C    3
5 YAL008W    3
6 YAL010C  176

head(Kdata)
  orfname   Ka   Ks Omega
1 YAL005C 0.01 0.12  0.07
2 YAL007C 0.01 0.23  0.06
3 YAL008W   NA   NA    NA
4 YAL009W 0.02 0.22  0.09
5 YAL010C 0.05 0.26  0.18
6 YAL012W 0.01 0.16  0.04

 m = lm(Kdata$Ka ~ Kdata$gDegree);
 > summary(m)
 Multiple R-squared: 0.006554
 p-value: 0.0005679

plot( Kdata$gDegree ~ Kdata$Ka, main='genetic interaction ~ Ka')
abline(m, col="blue" );


# 

