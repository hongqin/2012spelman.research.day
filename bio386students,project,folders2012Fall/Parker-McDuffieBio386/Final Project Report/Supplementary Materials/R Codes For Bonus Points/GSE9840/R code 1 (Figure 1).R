rm(list=ls())

setwd("C:/Users/Jasmine/Dropbox/Bio386_Jasmine")
list.files()
#Answers
#[1] "desktop.ini"                                    
#[2] "growth.fitness.hom.csv"                         
#[3] "GSE9840_log2CV.csv                                                                             
#[5] "lifespan.csv"                                   
#[6] "pairs.csv"                                      
#[7] "protein half life.csv"                          
#[8] "Sce.Spa.KaKs.csv"                               
#[9] "sgadata_costanzo2009_stringentCutoff_101120.csv"

list.files(, pattern="GSE9840")
#Answers
#"GSE98405_log2CV.csv"


CV.tb = read.csv("GSE9840_log2CV.csv")  #CV table
CV.tb$ORF = as.character( CV.tb$ORF )
str(CV.tb)
#Answers
#'data.frame':  5667 obs. of  4 variables:
# $ ORF     : chr  "YKR009C" "YDL157C" "YER059W" "YGL059W" ...
#$ myStddev: num  0.0285 0.0639 0.0583 0.0668 0.0282 ...
#$ myMean  : num  3.63 3.9 3.7 3.66 3.85 ...
#$ myCV    : num  0.00785 0.01638 0.01574 0.01826 0.00734 ...

#RLS table
RLS.tb = read.csv("lifespan.csv")
RLS.tb$ORF = as.character(RLS.tb$ORF)
str(RLS.tb)
#Answers
#'data.frame':  584 obs. of  17 variables:
# $ ORF           : chr  "YMR056C" "YCR107W" "YHR047C" "YAL054C" ...
#$ Gene          : Factor w/ 374 levels "","AAC1","AAD3",..: 2 3 4 5 6 7 8 9 10 11 ...
#$ RLS_Del_alpha : num  25.3 27.4 15.8 24.5 25.7 23.4 23.4 28.4 14.7 29.4 ...
#$ N             : int  -105 -20 -5 -10 -25 -10 -10 -15 -25 -10 ...
#$ BY4742        : num  23.9 23.8 22.4 26.1 28.4 22.7 29.6 25 25.8 30 ...
#$ N.1           : int  -135 -20 -5 -10 -25 -10 -10 -25 -25 -10 ...
#$ p.value       : num  2.2e-01 3.7e-01 1.9e-01 5.7e-01 4.0e-01 6.5e-01 8.1e-02 5.4e-01 9.5e-06 9.1e-01 ...
#$ RLS_Del_a     : num  28.7 0 0 0 0 0 0 0 0 0 ...
#$ N.2           : int  -100 0 0 0 0 0 0 0 0 0 ...
#$ BY4741        : num  26.8 0 0 0 0 0 0 0 0 0 ...
#$ N.3           : int  -120 0 0 0 0 0 0 0 0 0 ...
#$ p.value.1     : num  0.13 0 0 0 0 0 0 0 0 0 ...
#$ RLS_Del.pooled: num  26.9 27.4 15.8 24.5 25.7 23.4 23.4 28.4 14.7 29.4 ...
#$ N.4           : int  -205 -20 -5 -10 -25 -10 -10 -15 -25 -10 ...
#$ Pooled        : num  25.3 23.8 22.4 26.1 28.4 22.7 29.6 25 25.8 30 ...
#$ N.5           : int  -255 -20 -5 -10 -25 -10 -10 -25 -25 -10 ...
#$ p.value.2     : Factor w/ 176 levels "","1.00E-001",..: 49 70 29 108 77 120 147 103 169 163 ...

summary(RLS.tb)
#Answers
#ORF                 Gene     RLS_Del_alpha         N                BY4742     
#Length:584                :210   Min.   : 4.00   Min.   :-1810.00   Min.   :22.40  
#Class :character   NFT1   :  2   1st Qu.:21.80   1st Qu.:  -25.00   1st Qu.:24.70  
#Mode  :character   AAC1   :  1   Median :24.40   Median :  -15.00   Median :26.90  
#AAD3   :  1   Mean   :24.58   Mean   :  -24.06   Mean   :26.75  
#AAP1   :  1   3rd Qu.:27.80   3rd Qu.:   -5.00   3rd Qu.:29.00  
#ACS1   :  1   Max.   :38.80   Max.   :   -5.00   Max.   :32.80  
#(Other):368   NA's   :20      NA's   :20         NA's   :20     
#     N.1              p.value         RLS_Del_a          N.2               BY4741    
#Min.   :-1790.00   Min.   :0.0000   Min.   : 0.00   Min.   :-180.000   Min.   : 0.0  
#1st Qu.:  -30.00   1st Qu.:0.1100   1st Qu.: 0.00   1st Qu.:   0.000   1st Qu.: 0.0  
#Median :  -15.00   Median :0.3100   Median : 0.00   Median :   0.000   Median : 0.0  
#Mean   :  -29.25   Mean   :0.3729   Mean   : 2.93   Mean   :  -4.894   Mean   : 2.8  
#3rd Qu.:   -5.00   3rd Qu.:0.6025   3rd Qu.: 0.00   3rd Qu.:   0.000   3rd Qu.: 0.0  
#Max.   :   -5.00   Max.   :1.0000   Max.   :37.60   Max.   : 120.000   Max.   :31.0  
#NA's   :20         NA's   :20       NA's   :20      NA's   :20         NA's   :20    
#     N.3             p.value.1       RLS_Del.pooled       N.4            Pooled     
#Min.   :-220.000   Min.   :0.00000   Min.   : 4.00   Min.   :-1950   Min.   :22.40  
#1st Qu.:   0.000   1st Qu.:0.00000   1st Qu.:21.80   1st Qu.:  -25   1st Qu.:24.70  
#Median :   0.000   Median :0.00000   Median :24.30   Median :  -15   Median :26.90  
#Mean   :  -6.348   Mean   :0.02939   Mean   :24.25   Mean   :  -27   Mean   :26.73  
#3rd Qu.:   0.000   3rd Qu.:0.00000   3rd Qu.:27.50   3rd Qu.:   -5   3rd Qu.:29.00  
#Max.   : 100.000   Max.   :0.99000   Max.   :37.10   Max.   :  575   Max.   :32.80  
#NA's   :20         NA's   :20        NA's   :20      NA's   :20      NA's   :20     
#     N.5               p.value.2  
#Min.   :-1890.00            : 20  
#1st Qu.:  -35.00   3.10E-001: 12  
#Median :  -15.00   6.90E-001: 12  
#Mean   :  -33.32   1.00E+000: 11  
#3rd Qu.:   -5.00   2.90E-001: 11  
#Max.   :  555.00   4.30E-001: 11  
#NA's   :20         (Other)  :507  

#fitness table
list.files(, pattern="fitness")
fitness.tb = read.csv("growth.fitness.hoM.csv")
#fit.tb = read.table("Regression_Tc1_hom.txt", sep='\t')
#fit.tb = read.table("Regression_Tc1_hom_2012Nov13.txt", sep="\t", header=T)
fitness.tb$ORF = as.character(fitness.tb$orf)
str(fitness.tb)
#Answers
#'data.frame':  2706 obs. of  18 variables:
#  $ X             : int  1 2 3 4 5 6 7 8 9 10 ...
#$ gene          : Factor w/ 2706 levels "AAC1","aad14",..: 1426 1540 346 1556 1354 374 425 1332 678 1557 ...
#$ orf           : Factor w/ 2706 levels "YAL001C","YAL002W",..: 1 2 3 4 5 6 7 8 9 10 ...
#$ no..of.probes : int  2 2 0 2 2 2 2 2 2 2 ...
#$ measurments   : int  38 38 0 38 38 38 38 38 20 38 ...
#$ F.statistic   : num  15.51 28.07 NA 6.73 5.85 ...
#$ X.log.p.value.: num  7.54 10.43 NA 4.12 3.64 ...
#$ difference    : num  0.1145 0.0809 NA 0.0643 0.0488 ...
#$ YPD           : num  0.96 0.978 NA 0.967 1.03 ...
#$ YPDGE         : num  0.958 0.968 NA 0.98 1.031 ...
#$ YPG           : num  0.903 0.923 NA 1.002 1.024 ...
#$ YPE           : num  0.989 0.914 NA 1.014 1.054 ...
#$ YPL           : num  0.955 0.992 NA 1.001 1.01 ...
#$ YPD0.1DNP     : num  0.957 0.995 NA 1.015 1.021 ...
#$ YPD0.25DNP    : num  0.95 0.982 NA 1.031 1.047 ...
#$ YPG0.25DNP    : num  0.9 0.947 NA 1.007 1.034 ...
#$ YPG0.1DNP     : num  0.875 0.917 NA 1.001 1.005 ...
#$ ORF           : chr  "YAL001C" "YAL002W" "YAL003W"

fitness.tb[fitness.tb$ORF=="YDR001C",] #test
#Answer
#[1] X              gene           orf           
#[4] no..of.probes  measurments    F.statistic   
#[7] X.log.p.value. difference     YPD           
#[10] YPDGE          YPG          YPE           
#[13] YPL            YPD0.1DNP      YPD0.25DNP    
#[16] YPG0.25DNP     YPG0.1DNP      ORF           
#<0 rows> (or 0-length row.names)

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
#Answer
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

### RLS, CV
RLS.tb2 = merge(RLS.tb, CV.tb)

summary(lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myStddev))  #p=0.0332
#Answer
#Call:
#  lm(formula = RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myStddev)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-20.1839  -2.7656  -0.0109   3.1163  14.6760 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)        23.896      0.370  64.582   <2e-16 ***
# RLS.tb2$myStddev    1.020      0.542   1.882   0.0604 .  
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.099 on 495 degrees of freedom
#Multiple R-squared: 0.007104,  Adjusted R-squared: 0.005098 
#F-statistic: 3.542 on 1 and 495 DF,  p-value: 0.06043 



summary(lm(RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myMean)) #p=0.66
#Answer
#Call:
# lm(formula = RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myMean)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-20.3904  -2.6456  -0.1205   3.1921  14.3310 

#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    24.17278    1.28776  18.771   <2e-16 ***
# RLS.tb2$myMean  0.03003    0.14036   0.214    0.831    
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.117 on 495 degrees of freedom
#Multiple R-squared: 9.245e-05,  Adjusted R-squared: -0.001928 
#F-statistic: 0.04577 on 1 and 495 DF,  p-value: 0.8307 



plot( RLS.tb2$RLS_Del_alpha ~ RLS.tb2$myCV  )
abline(m, col='red')
large CV or Stddev -> more noisy -> less robust -> long RLS (?)

m = lm(RLS.tb2$RLS_Del.pooled  ~ RLS.tb2$myCV) #p=0.15
summary(m)
#Answer
#Call:
# lm(formula = RLS.tb2$RLS_Del.pooled ~ RLS.tb2$myCV)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-20.0112  -2.4591   0.2057   2.9851  12.9853 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    23.859      0.283  84.309   <2e-16 ***
# RLS.tb2$myCV    3.924      2.908   1.349    0.178    
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 4.663 on 495 degrees of freedom
#Multiple R-squared: 0.003665,  Adjusted R-squared: 0.001653 
#F-statistic: 1.821 on 1 and 495 DF,  p-value: 0.1778 


> ### CV, fitness
   fitness.tb2 = merge(CV.tb, fitness.tb) 
 m=lm(fitness.tb2$YPD ~ fitness.tb2$myCV); summary(m)
#Answer
#Call:
# lm(formula = fitness.tb2$YPD ~ fitness.tb2$myCV)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-0.89792  0.00269  0.03501  0.05034  0.12304 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      0.956793   0.003207 298.318  < 2e-16 ***
# fitness.tb2$myCV 0.260143   0.038039   6.839 1.07e-11 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 0.09846 on 1913 degrees of freedom
#(507 observations deleted due to missingness)
#Multiple R-squared: 0.02387,	Adjusted R-squared: 0.02336 
#F-statistic: 46.77 on 1 and 1913 DF,  p-value: 1.07e-11

plot( YPD ~ myCV, data=fitness.tb2)
abline(m, col="red")
#weak posiive, large CV -> more noise, less robust -> hight fitness? 

### PIN
#PIN.pairs.tb = read.csv("pairs.csv", colClass=c("character","character")) #Protein network in pairs
str(PIN.pairs.tb)
#Answer
#'data.frame':  13667 obs. of  2 variables:
# $ ORF1: chr  "YMR047C" "YOR128C" "YLL028W" "YMR292W" ...
#$ ORF2: chr  "YGR218W" "YCR067C" "YHR188C" "YDR414C" ...

longPids = c(PIN.pairs.tb$ORF1, PIN.pairs.tb$ORF2)
degree = table( longPids );
sum(degree); #check the counting result, the length of ids
#Answer
#27334

length(longPids)
#Answer
#27334

PIN = data.frame(degree);
PIN$ORF = as.character( PIN$longPids); #make sure gene names are "characters"
PIN = PIN[, c(3,2)]; #remove the first column of unecessary data
str(PIN)
#Answer 
'data.frame':  4478 obs. of  2 variables:
  $ ORF : chr  "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
$ Freq: int  1 3 3 2 7 24 5 4 3 5 ...

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
#Answer
#'data.frame':  2302 obs. of  3 variables:
#  $ ORF         : chr  "YAL002W" "YAL008W" "YAL010C" "YAL011W" ...
#$ plusGDegree : int  39 1 23 20 5 78 2 1 2 1 ...
#$ minusGDegree: int  183 2 153 169 25 368 6 5 25 85 ...
str(PIN)
#Answer
#'data.frame':  4478 obs. of  2 variables:
# $ ORF : chr  "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
#$ Freq: int  1 3 3 2 7 24 5 4 3 5 ...


RLS.tb2$pDegree =  PIN$Freq[ match(RLS.tb2$ORF, PIN$ORF) ]
RLS.tb2 = merge( RLS.tb2, GIN)
RLS.tb2 = merge(RLS.tb2, fitness.tb)

summary(lm(RLS_Del_alpha ~ pDegree, data=RLS.tb2))
#Answer
#Residuals:
# Min       1Q   Median       3Q      Max 
#-15.3851  -2.4257   0.1337   3.2026  14.2639 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept) 25.43987    0.65921  38.591   <2e-16
#pDegree     -0.05472    0.03625  -1.509    0.135

#(Intercept) ***
# pDegree        
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.378 on 85 degrees of freedom
#(26 observations deleted due to missingness)
#Multiple R-squared: 0.02611,  Adjusted R-squared: 0.01465 
#F-statistic: 2.278 on 1 and 85 DF,  p-value: 0.1349

summary(lm(RLS_Del_alpha ~ plusGDegree, data=RLS.tb2))
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ plusGDegree, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-13.9779  -2.8428   0.3526   3.0526  13.4758 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept) 25.65202    0.54741  46.861   <2e-16
#plusGDegree -0.10463    0.04884  -2.142   0.0344

#(Intercept) ***
# plusGDegree *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.178 on 111 degrees of freedom
#Multiple R-squared: 0.0397,  Adjusted R-squared: 0.03105 
#F-statistic: 4.589 on 1 and 111 DF,  p-value: 0.03436

summary(lm(RLS_Del_alpha ~ minusGDegree, data=RLS.tb2))
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ minusGDegree, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-13.6979  -2.8531   0.0536   3.3326  14.1572 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)  25.72689    0.64032  40.178   <2e-16
#minusGDegree -0.01610    0.01081  -1.489    0.139

#(Intercept)  ***
# minusGDegree    
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.232 on 111 degrees of freedom
#Multiple R-squared: 0.01959,  Adjusted R-squared: 0.01076 
#F-statistic: 2.218 on 1 and 111 DF,  p-value: 0.1392

summary(lm(RLS_Del_alpha ~ myCV, data=RLS.tb2)) # p = 0.006, R2=0.03
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ myCV, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-15.2573  -2.3414   0.0366   3.2313  13.3701 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)   23.647      0.982  24.080   <2e-16
#myCV          87.014     50.361   1.728   0.0868

#(Intercept) ***
# myCV        .  
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.214 on 111 degrees of freedom
#Multiple R-squared: 0.02619,  Adjusted R-squared: 0.01742 
#F-statistic: 2.985 on 1 and 111 DF,  p-value: 0.0868

summary(lm(RLS_Del_alpha ~ myCV + pDegree + plusGDegree + minusGDegree, data=RLS.tb2)) #p_myCV = 0.018
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ myCV + pDegree + plusGDegree + minusGDegree, 
#    data = RLS.tb2)

#Residuals:
# Min      1Q  Median      3Q     Max 
#-14.856  -2.318  -0.294   3.565  13.245 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)  24.78629    1.36663  18.137   <2e-16
#myCV         64.10426   62.01155   1.034    0.304
#pDegree      -0.05285    0.03632  -1.455    0.149
#plusGDegree  -0.09146    0.06366  -1.437    0.155
#minusGDegree  0.00317    0.01467   0.216    0.829

#(Intercept)  ***
# myCV            
#pDegree         
#plusGDegree     
#minusGDegree    
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.343 on 82 degrees of freedom
#(26 observations deleted due to missingness)
#Multiple R-squared: 0.0725,  Adjusted R-squared: 0.02726 
#F-statistic: 1.602 on 4 and 82 DF,  p-value: 0.1815 

summary(lm( myCV ~ RLS_Del_alpha + pDegree + plusGDegree + minusGDegree, data=RLS.tb2)) #
#Answers
#Call:
# lm(formula = myCV ~ RLS_Del_alpha + pDegree + plusGDegree + minusGDegree, 
#    data = RLS.tb2)

#Residuals:
# Min        1Q    Median        3Q       Max 
#-0.012123 -0.006048 -0.003277  0.003927  0.033963 

#Coefficients:
# Estimate Std. Error t value
#(Intercept)    1.230e-02  5.240e-03   2.347
#RLS_Del_alpha  2.007e-04  1.941e-04   1.034
#pDegree       -5.946e-05  6.475e-05  -0.918
#plusGDegree   -6.868e-05  1.138e-04  -0.604
#minusGDegree  -1.072e-05  2.593e-05  -0.413
#Pr(>|t|)  
#(Intercept)     0.0213 *
# RLS_Del_alpha   0.3043  
#pDegree         0.3612  
#plusGDegree     0.5478  
#minusGDegree    0.6805  
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 0.009454 on 82 degrees of freedom
#(26 observations deleted due to missingness)
#Multiple R-squared: 0.04469,  Adjusted R-squared: -0.001907 
#F-statistic: 0.9591 on 4 and 82 DF,  p-value: 0.4345

summary(lm( RLS_Del_alpha ~ myCV + pDegree + YPE + + plusGDegree + minusGDegree, data=RLS.tb2))   
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ myCV + pDegree + YPE + +plusGDegree + 
# minusGDegree, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-12.6224  -2.7163  -0.4115   3.3471  15.8481 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  12.70272    5.70746   2.226   0.0289 *
# myCV         61.94357   61.21760   1.012   0.3147  
#pDegree      -0.02462    0.03808  -0.647   0.5197  
#YPE          11.97489    5.47129   2.189   0.0315 *
# plusGDegree  -0.06851    0.06351  -1.079   0.2839  
#minusGDegree  0.01081    0.01486   0.728   0.4689  
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.253 on 80 degrees of freedom
#(27 observations deleted due to missingness)
#Multiple R-squared: 0.1251,  Adjusted R-squared: 0.07046 
#F-statistic: 2.289 on 5 and 80 DF,  p-value: 0.05356

summary(lm( RLS_Del_alpha ~ YPD, data=RLS.tb2))
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ YPD, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-15.1148  -2.9337   0.3343   3.2409  15.9587 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)   18.188      5.208   3.492  0.00069
#YPD            7.118      5.318   1.339  0.18345

#(Intercept) ***
# YPD            
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.264 on 110 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared: 0.01603,  Adjusted R-squared: 0.007084 
#F-statistic: 1.792 on 1 and 110 DF,  p-value: 0.1835 

summary(lm( RLS_Del_alpha ~ YPG, data=RLS.tb2))   
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ YPG, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-12.6672  -2.9525   0.0566   3.4112  13.9236 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
#(Intercept)   12.594      4.483   2.809  0.00588 **
# YPG           12.977      4.615   2.812  0.00583 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.125 on 110 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared: 0.06707,  Adjusted R-squared: 0.05859 
#F-statistic: 7.908 on 1 and 110 DF,  p-value: 0.005829 

summary(lm( RLS_Del_alpha ~ YPE, data=RLS.tb2)) 
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ YPE, data = RLS.tb2)

#Residuals:
# Min       1Q   Median       3Q      Max 
#-12.9610  -2.7789  -0.2885   3.1060  16.8448 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)   11.129      4.099   2.715 0.007699
#YPE           14.522      4.224   3.438 0.000828

#(Intercept) ** 
# YPE         ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.042 on 110 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared: 0.09704,  Adjusted R-squared: 0.08883 
#F-statistic: 11.82 on 1 and 110 DF,  p-value: 0.0008277 

summary(lm( RLS_Del_alpha ~ myCV  + YPD + YPE + YPG, data=RLS.tb2)) 
#Answers
#Call:
# lm(formula = RLS_Del_alpha ~ myCV + YPD + YPE + YPG, data = RLS.tb2)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-12.5549  -3.3928  -0.2098   3.2868  15.4849 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
#(Intercept)   16.170      5.477   2.953  0.00387 **
# myCV          64.526     49.141   1.313  0.19197   
#YPD           -9.564      6.988  -1.369  0.17395   
#YPE           25.701     10.922   2.353  0.02045 * 
#  YPG           -7.856     10.440  -0.753  0.45336   
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

#Residual standard error: 5.022 on 107 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared: 0.1288,  Adjusted R-squared: 0.09627 
#F-statistic: 3.956 on 4 and 107 DF,  p-value: 0.004927 

summary(lm( myCV ~ RLS_Del_alpha  + YPD + YPE + YPG, data=RLS.tb2))   
#Answers
#Call:
# lm(formula = myCV ~ RLS_Del_alpha + YPD + YPE + YPG, data = RLS.tb2)

#Residuals:
# Min        1Q    Median        3Q       Max 
#-0.014244 -0.006651 -0.003133  0.004618  0.033261 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.0008009  0.0111149   0.072    0.943
#RLS_Del_alpha 0.0002458  0.0001872   1.313    0.192
#YPD           0.0015226  0.0137556   0.111    0.912
#YPE           0.0022343  0.0218601   0.102    0.919
#YPG           0.0065718  0.0204182   0.322    0.748

#Residual standard error: 0.0098 on 107 degrees of freedom
#(1 observation deleted due to missingness)
#Multiple R-squared: 0.03587,  Adjusted R-squared: -0.0001733 
#F-statistic: 0.9952 on 4 and 107 DF,  p-value: 0.4135 
