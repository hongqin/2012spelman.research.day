
rm(list=ls())

debug = 0; 
require(flowCore);  require(flowClust); require(flowViz)

mydate = ''
mypatterns = c('BY', 'M5'); 
#mypath='2012JUNE28DHRDHE'
#mypath='2012JUNE28-BLANKS'
mypath='2012JUNE28-DHE'
mypath='2012JUNE28-DHR'

for( myp in mypatterns ) {
 fs = read.flowSet(path=mypath, pattern=myp)
 print(myp)
 pData(phenoData(fs))
# xyplot(`SSC-H` ~ `FSC-H`, data=fs);

 tf = transformList(from=colnames(fs), tfun=log10)
 fs2 = tf %on% fs
 rgate1 = rectangleGate('FSC-H'=c(2, 4), 'SSC-H'=c(2,4))
 my.filter = filter(fs2, rgate1)
 fs3 = Subset(fs2, my.filter)
 print(fs3)
 
 pdf(paste(mypath, myp, mydate, 'FL1-3scatter', 'pdf', sep='.'))
 xyplot(`FL3-H` ~ `FL1-H`, data=fs3, filter=rgate1);
 dev.off()

 pdf(paste(mypath, myp, mydate, 'FL1-3marginal', 'pdf', sep='.'))
 #densityplot(~ `FL1-H`, data=fs3, filter=curv1Filter("FL1-H")) 
 densityplot(~ ., fs3, channels=c("FL1-H", "FL3-H"), filter=list(curv1Filter("FL1-H"),  curv1Filter("FL3-H")))
 dev.off()
 
 pdf(paste(mypath, myp, mydate, 'SSC-FSC', 'pdf', sep='.'));
 xyplot(`SSC-H` ~ `FSC-H`, data=fs3, filter=rgate1);
 dev.off()
 
}


#quit("no")
