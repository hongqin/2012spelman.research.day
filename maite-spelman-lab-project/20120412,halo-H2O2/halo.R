tb = read.csv("halo-Ecoli-20120412.csv")
t.test( tb[,2], tb[,3], pair=T, al='gr')

pdf("halo-H2O2-Ecoli20120412.pdf",width=5,height=5)
plot( tb$D0.4 ~ tb$H2O2, col='red', type='p', pch=19, log='x', ylim = c(0,50), xlab="H2O2", ylab="Diameter")
lines( tb$D0.4 ~ tb$H2O2, col='red', lty=2)
points( tb$D0.1 ~ tb$H2O2, col='blue', pch=18)
lines( tb$D0.1 ~ tb$H2O2, col="blue", lty=2)
legend(0.01, 45, legend=c('0.4% glucose', '0.1% glucose'), col=c("red","blue"),  pch=c(19,18))
text(0.1, 5, 'p=0.005, paired t-test')
dev.off()

