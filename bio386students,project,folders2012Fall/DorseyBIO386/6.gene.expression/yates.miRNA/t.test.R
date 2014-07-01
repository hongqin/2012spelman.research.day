# pairwise t.test on a set of microRNA data

 #read in data
 #these values are log2 transformed
 data=read.csv("yates.miRNA.csv",row.names=1); 

 #find out sample names
 samples = colnames( data );
 samples [grep( "RC.77",samples) ]

 #pick only cell-line "RC.77"
 sub77 = data[, grep( "RC.77",samples) ];

 #reorder the columns by "Normal" and "Tumor"
 sample77names = colnames( sub77 );
 sub77 = sub77[, c(grep("Normal",sample77names), grep("Tumor", sample77names) )]
 sub77 = data.frame(sub77);

 
 # t.test on every probe through a for-loop
 for( i in 1:length(sub77[,1]) ) {
    wt    = t( sub77[i, 1:3] );
    tumor = t( sub77[i, 4:6] );
    tt = t.test( wt, tumor); 
    sub77$p[i] = tt$p.value;
	sub77$ratio[i] = mean(tumor) - mean(wt); 
 }

 ### probes with differetial expressions
 out001 = sub77[sub77$p<0.01, ] # 0.01 cutoff

 #volcano plot
 plot( -log10(sub77$p) ~ sub77$ratio ); 
 points( -log10(out001$p) ~ out001$ratio, col="red" ); #color significant ones
 
 #output your report
 write.csv(out001, "001.csv"); 
 
################################################# 
###              END                          ###
################################################# 
