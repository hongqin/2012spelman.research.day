# You are given the genetic interaction data in pairwise format, Costanzo 2009 Science.
# Please use following codes to load the interaction data
tb = read.delim("sgadata_costanzo2009.txt", header=F,)

names(tb) = c("ORF1", "name1", "ORF2", "name2", NA, NA, NA)

# Now, please calculate the number of interactions per gene.


