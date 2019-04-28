# Hedenfalk et al. data
# Requires package LBE

# source("https://bioconductor.org/biocLite.R")
#  biocLite("LBE")

library(LBE)
data("hedenfalk.pval")
bc.bonf = bonffun(hedenfalk.pval)
bc.by = byfun(hedenfalk.pval)
bc.newbh = wavebhfun(hedenfalk.pval)
bc.newabh = waveadabhfun(hedenfalk.pval)


data("golub.pval")
all.bonf = bonffun(golub.pval)
all.by = byfun(golub.pval)
all.newbh = wavebhfun(golub.pval)
all.newabh = waveadabhfun(golub.pval)

# Compare with regular Benjamini-Hochberg
# 
# Requires qvalue package, which can be called by
# source("https://bioconductor.org/biocLite.R")
#  biocLite("qvalue")

qvalue(hedenfalk.pval,fdr.level=0.05)
qvalue(golub.pval,fdr.level=0.05)
