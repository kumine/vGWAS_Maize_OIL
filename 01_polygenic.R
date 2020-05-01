
#library(GenABEL)
library(MASS)

# setwd("")

# Because package ‘GenABEL’ was removed from the CRAN repository.
# So the code of 'polygenic' was extracted from https://github.com/cran/GenABEL
source("./src//polygenic.R")

pfile = "./data/sim_pheno.txt"
kfile = "./data/IBS_kinship.txt"

pheno <- read.table(pfile, head = T, check = F, row = 1)
k <- read.table(kfile, row = 1)
colnames(k) <- rownames(k)

inds <- intersect(colnames(k), rownames(pheno))


pheno <- pheno[inds, ]
k <- k[inds, inds]

out <- sapply(colnames(pheno), function(i) {
  res <- polygenic(y,
    kinship.matrix = k,
    data = data.frame(y = pheno[, i])
  )
  c(res$grresidualY)
})

dimnames(out) = dimnames(pheno)
out = data.frame(iid = rownames(out),out)

write.table(out,file = "./data/grammar_residuals.txt",
            sep = "\t",row.name = F,quote = F)











