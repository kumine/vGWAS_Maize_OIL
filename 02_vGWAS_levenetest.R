library(data.table)

gfile <- "./data/geno.txt"
pfile <- "./data/grammar_residuals.txt"
nm <- "sim1"
ofile <- "./data/sim1"

source("src/levene.test.R")


## ------------------------start code-----------------------------------------------------------
data <- fread(gfile, sep = "\t", head = T, data.table = F)

map <- data[, 1:6]
geno <- as.matrix(data[, -c(1:6)])
# set het to missing
geno[geno == 1] <- NA


pheno <- read.table(pfile, head = T, row.names = 1)
id <- intersect(colnames(geno), rownames(pheno))

geno <- geno[, id]
pheno <- pheno[id, ]

##
y <- pheno[, nm]

p.value <- apply(geno, 1, function(x) {
  res <- try(levene.test(y, x), T)
  if (class(res) == "try-error") {
    NA
  } else {
    res
  }
})

res <- data.frame(map, p.value, stringsAsFactors = F)
res_gc <- vGWAS.gc(res)
write.table(res_gc, file = paste(ofile, "vGWAS.txt"), quote = F, sep = "\t", row.names = F)
