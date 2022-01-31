library(edgeR)

GSE95224 <- read.delim("E:/Arquivos GEO/GSE95224_raw_counts/GSE95224.csv",header = TRUE,row.names = 1)
GSE95224 <- GSE95224[,2:5]
group <- factor(c("H","WT","H", "WT"))
y <- DGEList(counts = GSE95224, group = group)

y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmFit(y, design)
 lrt <- glmLRT(fit)
 topTags(lrt)

pvalue <- exactTest(y)
pvalue$table

write.csv(pvalue$table,"edgeR_95224_p-valor.csv")
write.csv(y,"edgeR_95224.csv")

GSE111906_new <- read.delim("~/Documents/GEO/GSE111906_new.csv", header = TRUE, sep = ",",row.names = 1)
View(GSE111906_new)

GSE111906_new <- GSE111906_new[,11:22]
View(GSE111906_new)


groupie <- factor(c("Cntl", "Cntl", "Cntl","S", "S", "S","HS","HS","HS","H","H","H"))
y <- DGEList(counts = GSE111906_new, group = groupie)
y <- calcNormFactors(y)
design <- model.matrix(~groupie)
y <- estimateDisp(y,design)


pvalue <- exactTest(y)


write.csv(pvalue$table,"edgeR_111906_p-valor.csv")
write.table(y,"edgeR_GSE111906.txt",na='NA',col.names = FALSE,row.names = FALSE)
