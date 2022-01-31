#Pipeline para análise de microarranjo affymetrix, com normalização RMA, com base em infos disponíveis em:
#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#https://www.bioconductor.org/packages/devel/bioc/vignettes/affy/inst/doc/affy.pdf

#Dados do exemplo em:
#http://bioinf.wehi.edu.au/limma/data/ecoli-lrp.zip


library(affy)
library(limma)
library(statmod)
library(biomaRt)

#Settar o working directory para um pasta com todos os arquivos .CEL, o mapeamentos das amostras
#e o mapeamento de genes

#Ler o arquivo de mapeamento das amostras
targets<-readTargets("targets.txt")

#Criar um objeto affybatch com os arquivos na pasta
data<-ReadAffy()
#Executar a normalização RMA
eset<-rma(data)

#Atribuir as colunas os nomes das amostras
colnames(eset)<-row.names(targets)
head(exprs(eset))


#Selecionar no mapfile "target.txt" qual é a categoria que se deseja fazer
# a análise de DEGs, nesse caso será feita para experimentos e cepas, juntamente

Cond <- factor(targets$Condition, levels=c("Control","Mutant"))

#Construir a matrix de design com experimentos como var dependente da var independente, cepas
design <- model.matrix(~0+Cond)

#pipeline de DEGs Analysis do Limma
fit <- lmFit(eset, design)
fit <- treat(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit,p.value = 0.05, lfc = 1)

#Armazenar os resultados do objeto MArrayLM em uma tabela, sem limite de genes,
#com cutoff de p.value<0.05 e |log2-fold-change|1>
table<-topTreat(fit, number=Inf, adjust.method="BH",p.value=0.05, lfc=1)

ensembl = useMart(biomart= "ensembl", dataset = "mmusculus_gene_ensembl")
affy_ensembl= c("affy_mouse430a_2", "ensembl_gene_id")
probelist<-getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)
table$probes<-row.names(table)
table$EnsemblID <- probelist$ensembl_gene_id[match(table$probes, probelist$affy_mouse430a_2)]

table$EnsemblID<-na.omit(table$EnsemblID)

#Salvar os resultados em tabela por tabulação
write.table(table,"Tabelatopzera.txt",sep="\t")

