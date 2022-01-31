#Tutoriais consultados:
#https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html
#https://www.bioconductor.org/packages//2.7/bioc/vignettes/oligo/inst/doc/V5ExonGene.pdf


library(limma)
library(statmod)
library(biomaRt)
library(oligo)

#Settar o working directory para um pasta com todos os arquivos .CEL, o mapeamentos das amostras
#e o mapeamento de genes

#Ler o arquivo de mapeamento das amostras
targets<-readTargets("targets.txt")


#Então a leitura dos arquivos na pasta oligo
celFiles <- list.celfiles(getwd(), full.names=TRUE)
rawData <- read.celfiles(celFiles)

#Executar a normalização RMA
eset<-oligo::rma(rawData)

#Remover do environment os dados brutos
rm(rawData)

#Coletar os dados de anotação da Affymetrix
featureData(eset) <- getNetAffx(eset, type="transcript")


#Atribuir as colunas os nomes das amostras
rownames(targets) <- sampleNames(eset)
pData(eset)<-targets[2:ncol(targets)]
pData(eset)


#Remover as sondas que não tiverem anotação genênica
HasSymbol <- !is.na(eset@featureData@data$geneassignment)
eset <- eset[HasSymbol,]
dim(eset)

write.exprs(eset, file="exprs.txt", sep="\t")

#Selecionar no mapfile "target.txt" qual é a categoria que se deseja fazer
Exp <- factor(targets$Condition)

#Construir a matrix de design com experimentos como var dependente da var independente, cepas
design <- model.matrix(~0+Exp)
#***Não se esquecer do til ali em cima***

#pipeline de DEGs Analysis do Limma
fit <- lmFit(eset, design)
fit <- treat(fit,lfc=1, trend=TRUE, robust=TRUE)
results <- decideTests(fit,p.value = 0.05)

#Armazenar os resultados do objeto MArrayLM em uma tabela, sem limite de genes,
#com cutoff de p.value<0.05 e |log2-fold-change|1>
table<-topTreat(fit, number=Inf, adjust.method="BH",p.value=0.05,lfc=1)

#Selecionar só algumas colunas necessárias
table<-table[,c(8,18,21,22)]
colnames(table)<-c("GeneData","logFC","PValue","FDR")

#Pegar o EnsemblID de dentro da coluna gene assignment
ensemblIDFinder<-function(linha){
  geneTerms<-unlist(strsplit(linha[1]," //"))
  ensemblID<-grep("ENST",geneTerms,value=TRUE)[1]
  return(ensemblID)
}
table$EnsemblID<-apply(table,1,ensemblIDFinder)

#Remover a coluna original de gene assignment
table$GeneData<-NULL
row.names(table)<-table$EnsemblID
head(table)

#Salvar os resultados em tabela por tabulação
write.table(table,"Tabelatopzera.txt",sep="\t")

