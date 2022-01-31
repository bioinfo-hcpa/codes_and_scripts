#Tutorial para análise de DEGs de chip Agilent 4x44k, seguindo 
#pipeline disponibilizado no Limma UserGuide:
#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#Arquivo de anotação baixado na home da plataforma GPL13605 no GEO:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13605


library(limma)

#Leitura do targetfile, aceita apenas 1 filename por arquivo, portanto, deve-se garantir
#que o mesmo canal (cor) esteja em coerência com a condição da amostra

SDRF <- read.delim("targets.csv",sep = ",")

#Leitura dos dados brutos no working directory
x<-read.maimages(SDRF[,"Original_File"],source = "agilent")

#Normalização com RMA
y<-backgroundCorrect(x, normexp.method = "rma")
y<-normalizeBetweenArrays(y, method = "quantile")


#Fazer análise de DEGs
design <- modelMatrix(SDRF, ref="Cntl")
fit <- lmFit(y,design)
fit <- na.omit(fit)

#Ver a porcentagem de NA nas covariáveis pra ver se não influencia significativamente no resultado final
1-(length(na.omit(fit$Amean))/ length(fit$Amean))

#Atribuir as variaveis NA, o menor valor do dataset
fit$Amean[is.na(fit$Amean)] <- min(fit$Amean, na.rm = TRUE)

#Função treat do Limma com um na.omit lá no meio para dar certo
#Colocar os parametros aqui:
robust=TRUE
lfc<-1
trend=TRUE

###################################################################
sv<-squeezeVar(fit$sigma^2,fit$df.residual,covariate = fit$Amean,robust = robust, winsor.tail.p=c(0.05,0.1))
fit$df.prior <- sv$df.prior
fit$s2.prior <- sv$var.prior
fit$s2.post <- sv$var.post
df.total <- fit$df.residual + sv$df.prior
df.pooled <- sum(fit$df.residual,na.rm=TRUE)
df.total <- pmin(df.total,df.pooled)
fit$df.total <- df.total
lfc <- abs(lfc)
coefficients<-as.matrix(fit$coefficients)
acoef <- abs(coefficients)
se <- fit$stdev.unscaled*sqrt(fit$s2.post)
tstat.right <- (acoef-lfc)/se
tstat.left <- (acoef+lfc)/se
fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
fit$p.value <- pt(tstat.right, df=df.total,lower.tail=FALSE) + pt(tstat.left,df=df.total,lower.tail=FALSE)
tstat.right <- pmax(tstat.right,0)
####NA OMIT que faltava e cagava TUDO#####
coefficients<-na.omit(coefficients)
##########################################
fc.up <- (coefficients > lfc)
fc.down <- (coefficients < -lfc)
fit$t[fc.up] <- tstat.right[fc.up]
fit$t[fc.down] <- -tstat.right[fc.down]
fit$treat.lfc <- lfc
#######################################################
#Limpar o env
rm(se)
rm(sv)
rm(tstat.left)
rm(tstat.right)
rm(x)
rm(y)
rm(fc.down)
rm(fc.up)
rm(coefficients)
rm(acoef)

#Terminar pipeline com a função topTreat
table<-topTreat(fit, number=Inf, adjust.method="BH",p.value=0.05,lfc=1)


#Arquivo do Biomart contendo os identificadores de mRNA do Ensembl, EMA, RefSeq, e nome do gene
AnnotationFile <- read.delim("~/Downloads/GSE78889/martexport.txt",stringsAsFactors = FALSE)
colnames(AnnotationFile)
#[1] "Gene.stable.ID"                 "European.Nucleotide.Archive.ID" "RefSeq.mRNA.predicted.ID"      
#[4] "RefSeq.mRNA.ID"                 "Transcript.stable.ID"           "Gene.name"  

#Criar um vetor Booleano para ir colocando as linhas vazias
table$EnsemblId<-as.character.default(NA)
noEnsID<-is.na(table$EnsemblId)

#Matches de Transcript.stable.ID
table[noEnsID,]$EnsemblId <- AnnotationFile$Gene.stable.ID[match(table[noEnsID,]$SystematicName, AnnotationFile$Transcript.stable.ID)]

#Atualizar vetor
noEnsID<-is.na(table$EnsemblId)

#Matches de RefSeq.mRNA.ID
table[noEnsID,]$EnsemblId <- AnnotationFile$Gene.stable.ID[match(table[noEnsID,]$SystematicName, AnnotationFile$RefSeq.mRNA.ID)]

#Atualizar vetor
noEnsID<-is.na(table$EnsemblId)

#Matches de RefSeq.mRNA.predicted.ID
table[noEnsID,]$EnsemblId <- AnnotationFile$Gene.stable.ID[match(table[noEnsID,]$SystematicName, AnnotationFile$RefSeq.mRNA.predicted.ID)]

#Atualizar vetor
noEnsID<-is.na(table$EnsemblId)

#Matches de European.Nucleotide.Archive.ID
table[noEnsID,]$EnsemblId <- AnnotationFile$Gene.stable.ID[match(table[noEnsID,]$SystematicName, AnnotationFile$European.Nucleotide.Archive.ID)]

#Atualizar vetor
noEnsID<-is.na(table$EnsemblId)

#Matches de Gene.name
table[noEnsID,]$EnsemblId <- AnnotationFile$Gene.stable.ID[match(table[noEnsID,]$GeneName, AnnotationFile$Gene.name)]

#Atualizar vetor
noEnsID<-is.na(table$EnsemblId)

#Remover aquelas sondas sem mapeamento no Enseml E descrição da sonda e sondas de BG
Probename<-(table$ProbeName==table$GeneName)
unknownDes <- grepl("Unknown",table$Description)
emptyDes <- is.na(table$Description)
BGprobe1 <- grepl("DN",table$SystematicName)
BGprobe2 <- grepl("BU",table$SystematicName)
BGprobe3 <- grepl("TC",table$SystematicName)
LOCs <- grepl("LOC",table$GeneName)
libs <- grepl("Library",table$Description)
COs <- grepl("CO[[:digit:]]",table$GeneName,perl = TRUE)
CXs <- grepl("CX[[:digit:]]",table$GeneName,perl = TRUE)



table<-table[!(unknownDes | emptyDes| BGprobe1 | BGprobe2 | BGprobe3 | LOCs | libs | COs |CXs | Probename) && !is.na(table$EnsemblId),]
noEnsID<-is.na(table$EnsemblId)

#Cortar umas colunas...
table<-table[,c(7:16)]

#Salvar a tabela 
write.table(table,"Tabelatopzera.txt",sep="\t")

#Dar uma conferida nas linhas que restaram pra ver se tem algum gene symbol s/ match
View(table[noEnsID,])
