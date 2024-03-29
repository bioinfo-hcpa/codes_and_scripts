#Tutorial para análise de DEGs de chip Agilent, seguindo 
#pipeline disponibilizado no Limma UserGuide:
#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#Arquivo de anotação baixado na home da plataforma GPL13912 no GEO:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13912


library(limma)
library(lumiMouseAll.db)


Annotation <- read.delim("~/Downloads/GSE30657/MouseRef-8_V2_0_R3_11278551_A.txt")

#Passar a coluna de probes negatias
GSE30657 <- read.delim("~/Downloads/GSE30657/GSE30657_Full.txt", row.names = NULL)
GSE30657$EntrezID <-Annotation$Entrez_Gene_ID[match(GSE30657$ID_REF, Annotation$Array_Address_Id)]
write.table(GSE30657,"GSE30657_full_w_Ensembl.txt",sep = "\t")

#Leitura do targetfile
SDRF <- readTargets("targets.txt",sep = "\t", row.names = NULL)

#Leitura dos dados brutos no working directory
x<-read.ilmn(files = "Dilatation_cond.csv")


#Remoção de BG
y<-backgroundCorrect(x, normexp.method = "saddle")
y<-normalizeBetweenArrays(x, method = "quantile")



#Fazer análise de DEGs
Exp <- factor(SDRF$Condition)
design <- model.matrix(~0+Exp)

fit <- lmFit(y,design)


#Ver a porcentagem de NA nas covariáveis pra ver se não influencia significativamente no resultado final
napercent<-(1-(length(na.omit(fit$Amean))/ length(fit$Amean)))
napercent

#Ver se valores infinitos negativos e positivos

isinfneg<-is.infinite(min(fit$Amean))
isinfneg
neginfpercent<- (length(fit$Amean[is.infinite(fit$Amean)])/length(fit$Amean))
neginfpercent

isinfmax<-is.infinite(max(fit$Amean))

if(isinfneg && !isinfmax){
  minNoInfinite<-min(fit$Amean[is.finite(fit$Amean)],na.rm = TRUE)
  fit$Amean[is.infinite(fit$Amean)]<-minNoInfinite
  }

if(napercent==0 && (length(fit$coefficients[is.na(fit$coefficients)])==0)){
  #Caso não tenha NAs, segue a função normal treat do pipeline do Limma
  fit<-treat(fit,robust=TRUE,lfc=1.5,trend=TRUE)
  } else{
  fit <- na.omit(fit)
  
  
  #Atribuir as variaveis NA, o menor valor do dataset
  fit$Amean[is.na(fit$Amean)] <- min(fit$Amean, na.rm = TRUE)
  
  #Função treat do Limma com um na.omit lá no meio para dar certo
  #Colocar os parametros aqui:
  robust=TRUE
  lfc<-1.5
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
  }


table<-topTreat(fit, number=Inf, adjust.method="BH",p.value=0.05,lfc=1.5)


#Matches de Gene.name
Mart <- read.delim("~/Downloads/GSE30657/Entrez+Ensembl.txt")
table$EnsemblId <- Mart$Gene.stable.ID[match(table$TargetID, Mart$NCBI.gene.ID)]

#Tirar as probes sem EnsemblID
noEnsembl<-is.na(table$EnsemblId)
table<-table[!noEnsembl,]


#Salvar a tabela 
write.table(table,"Tabelatopzera.txt",sep="\t")
# 
