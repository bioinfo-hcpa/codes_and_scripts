#Analise basica de microbioma

#Carregando biblioteca
library(dada2); packageVersion("dada2")

#Sinalizando o caminho
path <- "caminho" #Diretorio dos fastq apos descompactacao
list.files(path)

#Amostras devem estar no formato "_R1_001.fastq" e "_R2_001.fastq"
#O que estiver a esquerda sera o nome da amostra
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Extraindo o nome das amostras
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Plotando e inspecionando a qualidade
#Fwd - ajustar o intervalo de amostra que deseja inspecionar
plotQualityProfile(fnFs[1:2])
#Rev - ajustar o intervalo de amostra que deseja inspecionar
plotQualityProfile(fnRs[1:2])

#Filtrar e aparar
#Sinalizando o caminho para os arquivos filtrados
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Filtrar e aparar
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,210), 
                     #trimLeft = c(FWDPRIMERLEN,REVPRIMERLEN), se for necessario remover primers
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # No Windows ajuste multithread=FALSE
head(out)

#Aprendendo as taxas de erro
#Fwd
errF <- learnErrors(filtFs, multithread=TRUE)
#Rev
errR <- learnErrors(filtRs, multithread=TRUE)

#Plotando erros
#Fwd
plotErrors(errF, nominalQ=TRUE)
#Rev
plotErrors(errR, nominalQ=TRUE)

#Dereplicacao
#Fwd
derepFs <- derepFastq(filtFs, verbose=TRUE)
#Rev
derepRs <- derepFastq(filtRs, verbose=TRUE)

#Nomeando as reads dereplicadas
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Inferencia das sequencias
#Fwd
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
#Rev
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#Inspecionando o objeto dada retornado
dadaFs[[1]]
dadaRs[[1]]

#Agrupando as sequencias pareadas
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]]) #Ajustar o intervalo de amostras a inspecionar

#Construindo a tabela de sequencias
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspecionando a distribuicao de tamanho das sequencias
table(nchar(getSequences(seqtab)))

#Removendo chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Inspecionando o percentual de sequencias restantes apos a remocao de chimeras
sum(seqtab.nochim)/sum(seqtab)

#Acompanhamento das reads atraves do pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# Se estiver processando uma unica amostra, remova o sapply: exemplo - substitua sapply(dadaFs, getN) por getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy
#Substituir o caminho pelo diretorio do arquivo do banco de dados
taxa <- assignTaxonomy(seqtab.nochim, "/home/Otavio/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) 

#Adicionar especies (se for o caso)
##Substituir o caminho pelo diretorio do arquivo do banco de dados
taxa <- addSpecies(taxa, "caminho")

#Inspecionando as assinaturas taxonomicas
taxa.print <- taxa # removendo o nome das sequencias para apresentacao
rownames(taxa.print) <- NULL
head(taxa.print)

#Sequencia da analise com phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

# Construindo um dataframe a partir da informacao dos arquivos do dada2
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(samples.out)
rownames(samdf) <- samples.out

# Criando uma 'ASV Table'
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# Renomeando as sequencias para melhor tabulacao
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# Convertendo o objeto PS em dataframe
library(metagMisc)
ASV_Table <- phyloseq_to_df(ps, addtax = T, addtot = T, addmaxrank = F, sorting = "abundance")

# Gerando arquivo CSV
write.table(ASV_Table, file = "caminho/nomedoarquivo",row.names = TRUE, dec = ",", sep = ";", quote = FALSE)

#Gerando legenda de sequencias
Legendadesequencias <- refseq(ps)
write.table(Legendadesequencias, file = "caminho/nomedoarquivo")

