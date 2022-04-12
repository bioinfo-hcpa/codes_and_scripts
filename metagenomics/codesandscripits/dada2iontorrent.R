## Carregar o pacote dada2.
library(dada2); packageVersion("dada2")

## Definir a variável path
path <- "caminho_diretorio" # Mudar o caminho para o diretório contendo os arquivos fastq individuais.
list.files(path)

## Agora vamos ler os nomes dos arquivos fastq e realizar algumas maniputações de strings
#Corrigir a string para o formato adequado de nome de arquivo ou renomear os arquivos no diretorio
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE)) 
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspecionar a qualidade das reads
plotQualityProfile(fnFs[1:66]) # aqui ajusta o intervalo de amostras que devem ser avaliadas quanto a qualidade
## É interessante trimar alguns nucleotídeos finais adicionais aos indicados pelo plot de qualidade

## Filter and trim
## Designar os nomes dos arquivos fastq.gz filtrados
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # coloca os arquivos filtrados na pasta ‘filtered’
names(filtFs) <- sample.names # nomeia os filtrados

## Filtragem - truncLen - tamanho das sequências
out <- filterAndTrim(fnFs, filtFs, truncLen=260,maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, trimLeft=34)
head(out)

## Aprendendo as taxas de erro
errF <- learnErrors(filtFs, multithread=TRUE)

## Visualizando as estimativas das taxas de erro
plotErrors(errF, nominalQ=TRUE)

## Inferência de amostras
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, pool=TRUE)

## Inspecionando o objeto dada
dadaFs[[1]] # indica amostra ou amostras

## Agora vamos converter o dado como se ele estivesse sendo agrupado
mergers <- dadaFs

## Inspecionando o merger dataframe das amostras
head(mergers[1:4]) # ajustar o intervalo das amostras que deseja inspecionar

## Construindo a tabela ASV
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Inspecionar distribuição das sequence lengths
table(nchar(getSequences(seqtab)))

## Remover chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## Verificando o percentual de chimeras removido
sum(seqtab.nochim)/sum(seqtab)

## Revisando as reads através do pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

## Assign taxonomy
## Para continuar, arrume o caminho do diretório no script
taxa <- assignTaxonomy(seqtab.nochim, "diretorio/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, verbose=TRUE)

## Inspecionando o assign taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

## Sequencia da análise com phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

## Construindo um dataframe a partir da informação dos arquivos do dada2
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(samples.out)
rownames(samdf) <- samples.out

## Criando uma 'ASV Table'
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

## Renomeando as sequencias para melhor tabulação
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

## Convertendo o objeto PS em dataframe
library(metagMisc)
ASV_Table <- phyloseq_to_df(ps, addtax = T, addtot = T, addmaxrank = F, sorting = "abundance")

## Gerando arquivo CSV
write.table(ASV_Table, file = "diretorio/ASV_Table.csv", row.names = TRUE, dec = ",", sep = ";", quote = FALSE)

##Gerando legenda de sequencias
Legendadesequencias <- refseq(ps)
write.table(Legendadesequencias, file = "diretorio/Legendadesequencias.csv")
