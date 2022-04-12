#Removendo primers
#Usamos cutadapt para esta funcao. Download e instrucoes em: http://cutadapt.readthedocs.io/en/stable/index.html
#Esse script foi especialmente desenvolvido para uso apos a execucao da etapa de deteccao de primers
#Ajuste o caminho para o cutadapt na sua maquina
cutadapt <- "diretorio/cutadapt"
system2(cutadapt, args = "--version") #Executa comandos shell em R

#Agora criaremos os outputs com primers removidos
#Os primers precisam estar na orientacao correta (Fwd/Rev)
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
#Trima FWD e o complemento reverso de REV das R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
#Trima REV e o complemento reverso de FWD das R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
#Executa Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output
                             fnFs.filtN[i], fnRs.filtN[i])) # input
}