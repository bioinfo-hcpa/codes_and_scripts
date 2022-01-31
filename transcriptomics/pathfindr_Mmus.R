library(pathfindR)


#Get gene ID, LFC and FDR
colnames(Corrected_signal)
Corrected_signal<-Corrected_signal[,c(5,2,4)]

#Get Gene Ontologies as dataFrame
colnames(GOtable)<-c("GO","EnsemblGenes")
GOtable<-as.data.frame(GOtable)

#Convert the table into a Large list of named char vector
listGOtable=list()
for (i in (1:nrow(GOtable))){
 genes=GOtable[i,2]
 myVector<-strsplit(genes,", ")
 names(myVector)<-GOtable[i,1]
 listGOtable<-base::append(listGOtable,myVector)
}


#Get GO descritpions
go_terms <- read.delim("~/Downloads/Biomart_geneSymbol/go_terms.mgi", header=FALSE, row.names = 2, stringsAsFactors = FALSE)
go_terms$V1<-NULL
Vecgo_terms <- setNames(as.character(go_terms$V3), rownames(go_terms))


path2SIF<-("~/Downloads/Biomart_geneSymbol/musculusPIN_noquotes.sif")
Corrected_signal<-na.omit(Corrected_signal)

output <- run_pathfindR(input = Corrected_signal,
                                gene_sets = "Custom",
                                custom_genes = listGOtable,
                                custom_descriptions = Vecgo_terms,
                                pin_name_path = path2SIF,
                                convert2alias=FALSE,
                                max_gset_size = Inf,
                                adj_method="none", enrichment_threshold=0.8, iterations=1,
                                score_quan_thr =0.5, sig_gene_thr = 0.5)

