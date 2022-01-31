    library(clusterProfiler)
    library(dplyr)
    library(KEGGprofile)
    library(ggplot2)
    library(data.table)
    library(gplots)
    library(limma)
      
      
      #Input datafile
      Corrected_signal <- read.delim("./Corrected_signal.txt", stringsAsFactors=FALSE, header=TRUE)
      
      #OPT. Type the identifier of the dataset, additionaly, the default it's the name of the working directory
      CompName<-''
      
      species<-'rno'
      dotplotgenes<-30
      heatmapgenes<-50
      fontSize<-8
      
      ########FROM HERE, nothing else need to be checked
      if(CompName==''){
        CompName=tail(unlist(strsplit(getwd(),"/"))[-2],n='1')
      }
      
      #Table conataing or similar from biomaRt package:
      #"Gene.stable.ID"   "Gene.name"   "NCBI.gene.ID"              
    
      if(c('mmu')==c(species)){
        BioMart <- read.delim('https://github.com/ldiass/MPSbase/raw/master/BioMart/Mmusculus_geneID_Biomart.txt', stringsAsFactors=FALSE )
        library(org.Mm.eg.db)
        myOrgDb<-org.Mm.eg.db
      }else if (c('hsa')==c(species)){
        if(c(species)==c('hsa')){
          BioMart <- read.delim('https://github.com/ldiass/MPSbase/raw/master/BioMart/hsapiens_geneID_biomart.txt')
        library(org.Hs.eg.db)
        myOrgDb<-org.Hs.eg.db
        }
        }else if(c(species)==c('rno')){
            BioMart <- read.delim('https://github.com/ldiass/MPSbase/raw/master/BioMart/rnovergicus.txt')
            library(org.Rn.eg.db)
            myOrgDb<-org.Rn.eg.db
            bioCNames<-colnames(BioMart)
            bioCNames[4]<-'NCBI.gene.ID'
            colnames(BioMart)<-bioCNames
        }else{
          warning("Species not found")
        }
      
      
      #Check the caps of the column EnsemblId
      if(length(grep("EnsemblId", colnames(Corrected_signal)))>0){
        EnsemblPos=as.numeric(grep("EnsemblId",colnames(Corrected_signal)))
      }else{
        EnsemblPos=as.numeric(grep("EnsemblID",colnames(Corrected_signal)))
      }
      
      #Get the Entrez ID
      Corrected_signal$Entrez<-BioMart$NCBI.gene.ID[match(Corrected_signal[,EnsemblPos], BioMart$Gene.stable.ID)]
      
      #Create a vector of Entrez ID values
      genes<-Corrected_signal$Entrez[!(is.na(Corrected_signal$Entrez))]
      
      #Order the table by pvalue and cuted the repeated entrez ID
      Corrected_signal<-Corrected_signal[order(Corrected_signal$adj.P.Val),]
      genes<-as.data.frame(genes) %>% distinct(.keep_all = TRUE)
      
      #Run Keggprofile
      keggResult<-find_enriched_pathway(as.matrix(genes), species = species, download_latest = TRUE)
      keggResult$stastic<-keggResult$stastic[order(keggResult$stastic$pvalueAdj),]
      
      #Adding species to the ID
      keggResult$stastic$KeggID<-unlist(lapply(rownames(keggResult$stastic),function(x){paste(species,x,sep = "")}))
      
      #Add one blank column to allow the match of the columns with the enrichGO table
      keggResult$stastic$Blank1<- 'empty'
      
      #Create the column where the geneList will be stored
      keggResult$stastic$GeneList<- 'empty'
      
      #Add the gene list to the static table
      GOsets<-keggResult$detail
      for(name in names(GOsets)){
      keggResult$stastic$GeneList[grep(name,keggResult$stastic$KeggID)]<-(paste(GOsets[[name]],collapse="/ "))}
      
      #Changing the position of the columns to allow them match with the enrichGo result
      keggResult$stastic<-keggResult$stastic[, colnames(keggResult$stastic)[c(1:5,7,6,8,9)]]
      
      #Write the table
      write.table(keggResult$stastic, file = paste('../',CompName,"_keggGO.txt",sep = ""), sep = "\t")
      
      #Ordering the result by Padj
      mydata<-keggResult$stastic[order(keggResult$stastic$pvalueAdj),]
      
      #Converting the Padj to -log(Padj)
      mydata$TransP<-as.numeric(lapply(mydata$pvalueAdj,function(x){-log10(x)}))
      
      #Setting the infinite numbers as the other minimals
      InfiniteValues<-is.infinite(mydata$TransP)
      maxtransp<-max(mydata$TransP[!InfiniteValues])
      mydata$TransP[InfiniteValues]<-maxtransp
      
      #Write the dotplot
      pdfname<-paste(CompName,"_kegg.pdf",sep = "")
      pdf(pdfname)
      
      mydata$Pathway_Name <- factor(mydata$Pathway_Name, levels = mydata$Pathway_Name[order(mydata$TransP)])
      ggplot(mydata[1:dotplotgenes,], aes_string(x='Percentage', y="Pathway_Name", size='Gene_Found', color='TransP')) +
      geom_point() +
      scale_color_continuous(low="red", high="blue", name = '-log(FDR)', guide=guide_colorbar(reverse=TRUE)) +
      ylab(NULL) + theme(axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"))   + scale_size(range=c(3, 8))
      
      dev.off()
      
    ##############################################################################################
      
      Corrected_signal <- na.omit(Corrected_signal)
      
      ego <- enrichGO(gene = Corrected_signal[,EnsemblPos],
                      OrgDb = myOrgDb,
                      keyType = 'ENSEMBL',
                      ont = 'ALL',
                      pAdjustMethod = 'BH',
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.1)
      
      write.table(ego, file =paste("../",CompName, "_GO.csv",sep = ""), sep = "\t")
      
      pdf(paste(CompName,"_GO_ALL.pdf",sep = ""))
      clusterProfiler::dotplot(ego, showCategory = dotplotgenes, font.size = 8)
      dev.off()
      
      
        #Write the dot_plot for each GO classification
        egobp <- enrichGO(gene = Corrected_signal[,EnsemblPos],
                          OrgDb = myOrgDb,
                          keyType = 'ENSEMBL',
                          ont = 'BP',
                          pAdjustMethod = 'BH',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.1)
        
        pdf(paste(CompName,"_bp.pdf",sep = ""))
        dotplot(egobp, showCategory = dotplotgenes, font.size = 8)
        dev.off()
        
        
        egomf <- enrichGO(gene = Corrected_signal[,EnsemblPos],
                          OrgDb = myOrgDb,
                          keyType = 'ENSEMBL',
                          ont = 'MF',
                          pAdjustMethod = 'BH',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.1)
        
        pdf(paste(CompName,"_mf.pdf",sep = ""))
        dotplot(egomf, showCategory = dotplotgenes, font.size = 8)
        dev.off()
        
        
        egocc <- enrichGO(gene = Corrected_signal[,EnsemblPos],
                          OrgDb = myOrgDb,
                          keyType = 'ENSEMBL',
                          ont = 'CC',
                          pAdjustMethod = 'BH',
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.1)
        
      pdf(paste(CompName,"_cc.pdf",sep = ""))
      dotplot(egocc, showCategory = dotplotgenes, font.size = 8)
      dev.off()
      
    
      DEGfile <- read.delim( "Corrected_signal.txt")
      targets<-readTargets("targets.txt")
      expressfile<- read.delim("NormalizedExpression.csv" )
      NumberOfGene<-heatmapgenes
      
      head(expressfile)
      #Which column the data numbers start
      ColumnStart<-1
      EnsSamples<-(nrow(targets)+ColumnStart-1)
      
      DEGfile<-DEGfile[order(abs(DEGfile$logFC), decreasing = TRUE),]
      DEGfile<-DEGfile[1:(NumberOfGene*2),]
      
      exprs<-(expressfile[match(DEGfile$EnsemblId,expressfile$EnsemblID),])
      
      exprs$geneSymbol<-BioMart$Gene.name[match(exprs$EnsemblID,BioMart$Gene.stable.ID)]
      exprs<-exprs[!duplicated(exprs$geneSymbol),]
      exprs<-exprs[!is.na(exprs$geneSymbol),]
      rownames(exprs)<-exprs$geneSymbol
      
      exprs<-as.matrix(exprs[1:NumberOfGene,ColumnStart:EnsSamples])
      
      
      my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
      
      
      # creates a 5 x 5 inch image
      png("heatmaps.png",    # create PNG for the heat map        
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8) 
      
      heatmap.2(exprs,scale = "row",col=my_palette, trace = "none", density.info = "none", labCol = targets$Condition,key.title = "NA",
                keysize = 1.0, cexCol = 1.0, offsetCol = 0.3)
      
      
      dev.off()
      
      #Clear environment
      #rm(list = ls(all.names = TRUE))
        
        
        
