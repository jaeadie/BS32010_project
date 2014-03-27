#load required packages
require(DESeq2)
require(RCurl)
require(biomaRt)
require(stringr)
#read in the data file and assign it to a variable
geneCounts<-read.delim("RNAseqCounts.txt",
                       head=T,sep="\t",skip=1, row.names=1)


#count the expressions that are not zero
nonZeroCounts<-geneCounts[rowSums(geneCounts[,6:28])>0,6:28]

# samples are either plus or minus. 
treatments <- as.factor(str_sub(colnames(nonZeroCounts),-9,-5))
#perform DESe1 analysis on the matrix
dds <- DESeqDataSetFromMatrix(as.matrix(nonZeroCounts),
                              as.data.frame(treatments),
                              design=~treatments)

dds$treatments <- relevel(dds$treatments,"minus")
#plot results in a file
dds <- DESeq(dds)
plotMA(dds, main="DESeq2")
DEresults <- results(dds)
write.csv(as.data.frame(DEresults), file="treatments_treated_results.csv")

# get annotations
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                host="www.ensembl.org", 
                path="/biomart/martservice")
#select which ensembl genome to use
gg4 <- useDataset("ggallus_gene_ensembl",mart=mart)
#select the appropriate annotations 
annot <- getBM(attributes=c("ensembl_gene_id",
                            "external_gene_id",
                            "affy_chicken"),
               filter="ensembl_gene_id",
               values=row.names(DEresults),
               mart=gg4)
#annotate the results using annot function
annotResults <- merge(annot,DEresults,by.x="ensembl_gene_id",by.y="row.names")

head(annotResults)
write.table(annotResults, "annotResults.txt", sep="\t", quote=FALSE)
return(annotResults)
#sort the results and display the top expressed genes based on the adjusted P value.
sortLimma<-head(annotResults[order(annotResults$padj),])
