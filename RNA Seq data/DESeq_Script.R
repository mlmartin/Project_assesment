require(DESeq2)
require(RCurl)
require(stringr)
library(limma)
dat<-read.delim("/usr/local/share/BS32010/expression/data/RNAseqCounts.txt", sep="\t", skip=1, head=T, row.names=1)

nonZerodat<-dat[rowSums(dat[,7:25])>0,7:25]

treatments<-as.factor(str_sub(colnames(nonZerodat),-9,-5))

dds<-DESeqDataSetFromMatrix(as.matrix(nonZerodat), as.data.frame(treatments), design=~treatments) 

dds$treatments<-relevel(dds$treatments, "minus")

dds<-DESeq(dds)

require(biomaRt)

DEresults<-results(dds)

mart<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", path="/biomart/martservice")

gg4<-useDataset("ggallus_gene_ensembl", mart=mart)

annot<-getBM(attributes=c("ensembl_gene_id", "external_gene_id", "affy_chicken"), filters="ensembl_gene_id", values=rownames(DEresults), mart=gg4, uniqueRows=TRUE)

annotResults<-merge(annot, DEresults, by.x="ensembl_gene_id", by.y="row.names")

sortLimma<-head(annotResults[order(annotResults$padj),])

write.table(sortLimma, "DESeqresults.txt", sep="\t", quote=FALSE)
