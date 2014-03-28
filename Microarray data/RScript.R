library(GEOquery)
library(affy)
library(chicken.db)
library(annotate)
library(simpleaffy)
library(limma)
samples<-list.files(pattern="*.CEL")
filenames<-c("ROS1-_9.CEL", "ROS1+_10.CEL", "ROS2-_11.CEL", "ROS2+_12.CEL", "ROS3-_13.CEL", "ROS3+_14.CEL", "ROS4-_15.CEL", "ROS4+_16.CEL")
samplenames<-c("ROS1-_9", "ROS1+_10", "ROS2-_11", "ROS2+_12", "ROS3-_13", "ROS3+_14", "ROS4-_15", "ROS4+_16")
state<-c("Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated", "Untreated", "Treated")
phenodata<-as.data.frame(cbind(filenames, samplenames, state))
mydata<-ReadAffy(phenoData=phenodata)

mydata.rma<-rma(mydata)


#Checking differences between normalised and unnormalised data
library(RColorBrewer)
library(affyPLM)
cols<-brewer.pal(8, "Set1")

boxplot(mydata, col=cols)
boxplot(mydata.rma, col=cols)

hist(mydata, col=cols)
hist(mydata.rma, col=cols)

#Performing a quality control of the data using affyPLM package to check each sample individually
mydata.qc<-fitPLM(mydata)
image(mydata.qc, which=1, add.legend=TRUE)

#Removing control probesets, genes with low expression variance and/or uniformly expressed close to the background detection levels.
mydata.filtered<-nsFilter(mydata.rma, require.entrez=FALSE, remove.dupEntrez=FALSE)

#Checking what got removed
mydata.filtered$filter_log

#Nothing got removed

#Looking for differentially expressed genes using Limma
samples<-mydata.rma$state
samples<-as.factor(samples)
samples

design<-model.matrix(~0 + samples)
colnames(design)<-c("Untreated", "Treated")
eset<-exprs(mydata.rma)
fit<-lmFit(exprs(mydata.filtered$eset), design)
contrast.matrix<-makeContrasts(treated_untreated = Treated - Untreated, levels=design)
treated_fits<-contrasts.fit(fit, contrast.matrix)
treated_ebFit<-eBayes(treated_fits)
topTable(treated_ebFit, number=10, coef=1)
nrow(topTable(treated_ebFit, coef=1, number=1000, lfc=2))
probeset.list<-topTable(treated_ebFit, coef=1, number=20, lfc=2)

annot<-as.data.frame(select(chicken.db, rownames(probeset.list), c("ENSEMBL", "SYMBOL")))
colnames(annot)<-c("PROBEID", "ENSEMBLID", "SYMBOL")

results<-merge(annot, probeset.list, by.x="PROBEID", by.y="row.names")
head(results)
write.table(results, "results.txt", sep="\t", quote=FALSE)
return(results)

