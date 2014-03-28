setwd("~/BS32010/project/data")
files<-dir(pattern="*.CEL")
library(affy)
library(chicken.db)
library(annotate)
library(limma)
mydata<-ReadAffy()
eset_rma<-rma(mydata)
boxplot(mydata)
boxplot(exprs(eset_rma))

#Diff expression

design<-model.matrix(~0 + sampleType)