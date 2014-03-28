#We first require all the packages we will use.
require(ape)
require(phangorn)
require(picante)

#Create an object that contains our data file.
x<-read.dna("CLUSTALOMEGA_alignment_mRNA.fasta", format="fasta")

#Create a first matrix using the default K80 (Kimura 1980, a Jukes-Cantor model variation)
d<-dist.dna(x)
write.table(as.matrix(d), "distances_unfit.csv")

#Tree construction and plotting using the neighbour-joining, FastME and BIOINJ method.
tr.nj<-nj(d)
plot(tr.nj)


#Calculating tree distorsion for NJ
dt.nj<-cophenetic(tr.nj)
dmat<-as.matrix(d)
nms<-rownames(dmat)
dt.nj<-dt.nj[nms, nms]
dt.nj<-as.dist(dt.nj)
plot(dt.nj-d,ylab="Residuals", cex=0.5, main="NJ")
abline(h=0, col="red", lty=3)

#Now we will test which substitution method is better for our data.
mt<-modelTest(as.phyDat(x))

#Fitting the model
fittree<-pml(tr.nj, phyDat(x), k=4,inv=0.2)
fittree=optim.pml(fittree, k=4,inv=0.2)
plot(fittree)
set.seed(8)
bs<-bootstrap.pml(fittree,bs=100,optNni=T)
treeBS<-root(plotBS(fittree$tree, bs), 1)
