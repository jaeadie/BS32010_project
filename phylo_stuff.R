#install and load the required packages
install.packages("ape")
install.packages("phangorn")
install.packages("picante")
require(ape)
require(phangorn)
require(picante)
#read in the fasta file and assign it to the variable x
x <- read.dna("filename.fasta", format="fasta")
#calculate the distances for each of the sequences in the file using the dist.dna command and store the results in the variable d
d <- dist.dna(x)
#convert the results and store them in a matrix, then write this matrix to the distance.csv file using the write.table command, for ease of viewing 
write.table(as.matrix(d),"distance.csv")
#recalculate the distances using an appropriate mutation model, in this case the K80 model
d2<-dist.dna(x, model="K80")
#use the neighbour join algorithm to calculate the tree positions for each sequence in the file
tr.nj<-nj(d)
#plot the tree from the neighbour join algorithm 
plot(tr.nj)

dt.nj<-cophenetic(tr.nj)

dmat<-as.matrix(d)

nms<-rownames(dmat)

dt.nj<-dt.nj[nms, nms]

dt.nj<-as.dist(dt.nj)
#plot the residuals from the nj analysis
plot(dt.nj-d,ylab="residuals", cex=0.5,main="NJ")

abline(h=0,lty=3)

fit <- pml(tr.nj,as.phyDat(x))

fit=optim.pml(fit,T)

plot(fit)

set.seed(8)
#bootstrap the results 100 times to improve confidence in the position of each node in the tree
bs<-bootstrap.pml(fit,bs=100,optNni=T)
#plot the boostrapped tree in prefferred format, in this case the p format
treeBS<-plotBS(fit$tree, type"p", bs)

mt<-modelTest(as.phyDat(x),G=F,I=F)

