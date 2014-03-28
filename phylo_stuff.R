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
#calculate the distortion
dt.nj<-cophenetic(tr.nj)
#pull out the taxa from (d) and assign them to dmat
dmat<-as.matrix(d)
#pull out the rownames of dmat and assign it to nms
nms<-rownames(dmat)
#update the cophenetic analysis with the rownames
dt.nj<-dt.nj[nms, nms]
#calcualte the new distances using as.dist
dt.nj<-as.dist(dt.nj)
#plot the residuals from the new distances with the originals subtracted to get the residuals
plot(dt.nj-d,ylab="residuals", cex=0.5,main="NJ")
abline(h=0,lty=3)
#fit the data to the tree using pml
fit <- pml(tr.nj,as.phyDat(x))
#optimise this, then plot and set a random seed for the boostrap
fit=optim.pml(fit,T)
plot(fit)
set.seed(8)

#bootstrap the results 100 times to improve confidence in the position of each node in the tree
bs<-bootstrap.pml(fit,bs=100,optNni=T)
#plot the boostrapped tree in prefferred format, in this case the p format
treeBS<-plotBS(fit$tree, type"p", bs)
#perform a model test to find the most appropriate mutation model
mt<-modelTest(as.phyDat(x),G=F,I=F)

