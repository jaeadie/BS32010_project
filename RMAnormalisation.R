source("http://bioconductor.org/biocLite.R")
library(affy)
library(affyPLM)
library(annotate)
library(chicken.db)
library(limma)

data=ReadAffy()
ed=exprs(data)
samp=sampleNames(data)
probes=featureNames(data)
#graphs for non-normalised data
png("ROS1.png", width = 1920, height = 1080, units = "px")
plot(ed[,"ROS1-_9.CEL"], ed[,"ROS1+_10.CEL"], xlab='ROS1-_9', ylab='ROS1+_10', main='ROS1')
dev.off()
png("ROS2.png", width = 1920, height = 1080, units = "px")
plot(ed[,"ROS2-_11.CEL"], ed[,"ROS2+_12.CEL"], xlab='ROS2-_11', ylab='ROS2+_12', main='ROS2')
dev.off()
png("ROS3.png", width = 1920, height = 1080, units = "px")
plot(ed[,"ROS3-_13.CEL"], ed[,"ROS3+_14.CEL"], xlab='ROS3-_13', ylab='ROS3+_14', main='ROS3')
dev.off()
png("ROS4.png", width = 1920, height = 1080, units = "px")
plot(ed[,"ROS4-_15.CEL"], ed[,"ROS4+_16.CEL"], xlab='ROS4-_15', ylab='ROS4+_16', main='ROS4')
dev.off()

#normalisation
nvals=rma(data)
ned=exprs(nvals)
nsamp=sampleNames(nvals)
nprobes=featureNames(nvals)
png("ROS1Normalised.png", width = 1920, height = 1080, units = "px")
plot(ned[,"ROS1-_9.CEL"], ned[,"ROS1+_10.CEL"], xlab='ROS1-_9', ylab='ROS1+_10', main='ROS1 Normalised')
dev.off()
png("ROS2Normalised.png", width = 1920, height = 1080, units = "px")
plot(ned[,"ROS2-_11.CEL"], ned[,"ROS2+_12.CEL"], xlab='ROS2-_11', ylab='ROS2+_12', main='ROS2 Normalised')
dev.off()
png("ROS3Normalised.png", width = 1920, height = 1080, units = "px")
plot(ned[,"ROS3-_13.CEL"], ned[,"ROS3+_14.CEL"], xlab='ROS3-_13', ylab='ROS3+_14', main='ROS3 Normalised')
dev.off()
png("ROS4Normalised.png", width = 1920, height = 1080, units = "px")
plot(ned[,"ROS4-_15.CEL"], ned[,"ROS4+_16.CEL"], xlab='ROS4-_15', ylab='ROS4+_16', main='ROS4 Normalised')
dev.off()