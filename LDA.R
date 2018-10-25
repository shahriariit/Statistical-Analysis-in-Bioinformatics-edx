source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(biomaRt)
peptides.txt <- read.table("peptidefrags.txt", header=FALSE)
peptides <-as.vector(peptides.txt$V1)
hist(peptides,breaks=400) 

library(gplots)

mascot.txt <- read.table("mascot.txt", header=FALSE)
mascot <-as.vector(mascot.txt$V1)
xtandem.txt <- read.table("xtandem.txt", header=FALSE)
xtandem <-as.vector(xtandem.txt$V1)
protpro.txt <- read.table("protpro.txt", header=FALSE)
protpro <-as.vector(protpro.txt$V1)

combinedMSdata <- list(Mascot=mascot, XTandem=xtandem, ProtPro=protpro)
venn(combinedMSdata)

library(timeSeries)
library(MASS)
library(rgl)
library(ggplot2)

Dataset <- read.csv("ms.csv", header=TRUE, na.strings="NA", dec=".", strip.white=TRUE)
RawData <- Dataset[,2:14]
filledcols = colSds(RawData) != 0.0
RawData <- RawData[,filledcols]

test1.lda <- lda(Dataset$X1 ~ . , data=Dataset)
test1.lda.values <- predict(test1.lda,Dataset)

x <- test1.lda.values$x[,1]
y <- test1.lda.values$x[,2]

class <- Dataset$X1

plotdata <- data.frame(class, x, y)
centroids <- aggregate(cbind(x,y)~class,plotdata,mean)

CentroidDistances <- dist(centroids, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
plot1 <- ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)
plot1
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7)
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black")
ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") + geom_text(aes(label=Dataset$X1),hjust=0, vjust=0, colour="black")
plot1 <- ggplot(plotdata,aes(x,y,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") + geom_text(aes(label=Dataset$X1),hjust=0, vjust=0, colour="black")
ggsave(filename="plot1.pdf")
write.csv(as.matrix(CentroidDistances, file = "centroiddistances.csv"))