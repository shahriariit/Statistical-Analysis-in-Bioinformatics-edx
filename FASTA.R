library(seqinr)
cox1 <- read.fasta(file = "cox1.fasta",seqtype = "AA")
length(cox1) 
cox1[1]

library(ape)
AB003468 <- read.GenBank("AB003468", as.character = "TRUE")

write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

install.packages("rentrez")
library(rentrez)

CloningVector <- AB003468[[1]]
count <- count(CloningVector,1) 
GC <- GC(CloningVector)
GCwindow <- seq(1, length(CloningVector)-200, by = 200)
n <- length(GCwindow)
Chunks <- numeric(n)

for(i in 1:n){
  Chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
  ChunkGC <- GC(Chunk)
  print(ChunkGC)
  Chunks[i]<-ChunkGC
}

plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content")


slidingwindowGCplot <- function(windowsize,inputseq){
  GCwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(GCwindow)
  Chunks <- numeric(n)
  for(i in 1:n){
    Chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
    ChunkGC <- GC(Chunk)
    print(ChunkGC)
    Chunks[i]<-ChunkGC
  }
  plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize ", windowsize))
}


# Protein seq statistics

library(Peptides)
aaComp(cox1[1])

charge(cox1)
charge(seq="FLPVLAG", pH=7, pKscale="EMBOSS")

hydrophobicity(cox1[1])

slidingwindowGCplot <- function(windowsize,inputseq){
  GCwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(GCwindow)
  Chunks <- numeric(n)
  for(i in 1:n){
    Chunk <- CloningVector[GCwindow[i]:(GCwindow[i]+199)]
    ChunkGC <- hydrophobicity(Chunk)
    print(ChunkGC)
    Chunks[i]<-ChunkGC
  }
  plot(GCwindow,Chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize ", windowsize))
}

