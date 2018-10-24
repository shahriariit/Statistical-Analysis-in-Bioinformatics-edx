source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

biocLite("msa")


library(Biostrings)
library(seqinr)
library(msa)

prokaryotes <- read.fasta(file = "prok.fasta", seqtype = "DNA")

seq1=as.character(prokaryotes[[1]])
seq1=paste(seq1,collapse ="")
seq2=as.character(prokaryotes[[2]])
seq2=paste(seq1,collapse ="")

pairalign <- pairwiseAlignment(pattern = seq2, subject = seq1)

pairalignString = BStringSet( c( toString( subject(pairalign) ), toString(pattern(pairalign))))

writeXStringSet(pairalignString, "aligned.txt", format="FASTA")

coxgenes <- read.fasta(file = "cox1multi.fasta", seqtype="AA")

cox1 <- as.character(coxgenes[[1]])

cox2 <- as.character(coxgenes[[2]])

coxAA <- readAAStringSet("cox1multi.fasta")
prokDNA <- readDNAStringSet("prok.fasta")

coxAligned <- msa(coxAA)
prokAligned <- msa(prokDNA)

print(coxAligned, show="complete")
print(prokAligned, show="complete")

msa(prokDNA, "ClustalW")
msa(prokDNA, "ClustalOmega")
msa(prokDNA, "Muscle")

msa(coxAA,cluster="upgma")

prokAlignStr = as(prokAligned, "DNAStringSet")
writeXStringSet(prokAlignStr, file="prokAligned.fasta")

coxAlignStr = as(coxAligned, "AAStringSet")
writeXStringSet(coxAlignStr, file="coxAligned.fasta")

write.phylip(coxAligned, "coxAligned.phylip")

prokAligned2 <- msaConvert(prokAligned, type="seqinr::alignment")
prokdist <- dist.alignment(prokAligned2, "identity")

library(ape)
prokTree <- nj(prokdist)

plot(prokTree)

library(phangorn)
prokAligned3 <- msaConvert(prokAligned, type="phangorn::phyDat")
ParsTree <- pratchet(prokAligned3)
plot(ParsTree)

fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
plot(fitJC)

bootstrapped <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=FALSE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bootstrapped, p = 50, type="p")

