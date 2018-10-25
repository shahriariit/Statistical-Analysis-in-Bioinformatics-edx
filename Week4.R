library(Biostrings)
library(seqinr)

AB003468 <- readDNAStringSet("AB003468.fasta")
AB003468 <- as.character(AB003468)

matchPattern("ATG", AB003468)

sequence <- AB003468

start_codon <- "ATG"
stop_codons <- c("TGA","TAA","TAG")

start_pos <- c()
revstart_pos <- c()
stop_pos <- c()
revstop_pos <- c()

matches <- matchPattern(start_codon,sequence)
start_pos <- c(start_pos, start(matches))

revmatches <- matchPattern(reverseComplement(DNAString(start_codon)), sequence)
revstart_pos <- c(revstart_pos, start(revmatches))

start_pos <- sort(start_pos)
revstart_pos <- sort(revstart_pos, decreasing = TRUE)

for (codon in stop_codons) {
  matches <- matchPattern(codon, sequence)
  stop_pos <- c(stop_pos, start(matches))
  revmatches <- matchPattern(reverseComplement(DNAString(codon)),sequence)
  revstop_pos <- c(revstop_pos, start(revmatches))
}

stop_pos <- sort(stop_pos)
revstop_pos <- sort(revstop_pos, decreasing = TRUE)

k <- 150
stop_pointers <- c(0,0,0)
count <- 0

for (current_start in start_pos) {
  frame <- (current_start%%3) + 1
  stop_pointer <- stop_pointers[frame]
  if (stop_pointer <= length(stop_pos) && (stop_pointer == 0 || stop_pos[stop_pointer] < current_start)) {
     stop_pointer <- stop_pointer + 1
     
     while((stop_pointer <= length(stop_pos)) && (stop_pos[stop_pointer] < current_start)||(((stop_pos[stop_pointer]%%3)+1)!=frame)){
       stop_pointer <- stop_pointer + 1
     }
     stop_pointers[frame]=stop_pointer
     
     if(stop_pointer <= length(stop_pos)){
       if ((stop_pos[stop_pointer] + 2 - current_start + 1) > k ) {
         count = count +1
         print(count)
         print("Frame")
         print(frame)
         print("Start:")
         print(current_start)
         print("Stop:")
         print(stop_pos[stop_pointer])
         print("Length:")
         lengths<-c(lengths,(stop_pos[stop_pointer] + 2 - current_start + 1))
         print(stop_pos[stop_pointer] + 2 - current_start + 1)
         print("Sequence:")
         print(subseq(sequence,current_start,stop_pos[stop_pointer]+2))
       }
     }
     
  }
}

#Reverse 3 frames

k <- 150
revstop_pointers <- c(0,0,0)
count <- 0

for (current_revstart in revstart_pos) {
  current_revstart <- current_revstart + 2
  frame <- (current_revstart%%3) + 1
  revstop_pointer <- revstop_pointers[frame]
  if (revstop_pointer <= length(revstop_pos) && (revstop_pointer == 0 || revstop_pos[revstop_pointer])) {
    revstop_pointer <- revstop_pointer + 1
    
    while((revstop_pointer <= length(revstop_pos)) && ((revstop_pos[revstop_pointer]+2 < current_revstart)||((((revstop_pos[revstop_pointer]+2)%%3)+1)!=frame))){
      revstop_pointer <- revstop_pointer + 1
    }
    revstop_pointers[frame]=revstop_pointer
    
    if(revstop_pointer <= length(revstop_pos)){
      if ((current_revstart-revstop_pos[revstop_pointer]) + 1 > k ) {
        count = count +1
        print(count)
        print("Frame")
        print(frame)
        print("Start:")
        print(current_revstart)
        print("Stop:")
        print(revstop_pos[revstop_pointer])
        print("Length:")
        lengths<-c(lengths,(current_revstart-revstop_pos[revstop_pointer]))
        print(current_revstart-revstop_pos[revstop_pointer])
        print("Sequence:")
        print(subseq(sequence,revstop_pos[revstop_pointer],current_revstart))
      }
    }
    
  }
}

lengths <- c(lengths, (stop_pos[stop_pointer]+2-current_start+1))
print("Length:"):
lengths <- c(lengths, (current_revstart - revstop_pos[revstop_pointer]))
lengths <- vector(mode="numeric")

lengths <- sort(lengths)
barplot(lengths)
plot(density(lengths))

bins <- seq(0,1000,50)
hist(lengths, breaks=bins, col="red", xlim=c(0,1000))


