#############################################################
#
#  Honours project  
# 
#  Code to align the barley sequences for CCR2 and KNAT
#
#############################################################

########## Clean-up the memory and start a new session #################

rm(list=ls())
dev.off()

########## Libraries required ##########################################

library('DECIPHER')
library('seqinr')
library('msa')

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of 
#the code
sessionInfo()

#set users working directory  
setwd("~/Documents/4th Year/Honours Project/Data/barley genes/seq FASTAS")

########## Preparing the files for alignment ###########################

setwd("~/Documents/4th Year/Honours Project/Data/barley genes/CMCF002-/CMCF002")

fileNames <- c('L2 KNAT 4 F', 'R2 KNAT 4 F', 'A3 KNAT 4 F', 'B3 KNAT 4 F', 'F3 KNAT 4 F', 'B4 KNAT 4 F', 
               'S4 KNAT 4 F', 'SW4 KNAT 4 F', 'L2 KNAT 4 R', 'R2 KNAT 4 R', 'A3 KNAT 4 R', 'B3 KNAT 4 R', 
               'F3 KNAT 4 R', 'B4 KNAT 4 R', 'S4 KNAT 4 R', 'SW4 KNAT 4 R', 'A5 KNAT 4 F', 'S5 KNAT 4 F', 
               'T5 KNAT 4 F', 'A6 KNAT 4 F', 'I6 KNAT 4 F', 'A1 KNAT 5 F', 'A2 KNAT 5 F', 'L2 KNAT 5 F', 
               'A5 KNAT 4 R', 'S5 KNAT 4 R', 'T5 KNAT 4 R', 'A6 KNAT 4 R', 'I6 KNAT 4 R', 'A1 KNAT 5 R',
               'A2 KNAT 5 R', 'L2 KNAT 5 R', 'R2 KNAT 5 F', 'A3 KNAT 5 F', 'B3 KNAT 5 F', 'F3 KNAT 5 F',
               'B4 KNAT 5 F', 'S4 KNAT 5 F', 'SW4 KNAT 5 F', 'A5 KNAT 5 F', 'R2 KNAT 5 R', 'A3 KNAT 5 R',
               'B3 KNAT 5 R', 'F3 KNAT 5 R', 'B4 KNAT 5 R', 'S4 KNAT 5 R', 'SW4 KNAT 5 R', 'A5 KNAT 5 R')

Files <- list.files(pattern = '.ab1')

for (counter in 1:96) {
  file.rename(from = Files[counter], to = paste(fileNames[counter],'.ab1', sep = ''))
}

########## aligning the sequences ######################################

#Creating lists of the forward and reveres sequences
forwardSeqFiles <- list.files(pattern = 'F.fasta')
reverseSeqFiles <- list.files(pattern = 'R.fasta')

# load the sequences from the file
forwardSeqs <- lapply(forwardSeqFiles, readDNAStringSet)
write.fasta(forwardSeqs,  names = lapply(forwardSeqs, names), file.out ='forwardSeqs.fasta')
forwardSeqs <- readDNAStringSet('forwardSeqs.fasta')

reverseSeqs <- lapply(reverseSeqFiles, readDNAStringSet)
write.fasta(reverseSeqs, names = lapply(reverseSeqs, names), file.out ='reverseSeqs.fasta')
reverseSeqs <- readDNAStringSet('reverseSeqs.fasta')
  
#importing the referece files
KNAT <- readDNAStringSet("~/Documents/4th Year/Honours Project/Data/barley genes/KNATGeneSequence.fasta")
CCR2 <- readDNAStringSet("~/Documents/4th Year/Honours Project/Data/barley genes/CCR2GeneSequence.fasta")

#nucleotide sequences need to be in the same orientation. If they are not, then they can be reoriented 
orintatedSeqs <- reverseComplement(reverseSeqs)

#subsetting for the alignment
LCCR2 <- c(forwardSeqs[grep('L C', names(forwardSeqs))], orintatedSeqs[grep('L C', names(orintatedSeqs))])
LKNAT <- c(forwardSeqs[grep('L K', names(forwardSeqs))], orintatedSeqs[grep('L K', names(orintatedSeqs))])
TCCR2 <- c(forwardSeqs[grep('T C', names(forwardSeqs))], orintatedSeqs[grep('T C', names(orintatedSeqs))])
TKNAT <- c(forwardSeqs[grep('T K', names(forwardSeqs))], orintatedSeqs[grep('T K', names(orintatedSeqs))])
A1CCR2 <- c(forwardSeqs[grep('A1 C', names(forwardSeqs))], orintatedSeqs[grep('A1 C', names(orintatedSeqs))])
A2CCR2 <- c(forwardSeqs[grep('A2 C', names(forwardSeqs))], orintatedSeqs[grep('A2 C', names(orintatedSeqs))])
A3CCR2 <- c(forwardSeqs[grep('A3 C', names(forwardSeqs))], orintatedSeqs[grep('A3 C', names(orintatedSeqs))])
A5CCR2 <- c(forwardSeqs[grep('A5 C', names(forwardSeqs))], orintatedSeqs[grep('A5 C', names(orintatedSeqs))])
A6CCR2 <- c(forwardSeqs[grep('A6 C', names(forwardSeqs))], orintatedSeqs[grep('A6 C', names(orintatedSeqs))])
B3CCR2 <- c(forwardSeqs[grep('B3 C', names(forwardSeqs))], orintatedSeqs[grep('B3 C', names(orintatedSeqs))])
B4CCR2 <- c(forwardSeqs[grep('B4 C', names(forwardSeqs))], orintatedSeqs[grep('B4 C', names(orintatedSeqs))])
F3CCR2 <- c(forwardSeqs[grep('F3 C', names(forwardSeqs))], orintatedSeqs[grep('F3 C', names(orintatedSeqs))])
I6CCR2 <- c(forwardSeqs[grep('I6 C', names(forwardSeqs))], orintatedSeqs[grep('I6 C', names(orintatedSeqs))])
L2CCR2 <- c(forwardSeqs[grep('L2 C', names(forwardSeqs))], orintatedSeqs[grep('L2 C', names(orintatedSeqs))])
R2CCR2 <- c(forwardSeqs[grep('R2 C', names(forwardSeqs))], orintatedSeqs[grep('R2 C', names(orintatedSeqs))])
S4CCR2 <- c(forwardSeqs[grep('S4 C', names(forwardSeqs))], orintatedSeqs[grep('S4 C', names(orintatedSeqs))])
S5CCR2 <- c(forwardSeqs[grep('S5 C', names(forwardSeqs))], orintatedSeqs[grep('S5 C', names(orintatedSeqs))])
SW4CCR2 <- c(forwardSeqs[grep('SW4 C', names(forwardSeqs))], orintatedSeqs[grep('SW4 C', names(orintatedSeqs))])
T5CCR2 <- c(forwardSeqs[grep('T5 C', names(forwardSeqs))], orintatedSeqs[grep('T5 C', names(orintatedSeqs))])
A1KNAT <- c(forwardSeqs[grep('A1 K', names(forwardSeqs))], orintatedSeqs[grep('A1 K', names(orintatedSeqs))])
A2KNAT <- c(forwardSeqs[grep('A2 K', names(forwardSeqs))], orintatedSeqs[grep('A2 K', names(orintatedSeqs))])
A3KNAT <- c(forwardSeqs[grep('A3 K', names(forwardSeqs))], orintatedSeqs[grep('A3 K', names(orintatedSeqs))])
A5KNAT <- c(forwardSeqs[grep('A5 K', names(forwardSeqs))], orintatedSeqs[grep('A5 K', names(orintatedSeqs))])
A6KNAT <- c(forwardSeqs[grep('A6 K', names(forwardSeqs))], orintatedSeqs[grep('A6 K', names(orintatedSeqs))])
B3KNAT <- c(forwardSeqs[grep('B3 K', names(forwardSeqs))], orintatedSeqs[grep('B3 K', names(orintatedSeqs))])
B4KNAT <- c(forwardSeqs[grep('B4 K', names(forwardSeqs))], orintatedSeqs[grep('B4 K', names(orintatedSeqs))])
F3KNAT <- c(forwardSeqs[grep('F3 K', names(forwardSeqs))], orintatedSeqs[grep('F3 K', names(orintatedSeqs))])
I6KNAT <- c(forwardSeqs[grep('I6 K', names(forwardSeqs))], orintatedSeqs[grep('I6 K', names(orintatedSeqs))])
L2KNAT <- c(forwardSeqs[grep('L2 K', names(forwardSeqs))], orintatedSeqs[grep('L2 K', names(orintatedSeqs))])
R2KNAT <- c(forwardSeqs[grep('R2 K', names(forwardSeqs))], orintatedSeqs[grep('R2 K', names(orintatedSeqs))])
S4KNAT <- c(forwardSeqs[grep('S4 K', names(forwardSeqs))], orintatedSeqs[grep('S4 K', names(orintatedSeqs))])
S5KNAT <- c(forwardSeqs[grep('S5 K', names(forwardSeqs))], orintatedSeqs[grep('S5 K', names(orintatedSeqs))])
SW4KNAT <- c(forwardSeqs[grep('SW4 K', names(forwardSeqs))], orintatedSeqs[grep('SW4 K', names(orintatedSeqs))])
T5KNAT <- c(forwardSeqs[grep('T5 K', names(forwardSeqs))], orintatedSeqs[grep('T5 K', names(orintatedSeqs))])


#perform the alignment 
alignedHaplotype1CCR2 <- DNAStringSet(msaClustalOmega( c(CCR2, LCCR2, TCCR2, A1CCR2) ,order = "input"))
alignedHaplotype1KNAT <- DNAStringSet(msaClustalOmega( c(KNAT, LKNAT, TKNAT, A1KNAT) ,order = "input"))
alignedHaplotype2CCR2 <- DNAStringSet(msaClustalOmega(c(CCR2, A2CCR2, L2CCR2, R2CCR2), order = "input")) 
alignedHaplotype2KNAT <- DNAStringSet(msaClustalOmega(c(KNAT, A2KNAT, L2KNAT, R2KNAT), order = "input"))
alignedHaplotype3CCR2 <- DNAStringSet(msaClustalOmega(c(CCR2, A3CCR2, B3CCR2, F3CCR2), order = "input"))
alignedHaplotype3KNAT <- DNAStringSet(msaClustalOmega(c(KNAT, A3KNAT, B3KNAT, F3KNAT), order = "input"))
alignedHaplotype4CCR2 <- DNAStringSet(msaClustalOmega(c(CCR2, B4CCR2, S4CCR2, SW4CCR2), order = "input"))
alignedHaplotype4KNAT <- DNAStringSet(msaClustalOmega(c(KNAT, B4KNAT, S4KNAT, SW4KNAT), order = "input"))
alignedHaplotype5CCR2 <- DNAStringSet(msaClustalOmega(c(CCR2, A5CCR2, S5CCR2, T5CCR2), order = "input"))
alignedHaplotype5KNAT <- DNAStringSet(msaClustalOmega(c(KNAT, A5KNAT, S5KNAT, T5KNAT), order = "input"))
alignedHaplotype6CCR2 <- DNAStringSet(msaClustalOmega(c(CCR2, A6CCR2, I6CCR2), order = "input"))
alignedHaplotype6KNAT <- DNAStringSet(msaClustalOmega(c(KNAT, A6KNAT, I6KNAT), order = "input"))


#visualising the alignments 
BrowseSeqs(alignedHaplotype1CCR2, highlight = 1)

