#############################################################
#
#  Honours project  
# 
#  Code to create SNP files for VEP analysis and
#  remove the unwanted data form VEP output
#
#############################################################

########## Clean-up the memory and start a new session ######
rm(list=ls())
dev.off()

########## Libraries required ###############################

#required packages 

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of 
#the code
sessionInfo()

#set users working directory  
setwd("~/Documents/4th Year/Honours Project/Data/rice genes/snp data")

########## Data cleaning for VEP ############################

#Search a directory and create a list of all the files present in it 
SNPFileNames <- list.files("~/Documents/4th Year/Honours Project/Data/rice genes/snp data")

for(File in SNPFileNames) {
 
  #import the csv file
  SNPdata <- read.table(File, header = T, sep = '\t')
  
  #The first step is to prune the data and sub set just for the data we need
  trimmedSNPdata <- SNPdata[, c( 1, 2, 3, 4, 6)]
  
  #Next step is to restructure the data to how we need it 
  trimmedSNPdata$allele <- paste(trimmedSNPdata$Ref_nt, '/', trimmedSNPdata$Alt_nt, sep = "")
  trimmedSNPdata$start <- trimmedSNPdata$Positon
  trimmedSNPdata$end <- trimmedSNPdata$Positon
  
  trimmedSNPDatafile2 <- trimmedSNPdata[, c( 1, 7, 8, 6)]
  strand <- readline(prompt=paste(File, "Strand direction (+/-): "))
  trimmedSNPDatafile2$direction <- strand
  names(trimmedSNPDatafile2) <- NULL
  
  #after this we need to export the file
  writePath <- paste("~/Documents/4th Year/Honours\ Project/Data/rice\ genes/Vep input/", File,
                     '.bed', sep = "")
  write.table(trimmedSNPDatafile2, writePath, sep="\t", row.names = F, quote = F)
  
}


######### removing the unwanted data form VEP output ########

#set users working directory 
setwd("~/Documents/4th Year/Honours Project/Data/rice genes/vep output")

#Search a directory and create a list of all the files present in it 
vepInputFileNames <- list.files("~/Documents/4th Year/Honours Project/Data/rice genes/Vep input")
vepOutputFileNames <- list.files("~/Documents/4th Year/Honours Project/Data/rice genes/vep output")
counter <- 1

for(File in vepOutputFileNames) {
  
  #creating the file path for the csv files
  vepInputFileName <- paste("~/Documents/4th Year/Honours Project/Data/rice genes/Vep input/", 
                            vepInputFileNames[counter], sep = "")
  
  #import the csv file
  Vepinputdata <- read.table(vepInputFileName, header = F, sep = '\t')
  Vepdata <- read.table(File, header = F, sep = '\t')

  #the next step is to create a series of if statements to deal with the possible strands of the various 
  #files 
  
  if (levels(Vepinputdata$V5) == '+'){
    trimmedVepOutput <- rownames(Vepdata)[which(Vepdata$V16 == "1")]
    snpsOfIntrests <- Vepdata[trimmedVepOutput, ]
    
    #assigning the new table column names
    colnames(snpsOfIntrests) <- c('#Uploaded_variation',	'Location',	'Allele',	'Gene',	'Feature',	
                                  'Feature_type',	'Consequence',	'cDNA_position',	'CDS_position',	
                                  'Protein_position',	'Amino_acids',	'Codons',	'Existing_variation',
                                  'IMPACT',	'DISTANCE',	'STRAND',	'FLAGS',	'SIFT')
    
    #after this we need to export the file
    writePath <- paste("~/Documents/4th Year/Honours\ Project/Data/rice\ genes/SNPs\ of\ intrest/", File, 
                       sep = "")
    write.table(snpsOfIntrests, writePath, sep="\t", row.names = F, quote = F)
  }
  
  if (levels(Vepinputdata$V5) == '-'){
    trimmedVepOutput <- rownames(Vepdata)[which(Vepdata$V16 == "-1")]
    snpsOfIntrests <- Vepdata[trimmedVepOutput, ]
    
    #assigning the new table column names
    colnames(snpsOfIntrests) <- c('#Uploaded_variation',	'Location',	'Allele',	'Gene',	'Feature',	
                                  'Feature_type',	'Consequence',	'cDNA_position',	'CDS_position',	
                                  'Protein_position',	'Amino_acids',	'Codons',	'Existing_variation',
                                  'IMPACT',	'DISTANCE',	'STRAND',	'FLAGS',	'SIFT')
    
    #after this we need to export the file
    writePath <- paste("~/Documents/4th Year/Honours\ Project/Data/rice\ genes/SNPs\ of\ intrest/", File, 
                       sep = "")
    write.table(snpsOfIntrests, writePath, sep="\t", row.names = F, quote = F)
  }

  counter <- counter + 1
  
}





