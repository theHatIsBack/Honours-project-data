#############################################################
#
#  Honours project  
# 
#  Code to create the gene model for rice figure
#  
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
setwd("~/Documents/4th Year/Honours Project/Data/rice genes")

########## creating the gene ids ############################

files <- list.files('~/Documents/4th Year/Honours Project/Data/rice genes/snp data')
genes <- c(unlist(strsplit(files, '.', fixed = T)))
genes <- genes[grep(genes, pattern = 'Os')]
geneIDs <- paste('LOC', genes, sep = '_')

#adding it them to the table
geneModel <- as.data.frame(geneIDs)


## creating the strand direction, start/end and width columns ####

counter <- 1
widths <- c()
ends <- c()
starts <- c()
directions <- c()

for(gene in geneIDs) {
  
  #strand direction
  strand <- readline(prompt=paste(gene, "Strand direction (+/-): "))
  directions[counter] <- strand
  
  #start position
  start <- readline(prompt=paste(gene, "Start Position (int): "))
  starts[counter] <- start
  
  #end position
  end <- readline(prompt=paste(gene, "End Position (int): "))
  ends[counter] <- end
  
  #widths
  widths[counter] <- as.integer(end) - as.integer(start)
  
  counter <- counter + 1
}

#putting it all together

geneModel$widths <- widths
geneModel$ends <- ends
geneModel$starts <- starts
geneModel$directions <- directions
geneModel$chr <- 'chr2'

#reordering the columns
geneModel <- geneModel[, c("chr", "starts", "ends", "widths", "directions", "geneIDs")]
colnames(geneModel) <- c("chromosome", "start", "end", "width", "strand", "gene")

#after this we need to export the model
writePath <- paste("~/Documents/4th Year/Honours\ Project/Data/rice\ genes/", 'rGeneModel',
                   '.txt', sep = "")
write.table(geneModel, writePath, sep="\t", row.names = F, quote = F)