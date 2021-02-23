#############################################################
#
#  Honours project  
# 
#  Code to create the rice figure
#  
#
#############################################################

####### Clean-up the memory and start a new session #########
rm(list=ls())
dev.off()

########## Libraries required ###############################

#required packages
library(readxl)
library(Gviz)
library(magrittr)
library(dplyr)
library(RColorBrewer)
library(gt)


#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of 
#the code
sessionInfo()

#set users working directory  
setwd("~/Documents/4th Year/Honours Project/Data/rice genes")

########## extracting the snps of intrest ###################

riceGenes <- read_excel("rice genes.xlsx", sheet = "SNPs", col_names = T, skip = 2)
riceGenesOfInterest <- rbind(riceGenes[(rownames(riceGenes)[grep(riceGenes$SIFT, pattern = 'del')]), ], 
                             riceGenes[(rownames(riceGenes)[grep(riceGenes$IMPACT, pattern = 'H')]), ], 
                             deparse.level = 0)
riceGenesOfInterest <- riceGenesOfInterest[!duplicated(riceGenesOfInterest), ]


############### creating the ideogram ######################

#importing the data  for the ideogram
IdeogramData <- read.table(file = "chrosome 2", header = T, sep = "\t", stringsAsFactors = FALSE)

#generating the ideogram
itrack <- IdeogramTrack(chromosome = "chr2", genome = "S", name=NULL, bands=IdeogramData, bevel = 0.125)

axis <- GenomeAxisTrack()


############### ploting the snps ######################

#the first step is to prepare the data for the annotation track
SNPPosList <- c(unlist(strsplit(riceGenesOfInterest$`#Uploaded_variation`, "_", fixed = T)))
SNPPosList2 <- SNPPosList[grep(SNPPosList, pattern = '4')] 
SNPPosList <- paste(SNPPosList2, '-', SNPPosList2, sep = '')

#creating a granges class of the data to feed into the annotation track function
SNPClass <- GRanges(seqnames = 'chr2', ranges = SNPPosList , strand = '+')

#creating the snps annotation track
SNPs <- AnnotationTrack(SNPClass, name = "SNPs", stacking="dense", stackHeight=0.25)


############### creating the gene model ######################

geneModel <- read.table('rGeneModel.txt', header = T, sep = '\t')

IDs <- c()
counter <- 1
for (gene in geneModel$gene) {
  for (pos in SNPPosList2) {
    if( pos >= geneModel$start[counter] & pos <= geneModel$end[counter]){
      IDs[counter] <- gene
    } 
  }
  counter <- counter + 1 
}
IDs <- IDs[grep(IDs, pattern = 'LOC')]

geneModel2 <- as.data.frame(c())
for(x in IDs){ 
  geneModel3 <- geneModel[(rownames(geneModel)[grep(geneModel$gene, pattern = x )]), ]
  geneModel2 <- rbind(geneModel2, geneModel3)
}


grtrack <- GeneRegionTrack(geneModel2, genome = 's', chromosome = 2, name = "Rice Gene Model", 
                           stackHeight=0.25)


############### plotting the expression data ######################

exspresionData <- read.table('~/Documents/4th Year/Honours Project/Data/rice genes/gene exspresion data/final_rice_v7_expression_matrix_48columns.txt',
                             header = T, sep = '\t', stringsAsFactors = F)
gOIExepresionData <- as.data.frame(c())
 
for(x in IDs){ 
  GED <- exspresionData[(rownames(exspresionData)[grep(exspresionData$Locus_id, pattern = x )]), ]
  gOIExepresionData <- rbind(gOIExepresionData, GED)
}

gOIExepresionData <- gOIExepresionData[, c("Locus_id", "SRR074151", "SRR074171", "SRR074146")]

gOIExepresionData2 <- gOIExepresionData %>% 
  mutate(SRR074151 = replace(SRR074151, SRR074151 == "yes", as.integer('1'))) %>%
  mutate(SRR074151 = replace(SRR074151, SRR074151 == "no", as.integer('0'))) %>%
  mutate(SRR074171 = replace(SRR074171, SRR074171 == "yes", as.integer('1'))) %>%
  mutate(SRR074171 = replace(SRR074171, SRR074171 == "no", as.integer('0'))) %>%
  mutate(SRR074146 = replace(SRR074146, SRR074146 == "yes", as.integer('1'))) %>%
  mutate(SRR074146 = replace(SRR074146, SRR074146 == "no", as.integer('0'))) 

colnames(gOIExepresionData2) <- c("Locus_id", "Stems", "Roots", "Leaves")

exspresionModel <- GRanges(seqnames = "chr2", strand = '*', 
                           ranges = IRanges(start = geneModel2$start, end = geneModel2$end), 
                           gOIExepresionData2 )

dTrack <- DataTrack(exspresionModel, name = "Exspresion Data")


############### plotting the different tracks ######################

plotTracks(c(itrack, axis, SNPs, grtrack, dTrack), from=4407800, to=4776300, transcriptAnnotation = "gene",
          type = 'heatmap', extend.left = 0.15, extend.right = 0.05,showSampleNames = TRUE, 
          cex.sampleNames = 0.6, separator = 4 , gradient = c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", 
                                                              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"))

#groupAnnotation = "group", just.group = "above"

############### creating a table of Rice genes and SNPs ######################

SNPPosList <- c(unlist(strsplit(riceGenesOfInterest$`#Uploaded_variation`, "_", fixed = T)))
SNPList <- SNPPosList[grep(SNPPosList, pattern = '/')]

geneIDs <- as.data.frame(c())
SNP <- c()
posistion <- c()

for (pos in 1:142) {
  for (count in 1:41) {
    if( SNPPosList2[pos] >= geneModel$start[count] & SNPPosList2[pos] <= geneModel$end[count]){
      x <- as.character(geneModel$gene[count])
      gene <- geneModel[(rownames(geneModel)[grep(geneModel$gene, pattern = x)]), ]
      SNP[pos] <- SNPList[pos]
      posistion[pos] <- SNPPosList2[pos]
      geneIDs <- rbind( geneIDs, gene)
      break
    } 
  }
}

geneIDs$snps <- na.omit(SNP)
geneIDs$position <- na.omit(posistion)
gOIExepresionData <- as.data.frame(c())

for(x in IDs){ 
  GED <- exspresionData[(rownames(exspresionData)[grep(exspresionData$Locus_id, pattern = x )]), ]
  gOIExepresionData <- rbind(gOIExepresionData, GED)
}

gOIExepresionData <- gOIExepresionData[(rownames(gOIExepresionData)[grep(gOIExepresionData$"SRR074151", 
                                                                         pattern = 'y')]), ] 
genesOfIntrest <- as.data.frame(c())

for(x in gOIExepresionData$Locus_id){ 
  genes <- geneIDs[(rownames(geneIDs)[grep(geneIDs$gene, pattern = x )]), ]
  genesOfIntrest <- rbind(genesOfIntrest, genes)
}

genesOfIntrest <- genesOfIntrest[,  c("chromosome", "gene", "position", "snps")]


genesOfIntrest <- rbind(genesOfIntrest[(rownames(genesOfIntrest)[grep(genesOfIntrest$gene, pattern = 'LOC_Os02g08420')]), ],
                         genesOfIntrest[(rownames(genesOfIntrest)[grep(genesOfIntrest$gene, pattern = 'LOC_Os02g09260')]), ])

functions <- c('cinnamoyl CoA reductase', 'cinnamoyl CoA reductase', 'cytochrome P450 71A1', 
               'cytochrome P450 71A1', 'cytochrome P450 71A1', 'cytochrome P450 71A1','cytochrome P450 71A1',
               'cytochrome P450 71A1', 'cytochrome P450 71A1', 'cytochrome P450 71A1', 'cytochrome P450 71A1',
               'cytochrome P450 71A1','cytochrome P450 71A1','cytochrome P450 71A1','cytochrome P450 71A1',
               'cytochrome P450 71A1','cytochrome P450 71A1','cytochrome P450 71A1')

genesOfIntrest$'Gene\ Annotation' <- functions

colnames(genesOfIntrest) <- c("Chromosome", "Gene", "Position", "Alleles", "Gene Annotation")

genesOfIntrest %>% gt()


