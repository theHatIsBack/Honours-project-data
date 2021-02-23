#############################################################
#
#  Honours project  
# 
#  code to analysis barley haplotype data
#
#############################################################

########## Clean-up the memory and start a new session #################

rm(list=ls())
dev.off()

########## Libraries required ##########################################

library(PMCMR)
library(magrittr)
library(dplyr)

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of 
#the code
sessionInfo()

#set users working directory  
setwd("~/Documents/4th Year/Honours Project/Data/barley genes")

##################### Plotting the raw data ############################

#importing the data
haplotypes <- read.table('haplotypeExspresion', header = T, sep = ',')

#boxplot of raw data
boxplot(haplotypes$X1~haplotypes$haplotype., xlab = "", ylab = "FT-IR Lignin Content (%)")

#Shapiro-wilk test for normal distribution:
shapiro.test(haplotypes$X1)
#As the P value is less than 0.05, distribution is not normal, therefore an analysis of variance can not
#be used

#a kruskal walis test was used to assess the variance in the data
kruskal.test(haplotypes$X1~haplotypes$haplotype.)

#Post-hoc Dunn's test, to look at this variance pairwise
posthoc.kruskal.dunn.test (x=haplotypes$X1, g=haplotypes$haplotype., p.adjust.method="BH")


##################### combining the haplotypes into low and high groups ############################
haplotypes$haplotype. <- as.character(haplotypes$haplotype.)

ligninContent <-  haplotypes %>% 
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 2", as.character('Low'))) %>%
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 5", as.character('Low'))) %>%
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 6", as.character('Low'))) %>%
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 1", as.character('High'))) %>%
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 3", as.character('High'))) %>%
  mutate(haplotype. = replace(haplotype., haplotype. == "Haplotype 4", as.character('High'))) %>%
  cbind(haplotypes[, 1])

names(ligninContent) <- c('lignin level', 'lignin content(%)', 'haplotype')

#boxplot of the data
boxplot(ligninContent$`lignin content(%)`~ligninContent$`lignin level`, xlab = "Lignin Content", 
        ylab = "FT-IR Lignin Content (%)")

#a kruskal walis test was used to assess the variance in the data
kruskal.test(ligninContent$`lignin content(%)`~ligninContent$`lignin level`)

# ***NOTE*** play around with factoring in lignin level and haplotype 


##################### quasi binomial approach ############################

#example of another possible approach 
modeltest <- with(ligninContent, glm((`lignin content(%)`/100)~`lignin level`, quasibinomial))

summary(modeltest)
