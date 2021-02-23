
# librarys required 
library(magrittr)
library(gt)
library(readr)
library(qpcR)

#setting the working environment
setwd("~/Documents/4th Year/Honours Project/Data/barley genes")

#reading in the files
CCR2 <- read_delim("data", "\t", escape_double = FALSE, 
                   col_types = cols(position = col_character(), 
                                    position_1 = col_character(), SNP_2 = col_character(), 
                                    position_2 = col_character(), position_3 = col_character(), 
                                    position_4 = col_character()), na = "empty", 
                   trim_ws = TRUE, skip = 1)

KNAT <- read_delim("data2", "\t", escape_double = FALSE, 
                   col_types = cols(position = col_character(), 
                                    position_1 = col_character(), SNP_2 = col_character(), 
                                    position_2 = col_character(), position_3 = col_character(), 
                                    position_4 = col_character(), position_5 = col_character()), 
                   na = "empty", trim_ws = TRUE, skip = 1)

CCR2LigninContent <- read_delim("ccr2 lignin content genes", 
                                "\t", escape_double = FALSE, 
                                col_types = cols(position = col_character(), position_1 = col_character()),
                                na = "empty", trim_ws = TRUE, skip = 1)
  
KNATLigninContent <- read_delim("knat lignin content genes", 
                                "\t", escape_double = FALSE,
                                col_types = cols(position = col_character(), position_1 = col_character()),
                                na = "empty", trim_ws = TRUE, skip = 1)

#creating the highlignin content data frame 
ligninlevel <- qpcR:::cbind.na(CCR2LigninContent[, c("SNP_1", "position_1")],
                               CCR2LigninContent[, c("SNP", "position")],
                               KNATLigninContent[, c("SNP_1", "position_1")], 
                               KNATLigninContent[, c("SNP", "position")]) 
# using the cbind.fill function to join to columns with different lengths
ligninlevel[is.na(ligninlevel)] <- '' 
# replacing NAs with an emmpty cell
names(ligninlevel) <- c('SNP', 'position', 'SNP_1', 'position_1', 'SNP_2', 'position_2', 'SNP_3', 'position_3')
#renaming the columns 


#plotting the tables
CCR2%>% gt() %>%
  tab_spanner(
    label = "Haplotype 2",
    columns = vars('position_1','SNP_1')
    ) %>%
  tab_spanner(
    label = "Haplotype 1",
    columns = vars('position','SNP')
  ) %>%
  tab_spanner(
    label = "Haplotype 3",
    columns = vars('position_2','SNP_2')
  ) %>%
  tab_spanner(
    label = "Haplotype 4",
    columns = vars('position_3','SNP_3')
  ) %>%
  tab_spanner(
    label = "Haplotype 5",
    columns = vars('position_4','SNP_4')
  ) %>%
  tab_spanner(
    label = "Haplotype 6",
    columns = vars('position_5','SNP_5')
  ) %>%
  cols_label(
    position = "position",
    SNP = "SNP",
    position_1 = "position",
    SNP_1 = "SNP",
    position_2 = "position",
    SNP_2 = "SNP",
    position_3 = "position",
    SNP_3 = "SNP",
    position_4 = "position",
    SNP_4 = "SNP",
    position_5 = "position",
    SNP_5 = "SNP"
    )

KNAT%>% gt() %>%
  tab_spanner(
    label = "Haplotype 2",
    columns = vars('position_1','SNP_1')
  ) %>%
  tab_spanner(
    label = "Haplotype 1",
    columns = vars('position','SNP')
  ) %>%
  tab_spanner(
    label = "Haplotype 3",
    columns = vars('position_2','SNP_2')
  ) %>%
  tab_spanner(
    label = "Haplotype 4",
    columns = vars('position_3','SNP_3')
  ) %>%
  tab_spanner(
    label = "Haplotype 5",
    columns = vars('position_4','SNP_4')
  ) %>%
  tab_spanner(
    label = "Haplotype 6",
    columns = vars('position_5','SNP_5')
  ) %>%
  cols_label(
    position = "position",
    SNP = "SNP",
    position_1 = "position",
    SNP_1 = "SNP",
    position_2 = "position",
    SNP_2 = "SNP",
    position_3 = "position",
    SNP_3 = "SNP",
    position_4 = "position",
    SNP_4 = "SNP",
    position_5 = "position",
    SNP_5 = "SNP"
  ) 

ligninlevel%>% gt() %>%
  tab_spanner(
    label = html("Low Lignin <br> CCR2 SNPs"),
    columns = vars('position_1','SNP_1')
  ) %>%
  tab_spanner(
    label = html("High Lignin <br> CCR2 SNPs"),
    columns = vars('position','SNP')
  ) %>%
  tab_spanner(
    label = html("Low Lignin <br> KNAT SNPs"),
    columns = vars('position_3','SNP_3')
  ) %>%
  tab_spanner(
    label = html("High Lignin <br> KNAT SNPs"),
    columns = vars('position_2','SNP_2')
  ) %>%
  cols_label(
    position = "position",
    SNP = "SNP",
    position_1 = "position",
    SNP_1 = "SNP",
    position_2 = "position",
    SNP_2 = "SNP",
    position_3 = "position",
    SNP_3 = "SNP"
  ) 

