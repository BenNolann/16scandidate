library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(magrittr)

#dont need this but took 7 years to download so dont want to get rid of it.
silva <- read.csv("~/Desktop/16scandidate/silva/SILVA_132_LSUParc_tax_silva.fasta", sep = ";", header = FALSE, stringsAsFactors = FALSE) 
str(silva)

#silva taxonomy and accessions are in the same file
silvatax <- read.csv("~/Desktop/16scandidate/silva/taxmap_embl_lsu_parc_132.txt", sep = ";", header = FALSE, stringsAsFactors = FALSE)
str(silvatax)

#substitute tabs
silvatax$V1 <- gsub("\t", " ", silvatax$V1)
silvatax$V6 <- gsub("\t", " ", silvatax$V6)
silvatax$V7 <- gsub("\t", " ", silvatax$V7)

           
#split first column
silvatax$accession <- gsub("(\\d*)\\s(.*)", "\\1", silvatax$V1)
silvatax$kingdom <- gsub("(\\S*)\\s(\\S*)\\s(\\S*)\\s(\\S*)", "\\4", silvatax$V1)
silvatax$V1 <- NULL


#headers
headers <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession", "Kingdom")
ordered_headers <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Accession")
colnames(silvatax) <- headers
silvatax <- silvatax[, ordered_headers]

#filter out everything except for bacteria
 silvataxclean <- silvatax %>%
                filter(Kingdom == "Bacteria") %>%
                 as.data.frame()
              

str(silvataxclean)

#Optional strain
silvataxclean$Genus <-  gsub("(\\S+)\\s(\\S+).*", "\\1", silvataxclean$Genus)
silvataxclean$Species <-  gsub("(\\S+)\\s(\\S+).*", "\\2", silvataxclean$Species)

#Make gs column ******FIX******
silvataxclean$GenusSpecies <- paste(silvataxclean$Genus, silvataxclean$Species, sep = " ")


#n for a given genus and species
silvataxclean <- group_by(gs) %>%
                mutate(silva_NGS = n()) %>%
                 as.data.frame()



#filter any non accessions
#silvataxbacc <- silvataxbac %>%
#               group_by(accession) %>%
#                filter(accession="\\w")
           
              
