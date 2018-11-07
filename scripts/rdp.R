library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
#Using bacteria and archaea header file to categorise, similar to silva and greengenes.
#Read in rdp bacteria header file with tabs
#rdpbacteria2 <- read.csv("~/Desktop/16scandidate/rdp/bacteriaheader", header = FALSE, sep = ";", stringsAsFactors = FALSE)

#Read in rdp bacteria header file with tabs
rdpbacteria <- read.csv("~/Desktop/16scandidate/rdp/bacteriaheader", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#rdpbacteria <- rdpbacteria[sample(nrow(rdpbacteria), 10000), ]
str(rdpbacteria)

# filter out uncultured bacterium and uncultured organisms and any other outliers
rdpbacteria <- rdpbacteria[!grepl(".*uncultured bacterium.*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured organism*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured soil*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured prokaryote*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*endophytic bacterium*", rdpbacteria$V1), ]
rdpbacteria$V1 <- gsub("(uncultured)\\s", "", rdpbacteria$V1)
rdpbacteria$V1 <- gsub("Candidatus ", "", rdpbacteria$V1)
rdpbacteria$V2 <- gsub("Candidatus ", "", rdpbacteria$V2)
rdpbacteria$acc <- gsub(">(.*?) .*", "\\1", rdpbacteria$V1) 
rdpbacteria$g <- gsub(">(.*?) +(.*?) .*", "\\2", rdpbacteria$V1) 
rdpbacteria$s <- gsub(">(.*?) +(.*?) (.*)", "\\3", rdpbacteria$V1) 
rdpbacteria$gs <- gsub(">(.*?) +(.*?) (.*)", "\\2 \\3", rdpbacteria$V1)


#remove subclass and suborder
#rdpbacteria$V2 <- gsub("\\s(suborder)", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("order;.*;suborder", "order", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("class;.*;subclass", "class", rdpbacteria$V2)


#remove titles of taxonomy
rdpbacteria$V2 <- gsub("domain;|phylum;|class;|order;|family;|genus;|sub;|genus", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("1", "", rdpbacteria$V2)

#remove other unhelpful text
rdpbacteria$V2 <- gsub("Lineage=Root;", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("rootrank;", "", rdpbacteria$V2)
#rdpbacteria$V2 <- gsub(" +", " ", rdpbacteria$V2)
#rdpbacteria$V2 <- gsub("^;+", "", rdpbacteria$V2)
#rdpbacteria$V3 <- sapply(rdpbacteria$V2, function(x) {length(strsplit(x, " ")[[1]])})
rdpbacteria$V3 <- str_count(rdpbacteria$V2, ";")

#histogram showing number of words in lineage 
hist(rdpbacteria$V3)

#lineage must contain 6 words
rdpbacteria <- rdpbacteria %>%
  filter(V3 =="6") %>%
  as.data.frame()

#Splitting up V2
rdpbacteria$Domain <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\1", rdpbacteria$V2)
rdpbacteria$Phylum <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\2", rdpbacteria$V2)
rdpbacteria$Class <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\3", rdpbacteria$V2)
rdpbacteria$Order <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\4", rdpbacteria$V2)
rdpbacteria$Family <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\5", rdpbacteria$V2)
rdpbacteria$Genus <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\6", rdpbacteria$V2)
rdpbacteria$V2 <- NULL
rdpbacteria$V1 <- NULL
rdpbacteria$V3 <- NULL

#Headers and ordered
headers <- c("Accessions", "g", "s", "gs", "Domain", "Phylum", "Class", "Order", "Family", "Genus")
ordered_headers <- c("Accessions", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "g", "s", "gs")
colnames(rdpbacteria) <- headers
rdpbacteria <- rdpbacteria[, ordered_headers]

# Separate strain from gs, however number of strings is inconsistent and strain isn't always just the last string
# rdpbacteria$gs <- gsub("(*)\\s(\\S+$)", "\\1", rdpbacteria$gs)
# rdpbacteria$g <-  gsub("(\\S+)\\s(\\S+)", "\\1", rdpbacteria$gs)
# rdpbacteria$s <-  gsub("(\\S+)\\s(\\S+)", "\\2", rdpbacteria$gs)

#if 'Genus' is not equal to 'g', get rid
table(rdpbacteria$Genus == rdpbacteria$g)
rdpbacteria <- rdpbacteria[rdpbacteria$Genus == rdpbacteria$g, ]

rdpbacteria <- rdpbacteria[!grepl("^sp.", rdpbacteria$s), ]

rdpbacteria$species <- gsub("([[:lower:]]*).*", "\\1", rdpbacteria$s)

#Read in archaea
rdparchaea <- read.csv("~/Desktop/16scandidate/rdp/archaeaheader", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# take care of those lacking classes\
rdparchaea$V2 <-ifelse(
  !grepl("class;", rdparchaea$V2),
  gsub("phylum;", "phylum;incertae_sedis;class;", rdparchaea$V2),
  rdparchaea$V2   
  )

rdparchaea$V1 <- gsub("Candidatus ", "", rdparchaea$V1)
rdparchaea$V2 <- gsub("Candidatus ", "", rdparchaea$V2)



#Remove tax titles
rdparchaea$V2 <- gsub("domain;|phylum;|class;|order;|family;|genus", "", rdparchaea$V2)

#Remove unhelpful text and count strings
rdparchaea$V2 <- gsub("Lineage=Root;", "", rdparchaea$V2)
rdparchaea$V2 <- gsub("rootrank;", "", rdparchaea$V2)
rdparchaea$V3 <- str_count(rdparchaea$V2, ";")

#histogram to show string number
hist(rdparchaea$V3)

#Split up V1
rdparchaea <- rdparchaea[!grepl(".*uncultured archaeon.*", rdparchaea$V1), ]
rdparchaea <- rdparchaea[!grepl(".*uncultured bacterium.*", rdparchaea$V1), ]
rdparchaea$V1 <- gsub("(uncultured)\\s", "", rdparchaea$V1)
rdparchaea$acc <- gsub(">(.*?) .*", "\\1", rdparchaea$V1) 
rdparchaea$g <- gsub(">(.*?) (.*?) .*", "\\2", rdparchaea$V1) 
rdparchaea$s <- gsub(">(.*?) (.*?) (.*)", "\\3", rdparchaea$V1) 
rdparchaea$gs <- gsub(">(.*?) (.*?) (.*)", "\\2 \\3", rdparchaea$V1)


#Split up V2. About a third contain 5 strings for tax, lacking class. Added in latin for unknown in class.
rdparchaea$Domain <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\1", rdparchaea$V2)
rdparchaea$Phylum <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\2", rdparchaea$V2)
rdparchaea$Class <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\3", rdparchaea$V2)
rdparchaea$Order <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\4", rdparchaea$V2)
rdparchaea$Family <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\5", rdparchaea$V2)
rdparchaea$Genus <-  gsub("(.*?);(.*?);(.*?);(.*?);(.*?);(.*?);", "\\6", rdparchaea$V2)
rdparchaea$V2 <- NULL
rdparchaea$V1 <- NULL
rdparchaea$V3 <- NULL

#lineage must contain 6 words
rdparchaea <- rdparchaea %>%
  filter(V3 =="6") %>%
  as.data.frame()

#If Genus is not equal to g, get rid
table(rdparchaea$Genus == rdparchaea$g)
rdparchaea <- rdparchaea[rdparchaea$Genus == rdparchaea$g, ]
rdparchaea <- rdparchaea[!grepl("^sp.", rdparchaea$s), ]
rdparchaea$species <- gsub("([[:lower:]]*).*", "\\1", rdparchaea$s)

#MERGE
