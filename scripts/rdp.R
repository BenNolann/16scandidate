library(stringr)

#Using bacteria and archaea header file to categorise, similar to silva and greengenes.
#Read in rdp bacteria header file with tabs
#rdpbacteria2 <- read.csv("~/Desktop/16scandidate/rdp/bacteriaheader", header = FALSE, sep = ";", stringsAsFactors = FALSE)

#Read in rdp bacteria header file with tabs
rdpbacteria <- read.csv("~/Desktop/16scandidate/rdp/bacteriaheader", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#rdpbacteria <- rdpbacteria[sample(nrow(rdpbacteria), 10000), ]
str(rdpbacteria)

rdpbacteria <- rdpbacteria[!grepl(".*uncultured bacterium.*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured organism*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured soil*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*uncultured prokaryote*", rdpbacteria$V1), ]
rdpbacteria <- rdpbacteria[!grepl(".*endophytic bacterium*", rdpbacteria$V1), ]
rdpbacteria$acc <- gsub(">(.*?) .*", "\\1", rdpbacteria$V1) 
rdpbacteria$g <- gsub(">(.*?) (.*?) .*", "\\2", rdpbacteria$V1) 
rdpbacteria$s <- gsub(">(.*?) (.*?) (.*)", "\\3", rdpbacteria$V1) 
rdpbacteria$gs <- gsub(">(.*?) (.*?) (.*)", "\\2 \\3", rdpbacteria$V1) 
# filter out uncultured bacterium and uncultured organisms


#rdpbacteria <- rdpbacteria[!grepl("uncultured", rdpbacteria$V1), ]
#replace ; with spaces
rdpbacteria$V2 <- gsub(";", " ", rdpbacteria$V2)
rdpbacteria$V1 <- gsub(";", " ", rdpbacteria$V1)

#remove subclass and suborder
rdpbacteria$V2 <- gsub("\\s(suborder)", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("(\\S+)\\s(suborder)", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("(\\S+)\\s(subclass)", "", rdpbacteria$V2)


#remove titles of taxonomy
rdpbacteria$V2 <- gsub("domain|phylum|class|order|family|genus|sub", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("", "", rdpbacteria$V2)
# rdpbacteria$V2 <- gsub("1", "", rdpbacteria$V2)

#remove other unhelpful text
rdpbacteria$V2 <- gsub("Lineage=Root", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("rootrank", "", rdpbacteria$V2)
rdpbacteria$V2 <- gsub(" +", " ", rdpbacteria$V2)
rdpbacteria$V2 <- gsub("^ ", "", rdpbacteria$V2)
#rdpbacteria$V3 <- sapply(rdpbacteria$V2, function(x) {length(strsplit(x, " ")[[1]])})
rdpbacteria$V3 <- str_count(rdpbacteria$V2, " ")


rdpbacteria$Domain <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\1", rdpbacteria$V2)
rdpbacteria$Phylum <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\2", rdpbacteria$V2)
rdpbacteria$Class <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\3", rdpbacteria$V2)
rdpbacteria$Order <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\4", rdpbacteria$V2)
rdpbacteria$Family <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\5", rdpbacteria$V2)
rdpbacteria$Genus <- gsub("(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)\\s(\\S+)", "\\6", rdpbacteria$V2)




#Read in rdp archaea header file with tabs
rdparchaea <-  read.csv("~/Desktop/16scandidate/rdp/archaeaheader", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#Read in rdp archaea header file with ;
#rdparchaea2 <-  read.csv("~/Desktop/16scandidate/rdp/archaeaheader", header = FALSE, sep = ";", stringsAsFactors = FALSE)
