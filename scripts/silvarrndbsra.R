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
rawsilvatax <-  read.csv("~/Desktop/16scandidate/silva/taxmap_embl_lsu_parc_132.txt", sep = ";", header = FALSE, stringsAsFactors = FALSE)
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

#Make gs column 
silvataxclean$GenusSpecies <- paste(silvataxclean$Genus, silvataxclean$Species, sep = " ")


#Get rid of empty GenusSpecies and empty accessions
silvataxclean <- silvataxclean %>%
                filter(Accession !=" ") %>%
                filter(GenusSpecies !=" ") %>%
                as.data.frame()

#n for a given genus and species
silvataxclean <- silvataxclean %>% 
  group_by(GenusSpecies) %>%
  mutate(silva_NGS = n()) %>%
  as.data.frame()

#Table containing just genuspecies and accession.
silvaclean <- silvataxclean[c("Accession", "GenusSpecies")]
write.table(silvaclean, file = "~/Desktop/16scandidate/silva/silvaclean", sep = "\t", col.names = FALSE, row.names = FALSE)

#rrndb
rawrrndb <- read.csv("~/Desktop/16scandidate/rrnDB/rrnDB-5.5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#separate genus and species
rawrrndb$genus <- gsub("(\\S+)\\s(\\S+).*", "\\1", rawrrndb$NCBI.scientific.name)
rawrrndb$species <- gsub("(\\S+)\\s(\\S+).*", "\\2", rawrrndb$NCBI.scientific.name)
rawrrndb$GenusSpecies <- paste(rawrrndb$genus, rawrrndb$species)              

#how many species in both silva and rrndb
table(silvataxclean$GenusSpecies %in% rawrrndb$GenusSpecies)
table(silvataxclean$Species %in% rawrrndb$species)

# Get mean and sddev for rRNA counts at the genus level
rawrrndb <- rawrrndb %>%
  group_by(genus) %>%
  mutate(mean16s=mean(X16S.gene.count, na.rm = T)) %>%
  mutate(stddev16s=sd(X16S.gene.count, na.rm = T)) %>%
  as.data.frame()

rawrrndb_unique <- rawrrndb[, c("genus", "mean16s", "stddev16s") ]
rawrrndb_unique <- rawrrndb_unique[!duplicated(rawrrndb_unique$genus), ]

table(silvataxclean$Genus %in% rawrrndb_unique$genus)

#merge
silvarrndb <- inner_join(silvataxclean, rawrrndb_unique, by=c("Genus"="genus"))

#sraHits
rawsrahits <- read.csv("~/Desktop/16scandidate/srafind/sraFind-Contig-biosample-with-SRA-hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#Optional strain
rawsrahits$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", rawsrahits$organism_ScientificName)

#Remove empty gs, remove long read sequencing platforms
srahits <- rawsrahits %>%
  filter(gs != "NA") %>%
  filter(platform != "OXFORD_NANOPORE") %>%
  filter(platform != "PACBIO_SMRT") %>% 
  as.data.frame()

# Add in extra column for n at a given genus and species
srahits <- srahits %>%
  group_by(gs) %>%
  mutate(sra_n_gs = n()) %>%
  as.data.frame()

#Unique at gs
srahits_unique_at_gs <- srahits[!duplicated(srahits$gs), ]

#merge to compare at gs level between silva, rrndb and sra
silvarrndbsra <- inner_join(silvarrndb, srahits_unique_at_gs, by=c("GenusSpecies"="gs"))

# Sorting by mean 16S gene counts
sixteensilvarrndbsra <- silvarrndbsra %>%
  filter(7 >= mean16s) %>%
  filter(3 <= mean16s) %>%
  filter(1 >= stddev16s | is.na(stddev16s) ) %>%
  as.data.frame()

srapure <- srahits[,c("run_SRAs","gs")]
write.table(srapure, file = "~/Desktop/16scandidate/srafind/srapure", sep = "\t", col.names = FALSE, row.names = FALSE) 

#Select rows of interest for plots etc.
sixteenclean <- sixteensilvarrndbsra[, c("Kingdom", "Phylum", "Class", "Order", "Family" , "Genus", "Species", "GenusSpecies", "silva_NGS", "sra_n_gs", "platform" )]
sixteenclean_unique <- sixteenclean[!duplicated(sixteenclean$GenusSpecies),]

# Add column for the percentile of silva_NGS and sra_n_gs in silva,rrndb,sra merge
perc_sra <- data.frame(value=quantile(unique(sixteensilvarrndbsra$sra_n_gs),  probs = seq(0,1, .01)))
perc_silva <- data.frame(value=quantile(unique(sixteensilvarrndbsra$silva_NGS), probs = seq(0,1, .01)))
perc_sra$percentile <-rownames(perc_sra)
perc_silva$percentile <-rownames(perc_silva)

sixteensilvarrndbsra$perc_sra <- NA
sixteensilvarrndbsra$perc_silva <- NA

#sra
for (perc in 1:nrow(perc_sra)){
  if(perc == 1){
    next()
  } else {
    for (r in 1:nrow(sixteenclean_unique)){
      if(
        sixteenclean_unique[r, "sra_n_gs"] >= perc_sra[perc-1, "value"] & 
        sixteenclean_unique[r, "sra_n_gs"] <= perc_sra[perc, "value"]
      ){
        
        sixteenclean_unique[r, "perc_sra"] <- perc_sra[perc-1, "percentile"]
      }
    }
  }
}

#silva
for (perc in 1:nrow(perc_silva)){
  if(perc == 1){
    next()
  } else {
    for (r in 1:nrow(sixteenclean_unique)){
      if(
        sixteenclean_unique[r, "silva_NGS"] >= perc_silva[perc-1, "value"] & 
        sixteenclean_unique[r, "silva_NGS"] <= perc_silva[perc, "value"]
      ){
        
        sixteenclean_unique[r, "perc_silva"] <- perc_silva[perc-1, "percentile"]
      }
    }
  }
}

#Convert percentile to integers
sixteenclean_unique$perc_silva <- as.numeric(gsub("%", "", sixteenclean_unique$perc_silva))
sixteenclean_unique$perc_sra <-  as.numeric(gsub("%", "", sixteenclean_unique$perc_sra))

#To label the plot
#sixteenclean_unique$labelGenusSpecies <- ifelse(sixteenclean_unique$labelGenusSpecies == , sixteenclean_unique$labelGenusSpecies, "")
#Lets make a plot- basis
p <- ggplot(sixteenclean_unique, aes(x=sra_n_gs, y=perc_silva, label=GenusSpecies)) +
       scale_x_log10() +
         labs(x="sraHits",
              y="silva",
              title="Count of organisms that are present in sraHits, SILVA and rrndb") +
         annotate("rect", xmin = 1, xmax = 75, ymin = 20, ymax = 60,
                  alpha = .2)

#-Actual plot
(ppoint <- p +
          geom_point() +
          theme(axis.line.x = element_line(colour="black", size=1),
                axis.line.y = element_line(colour="black", size=1),
                panel.background=element_rect(fill="transparent"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) +
           geom_label_repel())
                
      


