library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(magrittr)

#dont need this but took 7 years to download so dont want to get rid of it.
#silva <- read.csv("~/Desktop/16scandidate/silva/SILVA_132_LSUParc_tax_silva.fasta", sep = ";", header = FALSE, stringsAsFactors = FALSE) 
#str(silva)

#silva taxonomy and accessions are in the same file
#rawsilvatax <-  read.csv("~/Downloads/taxmap_ncbi_ssu_ref_nr99_132.txt", sep = ";", header = FALSE, stringsAsFactors = FALSE)
rawsilvatax <-  read.csv("~/Downloads/taxmap_ncbi_ssu_ref_nr99_132.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
silvatax <- read.csv("~/Downloads/taxmap_ncbi_ssu_ref_nr99_132.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
str(rawsilvatax)

#substitute tabs
# silvatax$V1 <- gsub("\t", " ", silvatax$V1)
# silvatax$V6 <- gsub("\t", " ", silvatax$V6)
# silvatax$V7 <- gsub("\t", " ", silvatax$V7)

# remove uncultured
silvatax <- silvatax[silvatax$V5 != "uncultured bacterium", ]
# select prokaryotes
silvatax <- silvatax[grepl("prokaryotes",silvatax$V4), ]

# remove Candidus, etc
silvatax$V5 <- gsub("Candidatus", "", silvatax$V5)
# remove unclassified, etc
silvatax <-silvatax[!grepl("[Uu]nclassified",silvatax$V4), ]

# remove unidentified
silvataxclean <- silvatax[!grepl("[Uu]nidentified",silvatax$V4), ]

#Optional strain
silvataxclean$g <- gsub("(\\S+)\\s(\\S+).*", "\\1", silvataxclean$V5)
silvataxclean$s <- gsub("(\\S+)\\s(\\S+).*", "\\2", silvataxclean$V5)
silvataxclean$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", silvataxclean$V5)
table(silvataxclean$s=="sp.")


#rename V1
silvataxclean$accession <- silvataxclean$V1
silvataxclean$V1 <- NULL

# #filter out everything except for bacteria
# silvataxclean <- silvatax %>%
#   filter(kingdom == "Bacteria") %>%
#   as.data.frame()


  

str(silvataxclean)

#Get rid of empty GenusSpecies and empty accessions
silvataxclean <- silvataxclean %>%
                filter(accession !=" ") %>%
                filter(gs !=" ") %>%
                filter(s != "sp.")
                as.data.frame()

#n for a given genus and species
silvataxclean <- silvataxclean %>% 
  group_by(gs) %>%
  mutate(silva_ngs = n()) %>%
  as.data.frame()

#Table containing just genuspecies and accession.
silvaclean <- silvataxclean[c("accession", "gs")]
write.table(silvaclean, file = "~/Desktop/16scandidate/silva/silvaclean", sep = "\t", col.names = FALSE, row.names = FALSE)

#rrndb
rawrrndb <- read.csv("~/Desktop/16scandidate/rrnDB/rrnDB-5.5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#separate genus and species
rawrrndb$genus <- gsub("(\\S+)\\s(\\S+).*", "\\1", rawrrndb$NCBI.scientific.name)
rawrrndb$species <- gsub("(\\S+)\\s(\\S+).*", "\\2", rawrrndb$NCBI.scientific.name)
rawrrndb$GenusSpecies <- paste(rawrrndb$genus, rawrrndb$species)              

#how many species in both silva and rrndb
table(silvataxclean$gs %in% rawrrndb$GenusSpecies)
table(silvataxclean$Species %in% rawrrndb$species)

# Get mean and sddev for rRNA counts at the genus level
rawrrndb <- rawrrndb %>%
  group_by(genus) %>%
  mutate(mean16s=mean(X16S.gene.count, na.rm = T)) %>%
  mutate(stddev16s=sd(X16S.gene.count, na.rm = T)) %>%
  as.data.frame()

rawrrndb_unique <- rawrrndb[, c("genus", "mean16s", "stddev16s") ]
rawrrndb_unique <- rawrrndb_unique[!duplicated(rawrrndb_unique$genus), ]

table(silvataxclean$g %in% rawrrndb_unique$genus)

#merge
silvarrndb <- inner_join(silvataxclean, rawrrndb_unique, by=c("g"="genus"))

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
silvarrndbsra <- inner_join(silvarrndb, srahits_unique_at_gs, by=c("gs"))

# Sorting by mean 16S gene counts
sixteensilvarrndbsra <- silvarrndbsra %>%
  filter(7 >= mean16s) %>%
  filter(3 <= mean16s) %>%
  filter(1 >= stddev16s | is.na(stddev16s) ) %>%
  as.data.frame()

#file containing just sras.
srapure <- srahits[,c("run_SRAs","gs")]
write.table(srapure, file = "~/Desktop/16scandidate/srafind/srapure", sep = "\t", col.names = FALSE, row.names = FALSE) 

#Select rows of interest for plots etc.
sixteenclean <- sixteensilvarrndbsra[, c("g", "s", "gs", "silva_ngs", "sra_n_gs", "platform" )]
sixteenclean_unique <- sixteenclean[!duplicated(sixteenclean$gs),]

# Add column for the percentile of silva_ngs and sra_n_gs in silva,rrndb,sra merge unique and clean
perc_sra <- data.frame(value=quantile(unique(sixteenclean_unique$sra_n_gs),  probs = seq(0,1, .01)))
perc_silva <- data.frame(value=quantile(unique(sixteenclean_unique$silva_ngs), probs = seq(0,1, .01)))
perc_sra$percentile <-rownames(perc_sra)
perc_silva$percentile <-rownames(perc_silva)

sixteenclean_unique$perc_sra <- NA
sixteenclean_unique$perc_silva <- NA

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
        sixteenclean_unique[r, "silva_ngs"] >= perc_silva[perc-1, "value"] & 
        sixteenclean_unique[r, "silva_ngs"] <= perc_silva[perc, "value"]
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
p <- ggplot(sixteenclean_unique, aes(x=sra_n_gs, y=perc_silva, label=gs)) +
       scale_x_log10() +
         labs(x="sraHits",
              y="silva",
              title="Count of organisms that are present in sraHits, SILVA and rrndb") +
         annotate("rect", xmin = 20, xmax = 65, ymin = 40, ymax = 60,
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
                
      


