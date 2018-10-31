#install.packages("tidyverse")
#install.packages("data.table")
#library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

# Read in greengenes taxa.
rawtax <- read.csv("~/Desktop/16scandidate/greengenes/gg_13_5_taxonomy.txt", sep=";", stringsAsFactors = FALSE, header = FALSE)
# split fist column  so that id is everything before tab and k is everything after
rawtax$id<- gsub("(\\d*)\\t(.*)", "\\1", rawtax$V1)
rawtax$k<-  gsub("(\\d*)\\t(.*)"," \\2", rawtax$V1)
rawtax$V1 <- NULL
headers <- c("p", "c", "o", "f", "g", "s", "id", "k")
ordered_headers <- c("k", "p", "c", "o", "f", "g", "s", "id")
colnames(rawtax) <- headers
tax <- rawtax[, ordered_headers]

# replace underscores in rows with no space.
for (column in ordered_headers){
 tax[, column]<- gsub(" .__", "",tax[, column])
}
tax$id <- as.numeric(tax$id)

# read in accessions from gg 
rawacc <- read.csv("~/Desktop/16scandidate/greengenes/gg_13_5_accessions.txt", sep ="\t", stringsAsFactors = FALSE, header = TRUE)
str(rawacc)
table(rawacc$accession_type)

# check all tax ids in accessions
table(tax$id %in% rawacc$X.gg_id)

# merge taxa and accessions by id.
both_acc_and_tax <- merge(tax, rawacc, by.x = "id", by.y = "X.gg_id")
colnames(both_acc_and_tax)

#combine g and s
both_acc_and_tax$gs <- paste(both_acc_and_tax$g, both_acc_and_tax$s, sep = " ") 
# Optional strain
both_acc_and_tax$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", both_acc_and_tax$gs) 
#for use in a bash loop to subset the fasta file
ggped <- both_acc_and_tax[both_acc_and_tax$gs == "Pediococcus acidilactici", c("id", "accession")] 

write.table(ggped, file = "~/Desktop/16scandidate/results/pedidacc", sep = "\t", col.names = FALSE, row.names = FALSE)



write.table(greenidacc, file = "~/Desktop/16scandidate/greengenes/greenidacc", sep = "\t", col.names = FALSE, row.names = FALSE)


# add column for n at a given genus and species
both_acc_and_tax <- both_acc_and_tax %>%
  group_by(gs)%>%
  mutate(gg_n_gs = n()) %>%
  as.data.frame()

#  filter out those lacking species and genus
both_acc_and_tax_filtered <- both_acc_and_tax %>%
  filter(g != "") %>%
  filter(s != "") %>%
  as.data.frame()

# All greengenes id, accession, gs
greenclean <- both_acc_and_tax_filtered[c("id", "accession", "gs")]

write.table(greenclean, file = "~/Desktop/16scandidate/greengenes/greenidacc", sep = "\t", col.names = FALSE, row.names = FALSE)

#ggplot(both_acc_and_tax_filtered, aes(x=gg_n_gs)) + geom_histogram(bins = 100)
#ggplot(both_acc_and_tax_filtered, aes(y=gg_n_gs)) + geom_boxplot() +geom_point(aes(x=0))
#ggplot(both_acc_and_tax_filtered, aes(y=gg_n_gs)) + geom_violin(aes(x=0))

head(both_acc_and_tax_filtered$gs)


#rrnDB_tidy.R was moved here
rawrrndb <- read.csv("~/Desktop/16scandidate/rrnDB/rrnDB-5.5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(rawrrndb)
colnames(rawrrndb)
colnames(both_acc_and_tax_filtered)
rownames(rawrrndb)


#   genus  species [optional strain]
# "(\\S+)\\s(\\S+).*"
rawrrndb$genus <- gsub("(\\S+)\\s(\\S+).*", "\\1", rawrrndb$`NCBI.scientific.name`)
rawrrndb$species <- gsub("(\\S+)\\s(\\S+).*", "\\2", rawrrndb$`NCBI.scientific.name`)
rawrrndb$gs <- paste(rawrrndb$genus, rawrrndb$species)

# Get mean and sddev for rRNA counts at the genus level
rawrrndb <- rawrrndb %>%
  group_by(genus) %>%
  mutate(mean16s=mean(X16S.gene.count, na.rm = T)) %>%
  mutate(stddev16s=sd(X16S.gene.count, na.rm = T)) %>%
  as.data.frame()


table(rawrrndb$genus %in% both_acc_and_tax_filtered$g)
#table(rawrrndb$gs %in% both_acc_and_tax_filtered$gs)
# rawrrndb$gs[!rawrrndb$gs %in% both_acc_and_tax_filtered$gs]
# both_acc_and_tax_filtered$gs[!both_acc_and_tax_filtered$gs %in% rawrrndb$gs]
table(both_acc_and_tax_filtered$gs %in% rawrrndb$gs)

rawrrndb_unique <- rawrrndb[, c("genus", "mean16s", "stddev16s") ]
rawrrndb_unique <- rawrrndb_unique[!duplicated(rawrrndb_unique$genus), ]
table(both_acc_and_tax_filtered$g %in% rawrrndb_unique$genus)
rrndb_gg_merge <- inner_join(
  both_acc_and_tax_filtered, 
  rawrrndb_unique,
 by=c("g" = "genus"))

# Read in sraFind file
rawsrahits <- read.csv("~/Desktop/16scandidate/srafind/sraFind-Contig-biosample-with-SRA-hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rawsrahits$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", rawsrahits$organism_ScientificName)


srahits <- read.csv("~/Desktop/16scandidate/srafind/sraFind-Contig-biosample-with-SRA-hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
srahits$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", srahits$organism_ScientificName)
srahits <- srahits %>%
  filter(gs != "NA") %>%
  filter(platform != "OXFORD_NANOPORE") %>%
  filter(platform != "PACBIO_SMRT") %>% 
  as.data.frame()

# Add in extra column for n at a given genus and species
srahits <- srahits %>%
    group_by(gs) %>%
    mutate(sra_n_gs = n()) %>%
    as.data.frame()
  

srahits_unique_at_gs <- srahits[!duplicated(srahits$gs), ]
#unique(srahits$gs)

#To compare n for given species on sra and gg
sra_rrndb_gg_merge <- inner_join(
  rrndb_gg_merge, 
  srahits_unique_at_gs, by = "gs")
table(rrndb_gg_merge$gs %in% srahits_unique_at_gs$gs)


# Sorting by mean 16S gene counts
srarrndbgg16 <- sra_rrndb_gg_merge %>%
  filter(7 >= mean16s) %>%
  filter(3 <= mean16s) %>%
  filter(1 >= stddev16s | is.na(stddev16s) ) %>%
  as.data.frame()
  
# select rows of interest for plotting, etc
srarrndbgg16_unique <- srarrndbgg16[, c("k", "p", "c", "o", "f" , "g", "s", "gs", "gg_n_gs", "sra_n_gs", "platform" ) ]
srarrndbgg16_unique <- srarrndbgg16_unique[!duplicated(srarrndbgg16_unique$gs),]

# Add column for the percentile of sra_n_gs and gg_n_gs in sra_gg_merge.
perc_df_sra <- data.frame(value=quantile(unique(srarrndbgg16$sra_n_gs),  probs = seq(0,1, .01)))
perc_df_gg <- data.frame(value=quantile(unique(srarrndbgg16$gg_n_gs), probs = seq(0,1, .01)))
perc_df_sra$percentile <-rownames(perc_df_sra)
perc_df_gg$percentile <-rownames(perc_df_gg)
srarrndbgg16$perc_sra <- NA
srarrndbgg16$perc_gg <- NA
# for sra
for (perc in 1:nrow(perc_df_sra)){
  if(perc == 1){
    next()
  } else {
    for (r in 1:nrow(srarrndbgg16_unique)){
      if(
        srarrndbgg16_unique[r, "sra_n_gs"] >= perc_df_sra[perc-1, "value"] & 
        srarrndbgg16_unique[r, "sra_n_gs"] <= perc_df_sra[perc, "value"]
          ){
        
        srarrndbgg16_unique[r, "perc_sra"] <- perc_df_sra[perc-1, "percentile"]
      }
    }
  }
}

# and greengenes
for (perc in 1:nrow(perc_df_gg)){
  if(perc == 1){
    next()
  } else {
    for (r in 1:nrow(srarrndbgg16_unique)){
      if(
        srarrndbgg16_unique[r, "gg_n_gs"] >= perc_df_gg[perc-1, "value"] & 
        srarrndbgg16_unique[r, "gg_n_gs"] <= perc_df_gg[perc, "value"]
      ){
        
        srarrndbgg16_unique[r, "perc_gg"] <- perc_df_gg[perc-1, "percentile"]
      }
    }
  }
}


# convert percentile to integers
srarrndbgg16_unique$perc_gg <- as.numeric(gsub("%", "", srarrndbgg16_unique$perc_gg))
srarrndbgg16_unique$perc_sra <- as.numeric(gsub("%", "", srarrndbgg16_unique$perc_sra))
# install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
# For the label on the plot
srarrndbgg16_unique$labelgs <- ifelse(srarrndbgg16_unique$gs == "Pediococcus acidilactici", srarrndbgg16_unique$gs, "")
g <- ggplot(srarrndbgg16_unique, aes(x=sra_n_gs, y=perc_gg, label=labelgs)) +
 # g <- ggplot(srarrndbgg16_unique, aes(x=sra_n_gs, y=gg_n_gs, label=labelgs)) +
  scale_x_log10() +
  #scale_y_log10() +
  labs(x="sraFind hit",
           y="Greengenes Percentile",
           colour="Genus",
           title="Count distribution of species",
           subtitle="Comparing counts of species based on Greengenes and sraFind \n in order to find a candidate.") +
        annotate("rect", xmin = 1, xmax = 50, ymin = 40, ymax = 60,
                 alpha = .2) 


(glight <- g + 
    geom_point() + 
    theme(axis.line.x = element_line(colour="black", size=1),
           axis.line.y = element_line(colour="black", size=1),
           panel.background=element_rect(fill="transparent"),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank()) +
       geom_label_repel() )

(gdark <- g + 
    geom_point(color="white", size=2) + 
    theme(text = element_text(colour = "white"),
          rect = element_rect(fill = ""),
                    axis.text = element_text(colour = "white"),
                    axis.ticks = element_line(colour = "white"),
                    axis.line.x = element_line(colour="white", size=1),
                    axis.line.y = element_line(colour="white", size=1),
                    panel.background=element_rect(fill="transparent"),
              plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank()) +
    scale_color_manual(values = c("white")) +
  geom_label_repel(segment.colour = "white", color="white", fill="transparent") )

# pdf(file = "./Desktop/tmp_dark.pdf", width = 7, height = 4, bg = "")
# print(gdark)
# dev.off()
ggsave(glight, filename = "./Desktop/tmp_light.pdf",  bg = "transparent")
ggsave(gdark, filename = "./Desktop/tmp_dark.png",  bg = "transparent",width = 7, height = 4.5)

 


