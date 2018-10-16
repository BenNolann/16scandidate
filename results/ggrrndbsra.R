#install.packages("tidyverse")
#install.packages("data.table")
#library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

rawtax <- read.csv("~/Desktop/16scandidate/greengenes/gg_13_5_taxonomy.txt", sep=";", stringsAsFactors = FALSE, header = FALSE)
rawtax$id<- gsub("(\\d*)\\t(.*)", "\\1", rawtax$V1)
rawtax$k<-  gsub("(\\d*)\\t(.*)"," \\2", rawtax$V1)
rawtax$V1 <- NULL
headers <- c("p", "c", "o", "f", "g", "s", "id", "k")
ordered_headers <- c("k", "p", "c", "o", "f", "g", "s", "id")
colnames(rawtax) <- headers
tax <- rawtax[, ordered_headers]

#rawtax <- read.csv("~/Desktop/16scandidate/greengenes/gg_13_5_taxonomy.txt", sep=";", stringsAsFactors = FALSE, header = FALSE)
#tax <- head(rawtax, 2000)

head(rawtax)
# split fist column  so that id is everything before tab and k is everything after

for (column in ordered_headers){
 tax[, column]<- gsub(" .__", "",tax[, column])
}
tax$id <- as.numeric(tax$id)

#########  exploring the data
#results_df <- starting_dataframe %>% DO_SOME_STUFF()
lok <- tax %>% 
  filter(str_detect(g, "Loktanella")) 

counts <- tax %>%
  group_by(k, p, c, o, f,  g, s) %>%
  mutate(n=n()) %>%
  distinct( k, p, n, c, o, f,  g, s)%>%
  as.data.frame()

#counts <- counts[!duplicated(counts),]
counts$name <- paste(counts$k, counts$p, counts$c, counts$o, counts$f, counts$g, counts$s)

 tmp <-tax %>%
  filter(s!="") %>%
   filter( g!="") %>%
  group_by(k, p, c, o, f,  g, s) %>%
  mutate(n=n()) %>%
  unite(name, k,p,c,o,f,g,s,n, remove = F) %>%
   distinct( name, n)%>%
   as.data.frame() 
ggplot(tmp, aes(x=reorder(s, -n), y=n, fill=k)) + geom_bar(stat="identity", color="black") +
    scale_y_continuous() + facet_grid(~k, scales = "free") +
    theme(
      axis.text.x = element_text(angle=45, hjust = 1)
  )


counts[counts$n == max(counts$n), ]
counts %>%
  group_by(k)%>%
  filter(s!="") %>%
  filter(n==max(n))

counts %>%
  filter(p=="") %>%
  filter(n==max(n))


table(counts$k)
#tax[tax$id==228054 & , "k"]


#gsub("\t"," ", tax)
rawacc <- read.csv("~/Desktop/16scandidate/greengenes/gg_13_5_accessions.txt", sep ="\t", stringsAsFactors = FALSE, header = TRUE)
str(rawacc)
table(rawacc$accession_type)

# check all tax ids in accessions
table(
  tax$id %in% rawacc$X.gg_id
)
# find which line taxid 4 is
grep("^755605$", rawacc$X.gg_id)


#  merging data
#a <- head(tax, 1999)
#b <- head(rawacc, 1999)
#colnames(a)
#colnames(b)
both_acc_and_tax <- merge(tax, rawacc, by.x = "id", by.y = "X.gg_id")
colnames(both_acc_and_tax)

both_acc_and_tax$gs <- paste(both_acc_and_tax$g, both_acc_and_tax$s, sep = " ") #combine g and s
both_acc_and_tax$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", both_acc_and_tax$gs) # Optional strain
# add column for n at a given genus and species
both_acc_and_tax <- both_acc_and_tax %>%
  group_by(gs)%>%
  mutate(n_gs = n()) %>%
  as.data.frame()

#####  filter out those lacking species and genus
both_acc_and_tax_filtered <- both_acc_and_tax %>%
  filter(g != "") %>%
  filter(s != "") %>%
  as.data.frame()


ggplot(both_acc_and_tax_filtered, aes(x=n_gs)) + geom_histogram(bins = 100)
ggplot(both_acc_and_tax_filtered, aes(y=n_gs)) + geom_boxplot() +geom_point(aes(x=0))
ggplot(both_acc_and_tax_filtered, aes(y=n_gs)) + geom_violin(aes(x=0))

head(both_acc_and_tax_filtered$gs)
#rrnDB_tidy.R was moved here
rawrrndb <- read.csv("~/Desktop/16scandidate/rrnDB/rrnDB-5.5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(rawrrndb)
colnames(rawrrndb)
colnames(both_acc_and_tax_filtered)
rownames(rawrrndb)



#   genus  species [optional strain]
# "(\\S+)\\s(\\S+).*"
rawrrndb$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", rawrrndb$`NCBI.scientific.name`)


table(rawrrndb$gs %in% both_acc_and_tax_filtered$gs)
table(both_acc_and_tax_filtered$gs %in% rawrrndb$gs)


rrndb_gg_merge <- inner_join(rawrrndb[!duplicated(rawrrndb$gs), ],
                              both_acc_and_tax_filtered[!duplicated(both_acc_and_tax_filtered$gs), ], 
                              by = "gs" )

#srahits_tidy.R

rawsrahits <- read.csv("~/Desktop/16scandidate/srafind/sraFind-Contig-biosample-with-SRA-hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rawsrahits$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", rawsrahits$organism_ScientificName)


srahits <- read.csv("~/Desktop/16scandidate/srafind/sraFind-Contig-biosample-with-SRA-hits.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
srahits$gs <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", srahits$organism_ScientificName)
srahits <- srahits %>%
  filter(gs != "NA") %>%
  as.data.frame()

# Add in extra column for n at a given genus and species
rawsrahits <- rawsrahits %>%
    group_by(gs) %>%
      mutate(n_gs = n()) %>%
        as.data.frame()
  


sra_rrndb_gg_merge <- merge(rrndb_gg_merge, srahits[!duplicated(srahits$gs),], by = "gs")


tmp <- sra_rrndb_gg_merge %>%
  filter(7 >= X16S.gene.count) %>%
  filter(4 <= X16S.gene.count) %>%
         as.data.frame()
  
###### Calculating the percentile of the n_gs in sra_rrndb_gg_merge
quantile(sra_rrndb_gg_merge$n_gs.y, probs = c(.1 ,.2, .3, .4, .5, .6, .7, .8 ,.9, 1))
quantile(sra_rrndb_gg_merge$n_gs.x, probs = c(.1 ,.2, .3, .4, .5, .6, .7, .8 ,.9, 1))




#sra_rrndb_gg_merge$percentilegg <- quantile(sra_rrndb_gg_merge$n_gs.x)
#sra_rrndb_gg_merge$percentilesra <- quantile(sra_rrndb_gg_merge$n_gs.y)
