#install.packages("tidyverse")
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
dir.create(camp_find_gg)
tax <- read.csv("~/Downloads/gg_13_5_taxonomy.txt", sep=";", stringsAsFactors = FALSE, header = FALSE)
#tax <- head(rawtax, 2000) 
head(tax)
# split fist column  so that id is everything before tab and k is everything after
tax$id<- gsub("(\\d*)\\t(.*)", "\\1", tax$V1)
tax$k<-  gsub("(\\d*)\\t(.*)"," \\2", tax$V1)
tax$V1 <- NULL

headers <- c("p", "c", "o", "f", "g", "s", "id", "k")
ordered_headers <- c("k", "p", "c", "o", "f", "g", "s", "id")
colnames(tax) <- headers
tax <- tax[, ordered_headers]
for (column in ordered_headers){
 tax[,column]<- gsub(" .__", "", tax[,column])
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
acc <- read.csv("~/Downloads/gg_13_5_accessions.txt", sep ="\t", stringsAsFactors = FALSE, header = TRUE)
str(acc)
table(rawacc$accession_type)

# check all tax ids in accessions
table(
  tax$id %in% rawacc$X.gg_id
)
# find which line taxid 4 is
grep("^755605$", rawacc$X.gg_id)


#  merging data
a <- head(tax, 1999)
b <- head(acc, 1999)
colnames(a)
colnames(b)
both_acc_and_tax <- merge(a, b, by.x = "id", by.y = "X.gg_id")
