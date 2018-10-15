rawrrndb <- read.csv("~/Desktop/16scandidate/rrnDB/rrnDB-5.5.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(rawrrndb)
colnames(rawrrndb)
colnames(both_acc_and_tax)
rownames(rawrrndb)
headrrndb <- head(rawrrndb, 2000)
headboth_acc_and_tax <- head(both_acc_and_tax, 2000)
rrndb_gg_merge <- merge(rawrrndb, both_acc_and_tax, by.x = "", by.y = "") #Note: NCBI.tax.id in rrnDB and id in greengenes are different
# What will be the link to merge the two tables? The genus and species name 
sum(rrndb$BioSample !="")

#   genus  species [optional strain]
# "(\\S+)\\s(\\S+).*"
rawrrndb$species_species <- gsub("(\\S+)\\s(\\S+).*", "\\1 \\2", rawrrndb$NCBI.scientific.name)

