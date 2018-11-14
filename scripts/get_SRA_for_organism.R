args = commandArgs(T)

# test args
# args=c("morganii", "Morganella morganii")
# setwd("~/Desktop/16scandidate/scripts/")
print(args)
USAGE <- 'USAGE: Rscript get_sra_for_organism.R  /path/to/output/dir/ "organism name"'
OUTPUT_DIR <- args[1]
ORGANISM <- args[2]
if (length(args) != 2){
  print(USAGE)
  stop("Not enough arguments")
}
if (dir.exists(OUTPUT_DIR)){
  print("Error: Directory exists!")
  stop()  
}else{
  dir.create(OUTPUT_DIR)
} 


#Download reads for a specific organism from NCBI based on results of .R
#INPUT_FILE containing accessions obtained from ggrrndbsra.R result for chosen organism

if (!file.exists("./srapure")){
  stop("srapure not found")
}

filtered_sra_results <- file.path(OUTPUT_DIR, "srafind")
print(paste0("grep '", ORGANISM, "' srapure > ", filtered_sra_results))
system(paste0("grep '", ORGANISM, "' srapure > ", filtered_sra_results)) 
lines <- readLines(filtered_sra_results)
read_set <- 1
for (line in lines){
  dir.create(file.path(OUTPUT_DIR, read_set))
  SUBOUTPUT_DIR_RAW <- file.path(OUTPUT_DIR, read_set, "raw")
  SUBOUTPUT_DIR_DOWN <- file.path(OUTPUT_DIR, read_set, "downsampled")
  linesplit <- gsub("\"", "", line, fixed=T)
  linesplit <- strsplit(linesplit, "\t")[[1]]
  dir.create(SUBOUTPUT_DIR_RAW)
  dir.create(SUBOUTPUT_DIR_DOWN)
  cmd <- paste0("fastq-dump --split-files ", linesplit[1], " -O ", SUBOUTPUT_DIR_RAW)
  print(cmd)
  system(cmd)
  read_number <- 1
  for (file in dir(SUBOUTPUT_DIR_RAW)){
    downcmd <- paste0("seqtk sample -s100 ", file, "1000000 > ", 
                      file.path(SUBOUTPUT_DIR_DOWN, "reads"), read_number, ".fq")
    read_number = read_number + 1
  }
  # downsample reads to 1,000,000

  
  read_set <- read_set + 1
}
