args = commandArgs(T)

# test args
# args=c("morganii", "/path/to/riboseed/git/")
#setwd("~/Desktop/16scandidate/scripts/")
print(args)
USAGE <- 'USAGE: Rscript run_pyani_on_SRAs.R  /path/to/output/dir/ /path/to/riboseed/git/'
INPUT_DIR <- args[1]
RIBOSEED_DIR <- args[2]

# list subdirectories
subdirectories <- list.dirs(INPUT_DIR, recursive = F )
 
for (d in subdirectories){
  name <- paste0(basename(INPUT_DIR), "_", basename(d))
  OUTPUT_DIR <- file.path(INPUT_DIR, d,  "pyani")
  reads <- list.files(d)
  if (length(reads) == 2){
    cwd <- getwd()
    cmd <- paste("./run_ANI.sh -e", name, "-f", reads[1], "-r", reads[2], "-o", file.path(cwd, OUTPUT_DIR ))
    print(cmd)
    setwd(file.path(RIBOSEED_DIR, "scripts"))
    system(cmd)
    setwd(cwd)
  }
}

