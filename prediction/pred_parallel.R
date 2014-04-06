#### Predicts nucleossome positioning for all fasta files in <pathToFolder>

args <- commandArgs(TRUE)

print_info = function () {
	cat(sprintf("concatenateCoverageMean.R
Usage: Rscript concatenateCoverageMean.R <pathToFolder> <outputDir>
Predicts nucleossome positioning for all fasta files in <pathToFolder>
"))
}

if (length(args) == 2) {
	path <- args[1]
	outputDir <- args[2]
} else {
	cat(sprintf("ERROR: Wrong number of arguments provided.
"))
	print_info()
}

# Load libraries
library(foreach)
library(doMC)
registerDoMC(24) # number of cores to use
if ("rtracklayer" %in% rownames(installed.packages())) {
	library("NuPoP")
} else {
	source("http://bioconductor.org/biocLite.R")
	biocLite("NuPoP")
}

files = list.files(path)
fastafiles = files[grep(pattern="*.fasta", x=files)]
len = length(fastafiles)

foreach(i=1:len) %dopar% {
	predName = paste(fastafiles[i], "Prediction4.txt", sep="_")
	if (file.exists(predName)){
		cat(sprintf("Prediction exists, skipping...
"))
	} else {
		predNuPoP(fastafiles[i], species=0, model=4)
		predLength = length(readLines(predName))
		pred = readNuPoP(predName, start = 1, end = predLength)
		scaffoldName = sub("^([^.]*).*", "\\1", fastafiles[i])
		write.table(pred, paste(outputDir, scaffoldName, ".nucl", sep=""), sep="\t", quote=FALSE)	
	}
}