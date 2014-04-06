#### Concatenates files in the same directory and outputs values for all genes in every feature

args <- commandArgs(TRUE)

print_info = function () {
	cat(sprintf("concatenateCoverageMean.R
Usage: Rscript concatenateCoverageMean.R <pathToFolder> <outputFile>
Concatenates files in the same directory and averages values for all features.
"))
}

if (length(args) == 2) {
	path <- args[1]
	outputFile <- args[2]
} else {
	cat(sprintf("ERROR: Wrong number of arguments provided.
"))
	print_info()
}

# Load libraries
require('gtools')
require('stringr')

filenames = mixedsort(list.files(path)[grep(x=list.files(path), pattern=".bed")])

# find out the length of the files
dataRows <- nrow(read.table(paste(path, filenames[1], sep="/"), sep="\t", header=FALSE))
dataCols <- length(mixedsort(filenames))

concate = as.data.frame(matrix(NA, ncol=dataCols+1, nrow=dataRows))

### add coverage values for remaining bins to table
for(bin in 1:length(filenames)) {
	# read sample bin
	binCur <- read.delim(paste(path, filenames[bin], sep="/"),  sep="\t", header=FALSE)
	colnames(binCur) <- c("chr", "start", "end", "ID", "count")

	binCurName = gsub(x=filenames[bin], pattern="*.bed",replacement="")
	binCurName = gsub(x=binCurName, pattern="win_*",replacement="")

	if (bin == 1) {
		concate[,bin] = as.character(binCur$ID)
		names(concate)[1] <- "gene"
	}

	#Add data to concate and column names (bin identity)
	concate[, bin + 1] <- as.numeric(binCur$count)
	colnames(concate)[bin + 1] <- binCurName
}
# save
save(concate, file = outputFile)