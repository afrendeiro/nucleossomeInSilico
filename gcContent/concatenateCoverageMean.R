#### Concatenates files in the same directory and averages values for all features

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

concateMean = as.data.frame(matrix(NA, ncol=dataCols, nrow=1))

### add coverage values for remaining bins to table
for(bin in 1:length(filenames)) {
	# read sample bin
	binCur <- read.delim(paste(path, filenames[bin], sep="/"),  sep="\t", header=FALSE)
	colnames(binCur) <- c("chr", "start", "end", "ID", "gcContent")

	binCurName = gsub(x=filenames[bin], pattern="*.bed",replacement="")
	binCurName = gsub(x=binCurName, pattern="win_*",replacement="")

	# Add data to concateMean and column names (bin identity)
	concateMean[,bin] <- mean(as.numeric(binCur$gcContent))
	colnames(concateMean)[bin] <- binCurName
}
# save
save(concateMean, file=outputFile)