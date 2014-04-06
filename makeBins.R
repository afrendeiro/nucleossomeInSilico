##### Makes bins around TSS/TES with the required width/length/step

print_info = function () {
	cat(sprintf("
makeBins.R
Usage: Rscript makeBins.R <TSS/TES> <distance> <width> <step> <outputFolder>
Makes bins around TSS/TES with the required width/length/step.
"))
	quit(save="no", status=1)
}

# Parse command-line arguments
args <- commandArgs(TRUE)

if (length(args) == 5) {
	features <- c("TSS", "TES")
	feature = args[1]
	if (!feature %in% features){
		cat(sprintf("ERROR: Unrecognized feature provided."))
		print_info()
	}
	distance = as.numeric(args[2])
	width = as.numeric(args[3])
	step = as.numeric(args[4])
	if (!distance %% 1 == 0 | !width %% 1 == 0 | !step %% 1 == 0) {
		cat(sprintf("ERROR: Provided window values are not integers."))
		print_info()
	}
	outputFolder = args[5]
	dir.create(outputFolder, recursive=TRUE)
} else {
	cat(sprintf("ERROR: Wrong number of arguments provided."))
	print_info()
}

# Load libraries
library(foreach)
library(doMC)
registerDoMC(24) # number of cores to use
options(scipen=999) # avoid output of scientific notation
if ("rtracklayer" %in% rownames(installed.packages())) {
	library("rtracklayer")	
} else {
	source("http://bioconductor.org/biocLite.R")
	biocLite("rtracklayer")
}

makeBins = function(feature, bin, left, right) {
	# Initialize variables
	chr2 = character(0)
	start = integer(0)
	end = integer(0)
	ID = character(0)

	# Do the plus strand genes
	genes.positive = genes[strand(genes) == "+" , ]
	chrs = unique(as.character(chrom(genes.positive)))

	# Go by chromossomes
	for (chr in chrs){
		# Select genes in this chromossome
		genes.chr = genes.positive[chrom(genes.positive) == chr , ]
		# Get chromossome size
		chr.size = chrSizes[which(chr == chrSizes$chr),2]

		for (i in 1:length(genes.chr)) {
			# Discard genes too close to chromossome start/end
			if (start(genes.chr)[i] > chrBuffer && end(genes.chr)[i] < chr.size - chrBuffer) {
				if (feature == "TSS"){
					T = start(genes.chr)[i]
				} else {
					T = end(genes.chr)[i]	
				}
				chr2 <- c(chr2, chr)
				start <- c(start, T + left)
				end <- c(end, T + right)
				ID <- c(ID, as.character(genes.chr$ID[i]))
			}
		}
	}

	# Do the minus strand genes
	genes.negative = genes[strand(genes) == "-" , ]
	chrs = unique(as.character(chrom(genes.negative)))
	
	for (chr in chrs){
		# Select genes in this chromossome
		genes.chr = genes.negative[chrom(genes.negative) == chr , ]
		# Get chromossome size
		chr.size = chrSizes[which(chr == chrSizes$chr),2]

		for (i in 1:length(genes.chr)) {
			# Discard genes too close to chromossome start/end
			if (start(genes.chr)[i] > chrBuffer && end(genes.chr)[i] < chr.size - chrBuffer) {
				if (feature == "TSS"){
					T = end(genes.chr)[i]
				} else {
					T = start(genes.chr)[i]	
				}
				chr2 <- c(chr2, chr)
				start <- c(start, T - right)
				end <- c(end, T - left)
				ID <- c(ID, as.character(genes.chr$ID[i]))
			}
		}
	}

	outputFile = paste(outputFolder, "/win_", bin, ".bed", sep="")
	output = data.frame(chr2, start, end, ID)

	write.table(output, file=outputFile, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)

	cat(sprintf(paste("Finished bin ",bin, ". Saved in ", outputFile, "
", sep="")))
}

bins = seq(-distance, distance, step)
len = length(bins)

# load gene annotation
gff = import("~/data/oikopleura/annotation/Oikopleura_gene_models.reference.gff3")
genes <<- gff[gff$type == "gene", ]

# read chromossome length info
chrSizes <<- read.delim("/home/s3/afr/data/oikopleura/assembly/Oikopleura_reference_chrSizes.tsv", sep="\t",header=F)
colnames(chrSizes) <- c("chr", "size")

# Buffer distance from chromossome ends
chrBuffer <<- 1101

foreach(i=1:len) %dopar% {
	left = bins[i] - width / 2
	right = bins[i] + width / 2
	makeBins(feature, bins[i], left, right)
}