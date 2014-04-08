library(rtracklayer)

operons <- import("~/data/oikopleura/annotation/operons_GGB.gff")

gene <- character(0)
OperonPosition <- integer(0)
OperonLength <- integer(0)

for (operon in 1:length(operons)) {
	# remove "operon " and split string using by ","
	op = unlist(strsplit(x = sub(x = as.character(operons$group[operon]), pattern = "operon ", replacement = ""), split = ","))
	for (n in 1:length(op)) {
		# extract geneID, position in operon, number of genes in operon
		gene <- c(gene, sub(x=op[n], pattern="GSOIDT", replacement="GSOIDG"))
		OperonPosition <- c(OperonPosition, n)
		OperonLength <- c(OperonLength, length(op))
	}
}
splitOperons = data.frame(gene, OperonPosition, OperonLength)

save(splitOperons, file="~/data/oikopleura/annotation/operon_gene_position.split.Rdata")