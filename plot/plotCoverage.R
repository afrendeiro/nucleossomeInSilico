nucleossomeDir = "/data/oikopleura/nucleossome/"

# read nuc pos info
load(paste(nucleossomeDir, "nucleossome.TSS_coverage.Rdata", sep="")
NucleossomesonTSS = concate
load(paste(nucleossomeDir, "nucleossome.TES_coverage.Rdata", sep="")
NucleossomesonTES = concate

# Build dataframes with 'factors'
nucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(NucleossomesonTSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "All"
)
nucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(NucleossomesonTES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "All"
)

gcDir = "/data/oikopleura/nucleossome/"
# read GC content info
load(paste(gcDir, "gcContent.TSS_coverage.mean.Rdata"),sep="")
gcContentTSSs = concateMean
load(paste(gcDir, "gcContent.TES_coverage.mean.Rdata"),sep="")
gcContentTESs = concateMean

gcTSS <- data.frame(
	x = unlist(as.numeric(colnames(gcContentTSSs))),
	y = unlist(gcContentTSSs),
	type = "GC",
	region = "TSS",
	sub = "All"
)
gcTES <- data.frame(
	x = unlist(as.numeric(colnames(gcContentTESs))),
	y = unlist(gcContentTESs),
	type = "GC",
	region = "TES",
	sub = "All"
)

#### TRANSPLICED GENES
# read SL gene list
SLgenes <- "~/data/oikopleura/trans-splicing/"
SL<-read.delim(paste(SLgenes, "SL_AG_genes+IDs.bed",sep=""), header=F)
SLgenes<-unique(SL$V5)
SL_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% SLgenes), ]
SL_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% SLgenes), ]


SLnucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(SL_nuc_TSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "SL"
)
SLnucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(SL_nuc_TES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "SL"
)
SLgcTSS <- data.frame(
	x = unlist(as.numeric(colnames(gcContentTSSs))),
	y = unlist(gcContentTSSs),
	type = "GC",
	region = "TSS",
	sub = "SL"
)
SLgcTES <- data.frame(
	x = unlist(as.numeric(colnames(gcContentTESs))),
	y = unlist(gcContentTESs),
	type = "GC",
	region = "TES",
	sub = "SL"
)
#### Operons
op <- "~/data/oikopleura/annotation/"
load(paste(op, "operon_gene_position.split.Rdata"), sep="")
# this is called splitOperons

# grab all genes in operons
OP_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% splitOperons$gene), ]
OP_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% splitOperons$gene), ]

# grab genes that are in position 1 of operon
OP1 = splitOperons[splitOperons$OperonPosition == 1, ]
OP1_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% OP1$gene), ]
OP1_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% OP1$gene), ]

# grab genes that are in position 2 or MORE of operon
OP2m = splitOperons[splitOperons$OperonPosition >= 2 & splitOperons$OperonPosition < splitOperons$OperonLength, ]
OP2m_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% OP2m$gene), ]
OP2m_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% OP2m$gene), ]

# grab genes that are in the LAST position of a operon
OPL = splitOperons[splitOperons$OperonPosition == splitOperons$OperonLength, ]
OPL_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% OPL$gene), ]
OPL_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% OPL$gene), ]

OPnucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(OP_nuc_TSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "Operon-All"
)
OPnucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(OP_nuc_TES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "Operon-All"
)

OP1nucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(OP1_nuc_TSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "Operon-First"
)
OP1nucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(OP1_nuc_TES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "Operon-First"
)

OP2MnucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(OP2m_nuc_TSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "Operon-Middle"
)
OP2MnucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(OP2m_nuc_TES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "Operon-Middle"
)

OPLnucleoTSS <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTSS[,-1]))),
	y = unlist(colMeans(OPL_nuc_TSS[,-1])),
	type = "Nucleossomes",
	region = "TSS",
	sub = "Operon-Last"
)
OPLnucleoTES <- data.frame(
	x = unlist(as.numeric(colnames(NucleossomesonTES[,-1]))),
	y = unlist(colMeans(OPL_nuc_TES[,-1])),
	type = "Nucleossomes",
	region = "TES",
	sub = "Operon-Last"
)

# Join
NucleossomeGC <- rbind(nucleoTSS, nucleoTES, gcTSS, gcTES, SLnucleoTSS, SLnucleoTES, SLgcTSS, SLgcTES, OPnucleoTSS, OPnucleoTES, OP1nucleoTSS, OP1nucleoTES, OP2MnucleoTSS, OP2MnucleoTES, OPLnucleoTSS, OPLnucleoTES)

# Save
save(NucleossomeGC, file="~/data/oikopleura/nucleossome/nucleossome.Rdata")

##### Plot
library(ggplot2)

# Plot with 2 factors (type & region)
pdf("~/data/oikopleura/nucleossome/plots/Nucleossomes+GC_TSS_TES_allgenes.40bpWind+overlap.pdf")
p <- ggplot(data = NucleossomeGC, mapping = aes(x = x, y = y)) +
	geom_line() +
	facet_grid(type ~ region, scale = "free") +
	xlab("distance to transcription start / end (bp)") +
	ylab("nucleossome occupancy / gc content") +
	theme_bw()
	p
dev.off()

library(ggplot2)
# Plot with 2 factors (region, SL/notSL/Operon...)
plotNuc = NucleossomeGC[NucleossomeGC$type == 'Nucleossomes',]
pdf("~/data/oikopleura/nucleossome/plots/Nucleossomes_TSS_TES_allgenes+SLgenes+Operons.40bpWind+overlap.pdf")
p <- ggplot(data = plotNuc, mapping = aes(x = x, y = y)) +
	geom_line() +
	facet_wrap(sub ~ region, ncol = 2) +
	xlab("distance to transcription start / end (bp)") +
	ylab("nucleossome occupancy / gc content") +
	theme_bw()
	p
dev.off()

library(ggplot2)
# Plot with 3 factors (type & region, SL/notSL/Operon...)
pdf("~/data/oikopleura/nucleossome/plots/Nucleossomes+GC_TSS_TES_allgenes+SLgenes+Operons.40bpWind+overlap.pdf")
p <- ggplot(data = NucleossomeGC, mapping = aes(x = x, y = y)) +
	geom_line() +
	facet_grid(type ~ region ~ sub, scale = "free") +
	xlab("distance to transcription start / end (bp)") +
	ylab("nucleossome occupancy / gc content") +
	theme_bw()
	p
dev.off()