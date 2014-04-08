# read nuc pos info
load("~/data/oikopleura/nucleossome/nucleossome.TSS_coverage.Rdata")
NucleossomesonTSS = concate
load("~/data/oikopleura/nucleossome/nucleossome.TES_coverage.Rdata")
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

# read GC content info
load("~/data/oikopleura/gcContent/gcContent.TSS_coverage.mean.Rdata")
gcContentTSSs = concateMean
load("~/data/oikopleura/gcContent/gcContent.TES_coverage.mean.Rdata")
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
SL<-read.delim("~/data/oikopleura/trans-splicing/SL_AG_genes+IDs.bed",header=F)
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

load("~/data/oikopleura/annotation/operon_gene_position.split.Rdata")
# this is called splitOperons

# grab all genes in operons
OP_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% splitOperons$gene), ]
OP_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% splitOperons$gene), ]

# grab genes that are in position 1 of operon
OP1 = splitOperons[splitOperons$OperonPosition == 1, ]
OP1_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% OP1$gene), ]
OP1_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% OP1$gene), ]

# grab genes that are in position 2 or MORE of operon
OP2m = splitOperons[splitOperons$OperonPosition >= 2, ]
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