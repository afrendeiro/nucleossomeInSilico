# read nuc pos info
load("~/data/oikopleura/nucleossome/nucleossome.TSS_coverage.Rdata")
NucleossomesonTSS = concate
load("~/data/oikopleura/nucleossome/nucleossome.TES_coverage.Rdata")
NucleossomesonTES = concate

# read GC content info
load("~/data/oikopleura/gcContent/gcContent.TSS_coverage.mean.Rdata")
gcContentTSSs = concateMean
load("~/data/oikopleura/gcContent/gcContent.TES_coverage.mean.Rdata")
gcContentTESs = concateMean

#### TRANSPLICED GENES
# read SL gene list
SL<-read.delim("~/data/oikopleura/trans-splicing/SL_AG_genes+IDs.bed",header=F)
SLgenes<-unique(SL$V5)
SL_nuc_TSS <- NucleossomesonTSS[which(NucleossomesonTSS$gene %in% SLgenes), ]
SL_nuc_TES <- NucleossomesonTES[which(NucleossomesonTES$gene %in% SLgenes), ]


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
# Join
NucleossomeGC <- rbind(nucleoTSS, nucleoTES, gcTSS, gcTES, SLnucleoTSS, SLnucleoTES, SLgcTSS, SLgcTES)

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

# Plot with 3 factors (type & region, SL/notSL)
pdf("~/data/oikopleura/nucleossome/plots/Nucleossomes+GC_TSS_TES_allgenes+SLgenes.40bpWind+overlap.pdf")
p <- ggplot(data = NucleossomeGC, mapping = aes(x = x, y = y)) +
	geom_line() +
	facet_grid(type ~ region ~ sub, scale = "free") +
	xlab("distance to transcription start / end (bp)") +
	ylab("nucleossome occupancy / gc content") +
	theme_bw()
	p
dev.off()