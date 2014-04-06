
# Make nucleossome prediction
python ./prediction/separate_scaffolds.py
Rscript ./prediction/pred_parallel.R ~/data/oikopleura/assembly/scaffolds/ ~/data/oikopleura/nucleossome/scaffolds/ 
sh ./prediction/makeBed.sh
sh ./prediction/makeBedGraphs.sh

# Compute GC content in features of interest (TSS, TES)
sh ./gcContent/gcContent.sh

# Compute GC content in features of interest (TSS, TES)
sh ./nucleossomeEnrichment/nucleossomeCoverage.sh

# Plot
R ./plot/plotCoverage.R