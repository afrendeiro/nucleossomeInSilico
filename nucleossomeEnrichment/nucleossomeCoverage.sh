### Make bins around features of interest (TSSs, TESs...)

for FEATURE in TSS TES
do
    OUTDIR=~/data/oikopleura/windows/${FEATURE}
    mkdir -p $OUTDIR
    # Call script "Rscript makeBins.R" to see arguments
    Rscript ../makeBins.R $FEATURE 500 200 20 $OUTDIR
done

### Compute Nucleossome enrichment in windows around features of interest

NUCLEOSSOME=~/data/oikopleura/nucleossome/bound.bed

for FEATURE in TSS TES
do
    OUTDIR=/home/s3/afr/data/oikopleura/nucleossome/$FEATURE
    mkdir -p $OUTDIR

    cd ~/data/oikopleura/windows/$FEATURE

    # Calculate coverage in those windows
    for WINDOW in `ls `
    do
        WINDOWNAME=`basename ${WINDOW/.bed/}`
        echo "doing window $WINDOWNAME on feature $FEATURE"
        intersectBed -wa -c -a $WINDOW -b $NUCLEOSSOME > $OUTDIR/$WINDOWNAME.bed
    done
done

### Concatenate window values
for FEATURE in TSS TES
do
    Rscript concatenateCoverage.R ~/data/oikopleura/nucleossome/$FEATURE ~/data/oikopleura/nucleossome/nucleossome.${FEATURE}_coverage.Rdata
done