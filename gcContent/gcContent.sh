### Make bins around features of interest (TSSs, TESs)

for FEATURE in TSS TES
do
    OUTDIR=~/data/oikopleura/windows/${FEATURE}
    mkdir -p $OUTDIR
    # Call script "Rscript makeBins.R" to see arguments
    Rscript ../makeBins.R $FEATURE 500 200 20 $OUTDIR
done

### Compute GC content in windows around features of interest

GENOME=~/data/oikopleura/assembly/Oikopleura_reference_masked_v3.0.fa

for FEATURE in TSS TES
do
    OUTDIR=/home/s3/afr/data/oikopleura/gcContent/$FEATURE
    mkdir -p $OUTDIR

    cd /home/s3/afr/data/oikopleura/windows/$FEATURE
    for WINDOW in `ls `
    do
        WINDOWNAME=`basename ${WINDOW/.bed/}`
        echo "doing window $WINDOWNAME on feature $FEATURE"

        bedtools nuc -fi $GENOME -bed $WINDOW > $OUTDIR/$WINDOWNAME.bed
    done

    cd $OUTDIR
    for file in `ls`; do tail -n +2 $file | cut -f 1,2,3,4,6 > tmp; mv tmp $file; done
done

### Concatenate window values
for FEATURE in TSS TES
do
    Rscript concatenateCoverageMean.R ~/data/oikopleura/gcContent/$FEATURE ~/data/oikopleura/gcContent/gcContent.${FEATURE}_coverage.mean.Rdata
done
