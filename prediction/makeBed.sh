# Bound or not (ON/OFF)
cd ~/data/oikopleura/nucleossome/scaffolds

for FILE in `find -name '*.nucl'`
do
SCAFFOLD=`basename $FILE .nucl`
# filter bases not nucleossome bound
tail -n +2 $FILE | awk -v scaffold=$SCAFFOLD '$5 == 1 {print scaffold, $1, ($1+1)}' OFS='\t' > ../bound/$SCAFFOLD.bound.bed
done

## join files and sort
cd ~/data/oikopleura/nucleossome/bound

for FILE in `find -name '*.bed'`
do
cat $FILE >> ../bound.bed
done

cd ~/data/oikopleura/nucleossome

# merge overlapping in one entry (nucleossome position)
bedtools merge -i bound.bed | bedtools sort > boundtmp
mv boundtmp bound.bed
