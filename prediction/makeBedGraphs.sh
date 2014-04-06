

# Start probability

for FILE in `find -name '*.nucl'`
do
SCAFFOLD=`basename $FILE .nucl`
tail -n +2 $FILE | awk -v scaffold=$SCAFFOLD '{print scaffold, ($1-1), ($1), $3}' OFS='\t' > ../startProb/$SCAFFOLD.startProb.bedgraph
done

## join files and sort
cd ../startProb/

for FILE in `find -name '*.bedgraph'`
do
cat $FILE >> ../startProb.bedgraph 
done

bedtools sort -i ../startProb.bedgraph > start

echo "track type=bedGraph name='P-NucStart' description='Probability of Nucleossome Start' visibility=full color=200,100,0 altColor=0,100,200 priority=20" > ../startProb.bedgraph
cat start >> ../startProb.bedgraph

bedGraphToBigWig ../startProb.bedgraph ~/data/oikopleura/assembly/Oikopleura_chrSizes.txt ../startProb.bw

# Occupancy probability

for FILE in `find -name '*.nucl'`
do
SCAFFOLD=`basename $FILE .nucl`
tail -n +2 $FILE | awk -v scaffold=$SCAFFOLD '{print scaffold, $1, ($1+1), $4}' OFS='\t' > ../occupProb/$SCAFFOLD.occupProb.bedgraph
done

## join files and sort
cd ../occupProb/

for FILE in `find -name '*.bedgraph'`
do
cat $FILE >> ../occupProb.bedgraph 
done

bedtools sort -i ../occupProb.bedgraph > occup

echo "track type=bedGraph name='P-Occupancy' description='Probability of Nucleossome Occupancy' visibility=full color=200,100,0 altColor=0,100,200 priority=20" > ../occupProb.bedgraph
cat occup >> ../occupProb.bedgraph

bedGraphToBigWig ../occupProb.bedgraph ~/data/oikopleura/assembly/Oikopleura_chrSizes.txt ../occupProb.bw

# Affinity

for FILE in `find -name '*.nucl'`
do
SCAFFOLD=`basename $FILE .nucl`
# filter bases without affinity value
tail -n +2 $FILE | awk -v scaffold=$SCAFFOLD '$6 != "NA" {print scaffold, $1, ($1+1), $6}' OFS='\t' > ../aff/$SCAFFOLD.aff.bedgraph
done

## join files and sort
cd ../aff/

for FILE in `find -name '*.bedgraph'`
do
cat $FILE >> ../aff.bedgraph 
done

bedtools sort -i ../aff.bedgraph > aff

echo "track type=bedGraph name='Affinity' description='Measure of Nucleossome Affinity' visibility=full color=200,100,0 altColor=0,100,200 priority=20" > ../aff.bedgraph
cat aff >> ../aff.bedgraph

bedGraphToBigWig ../aff.bedgraph ~/data/oikopleura/assembly/Oikopleura_chrSizes.txt ../aff.bw