import csv
from Bio import SeqIO
 
fastagenome = "~/data/oikopleura/assembly/Oikopleura_reference_masked_v3.0.fa.msk"
outpath = "~/data/oikopleura/assembly/scaffolds/"

for seq_record in SeqIO.parse(fastagenome, "fasta"):
	filename = outpath + str(seq_record.id) + ".fasta"
	output_handle = open(filename, "w")
	SeqIO.write(seq_record, output_handle, "fasta")
	output_handle.close()
