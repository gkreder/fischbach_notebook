
######################################################################
#                                                                 
#	Convert genbank files to fasta files			         
#	Gabe Reder
#	May 2017
#                                                        
######################################################################
#!/bin/bash
# Takes all gb files in ./HHV_gb and converts to fasta files by calling
# python script in folder genbank_to_fasta_v1.1
# Can change the tags that get included in the resulting fasta files
# by modifying the -q flag

cd ./data/HHV_gb
array=(*)
cd ..
top_dir="$PWD"

for i in "${array[@]}"
do
	in=$top_dir/HHV_gb/$i
	out=$top_dir/HHV_fasta/${i/.gb/.fasta}
	python ./genbank_to_fasta_v1.1/genbank_to_fasta.py -i $in  -o $out -d 'pipe' -q "locus_tag,gene,protein_id,product,location"
done
