
######################################################################
#                                                                 
#	Runs all-vs-all blast
#	Gabe Reder
#	May 2017
#                                                        
######################################################################
#!/bin/bash
# Makes blastdb by calling makeblast db. Then runs all vs all blast 
# with the -outfmt 6 flag to get tabular output. The rest of the 
# blastp params are default 

cd data
top_dir="$PWD"
# array=(*)
# cd ..
# top_dir="$PWD"
# IFS='_'
cd orthomcl_files


echo "----------------------------------------"
echo "Making BLAST db"
echo "----------------------------------------"
makeblastdb -in goodProteins.fasta -dbtype prot -out ortho_blast_db


echo "----------------------------------------"
echo "Running BLAST all-vs-all"
date
echo "----------------------------------------"
blastp -db ortho_blast_db -query goodProteins.fasta -outfmt 6 -out all-vs-all.out
