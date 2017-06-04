
######################################################################
#                                                                 
#	Parses Blast results
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
echo "Parsing BLAST results"
echo "----------------------------------------"
orthomclBlastParser all-vs-all.out ../HHV_compliant >> similarSequences.txt
