
######################################################################
#                                                                 
#	Filters orthomcl compliant fasta files
#	Gabe Reder
#	May 2017
#                                                        
######################################################################
#!/bin/bash
# produces a single file called goodProteins.fasta to run BLAST on. Any 
# sequences deemed poor quality are placed in poorProteins.fasta. 
# Filtering is based on length and percent stop codon (both values can
# be adjusted). Suggested value for minimum length is 10. Suggested
# value for percent stop codon is 20. 

# cd HHV_fasta
# array=(*)
# cd ..
# top_dir="$PWD"
# IFS='_'
MIN_LENGTH=10
PERCENT_STOP=20

orthomclFilterFasta "HHV_compliant" "$MIN_LENGTH" "$PERCENT_STOP" "../data/orthomcl_files/goodProteins.fasta" "../data/orthomcl_files/poorProteins.fasta"
