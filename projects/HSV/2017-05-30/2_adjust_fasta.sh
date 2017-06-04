
######################################################################
#                                                                 
#	Adjusts fasta files to be compliant with orthomcl		         
#	Gabe Reder
#	May 2017
#                                                        
######################################################################
#!/bin/bash
# Adjusts converted fasta files so that they are compliant with 
# orthomcl pipeline

cd ./data/HHV_fasta
array=(*)
cd ..
top_dir="$PWD"
IFS='_'


for i in "${array[@]}"
do
	read -ra split_string <<< "$i"
	HHV_NAME=$split_string
	HHV_NAME=${HHV_NAME/-/_}	
	in=$top_dir/HHV_fasta/$i
	out=$top_dir/HHV_compliant/$HHV_NAME.fasta
	/usr/local/orthomcl/bin/orthomclAdjustFasta "$HHV_NAME" "$in" "3"
	mv "$HHV_NAME.fasta" "HHV_compliant/"
done
