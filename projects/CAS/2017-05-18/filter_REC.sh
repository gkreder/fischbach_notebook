#!/bin/bash


cd all_pfam
array=(*.pfd)

for i in "${array[@]}"
do
	if !($(grep "-q" "Cas9_REC" "$i"))
		then
			rm "${i/.pfd/.pfs}"
			rm $i
		fi
done