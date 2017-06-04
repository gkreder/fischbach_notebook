array=(1,2,3)
# echo $array

rm -r ./data/orthomcl_files/mcl_logs
mkdir ./data/orthomcl_files/mcl_logs
rm -r ./data/orthomcl_files/groups
mkdir ./data/orthomcl_files/groups

for i in `seq 1 .1 10`
do
	echo "----------------------------------------"
	echo "$i"
	echo "----------------------------------------"
	I=$i
	PREF=${i/./_}
	mcl ./data/orthomcl_files/mclInput --abc -I $I -o ./data/orthomcl_files/mclOutput &> ./data/orthomcl_files/mcl_logs/mcl_$PREF.log
	# "mcl ./data/orthomcl_files/mclInput --abc -I $I -o ./data/orthomcl_files/mclOutput"
	orthomclMclToGroups INFL_$PREF 1000 < ./data/orthomcl_files/mclOutput > ./data/orthomcl_files/groups/group_$PREF.txt

done


# mcl ./data/orthomcl_files/mclInput --abc -I $I -o ./data/orthomcl_files/mclOutput > ./data/orthomcl_files/mcl_logs/mcl_$PREF.log
# mcl ./data/orthomcl_files/mclInput --abc -I 1.5 -o ./data/orthomcl_files/mclOutput
# orthomclMclToGroups my_prefix 1000 < ./data/orthomcl_files/mclOutput > ./data/orthomcl_files/groups.txt