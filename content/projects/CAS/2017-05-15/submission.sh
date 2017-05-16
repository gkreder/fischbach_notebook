#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -o /netapp/home/gkreder/Fischbach/std_out
#$ -e /netapp/home/gkreder/Fischbach/std_err
#$ -l arch=linux-x64
#$ -l netapp=8G,scratch=8G
#$ -l h_rt=48:00:00
#$ -t 1-229



cd all_gb
array=(*)
index=$(($SGE_TASK_ID-1))
cd ..
python cas9_pfam_cluster.py ${array[$index]}
