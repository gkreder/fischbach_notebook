# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# repr() to turn list to string and eval() to get back to list
# http://sebastianraschka.com/Articles/2013_sqlite_database.html
# import sqlite3
import time
from parse_clusters import *
import sys
import os
import pickle as pkl

FNAME_DICT = {
'HHV_1' : 'HHV-1_Strain_17_NC_001806.59968.fasta',
'HHV_2' : 'HHV-2_NC_001798.75822.fasta',
'HHV_3' : 'HHV-3_VZV_NC_001348.13322.fasta',
'HHV_4' : 'HHV-4_EBV_Reference_NC_007605.34185.fasta',
'HHV_5' : 'HHV-5_HCMV_NC_006273.43120.fasta',
'HHV_6A' : 'HHV-6A_NC_001664.22003.fasta',
'HHV_6B' : 'HHV-6B_NC_000898.82020.fasta',
'HHV_7' : 'HHV-7_NC_001716.84724.fasta',
'HHV_8' : 'HHV-8_NC_009333.81534.fasta',
'CeHV_1' : 'CeHV-1_E2490_AF533768.fasta',
'MuHV_68' : 'MuHV-68_NC_001826.fasta'
}

DEFAULT_FILE = 'group_1.txt'
# PATH_TO_ORTHO = '/Users/student/Dropbox/UCSF/Fischbach/HSV/orthomcl_fischbach/'
PATH_TO_ORTHO = '/Users/student/Dropbox/UCSF/Fischbach/HSV/data/orthomcl_files/'
# PATH_TO_DATA_FOLDER = '/Users/student/Dropbox/UCSF/Fischbach/data/HSV/'
PATH_TO_DATA_FOLDER = '/Users/student/Dropbox/UCSF/Fischbach/HSV/data/orthomcl_files/groups/'
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# 
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
start_time = time.clock()

if len(sys.argv) > 1:
    file = PATH_TO_DATA_FOLDER + sys.argv[1]
else:
    file = PATH_TO_DATA_FOLDER + DEFAULT_FILE


# db_name = file.split('/')[-1].replace('.txt', '')
# db_name = 'HHV_ortho_' + db_name

df_name = file.split('/')[-1].replace('.txt', '')
df_name_groups = 'HHV_ortho_groups_1.df'
df_name_counts = 'HHV_ortho_counts_1.df'
lines = 0
lst = list()
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# 
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
def line_to_dict(line):
    d = {}
    for key in FNAME_DICT.keys():
        d[key] = []
    
    for HHV_genes in line:
        gene = HHV_genes.split('|')
        current_HHV = gene[0]
        gene = gene[1]
        d[current_HHV].append(gene)
        # print(current_HHV, gene)
    return d

# String Version
def line_to_dict_string(line):
    d = {}
    for key in FNAME_DICT.keys():
        d[key] = ''

    for HHV_genes in line:
        gene = HHV_genes.split('|')
        current_HHV = gene[0]
        gene = gene[1]
        if len(d[current_HHV]) == 0:
            d[current_HHV] = gene
        else:
            d[current_HHV] = d[current_HHV] + '|' + gene

    for key in FNAME_DICT.keys():
        d[key] = [d[key]]

    return d

def line_to_dict_counts(line):
    d = {}
    for key in FNAME_DICT.keys():
        d[key] = 0

    for HHV_genes in line:
        gene = HHV_genes.split('|')
        current_HHV = gene[0]
        gene = gene[1]
        d[current_HHV] += 1

    for key in FNAME_DICT.keys():
        d[key] = [d[key]]

    return d
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
# 
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #
df_columns = list(FNAME_DICT.keys())
df_columns = df_columns
df_groups = pd.DataFrame(columns = df_columns)
df_counts = pd.DataFrame(columns = df_columns)
with open(file, "rb") as groups:
    for i,line in enumerate(groups):
        line = line.strip()
        line = line.decode('UTF-8').split()
        line = line[1:]
        line.sort()
        d_groups = line_to_dict_string(line)
        d_counts = line_to_dict_counts(line)

        df_groups = df_groups.append(pd.DataFrame(d_groups), ignore_index = True)
        df_counts = df_counts.append(pd.DataFrame(d_counts), ignore_index = True)

        lines += 1


pkl.dump(df_groups, open(PATH_TO_ORTHO + df_name_groups, 'wb'))
pkl.dump(df_counts, open(PATH_TO_ORTHO + df_name_counts, 'wb'))
