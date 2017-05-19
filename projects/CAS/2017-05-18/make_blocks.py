import os
import pandas as pd
import pickle as pkl
import sys
from pfam_object import pfam

DATA_DIR = './done/'
COLUMNS = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
WINDOW_HALF_SIZE = 25

# pfam_term = 'Cas9_REC'

if len(sys.argv) < 2:
	print('Error: No PFAM argument term given')
	sys.exit(1)

pfam_term = sys.argv[1]


def get_surrounding_pfams(file, term_location_start, term_location_end):
	with open(file, 'r') as f:
		lines = f.readlines()

		for line in lines:
			line_elements = line.strip().split('\t')




top_dir = os.getcwd()
data_dir = top_dir + '/data/' + pfam_term
pfam_dir = data_dir + '/pfam/'

if not os.path.exists(data_dir):
	print("Error: Couldn't find PFAM directory for given term")
	print("\t Searched for directory: " + data_dir)
	sys.exit(1)

os.chdir(data_dir)
pfam_dicts = []

if not os.path.exists(pfam_dir):
	print("Error: Couldn't find PFAM directory for given term")
	print("\t Searched for directory: " + pfam_dir)
	sys.exit(1)

os.chdir(pfam_dir)

df = pd.DataFrame()
for file in os.listdir():
	if '.pfd' in file:
		pfam_dicts.append(file)



term_lines = []
before_blocks = []
after_blocks = []
term_blocks = []
for pfd_file in pfam_dicts:
	pfd = open(pfd_file, 'r')
	lines = pfd.readlines()
	for line_number, line in enumerate(lines):
		line_elements = line.strip().split('\t')		
		pfam_id = line_elements[-2]
		if pfam_term in pfam_id:
			term_lines.append(line_elements)
			before_block = []
			after_block = []
			for i in range(max(line_number - WINDOW_HALF_SIZE, 0), line_number):
				pfam_temp = pfam(lines[i].strip().split('\t'))
				before_block.append(pfam_temp)
			for i in range(line_number + 1, min(line_number + WINDOW_HALF_SIZE + 1, len(lines) - 1)):
				pfam_temp = pfam(lines[i].strip().split('\t'))
				after_block.append(pfam_temp)
			pfam_temp = pfam(lines[line_number].strip().split('\t'))
			term_blocks.append(pfam_temp)
			before_blocks.append(before_block)
			after_blocks.append(after_block)


print(after_blocks[0][0].sequence)
os.chdir(data_dir)
file_name = pfam_term + '_blocks.out'
pkl.dump((term_blocks, before_blocks, after_blocks), open(file_name, 'wb'))

