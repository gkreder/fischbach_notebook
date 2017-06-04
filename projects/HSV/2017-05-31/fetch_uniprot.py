import urllib
from bioservices import UniProt
from Bio import SeqIO
import xml.etree.ElementTree as ET
import pandas as pd
import time
# url = 'http://www.uniprot.org/uniprot/'
# params = {
# 'query' : 'YP_401665.1',
# 'sort' : 'score'
# }


GB_DIR = '/Users/student/Dropbox/UCSF/Fischbach/data/HSV/GB/'
DATA_DIR = '/Users/student/Dropbox/UCSF/Fischbach/data/HSV/'
OUT_FILE = 'HHV_id_annotations.df'


def parse_results(res):
	lines = res.split('\n')	
	headers = lines[0]
	headers = headers.split('\t')
	# headers = [x for x in headers if x != '']
	results_parsed = []
	result_lines = lines[1 : ]
	keys = []


	for header in headers:
		keys.append(header)

	for line_index, line in enumerate(result_lines):
		current_result = {}
		entries = line.split('\t')
		for entry_index, entry in enumerate(entries):
			key = keys[entry_index]
			current_result[key] = entry
		current_result['rank'] = line_index
		results_parsed.append(current_result)
	
	return(results_parsed)


def create_record_dict(id_string, u):
	res = u.search(id_string)
	results = parse_results(res)

	if len(results) > 0:
		best_result = results[0]['Entry']
		u_entry = u.retrieve(best_result, frmt = 'xml')
		out_file = open('temp.xml', 'w')
		out_file.write(str(u_entry))
		out_file.close()

		# in_file = open('temp.xml', 'r')
		# xml_lines = in_file.read()
		record = SeqIO.read('temp.xml', 'uniprot-xml')
		record_dict = {}
		db_keys = []
		for key in dir(record):
			if key[0] != '_':
				attr = getattr(record, key)
				# print('------------------------------------------------------------------------')
				# print(key)
				# print('------------------------------------------------------------------------')
				# print(attr)
				# print('\n\n')

				if key == 'annotations':
					annotations = attr
					for annotation in annotations:
						annotation_attr = annotations[annotation]
						if type(annotation_attr) is str:
							annotation_attr = annotation_attr.strip()
						record_dict[annotation] = annotation_attr
				else:
					if type(attr) is str:
						attr = attr.strip()
					record_dict[key] = attr

	else:
		print('----------------------------')
		print("No RESULTS FOR: " + id_string)
		print('----------------------------')
		record_dict = {}


	return(record_dict)



# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# df = pd.DataFrame()
# u = UniProt(verbose = False)
# id_string = 'YP_401665.1'
# record_dict = create_record_dict(id_string)

# # record_df = pd.DataFrame(record_dict)
# # df = df.append(record_df)
# # print(list(df))
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------



