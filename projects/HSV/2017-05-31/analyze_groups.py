######################################################################
#                                                            
# HHV Groups Analysis
# Gabe Reder
# May 2017
#                                       
######################################################################
import os
import sys
from groups_object import cluster, data_set
from bioservices import UniProt
import pickle as pkl
from fetch_uniprot import create_record_dict

os.chdir('../data/')
data_dir = os.getcwd()
orthomcl_files_dir = data_dir + '/orthomcl_files/'
groups_dir = orthomcl_files_dir + 'groups/'


def grab_data():
	# Scrape groups data and store in data (list)
	if not os.path.exists(orthomcl_files_dir + 'all_groups_data.pkl'):
		print('Creating data from files and dumping...')
		groups_files = os.listdir(groups_dir)
		data = []
		for file in groups_files:
			clusters = []
			inflation_val = file.replace('group_', '').replace('.txt', '').replace('_', '.')
			with open(groups_dir + file, 'r') as f:
				lines = f.read().splitlines()
			for line in lines:
				cluster_ob = cluster(line)
				clusters.append(cluster_ob)

			current_data_set = data_set(clusters, inflation_val)
			data.append(current_data_set)
		pkl.dump(data, open(orthomcl_files_dir + 'all_groups_data_no_records.pkl', 'wb'))

		gene_records = {}
		u = UniProt(verbose = False)
		for ds_number, ds_obj in enumerate(data):
			print('Dataset ' + str(ds_number) + ' of ' + str(len(data)))
			for cluster in ds_obj.clusters:
				for gene in cluster.genes:
					if gene not in gene_records:
						record = create_record_dict(gene, u)
						gene_records[gene] = record
					record = gene_records[gene]
					cluster.add_gene_record(record)

		pkl.dump(data, open(orthomcl_files_dir + 'all_groups_data.pkl', 'wb'))

	else:
		print('Grabbing from .pkl file...')
		data = pkl.load(open(orthomcl_files_dir + 'all_groups_data.pkl', 'rb'))

	return data




data = grab_data()



# # test_gene = data[0].clusters[0].genes[3]
# test_gene = 'YP_401705.1'
# u = UniProt(verbose = False)
# record_dict = create_record_dict(test_gene, u)
# print(record_dict)

		