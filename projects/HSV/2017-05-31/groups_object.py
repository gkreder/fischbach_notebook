from fetch_uniprot import create_record_dict
from bioservices import UniProt

class cluster(object):
	"""
		Attributes:
			text
			species
			genes
			species_genes
	"""
	def __init__(self, line):
		self.text = line
		line = line.split(':')[1]
		elements = [x for x in line.split(' ') if x != '']
		species_list = []
		genes_list = []
		record_dicts = []
		species_genes_dict = {}
		for element in elements:
			species = element.split('|')[0]
			gene = element.split('|')[1]
			genes_list.append(gene)
			if species not in species_list:
				species_list.append(species)
				species_genes_dict[species] = [gene]
			else:
				species_genes_dict[species].append(gene)


		self.species = species_list
		self.genes = genes_list
		self.species_genes = species_genes_dict
		self.gene_records = []

	def add_gene_record(self, gene_record):
		self.gene_records.append(gene_record)

class data_set(object):
	"""
		Attributes:
			i_val: inflation value used to generate the data_set
			clusters: 
			

	"""
	def __init__(self, clusters, i_val):
		self.i_val = i_val
		self.clusters = clusters

	def add_gene_records(self, gene_records):
		for cluster in clusters:
			for gene in cluster.genes:
				record = gene_records[gene]
				cluster.add_gene_record(record)

	def get_species_clusters(self, species_list):
		"""
		For 
		"""

