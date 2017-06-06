import pickle as pkl
import random
import pandas as pd
from bokeh.plotting import ColumnDataSource
import copy



def get_source(source_name = 'data_source.pkl'):
    DATA_DIR = '/Users/student/Dropbox/UCSF/Fischbach/HSV/data/bokeh_files/'
    # SOURCE_NAME = 'data_source.pkl'
    SOURCE_NAME = source_name
    df_counts, HHV_MAP = pkl.load(open(DATA_DIR + SOURCE_NAME, 'rb'))
    source = ColumnDataSource(df_counts)
    return source, HHV_MAP

def dump_source():
    DATA_DIR = '/Users/student/Dropbox/UCSF/Fischbach/HSV/data/bokeh_files/'
    ORTHO_DIR = '/Users/student/Dropbox/UCSF/Fischbach/HSV/data/orthomcl_files/'
    SOURCE_NAME = 'data_source.pkl'
    DB_FNAME = 'HHV_ortho_groups.db'
    DF_FNAME_GROUPS = 'HHV_ortho_groups_80.df'
    DF_FNAME_COUNTS = 'HHV_ortho_counts_80.df'
    DF_ANNOTATIONS = 'HHV_id_annotations.df'

    df_counts = pkl.load(open(ORTHO_DIR + DF_FNAME_COUNTS, 'rb'))
    df_groups = pkl.load(open(ORTHO_DIR + DF_FNAME_GROUPS, 'rb'))
    df_annotations = pkl.load(open(DATA_DIR + DF_ANNOTATIONS, 'rb'))


    max_count = df_counts.max().max()
    HHV_names = list(df_counts.columns.values)
    HHV_MAP = {}
    for i,HHV_name in enumerate(HHV_names):
        HHV_MAP[HHV_name] = i

    cluster_numbers = list(df_counts.index)
    color = []
    counts = []
    # HHV_names = []
    HHV_names_mapped = []
    cluster_numbers_mapped = []
    HHV_names_raw = []
    genes = []
    products = []
    notes = []
    x_positions = []
    y_positions = []



    for c, cluster_number in enumerate(cluster_numbers):
        for h, HHV_name in enumerate(HHV_names):
            count = df_counts.loc[cluster_number][HHV_name]
            # print(count, cluster_number, HHV_name)
            counts.append(count)
            color.append("#%02x%02x%02x" % (int(255), int(255 - (count / max_count) * 255.0), int(255 - (count / max_count) * 255.0)))
            cluster_numbers_mapped.append(cluster_number)
            HHV_names_mapped.append(HHV_MAP[HHV_name])
            HHV_names_raw.append(HHV_name)
            gene_list = df_groups.loc[cluster_number][HHV_name]
            genes.append(gene_list)
            x_positions.append(100 * random.randrange(100))
            y_positions.append(100 * random.randrange(100))


            product = ''
            note = ''
            for gene in gene_list.split("|"):
                df_match = df_annotations.loc[df_annotations['id'] == gene]
                if not df_match.empty:
                    product = product + df_match['product'].values[0] + '|'
                    note = note + df_match['note'].values[0] + '|'

            product = product[0 : len(product) - 1]
            note = note[0 : len(note) - 1]
            products.append(product)
            notes.append(note)


    cluster_numbers_mapped = [20 * i for i in cluster_numbers_mapped]

    data_dict = {'cluster_number' : cluster_numbers_mapped, 'color' : color, 
                 'HHV_name_mapped' : HHV_names_mapped, 'gene' : genes,
                 'product' : products, 'note' : notes, 
                 'HHV_name_raw' : HHV_names_raw, 'x' : x_positions,
                 'y' : y_positions}

    df_counts = pd.DataFrame(data_dict, columns = list(data_dict.keys()))
    # source = ColumnDataSource(df_counts)

    pkl.dump((df_counts, HHV_MAP), open(DATA_DIR + SOURCE_NAME, 'wb'))

    return data_dict

def get_groups_data():
    data_dict = dump_source()
    data_dict_2 = copy.deepcopy(data_dict)
    data_dict_3 = copy.deepcopy(data_dict)

    colors_2 = []
    colors_3 = []

    for i in range(len(data_dict['color'])):
        colors_2.append("#%02x%02x%02x" % (int(255), int(255), int(255)))
        colors_3.append("#%02x%02x%02x" % (int(0), int(0), int(0)))

    data_dict_2['color'] = colors_2
    data_dict_3['color'] = colors_3


    df_counts_1 = pd.DataFrame(data_dict, columns = list(data_dict.keys()))
    df_counts_2 = pd.DataFrame(data_dict_2, columns = list(data_dict_2.keys()))
    df_counts_3 = pd.DataFrame(data_dict_3, columns = list(data_dict_3.keys()))

    source_1 = ColumnDataSource(df_counts_1)
    source_2 = ColumnDataSource(df_counts_2)
    source_3 = ColumnDataSource(df_counts_3)

    # sources = [source_1, source_2, source_3]
    groups_data = {'groups_1' : [source_1], 'groups_2': [source_2], 'groups_3': [source_3]}

    return groups_data



