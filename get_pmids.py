#Copyright (C) Systems & Technology Research
#Use of this software is subject to the restrictions in license.txt
#Distribution A: Approved for Public Release, Distribution Unlimited
import numpy as np
import pandas as pd
import json
from tqdm import tqdm

'''
Takes in a list of json file paths and converts jsons to contain pmids for each relationship
Highly recommended that if you are planning to convert many json files, only call this function once, 
because it takes extra time to load each of the csvs to dataframes at the beginning

Creates new files using {json_path}_with_pmids.json naming
'''
def get_pmids(json_path_list, knowledge_graph_folder_path):
    
    KG_path = knowledge_graph_folder_path

    print('Loading BLENDER Knowledge Graph CSVs to Pandas Dataframe')
    chemicals_df = pd.read_csv(f'{KG_path}/chemicals.csv', delimiter = "	")
    diseases_df = pd.read_csv(f'{KG_path}/diseases.csv', delimiter = "	")
    genes_df = pd.read_csv(f'{KG_path}/genes.csv', delimiter = "	")

    chem_to_gene_df = pd.read_csv(f'{KG_path}/chem_gene_ixns_relation.csv', delimiter = "	")
    chem_to_disease_df = pd.read_csv(f'{KG_path}/chemicals_diseases_relation.csv', delimiter = "	")
    gene_to_disease_df = pd.read_csv(f'{KG_path}/genes_diseases_relation.csv', delimiter = "	")
    
    for json_path in tqdm(json_path_list):
        jd = None
        with open(f'{json_path}') as f:
            jd = f.read()

        sample_subgraph = json.loads(jd)

        # Load a categories array (this assues that ids are presented in numerical order)

        category_arr = list()
        for category in sample_subgraph['categories']:
            category_arr.append(category['name'])

        # Load a nodes array (same assumption as above)
        nodes_arr = list()
        for node in sample_subgraph['nodes']:
            nodes_arr.append({'name': node['name'], 'category': category_arr[node['category']]})



        print('-'*20)
        print('-'*20)
        print('Looping through edges and collecting all PMIDS')

        # Loop through links and collect desired pmids
        links_arr = list()
        for link in tqdm(sample_subgraph['links']):
            src_node = nodes_arr[link['source']]['name']
            src_category = nodes_arr[link['source']]['category']
            trg_node = nodes_arr[link['target']]['name']
            trg_category = nodes_arr[link['target']]['category']

            
            if src_category == "Chemical" and trg_category == "Disease":
            
            # Find src_id in chemicals_df
                src_id = list(chemicals_df[chemicals_df['ChemicalName'] == src_node]['ChemicalID'])[0]

            # Find target_id in diseases_df
                trg_id = list(diseases_df[diseases_df['DiseaseName'] == trg_node]['DiseaseID'])[0]

            # Search chem_to_disease_df

                # Search Chemicals that match src_id
                all_chems = chem_to_disease_df[chem_to_disease_df['ChemicalID'] == src_id]
                
                # Search Diseases that match trg_id
                pmids = list(all_chems[all_chems['DiseaseID'] == trg_id]['pmids'].dropna())
                pmids_list = [x for xs in pmids for x in xs.split('|')]
                
                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string
                

            elif src_category == "Disease" and trg_category == "Chemical":
            # Find src_id in diseases_df
                src_id = list(diseases_df[diseases_df['DiseaseName'] == src_node]['DiseaseID'])[0]

            # Find target_id in chemicals_df
                trg_id = list(chemicals_df[chemicals_df['ChemicalName'] == trg_node]['ChemicalID'])[0]

            # Search chem_to_disease_df
                # Search Chemicals that match src_id
                all_chems = chem_to_disease_df[chem_to_disease_df['ChemicalID'] == trg_id]
                
                # Search Diseases that match trg_id
                pmids = list(all_chems[all_chems['DiseaseID'] == src_id]['pmids'].dropna())
                pmids_list = [x for xs in pmids for x in xs.split('|')]
                
                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string


            elif src_category == "Gene" and trg_category == "Chemical":

            # Find src_id in genes_df
                src_id = list(genes_df[genes_df['GeneName'] == src_node]['GeneID'])[0]
                
            # Find target_id in chemicals_df
                trg_id = list(chemicals_df[chemicals_df['ChemicalName'] == trg_node]['ChemicalID'])[0]

            # Search chem_to_gene_df
                
                # Search for Chemical in ChemicalID
                all_chems = chem_to_gene_df[chem_to_gene_df['ChemicalID'] == trg_id]

                # Search for Gene in GeneID
                pmids = list(all_chems[all_chems['GeneID'] == src_id]['pmids'].dropna())
                
                pmids_list = [x for xs in pmids for x in xs.split('|')]

                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string

            elif src_category == "Chemical" and trg_category == "Gene":

            # Find src_id in chemicals_df
                src_id = list(chemicals_df[chemicals_df['ChemicalName'] == src_node]['ChemicalID'])[0]
                
            # Find target_id in genes_df
                trg_id = list(genes_df[genes_df['GeneName'] == trg_node]['GeneID'])[0]

            # Search chem_to_gene_df
                # Search for Chemical in ChemicalID
                all_chems = chem_to_gene_df[chem_to_gene_df['ChemicalID'] == src_id]

                # Search for Gene in GeneID
                pmids = list(all_chems[all_chems['GeneID'] == trg_id]['pmids'].dropna())
                
                pmids_list = [x for xs in pmids for x in xs.split('|')]
                
                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string


            elif src_category == "Disease" and trg_category == "Gene":
            # Find src_id in diseases_df
                src_id = list(diseases_df[diseases_df['DiseaseName'] == src_node]['DiseaseID'])[0]

            # Find target_id in genes_df
                trg_id = list(genes_df[genes_df['GeneName'] == trg_node]['GeneID'])[0]

            # Search gene_to_disease_df
                all_genes = gene_to_disease_df[gene_to_disease_df['GeneID'] == trg_id]
                pmids = list(all_genes[all_genes['DiseaseID'] == src_id]['pmids'].dropna())

                pmids_list = [x for xs in pmids for x in xs.split('|')]

                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string

            elif src_category == "Gene" and trg_category == "Disease":

            # Find src_id in genes_df
                src_id = list(genes_df[genes_df['GeneName'] == src_node]['GeneID'])[0]

            # Find target_id in diseases_df
                trg_id = list(diseases_df[diseases_df['DiseaseName'] == trg_node]['DiseaseID'])[0]

            # Search gene_to_disease_df
                all_genes = gene_to_disease_df[gene_to_disease_df['GeneID'] == src_id]
                pmids = list(all_genes[all_genes['DiseaseID'] == trg_id]['pmids'].dropna())
                
                pmids_list = [x for xs in pmids for x in xs.split('|')]
                
                edgetype_string = link['edgetype']
                new_edgetype_string = edgetype_string + '\nsource:'
                for pmid in pmids_list:
                    new_edgetype_string += pmid + ","
                link['edgetype'] = new_edgetype_string
    
        # Save sample_subgraph to a new json file
        # Get just the path without .json ending
        new_json_path = json_path.replace('.json', '_with_pmids.json')
        with open(f'{new_json_path}', "w") as f:
            f.write(json.dumps(sample_subgraph))


