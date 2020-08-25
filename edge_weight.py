#Copyright (C) Systems & Technology Research
#Use of this software is subject to the restrictions in license.txt
#Distribution A: Approved for Public Release, Distribution Unlimited
import json
import pandas as pd
from tqdm import tqdm
from datetime import datetime
import numpy as np
import os

from Bio import Entrez
from Bio.Entrez import efetch, read
from Bio import Medline

import networkx as nx
import h_index_finder.h_index_finder as h_index_finder

import holoviews as hv
from holoviews import opts
import matplotlib.pyplot as plt
import datashader as ds
import datashader.transfer_functions as tf
from datashader.layout import random_layout, circular_layout, forceatlas2_layout
from datashader.bundling import connect_edges, hammer_bundle

from fa2 import ForceAtlas2

import seaborn as sns


#from sqlalchemy import create_engine
hv.extension('bokeh')

from bokeh.models import HoverTool


# Global Variables
id_2_name_dict = None
id_2_category_dict = None
all_pmids = None 
json_data = None 
pmid_dict = None
rel_name = None
api_key = None


# Blotzman Function
def boltzman(x, xmid, tau):
    return 1. / (1. + np.exp(-(x-xmid)/tau))


# Collect required information from NCBI database to calculate all different ways of edge weights
'''
Following information are collected:
Author names
Authors' affiliations
Publication date of paper/pmid-pmcid
All pmids citing this paper
Grants used to publish this paper

Note: If a pmid refers to pmcid, use pmcid to collect information since it has more details (e.g., full author 
names rather than first initial and last name)
'''

def collect_NCBI():
    global all_pmids
    global pmid_dict

    if os.path.exists(f'./{rel_name}/{rel_name}_pmid_dict.json'):
        with open(f'./{rel_name}/{rel_name}_pmid_dict.json', 'r') as f:
            jd = f.read()
            temp_dict = json.loads(jd)
        pmid_dict.update(temp_dict)
        return pmid_dict
    
    for idx in tqdm(range(len(all_pmids))):
        pmid = all_pmids[idx]
        # get records for each pmid
        fetch_records_handle1 = efetch(db="pubmed", id=str(pmid), rettype="medline",retmode="text")
        # parse fetched records
        records1 = Medline.parse(fetch_records_handle1)

        # Need to iterate over records to extract information
        for record1 in records1:
            # try except check to be sure that NCBI is not returning empty result
            try:
                # let's get pmcid if exists
                id2 = record1['PMC'][3:]
                #print('PMC',id2)

                # get records for pmcid
                fetch_records_handle2 = efetch(db="pubmed", id=str(id2),
                                               rettype="medline",retmode="text")    
                # parse records for pmcid
                records2 = Medline.parse(fetch_records_handle2)

                # Need to iterate over records to extract information
                '''
                Collect following information: authors, authors' affiliations, publication date, citations, grants
                Store all these information in an dictionary (pmid_dict)
                '''
                for record2 in records2:
                    authors = record2['FAU']
                    affiliations = record2['AD']
                    pub_date = record2['DCOM']
                    citations = get_links_id(pmid)
                    grants = record2['GR']
                    pmid_dict[pmid] = {'pmcid_number':id2,'pmcid':True,'authors':authors,'affiliations':affiliations,'grants':grants,'pub_date':pub_date,'citations':citations}
            except:
                authors = record1['FAU']
                try:
                    affiliations = record1['AD']
                except:
                    affiliations = ''
                try:
                    pub_date = record1['DCOM']
                except:
                    pub_date = ''
                try:                
                    citations = get_links_id(pmid)
                except:
                    citations = ''
                try:
                    grants = record1['GR']
                except:
                    grants = ''
                pmid_dict[pmid] = {'pmcid_number':'','pmcid':False,'authors':authors,'affiliations':affiliations,'grants':grants,'pub_date':pub_date,'citations':citations}


    with open(f'./{rel_name}/{rel_name}_pmid_dict.json', 'w') as output:
        output.write(json.dumps(pmid_dict))                
    
    return pmid_dict

# Let's get citations from pmid
'''
This function gets all pmids citing a given pmid
'''
def get_links_id(pmid):
    link_list = []
    links = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pmc_refs")
    
    record = Entrez.read(links)
    records = record[0][u'LinkSetDb'][0][u'Link']
    
    for link in records:
        link_list.append(link[u'Id'])

    return link_list


# Edge Weight Based on Citations
'''
Get all citation information from pmid_dict
'''

def get_citation_edge_weights():
    global pmid_dict
    global json_data

    all_citations = []
    # let's track the max # of citations
    max_length_citations = 0
    for key in pmid_dict.keys():
        citations = pmid_dict[key]['citations']
        if len(citations) > 1: # double check for length of number of citations
            if len(citations) > max_length_citations:
                # update length of max. # of citations
                max_length_citations = len(citations)

            # collect all citations in a lit
            for citation in citations:
                all_citations.append(int(citation))
        elif len(citations) == 1:
            all_citations.append(int(citations[0]))

    citation_weights = dict()
    
    for idx in range(len(json_data['links'])):
        cnt = 0 #TODO, maybe this should be 1? No PMIDS = 0 edgeweight right now

        # get all pmids from original json file (edges have more than one pmids)
        if json_data['links'][idx]['edgetype'].split(':')[1] != '':
            pmids = json_data['links'][idx]['edgetype'].split(':')[1].split(',')

            for pmid in pmids:
                if pmid == '':
                    continue
                cnt += len(pmid_dict[pmid]['citations'])

        # calculate edge weight based on # of citations normalized by max # of citations
        w = cnt/(max_length_citations+1)
        
        t_str = json_data['nodes'][json_data['links'][idx]['source']]['name']+'-->'+json_data['nodes'][json_data['links'][idx]['target']]['name']
        citation_weights[t_str] = w
    
    return citation_weights

# Edge Weight Based on Paper Publication Date

def get_publication_edge_weights():
    global pmid_dict
    global json_data

    # Parse publication date as datetime (it should be added into the loop where we create pmid_dict to hold all information)
    pmid_dict_year = dict()
    for pmid in pmid_dict.keys():
        if pmid_dict[pmid]['pub_date'] != '':
            pmid_dict_year[pmid] = datetime.strptime(pmid_dict[pmid]['pub_date'],'%Y%m%d').year    
        else:
            pmid_dict_year[pmid] = -np.infty

    # Calculate publication edge weight
    edge_weights_from_pub_year = dict()
    for idx in range(len(json_data['links'])):

        # let's find the oldest publication in an edge for all given edge-pmids 
        oldest_pub_year = -np.infty

        if json_data['links'][idx]['edgetype'].split(':')[1] != '':
            pmids = json_data['links'][idx]['edgetype'].split(':')[1].split(',')
            for pmid in pmids:
                if pmid == '':
                    continue
                if pmid_dict_year[pmid] > oldest_pub_year:
                    oldest_pub_year = pmid_dict_year[pmid]

            # calculate edge weight assume that mid and tau are 10 and 2, respectively 
            w = (1-boltzman(2020-oldest_pub_year, 10, 3))/(1-boltzman(0, 10, 3))
        else:
            w = 0 # there is no pmids so this should be edge weight of 0
        t_str = json_data['nodes'][json_data['links'][idx]['source']]['name']+'-->'+json_data['nodes'][json_data['links'][idx]['target']]['name']
        edge_weights_from_pub_year[t_str] = w

    return edge_weights_from_pub_year

# Let create an edge weight based on fusing publication year and # of citations edge weights (equally weight)
'''
w = 0.5 * weight from pub_year + 0.5 * weight from # of citations

it might be good idea to provide user adjustable weights
'''

def get_publication_citation_weights():
    global pmid_dict
    global json_data

    # Parse publication date as datetime (it should be added into the loop where we create pmid_dict to hold all information)
    pmid_dict_year = dict()
    for pmid in pmid_dict.keys():
        if pmid_dict[pmid]['pub_date'] != '':
            pmid_dict_year[pmid] = datetime.strptime(pmid_dict[pmid]['pub_date'],'%Y%m%d').year    
        else:
            pmid_dict_year[pmid] = -np.infty
    

    edge_weights_from_pub_year_and_citation = dict()
    for idx in range(len(json_data['links'])):


        if json_data['links'][idx]['edgetype'].split(':')[1] != '':
            # find the oldest publication for a given edge
            oldest_pub_year = -np.infty


            pmids = json_data['links'][idx]['edgetype'].split(':')[1].split(',')
            for pmid in pmids:
                if pmid == '':
                    continue
                if pmid_dict_year[pmid] > oldest_pub_year:
                    oldest_pub_year = pmid_dict_year[pmid]

            # find # of citations for a given edge
            cnt = 1


            pmids = json_data['links'][idx]['edgetype'].split(':')[1].split(',')
            for pmid in pmids:
                if pmid == '':
                    continue
                cnt += len(pmid_dict[pmid]['citations'])

            # calculte edge weight
            w = (1-boltzman(2020-oldest_pub_year, 10, 3))/(1-boltzman(0, 10, 3))*0.5+(1-boltzman(cnt, 15, 3))*0.5
        else:
            w = 0 # no pmids = edge of 0
        t_str = json_data['nodes'][json_data['links'][idx]['source']]['name']+'-->'+json_data['nodes'][json_data['links'][idx]['target']]['name']
        edge_weights_from_pub_year_and_citation[t_str] = w

    return edge_weights_from_pub_year_and_citation


def get_h_index_weights(email):
    global json_data
    global pmid_dict
    global id_2_category_dict
    global id_2_name_dict
    global api_key

    # Collect H-Index for all authors for a given subgraph

    print('Loading all Authors into List...')

    author_2_h_index = dict()
    authors = list()
    for key in tqdm(pmid_dict.keys()):
        for item in pmid_dict[key]['authors']:
            authors.append(item)


    print('Collecting H-Index for All Authors using h_index_finder.py')

    temp_author_dict = h_index_finder.find_h_index(email, authors, api_key)
    author_2_h_index.update(temp_author_dict)

    # Let's create an edge weight based on h-index
    '''
    The idea is that the highest h-index will be used for edge-weight calculation.  
    '''
    G_h_index = nx.Graph()
    source = []
    target = []
    edge_weights_from_h_index = dict()
    for idx in range(len(json_data['links'])):
        # get source(s) and target(t) nodes
        s = json_data['links'][idx]['source']
        t = json_data['links'][idx]['target']
        
        # create an edge
        e = (s,t)
        
        # add nodes into the graph
        G_h_index.add_node(s,name = id_2_name_dict[s], 
                        categories=id_2_category_dict[json_data['links'][idx]['source']])
        G_h_index.add_node(t,name=id_2_name_dict[t], 
                        categories=id_2_category_dict[json_data['links'][idx]['target']])
        
        max_h_index = -1
        # find the highest h-index for a given edge
        if json_data['links'][idx]['edgetype'].split(':')[1] != '':
            pmids = json_data['links'][idx]['edgetype'].split(':')[1].split(',')
            for pmid in pmids:
                if pmid == '':
                    continue
                authors = pmid_dict[pmid]['authors']
                max_h_index = -1
                for author in authors:
                    if author_2_h_index[author] > max_h_index:
                        max_h_index = author_2_h_index[author]
                        
            # calculte edge weight by mapping maximum h-index to Boltzmann function
            w = boltzman(max_h_index,20,2)
        else:
            w = 0 # no pmids = edgeweight of 0
        
        t_str = json_data['nodes'][json_data['links'][idx]['source']]['name']+'-->'+json_data['nodes'][json_data['links'][idx]['target']]['name']
        edge_weights_from_h_index[t_str] = w

        # create an edge
        G_h_index.add_edge(*e,value=w)
        source.append(s)
        target.append(t)

    # edge_weights_from_h_index = pd.DataFrame({'source':source,'target':target})
    G_h_index.number_of_nodes(),G_h_index.number_of_edges()

    return edge_weights_from_h_index

def get_original_edge_weights():
    global json_data
    global id_2_category_dict
    global id_2_name_dict

    # Get Original Edge weights
    G = nx.Graph()

    # extract source nodes, target nodes, and edges from json data
    source = []
    target = []
    edge_weights_from_original = dict()
    for idx in range(len(json_data['links'])):
        s = json_data['links'][idx]['source']
        t = json_data['links'][idx]['target']
        e = (s,t)
        G.add_node(s,name = id_2_name_dict[s], categories=id_2_category_dict[json_data['links'][idx]['source']])
        G.add_node(t,name=id_2_name_dict[t], categories=id_2_category_dict[json_data['links'][idx]['target']])
        w = 0
        if json_data['links'][idx]['edgetype'].split(':')[1] != '':
            w = len(json_data['links'][idx]['edgetype'].split(':')[1].split(','))/len(all_pmids)
        t_str = json_data['nodes'][json_data['links'][idx]['source']]['name']+'-->'+json_data['nodes'][json_data['links'][idx]['target']]['name']
        edge_weights_from_original[t_str] = w
        G.add_edge(*e,value=w)
        source.append(s)
        target.append(t)
    edges = pd.DataFrame({'source':source,'target':target})
    G.number_of_nodes(),G.number_of_edges()

    # initiate forceatlas algorithm to calculate node positions
    forceatlas2 = ForceAtlas2(
                            # Behavior alternatives
                            outboundAttractionDistribution=True,  # Dissuade hubs
                            linLogMode=False,  # NOT IMPLEMENTED
                            adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                            edgeWeightInfluence=10.0,

                            # Performance
                            jitterTolerance=10.0,  # Tolerance
                            barnesHutOptimize=True,
                            barnesHutTheta=1.2,
                            multiThreaded=False,  # NOT IMPLEMENTED

                            # Tuning
                            scalingRatio=15.0,
                            strongGravityMode=False,
                            gravity=50,

                            # Log
                            verbose=True)

    # calculate nodes' positions
    positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=10000)

    # let's plot the same graph via holoviews
    '''
    click on nodes to see the details
    '''
    hv.Graph.from_networkx(G,positions,kwargs={'weight':'weight'}).opts(directed=True,node_size=10,arrowhead_length=0.009,tools=['hover']).options(color_index='categories', cmap='Category10')

    # return original edge weights
    return edge_weights_from_original




# START HERE
'''
This is the public function to calculate edge weights

Args:
sub_graph_json: a json object already read in from a file
email: a valid email for api calls
json_name: name of json file or relationship

Returns:
dataframe: dataframe containing each relationship and the various edge weights
Also creates a csv from the dataframe save to a folder called {json_name}
'''

def calculate_edge_weights(sub_graph_json, email, json_name, apikey=None):
    global id_2_name_dict
    global id_2_category_dict
    global all_pmids
    global json_data
    global pmid_dict
    global rel_name
    global api_key

    rel_name = json_name
    api_key = apikey

    id_2_name_dict = dict()
    id_2_category_dict = dict()
    all_pmids = []
    json_data = dict()
    pmid_dict = dict()


    json_data = json.loads(sub_graph_json)
    Entrez.email = email

    if not (os.path.exists(f'./{rel_name}')):
        os.mkdir(f'./{rel_name}')
    
    # Book-Keeping of Node Attributes

    for idx in range(len(json_data['nodes'])):
        id_2_name_dict[json_data['nodes'][idx]['id']] = json_data['nodes'][idx]['name']
        id_2_category_dict[json_data['nodes'][idx]['id']] = json_data['nodes'][idx]['category']

    # List all pmids from json file


    print('Creating a List of All PMIDS in Subgraph')

    for idx1 in range(len(json_data['links'])):
        if json_data['links'][idx1]['edgetype'].split(':')[1] != '':
            pmids = json_data['links'][idx1]['edgetype'].split(':')[1].split(',')
            for idx2 in range(len(pmids)):
                pmid = pmids[idx2]
                if pmid not in all_pmids and (len(pmid) > 1):
                    all_pmids.append(pmid)

    
    # Collect required information from NCBI database

    print('Collecting PMID Metadata from NCBI')

    collect_NCBI()

    # Edge weight based on original calculation
    print('-'*20)
    print('Calculating Original Edge Weights...')

    edge_weight_original = get_original_edge_weights()

    # Edge weight based on # of citations (Boltzmann Citation Weight)
    print('-'*20) 
    print('Calculating Citation Edge Weights...')

    edge_weight_citations = get_citation_edge_weights()
    
    # Edge weight based on paper publication date (Publication Year Weight)
    print('-'*20)
    print('Calculating Publication-Data Edge Weights...')

    edge_weight_publications = get_publication_edge_weights()
    
    # Combine above methods (Publication Year and Citation Weight)
    print('-'*20)
    print('Calculating Combo Citation/Publication Edge Weights...')

    edge_weight_citation_and_pub = get_publication_citation_weights()
    
    # Edge weight based on h-index of author (H-Index Weight)
    print('-'*20)
    print('Calculating Author H-Index Edge Weights...')

    edge_weight_h_index = get_h_index_weights(email)
    
    # Dump edge weights into CSV file
    print('-'*20)
    print('Dumping Edge Weights to CSV File...')

    keys = edge_weight_citation_and_pub.keys()
    edge,w_original, w_citation,w_pub_year,w_citation_pub_year,w_h_index = [],[],[],[],[],[]
    for key in keys:
        edge.append(key)
        w_original.append(float(edge_weight_original[key]))
        w_citation.append(float(edge_weight_citations[key]))
        w_pub_year.append(float(edge_weight_publications[key]))
        w_citation_pub_year.append(float(edge_weight_citation_and_pub[key]))
        w_h_index.append(float(edge_weight_h_index[key]))
    
    df = pd.DataFrame({'Edge':edge, 'Original_Weight':w_original, 'Boltzmann_Citation_Weight':w_citation,
                   'Publication_Year_Weight':w_pub_year,'Publication_Year_and_Citation_Weight':w_citation_pub_year, 'H_Index_Weight': w_h_index})
    
    df.to_csv(f'./{rel_name}/{rel_name}.csv',index=False)
    
    return df
