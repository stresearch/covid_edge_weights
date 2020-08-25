#Copyright (C) Systems & Technology Research
#Use of this software is subject to the restrictions in license.txt
#Distribution A: Approved for Public Release, Distribution Unlimited
'''
DEMO FOR SITE VISIT:

Should be run from ./covid_edge_weights folder

Drag and drop ./sample/2hop_candesartan_COVID-19_dynamic and ./sample/2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic
into the ./covid_edge_weghts folder so that we can utilize the already calculated data (speedy)


Make sure KG is locally downloaded and you alter the code to identify its location on your computer
'''



'''
Demo Run 1: 2hop_candesartan_COVID-19_dynamic
'''

# Don't need to collect pmids for this one!

import edge_weight
import plot_subgraph.plot_subgraph as ps

print('-'*70)
print('-'*70)
print('-'*8 + '2hop_candesartan_COVID-19_dynamic Edge Weight Pipeline' + '-'*8)
print('-'*70)
print('-'*70)


# Load JSON File
jd = None
with open('2hop_candesartan_COVID-19_dynamic/2hop_candesartan_COVID-19_dynamic.json') as f:
    jd = f.read()

print('-'*30 + 'edge_weight.py:' + '-'*30)
# Calculate Edge Weights 
dataframe = edge_weight.calculate_edge_weights(jd,'lauren.lederer@stresearch.com','2hop_candesartan_COVID-19_dynamic', '8c764b8875ee55edc5e399779f09470b6e08')

print('-'*30 + 'plot_subgraph.py:' + '-'*30)
# Plot Subgraph
ps.plot_subgraph('2hop_candesartan_COVID-19_dynamic/2hop_candesartan_COVID-19_dynamic.csv', '2hop_candesartan_COVID-19_dynamic')

print('-'*30 + 'Complete!'+ '-'*30)
'''
 Demo Run 2: 2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic
'''
import get_pmids
import edge_weight
import plot_subgraph.plot_subgraph as ps

print('-'*70)
print('-'*70)

print('-2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic Edge Weight Pipeline-')
print('-'*70)
print('-'*70)


# If PMIDS are not found in json subgraph, use get_pmids to create a json with pmids
# Requires locally downloaded Knowledge Graph folder from http://blender.cs.illinois.edu/covid19/
# If you are planning on converting several jsons, do so with one call of this function for speed
json_path_list = ['./2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic/2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic.json']
local_knowledge_graph_path = '../KG'

print('-'*30 + 'get_pmids.py:'+ '-'*30)

get_pmids.get_pmids(json_path_list, local_knowledge_graph_path)
# New jsons will be saved with a file name ending of _with_pmids.json

jd = None
with open('./2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic/2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic_with_pmids.json') as f:
    jd = f.read()

print('-'*30 + 'edge_weight.py:'+ '-'*30)

dataframe = edge_weight.calculate_edge_weights(jd,'lauren.lederer@stresearch.com' ,'2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic', '8c764b8875ee55edc5e399779f09470b6e08')

print('-'*30 + 'plot_subgraph.py:'+ '-'*30)

ps.plot_subgraph('2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic/2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic.csv', '2hop_Vancomycin_cathepsin-S-L-homeolog_dynamic')

print('-'*30 + 'Complete!'+ '-'*30)




