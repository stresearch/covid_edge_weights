### Functionality

#### edge_weight.py

Calculates edge weights using:
- Original Weight
- Boltzman Citation Weight
- Publication Year Weight
- Combination of Citation and Publication
- H-Index of Author

Returns a pandas dataframe with edges as rows and weights as columns and create a .csv


#### plot_subgraph.py

Creates a subgraph image and html using the weights calculated in edge_weight.py (requires output .csv)


#### get_pmids.py

Newer Blender subgraphs do not contain the pmids with evidence for graph edges. This code will take a
Blender subgraph and use the Knowledge Graph folder to add the pmids as expected by edge_weight.py

#### Usage

##### Usage

```python
import get_pmids
import edge_weight
import plot_subgraph.plot_subgraph as ps

# If PMIDS are not found in json subgraph, use get_pmids to create a json with pmids
# Requires locally downloaded Knowledge Graph folder from http://blender.cs.illinois.edu/covid19/
# If you are planning on converting several jsons, do so with one call of this function for speed
json_path_list = ['./sample_subgraph1.json', './sample_subgraph2.json', './desired_subgraph.json']
local_knowledge_graph_path = './KG'

get_pmids.get_pmids(json_path_list, local_knowledge_graph_path)
# New jsons will be saved with a file name ending of _with_pmids.json

jd = None
with open('desired_subgraph_with_pmids.json') as f:
    jd = f.read()
    
dataframe = edge_weight.calculate_edge_weights(jd,'john.doe@email.com','desired_subgraph_name', 'api key from NCBI (optional)')

ps.plot_subgraph('desired_subgraph_name/desired_subgraph.csv', 'desired_subgraph_name')


```