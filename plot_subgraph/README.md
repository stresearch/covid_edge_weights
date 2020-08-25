# Generate Interactive Visualization of Subgraphs
Read in a subgraph file as a csv and generate a static html with a dropdown menu of different weight types
### Setup Environment
>set up conda environment:
```bash
conda create -n covid python=3.7
conda activate covid
conda install -c pyviz holoviews pandas networkx
```
>running a notebook
```bash
conda install jupyterlab
jupyter labextension install @pyviz/jupyterlab_pyviz
jupyter lab
```



```python
import holoviews as hv
from holoviews.operation.datashader import bundle_graph
from colorcet import kb, kr
import numpy as np
import networkx as nx
import pandas as pd
from IPython import get_ipython
IS_NOTEBOOK = "IPKernelApp"  in get_ipython().config
```

## Read in the data
I read in a csv of edges with different edge weights for every edge. Ideally, Cem will also write each node (left-edge->right) as a seperate column.  
### Create a nodes dataframe


```python
data = pd.read_csv("https://855da60d-505b-4eee-942c-e19fb87dcc5f.s3.amazonaws.com/covid/2hop_benazepril_COVID-19.csv")

#this is since I have CSV that does not explicitly specify node types
data["left"] = data.Edge.apply(lambda a: a.split("-->")[0])
data["right"] = data.Edge.apply(lambda a: a.split("-->")[1])

#get all unique nodes
nodes = pd.concat([data.left,data.right],ignore_index=True).drop_duplicates().reset_index(drop=True).to_frame("name").reset_index()
data["left_index"] = nodes.set_index("name").loc[data.left,"index"].values
data["right_index"] = nodes.set_index("name").loc[data.right,"index"].values

#also temporary type should be in the csv
nodes["type"] = nodes.name.apply(lambda a: "Chem/Gene" if a in data.left.values else "Disease")
```

### Create NetworkX Graph
Make a directed graph from a list of edges.  
Also let's compute node layout as well... layout is just position of every node on the canvas


```python
graph = nx.DiGraph()
graph.add_edges_from(data.loc[:,["left_index","right_index"]].values)
layout = nx.layout.fruchterman_reingold_layout(graph.to_undirected(),.5)
layout = pd.DataFrame(layout,index=["x","y"]).T.sort_index()
nodes = pd.concat([nodes,layout],axis=1)
```

## Let's Plot Some Graphs
We have the following dataframes


```python
data.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Edge</th>
      <th>Original_Weight</th>
      <th>Boltzmann_Citation_Weight</th>
      <th>Publication_Year_Weight</th>
      <th>Publication_Year_and_Citation_Weight</th>
      <th>left</th>
      <th>right</th>
      <th>left_index</th>
      <th>right_index</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>gap junction protein alpha 1--&gt;benazepril</td>
      <td>0.0058</td>
      <td>9.820138e-01</td>
      <td>0.603353</td>
      <td>0.792683</td>
      <td>gap junction protein alpha 1</td>
      <td>benazepril</td>
      <td>0</td>
      <td>29</td>
    </tr>
    <tr>
      <th>1</th>
      <td>gap junction protein alpha 1--&gt;Tetradecanoylph...</td>
      <td>0.0095</td>
      <td>5.268630e-11</td>
      <td>0.123455</td>
      <td>0.061728</td>
      <td>gap junction protein alpha 1</td>
      <td>Tetradecanoylphorbol Acetate</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Tetradecanoylphorbol Acetate--&gt;COVID-19</td>
      <td>0.0058</td>
      <td>0.000000e+00</td>
      <td>1.000000</td>
      <td>0.500000</td>
      <td>Tetradecanoylphorbol Acetate</td>
      <td>COVID-19</td>
      <td>1</td>
      <td>35</td>
    </tr>
    <tr>
      <th>3</th>
      <td>angiotensin II receptor type 1--&gt;benazepril</td>
      <td>0.0064</td>
      <td>9.906840e-01</td>
      <td>0.164537</td>
      <td>0.577610</td>
      <td>angiotensin II receptor type 1</td>
      <td>benazepril</td>
      <td>2</td>
      <td>29</td>
    </tr>
    <tr>
      <th>4</th>
      <td>angiotensin II receptor type 1--&gt;Fadrozole</td>
      <td>0.0064</td>
      <td>4.742587e-02</td>
      <td>0.603353</td>
      <td>0.325389</td>
      <td>angiotensin II receptor type 1</td>
      <td>Fadrozole</td>
      <td>2</td>
      <td>3</td>
    </tr>
  </tbody>
</table>
</div>




```python
nodes.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>index</th>
      <th>name</th>
      <th>type</th>
      <th>x</th>
      <th>y</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>gap junction protein alpha 1</td>
      <td>Chem/Gene</td>
      <td>0.042898</td>
      <td>-0.577032</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>Tetradecanoylphorbol Acetate</td>
      <td>Chem/Gene</td>
      <td>-0.416220</td>
      <td>-0.918435</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>angiotensin II receptor type 1</td>
      <td>Chem/Gene</td>
      <td>-0.006889</td>
      <td>0.122234</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>Fadrozole</td>
      <td>Chem/Gene</td>
      <td>-0.352456</td>
      <td>0.410615</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>adrenoceptor beta 2</td>
      <td>Chem/Gene</td>
      <td>0.530451</td>
      <td>0.114279</td>
    </tr>
  </tbody>
</table>
</div>



First let's make some curvy edge paths (very necessary)


```python
def make_paths(nodes,data,curvature=.1):
    x = nodes.loc[data.left_index,["x","y"]].values
    y = nodes.loc[data.right_index,["x","y"]].values
    mid = .5*x + .5*y
    d = y-x
    dnorm = d/np.linalg.norm(d,axis=1, keepdims=True)
    p = x - (x*dnorm).sum(-1,keepdims=True)*dnorm
    pnorm = p/np.linalg.norm(p,axis=1,keepdims=True)
    control = mid + curvature*pnorm
    steps=np.linspace(0, 1, 100)[None,:,None]
    paths = (1-steps)**2*x[:,None,:] + 2*(1-steps)*steps*control[:,None,:]+steps**2*y[:,None,:]
    return list(paths)

paths = make_paths(nodes,data)
```

## Do the plotting
You can read more here: https://holoviews.org/user_guide/Network_Graphs.html


```python
hv.extension("bokeh")
```


```python
nodes_plot = hv.Nodes(nodes,kdims=["x","y","index"],vdims=["name","type"]).opts(show_legend=True)
edges_plot = hv.Graph((data,nodes_plot,paths),kdims=["left_index","right_index"],vdims=["Original_Weight","Boltzmann_Citation_Weight","Publication_Year_Weight","Publication_Year_and_Citation_Weight"])
# edges_plot = bundle_graph(edges_plot)


labels = hv.Labels(nodes,kdims=["x","y"],vdims=["name"]).opts(text_align="left", text_font_size="8pt",xoffset=.05,text_alpha = .5, text_color="white")

weight_field = "Boltzmann_Citation_Weight"

weights = ["Original_Weight","Boltzmann_Citation_Weight","Publication_Year_Weight","Publication_Year_and_Citation_Weight"]

# need to change this
node_type_size_dict = {"Chem/Gene":15,"Disease":30}

def make_plot(weight_field):
    return edges_plot.opts(clone = True, node_color = "type", node_line_width=0, node_size = hv.dim("type").categorize(node_type_size_dict),
                node_cmap="Category10_r",
                edge_color=weight_field,edge_alpha=.9,
                edge_line_width = 1.5*hv.dim(weight_field).norm()+.3,
                edge_cmap="Reds", xaxis=None, yaxis=None, width=1000, height=800,
                colorbar=True,  inspection_policy="edges" ,bgcolor="#282829")*labels



out = hv.HoloMap({w: make_plot(w) for w in weights},kdims = "weight_type")
out
```


```python
# from IPython.display import HTML
# HTML(open("subgraph.html").read())
```

save it


```python
hv.save(out,"subgraph.html")
if not IS_NOTEBOOK:
    hv.save(out,"img1.gif")
```

![Subgraph](imgs/img1.gif)

We can also see all of them side by side


```python
out2 = out.opts(width=600,height=600).layout("weight_type").cols(2)
if not IS_NOTEBOOK:
    hv.save(out2,"imgs/img2.png")
hv.save(out2,"subgraph_side_by_side.html")
out2
```

![Side by Side](imgs/img2.png)
