#!/usr/bin/env python
# coding: utf-8
#Copyright (C) Systems & Technology Research
#Use of this software is subject to the restrictions in license.txt
#Distribution A: Approved for Public Release, Distribution Unlimited

import holoviews as hv
from holoviews.operation.datashader import bundle_graph
from colorcet import kb, kr
import numpy as np
import networkx as nx
import pandas as pd
import os

# also requires selenium
# conda install -c conda-forge firefox geckodriver

edges_plot = None
nodes_plot = None
node_type_size_dict = None
labels = None



# First let's make some curvy edge paths (very necessary)

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


def make_plot(weight_field):
    return edges_plot.opts(clone = True, node_color = "type", node_line_width=0, node_size = hv.dim("type").categorize(node_type_size_dict),
                node_cmap="Category10_r",
                edge_color=weight_field,edge_alpha=.9,
                edge_line_width = 1.5*hv.dim(weight_field).norm()+.3,
                edge_cmap="Reds", xaxis=None, yaxis=None, width=1000, height=800,
                colorbar=True,  inspection_policy="edges" ,bgcolor="#282829")*labels

# ## Read in the data
# I read in a csv of edges with different edge weights for every edge. Ideally, Cem will also write each node (left-edge->right) as a seperate column.  
# ### Create a nodes dataframe

def plot_subgraph(csv_path, rel_name):
    global edges_plot
    global nodes_plot
    global node_type_size_dict
    global labels

    with open(csv_path, 'r') as w:
        data = pd.read_csv(w)

    if not (os.path.exists(f'./{rel_name}')):
        os.mkdir(f'./{rel_name}')
    
    if not (os.path.exists(f'./{rel_name}/subgraphs')):
        os.mkdir(f'./{rel_name}/subgraphs')
    
    #this is since I have CSV that does not explicitly specify node types
    data["left"] = data.Edge.apply(lambda a: a.split("-->")[0])
    data["right"] = data.Edge.apply(lambda a: a.split("-->")[1])

    #get all unique nodes
    nodes = pd.concat([data.left,data.right],ignore_index=True).drop_duplicates().reset_index(drop=True).to_frame("name").reset_index()
    data["left_index"] = nodes.set_index("name").loc[data.left,"index"].values
    data["right_index"] = nodes.set_index("name").loc[data.right,"index"].values

    #also temporary type should be in the csv
    nodes["type"] = nodes.name.apply(lambda a: "Chem/Gene" if a in data.left.values else "Disease")


    # ### Create NetworkX Graph
    # Make a directed graph from a list of edges.  
    # Also let's compute node layout as well... layout is just position of every node on the canvas

    graph = nx.DiGraph()
    graph.add_edges_from(data.loc[:,["left_index","right_index"]].values)
    layout = nx.layout.fruchterman_reingold_layout(graph.to_undirected(),.5)
    layout = pd.DataFrame(layout,index=["x","y"]).T.sort_index()
    nodes = pd.concat([nodes,layout],axis=1)


    # ## Let's Plot Some Graphs
    paths = make_paths(nodes,data)

    hv.extension("bokeh")


    nodes_plot = hv.Nodes(nodes,kdims=["x","y","index"],vdims=["name","type"]).opts(show_legend=True)
    edges_plot = hv.Graph((data,nodes_plot,paths),kdims=["left_index","right_index"],vdims=["Original_Weight","Boltzmann_Citation_Weight","Publication_Year_Weight","Publication_Year_and_Citation_Weight", "H_Index_Weight"])
    # edges_plot = bundle_graph(edges_plot)


    labels = hv.Labels(nodes,kdims=["x","y"],vdims=["name"]).opts(text_align="left", text_font_size="8pt",xoffset=.05,text_alpha = .5, text_color="white")

    weight_field = "Boltzmann_Citation_Weight"

    weights = ["Original_Weight","Boltzmann_Citation_Weight","Publication_Year_Weight","Publication_Year_and_Citation_Weight", "H_Index_Weight"]

    # need to change this
    node_type_size_dict = {"Chem/Gene":15,"Disease":30}

    out = hv.HoloMap({w: make_plot(w) for w in weights},kdims = "weight_type")


    print('Saving Subgraphs to file...')

    hv.save(out,f"./{rel_name}/subgraphs/{rel_name}_subgraph.html")

    hv.save(out,f"./{rel_name}/subgraphs/{rel_name}_img1.gif")


    # ![Subgraph](imgs/img1.gif)

    # We can also see all of them side by side


    out2 = out.opts(width=600,height=600).layout("weight_type").cols(2)

    hv.save(out2,f"./{rel_name}/subgraphs/{rel_name}_img2.png")
    hv.save(out2,f"./{rel_name}/subgraphs/{rel_name}_subgraph_side_by_side.html")

    # ![Side by Side](imgs/img2.png)
