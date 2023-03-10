#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 10:30:03 2023

@author: robertbenassai
"""

import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
bcn_edges = gpd.read_file(path_bcn, crs="EPSG:4326")

bcn_nodes = gpd.read_file('/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/vertices.shp',
                          crs="EPSG:4326")

missing_edges = pd.read_csv(
    '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/missing_edges.csv')
bcn_distr = gpd.read_file('/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp',
                          crs="EPSG:4326")
nodes_w_distr = gpd.sjoin(bcn_nodes, bcn_distr,
                        how="left", predicate='within')
distr_df = nodes_w_distr.loc[nodes_w_distr['DISTRICTE'] == 
                             "{:02d}".format(1)]
 # Construct the primal graph
#%%
from shapely.geometry import LineString

merged = pd.merge(bcn_nodes, missing_edges, how = 'right', left_on='node_id', 
                  right_on='FROM_NODE')
merged_to = pd.merge(bcn_nodes, missing_edges, how = 'right', left_on='node_id', 
                  right_on='TO_NODE')
merged['geometry_to'] = merged_to.geometry

missing_edge_gpd = gpd.GeoDataFrame()
missing_edge_gpd['FROM_NODE'] = merged.FROM_NODE
missing_edge_gpd['TO_NODE'] = merged.TO_NODE
missing_edge_gpd['geometry'] = merged.apply(lambda row: LineString([row['geometry'], row['geometry_to']]), axis=1) #Create a linestring column
missing_edge_gpd['length (m)'] = missing_edge_gpd.apply(lambda row: row['geometry'].length, axis=1) #Create a linestring column

eid = np.linspace(38625, 38625+49, 50)

missing_edge_gpd['eid'] = eid

missing_edge_gpd.to_file('/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/missing_edges.shp')
#%%
import momepy
distr = nx.DiGraph()
distr.add_nodes_from(distr_df['uid'])
ind_to_pos = {nodes_w_distr.loc[i,'uid']:(nodes_w_distr.loc[i,'geometry'].x,
                                      nodes_w_distr.loc[i,'geometry'].y) for i in 
              range(len(nodes_w_distr))}
#nx.set_node_attributes(bcn_graph, ind_to_pos, 'pos')

for edge in bcn_edges['to_from']:
    sep = edge.split(',')
    if len(sep) == 2:
        if int(sep[0]) and int(sep[1]) in distr.nodes:
            distr.add_edge(int(sep[0]), int(sep[1]))
            
for i,j in zip(missing_edges['FROM_NODE'], missing_edges['TO_NODE']):
    if (i and j) in distr.nodes:
        distr.add_edge(i, j)
        
#%%
import matplotlib.pyplot as plt
plt.figure()
nx.draw_networkx(distr, with_labels=False, pos = ind_to_pos, node_size = 3)


