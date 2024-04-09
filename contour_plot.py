#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:02:05 2024

@author: robertbenassai
"""

import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.tri as tri
import networkx as nx
import geopandas as gpd
import pandas as pd
import contextily as ctx
from shapely.geometry import Point, LineString, Polygon

def distr_to_nx(distr_ind:int, path_edges: str, path_distr: str, edge_attr 
                = np.empty(shape = (1)), name :str = 'length', 
                path_nodes: str = '/Users/robertbenassai/Documents/UOC/'+
                'project_HHD/xarxaneta/data/paris/shp/nodes/nodes.shp'):
    
    # edges and nodes of the whole graph
    bcn_edges = gpd.read_file(path_edges, crs = "ESPG:2154") #crs="EPSG:25831")


    bcn_nodes = gpd.read_file(path_nodes, crs = "ESPG:2154")
                          #crs="EPSG:25831")
    
    if distr_ind == 0: # graph of the entire city
    
        # initialize graph of the district
        distr = nx.DiGraph()
        distr.add_nodes_from(bcn_nodes['uid'])
        ind_to_pos = {int(bcn_nodes.loc[i,'uid']):(bcn_nodes.loc[i,'geometry'].x,
                                              bcn_nodes.loc[i,'geometry'].y) 
                      for i in range(len(bcn_nodes))}
    
        for i, j in zip(bcn_edges['i'], bcn_edges['j']):
            
            distr.add_edge(int(i), int(j), length =\
                            np.linalg.norm(np.array(ind_to_pos[int(i)]) - 
                                          np.array(ind_to_pos[int(j)])))

        # for e in bcn_edges['to_from']:
        #     if e == None:
        #         continue
        #     ind = ''
        #     for a in e:
                
        #         if a.isdigit():
        #             ind += a
        #         if a == ',':
        #             i = int(ind)
        #             ind = ''
        #     j = int(ind)
            
        #     if i == '' or j == '': #sometimes edges only have 1 node
        #         continue
            
        #     distr.add_edge(int(i), int(j), length =\
        #                     np.linalg.norm(np.array(ind_to_pos[int(i)]) - 
        #                                   np.array(ind_to_pos[int(j)])))
        #     i = j = ''
        
    else: # get the graph of a given district

        # districts
        bcn_distr = gpd.read_file(path_distr,
                                  crs="EPSG:25831")
        
        # joined dataframes (districts and nodes)
        nodes_w_distr = gpd.sjoin(bcn_nodes, bcn_distr,
                                how="left", predicate='within')
        
        # select only the nodes of the given district
        distr_df = nodes_w_distr.loc[nodes_w_distr['DISTRICTE'] == 
                                     "{:02d}".format(distr_ind)]
        
        # initialize graph of the district
        distr = nx.DiGraph()
        distr.add_nodes_from(distr_df['uid'])
        ind_to_pos = {nodes_w_distr.loc[i,'uid']:(nodes_w_distr.loc[i,'geometry'].x,
                                              nodes_w_distr.loc[i,'geometry'].y) for i in 
                      range(len(nodes_w_distr))}
    
    
        for i, j in zip(bcn_edges['i'], bcn_edges['j']):
            if i and j in distr.nodes:
                distr.add_edge(int(i), int(j), length =\
                           np.linalg.norm(np.array(ind_to_pos[int(i)]) - 
                                          np.array(ind_to_pos[int(j)])))
                    
        if np.size(edge_attr) != 0:
            
            new_attr = {}
            
            for i, j in zip(bcn_edges['i'], bcn_edges['j']):
                if i and j in distr.nodes:
                    new_attr[(int(i),int(j))] = edge_attr[int(i), int(j)]
            
            nx.set_edge_attributes(distr, new_attr, name = name)
            
    #checking for self loops
    self_loops = []
    print('checking for self loops...')
    for edge in distr.edges:
        if edge[0] == edge[1]:
            print('self loop found', edge)
            self_loops.append(edge)
    distr.remove_edges_from(self_loops)
    
    print('checking for disconnected components...')
    remove_nodes = list(nx.connected_components(distr.to_undirected()))
    
    print('total number of connected components', len(remove_nodes))
        
    
    remove_nodes.remove(max(remove_nodes, key = len))
    for node_set in remove_nodes:
        distr.remove_nodes_from(node_set)
    
    #remove nodes of degree 1:
    while True:
        
        deg_1_nodes = [i for i,deg in dict(nx.degree(distr)).items() if deg == 1]
        
        if len(deg_1_nodes) == 0:
            break
        distr.remove_nodes_from(deg_1_nodes)

    nx.set_node_attributes(distr, ind_to_pos, name = 'pos')
    
    #relabeling nodes to ascending order
    
    dict_relabel = {node: i for i, node in enumerate(sorted(distr.nodes))}
    nx.relabel_nodes(distr, dict_relabel, copy = False)    
    ind_to_pos_upd = {dict_relabel[node]: ind_to_pos[node] for node in 
                      dict_relabel.keys()}
    
    # print('plotting final graph...')

    f, ax = plt.subplots(figsize=(12, 6))
    nx.draw(distr, pos = ind_to_pos_upd, node_size = 2, ax = ax, with_labels=False)
    print('total number of nodes: ', len(distr.nodes))
    return(distr, ind_to_pos_upd)



distr_ind = 0

path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'
# path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/data/paris/shp/edges/edges.shp'
p_nodes = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/nodes_clean_net_willum.shp'

bcn_graph, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr, path_nodes = p_nodes)

x = []
y = []


for i in range(len(bcn_graph.nodes)):
    xy = ind_to_pos[i]
    x.append(xy[0])
    y.append(xy[1])

with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/bcn/pot.csv', 'r') as tot_f:
    reader = csv.reader(tot_f)
    z = [-float(row[1]) for row in reader]
    
with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/bcn/w_g.csv', 'r') as tot_f:
    reader = csv.reader(tot_f)
    gcomp = {(int(row[0]), int(row[1])) if float(row[2]) > 0 else (int(row[1]), int(row[0])) : np.abs(float(row[2]))  for row in reader}

#%% POTENTIAL CONTOUR PLOT

#getting gpd dataframe
# Convert nodes to GeoDataFrame
nodes_df = pd.DataFrame(list(ind_to_pos.items()), columns=['node', 'pos'])
geometry = gpd.points_from_xy(nodes_df['pos'].apply(lambda x: x[0]), nodes_df['pos'].apply(lambda x: x[1]))
nodes_gdf = gpd.GeoDataFrame(nodes_df, geometry=geometry)



fig, ax2 = plt.subplots(nrows=1, ncols = 1, figsize = (14, 10))
ngridx = 100

ngridy = 200
# -----------------------
# Interpolation on a grid
# -----------------------
# A contour plot of irregularly spaced data coordinates
# via interpolation on a grid.
z = np.array(z)

z_pos = z[z>0]
z_neg = z[z<0]
print(np.max(z_neg))
print(len(z_neg), len(z_pos))
percentile_pot_pos = np.percentile(z_pos, 90)
percentile_pot_neg = np.percentile(z_neg, 10)
print(percentile_pot_neg, percentile_pot_pos)
# Create grid values first.
grid_xlim_up = np.max([np.min(x), np.max(x)])
grid_xlim_lo = np.min([np.min(x), np.max(x)])

grid_ylim_up = np.max([np.min(y), np.max(y)])
grid_ylim_lo = np.min([np.min(y), np.max(y)])

xi = np.linspace(grid_xlim_lo, grid_xlim_up, ngridx)
yi = np.linspace(grid_ylim_lo, grid_ylim_up, ngridy)

# Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
triang = tri.Triangulation(x, y)
interpolator = tri.LinearTriInterpolator(triang, z)
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)

# Note that scipy.interpolate provides means to interpolate data on a grid
# as well. The following would be an alternative to the four lines above:
# from scipy.interpolate import griddata
# zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear')
cmap= "RdBu_r"
# cmap = "PuOr_r"

ax2.contour(xi, yi, zi, levels=17, vmin = percentile_pot_neg, 
                     vmax = percentile_pot_pos, cmap = cmap, alpha=0.6)
ax2.contourf(xi, yi, zi, levels=17, vmin = percentile_pot_neg, 
                     vmax = percentile_pot_pos, cmap = cmap, alpha=0.1)


sm = plt.cm.ScalarMappable(cmap=cmap, 
                            norm=plt.Normalize(vmin=percentile_pot_neg, vmax=percentile_pot_pos))
sm._A = []
cbar = fig.colorbar(sm, ax=ax2)
cbar.set_label(r'node potential', fontsize = 18)
cbar.ax.tick_params(labelsize=18)


providers = ctx.providers.flatten()
selection = ['OpenStreetMap.Mapnik',
             'OpenStreetMap.HOT',
             'CartoDB.Positron',
             'Stamen.Terrain',
             'CartoDB.VoyagerNoLabels',
             'CartoDB.DarkMatter',
             'NASAGIBS.ASTER_GDEM_Greyscale_Shaded_Relief',
             'NASAGIBS.ViirsEarthAtNight2012',
             'NASAGIBS.BlueMarble',
             'Thunderforest.MobileAtlas'
            ]

ctx.add_basemap(ax2, crs="EPSG:25831", source=providers[selection[4]])
# ctx.add_basemap(ax2, crs="EPSG:2154", source=providers[selection[4]])

nodes_gdf.plot(ax=ax2, markersize=0.5, alpha=0.3, label='Nodes', color = 'gray')
ax2.axis('off')

ax2.set_xlim(grid_xlim_lo, grid_xlim_up)
ax2.set_ylim(grid_ylim_lo, grid_ylim_up)

plt.tight_layout()
plt.show()

#%% gradient comp plot

def average_field(G, gcomp, pos):
    g_edges = list(gcomp.keys())
    
    avg_field = {}
    for node in G.nodes():
        #first check out edges
        x = pos[node][0]
        y = pos[node][1]
        
        x_ls = []
        y_ls = []
        norm_ls = []
        for out_e in G.out_edges(node):
            
            if out_e in g_edges:
                if gcomp[out_e] > 0:
                    
                    xi = pos[out_e[1]][0]
                    yi = pos[out_e[1]][1]
                    
                    norm = np.sqrt((xi-x)**2+(yi-y)**2)
                    alph = np.arctan((yi-y)/(xi-x))
                    
                    xf = np.cos(alph) * gcomp[out_e]
                    yf = np.sin(alph) * gcomp[out_e]
                                         
                    x_ls.append(xf)
                    y_ls.append(yf)
                    norm_ls.append(norm)
                    
        for in_e in G.in_edges(node):
            
            if in_e in g_edges:
                if gcomp[in_e] < 0:
                    xi = pos[in_e[0]][0]
                    yi = pos[in_e[0]][1]
                    
                    norm = np.sqrt((xi-x)**2+(yi-y)**2)
                    alph = np.arctan((yi-y)/(xi-x))
                    
                    xf = np.cos(alph) * gcomp[out_e]
                    yf = np.sin(alph) * gcomp[out_e]
            
                    x_ls.append(xf)
                    y_ls.append(yf)
                    norm_ls.append(norm)
        
        avg_norm = np.mean(norm_ls)
        
        new_norm = np.sqrt((np.mean(x_ls))**2+(np.mean(y_ls))**2)
        
        avg_x = x + np.mean(x_ls)*avg_norm/new_norm
        avg_y = y + np.mean(y_ls)*avg_norm/new_norm
        
        if np.isnan(avg_x) or np.isnan(avg_y):
            avg_x = x
            avg_y = y
            new_norm = 0
        avg_field[node] = (avg_x, avg_y, new_norm)
    
    return(avg_field)
        
#%%
avg_field = average_field(bcn_graph, gcomp, ind_to_pos)
#%%
# Convert edges to GeoDataFrame
# edges_df = pd.DataFrame(list(bcn_graph.edges()), columns=['node1', 'node2'])
# geometry_edges = [LineString([ind_to_pos[edge[0]], ind_to_pos[edge[1]]]) for edge in bcn_graph.edges()]
# edges_gdf = gpd.GeoDataFrame(edges_df, geometry=geometry_edges)

edges_df = pd.DataFrame()
geometry_edges = [LineString([ind_to_pos[node], (avg_field[node][0], avg_field[node][1])]) for node in bcn_graph.nodes()]
edges_gdf = gpd.GeoDataFrame(edges_df, geometry=geometry_edges)

cmap = plt.cm.get_cmap('inferno').copy()
cmap.set_under(color='gray', alpha = 0.5)    
edges_gdf['avg_grad_comp'] = {i: avg_field[node][2] for i, node in enumerate(bcn_graph.nodes)}

fig, ax = plt.subplots(nrows=1, ncols = 1, figsize = (10, 14))

max_val = np.percentile(list(gcomp.values()), 95)
min_val = np.percentile(list(gcomp.values()), 10)

edges_gdf.plot(column ='avg_grad_comp', ax=ax, markersize=0.5, alpha=1, 
               label='edges', cmap = cmap, vmin = min_val, vmax = max_val)
ctx.add_basemap(ax, crs="EPSG:2154", source=providers[selection[4]])

sm = plt.cm.ScalarMappable(cmap=cmap, 
                            norm=plt.Normalize(vmin= min_val, vmax=max_val))
sm._A = []
cbar = plt.colorbar(sm, ax = ax)
cbar.set_label('Gradient Component', fontsize = 18)
cbar.ax.tick_params(labelsize=18)


ax.axis('off')

ax.set_xlim(grid_xlim_lo, grid_xlim_up)
ax.set_ylim(grid_ylim_lo, grid_ylim_up)

plt.tight_layout()
plt.show()

#%%
# ----------
# Tricontour
# ----------
# Directly supply the unordered, irregularly spaced coordinates
# to tricontour.
fig, ax2 = plt.subplots(nrows=1, figsize = (12, 10))

percentile_pot = np.percentile(z, 99)
print(percentile_pot)
cmap= "RdBu_r"
ax2.tricontourf(x, y, z, levels=30,vmin = -percentile_pot, vmax = percentile_pot, cmap = cmap)


sm = plt.cm.ScalarMappable(cmap=cmap, 
                            norm=plt.Normalize(vmin=-percentile_pot, vmax=percentile_pot))
sm._A = []
cbar = fig.colorbar(sm, ax=ax2)
cbar.set_label(r'node potential', fontsize = 18)
cbar.ax.tick_params(labelsize=18)

ax2.set_xlim(grid_xlim_lo, grid_xlim_up)
ax2.set_ylim(grid_ylim_lo, grid_ylim_up)
ax2.plot(x, y, 'ko', ms=0.5, alpha = 0.3)
ax2.axis('off')

plt.show()

#%% hexbin of gradient component
import matplotlib.colors as colors
import xyzservices.providers as xyz

#first we find the midpoints of the edges
xpos = []
ypos = []
mag = []

for edge in gcomp.keys():
    
    x0, y0 = ind_to_pos[edge[0]]
    xf, yf = ind_to_pos[edge[1]]
    
    norm = np.sqrt((xf-x0)**2 + (yf-y0)**2)
    alph = np.arctan((yf-y0)/(xf-x0))
                    
    x_mid = np.cos(alph) * norm*0.5 + x0
    y_mid = np.sin(alph) * norm*0.5 + y0
    
    xpos.append(x_mid)
    ypos.append(y_mid)
    
    mag.append(np.abs(gcomp[edge]))


# --------
# PLOT
# --------

cmap = plt.colormaps['bone_r']
cmap.set_under(color='gray', alpha = 0)    

gridsize = 150
mincnt = 1

max_val = np.percentile(mag, 90)
c_norm = colors.Normalize(vmin=0, vmax=max_val)


fig, ax = plt.subplots(figsize=(8, 14))
hb = ax.hexbin(xpos, ypos, C=mag, gridsize=gridsize, cmap=cmap, norm=c_norm, mincnt=mincnt, alpha=0.5)

cmap2= "RdBu_r"

# ax.contour(xi, yi, zi, levels=17, vmin = percentile_pot_neg, 
#                      vmax = percentile_pot_pos, cmap = cmap2, alpha=1)
ctx.add_basemap(ax, crs="EPSG:25831", source=xyz.OpenStreetMap.BZH)
# providers[selection[1]])

offs = hb.get_offsets()

cbar = plt.colorbar(hb, orientation='vertical')
# cbar.set_label(r'diff 20% - 0% tolerance BW', fontsize = 18)
cbar.set_label(r'abs gradient component', fontsize = 18)
cbar.ax.tick_params(labelsize=16)

ax.axis('off')
# Show the plot
plt.tight_layout()
plt.show()
    