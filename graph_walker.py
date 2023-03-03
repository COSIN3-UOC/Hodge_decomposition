# -*- coding: utf-8 -*-
"""
Code author: Robert Benassai Dalmau

In this code, random walkers will be introduced to directed and undirected
graphs. The flows are going to be added to the edges and the Hodge
decomposition.
"""

import time
import json
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point, LineString
from libpysal import weights
from contextily import add_basemap
import momepy
import pandas as pd
import geopandas as gpd
import math
import scipy
from sympy import LeviCivita
from itertools import combinations, permutations
import EoN
import operator
from scipy.spatial import Delaunay
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np


'''
moves "n_walk" walkers/node in an undirected graph "G" "steps" times. The initial positions of the
walkers are randomly selected amongst the nodes. The function returns the graph
with the amount of times any walker has passed through each node as attributes
to the nodes

'''


def graph_walkers(G, steps, n_walk):
    path = []
    # let's add a list of n_walk random walkers
    initial_nodes = list(G.nodes)
    for i in range(1, n_walk):
        print(i)
        initial_nodes += initial_nodes
    # path stores the position of each walker in each time step
    path.append(initial_nodes)
    # weights counts the amount of times a walker has visited each node
    weights = {i_node: path[0].count(i_node) for i_node in list(G.nodes)}
    # now, let's move the walker to a connected node
    # first let's see the available nodes to go to
    for step in range(steps):
        # list of neighboring nodes of each walker
        neighbors = [[n for n in G.neighbors(walker)]
                                             for walker in initial_nodes]
        # now move the random walker to another node
        final_nodes = [random.choice(goto_node) for goto_node in neighbors]
        path.append(final_nodes)
        initial_nodes = final_nodes
        # count the occutpation of each node after the moves
        for i_node in list(G.nodes):
            weights[i_node] += final_nodes.count(i_node)
    # set node value as the number of visits of each value
    # print(weights)
    nx.set_node_attributes(G, weights, name='weights')
    return(G)


# %%STAR GRAPH
star = nx.star_graph(100)
walked = graph_walkers(star, 100, 1)
node_labels = nx.get_node_attributes(walked, 'weights')

color_p = list(node_labels.values())
colors = range(min(color_p), max(color_p))
cmap = plt.cm.Blues
vmin = min(colors)
vmax = 30
plt.title('100 walkers 10 steps')
pos = nx.spring_layout(walked)
nx.draw_networkx(walked, pos=pos, with_labels=False, node_color=color_p,
                 cmap=cmap, vmin=vmin, vmax=vmax)
nx.draw_networkx_labels(walked, pos, labels=node_labels)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label('# visits')
plt.show()
# %% DELAUNAY GRAPHS

'''
we will create 2 clusters of points normally distributed
'''
# x,y coords of points
nodes_cl1 = 100
nodes_cl2 = 100
clust_1 = np.random.normal(0, 1, size=(nodes_cl1, 2))
clust_2 = np.random.normal(5, 1, size=(nodes_cl2, 2))

clusters = np.concatenate((clust_1, clust_2), axis=0)
tri = Delaunay(clusters)
plt.figure()
plt.triplot(clusters[:, 0], clusters[:, 1], tri.simplices)
plt.plot(clust_1[:, 0], clust_1[:, 1], 'or', label='cluster 1')
plt.plot(clust_2[:, 0], clust_2[:, 1], 'ob', label='cluster 2')
plt.legend()
plt.show()

# create a set for edges that are indexes of the points
edges = set()
# for each Delaunay triangle
for n in range(tri.nsimplex):
    # for each edge of the triangle
    # sort the vertices
    # (sorting avoids duplicated edges being added to the set)
    # and add to the edges set
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 1]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 1], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))


# make a graph based on the Delaunay triangulation edges
graph = nx.Graph(list(edges))
# plt.figure()
# plot graph
# dictionary of node:position
pointIDXY = dict(zip(range(len(clusters)), clusters))
nx.set_node_attributes(graph, pointIDXY, name='pos')
# nx.draw(graph, pointIDXY)
# plt.show()


# %%
''' Now that we have created the Delaunay graph, we can put random walkers in it'''

n_walkers = 1
steps = 100
walked_D = graph_walkers(graph, steps, n_walkers)
# getting the labels (visits)
D_labels = nx.get_node_attributes(walked_D, 'weights')

# %% PLOT WITH COLOR GRADIENT
# plotting nodes with color gradient
color_p = list(D_labels.values())
colors = range(min(color_p), max(color_p))
fig, ax = plt.subplots()
cmap = plt.cm.Blues
vmin = min(colors)
vmax = max(colors)

nx.draw_networkx(walked_D, pointIDXY, with_labels=False, node_color=color_p,
                 cmap=cmap, vmin=vmin, vmax=vmax, node_size=100)

# nx.draw_networkx_labels(walked_D ,pointIDXY, labels = D_labels)
ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
plt.colorbar(sm)
plt.scatter([0, 5], [0, 5], s=100, color='r')
plt.show()

sorted_val = sorted(D_labels.items(), key=operator.itemgetter(1), reverse=True)
print(sorted_val)
# %% SIZE OF DOTS ACCORDING TO WALKERS
fig, ax = plt.subplots()

nx.draw_networkx(walked_D, pointIDXY, with_labels=False, node_color=color_p,
                 cmap=cmap, vmin=vmin, vmax=vmax,
                 node_size=[v * 1 for v in color_p])

# nx.draw_networkx_labels(walked_D ,pointIDXY, labels = D_labels)
ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
plt.colorbar(sm)
# plt.scatter([0,5], [0,5], s=100, color='r')
plt.show()

# %% node degree vs walkers in node
dict_list = []
num_walk_list = [100, 1000, 10000]
for n_walkers in num_walk_list:
    steps = 100
    walked_D = graph_walkers(graph, n_walkers, steps)
    # getting the labels (visits)
    D_labels = nx.get_node_attributes(walked_D, 'weights')
    deg_walk = [[walked_D.degree(i), D_labels[i]] for i in D_labels.keys()]

    degrees = sorted([i[1] for i in walked_D.degree()])
    degrees = [*set(degrees)]

    deg_dict = {}
    for deg in degrees:
        deg_list = []
        for i in deg_walk:
            if i[0] == deg:
                deg_list.append(i[1])
        deg_dict[deg] = np.array(deg_list)*100/((steps+1)*n_walkers)
        deg_list = []
    dict_list.append(deg_dict)
# %%  BOX PLOTS

fig, axs = plt.subplots(1, 3, figsize=(15, 5))
i = 0
for deg_dict, ax in zip(dict_list, axs.ravel()):
    ax.boxplot(deg_dict.values(), positions=list(deg_dict.keys()))
    ax.set_title('degree vs visits with ' +
                 str(num_walk_list[i])+' random walkers')
    # ax.set_xticklabels(deg_dict.keys())
    ax.set_xlabel('node degree')
    ax.set_ylabel('normalized visits per node (100%)')
    i += 1
plt.show()

# %%

'''VECTOR FIELDS USING RANDOM WALKER TRAJECTORIES'''

''' first we need to assign the distances to each edge'''

for node in graph.nodes:
# first neighbors of each node that have not been considered yet (avoid 12 21
# repetitions)
    first_n = [n for n in (list(graph.successors(node)) +
                           list(graph.predecessors)) if n > node]
    dist = {(node, n): np.linalg.norm(graph.nodes[n]['pos']-graph.nodes[node]['pos'])
            for n in first_n}
    nx.set_edge_attributes(graph, dist, name='dist')
print(len(nx.get_edge_attributes(graph, 'dist')), len(graph.edges))
# %%
a = 5
print(list(graph.neighbors(a)))
a_dists = np.array([graph[i[0]][i[1]]['dist'] for i in graph.edges(a)])
print(graph.edges(a))
print(np.where(a_dists/2 < 0.3))
print(graph.nodes[a]['pos'])
# %% NODE_WALKERS (walk during Dt) (time/n_wak = 25s for clusters of 100 points)
'''Function that puts n = len(list(G.nodes)) random walkers in a Digraph
transitable in both directions and moves them during Dt time. The walkers always
move until they are not able to go anywhere else with the time left'''

'''G: nx.DiGraph
   Dt: Maximum moving time
   v: velocity of the walkers
   pos: dictionary with positions of each node in the DiGraph {node:(posx,pos_y)}
   n_walk: number of realisations of the random walk during Dt max time'''


def node_walkers(G, Dt, v, pos, n_walk):
    # calculating the time intrinsec to each edge
    delta_t = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
            for edge in G.edges}
    nx.set_edge_attributes(G, delta_t, name='dt')
    # upper bound for steps
    min_edge_time = min(nx.get_edge_attributes(G, 'dt').values())
    max_steps = int(Dt/min_edge_time)
    print('maximum steps allowed each realisation', max_steps)
    # overall_node/edge_weights account for all the edge/node passings summed for
    # n_walk iterations
    # overall_node_weights= {node:0 for node in list(G.nodes)}
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # path = []
        initial_nodes = sorted(list(G.nodes))
        # path stores the position of each walker in each time step
        # path.append(initial_nodes)
        # weights/edge_weights counts the amount of times a walker has visited
        # each node/edge in 1 realisation
        # weights = {i_node:path[0].count(i_node) for i_node in list(G.nodes)}
        edge_weights = {i_edge: 0 for i_edge in list(G.edges)}
        # first let's difine the time used for each walker
        time = np.zeros(len(initial_nodes))
        last_time = np.copy(time)  # just for checking
        # will store the position
        final_nodes = np.copy(np.array(initial_nodes))
        # of the walkers after 1 step
        # 1 step moves (or tries to move according to time left) all the walkers
        # through a neighboring edge
        for step in range(max_steps):
            # list of neighboring nodes of each walker
            neighbors = [{n: (G[walker][n]['dt'] if (walker, n) in list(G.edges)
                          else G[n][walker]['dt']) for n in
                          nx.all_neighbors(G, walker)} for walker in initial_nodes]
#            print(neighbors)
            # now randomly move the random walker to another node only if time+edge_time
            # is < Dt
            ind = 0
            # checks how many walkers have nowhere left to go
            ended_walkers = 0
            # print(step)
            for goto_nodes in neighbors:
                # time of each possible edge for a given node

                if np.any(time[ind] + np.array(list(goto_nodes.values())) <= Dt):
                    possible_nodes = [k for k in goto_nodes.keys() if
                                      (goto_nodes[k]+time[ind]) <= Dt]
                    sel = random.choice(possible_nodes)
                    final_nodes[ind] = sel
                else:
                    # print(time[ind])
                    ended_walkers += 1
                ind += 1
            time = [t+G[initial_nodes[j]][final_nodes[j]]['dt'] if
                    (initial_nodes[j], final_nodes[j]) in list(G.edges) else
                    t+G[final_nodes[j]][initial_nodes[j]]['dt'] if
                    (final_nodes[j], initial_nodes[j]) in list(G.edges) else t
                    for j, t in enumerate(time)]
            # print(time)
            # path.append(np.copy(final_nodes))
            # counting edge visits according to direction of edge
            for node_i, node_f in zip(initial_nodes, final_nodes):
                if node_i != node_f:
                    if (node_i, node_f) in list(G.edges):
                        edge_weights[(node_i, node_f)] += 1
                    elif (node_f, node_i) in list(G.edges):
                        edge_weights[(node_f, node_i)] -= 1
                    else:
                        print('ini different than final but edge not in graph')
            # count the occutpation of each node after the moves
            # for node_i, node_f in zip(initial_nodes, final_nodes):
            #     if node_i != node_f:
            #         weights[node_i] += 1

            initial_nodes = np.copy(final_nodes)
            if ended_walkers == len(initial_nodes):
                print((np.array(time) == last_time).all())
                print('all walkers finished, steps needed for this realisation '
                      + str(step)+'/'+str(max_steps))
                break
            last_time = np.copy(np.array(time))
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
    #     for node in G.nodes:
    #         overall_node_weights[node] += weights[node]
    #     #set node value as the number of visits of each value
    #     # print(weights)
    # nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)


# %% OPTIMIZED NODE_WALKERS (time/n_wak = 11s for clusters of 100 points)
'''Function that puts n = len(list(G.nodes)) random walkers in a Digraph
transitable in both directions and moves them during Dt time. The walkers always
move until they are not able to go anywhere else with the time left'''

'''G: nx.DiGraph
   Dt: Maximum moving time
   v: velocity of the walkers
   pos: dictionary with positions of each node in the DiGraph {node:(posx,pos_y)}
   n_walk: number of realisations of the random walk during Dt max time'''


def node_walkers(G, Dt, v, pos, n_walk):
    # calculating the time intrinsec to each edge
    delta_t = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
            for edge in G.edges}
    nx.set_edge_attributes(G, delta_t, name='dt')
    # upper bound for steps
    min_edge_time = min(delta_t.values())
    max_steps = int(Dt/min_edge_time)
    print('maximum steps allowed each realisation', max_steps)
    # overall_edge_weights account for all the edge passings summed for
    # n_walk iterations
    # MSD_dict = {i:[] for i in G.nodes}
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {i: i for i in G.nodes}

        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        edge_weights = {i_edge: 0 for i_edge in list(G.edges)}
        # first let's difine the time used for each walker
        time = {node_ind: 0 for node_ind in G.nodes}
        active_walkers = {node_ind: True for node_ind in G.nodes}
        # 1 step moves (or tries to move according to time left) all the walkers
        # through a neighboring edge
        for step in range(max_steps):
            # list of active walkers indices
            active_ls = [walk_ind for walk_ind in active_walkers.keys()
                         if active_walkers[walk_ind]]
            if len(active_ls) == 0:  # all walkers have finished
                print('all walkers finished, steps needed for this realisation '
                      + str(step)+'/'+str(max_steps))
                break

            # list of neighboring nodes of each walker
            neighbors = [{n: (G[pos_walkers[walker]][n]['dt'] if (pos_walkers[walker], n)
                              in list(G.edges) else G[n][pos_walkers[walker]]['dt'])
                          for n in nx.all_neighbors(G, pos_walkers[walker])}
                         for walker in active_ls]
#            print(neighbors)
            # now randomly move the random walker to another node only if time+edge_time
            # is < Dt
            for walker, goto_nodes in zip(active_ls, neighbors):
                # time of each possible edge for a given node
                if np.any(time[walker] + np.array(list(goto_nodes.values())) <= Dt):
                    # list of the possible final nodes
                    possible_nodes = [k for k in goto_nodes.keys() if
                                      (goto_nodes[k]+time[walker]) <= Dt]
                    # random selection between the final nodes
                    sel = random.choice(possible_nodes)
                    # updating edge flow dict
                    if (pos_walkers[walker], sel) in list(G.edges):
                        edge_weights[(pos_walkers[walker], sel)] += 1
                    elif (sel, pos_walkers[walker]) in list(G.edges):
                        edge_weights[(sel, pos_walkers[walker])] -= 1
                    else:
                        print('ini different than final but edge not in graph')

                    # updating walker position and time
                    pos_walkers[walker] = sel
                    time[walker] += goto_nodes[sel]
                    # MSD_dict[walker].append([time[walker],np.array(pos[sel])])
                else:
                    # desactivate the walker that has finished
                    active_walkers[walker] = False
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)  # MSD_dict)


# %%
min_dist = min(nx.get_edge_attributes(graph, 'dist').values())
max_dist = max(nx.get_edge_attributes(graph, 'dist').values())
v = max_dist/1  # min dist walked in Delta t = 1
Dt = max_dist/(v)
ini, final, path = node_walkers(graph, Dt, v)
# %%
# path1 = [graph.nodes[path[i][0]]['pos'] for i in range(len(path))]
# path2 = [graph.nodes[path[i][1]]['pos'] for i in range(len(path))]
x, y = zip(*ini)
u, v = zip(*final)

x = np.array(x)
y = np.array(y)
u = np.array(u)
v = np.array(v)
plt.figure()
u = (u-x)/np.sqrt(max_dist)
v = (v-y)/np.sqrt(max_dist)
plt.triplot(clusters[:, 0], clusters[:, 1], tri.simplices)
plt.plot(clust_1[:, 0], clust_1[:, 1], 'or', label='cluster 1')
plt.plot(clust_2[:, 0], clust_2[:, 1], 'ob', label='cluster 2')
plt.legend()
# plt.plot(np.array(path1)[:,0],np.array(path1)[:,1], 'r','-')
# plt.scatter(x,y, s = 100, color ='y')
# plt.scatter(u[0]+x[0],v[0]+y[0], s = 100, color = 'y')
# plt.plot([x[0], u[0]+x[0]], [y[0], v[0]+y[0]])

# nx.draw_networkx(graph, pointIDXY, with_labels=False, node_size = 10, ax = ax)
plt.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=1.)
# ax.quiver(np.array(path1)[0,0], np.array(path1)[0,1], np.array(path1)[len(path1)-1,0],
           # np.array(path1)[len(path1)-1,1], angles='xy', scale_units='xy', scale=1)


# %% WALKERS ON A DIRECTED GRAPH (walk during amount of steps)

'''
moves "n_walk" walkers/node in a directed graph "G" "steps" times. The initial positions of the
walkers are randomly selected amongst the nodes. The function returns the graph
with the amount of times any walker has passed through each node as attributes
to the nodes

'''


def digraph_walkers(G, steps, n_walk):

    # let's add a list of n_walk random walkers
    overall_node_weights = {node: 0 for node in list(G.nodes)}
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        path = []
        initial_nodes = list(G.nodes)
        # initial_nodes += initial_nodes
        # path stores the position of each walker in each time step
        path.append(initial_nodes)
        # weights counts the amount of times a walker has visited each node
        weights = {i_node: path[0].count(i_node) for i_node in list(G.nodes)}
        edge_weights = {i_edge: 0 for i_edge in list(G.edges)}
        # now, let's move the walker to a connected node
        # first let's see the available nodes to go to
        for step in range(steps):
            # list of neighboring nodes (successors) of each walker
            neighbors = [[n for n in (list(G.neighbors(walker)) +
                                      list(G.predecessors(walker)))]
                         for walker in initial_nodes]
            # now move the random walker to another node
            final_nodes = [random.choice(goto_nodes)
                                         for goto_nodes in neighbors]

            path.append(final_nodes)

            # counting edge visits according to direction of edge
            for node_i, node_f in zip(initial_nodes, final_nodes):
                if node_i != node_f:
                    if (node_i, node_f) in list(G.edges):
                        edge_weights[(node_i, node_f)] += 1
                    else:
                        edge_weights[(node_f, node_i)] -= 1
            # count the occutpation of each node after the moves
            for node_i, node_f in zip(initial_nodes, final_nodes):
                if node_i != node_f:
                    weights[node_i] += 1

            initial_nodes = final_nodes
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
        for node in G.nodes:
            overall_node_weights[node] += weights[node]
        # set node value as the number of visits of each value
        # print(weights)
    nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)

# %%


'''Hodge decomposition of a tree directed graph'''

# creating the graph
G = nx.balanced_tree(2, 6)
pos = EoN.hierarchy_pos(G, 0, width=1)  # gets the positions for the
# graph to be seen as a tree. only works in non directed graps. Get pos on non directed
G2 = nx.Graph(G)
tree = nx.DiGraph(G)
out_edges = [edge for edge in tree.edges if edge[0]
    < edge[1]]  # removing all outward edges
tree.remove_edges_from(out_edges)
nx.draw_networkx(tree, pos=pos, with_labels=True)

# %% Walking the graph
steps = 10
n_walk = 1
tree_w = digraph_walkers(tree, steps, n_walk)
edge_labels = nx.get_edge_attributes(tree_w, 'edge_visits')
nx.draw_networkx(tree, pos=pos, with_labels=True)
nx.draw_networkx_edge_labels(tree_w, pos, edge_labels=edge_labels,
                             rotate=False, font_size=8)

# %%

''' DIVERGENCE'''

div = {node: 0 for node in tree.nodes}

for node in tree_w.nodes:
    for in_edge in tree_w.in_edges(node):
        div[node] += tree_w[in_edge[0]][in_edge[1]]['edge_visits']
    for out_edge in tree_w.out_edges(node):
        div[node] -= tree_w[out_edge[0]][out_edge[1]]['edge_visits']

nx.set_node_attributes(tree, div, name='div')

# plotting node div with color gradient

color_p = list(div.values())
colors = range(min(color_p), max(color_p))
fig, ax = plt.subplots()
cmap = plt.cm.Blues
vmin = min(colors)
vmax = max(colors)

nx.draw_networkx(tree_w, pos=pos, with_labels=False, node_color=color_p,
                 cmap=cmap, vmin=vmin, vmax=vmax, node_size=100)

nx.draw_networkx_labels(tree_w, pos, labels=div)

ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
plt.colorbar(sm)
plt.show()


# %%Original triangle finder function


def countTriangle(G, isDirected):

    g = nx.adjacency_matrix(G).toarray()
    nodes = len(g)
    count_Triangle = 0
    triangles = {}
   	# Consider every possible
   	# triplet of edges in graph
    for i in range(nodes):
        for j in range(nodes):
            for k in range(nodes):

				# check the triplet
				# if it satisfies the condition
                if (i != j and i != k and j != k and ((i, j) in G.edges or (j, i)
                    in G.edges) and ((j, k) in G.edges or (k, j) in G.edges)
                    and ((i, k) in G.edges or (k, i) in G.edges)):

                    if (i, j) in G.edges:
                        i1 = (i, j)
                    elif (j, i) in G.edges:
                        i1 = (j, i)

                    if (j, k) in G.edges:
                        i2 = (j, k)
                    elif (k, j) in G.edges:
                        i2 = (k, j)

                    if (i, k) in G.edges:
                        i3 = (i, k)
                    elif (k, i) in G.edges:
                        i3 = (k, i)
                        # print(list(permutations([i1, i2, i3])))
                    for tri in permutations([i1, i2, i3]):
                        if list(tri) in triangles.values():
                            rep = True
                            break
                        else:
                            rep = False
                    if rep == False:
                        x = list(set(i1+i2+i3))
                        triangles[tuple(sorted(x))]\
                            = [i1, i2, i3]
                        count_Triangle += 1
	# If graph is directed , division is done by 3
	# else division by 6 is done
    if isDirected:
        return count_Triangle, triangles
    else:
        print('Error: Graph is not directed')

# %% ALTERNATIVE WAY TO FND TRIANGLES (faster)


def find_triangles(G):
    """
    Returns a list of all triangles in the graph G.
    """
    triangles = {}
    count_Triangle = 0
    for node in G.nodes():
        # Get the neighbors of the node
        neighbors = set(list(G.neighbors(node)) + list(G.predecessors(node)))
        for neighbor in neighbors:
            # Get the common neighbors of the node and its neighbor
            common_neighbors = neighbors.intersection(set(list(G.neighbors(neighbor))
                                                          + list(G.predecessors(neighbor))))
            for common_neighbor in common_neighbors:
                # Add the triangle to the list
                if node < neighbor and neighbor < common_neighbor:
                    x = [node, neighbor, common_neighbor]
                    if (node, neighbor) in G.edges:
                        edge_1 = [node, neighbor]
                    else:
                        edge_1 = [neighbor, node]
                    if (node, common_neighbor) in G.edges:
                        edge_2 = [node, common_neighbor]
                    else:
                        edge_2 = [common_neighbor, node]

                    if (neighbor, common_neighbor) in G.edges:
                        edge_3 = [neighbor, common_neighbor]
                    else:
                        edge_3 = [common_neighbor, neighbor]

                    triangles[tuple(x)] = [edge_1, edge_2, edge_3]
                    count_Triangle += 1
    return count_Triangle, triangles

# %% HODGE DECOMPOSITION FUNCTION


def hodge_decomposition(G, attr_name):

# vector of the edge attribute considered where each component is the value of the
# attr of a given edge
    g_field = np.transpose(np.array([G[edge[0]][edge[1]][attr_name] for edge
                                     in G.edges]))

# Computing divergence
    div = {node: 0 for node in G.nodes}

    for node in G.nodes:
        for in_edge in G.in_edges(node):
            div[node] -= G[in_edge[0]][in_edge[1]][attr_name]
        for out_edge in G.out_edges(node):
            div[node] += G[out_edge[0]][out_edge[1]][attr_name]
    nx.set_node_attributes(G, div, name='div')

# GRADIENT OPERATOR
    grad_arr = []
    for edge in G.edges:
        row = []
        for node in G.nodes:
            if edge[1] == node:
                row.append(1)
            elif edge[0] == node:
                row.append(-1)
            else:
                row.append(0)
        grad_arr.append(row)
        row = []
    grad_arr = np.array(grad_arr)
    # print('gradient op\n',grad_arr)
# LAPLACIAN MATRIX
    lap = np.transpose(grad_arr).dot(grad_arr)
    # print('checking solvability',np.linalg.cond(lap))
    # OBTAINING THE GRADIENT COMPONENT
    # compute the pseudo inverse of the laplacian to solve the syst of equations
    # computationally intensive for large problems
    Apinv = np.linalg.pinv(lap)
    # apply the pseudo-inverse to the rhs vector to obtain the `solution'
    div_arr = np.transpose(np.array(list(div.values())))
    # print('divergence\n',div_arr)

    pot_field = np.squeeze(np.asarray(Apinv.dot(-div_arr)))
    # print('error in node potentials', lap.dot(pot_field)+div_arr)
    # gradient of the potential field
    # pot_field.insert(remove_index, 0)
    pot_nodes = {n: pot for n, pot in enumerate(pot_field)}
    # print('node potentials', pot_nodes)
    grad_comp_arr = np.transpose(grad_arr.dot(pot_field))
    grad_comp = {edge: grad_comp_arr[i] for i, edge in enumerate(G.edges)}
# SOLENOIDAL COMPONENT
# Calculating the edge-2 matrix or oriented-face incidence matrix which
# corresponds to the curl operator
    n_tri, tri_dict = find_triangles(G)
    print('number of triangles', n_tri)
    if n_tri != 0:
        curl_op = []
        # creating the curl operator:
        for triangle in tri_dict.keys():
            row = []
            for edge in G.edges:
#                print(triangle)
                if (triangle[1] == min(edge) and triangle[2] == max(edge)) or\
                (triangle[0] == min(edge) and triangle[1] == max(edge)):
                    if edge[0] < edge[1]:
                        row.append(1)
                    else:
                        row.append(-1)
                elif triangle[0] == min(edge) and triangle[2] == max(edge):
                    if edge[0] < edge[1]:
                        row.append(-1)
                    else:
                        row.append(1)
                else:
                    row.append(0)

#            print(row, tri_dict[triangle], list(G.edges))
            curl_op.append(row)

    # Caclulating the delta1 delta1* matrix (rot_arr) (using levi Chivita tensor
    # a la Hodge Tutorial)
        # rot_arr = []
        # for row_el in tri_dict.keys():
        #     row = []
        #     for col_el in tri_dict.keys():
        #         if col_el == row_el:
        #             row.append(3)
        #         elif len(set(row_el).symmetric_difference(col_el)) == 2:
        #             # returns all items in both sets but not the ones that are in both sets
        #             diff = list(set(row_el).symmetric_difference(col_el))
        #             c_mod = [diff[0] if x==diff[1] else diff[1] if x==diff[0]
        #                      else x for x in row_el]
        #             c = sorted(c_mod)

        #             row.append(int(LeviCivita(c_mod.index(c[0]),c_mod.index(c[1])
        #                                       ,c_mod.index(c[2]))))
        #         elif len(set(row_el).symmetric_difference(col_el)) > 2:
        #             row.append(0)
        #     rot_arr.append(row)
        # rot_arr = np.array(rot_arr)
        curl_op = np.array(curl_op)
        # print(rot_arr)
        rot_arr_inv = np.linalg.pinv(curl_op.dot(
            np.transpose(curl_op))).dot(curl_op)

#        print(curl_op)
        # computing the curl of the graph
        rot = curl_op.dot(g_field)
        # print('big rot arr\n',curl_op.dot(np.transpose(curl_op)))
        # print('curl_op',curl_op)
        # print('field',g_field)
        # print('solvability of curl comp', np.linalg.cond(rot_arr))
        # rot_pinv = np.linalg.pinv(rot_arr)

        # solving the system of equations
        tri_pot = rot_arr_inv.dot(g_field)
        # tri_pot = np.squeeze(np.asarray(rot_pinv.dot(rot)))
        # print('error curl component',
        #       curl_op.dot(np.transpose(curl_op)).dot(tri_pot)-rot)
        # solenoidal component is delta1* (transpose of curl op) of the potential:
        g_s = np.transpose(curl_op).dot(tri_pot)
        print('curl of graph',
              rot)
        sol_comp = {edge: comp for edge, comp in zip(G.edges, g_s)}
        # for edge, comp in zip(G.edges, g_s):
        #     sol_comp[edge] = comp
    else:
        g_s = np.transpose(np.zeros(np.shape(grad_comp_arr)))
        sol_comp = {edge: 0 for edge in G.edges}
        # for edge in G.edges:
        #     sol_comp[edge] = 0

# HARMONIC COMPONENT

    if n_tri != 0:
        big_arr = grad_arr.dot(np.transpose(grad_arr)) +\
            np.transpose(curl_op).dot(curl_op)
    else:
        big_arr = grad_arr.dot(np.transpose(grad_arr))
#        print(big_arr)

    g_har_vec = g_field - grad_comp_arr - g_s
    thrs = 10**(-10)
    if np.all(np.abs(big_arr.dot(g_har_vec)) < thrs):

        har_comp = {edge: har for edge, har in zip(G.edges, g_har_vec)}
        # for edge, har in zip(G.edges, g_har_vec):
        #     har_comp[edge] = har
    else:
        print('problem in the harmonic component')
        print('error of the harmonic component', big_arr.dot(g_har_vec))
        # print('divergence of the harmonic component',
        #       grad_arr.dot(np.transpose(grad_arr)).dot(g_har_vec))
        # print('curl of the harmonic component',
        #       curl_op.dot(g_har_vec))

        har_comp = {edge: har for edge, har in zip(G.edges, g_har_vec)}
        # for edge, har in zip(G.edges, g_har_vec):
        #     har_comp[edge] = har
    return grad_comp, sol_comp, har_comp, pot_nodes, div
# %% CUSTOM GRAPH


g = nx.DiGraph()

# SQUARE
# g.add_nodes_from([0, 1, 2, 3])
# g.add_edges_from([(0,1), (1,2), (2,3), (3,0)])
# attr = {(0,1): 2, (1,2): 2, (2,3): 2, (3,0):-2}

# 2 triangles
g.add_nodes_from([0, 1, 2, 3, 4])
g.add_edges_from([(0, 1), (2, 0), (1, 2), (1, 3), (2, 3), (3, 4), (0, 4)])
attr = {(0, 1): 3, (2, 0): 3, (1, 2): 2, (1, 3): 1, (2, 3): 1, (3, 4): -2, (0, 4): 1}
# line
# g.add_nodes_from([0, 1, 2, 3])
# g.add_edges_from([(0,1), (1,2), (2,3)])
# attr = {(0,1): 2,(1,2): 2, (2,3):2}


nx.set_edge_attributes(g, attr, name='edge_visits')
grad_comp, sol_comp, har_comp, div = hodge_decomposition(g, 'edge_visits')
grad_comp = {comp: round(grad_comp[comp], 4) for comp in grad_comp.keys()}
sol_comp = {comp: round(sol_comp[comp], 4) for comp in sol_comp.keys()}
har_comp = {comp: round(har_comp[comp], 4) for comp in har_comp.keys()}
# pos = EoN.hierarchy_pos(G, 0, width = 1)
# g = tree_w
pos = nx.planar_layout(g)
# pos = nx.spring_layout(g, k = 0.7)
# pos = {0: (0,0), 1: (1,0), 2: (1,1), 3: (0,1)}

nx.draw_networkx(g, pos=pos, with_labels=True, node_size=250, font_size=15,
                 node_color='#AED0EE')
nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=attr,
                             rotate=False, font_size=15)
plt.axis('off')
# %%

plt.subplots(2, 2, figsize=(10, 10))
plt.subplot(221)
plt.title('original')
nx.draw_networkx(g, pos=pos, with_labels=True, node_size=250, font_size=15,
                 node_color='#AED0EE')
nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=attr,
                             rotate=False, font_size=15)

plt.subplot(222)
plt.title('gradient')
nx.draw_networkx(g, pos=pos, with_labels=True, node_size=250, font_size=15,
                 node_color='#AED0EE')
nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=grad_comp,
                             rotate=False, font_size=15)

plt.subplot(223)
plt.title('solenoidal')
nx.draw_networkx(g, pos=pos, with_labels=True, node_size=250, font_size=15,
                 node_color='#AED0EE')
nx.draw_networkx_edge_labels(g, pos, edge_labels=sol_comp,
                             rotate=False, font_size=15)

plt.subplot(224)
plt.title('harmonic')
nx.draw_networkx(g, pos=pos, with_labels=True, node_size=250, font_size=15,
                 node_color='#AED0EE')
nx.draw_networkx_edge_labels(g, pos=pos, edge_labels=har_comp,
                             rotate=False, font_size=15)
plt.tight_layout()

# %% WHEEL

''' CIRCULAR GRAPH: 3 concentric wheels'''
# number of nodes
n_nodes = 60
# wheel graph
# wheel = nx.circular_ladder_graph(n_nodes)
# position of the nodes
wheel = nx.Graph()
wheel.add_nodes_from(list(range(n_nodes)))
wheel.add_edges_from([(i, i+1) if i <= (n_nodes/3-2) else (0, i)
                     for i in range(int(n_nodes/3))])
wheel.add_edges_from([(i, i+1) if i <= (2*n_nodes/3-2) else (n_nodes/3, i) for i
                      in range(int(n_nodes/3), 2*int(n_nodes/3))])
wheel.add_edges_from([(i, i+1) if i <= (n_nodes-2) else (int(2*n_nodes/3), i) for i in
                      range(2*int(n_nodes/3), int(n_nodes))])
wheel.add_edges_from([(i, i+20) for i in range(int(2*n_nodes/3))])
R = 1
R2 = 1.5
R3 = 2
dalph = 2*np.pi/(n_nodes/3)

pos_c = {}
for i in range(int(n_nodes/3)):
    pos_c[i] = (R*np.cos(dalph*i), R*np.sin(dalph*i))
    pos_c[i+n_nodes/3] = (R2*np.cos(dalph*i), R2*np.sin(dalph*i))
    pos_c[i+2*n_nodes/3] = (R3*np.cos(dalph*i), R3*np.sin(dalph*i))
# nx.draw_networkx(wheel, with_labels = True, pos = pos_c)

# transforming the wheel graph into a digraph and removing ill posing edges
wheel_d = nx.DiGraph(wheel)
out_edges = [edge for edge in wheel_d.edges if edge[0]
    > edge[1]]  # removing all outward edges
wheel_d.remove_edges_from(out_edges)

plt.figure(figsize=(8, 8))
nx.draw_networkx(wheel_d, pos=pos_c, with_labels=True)

Dt = 10
n_walk = 10
v = 1

walk_wh = node_walkers(wheel_d, Dt, v, pos_c, n_walk)
# walk_wh = digraph_walkers(wheel_d, steps, n_walk)

plt.figure(figsize=(8, 8))
nx.draw_networkx(walk_wh, pos=pos_c, with_labels=True)
nx.draw_networkx_edge_labels(walk_wh, pos=pos_c, edge_labels=nx.get_edge_attributes(walk_wh, 'edge_visits'),
                             rotate=False, font_size=8)
plt.tight_layout()
plt.show()
# %%
grad_wh, sol_wh, har_wh, pot_wh, div_wh = hodge_decomposition(
    walk_wh, 'edge_visits')
# %%

grad_wh = {comp: round(grad_wh[comp], 2) for comp in grad_wh.keys()}
sol_wh = {comp: round(sol_wh[comp], 2) for comp in sol_wh.keys()}
har_wh = {comp: round(har_wh[comp], 2) for comp in har_wh.keys()}

edge_wh = nx.get_edge_attributes(walk_wh, 'edge_visits')
plt.subplots(2, 2, figsize=(10, 15))
plt.subplot(221)
plt.title('original')
nx.draw_networkx(wheel_d, pos=pos_c, with_labels=False, node_size=10)
nx.draw_networkx_edge_labels(wheel_d, pos=pos_c, edge_labels=edge_wh,
                             rotate=False, font_size=8)

plt.subplot(222)
plt.title('gradient')
nx.draw_networkx(wheel_d, pos=pos_c, with_labels=False, node_size=10)
nx.draw_networkx_edge_labels(wheel_d, pos=pos_c, edge_labels=grad_wh,
                             rotate=False, font_size=8)

plt.subplot(223)
plt.title('solenoidal')
nx.draw_networkx(wheel_d, pos=pos_c, with_labels=False, node_size=10)
nx.draw_networkx_edge_labels(wheel_d, pos_c, edge_labels=sol_wh,
                             rotate=False, font_size=8)

plt.subplot(224)
plt.title('harmonic')
nx.draw_networkx(wheel_d, pos=pos_c, with_labels=False, node_size=10)
nx.draw_networkx_edge_labels(wheel_d, pos=pos_c, edge_labels=har_wh,
                             rotate=False, font_size=8)
plt.tight_layout()
plt.show()

# %%
# %%IMPORTANCE OF EACH COMPONENT

w = np.array(list(edge_wh.values()))
wg = np.array(list(grad_wh.values()))
ws = np.array(list(sol_wh.values()))
wh = np.array(list(har_wh.values()))
weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)

# %% 2 CLUSTER DELAUNAY
'''Delaunay'''

'''
we will create 2 clusters of points normally distributed
'''
# x,y coords of points
np.random.seed(1000)
nodes_cl1 = 200
nodes_cl2 = 200
clust_1 = np.random.normal(0, 1, size=(nodes_cl1, 2))
clust_2 = np.random.normal(5, 1, size=(nodes_cl2, 2))

# clust_1 = np.array([[0,0], [0,1], [1,0]])
# clust_2 = clust_1+2

# clust_1 = np.array([[0,0], [0.5,0.5], [1,0], [0.5,-0.5]])
# clust_2 = clust_1+2

clusters = np.concatenate((clust_1, clust_2), axis=0)
tri = Delaunay(clusters)
plt.figure()
plt.triplot(clusters[:, 0], clusters[:, 1], tri.simplices)
plt.plot(clust_1[:, 0], clust_1[:, 1], 'or', label='cluster 1')
plt.plot(clust_2[:, 0], clust_2[:, 1], 'ob', label='cluster 2')
plt.legend()
plt.show()

# create a set for edges that are indexes of the points
edges = set()
# for each Delaunay triangle
for n in range(tri.nsimplex):
    # for each edge of the triangle
    # sort the vertices
    # (sorting avoids duplicated edges being added to the set)
    # and add to the edges set
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 1]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 1], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))


# make a graph based on the Delaunay triangulation edges
ini_del = nx.Graph(list(edges))
# plt.figure()
# plot graph

# dictionary of node:position
pointIDXY = dict(zip(range(len(clusters)), clusters))
nx.set_node_attributes(ini_del, pointIDXY, name='pos')
# nx.draw(graph, pointIDXY)
# plt.show()


# DELAUNAY
del_dg = nx.DiGraph(ini_del)
# changing problematic edges (edge[0]> edge[1] changed to edge[0] < edge[1])
out_edges = [edge for edge in del_dg.edges if edge[1]
    > edge[0]]  # removing problematic edges
del_dg.remove_edges_from(out_edges)
# for edge in out_edges:
#     if (edge[1],edge[0]) not in del_dg.edges:
#         del_dg.add_edge(edge[1],edge[0])
# nx.draw_networkx(del_dg, pos = pointIDXY, with_labels = True)
Dt = 10
v = 1
n_walk = 10
walk_del = node_walkers(del_dg, Dt, v, pointIDXY, n_walk)


# %%
grad_del, sol_del, har_del, pot_del, div_del = hodge_decomposition(
    walk_del, 'edge_visits')
# %%
pos_c = pointIDXY

for node in walk_del.nodes:
    for in_edge in walk_del.in_edges(node):
        div[node] -= np.array([grad_del[in_edge],
                              float(walk_del[in_edge[0]][in_edge[1]]['edge_visits'])])
    for out_edge in walk_del.out_edges(node):
        div[node] += np.array([grad_del[out_edge],
                              float(walk_del[out_edge[0]][out_edge[1]]['edge_visits'])])
print('divergence of graph vs divergence of gradient component', div)

grad_del = {comp: round(grad_del[comp], 2) for comp in grad_del.keys()}
sol_del = {comp: round(sol_del[comp], 2) for comp in sol_del.keys()}
har_del = {comp: round(har_del[comp], 2) for comp in har_del.keys()}

edge_del = nx.get_edge_attributes(walk_del, 'edge_visits')


plt.subplots(2, 2, figsize=(15, 15))
plt.subplot(221)
plt.title('Original Graph')

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_div = list(div_del.values())
colors_div = range(int(min(color_div)), int(max(color_div)))
cmap_div = plt.cm.seismic
vmin_div = -max(colors_div)
vmax_div = max(colors_div)

nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
                       # node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                       # vmax=vmax_div, node_size = 10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_p,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega\right|$')

sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                              vmax=vmax_div))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node divergence')

plt.subplot(222)
color_g = np.abs(np.array(list(grad_del.values())))
# plotting edges with color gradient

color_pot = list(pot_del.values())
colors_pot = range(int(min(color_pot)), int(max(color_pot)))
cmap_pot = plt.cm.seismic
vmin_pot = min(colors_pot)
vmax_pot = max(colors_pot)

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)
plt.title('Gradient component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None,
                        node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_g,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_g\right|$')

sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                              vmax=vmax_pot))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node potentials')

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_s = np.abs(np.array(list(sol_del.values())))
plt.subplot(223)
plt.title('Solenoidal Component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_s,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
# nx.draw_networkx(walk_del,pos = pos_c, with_labels=False, node_size = 10,
#                  edge_color=color_s, edge_cmap=cmap, vmin=vmin, vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_s\right|$')


color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_h = np.abs(np.array(list(har_del.values())))
plt.subplot(224)
plt.title('Harmonic Component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_h,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_h\right|$')

plt.tight_layout()
plt.show()

# %%IMPORTANCE OF EACH COMPONENT

w = np.array(list(edge_del.values()))
wg = np.array(list(grad_del.values()))
ws = np.array(list(sol_del.values()))
wh = np.array(list(har_del.values()))
weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)
# %% FRUCHT GRAPH TO SEE HARMONIC SOLENOIDAL AND GRADIENT COMPONENTS AT THE SAME TIME

# frucht graph
frucht = nx.frucht_graph()
# position of the nodes

# transforming the frucht graph into a digraph
frucht_d = nx.DiGraph(frucht)
out_edges = [edge for edge in frucht_d.edges if edge[1]
    > edge[0]]  # removing problematic edges
frucht_d.remove_edges_from(out_edges)

pos_c = nx.planar_layout(frucht_d)
steps = 10
n_walk = 100
walk_fr = digraph_walkers(frucht_d, steps, n_walk)

plt.figure(figsize=(8, 8))
nx.draw_networkx(walk_fr, pos=pos_c, with_labels=True)
nx.draw_networkx_edge_labels(walk_fr, pos=pos_c, edge_labels=nx.get_edge_attributes(walk_fr, 'edge_visits'),
                             rotate=False, font_size=8)
plt.tight_layout()
plt.show()

grad_fr, sol_fr, har_fr, div = hodge_decomposition(walk_fr, 'edge_visits')


grad_fr = {comp: round(grad_fr[comp], 2) for comp in grad_fr.keys()}
sol_fr = {comp: round(sol_fr[comp], 2) for comp in sol_fr.keys()}
har_fr = {comp: round(har_fr[comp], 2) for comp in har_fr.keys()}


edge_fr = nx.get_edge_attributes(walk_fr, 'edge_visits')
plt.subplots(2, 2, figsize=(10, 15))
plt.subplot(221)
plt.title('original')
nx.draw_networkx(walk_fr, pos=pos_c, with_labels=False, node_size=10)
nx.draw_networkx_edge_labels(walk_fr, pos=pos_c, edge_labels=edge_fr,
                             rotate=False, font_size=8)

# plotting edges with color gradient
color_p = np.abs(np.array(list(edge_fr.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Blues
vmin = min(colors)
vmax = max(colors)


color_g = np.abs(np.array(list(grad_fr.values())))

plt.subplot(222)
plt.title('gradient')
nx.draw_networkx(walk_fr, pos=pos_c, with_labels=False, node_size=10,
                 edge_color=color_g, edge_cmap=cmap, vmin=vmin, vmax=vmax)
nx.draw_networkx_edge_labels(walk_fr, pos=pos_c, edge_labels=grad_fr,
                             rotate=False, font_size=8)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
plt.colorbar(sm)

color_s = np.abs(np.array(list(sol_fr.values())))

plt.subplot(223)
plt.title('solenoidal')
nx.draw_networkx(walk_fr, pos=pos_c, with_labels=False, node_size=10,
                 edge_color=color_s, edge_cmap=cmap, vmin=vmin, vmax=vmax)
nx.draw_networkx_edge_labels(walk_fr, pos_c, edge_labels=sol_fr,
                             rotate=False, font_size=8)
plt.colorbar(sm)

color_h = np.abs(np.array(list(har_fr.values())))
plt.subplot(224)
plt.title('harmonic')
nx.draw_networkx(walk_fr, pos=pos_c, with_labels=False, node_size=10,
                 edge_color=color_h, edge_cmap=cmap, vmin=vmin, vmax=vmax)
nx.draw_networkx_edge_labels(walk_fr, pos=pos_c, edge_labels=har_fr,
                             rotate=False, font_size=8)
plt.colorbar(sm)

plt.tight_layout()
plt.show()
# %% CIRCULATING DELAUNAYS
'''CIRCULATING DELAUNAYS'''


'''
we will create 4 clusters of points normally distributed
'''
# x,y coords of points
np.random.seed(1000)
nodes_cl = 100
clust_1 = np.random.normal((2, 0), 0.6, size=(nodes_cl, 2))
clust_2 = np.random.normal((-2, 0), 0.6, size=(nodes_cl, 2))
clust_3 = np.random.normal((0, -2), 0.6, size=(nodes_cl, 2))
clust_4 = np.random.normal((0, 2), 0.6, size=(nodes_cl, 2))


clusters = np.concatenate((clust_1, clust_2, clust_3, clust_4), axis=0)
tri = Delaunay(clusters)
plt.figure()
plt.triplot(clusters[:, 0], clusters[:, 1], tri.simplices)
plt.plot(clust_1[:, 0], clust_1[:, 1], 'or', label='cluster 1')
plt.plot(clust_2[:, 0], clust_2[:, 1], 'ob', label='cluster 2')
plt.plot(clust_3[:, 0], clust_3[:, 1], 'og', label='cluster 3')
plt.plot(clust_4[:, 0], clust_4[:, 1], 'oy', label='cluster 4')
plt.legend()
plt.show()

# create a set for edges that are indexes of the points
edges = set()
# for each Delaunay triangle
for n in range(tri.nsimplex):
    # for each edge of the triangle
    # sort the vertices
    # (sorting avoids duplicated edges being added to the set)
    # and add to the edges set
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 1]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 0], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))
    edge = sorted([tri.vertices[n, 1], tri.vertices[n, 2]])
    edges.add((edge[0], edge[1]))


# make a graph based on the Delaunay triangulation edges
ini_del = nx.Graph(list(edges))
# plt.figure()
# plot graph

# dictionary of node:position
pointIDXY = dict(zip(range(len(clusters)), clusters))
nx.set_node_attributes(ini_del, pointIDXY, name='pos')
# nx.draw(graph, pointIDXY)
# plt.show()


# CIRCULATING DELAUNAY
del_dg = nx.DiGraph(ini_del)

out_edges = [edge for edge in del_dg.edges if edge[1]
    > edge[0]]  # removing problematic edges
del_dg.remove_edges_from(out_edges)

remove = []
for edge in del_dg.edges:
    pos_1 = pointIDXY[edge[0]]
    pos_2 = pointIDXY[edge[1]]
    vec = pos_2-pos_1
    if (pos_1 in clust_1 and pos_2 in clust_1) or (pos_1 in clust_2 and pos_2
    in clust_2):
        if abs(vec[0])/abs(vec[1]) > 1.3:
            remove.append(edge)
    elif (pos_1 in clust_3 and pos_2 in clust_3) or (pos_1 in clust_4 and pos_2
    in clust_4):
        if abs(vec[0])/abs(vec[1]) < 0.7:
            remove.append(edge)
    elif (pos_1 in clust_1 and pos_2 in clust_2) or (pos_1 in clust_3 and pos_2
    in clust_4) or (pos_1 in clust_2 and pos_2 in clust_1) or (pos_1 in clust_4
                                                                and pos_2
                                                                in clust_3):
        remove.append(edge)
del_dg.remove_edges_from(remove)
plt.figure()
nx.draw_networkx(del_dg, pos=pointIDXY, with_labels=False, node_size=10)
# for node in del_dg.nodes:
#     print(node, 'neighbours',list(del_dg.predecessors(node))+list(del_dg.successors(node)))
#     if len(list(del_dg.predecessors(node))+list(del_dg.successors(node))) == 0:
#         print('no neigbours')
Dt = 100
v = 0.1
n_walk = 1
start_time = time.time()
walk_del = node_walkers(del_dg, Dt, v, pointIDXY, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
# MSD
# times_all = []
# MSD_all = []
# for walker in tracking_del.keys():
#     times = []
#     MSD = []
#     for step in tracking_del[walker]:
#         times.append(step[0])
#         MSD.append(np.linalg.norm(step[1]-pointIDXY[walker])**2)
#     times_all.append(times)
#     MSD_all.append(MSD)

# plt.figure()
# plt.xlabel('Time')
# plt.ylabel('Squared Displacement')
# for i, walker in enumerate(times_all):
#     plt.plot(walker, MSD_all[i], ls = '-',color = (0, 200/255, 202/255, 0.3))
# %%
grad_del, sol_del, har_del, pot_del, div_del = hodge_decomposition(
    walk_del, 'edge_visits')
# %%
pos_c = pointIDXY


grad_del = {comp: round(grad_del[comp], 2) for comp in grad_del.keys()}
sol_del = {comp: round(sol_del[comp], 2) for comp in sol_del.keys()}
har_del = {comp: round(har_del[comp], 2) for comp in har_del.keys()}

edge_del = nx.get_edge_attributes(walk_del, 'edge_visits')


plt.subplots(2, 2, figsize=(15, 15))
plt.subplot(221)
plt.title('Original Graph')

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_div = list(div_del.values())
colors_div = range(int(min(color_div)), int(max(color_div)))
cmap_div = plt.cm.seismic
vmin_div = -max(colors_div)
vmax_div = max(colors_div)

nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
                       # node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                       # vmax=vmax_div, node_size = 10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_p,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega\right|$')

sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                              vmax=vmax_div))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node divergence')

plt.subplot(222)
color_g = np.abs(np.array(list(grad_del.values())))
# plotting edges with color gradient

color_pot = list(pot_del.values())
colors_pot = range(int(min(color_pot)), int(max(color_pot)))
cmap_pot = plt.cm.seismic
vmin_pot = min(colors_pot)
vmax_pot = max(colors_pot)

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)
plt.title('Gradient component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None,
                        node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_g,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_g\right|$')

sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                              vmax=vmax_pot))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node potentials')

color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_s = np.abs(np.array(list(sol_del.values())))
plt.subplot(223)
plt.title('Solenoidal Component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_s,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
# nx.draw_networkx(walk_del,pos = pos_c, with_labels=False, node_size = 10,
#                  edge_color=color_s, edge_cmap=cmap, vmin=vmin, vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_s\right|$')


color_p = np.abs(np.array(list(edge_del.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_h = np.abs(np.array(list(har_del.values())))
plt.subplot(224)
plt.title('Harmonic Component')
nx.draw_networkx_nodes(walk_del, pos=pos_c, label=None, node_size=10)
nx.draw_networkx_edges(walk_del, pos=pos_c, label=None, edge_color=color_h,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_h\right|$')

plt.tight_layout()
plt.show()
# %%


def plot_hodge(walk_graph, grad_comp, sol_comp, har_comp, pot, div, pos):

    edge_graph = nx.get_edge_attributes(walk_graph, 'edge_visits')

    w = np.array(list(edge_graph.values()))
    wg = np.array(list(grad_comp.values()))
    ws = np.array(list(sol_comp.values()))
    wh = np.array(list(har_comp.values()))
    weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
    weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
    weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

    print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)

    plt.subplots(2, 2, figsize=(15, 15))
    plt.subplot(221)
    plt.title('Original Graph 100%')

    color_p = np.abs(np.array(list(edge_graph.values())))
    colors = np.linspace(0, np.max(color_p))
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_div = list(div.values())
    colors_div = range(int(min(color_div)), int(max(color_div)))
    cmap_div = plt.cm.seismic
    vmin_div = -max(colors_div)
    vmax_div = max(colors_div)

    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=1,
                            node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                            vmax=vmax_div)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_p,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega\right|$')

    sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                                  vmax=vmax_div))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Node divergence')

    plt.subplot(222)
    color_g = np.abs(np.array(list(grad_comp.values())))
    # plotting edges with color gradient

    color_pot = list(pot.values())
    colors_pot = range(int(min(color_pot)), int(max(color_pot)))
    cmap_pot = plt.cm.seismic
    vmin_pot = min(colors_pot)
    vmax_pot = max(colors_pot)

    color_p = np.abs(np.array(list(edge_graph.values())))
    colors = np.linspace(0, np.max(color_p))
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)
    plt.title('Gradient component ' + str(round(weight_g*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None,
                            node_size=1, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_g\right|$')

    sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                                  vmax=vmax_pot))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Node potentials')

    color_p = np.abs(np.array(list(edge_graph.values())))
    colors = np.linspace(0, np.max(color_p))
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_s = np.abs(np.array(list(sol_comp.values())))
    plt.subplot(223)
    plt.title('Solenoidal Component ' + str(round(weight_s*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=1,
                           node_color=np.array([[0, 88/255, 173/355, 0.5]]))
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_s,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
    # nx.draw_networkx(walk_graph,pos = pos, with_labels=False, node_size = 5,
    #                  edge_color=color_s, edge_cmap=cmap, vmin=vmin, vmax=vmax)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_s\right|$')

    color_p = np.abs(np.array(list(edge_graph.values())))
    colors = np.linspace(0, np.max(color_p))
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_h = np.abs(np.array(list(har_comp.values())))
    plt.subplot(224)
    plt.title('Harmonic Component ' + str(round(weight_h*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=1,
                           node_color=np.array([[0, 88/255, 173/255, 0.5]]))
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_h,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_h\right|$')

    plt.tight_layout()
    plt.show()


# %%
del_adj = nx.adjacency_matrix(walk_del, weight='edge_visits')
del_adj = del_adj.todense()
plt.figure()
plt.matshow(del_adj)
# %%IMPORTANCE OF EACH COMPONENT

w = np.array(list(edge_del.values()))
wg = np.array(list(grad_del.values()))
ws = np.array(list(sol_del.values()))
wh = np.array(list(har_del.values()))
weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)

# %%
plt.figure()
plt.hist(list(div_del.values()), bins=50, density=True)
mean = sum(list(div_del.values()))/len(list(div_del.values()))
std = np.std(np.array(list(div_del.values())))
x = np.linspace(-60, 20, 100)


def gauss(x, mu, std):
    return(np.exp(-0.5*((x-mu)/std)**2)/(std*np.sqrt(2*np.pi)))


plt.plot(x, gauss(x, mean, std))
# %% STAR GRAPH TO SEE THE POTENTIALS

star = nx.star_graph(30)
star_dg = nx.DiGraph(star)

out_edges = [edge for edge in star_dg.edges if edge[1]
    > edge[0]]  # removing problematic edges
star_dg.remove_edges_from(out_edges)
pos = nx.kamada_kawai_layout(star_dg)

Dt = 100
v = 1
n_walk = 20
walk_star = node_walkers(star_dg, Dt, v, pos, n_walk)


# %%
grad_star, sol_star, har_star, pot_star = hodge_decomposition(
    walk_star, 'edge_visits')
# %%
edge_star = nx.get_edge_attributes(walk_star, 'edge_visits')
grad_star = {comp: round(grad_star[comp], 1) for comp in grad_star.keys()}
har_star = {comp: round(har_star[comp], 2) for comp in har_star.keys()}

plt.subplots(2, 2, figsize=(15, 15))
plt.subplot(221)
plt.title('Original Graph')

color_p = np.abs(np.array(list(edge_star.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)
print(vmin, vmax)

nx.draw_networkx_nodes(walk_star, pos=pos, node_size=10)
nx.draw_networkx_edges(walk_star, pos=pos, label=None, edge_color=color_p,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
nx.draw_networkx_edge_labels(walk_star, pos=pos, edge_labels=edge_star,
                             rotate=False, font_size=8)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega\right|$')

plt.subplot(222)
color_g = np.abs(np.array(list(grad_star.values())))
# plotting edges with color gradient

color_pot = list(pot_star.values())
colors_pot = range(int(min(color_pot)), int(max(color_pot)))
cmap_pot = plt.cm.seismic
vmin_pot = min(colors_pot)
vmax_pot = max(colors_pot)

color_p = np.abs(np.array(list(edge_star.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)
print(vmin, vmax)
plt.title('Gradient component')
nx.draw_networkx_nodes(walk_star, pos=pos, label=None,
                       node_color=color_pot, cmap=cmap_pot, vmin=vmin_pot,
                       vmax=vmax_pot, node_size=10)
nx.draw_networkx_edges(walk_star, pos=pos, label=None, edge_color=color_g,
                       edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax)
nx.draw_networkx_edge_labels(walk_star, pos=pos, edge_labels=grad_star,
                             rotate=False, font_size=8)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_g\right|$')

sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                              vmax=vmax_pot))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node potentials')

color_p = np.abs(np.array(list(edge_star.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_s = np.abs(np.array(list(sol_star.values())))
plt.subplot(223)
plt.title('Solenoidal Component')
nx.draw_networkx(walk_star, pos=pos, with_labels=False, node_size=10,
                 edge_color=color_s, edge_cmap=cmap, vmin=vmin, vmax=vmax)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_s\right|$')


color_p = np.abs(np.array(list(edge_star.values())))
colors = np.linspace(0, np.max(color_p))
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

color_h = np.abs(np.array(list(har_star.values())))
plt.subplot(224)
plt.title('Harmonic Component')
nx.draw_networkx(walk_star, pos=pos, with_labels=False, node_size=10,
                 edge_color=color_h, edge_cmap=cmap, vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
cbar = plt.colorbar(sm)
cbar.set_label(r'$\left|\omega_h\right|$')

plt.tight_layout()
plt.show()

# %% MAPS, GEOPANDAS

bcn = gpd.read_file(
    '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp', crs="EPSG:4326")
# f, ax = plt.subplots(figsize=(10,10))
# bcn.plot(ax = ax, color = "lightgray")

# Construct the primal graph
G_primal = momepy.gdf_to_nx(bcn, approach="primal",
                            multigraph=False, directed=False)

# Plot
# f, ax = plt.subplots(1, 2, figsize=(12, 6), sharex=True, sharey=True)
# bcn.plot(color="k", ax=ax[0])
# for i, facet in enumerate(ax):
#     facet.set_title(("Streets", "Graph")[i])
#     facet.axis("off")
#     add_basemap(facet)
# nx.draw(
#     G_primal, {n: [n[0], n[1]] for n in list(G_primal.nodes)}, ax=ax[1], node_size=5
# )

# %% Districtes de Barcelona

bcn_distr = gpd.read_file(
    '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp', crs="EPSG:4326")
# f, ax = plt.subplots(figsize=(10,10))
# bcn_distr.plot(ax = ax, color = "black")
# %%

nodes_df = gpd.GeoDataFrame({'PointsX': list(np.array(G_primal.nodes)[:, 0]),
                             'PointsY': list(np.array(G_primal.nodes)[:, 1])})
nodes_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(nodes_df.PointsX,
                                                         nodes_df.PointsY))
# %%
nodes_w_distr = gpd.sjoin(nodes_gdf, bcn_distr, how="left", predicate='within')
# %%
f, ax = plt.subplots(figsize=(12, 6))
nodes_w_distr.plot(color="k", ax=ax)
# %% EIXAMPLE (02)
eixample_df = nodes_w_distr.loc[nodes_w_distr['CODI_UA'] == f"{2:02}"]
# Construct the primal graph
eixample = momepy.gdf_to_nx(eixample_df, approach="primal", multigraph=False,
                            directed=True)
eixample.remove_edges_from(list(eixample.edges))
edge_list = []
for edge in G_primal.edges:
    if (edge[0] and edge[1]) in eixample.nodes:
        eixample.add_edge(edge[0], edge[1])


# %%
f, ax = plt.subplots(figsize=(12, 6))
nx.draw(eixample, pos={n: [n[0], n[1]]
        for n in list(eixample.nodes)}, node_size=2, ax=ax)
# %%
pos_to_index = {node: i for (i, node) in enumerate(eixample.nodes)}
ind_to_pos = {i: node for (i, node) in enumerate(eixample.nodes)}
eixample = nx.relabel_nodes(eixample, pos_to_index, copy=False)
# %%
f, ax = plt.subplots(figsize=(12, 6))
nx.draw(eixample, pos=ind_to_pos, node_size=2, ax=ax)
nx.set_node_attributes(eixample, ind_to_pos, name='pos')

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in eixample.edges])
mean_dist = np.mean(dists)
# %%
start_time = time.time()
v = mean_dist/3
Dt = 15
n_walk = 20
walk_eix = node_walkers(eixample, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_eix, sol_eix, har_eix, pot_eix, div_eix = hodge_decomposition(
    walk_eix, 'edge_visits')
# %%
with open('eixample_dec.txt', 'w') as eix_file:
    eix_file.write(str(grad_eix)+'\n')
    eix_file.write(str(sol_eix)+'\n')
    eix_file.write(str(har_eix)+'\n')
    eix_file.write(str(pot_eix)+'\n')
    eix_file.write(str(div_eix)+'\n')
eix_file.close()
# %%
plot_hodge(walk_eix, grad_eix, sol_eix, har_eix, pot_eix, div_eix, ind_to_pos)

# %%
plt.figure()
deg_pot = []
for node in eixample.degree:
    deg_pot.append([int(node[1]), pot_eix[node[0]]])

plt.scatter(*zip(*deg_pot))
plt.xlabel('edge degree')
plt.ylabel('edge potential')
# %%
dict_ls = []
with open('eixample_dec.txt', 'r', newline='\n') as eix_file:
    data = eix_file.readlines()
    for i in data:
        dict_ls.append(eval(i))
# %%CIUTAT VELLA (01)

cv_df = nodes_w_distr.loc[nodes_w_distr['DISTRICTE'] == f"{1:02}"]
# Construct the primal graph
cv = momepy.gdf_to_nx(cv_df, approach="primal", multigraph=False,
                            directed=True)
cv.remove_edges_from(list(cv.edges))
edge_list = []
for edge in G_primal.edges:
    if (edge[0] and edge[1]) in cv.nodes:
        cv.add_edge(edge[0], edge[1])


# %%
f, ax = plt.subplots(figsize=(12, 6))
nx.draw(cv, pos={n: [n[0], n[1]] for n in list(cv.nodes)}, node_size=2, ax=ax)
# %%
pos_to_index = {node: i for (i, node) in enumerate(cv.nodes)}
ind_to_pos = {i: node for (i, node) in enumerate(cv.nodes)}
cv = nx.relabel_nodes(cv, pos_to_index, copy=False)
self_loops = []
for edge in cv.edges:
    if edge[0] == edge[1]:
        self_loops.append(edge)
cv.remove_edges_from(self_loops)
# %%
f, ax = plt.subplots(figsize=(12, 6))
nx.draw(cv, pos=ind_to_pos, node_size=2, ax=ax)
nx.set_node_attributes(cv, ind_to_pos, name='pos')

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in cv.edges])
mean_dist = np.mean(dists)
# %%
start_time = time.time()
v = mean_dist/3
Dt = 15
n_walk = 20
walk_cv = node_walkers(cv, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grav_cv, sol_cv, har_cv, pot_cv, div_cv = hodge_decomposition(
    walk_eix, 'edge_visits')
# %%
with open('cv_dec.txt', 'w') as file:
    file.write(str(grav_cv)+'\n')
    file.write(str(sol_cv)+'\n')
    file.write(str(har_cv)+'\n')
    file.write(str(pot_cv)+'\n')
    file.write(str(div_cv)+'\n')
file.close()
# %%
plot_hodge(walk_cv, grav_cv, sol_cv, har_cv, pot_cv, div_cv, ind_to_pos)

# %%
dict_ls = []
with open('cv_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

# %%


def distr_to_nx(distr_ind:int, path_bcn: str, path_distr: str):

    bcn = gpd.read_file(path_bcn, crs="EPSG:4326")
    # f, ax = plt.subplots(figsize=(10,10))
    # bcn.plot(ax = ax, color = "lightgray")

    # Construct the primal graph
    G_primal = momepy.gdf_to_nx(
        bcn, approach="primal", multigraph=False, directed=False)
    # Districtes de Barcelona

    bcn_distr = gpd.read_file(path_distr, crs="EPSG:4326")

    nodes_df = gpd.GeoDataFrame({'PointsX': list(np.array(G_primal.nodes)[:, 0]),
                                 'PointsY': list(np.array(G_primal.nodes)[:, 1])})
    nodes_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(nodes_df.PointsX,
                                                         nodes_df.PointsY))

    nodes_w_distr = gpd.sjoin(nodes_gdf, bcn_distr,
                            how="left", predicate='within')

    
    distr_df = nodes_w_distr.loc[nodes_w_distr['DISTRICTE'] == 
                                 "{:02d}".format(distr_ind)]
    # Construct the primal graph
    distr = momepy.gdf_to_nx(distr_df, approach="primal", multigraph = False,
                                directed = True)
    distr.remove_edges_from(list(distr.edges))
    edge_list = []
    for edge in G_primal.edges:
        if (edge[0] and edge[1]) in distr.nodes:
            distr.add_edge(edge[0], edge[1])
    #changing index labelling
    pos_to_index = {node: i for (i, node) in enumerate(distr.nodes)}
    ind_to_pos = {i: node for (i, node) in enumerate(distr.nodes)}
    distr = nx.relabel_nodes(distr, pos_to_index, copy=False)
    #checking for self loops
    self_loops = []
    print('checking for self loops...')
    for edge in distr.edges:
        if edge[0] == edge[1]:
            print('self loop found', edge)
            self_loops.append(edge)
    distr.remove_edges_from(self_loops)
    print('plotting final graoh...')
    f, ax = plt.subplots(figsize=(12, 6))
    nx.draw(distr, pos = ind_to_pos,node_size = 2, ax = ax)
    nx.set_node_attributes(distr, ind_to_pos, name = 'pos')
    return(distr)