#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 19:43:37 2023

@author: robertbenassai
In this code, random walkers will be introduced to directed and undirected
graphs. The flows are going to be added to the edges and the Hodge
decomposition.
"""

import time
import json
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point, LineString
from contextily import add_basemap
import momepy
import pandas as pd
import geopandas as gpd
import scipy
from sympy import LeviCivita
from itertools import combinations, permutations
from scipy.spatial import Delaunay
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np

# %% WALKERS ON A DIRECTED GRAPH (walk during amount of steps)

'''
moves "n_walk" walkers/node in a directed graph "G" "steps" times. The initial positions of the
walkers are randomly selected amongst the nodes. The function returns the graph
with the amount of times any walker has passed through each node as attributes
to the nodes

'''


def digraph_walkers(G, steps, n_walk):

    # let's add a list of n_walk random walkers
    # overall_node_weights = {node: 0 for node in list(G.nodes)}
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # path = []
        initial_nodes = list(G.nodes)
        # initial_nodes += initial_nodes
        # path stores the position of each walker in each time step
        # path.append(initial_nodes)
        # weights counts the amount of times a walker has visited each node
        # weights = {i_node: path[0].count(i_node) for i_node in list(G.nodes)}
        edge_weights = {i_edge: 0 for i_edge in list(G.edges)}
        # now, let's move the walker to a connected node
        # first let's see the available nodes to go to
        for step in range(steps):
            # print(step)
            # list of neighboring nodes (successors) of each walker
            neighbors = [[n for n in (list(G.neighbors(walker)) +
                                      list(G.predecessors(walker)))]
                         for walker in initial_nodes]
            # now move the random walker to another node
            final_nodes = [random.choice(goto_nodes)
                                         for goto_nodes in neighbors]

            # path.append(final_nodes)

            # counting edge visits according to direction of edge
            for node_i, node_f in zip(initial_nodes, final_nodes):
                if node_i != node_f:
                    if (node_i, node_f) in list(G.edges):
                        edge_weights[(node_i, node_f)] += 1
                    else:
                        edge_weights[(node_f, node_i)] -= 1
            # count the occutpation of each node after the moves
            # for node_i, node_f in zip(initial_nodes, final_nodes):
            #     if node_i != node_f:
            #         weights[node_i] += 1

            initial_nodes = final_nodes
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
        # for node in G.nodes:
        #     overall_node_weights[node] += weights[node]
        # set node value as the number of visits of each value
        # print(weights)
    # nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)
# %% ABSORBING RANDOM WALK ON A DIRECTED GRAPH (walk during amount of steps)

'''
moves walkers in a directed graph "G" "steps" times. The initial positions of the
walkers are randomly selected amongst the nodes. The function returns the graph
with the amount of times any walker has passed through each node as attributes
to the nodes. The absorbing node is every node of the graoh once such that all 
the source-target pairs are considered.

'''


def absorbing_walkers(G, n_walk):

    # let's add a list of n_walk random walkers
    # overall_node_weights = {node: 0 for node in list(G.nodes)}
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    total_nodes = len(list(G.nodes))
#    target is the absorbing node
    for target in sorted(list(G.nodes)):
        print('target '+str(target)+'/'+str(total_nodes-1))
        for w in range(n_walk):
    #        taking the initial nodes such that target>source to avoud repetitions
            initial_nodes = [node for node in G.nodes if node<target]
            # initial_nodes += initial_nodes
            # path stores the position of each walker in each time step
            #path.append(initial_nodes)
            # weights counts the amount of times a walker has visited each node
            #weights = {i_node: path[0].count(i_node) for i_node in list(G.nodes)}
            edge_weights = {i_edge: 0 for i_edge in G.edges}
            # now, let's move the walker to a connected node
            # first let's see the available nodes to go to
            while len(initial_nodes)!=0:
                # list of neighboring nodes (successors) of each walker
                neighbors = [[n for n in (list(G.neighbors(walker)) +
                                          list(G.predecessors(walker)))]
                             for walker in initial_nodes]
                # now move the random walker to another node
                final_nodes = [random.choice(goto_nodes)
                                             for goto_nodes in neighbors]
    
                #path.append(final_nodes)
    
                # counting edge visits according to direction of edge
                for node_i, node_f in zip(initial_nodes, final_nodes):
                    if node_i != node_f:
                        if (node_i, node_f) in list(G.edges):
                            edge_weights[(node_i, node_f)] += 1
                        else:
                            edge_weights[(node_f, node_i)] -= 1
                # count the occutpation of each node after the moves
                # for node_i, node_f in zip(initial_nodes, final_nodes):
                #     if node_i != node_f:
                #         weights[node_i] += 1
    #           checking if the walkers reached the target
                if target in final_nodes:
                    final_nodes.remove(target)
                initial_nodes = final_nodes
            for edge in G.edges:
                overall_edge_weights[edge] += edge_weights[edge]/n_walk
        # for node in G.nodes:
        #     overall_node_weights[node] += weights[node]
        # set node value as the number of visits of each value
        # print(weights)
    # nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)

# %% OPTIMIZED NODE_WALKERS (time/n_wak = 11s for clusters of 100 points)
from collections import Counter
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
    Dt = Dt
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
            if len(active_ls) <= 0:  # all walkers have finished
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
                    possible_nodes = [k for k in goto_nodes.keys()] #if
                                      #(goto_nodes[k]+time[walker]) <= Dt]
                    # random selection between the final nodes
                    sel = random.choice(possible_nodes)
                    #if the choice surpasses the allowed time we deactivate the walker
                    if time[walker] + goto_nodes[sel]>Dt:
                        # desactivate the walker that has finished
                        active_walkers[walker] = False
                        print(Counter(active_walkers.values()))
                        continue
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
                    print(Counter(active_walkers.values()))
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G)  # MSD_dict)

#%%
''' This funciton moves n_walk random walkers from each node in the graph 
during a time Dt in a periodic boundary condition square lattice.

G: an nx.DiGraph
Dt: total simulation time
v: velocity of the walker
pos: dictionary with node_id:position of the node
N: number of nodes of the original square lattice
follow_node: int, id of the node to follow the trajectory
'''
def periodic_walkers(G, Dt, v, pos, n_walk, N, follow_node):
    # calculating the time intrinsec to each edge
    delta_t = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
            for edge in G.edges}
    path_tot = []
    lattice_gap = 1
    for i in range(N):
        delta_t[(i, i+N*N-N)] = lattice_gap
    #    print((i, i+N*N-N) in list(PBC_lattice.edges))
        if i==0:
            delta_t[(i,i+N-1)] = lattice_gap 
    #        print((i,i+N-1) in list(PBC_lattice.edges))
        else:
            delta_t[(i*N, i*N+N-1)] = lattice_gap
    #        print((i*N, i*N+N-1) in list(PBC_lattice.edges))
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
        path = []
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {i: i for i in G.nodes}
        path.append(pos_walkers[follow_node])
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
            if len(active_ls) <= 0:  # all walkers have finished
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
                    possible_nodes = [k for k in goto_nodes.keys()] #if
                                      #(goto_nodes[k]+time[walker]) <= Dt]
                    # random selection between the final nodes
                    sel = random.choice(possible_nodes)
                    #if the choice surpasses the allowed time we deactivate the walker
                    if time[walker] + goto_nodes[sel]>Dt:
                        # desactivate the walker that has finished
                        active_walkers[walker] = False
                        # print(Counter(active_walkers.values()))
                        continue
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
                    # print(Counter(active_walkers.values()))
            path.append(pos_walkers[follow_node])
        path_tot.append(path)
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G, path_tot)  # MSD_dict)
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
    pot_nodes = {n: pot for n, pot in zip(G.nodes, pot_field)}
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
        print('max error', np.max(big_arr.dot(g_har_vec)))
        # print('divergence of the harmonic component',
        #       grad_arr.dot(np.transpose(grad_arr)).dot(g_har_vec))
        # print('curl of the harmonic component',
        #       curl_op.dot(g_har_vec))

        har_comp = {edge: har for edge, har in zip(G.edges, g_har_vec)}
        # for edge, har in zip(G.edges, g_har_vec):
        #     har_comp[edge] = har
    return grad_comp, sol_comp, har_comp, pot_nodes, div
# %%

def plot_hodge(walk_graph, grad_comp, sol_comp, har_comp, pot, div, pos):

    edge_graph = nx.get_edge_attributes(walk_graph, 'edge_visits')
    
    g_field = np.array([walk_graph[edge[0]][edge[1]]['edge_visits'] for edge
                                     in walk_graph.edges])
    percentile = np.percentile(g_field, 95)
    

    w = np.array(list(edge_graph.values()))
    wg = np.array(list(grad_comp.values()))
    ws = np.array(list(sol_comp.values()))
    wh = np.array(list(har_comp.values()))
    weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
    weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
    weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

    print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)
    # und_walk = walk_graph.to_undirected()
    # betweenness_cent = nx.betweenness_centrality(und_walk, weight = 
    #                                              'dt', normalized =False)
    deg = walk_graph.degree()
    deg_dict = {node[0]: node[1] for node in deg}
    plt.subplots(2, 2, figsize=(15, 15))
    plt.subplot(221)
    plt.title('Original Graph 100%')

    # color_p = np.abs(np.array(list(edge_graph.values())))
    color_p = np.array(list(edge_graph.values()))
    # colors = np.linspace(0, percentile)
    # cmap = plt.cm.Oranges
    colors = np.linspace(-percentile, percentile)
    cmap = plt.cm.seismic
    vmin = min(colors)
    vmax = max(colors)

    # color_div = list(div.values())
    # colors_div = range(int(min(color_div)), int(max(color_div)))
    # cmap_div = plt.cm.RdBu_r
    # vmin_div = min(colors_div)
    # vmax_div = round(max(color_div))

    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4, 
                           node_color='#D3D3D3')
                            # node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                            # vmax=vmax_div)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_p,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 4)

    sm = plt.cm.ScalarMappable(cmap=cmap, 
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega\right|$')

    # sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
    #                                                               vmax=vmax_div))
    # sm2._A = []
    # cbar2 = plt.colorbar(sm2, location='right')
    # cbar2.set_label(r'Node div')

    plt.subplot(222)
    color_g = np.abs(np.array(list(grad_comp.values())))
    # plotting edges with color gradient

    color_pot = list(pot.values())
    colors_pot = range(-int(max(color_pot)), int(max(color_pot)))
    cmap_pot = plt.cm.PRGn
    vmin_pot = np.min(color_pot)
    vmax_pot = np.max(color_pot)

    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)
    plt.title('Gradient component ' + str(round(weight_g*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None,
                            node_size=8, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 5, node_size = 8)

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

    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_s = np.abs(np.array(list(sol_comp.values())))
    plt.subplot(223)
    plt.title('Solenoidal Component ' + str(round(weight_s*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
                           node_color='#D3D3D3')
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_s,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 4)


    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_s\right|$')

    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_h = np.array(list(har_comp.values()))
    plt.subplot(224)
    plt.title('Harmonic Component ' + str(round(weight_h*100, 1))+'%')
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
                           node_color='#D3D3D3')
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_h,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 4)
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_h\right|$')

    plt.tight_layout()
    # plt.savefig('/Users/robertbenassai/Documents/UOC/figs/lattice_evolution/lattice_dt'+str(Dt)+'.png')
    plt.show()


# %% MAPS, GEOPANDAS
import time
import momepy
def distr_to_nx(distr_ind:int, path_bcn: str, path_distr: str):
    
    bcn_edges = gpd.read_file(path_bcn, crs="EPSG:25831")

    bcn_nodes = gpd.read_file('/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/vertices.shp',
                              crs="EPSG:25831")

    missing_edges = pd.read_csv(
        '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/missing_edges.csv')
    bcn_distr = gpd.read_file('/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp',
                              crs="EPSG:25831")
    nodes_w_distr = gpd.sjoin(bcn_nodes, bcn_distr,
                            how="left", predicate='within')
    distr_df = nodes_w_distr.loc[nodes_w_distr['DISTRICTE'] == 
                                 "{:02d}".format(distr_ind)]
    distr = nx.DiGraph()
    distr.add_nodes_from(distr_df['uid'])
    ind_to_pos = {nodes_w_distr.loc[i,'uid']:(nodes_w_distr.loc[i,'geometry'].x,
                                          nodes_w_distr.loc[i,'geometry'].y) for i in 
                  range(len(nodes_w_distr))}
    # pos_to_index = {(nodes_w_distr.loc[i,'geometry'].x, 
    #                  nodes_w_distr.loc[i,'geometry'].y):nodes_w_distr.loc[i,'uid']
    #                 for i in range(len(nodes_w_distr))}
    #nx.set_node_attributes(bcn_graph, ind_to_pos, 'pos')

    for edge in bcn_edges['to_from']:
        sep = edge.split(',')
        if len(sep) == 2:
            if int(sep[0]) and int(sep[1]) in distr.nodes:
                distr.add_edge(int(sep[0]), int(sep[1]))
                
    for i,j in zip(missing_edges['FROM_NODE'], missing_edges['TO_NODE']):
        if (i and j) in distr.nodes:
            distr.add_edge(i, j)
    
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
    plt.figure()
    for component in remove_nodes:
        subg = distr.subgraph(component)
        nx.draw_networkx(subg, pos = ind_to_pos, with_labels=False,
                         node_color = [[random.uniform(0, 1),
                                        random.uniform(0, 1),
                                        random.uniform(0, 1)]],
                         node_size = 3)
    remove_nodes.remove(max(remove_nodes, key = len))
    for node_set in remove_nodes:
        distr.remove_nodes_from(node_set)
    print('plotting final graph...')
    f, ax = plt.subplots(figsize=(12, 6))
    nx.draw(distr, pos = ind_to_pos,node_size = 2, ax = ax)
    nx.set_node_attributes(distr, ind_to_pos, name = 'pos')
    return(distr, ind_to_pos)

#%%Structural_ratios
#calculates the ratio between the dimension of each hodge component and the dim
#of the original graph
def structural_ratios(G):
    dim_orig = G.number_of_edges()
    n = G.number_of_nodes()
    con_comp = len(list(nx.connected_components(G.to_undirected())))
    dim_grad = n-con_comp
    
    dim_curl, _ = find_triangles(G) #dimension of curl is the number of linearly 
                                    #independent triangles
    return(dim_grad/dim_orig, dim_curl/dim_orig, 1-dim_grad/dim_orig-dim_curl/dim_orig)
    
    
# %% EIXAMPLE (02)
'''L'EIXAMPLE'''
distr_ind = 2
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

eixample, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)
# %%
f, ax = plt.subplots(figsize=(12, 6))
nx.draw(eixample, pos=ind_to_pos, node_size=2, ax=ax)
nx.set_node_attributes(eixample, ind_to_pos, name='pos')

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in eixample.edges])
mean_dist = np.mean(dists)
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(eixample)

print(st_ratio_g, st_ratio_s, st_ratio_h)
# %%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_eix = node_walkers(eixample, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_eix, sol_eix, har_eix, pot_eix, div_eix = hodge_decomposition(
    walk_eix, 'edge_visits')
# %%
with open('eixample_dec_dt_15_connected.txt', 'w') as eix_file:
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
'''CIUTAT VELLA'''

distr_ind = 1
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

cv, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
import scipy.stats as ss
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in cv.edges])
#MLE

P = ss.expon.fit(dists)

#plotting
rX = np.linspace(0,100, 100)
rP = ss.expon.pdf(rX, *P)
#Yup, just unpack P with *P, instead of scale=XX and shape=XX, etc.

#need to plot the normalized histogram with `normed=True`
plt.figure()
plt.plot(rX, rP, color = 'r', linestyle = '--')
plt.xlim(0,100)
plt.hist(dists, bins = 100,range=(0,100), density=True, color = 'orange')
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(cv)

print(st_ratio_g, st_ratio_s, st_ratio_h)
#%%
#mean_dist = np.mean(dists)
start_time = time.time()
v = 1.42
Dt = 13*60
n_walk = 20
walk_cv = node_walkers(cv, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
#%%
grad_cv, sol_cv, har_cv, pot_cv, div_cv = hodge_decomposition(
    walk_cv, 'edge_visits')
# %%
with open('cv_dec_dt_13_connected.txt', 'w') as file:
    file.write(str(grad_cv)+'\n')
    file.write(str(sol_cv)+'\n')
    file.write(str(har_cv)+'\n')
    file.write(str(pot_cv)+'\n')
    file.write(str(div_cv)+'\n')
file.close()

# %%
from scipy.stats import norm
import matplotlib.mlab as mlab
g_field = np.array([walk_cv[edge[0]][edge[1]]['edge_visits'] for edge
                                 in walk_cv.edges])

# best fit of data
(mu, sigma) = norm.fit(g_field)

# the histogram of the data
n, bins, _ = plt.hist(g_field, 60, facecolor='violet', alpha=0.75,density=True)

# add a 'best fit' line
y = norm.pdf(bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=2)

#plot
plt.xlabel('Edge flow')
plt.ylabel('Frequency')
plt.title(r'$\mu=%.3f,\ \sigma=%.3f$' %(mu, sigma))

#%%
pot_field = np.array(list(pot_cv.values()))
#plot
plt.xlabel('Node potential')
plt.ylabel('Frequency')
plt.hist(pot_field, bins = 50, edgecolor ='black')

#%%
plot_hodge(walk_cv, grad_cv, sol_cv, har_cv, pot_cv, div_cv, ind_to_pos)

#%%
dict_ls = []
with open('cv_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()
#%%
grad_w = [72.6,68.4, 63.6, 61.8, 59.1, 56.6, 52.5, 50.5, 51.1]
har_w = [26.2, 30.6, 34.7, 36.8, 39.1, 41.5, 45.7, 47.5, 46.9]
time_ls = [5, 7, 10, 13, 15, 20, 25, 30, 35]

plt.figure()
plt.plot(time_ls, grad_w, color = 'b' , linestyle = '-', marker = '.', label = 'gradient strength ratio')
plt.plot(time_ls, har_w, color = 'r', linestyle = '-', marker = '.', label = 'harmonic strength ratio')
plt.hlines(st_ratio_g*100, 0, 35, color = 'b' , linestyle = '--', label = 'gradient structural ratio')
plt.hlines(st_ratio_h*100, 0, 35, color = 'r', linestyle = '--', label = 'harmonic structural ratio')
plt.xlabel('Simulation time')
plt.ylabel('strength ratio')
plt.legend()

#%% discrete walk
start_time = time.time()
steps = 80
n_walk = 20
walk_cv = digraph_walkers(cv, steps, n_walk) #10, 20, 30, 40, 60, 70, 80
print("--- %s seconds ---" % (time.time() - start_time))
#%%
grad_cv, sol_cv, har_cv, pot_cv, div_cv = hodge_decomposition(
    walk_cv, 'edge_visits')
#%%
plot_hodge(walk_cv, grad_cv, sol_cv, har_cv, pot_cv, div_cv, ind_to_pos)
#%%
grad_w = [51.4, 38.7, 31.6, 27.6, 22.6, 22.0, 19.6, 20.8, 18.6]
har_w = [46.3, 59.1, 65.2, 69.1, 74.4, 73.8, 77.3, 76.4, 78.4]
time_ls = [10,20,30,40,50,60,70,80,100]

plt.figure()
plt.plot(time_ls, grad_w, color = 'b' , linestyle = '-', marker = '.', label = 'gradient strength ratio')
plt.plot(time_ls, har_w, color = 'r', linestyle = '-', marker = '.', label = 'harmonic strength ratio')
plt.hlines(st_ratio_g*100, 0, 100, color = 'b' , linestyle = '--', label = 'gradient structural ratio')
plt.hlines(st_ratio_h*100, 0, 100, color = 'r', linestyle = '--', label = 'harmonic structural ratio')
plt.xlabel('Simulation time')
plt.ylabel('strength ratio')
plt.legend()
#%% SANTS MONTJUÏC
'''SANTS - MONTJUÏC'''
distr_ind = 3
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

sts_mj, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in sts_mj.edges])

plt.xlim(0,100)
P =ss.pareto.fit(dists)
rX = np.linspace(0,100, 100)
rP = ss.pareto.pdf(rX, *P)
Pe =ss.expon.fit(dists)
rPe = ss.expon.pdf(rX, *Pe)
plt.plot(rX,rP, 'r--', label = 'pareto fit')
plt.plot(rX, rPe, 'b--', label = 'exponential fit')
plt.hist(dists, bins = 100,range=(0,100), density = True)
plt.legend()
print(min(dists))
#%%
mean_dist = np.mean(dists)
print(mean_dist)
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_sts_mj = node_walkers(sts_mj, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(sts_mj)

print(st_ratio_g, st_ratio_s, st_ratio_h)
# %%
grad_sts_mj, sol_sts_mj, har_sts_mj, pot_sts_mj, div_sts_mj = hodge_decomposition(
    walk_sts_mj, 'edge_visits')
# %%
with open('sts_mj_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_sts_mj)+'\n')
    file.write(str(sol_sts_mj)+'\n')
    file.write(str(har_sts_mj)+'\n')
    file.write(str(pot_sts_mj)+'\n')
    file.write(str(div_sts_mj)+'\n')
file.close()
# %%
plot_hodge(walk_sts_mj, grad_sts_mj, sol_sts_mj, har_sts_mj, pot_sts_mj, div_sts_mj, ind_to_pos)

# %%
dict_ls = []
with open('sts_mj_dec_dt_15_connected.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()
#%%

und_walk = walk_sts_mj.to_undirected()
betweenness_cent = nx.betweenness_centrality(und_walk, weight = 
                                             'dt', normalized =False)
pot_vs_cent = [[pot_sts_mj[node],betweenness_cent[node]] for node in pot_sts_mj.keys()]


#%% LES CORTS
'''LES CORTS'''
distr_ind = 4
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

corts, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in corts.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(corts)

print(st_ratio_g, st_ratio_s, st_ratio_h)
#%%
mean_dist = np.mean(dists)
print(mean_dist)
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_corts = node_walkers(corts, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_corts, sol_corts, har_corts, pot_corts, div_corts = hodge_decomposition(
    walk_corts, 'edge_visits')
# %%
with open('corts_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_corts)+'\n')
    file.write(str(sol_corts)+'\n')
    file.write(str(har_corts)+'\n')
    file.write(str(pot_corts)+'\n')
    file.write(str(div_corts)+'\n')
file.close()
# %%
plot_hodge(walk_corts, grad_corts, sol_corts, har_corts, pot_corts, div_corts, ind_to_pos)

# %%
dict_ls = []
with open('corts_dec_dt_15_connected.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()
#%%
walks = nx.get_edge_attributes(walk_cv, 'edge_visits')
dists_dict = {edge: np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in cv.edges}

ls = [[dists_dict[edge], abs(walks[edge])] for edge in walks.keys()]
plt.figure()

plt.scatter(*zip(*ls), s = 10)
plt.xlabel('edge distance')
plt.ylabel(r'$\left|\omega\right|$')
#%% SARRIÀ - ST. GERVASI
'''SARRIÀ - SANT GERVASI'''
distr_ind = 5
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

sarr, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in sarr.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(sarr)

print(st_ratio_g, st_ratio_s, st_ratio_h)
#%%
mean_dist = np.mean(dists)
print(mean_dist)
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_sarr = node_walkers(sarr, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_sarr, sol_sarr, har_sarr, pot_sarr, div_sarr = hodge_decomposition(
    walk_sarr, 'edge_visits')
# %%
with open('sarr_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_sarr)+'\n')
    file.write(str(sol_sarr)+'\n')
    file.write(str(har_sarr)+'\n')
    file.write(str(pot_sarr)+'\n')
    file.write(str(div_sarr)+'\n')
file.close()
# %%
plot_hodge(walk_sarr, grad_sarr, sol_sarr, har_sarr, pot_sarr, div_sarr, ind_to_pos)

# %%
dict_ls = []
with open('sarr_dec_dt_15_connected.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

#%% GRÀCIA
'''GRÀCIA'''
distr_ind = 6
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

gracia, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in gracia.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_gracia = node_walkers(gracia, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_gracia, sol_gracia, har_gracia, pot_gracia, div_gracia = hodge_decomposition(
    walk_gracia, 'edge_visits')
# %%
with open('gracia_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_gracia)+'\n')
    file.write(str(sol_gracia)+'\n')
    file.write(str(har_gracia)+'\n')
    file.write(str(pot_gracia)+'\n')
    file.write(str(div_gracia)+'\n')
file.close()
# %%
plot_hodge(walk_gracia, grad_gracia, sol_gracia, har_gracia, pot_gracia, div_gracia, ind_to_pos)

# %%
dict_ls = []
with open('gracia_dec_dt_15_connected.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

#%% HORTA-GUINARDÓ
'''HORTA-GUINARDÓ'''
distr_ind = 7
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

horta, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in horta.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_horta = node_walkers(horta, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_horta, sol_horta, har_horta, pot_horta, div_horta = hodge_decomposition(
    walk_horta, 'edge_visits')
# %%
with open('horta_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_horta)+'\n')
    file.write(str(sol_horta)+'\n')
    file.write(str(har_horta)+'\n')
    file.write(str(pot_horta)+'\n')
    file.write(str(div_horta)+'\n')
file.close()
# %%
plot_hodge(walk_horta, grad_horta, sol_horta, har_horta, pot_horta, div_horta, ind_to_pos)

# %%
dict_ls = []
with open('horta_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()
#%% NOU BARRIS
'''NOU BARRIS'''
distr_ind = 8
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

noub, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in noub.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
mean_dist = np.mean(dists)
print(mean_dist)
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_noub = node_walkers(noub, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_noub, sol_noub, har_noub, pot_noub, div_noub = hodge_decomposition(
    walk_noub, 'edge_visits')
# %%
with open('noub_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_noub)+'\n')
    file.write(str(sol_noub)+'\n')
    file.write(str(har_noub)+'\n')
    file.write(str(pot_noub)+'\n')
    file.write(str(div_noub)+'\n')
file.close()
# %%
plot_hodge(walk_noub, grad_noub, sol_noub, har_noub, pot_noub, div_noub, ind_to_pos)

# %%
dict_ls = []
with open('noub_dec_dt_15_connected.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

#%% ST. ANDREU
'''SANT ANDREU'''
distr_ind = 9
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

st_and, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in st_and.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_st_and = node_walkers(st_and, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_st_and, sol_st_and, har_st_and, pot_st_and, div_st_and = hodge_decomposition(
    walk_st_and, 'edge_visits')
# %%
with open('st_and_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_st_and)+'\n')
    file.write(str(sol_st_and)+'\n')
    file.write(str(har_st_and)+'\n')
    file.write(str(pot_st_and)+'\n')
    file.write(str(div_st_and)+'\n')
file.close()
# %%
plot_hodge(walk_st_and, grad_st_and, sol_st_and, har_st_and, pot_st_and, div_st_and, ind_to_pos)

# %%
dict_ls = []
with open('st_and_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

#%% ST. MARTI

'''SANT MARTÍ'''
distr_ind = 10
path_bcn = '/Users/robertbenassai/Documents/UOC/alertadadesconfidencialsdatosconfidencialesconfid/edges.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

st_mart, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

# %%
dists = np.array([np.linalg.norm(np.array(ind_to_pos[edge[1]])-np.array(ind_to_pos[edge[0]]))
         for edge in st_mart.edges])

n, bins, _ = plt.hist(dists, bins = 100,range=(0,100))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
#%%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
walk_st_mart = node_walkers(st_mart, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
grad_st_mart, sol_st_mart, har_st_mart, pot_st_mart, div_st_mart = hodge_decomposition(
    walk_st_mart, 'edge_visits')
# %%
with open('st_mart_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_st_mart)+'\n')
    file.write(str(sol_st_mart)+'\n')
    file.write(str(har_st_mart)+'\n')
    file.write(str(pot_st_mart)+'\n')
    file.write(str(div_st_mart)+'\n')
file.close()
# %%
plot_hodge(walk_st_mart, grad_st_mart, sol_st_mart, har_st_mart, pot_st_mart, div_st_mart, ind_to_pos)

# %%
dict_ls = []
with open('st_mart_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()

#%%BUILDING A LATTICE GRAPH WITH PERIODIC BOUNDARY CONDITIONS AND ADDING NODES
#AT THE CENTRE
N = 8
PBC_lattice =nx.grid_2d_graph(N,N,periodic=True)
new_edges = []
# Define the vertices of the original square
A = (3, 3)
B = (3, 4)
C = (4, 4)
D = (4, 3)

# Define the number of iterations (i.e., how many times to divide each square into four smaller squares)
num_iterations = 2

# Initialize the list of squares to divide with the original square
squares_to_divide = [(A, B, C, D)]

# Iterate through the specified number of iterations
for iteration in range(num_iterations):
    new_squares = []
    # Divide each square in the current list into four smaller squares
    for square in squares_to_divide:
        # Find the midpoints of each side of the square
        AB_midpoint = ((square[0][0] + square[1][0]) / 2, (square[0][1] + square[1][1]) / 2)
        BC_midpoint = ((square[1][0] + square[2][0]) / 2, (square[1][1] + square[2][1]) / 2)
        CD_midpoint = ((square[2][0] + square[3][0]) / 2, (square[2][1] + square[3][1]) / 2)
        DA_midpoint = ((square[3][0] + square[0][0]) / 2, (square[3][1] + square[0][1]) / 2)

        # Find the midpoint of the square
        square_midpoint = ((square[0][0] + square[2][0]) / 2, (square[0][1] + square[2][1]) / 2)

        # Find the midpoint of the diagonal of the square
        diagonal_midpoint = ((square[0][0] + square[2][0]) / 2, (square[1][1] + square[3][1]) / 2)

        # Define the four smaller squares
        square_1 = (square[0], AB_midpoint, square_midpoint, DA_midpoint)
        square_2 = (AB_midpoint, square[1], BC_midpoint, square_midpoint)
        square_3 = (square_midpoint, BC_midpoint, square[2], CD_midpoint)
        square_4 = (DA_midpoint, square_midpoint, CD_midpoint, square[3])

        # Add the four smaller squares to the list of new squares
        new_squares.extend([square_1, square_2, square_3, square_4])

    # Update the list of squares to divide with the new list of smaller squares
    squares_to_divide = new_squares

# Find the vertexes of the final squares excluding the ones from the original square
new_vertexes = []
for square in squares_to_divide:
    for vertex in square:
        if vertex not in [A, B, C, D] and vertex not in new_vertexes:
            new_vertexes.append(vertex)

# Print the new vertexes
print(new_vertexes)


#%%ADDING THE NEW EDGES TO THE GRAPH
for i in range(0, 4):
    new_edges += [((4,3+i/4), (4, 3+(i+1)/4))]
    new_edges += [((3+i/4,4), (3+(i+1)/4, 4))]
    for j in range(0, 4):
        new_edges += [((3+i/4, 3+j/4), (3+i/4, 3+(j+1)/4)), ((3+i/4, 3+j/4), (3+(i+1)/4, 3+j/4)), 
                      ((3+i/4, 3+j/4), (3+(i+1)/4, 3+(j+1)/4)), (((3+i/4, 3+(j+1)/4)), (3+(i+1)/4, 3+j/4))]
PBC_lattice.remove_edges_from([((3,3),(3,4)), ((3, 3), (4,3)), ((3,4), (4,4)),
                               ((4,3), (4,4))])
PBC_lattice.add_nodes_from(new_vertexes)
PBC_lattice.add_edges_from(new_edges)
#pos = {i * N + j:(i, j) for i, j in PBC_lattice.nodes()}
pos = {i: j for i,j in enumerate(PBC_lattice.nodes)}
#labels = dict( ((i, j), i * N + j) for i, j in PBC_lattice.nodes() )
labels = {i: k for k, i in enumerate(PBC_lattice.nodes)}
nx.relabel_nodes(PBC_lattice, labels, copy=False)

PBC_lattice = nx.DiGraph(PBC_lattice)
out_edges = [edge for edge in PBC_lattice.edges if edge[0]
    > edge[1]]  # removing all outward edges
PBC_lattice.remove_edges_from(out_edges)
nx.draw_networkx(PBC_lattice, pos=pos, with_labels=True, node_size = 10)
#%% PERIODIC BOUNDARY CONDITIONS CONTINUOUS RANDOM WALK
v = 1
Dt = 20
# for Dt in range(95,105, 5):
n_walk = 20
follow_node = 69
walked_PBC, path = periodic_walkers(PBC_lattice, Dt, v, pos, n_walk, N, follow_node)

grad, sol, har, pot, div = hodge_decomposition(
    walked_PBC, 'edge_visits')
print(structural_ratios(PBC_lattice))
plot_hodge(walked_PBC, grad, sol, har, pot, div, pos)
#%%
for k in range(n_walk):
    plt.figure()
    path_edges = [(path[k][i],path[k][i+1])  if path[k][i]<path[k][i+1] else
                  (path[k][i+1],path[k][i]) for i in range(len(path[k])-1)]
    edge_color_list = ['red' if edge in path_edges
                       else (210/255, 210/255, 210/255, 0.2) for edge in PBC_lattice.edges]
    final_nodes = [path[i][-1] for i, j in enumerate(path)]
    node_color_list = ['orange' if node in final_nodes else 'blue' for node in PBC_lattice.nodes]
    nx.draw_networkx(PBC_lattice, pos=pos, with_labels=True, node_size = 10, 
                      edge_color = edge_color_list, node_color=node_color_list)
    plt.title('Final nodes for node '+str(follow_node))
    plt.tight_layout()
print(structural_ratios(PBC_lattice))
#%%
n_w_ls = np.array([10, 20, 30, 40, 50, 100])*len(list(PBC_lattice.nodes))
grad_pc = [14.1, 23.6, 33.7, 37.4, 39.4, 58.4]

plt.figure()
plt.plot(n_w_ls, grad_pc, '.-')
plt.xlabel('total number of walkers')
plt.ylabel('gradient component (%)')
plt.title('Periodic random walk for Dt  = 15')
plt.tight_layout()

#%% absorbing walk on a lattice

walk_abs = absorbing_walkers(PBC_lattice, 20)
grad, sol, har, pot, div = hodge_decomposition(
    walk_abs, 'edge_visits')

plot_hodge(walk_abs, grad, sol, har, pot, div, pos)


#%%NEWMANN BETWENNESS FUNCTION

def newman_bet(A:np.array,lap:np.array):
    '''
    

    Parameters
    ----------
    A : np.array
        Adjacency matrix of the graph.
    lap : np.array
        Graph laplacian.

    Returns
    -------
    rw_edge_betw: array of edge flows according to Newmann's RW betweenness.

    '''
    #1 constructing D-A
    n = len(A)
    #rw_node_betw = np.zeros(n)
    rw_edge_btw = np.zeros(np.shape(lap))
    for target in range(n):
        print(target)
        M = np.copy(lap)
        #2 deleting row and column of absorbing node
        M = np.delete(M, target, 0)
        M = np.delete(M, target, 1)
        #3 invert matrix
        M_inv = np.linalg.pinv(M)
        M_inv = np.insert(M_inv, target, np.zeros(np.shape(M_inv)[0]), axis = 1)
        M_inv = np.insert(M_inv, target, np.zeros(np.shape(M_inv)[1]), axis = 0)
        for source in range(target):
            V = M_inv[:, source]-M_inv[:, target]
            #current = np.zeros(np.shape(V))
            for i, V_i in enumerate(V):
                #current[i] = np.sum(A[i,:].dot(np.abs(V[i]-V)))
                rw_edge_btw[i] += V_i-V
            #current *= 0.5
            #rw_node_betw += current
    rw_edge_btw = np.multiply(rw_edge_btw,A)
    # const = 1/(0.5*n*(n-1))
    # rw_edge_btw *= const
    # rw_node_betw *= const
    return(rw_edge_btw)#, rw_node_betw)
#%% running newmann
start_time = time.time()
A = nx.adjacency_matrix(nx.to_undirected(PBC_lattice), nodelist= sorted(PBC_lattice.nodes)).toarray()
lap = nx.laplacian_matrix(nx.to_undirected(PBC_lattice), nodelist= sorted(PBC_lattice.nodes)).toarray()
rw_edge_btw = newman_bet(A, lap)
print("--- %s seconds ---" % (time.time() - start_time))
#%% plot newmann
n = len(A)
# const = (0.5*n*(n-1))
# rw_edge_btw *= const
# rw_node_betw *= const
rw_bet_dict = {edge: rw_edge_btw[edge[0]][edge[1]] for edge in 
               PBC_lattice.edges}

nx.set_edge_attributes(PBC_lattice, rw_bet_dict, name = 'edge_visits')
grad_PBC, sol_PBC, har_PBC, pot_PBC, div_PBC = hodge_decomposition(
    PBC_lattice, 'edge_visits')
print(pot_PBC)
plot_hodge(PBC_lattice, grad_PBC, sol_PBC, har_PBC, pot_PBC, div_PBC, pos)
#%% Let's see if the gradient component is always the same and see the fluctuations
boxpl_ls_g = []
boxpl_ls_h = []
for i in range(5):
    walk_abs = absorbing_walkers(PBC_lattice, 50)
    grad, sol, har, pot, div = hodge_decomposition(
        walk_abs, 'edge_visits')

    diff_h = [har_PBC[edge] - har[edge] for edge in PBC_lattice.edges]
    diff_g = [grad_PBC[edge] - grad[edge] for edge in PBC_lattice.edges]
    boxpl_ls_g += diff_g
    boxpl_ls_h += diff_h
boxpl = []
boxpl.append(boxpl_ls_g)
boxpl.append(boxpl_ls_h)
#%%
plt.figure()
plt.boxplot(boxpl, labels =  ['Gradient discrepancy', 'Harmonic discrepancy'])
plt.ylabel('Edge discrepancy from predicted component')
#%%
A = nx.adjacency_matrix(nx.to_undirected(cv), nodelist= sorted(cv.nodes)).toarray()
lap = nx.laplacian_matrix(nx.to_undirected(cv), nodelist= sorted(cv.nodes)).toarray()
rw_edge_btw, rw_node_betw = newman_bet(A, lap)

#%% CONTINUUM NEWMANN WITHOUT ABSORBINGS

'''first we build the transition rates matrix'''

# building the transition rate matrix
#dict of matrix ind to node index
PBC_und = PBC_lattice.to_undirected()

#distances of the edges in the graph
dists = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))
        for edge in PBC_und.edges}
path_tot = []
lattice_gap = 1
for i in range(N):
    dists[(i, i+N*N-N)] = lattice_gap
    dists[(i+N*N-N, i)] = lattice_gap
    if i==0:
        dists[(i,i+N-1)] = lattice_gap 
        dists[(i+N-1, i)] = lattice_gap 
    else:
        dists[(i*N, i*N+N-1)] = lattice_gap
        dists[(i*N+N-1, i*N)] = lattice_gap

#mapping of each node to its index in the transition rate array
inode_to_iarr = {node:i for i, node in enumerate(PBC_und.nodes)}
iarr_to_inode = {inode_to_iarr[i]:i for i in inode_to_iarr.keys()}
#building the array
trans_rates = np.zeros((len(PBC_und.nodes), len(PBC_und.nodes)))
for node in PBC_und.nodes:
    #Degree of the departure node
    k = len(list(PBC_und.neighbors(node)))
    # k = np.sum(np.array([dists[(node, neigh)] if (node, neigh) in list(PBC_und.edges)
                         # else dists[(neigh, node)]
                         # for neigh in PBC_und.neighbors(node)]))
                         
    #for each neighbour we calculate the transition rate probability of going 
    #from node to final as v/(dist*k)
    for final in PBC_und.neighbors(node):
        i, j = inode_to_iarr[node], inode_to_iarr[final]
        dist = np.linalg.norm(np.array(pos[final])-np.array(pos[node]))
        trans_rates[i][j] = v/(k*dist)
        # trans_rates[i][j] = 1/k
# the diagonal is -sum of the off diag elements of the col
i,j = np.indices(trans_rates.shape)
diagonal = np.sum(trans_rates, axis=0)
trans_rates[i==j] = -diagonal
#%% FINDING THE EIGENVECS AND EIGENVALS. The eigenvector with eigenvalue 1 is 
#the stationary state in theory
eigval, eigvecs = scipy.linalg.eig(np.transpose(trans_rates))
#%%
print(np.where(np.real(np.round(eigval, 5)) == 0))
stationary_node_prob = np.real(eigvecs)[:,np.where(np.real(np.round(eigval, 2)) == 0)[0]]
print(Counter(np.transpose(stationary_node_prob)[0]))
#%% SYSTEM OF ODES
from scipy.integrate import odeint
v=1
def dp_dt(p, t, R, dummy):
    '''

    Parameters
    ----------
    p : np.array
        Array of probabilities of being at each node at time t.
    t : np.array
        Array (linspace) of integration times.
    R : np.array
        Transition rates matrix, each coefficient of the differential equations.
    dummy : None
        dummy vvariable for the odeint to work.

    Returns
    -------
    dy: dp.

    '''
    dy = R.dot(p)
    return(dy)
#array of flows through the edges
net_flow = np.zeros(np.shape(trans_rates))
#initial probabilities
p0 = np.zeros(len(PBC_lattice.nodes))
for initial in range(len(p0)):
    p0 = np.zeros(len(PBC_lattice.nodes))
    p0[initial] = 1
    
    t = np.linspace(start=0, stop=20,num=20000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    
    for j, p_j in enumerate(flows):
        
        for i in range(len(PBC_lattice.nodes)):
            net_flow[i,j] -= p_j*trans_rates[i][j] - flows[i]*trans_rates[j][i]

net_flow*= 20
cont_new = {edge: net_flow[inode_to_iarr[edge[0]]][inode_to_iarr[edge[1]]] 
            for edge in PBC_lattice.edges}
nx.set_edge_attributes(PBC_lattice, cont_new, name = 'edge_visits')
#%%
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = \
hodge_decomposition(PBC_lattice, 'edge_visits')
#pot_cont_new = {iarr_to_inode[i]: prob for i, prob in enumerate(sol[-1][:])}
plot_hodge(PBC_lattice, g_cont_new, s_cont_new, h_cont_new, pot_cont_new, 
           div_cont_new, pos)
#%%
plt.figure()
sol = pd.DataFrame(sol)
for i in range(len(PBC_lattice.nodes)):
    plt.plot(t,sol[i])
plt.show()
print(8/(2*len(list(PBC_lattice.edges))))
#%%
for i in pot.keys():
    print(pot[i], pot_cont_new[i])
    
#%% WHAT HAPPENS IF WE PERFORM A DISCRETE RW

Dt = 200
# for Dt in range(95,105, 5):
n_walk = 20

walked_PBC_discr = digraph_walkers(PBC_lattice, Dt, n_walk)

grad, sol, har, pot, div = hodge_decomposition(
    walked_PBC_discr, 'edge_visits')

plot_hodge(walked_PBC_discr, grad, sol, har, pot, div, pos)   
#%% NEWMANN RNADOM WALK BETWEENNESS USING NX INTEGRATED FUNCTION
edge_rw_centr = nx.edge_current_flow_betweenness_centrality(nx.to_undirected(PBC_lattice), 
                                                            normalized=False)

for edge in PBC_lattice.edges:
    if edge not in list(edge_rw_centr.keys()):
        edge_rw_centr[edge] = edge_rw_centr[(edge[1],edge[0])]
        del edge_rw_centr[(edge[1],edge[0])]

# for edge in edge_rw_centr.keys():
#     edge_rw_centr[edge] *= 200

nx.set_edge_attributes(PBC_lattice, edge_rw_centr, name = 'edge_visits')

grad, sol, har, pot, div = hodge_decomposition(
    PBC_lattice, 'edge_visits')

plot_hodge(PBC_lattice, grad, sol, har, pot, div, pos)   
#%%

edge_rw_centr = nx.edge_current_flow_betweenness_centrality(nx.to_undirected(cv), 
                                                            normalized=False)

for edge in cv.edges:
    if edge not in list(edge_rw_centr.keys()):
        edge_rw_centr[edge] = -edge_rw_centr[(edge[1],edge[0])]
        del edge_rw_centr[(edge[1],edge[0])]

# for edge in edge_rw_centr.keys():
#     edge_rw_centr[edge] *= 200

nx.set_edge_attributes(cv, edge_rw_centr, name = 'edge_visits')

grad, sol, har, pot, div = hodge_decomposition(
    cv, 'edge_visits')

plot_hodge(cv, grad, sol, har, pot, div, ind_to_pos)  
#%%DISCRETE NON ABSORBING NEWMANN: ANALYTIC DISCRETE RW

''' first we build the transition matrix'''
# building the transition rate matrix
#dict of matrix ind to node index
PBC_und = PBC_lattice.to_undirected()
deg_arr = []
#mapping of each node to its index in the transition rate array
inode_to_iarr = {node:i for i, node in enumerate(PBC_und.nodes)}
iarr_to_inode = {inode_to_iarr[i]:i for i in inode_to_iarr.keys()}
#building the array
trans_matrix = np.zeros((len(PBC_und.nodes), len(PBC_und.nodes)))
for node in PBC_und.nodes:
    #Degree of the departure node
    k = len(list(PBC_und.neighbors(node)))
    deg_arr.append(1/k)
    #for each neighbour we calculate the transition rate probability of going 
    #from node to final as v/(dist*k)
    for final in PBC_und.neighbors(node):
        i, j = inode_to_iarr[node], inode_to_iarr[final]
        trans_matrix[j][i] = 1/k
deg_arr = np.array(deg_arr)
#%% raising the array to the amount of steps
steps = 200
edge_discr_flow = dict.fromkeys(PBC_lattice.edges, 0)
trans_probabilities = np.zeros(np.shape(trans_matrix))
for n in range(0, steps):
# the probability of starting at i and endng at j after n steps is each element 
# of this matrix
    trans_probabilities += np.linalg.matrix_power(trans_matrix, n)
#now we sum all the probabilities of finishing at node j and multiply this value
# by the degree of j in order to find the expected value of transitions from j 
#to i at step t
visit_prob = trans_probabilities.dot(np.ones_like(deg_arr))
probs = np.copy(visit_prob)
prob_node = {iarr_to_inode[i]: prob for i, prob in enumerate(probs)}
visit_prob *= deg_arr
for edge in PBC_lattice.edges:
    i, j = inode_to_iarr[edge[0]], inode_to_iarr[edge[1]]
    edge_discr_flow[edge] -= (visit_prob[j] - visit_prob[i])*20
        
nx.set_edge_attributes(PBC_lattice, edge_discr_flow, name = 'edge_visits')
grad_discr, sol_discr, har_discr, pot_discr, div_discr = hodge_decomposition(
    PBC_lattice, 'edge_visits')
plot_hodge(PBC_lattice, grad_discr, sol_discr, har_discr, pot_discr, div_discr,
           pos)  

for i in pot.keys():
    print(pot[i], pot_discr[i])

#%% BOXPLOT OF POTENTIALS
Dt = 200
# for Dt in range(95,105, 5):
n_walk = 20
boxpl_ls_p = []
boxpl_ls_g = []
for i in range(20):
    walked_PBC_discr = digraph_walkers(PBC_lattice, Dt, n_walk)
    
    grad, sol, har, pot, div = hodge_decomposition(
        walked_PBC_discr, 'edge_visits')

    pot_diff = [(pot_discr[node]-pot[node]) for node in PBC_lattice.nodes]
    grad_diff = [(grad_discr[edge]-grad[edge]) for edge in PBC_lattice.edges]
    boxpl_ls_p += pot_diff
    boxpl_ls_g += grad_diff
boxpl = []
boxpl.append(boxpl_ls_p)
boxpl.append(boxpl_ls_g)

plt.figure()
plt.boxplot(boxpl, labels =  ['Potential discrepancy'])#, 'Gradient discrepancy'])
plt.ylabel('Relative discrepancy from predicted random walk')
plt.tight_layout()
for i in pot.keys():
    print(pot[i], pot_discr[i])
#%% HISTOGRAM OF THE DATA
from scipy.stats import norm
# Fit a normal distribution to the data:
mu, std = norm.fit(boxpl[0])
plt.xlabel('Discrepancy in the Potentials')
plt.title("Fitted discrepancy distribution: mu = %.2f,  std = %.2f" % (mu, std))
# Plot the histogram.
plt.hist(boxpl[0], bins=50, density=True, alpha=0.6, color='b')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)


plt.show()

#%%
v = 1
Dt = 20
# for Dt in range(95,105, 5):
n_walk = 20
follow_node = 69
for i in range(40):
    walked_PBC, path = periodic_walkers(PBC_lattice, Dt, v, pos, n_walk, N, follow_node)
    
    grad, sol, har, pot, div = hodge_decomposition(
        walked_PBC, 'edge_visits')
    # print(structural_ratios(PBC_lattice))
    # plot_hodge(walked_PBC, grad, sol, har, pot, div, pos)