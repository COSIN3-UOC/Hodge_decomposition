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
import hodge_dec.hodge_decomposition as hd

# %% WALKERS ON A DIRECTED GRAPH (walk during amount of steps)
from collections import Counter
'''
moves "n_walk" walkers/node in a directed graph "G" "steps" times. The initial positions of the
walkers are randomly selected amongst the nodes. The function returns the graph
with the amount of times any walker has passed through each node as attributes
to the nodes

'''


def digraph_walkers(G, steps, n_walk):

    # overall_node_weights = {node: 0 for node in list(G.nodes)}
    #total occupation of each node at each time step
    occupation_t = [dict.fromkeys(G.nodes, 0) for i in range(steps+1)]
    occupation_t[0] = dict.fromkeys(G.nodes, n_walk)
    #overall edge flow
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # path = []
        # let's add a list of n_walk random walkers
        initial_nodes = list(G.nodes)
        # initial_nodes += initial_nodes
        # path stores the position of each walker in each time step
        # path.append(initial_nodes)
        # weights counts the amount of times a walker has visited each node
        # weights = {i_node: path[0].count(i_node) for i_node in list(G.nodes)}
        edge_weights = {i_edge: 0 for i_edge in list(G.edges)}
        # now, let's move the walker to a connected node
        # first let's see the available nodes to go to
        for step in range(1, steps+1):
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
            
            #counting the occupation of each node at step 
            occ_dict = Counter(final_nodes)
            #we take the dict corresponding to step t and sum the occupation of
            #each node at time t
            total_occ_step = occupation_t[step]
            for node in occ_dict.keys():
                total_occ_step[node] += occ_dict[node]
            #overwrite the occupancy
            occupation_t[step] = total_occ_step
            
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
        # for node in G.nodes:
        #     overall_node_weights[node] += weights[node]
        # set node value as the number of visits of each value
        # print(weights)
    # nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G, occupation_t)
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

def node_walkers(G, Dt, v, pos, n_walk):
    '''
    Function that puts n = len(list(G.nodes)) random walkers in a Digraph
    transitable in both directions and moves them during Dt time. 

    Parameters
    ----------
    G : nx.DiGraph
        Input graph to walk on.
    Dt : Float
        Maximum moving time.
    v : Float
        velocity of the walkers.
    pos : Dictionary
        dictionary with positions of each node in the DiGraph {node:(posx,pos_y)}.
    n_walk : Integer
        number of realisations of the random walk during Dt max time'.

    Returns
    -------
    G: nx.DiGraph with the edge flows in an attribute called 'edge_visits'.

    '''
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
    overall_edge_weights = dict.fromkeys(G.edges, 0)
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {i: i for i in G.nodes}

        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        edge_weights = dict.fromkeys(G.edges, 0)
        # first let's difine the time used for each walker
        time = dict.fromkeys(G.nodes, 0)
        active_walkers = dict.fromkeys(G.nodes, True)
        # 1 step moves (or tries to move according to time left) all the walkers
        # through a neighboring edge
        for step in range(max_steps):
            # list of active walkers indices
            if not step%10:
                print(Counter(active_walkers.values()))
            active_ls = [walk_ind for walk_ind, value in active_walkers.items()
                         if value]
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
                    possible_nodes = [k for k in goto_nodes.keys()] #if
                                      #(goto_nodes[k]+time[walker]) <= Dt]
                    # random selection between the final nodes
                    sel = random.choice(possible_nodes)
                    #if the choice surpasses the allowed time we deactivate the walker
                    if time[walker] + goto_nodes[sel]>Dt:
                        # deactivate the walker that has finished with a 50% chance
                        # if random.choice([0,1]) == 1:
                        #     pass
                        # else:
                        active_walkers[walker] = False
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
        if i==0:
            delta_t[(i,i+N-1)] = lattice_gap 
        else:
            delta_t[(i*N, i*N+N-1)] = lattice_gap
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
            active_ls = [walk_ind for walk_ind, value in active_walkers.items()
                         if value]
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
                        # deactivate the walker that has finished with a 50% chance
                        # if random.choice([0,1]) == 1:
                        #     pass
                        # else:
                        active_walkers[walker] = False
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

# %% Node-centric random walk

def node_centric(G, Dt, v, pos, n_walk):
    '''
    Function that performs a node centric RW. The function puts n = 
    len(list(G.nodes)) random walkers in a Digraph transitable in both 
    directions and moves them during Dt time. 

    Parameters
    ----------
    G : nx.DiGraph
        Input graph to walk on.
    Dt : Float
        Maximum moving time.
    v : Float
        velocity of the walkers.
    pos : Dictionary
        dictionary with positions of each node in the DiGraph {node:(posx,pos_y)}.
    n_walk : Integer
        number of realisations of the random walk during Dt max time'.

    Returns
    -------
    G: nx.DiGraph with the edge flows in an attribute called 'edge_visits'.

    '''
    path_tot = []
    # calculating the time intrinsec to each edge
    delta_t = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
            for edge in G.edges}
    nx.set_edge_attributes(G, delta_t, name='dt')
    # dict of 1/lambda
    scales = {}
    for node in G.nodes:
        neigh_times = [delta_t[(node,neigh)] if (node, neigh) in list(G.edges) else 
                        delta_t[(neigh,node)] for neigh in nx.all_neighbors(G, node)]
        scales[node] = np.mean(neigh_times)
    # upper bound for steps
    min_edge_time = min(delta_t.values())
    max_steps = int(Dt/min_edge_time)
    print('maximum steps allowed each realisation', max_steps)
    # overall_edge_weights account for all the edge passings summed for
    # n_walk iterations
    overall_edge_weights = dict.fromkeys(G.edges, 0)
    for i in range(0, n_walk):
        path = []
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {i: i for i in G.nodes}
        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        edge_weights = dict.fromkeys(G.edges, 0)
        # first let's difine the time used for each walker
        time = dict.fromkeys(G.nodes, 0)
        path.append([list(pos_walkers.values()), list(time.values())])
        active_walkers = dict.fromkeys(G.nodes, True)
        # 1 step moves (or tries to move according to time left) all the walkers
        # through a neighboring edge
        for step in range(max_steps):
            # list of active walkers indices
            # if not step%10:
            #     print(Counter(active_walkers.values()))
            active_ls = [walk_ind for walk_ind, value in active_walkers.items()
                         if value]
            if len(active_ls) == 0:  # all walkers have finished
                print('all walkers finished, steps needed for this realisation '
                      + str(step)+'/'+str(max_steps))
                break

            # list of neighboring nodes of each walker
            neighbors = {walker:list(nx.all_neighbors(G, pos_walkers[walker]))
                         for walker in active_ls}
            # now randomly move the random walker to another node only if time+edge_time
            # is < Dt
            for walker, possible_nodes in neighbors.items():
                # time of each possible edge for a given node
                interval = np.random.exponential(scales[pos_walkers[walker]])
                if time[walker] + interval <= Dt:
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
                    time[walker] += interval
                else:
                    # desactivate the walker that has finished
                    active_walkers[walker] = False
                    # print(Counter(active_walkers.values()))
            path.append([list(pos_walkers.values()), list(time.values())])
        path_tot.append(path)
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G, path_tot)

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
    print('error in node potentials', lap.dot(pot_field)+div_arr)
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
        print('error curl component',
              curl_op.dot(np.transpose(curl_op)).dot(tri_pot)-rot)
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
    
    #mean neigh dists for each node
    
    dists = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))
            for edge in walk_graph.edges}
    
    nx.set_edge_attributes(walk_graph, dists, 'length')
    
    mean_dists = {}
    for node in walk_graph.nodes:
        adj_edg = list(walk_graph.in_edges(node)) + \
        list(walk_graph.out_edges(node))
        
        avg_dst = 0
        for e in adj_edg:
            avg_dst += walk_graph[e[0]][e[1]]['length']
        
        avg_dst /= len(adj_edg)
        mean_dists[node] = avg_dst
    
        
    
    
    percentile = np.percentile(g_field, 95)
    percentile_pot = np.percentile(list(pot.values()), 95)
    percentile_dist = np.percentile(list(mean_dists.values()), 95)

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
    plt.subplots(1, 2, figsize=(16, 5))
    plt.subplot(121)
    plt.title('Original Graph 100%', fontsize = 20)

    # color_p = np.abs(np.array(list(edge_graph.values())))
    color_p = np.array(list(edge_graph.values()))
    # colors = np.linspace(0, percentile)
    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_div = list(mean_dists.values())
    colors_div = range(int(min(color_div)), int(max(color_div)))
    cmap_div = plt.cm.PRGn
    # vmin_div = min(colors_div)
    # vmax_div = round(max(color_div))
    vmax_div = percentile_dist
    vmin_div = min(color_div)


    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=6, 
                           # node_color='#D3D3D3')
                            node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                            vmax=vmax_div)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_p,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 6)

    sm = plt.cm.ScalarMappable(cmap=cmap, 
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\omega_{tot}$', fontsize = 18)
    cbar.ax.tick_params(labelsize=18)


    sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                                  vmax=vmax_div))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Avg. surrounding edge lengths', fontsize = 18)
    cbar2.ax.tick_params(labelsize=18)
    plt.subplot(122)
    
    #GRADIENT COMPONENT
    
    color_g = np.abs(np.array(list(grad_comp.values())))
    # plotting edges with color gradient

    color_pot = list(pot.values())
    cmap_pot = plt.cm.PRGn
    # vmax_pot = max([int(np.max(color_pot)), int(np.abs(np.min(color_pot)))])
    vmax_pot = percentile_pot
    vmin_pot = -vmax_pot
    # vmax_pot = np.max(color_pot)
    # vmin_pot = np.min(color_pot)

    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)
    plt.title('Gradient component ' + str(round(weight_g*100, 1))+'%', fontsize = 20)
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None,
                            node_size=6, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 5, node_size = 6)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_g\right|$', fontsize = 18)
    cbar.ax.tick_params(labelsize=18)
    
    sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                                  vmax=vmax_pot))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Node potentials', fontsize = 18)
    cbar2.ax.tick_params(labelsize = 18)
    
    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_s = np.abs(np.array(list(sol_comp.values())))
    # plt.subplot(223)
    # plt.title('Solenoidal Component ' + str(round(weight_s*100, 1))+'%', fontsize = 20)
    # nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
    #                        node_color='#D3D3D3')
    # nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_s,
    #                        edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
    #                        arrowsize = 5, node_size = 4)


    # sm = plt.cm.ScalarMappable(
    #     cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    # sm._A = []
    # cbar = plt.colorbar(sm)
    # cbar.set_label(r'$\left|\omega_s\right|$', fontsize = 18)
    # cbar.ax.tick_params(labelsize=18)
    
    
    # colors = np.linspace(0, percentile)
    # cmap = plt.cm.Oranges
    # vmin = min(colors)
    # vmax = max(colors)

    # color_h = np.array(list(har_comp.values()))
    # plt.subplot(224)
    # plt.title('Harmonic Component ' + str(round(weight_h*100, 1))+'%', fontsize = 20)
    # nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
    #                        node_color='#D3D3D3')
    # nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_h,
    #                        edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
    #                        arrowsize = 5, node_size = 4)
    # sm = plt.cm.ScalarMappable(
    #     cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    # sm._A = []
    # cbar = plt.colorbar(sm)
    # cbar.set_label(r'$\left|\omega_h\right|$', fontsize = 18)
    # cbar.ax.tick_params(labelsize=18)
    plt.tight_layout()
    # plt.savefig('/Users/robertbenassai/Documents/UOC/figs/lattice_evolution/lattice_dt'+str(Dt)+'.png')
    plt.show()


# %% MAPS, GEOPANDAS
import time
import momepy
def distr_to_nx(distr_ind:int, path_bcn: str, path_distr: str):
    
    # edges and nodes of the whole graph
    bcn_edges = gpd.read_file(path_bcn, crs="EPSG:25831")


    bcn_nodes = gpd.read_file('/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/nodes_clean_net_willum.shp',
                          crs="EPSG:25831")
    
    if distr_ind == 0: # graph of the entire city
    
        # initialize graph of the district
        distr = nx.DiGraph()
        distr.add_nodes_from(bcn_nodes['uid'])
        ind_to_pos = {bcn_nodes.loc[i,'uid']:(bcn_nodes.loc[i,'geometry'].x,
                                              bcn_nodes.loc[i,'geometry'].y) 
                      for i in range(len(bcn_nodes))}
    
    
        for i, j in zip(bcn_edges['i'], bcn_edges['j']):
            
            distr.add_edge(int(i), int(j))
        
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
                distr.add_edge(int(i), int(j))
                
    
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
    
    #plot of components
    
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
    
    #relabeling nodes to ascending order
    dict_relabel = {node: i for i, node in enumerate(distr.nodes)}
    distr = nx.relabel_nodes(distr, dict_relabel, copy = True)    
    ind_to_pos_upd = {dict_relabel[node]: ind_to_pos[node] for node in 
                      dict_relabel.keys()}
    return(distr, ind_to_pos_upd)

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
    
#%%DISCRETE NON ABSORBING NEWMANN: ANALYTIC DISCRETE RW
def build_trans_matrix(G):
    '''
    Builds the transition  matrix or generator matrix for an unweighted 
    random walk in a topological graph Tij = 1/k_j

    Parameters
    ----------
    G : nx.DiGraph
        Digraph from which the transition rates will be computed.

    Returns
    -------
    trans_matrix : npumpy ndarray.
                  transition matrix

    '''
    G_und = G.to_undirected()
    # building the transition rate matrix
    #dict of matrix ind to node index
    deg_arr = []
    #mapping of each node to its index in the transition rate array
    inode_to_iarr = {node:i for i, node in enumerate(G_und.nodes)}
    #building the array
    trans_matrix = np.zeros((len(G_und.nodes), len(G_und.nodes)))
    for node in G_und.nodes:
        #Degree of the departure node
        k = len(list(G_und.neighbors(node)))
        deg_arr.append(1/k)
        #for each neighbour we calculate the transition rate probability of going 
        #from node to final as v/(dist*k)
        for final in G_und.neighbors(node):
            i, j = inode_to_iarr[node], inode_to_iarr[final]
            trans_matrix[j][i] = 1/k
    print(np.sum(trans_matrix, axis = 0))
    return(trans_matrix)
#%% raising the array to the amount of steps
def discrete_rw_edge_flow(G, trans_matrix, steps, n_walk):
    '''
    Calculates the edge flows of a discrete time random walk on a lattice.

    Parameters
    ----------
    G : nx.DiGraph
        Directed graph in which the random walker moves.
    trans_matrix : Numpy ndarray
        Transition matrix.
    steps : int
        total steps or jumps of the walker.
    n_walk : int
        number of repetitions of each walk.

    Returns
    -------
    G : the input graph with an edge attribute called 'Edge_visits' corresponding
    to the RW edge flow.

    '''
    edge_discr_flow = dict.fromkeys(G.edges, 0)
    trans_probabilities = np.zeros(np.shape(trans_matrix))
    prob_evo = []
    iarr_to_inode = {i:node for i, node in enumerate(G.nodes)}

    deg_arr =np.array([1/len(list(G.successors(node))+list(G.predecessors(node)))
                       for node in G.nodes])
    for n in range(0, steps+1):
    # the probability of starting at i and endng at j after n steps is each element 
    # of this matrix
        trans_probabilities += np.linalg.matrix_power(trans_matrix, n)
        probs = np.linalg.matrix_power(trans_matrix, n).dot(np.ones_like(deg_arr))
    #now we sum all the probabilities of finishing at node j and multiply this value
    # by the degree of j in order to find the expected value of transitions from j 
    #to i at step t
        prob_node = {iarr_to_inode[i]: prob for i, prob in enumerate(probs)}
        prob_evo.append(prob_node)
    visit_prob = trans_probabilities.dot(np.ones_like(deg_arr))
    inode_to_iarr = {node:i for i, node in enumerate(G.nodes)}
    # prob_node = {iarr_to_inode[i]: prob for i, prob in enumerate(probs)}
    visit_prob *= deg_arr
    for edge in G.edges:
        i, j = inode_to_iarr[edge[0]], inode_to_iarr[edge[1]]
        edge_discr_flow[edge] += (visit_prob[i] - visit_prob[j])*n_walk
            
    nx.set_edge_attributes(G, edge_discr_flow, name = 'edge_visits')
    # stat_prob = np.linalg.matrix_power(trans_matrix, steps)
    # ini_prob = np.zeros_like(deg_arr)
    # ini_prob[0] = 1
    # probs = stat_prob.dot(ini_prob)
    # prob_node = {iarr_to_inode[i]: prob for i, prob in enumerate(probs)}
    return(G, prob_evo)

#%% CONTINUUM NEWMANN WITHOUT ABSORBINGS

'''first we build the transition rates matrix'''

# building the transition rate matrix
#dict of matrix ind to node index
def build_trans_rates_matrix(G, pos, v):
    '''
    Builds the transition rates matrix or generator matrix for an unweighted 
    random walk in a topological graph Rij = v/(k_j*Dx_ij)

    Parameters
    ----------
    G : nx.DiGraph
        Digraph from which the transition rates will be computed.
    pos : Dict
        dictionary of node:(pos_x, pos_y).
    v : float
        Velocity of the random walker
    Returns
    -------
    trans_rates : npumpy ndarray.
                  transition rates matrix

    '''
    G_und = G.to_undirected()
    
    #distances of the edges in the graph
    dists = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
            for edge in G_und.edges}
    
    # ONLY FOR PBC_LATTICE ----------------------------------------------------
    # lattice_gap = 1
    # for i in range(N):
    #     dists[(i, i+N*N-N)] = lattice_gap
    #     dists[(i+N*N-N, i)] = lattice_gap
    #     if i==0:
    #         dists[(i,i+N-1)] = lattice_gap 
    #         dists[(i+N-1, i)] = lattice_gap 
    #     else:
    #         dists[(i*N, i*N+N-1)] = lattice_gap
    #         dists[(i*N+N-1, i*N)] = lattice_gap
    # -------------------------------------------------------------------------
    #mapping of each node to its index in the transition rate array
    inode_to_iarr = {node:i for i, node in enumerate(G_und.nodes)}
    #building the array
    trans_rates = np.zeros((len(G_und.nodes), len(G_und.nodes)))
    for node in G_und.nodes:
        #Degree of the departure node
        k_deg = len(list(G_und.neighbors(node)))
        # k = np.sum(np.array([dists[(node, neigh)] if (node, neigh) in list(G_und.edges)
        #                       else dists[(neigh, node)]
        #                       for neigh in G_und.neighbors(node)]))
                             
        #for each neighbour we calculate the transition rate probability of going 
        #from node to final as v/(dist*k)
        neigh_dists = np.array([dists[(node,final)] if (node,final) in dists.keys()
                                else dists[(final,node)] for final 
                                in G_und.neighbors(node)])
        mean_dist = np.mean(neigh_dists)
        for final in G_und.neighbors(node):
            i, j = inode_to_iarr[node], inode_to_iarr[final]
            # if (node,final) in dists.keys():
            #     dist = dists[(node, final)]
            # else:
            #     dist = dists[(final, node)]
            trans_rates[j][i] = 1/(k_deg*mean_dist)
    # the diagonal is -sum of the off diag elements of the col
    i,j = np.indices(trans_rates.shape)
    diagonal = np.sum(trans_rates, axis=0)
    trans_rates[i==j] = -diagonal
    return(trans_rates)


#%% SYSTEM OF ODES
from scipy.integrate import odeint

def dp_dt(p, t, R, dummy):
    '''
    System of linear differential equations for the probabilities
    
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

def solve_continuous_rw_flow(G:nx.DiGraph, trans_rates:np.array, Dt:float, n_walk:int):
    '''
    Solves the linear differential equations for the node probabilities and 
    finds the ende flows of continuous random walks on a topological graph.

    Parameters
    ----------
    G : nx.DiGraph
        Directed graph in which the RW takes place.
    trans_rates : Numpy ndarray
        Array of transition rates.
    Dt : float
        Total walking time.
    n_walk : int
        number of repetitions of the walks.

    Returns
    -------
    nx.DiGraph with an edge attribute called 'edge_visits' with the edge random
    walk flow.

    '''
    #mapping of each node to its index in the transition rate array
    G_und = nx.to_undirected(G)
    inode_to_iarr = {node:i for i, node in enumerate(G_und.nodes)}
    #array of flows through the edges
    net_flow = np.zeros(np.shape(trans_rates))
    #initial probabilities
    # for initial in range(len(list(G.nodes))):
    # p0 = np.ones(len(G.nodes))
    p0 = np.ones(len(G.nodes))
    
    t = np.linspace(start=0, stop=Dt,num=10000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    
    for j, p_j in enumerate(flows):
        
        for i in range(len(G.nodes)):
            net_flow[i,j] += flows[i]*trans_rates[j][i] - p_j*trans_rates[i][j]

    net_flow *= n_walk
    cont_new = {edge: net_flow[inode_to_iarr[edge[0]]][inode_to_iarr[edge[1]]] 
                for edge in G.edges}
    nx.set_edge_attributes(G, cont_new, name = 'edge_visits')
    return(G, sol)

# %% EIXAMPLE (02)
'''L'EIXAMPLE'''
distr_ind = 2
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

#THEO
trans_rates_eix = build_trans_rates_matrix(eixample, ind_to_pos, v)
walk_eix, _ = solve_continuous_rw_flow(eixample.copy(), trans_rates_eix, Dt, n_walk)

#%% simulation
# walk_eix = node_walkers(eixample, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_eix = hd.reverse_negative_edges(walk_eix)
grad_eix, sol_eix, har_eix, pot_eix, div_eix = hodge_decomposition(
    rev_eix, 'edge_visits')
# %%
with open('eixample_dec_dt_15_connected.txt', 'w') as eix_file:
    eix_file.write(str(grad_eix)+'\n')
    eix_file.write(str(sol_eix)+'\n')
    eix_file.write(str(har_eix)+'\n')
    eix_file.write(str(pot_eix)+'\n')
    eix_file.write(str(div_eix)+'\n')
eix_file.close()
# %%
plot_hodge(rev_eix, grad_eix, sol_eix, har_eix, pot_eix, div_eix, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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
trans_rates_cv = build_trans_rates_matrix(cv, ind_to_pos, v)
walk_cv, _ = solve_continuous_rw_flow(cv.copy(), trans_rates_cv, Dt, n_walk)
# walk_cv = node_walkers(cv, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
#%%
rev_cv = hd.reverse_negative_edges(walk_cv)
grad_cv, sol_cv, har_cv, pot_cv, div_cv = hd.hodge_decomposition(
    rev_cv, 'edge_visits')
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
#plot
plt.xlabel('Node potential')
plt.ylabel('Frequency')
plt.hist(pot_field, bins = 50, edgecolor ='black')

#%%

plot_hodge(rev_cv, grad_cv, sol_cv, har_cv, pot_cv, div_cv, ind_to_pos)

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
rev_cv = hd.reverse_negative_edges(walk_cv)
grad_cv, sol_cv, har_cv, pot_cv, div_cv = hodge_decomposition(
    rev_cv, 'edge_visits')
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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20
#THEO
trans_rates_sts = build_trans_rates_matrix(sts_mj, ind_to_pos, v)
walk_sts_mj, _ = solve_continuous_rw_flow(sts_mj.copy(), trans_rates_sts, Dt, n_walk)

# walk_sts_mj = node_walkers(sts_mj, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
#%%
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(sts_mj)

print(st_ratio_g, st_ratio_s, st_ratio_h)
# %%
rev_sts_mj = hd.reverse_negative_edges(walk_sts_mj)
grad_sts_mj, sol_sts_mj, har_sts_mj, pot_sts_mj, div_sts_mj = hodge_decomposition(
    rev_sts_mj, 'edge_visits')
# %%
with open('sts_mj_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_sts_mj)+'\n')
    file.write(str(sol_sts_mj)+'\n')
    file.write(str(har_sts_mj)+'\n')
    file.write(str(pot_sts_mj)+'\n')
    file.write(str(div_sts_mj)+'\n')
file.close()
# %%
plot_hodge(rev_sts_mj, grad_sts_mj, sol_sts_mj, har_sts_mj, pot_sts_mj, div_sts_mj, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20

#THEO
trans_rates_corts = build_trans_rates_matrix(corts, ind_to_pos, v)
walk_corts, _ = solve_continuous_rw_flow(corts.copy(), trans_rates_corts, Dt, n_walk)
# walk_corts = node_walkers(corts, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_corts = hd.reverse_negative_edges(walk_corts)
grad_corts, sol_corts, har_corts, pot_corts, div_corts = hodge_decomposition(
    rev_corts, 'edge_visits')
# %%
with open('corts_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_corts)+'\n')
    file.write(str(sol_corts)+'\n')
    file.write(str(har_corts)+'\n')
    file.write(str(pot_corts)+'\n')
    file.write(str(div_corts)+'\n')
file.close()
# %%
plot_hodge(rev_corts, grad_corts, sol_corts, har_corts, pot_corts, div_corts, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20

#THEO
trans_rates_sarr = build_trans_rates_matrix(sarr, ind_to_pos, v)
walk_sarr, _ = solve_continuous_rw_flow(sarr.copy(), trans_rates_sarr, Dt, n_walk)
# walk_sarr = node_walkers(sarr, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_sarr = hd.reverse_negative_edges(walk_sarr)
grad_sarr, sol_sarr, har_sarr, pot_sarr, div_sarr = hodge_decomposition(
    rev_sarr, 'edge_visits')
# %%
with open('sarr_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_sarr)+'\n')
    file.write(str(sol_sarr)+'\n')
    file.write(str(har_sarr)+'\n')
    file.write(str(pot_sarr)+'\n')
    file.write(str(div_sarr)+'\n')
file.close()
# %%
plot_hodge(rev_sarr, grad_sarr, sol_sarr, har_sarr, pot_sarr, div_sarr, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

#THEO
trans_rates_gracia = build_trans_rates_matrix(gracia, ind_to_pos, v)
walk_gracia, _ = solve_continuous_rw_flow(gracia.copy(), trans_rates_gracia, Dt, n_walk)

# walk_gracia = node_walkers(gracia, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_gracia = hd.reverse_negative_edges(walk_gracia)
grad_gracia, sol_gracia, har_gracia, pot_gracia, div_gracia = hodge_decomposition(
    rev_gracia, 'edge_visits')
# %%
with open('gracia_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_gracia)+'\n')
    file.write(str(sol_gracia)+'\n')
    file.write(str(har_gracia)+'\n')
    file.write(str(pot_gracia)+'\n')
    file.write(str(div_gracia)+'\n')
file.close()
# %%
plot_hodge(rev_gracia, grad_gracia, sol_gracia, har_gracia, pot_gracia, div_gracia, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

#THEO
trans_rates_horta = build_trans_rates_matrix(horta, ind_to_pos, v)
walk_horta, _ = solve_continuous_rw_flow(horta.copy(), trans_rates_horta, Dt, n_walk)

# walk_horta = node_walkers(horta, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_horta = hd.reverse_negative_edges(walk_horta)

grad_horta, sol_horta, har_horta, pot_horta, div_horta = hodge_decomposition(
    rev_horta, 'edge_visits')
# %%
with open('horta_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_horta)+'\n')
    file.write(str(sol_horta)+'\n')
    file.write(str(har_horta)+'\n')
    file.write(str(pot_horta)+'\n')
    file.write(str(div_horta)+'\n')
file.close()
# %%
plot_hodge(rev_horta, grad_horta, sol_horta, har_horta, pot_horta, div_horta, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20

#THEO
trans_rates_noub = build_trans_rates_matrix(noub, ind_to_pos, v)
walk_noub, _ = solve_continuous_rw_flow(noub.copy(), trans_rates_noub, Dt, n_walk)

# walk_noub = node_walkers(noub, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_noub = hd.reverse_negative_edges(walk_noub)

grad_noub, sol_noub, har_noub, pot_noub, div_noub = hodge_decomposition(
    rev_noub, 'edge_visits')
# %%
with open('noub_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_noub)+'\n')
    file.write(str(sol_noub)+'\n')
    file.write(str(har_noub)+'\n')
    file.write(str(pot_noub)+'\n')
    file.write(str(div_noub)+'\n')
file.close()
# %%
plot_hodge(rev_noub, grad_noub, sol_noub, har_noub, pot_noub, div_noub, ind_to_pos)

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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

#THEO
trans_rates_st_and = build_trans_rates_matrix(st_and, ind_to_pos, v)
walk_st_and, _ = solve_continuous_rw_flow(st_and.copy(), trans_rates_st_and, Dt
                                          , n_walk)

# walk_st_and = node_walkers(st_and, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_st_and = hd.reverse_negative_edges(walk_st_and)

grad_st_and, sol_st_and, har_st_and, pot_st_and, div_st_and = hodge_decomposition(
    rev_st_and, 'edge_visits')
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
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

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

#THEO
trans_rates_st_mart = build_trans_rates_matrix(st_mart, ind_to_pos, v)
walk_st_mart, _ = solve_continuous_rw_flow(st_mart.copy(), trans_rates_st_mart,
                                           Dt, n_walk)

# walk_st_mart = node_walkers(st_mart, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))
# %%
rev_st_mart = hd.reverse_negative_edges(walk_st_mart)

grad_st_mart, sol_st_mart, har_st_mart, pot_st_mart, div_st_mart = hodge_decomposition(
    rev_st_mart, 'edge_visits')
# %%
with open('st_mart_dec_dt_15_connected.txt', 'w') as file:
    file.write(str(grad_st_mart)+'\n')
    file.write(str(sol_st_mart)+'\n')
    file.write(str(har_st_mart)+'\n')
    file.write(str(pot_st_mart)+'\n')
    file.write(str(div_st_mart)+'\n')
file.close()
# %%
plot_hodge(rev_st_mart, grad_st_mart, sol_st_mart, har_st_mart, pot_st_mart, div_st_mart, ind_to_pos)

# %%
dict_ls = []
with open('st_mart_dec.txt', 'r', newline='\n') as file:
    data = file.readlines()
    for i in data:
        dict_ls.append(eval(i))
file.close()


#%% AL THE CITY

'ALL THE CITY'

distr_ind = 0
path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'

bcn_graph, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)

#%%
start_time = time.time()
v = 1.42
Dt = 15*60
n_walk = 20

#THEO
trans_rates_bcn = build_trans_rates_matrix(bcn_graph, ind_to_pos, v)
walk_bcn, _ = solve_continuous_rw_flow(bcn_graph.copy(), trans_rates_bcn,
                                           Dt, n_walk)

# walk_st_mart = node_walkers(st_mart, Dt, v, ind_to_pos, n_walk)
print("--- %s seconds ---" % (time.time() - start_time))

rev_bcn = hd.reverse_negative_edges(walk_bcn)

grad_bcn, sol_bcn, har_bcn, pot_bcnt, div_bcn = hodge_decomposition(
    rev_bcn, 'edge_visits')

plot_hodge(rev_bcn, grad_bcn, sol_bcn, har_bcn, pot_bcnt, div_bcn, ind_to_pos)

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
Dt = 100
# for Dt in range(95,105, 5):
n_walk = 20
follow_node = 69
walked_PBC, path = periodic_walkers(PBC_lattice.copy(), Dt, v, pos, n_walk, N, follow_node)

grad_exp_c, sol_exp_c, har_exp_c, pot_exp_c, div_exp_c = hodge_decomposition(
    walked_PBC, 'edge_visits')
print(structural_ratios(PBC_lattice))
plot_hodge(walked_PBC, grad_exp_c, sol_exp_c, har_exp_c, pot_exp_c, div_exp_c, pos)
#%% PLOTTING PATHS
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
plt.figure()
plt.boxplot(boxpl, labels =  ['Gradient discrepancy', 'Harmonic discrepancy'])
plt.ylabel('Edge discrepancy from predicted component')



#%% FINDING THE EIGENVECS AND EIGENVALS. The eigenvector with eigenvalue 1 is 
#the stationary state in theory
v = 1
trans_rates = build_trans_rates_matrix(PBC_lattice, pos, v)
#%%
eigval, eigvecs = scipy.linalg.eig(trans_rates)
print(np.where(np.real(np.round(eigval, 5)) == 0))
stationary_node_prob = np.real(eigvecs)[:,np.where(np.real(np.round(eigval, 2)) == 0)[0]]
print(Counter(np.transpose(stationary_node_prob)[0]))

#%%
Dt = 100
n_walk = 20
lattice_theo, solution_th = solve_continuous_rw_flow(PBC_lattice.copy(), trans_rates, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = \
hodge_decomposition(lattice_theo, 'edge_visits')
#pot_cont_new = {iarr_to_inode[i]: prob for i, prob in enumerate(sol[-1][:])}
#%%
plot_hodge(lattice_theo, g_cont_new, s_cont_new, h_cont_new, pot_cont_new, 
           div_cont_new, pos)

#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS CTRW
pot_list = [(pot_exp_c[i], pot_cont_new[i]) for i in pot_exp_c.keys()]
x, y = zip(*pot_list)
    
a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
pot_linsp = np.linspace(-40, 20, 50)
plt.figure()
plt.xlabel('Simulated Potential')
plt.ylabel('Analytical Potential')
plt.scatter(x, y, s = 10)
plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+'x +'+
         str(round(b, 2))+'\n'+r'$r^2 = $' +str(round(r_value**2, 2)))
plt.plot(x, np.array(x), c='r', label = 'y = x')
plt.legend()
plt.tight_layout()
# plt.figure()
# plt.hist(np.array(y)/np.array(x), bins = 50)
#%% Time evolution of node probabilities
plt.figure()
t = np.linspace(0, Dt, 10000)
sol = pd.DataFrame(solution_th/len(PBC_lattice.nodes))
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
for i in range(len(PBC_lattice.nodes)):
    plt.plot(t,sol[i])
plt.show()
print(8/(2*len(list(PBC_lattice.edges))), 5/(2*len(list(PBC_lattice.edges))), 
      4/(2*len(list(PBC_lattice.edges))))
plt.tight_layout()
#%% WHAT HAPPENS IF WE PERFORM A DISCRETE RW

Dt = 25
# for Dt in range(95,105, 5):
n_walk = 20

walked_PBC_discr, pr = digraph_walkers(PBC_lattice, Dt, n_walk)

grad, sol, har, pot, div = hodge_decomposition(
    walked_PBC_discr, 'edge_visits')
#%%
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

#%%
trans_matrix = build_trans_matrix(PBC_lattice)
steps = 25
n_walk = 20
PBC_lattice, pr = discrete_rw_edge_flow(PBC_lattice, trans_matrix, steps, n_walk)
grad_discr, sol_discr, har_discr, pot_discr, div_discr = hodge_decomposition(
    PBC_lattice, 'edge_visits')
#%%
# plt.hist(prob_node.values())
# n_edg = len(list(PBC_lattice.edges))
# print(4/(2*n_edg), 5/(2*n_edg), 8/(2*n_edg))
#%%
plot_hodge(PBC_lattice, grad_discr, sol_discr, har_discr, pot_discr, div_discr,
           pos)  

#%% RANDOM GEOMETRIC GRAPH DISCRETE WALK

# Use seed when creating the graph for reproducibility
G = nx.random_geometric_graph(200, 0.125, seed=896803)
# position is stored as node attribute data for random_geometric_graph
pos_c = nx.get_node_attributes(G, "pos")
# transforming the graph into a digraph
geom_d = nx.DiGraph(G)
out_edges = [edge for edge in geom_d.edges if edge[1]
    > edge[0]]  # removing problematic edges
geom_d.remove_edges_from(out_edges)

steps = 50
n_walk = 80
walk_geom = digraph_walkers(geom_d, steps, n_walk)
grad, sol, har, pot, div = hodge_decomposition(
    walk_geom, 'edge_visits')
plot_hodge(walk_geom, grad, sol, har, pot, div, pos_c)  
#analytical frucht
trans_matrix = build_trans_matrix(geom_d)
steps = 50
n_walk = 80
geom_d = discrete_rw_edge_flow(geom_d, trans_matrix, steps, n_walk)
grad_discr, sol_discr, har_discr, pot_discr, div_discr = hodge_decomposition(
    geom_d, 'edge_visits')
plot_hodge(geom_d, grad_discr, sol_discr, har_discr, pot_discr, div_discr,
           pos_c)  
#%% RANDOM GEOMETRIC GRAPH CONTINUOUS
# Use seed when creating the graph for reproducibility
G = nx.random_geometric_graph(200, 0.125, seed=896803)
# position is stored as node attribute data for random_geometric_graph
pos_c = nx.get_node_attributes(G, "pos")
# transforming the graph into a digraph
geom_d = nx.DiGraph(G)
v = 0.1
Dt = 20
n_walk = 20
#simulation 
walk_geom = node_walkers(geom_d, Dt, v, pos_c, n_walk)
grad, sol, har, pot, div = hodge_decomposition(
    walk_geom, 'edge_visits')
#analytical
trans_rates = build_trans_rates_matrix(geom_d, pos_c, v)
geom_d, solutions_geom = solve_continuous_rw_flow(geom_d, trans_rates, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = \
hodge_decomposition(geom_d, 'edge_visits')
#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS
pot_list = [(pot[i], pot_cont_new[i])for i in pot.keys()]
x, y = zip(*pot_list)
    
a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
pot_linsp = np.linspace(-10, 20, 50)
plt.figure()
plt.xlabel('Simulated Potential')
plt.ylabel('Analytical Potential')
plt.scatter(x, y, s = 10)
plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+'x +'+
         str(round(b, 2))+'\n'+r'$r^2 = $' +str(round(r_value**2, 2)))
plt.plot(x, np.array(x), c = 'red', label = 'y = x')
plt.legend()
plt.tight_layout()

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
#%% CUSTOM GRAPH TO STUDY PROBABILITIES
kite = nx.DiGraph()
kite_nodes = [i for i in range(7)]
kite_edges = [(0,1), (1,2), (2,3), (2,4), (2,5), (5,6), (3,6), (4,6), 
              (4,5), (3,4)]
kite.add_nodes_from(kite_nodes)
kite.add_edges_from(kite_edges)
pos_kite = {0:(0,0), 1:(0,1), 2:(0,2), 3:(0.2,3), 4:(0,3), 5:(-0.2, 3), 6:(0,3.5)}
nx.set_node_attributes(kite, pos_kite, 'pos')
nx.draw_networkx(kite, pos_kite)
#%%

kite = nx.DiGraph()
kite_nodes = [i for i in range(4)]
kite_edges = [(0,1), (0,2), (0,3)]
kite.add_nodes_from(kite_nodes)
kite.add_edges_from(kite_edges)
pos_kite = {0:(0,0), 1:(0,1), 2:(2,2), 3:(3,4)}
nx.set_node_attributes(kite, pos_kite, 'pos')
nx.draw_networkx(kite, pos_kite)
#%% SIMULATION DISCRETE
steps = 100
n_walk = 200
walk_kite, occupations = digraph_walkers(kite, steps, n_walk)
trans_matrix = build_trans_matrix(kite)
kite, theo_occupations = discrete_rw_edge_flow(kite, trans_matrix, steps, n_walk)
grad_discr, sol_discr, har_discr, pot_discr, div_discr = hodge_decomposition(
    kite, 'edge_visits')
#%%PLOT KITE
step = 19
theo_new = {key: round(value * n_walk) for key, value in 
            theo_occupations[step].items()}
plt.subplots(1, 2, figsize=(15, 15))
plt.subplot(121)
plt.title('Simulation')

color_prob = list(occupations[step].values())
colors_prob = range(int(min(color_prob)), int(max(color_prob)))
cmap_prob = plt.cm.Oranges
vmin_div = min(colors_prob)
vmax_div = round(max(colors_prob))

nx.draw_networkx_nodes(walk_kite, pos=pos_kite, label=None, node_size=100,
                       # node_color='#D3D3D3')
                        node_color=color_prob, cmap=cmap_prob, vmin=vmin_div,
                        vmax=vmax_div)
nx.draw_networkx_labels(walk_kite, pos = pos_kite, labels=occupations[step])
nx.draw_networkx_edges(walk_kite, pos=pos_kite, label=None, arrowsize = 5, 
                       node_size = 100)

sm2 = plt.cm.ScalarMappable(cmap=cmap_prob, norm=plt.Normalize(vmin=vmin_div,
                                                              vmax=vmax_div))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node occupation')

plt.subplot(122)

color_pot = list(theo_new.values())
cmap_pot = plt.cm.Oranges
vmax_pot = round(max(color_pot))
vmin_pot = round(min(color_pot))


plt.title('Analytical')
nx.draw_networkx_nodes(walk_kite, pos=pos_kite, label=None, node_size=100, 
                       # node_color='#D3D3D3')
                        node_color=color_pot, cmap=cmap_pot, vmin=vmin_pot,
                        vmax=vmax_div)
nx.draw_networkx_labels(walk_kite, pos = pos_kite, labels=theo_new)
nx.draw_networkx_edges(walk_kite, pos=pos_kite, label=None, arrowsize = 5, 
                       node_size = 100)

sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                              vmax=vmax_pot))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node potentials')
plt.tight_layout()
#%% CORRELATION OF OCCUPANCY

sim_occ = []
for element in occupations:
    sim_occ += list(element.values())
analyt_occ = []
for element in theo_occupations:
    analyt_occ += list(element.values())
a, b, r_value, p_value, std_err = scipy.stats.linregress(np.array(sim_occ),
                                                          np.array(analyt_occ)*n_walk)

# for element1, element2 in zip(occupations, theo_occupations[1:31]):
plt.figure()
plt.xlabel('Simulated occupancy')
plt.ylabel('Expected occupancy')
    # plt.scatter(element1.values(), np.array(list(element2.values()))*n_walk, s = 10)
plt.scatter(sim_occ, np.array(analyt_occ)*n_walk, s = 10)
plt.plot(sim_occ, np.array(sim_occ)*a+b, c = 'black', label = 'y = '+str(round(a,2))+'x +'+
              str(round(b, 2))+'\n'+r'$r^2 = $' +str(round(r_value**2, 2)))
plt.legend()
plt.tight_layout()
#%% EVOLUTION OF THE OCCUPATION
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='black', ls = '-.'),
                Line2D([0], [0], color = 'black', ls = '-')]

time = np.linspace(0, 101, 101)
colors = ['r', 'g', 'b', 'orange', 'yellow', 'black', 'purple']
plt.figure()
theo_evo = np.array([list(element.values()) for element in theo_occupations])
sim_evo = np.array([list(element.values()) for element in occupations])

for i in range(7):
    plt.plot(time, sim_evo[:,i], ls = '-.', c = colors[i])
    plt.plot(time, theo_evo[:,i]*n_walk, label = kite.degree[i], ls = '-', 
             c = colors[i])
plt.xlabel('time (steps)')
plt.ylabel('Occupation')
plt.legend(custom_lines, ['simulation', 'theoretical'], loc='upper center', 
           ncol = 2, bbox_to_anchor=(0.5, 1.1))
plt.tight_layout()

#%% SIMULATION CONTINUOUS
v = 1
Dt = 1.01
n_walk = 30
walk_kite = node_walkers(kite.copy(), Dt, v, pos_kite, n_walk)
grad_k, sol_k, har_k, pot_k, div_k = hodge_decomposition(walk_kite, 'edge_visits')
#%%
trans_rate_matrix = build_trans_rates_matrix(kite, pos_kite, v)
kite_th, solution_k = solve_continuous_rw_flow(kite, trans_rate_matrix, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = hodge_decomposition(kite_th, 'edge_visits')

#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS CTRW
pot_list = [(pot_k[i], pot_cont_new[i]) for i in pot_k.keys()]
x, y = zip(*pot_list)
    
a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
pot_linsp = np.linspace(-40, 20, 50)
plt.figure()
plt.xlabel('Simulated Potential')
plt.ylabel('Analytical Potential')
plt.scatter(x, y, s = 10)
plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+'x +'+
         str(round(b, 2))+'\n'+r'$r^2 = $' +str(round(r_value**2, 2)))
plt.plot(x, np.array(x), c='r', label = 'y = x')
plt.legend()
plt.tight_layout()
#%%
plot_hodge(kite_th, g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new,
           pos_kite)
#%%
plot_hodge(walk_kite, grad_k, sol_k, har_k, pot_k, div_k,
           pos_kite)
#%% Time evolution of node probabilities
plt.figure()
t = np.linspace(0, Dt, len(sol))
sol = pd.DataFrame(solution_k)
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
for i in range(len(kite_th.nodes)):
    plt.plot(t,sol[i], label = 'node'+str(i))
plt.legend()
plt.tight_layout()

#%% ERDOS_RENYI
Dt = 15
v = 1
n_walk = 200
erd_reny = nx.erdos_renyi_graph(50, 0.1, seed = 1000, directed=False)
erd_reny = erd_reny.to_directed()
out_edges = [edge for edge in erd_reny.edges if edge[1]
    < edge[0]]  # removing all outward edges
erd_reny.remove_edges_from(out_edges)
pos_ER = nx.spring_layout(erd_reny, seed = 1050)
nx.draw_networkx(erd_reny, with_labels = False, node_size = 20, pos = pos_ER)
#%%
walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
grad_ER, sol_ER, har_ER, pot_ER, div_ER = hodge_decomposition(walked_ER, 'edge_visits')

#%%

trans_rate_ER = build_trans_rates_matrix(erd_reny.copy(), pos_ER, v)
ER_th, solution_ER = solve_continuous_rw_flow(erd_reny, trans_rate_ER, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = hodge_decomposition(ER_th, 'edge_visits')

#%%
plot_hodge(ER_th, g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new,
           pos_ER)
#%%
plot_hodge(walked_ER, grad_ER, sol_ER, har_ER, pot_ER, div_ER,
           pos_ER)
#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS CTRW
pot_list = [(pot_ER[i], pot_cont_new[i]) for i in pot_ER.keys()]
x, y = zip(*pot_list)
    
a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
pot_linsp = np.linspace(-40, 20, 50)
plt.figure()
plt.xlabel('Simulated Potential')
plt.ylabel('Analytical Potential')
plt.scatter(x, y, s = 10)
plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+r'$\pm$'
         +str(round(std_err,2))+'x +'+str(round(b, 2))+'\n'+r'$r^2 = $' 
         +str(round(r_value**2, 2)))
plt.plot(x, np.array(x), c='r', label = 'y = x')
plt.legend()
plt.tight_layout()
#%%
grad_list = [(grad_ER[i], g_cont_new[i]) for i in grad_ER.keys()]
x, y = zip(*grad_list)
    
a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
pot_linsp = np.linspace(-40, 20, 50)
plt.figure()
plt.xlabel('Simulated gradient component')
plt.ylabel('Analytical gradient component')
plt.scatter(x, y, s = 10)
plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+'x +'+
         str(round(b, 2))+'\n'+r'$r^2 = $' +str(round(r_value**2, 2)))
plt.plot(x, np.array(x), c='r', label = 'y = x')
plt.legend()
plt.tight_layout()

#%%
import math
import pandas as pd
n_intervals = 500
intervals = list(np.linspace(0, Dt, n_intervals))
#dataframe where we will put the noed occupation at each time interval
occupation_df = pd.DataFrame({node:np.full_like(intervals, 0) for node 
                              in erd_reny.nodes})

#we take every walk from n_walk realizations
for walk in paths:
    # initialize a dataframe with the node where the walker jumps to and the 
    #time at which they happen of each walker
    occupation_evo_df = pd.DataFrame({node:np.full_like(intervals, None) for node 
                                  in erd_reny.nodes})
    #for every jump we find the interval correponding to the jump time
    for i, step in enumerate(walk):
        index =  0
        nodes = step[0]
        times = step[1]
        for node, t in zip(nodes, times):
            # node, t = item[0], item[1]
            if i == 0:
                j = 0
                occupation_evo_df.iloc[j, index] = node
            else:
                intervals.append(t)
                intervals = sorted(intervals)
                j = intervals.index(t)
                intervals.remove(t)
                occupation_evo_df.iloc[j, index] = node
            index += 1
    #filling the df with the current node at all times
    for column in occupation_evo_df.columns:
        last_node = occupation_evo_df.loc[0, column]
        for n in range(0,n_intervals):
            if math.isnan(occupation_evo_df.loc[n, column]):
                occupation_evo_df.loc[n, column] = last_node
            else:
                last_node = occupation_evo_df.loc[n, column]
    #now we count how many walkers are at each node at each time step
    for i_row in range(0, len(occupation_df)):
        row = occupation_evo_df.iloc[i_row]
        occupations = dict(Counter(row))
        for key, value in occupations.items():
            occupation_df.loc[i_row, key] += value
            # new_value = cell + value
            # occupation_df.loc[i_row, key] = new_value

occupation_df = occupation_df.divide(n_walk*len(list(erd_reny.nodes)))

#%%
plt.figure()
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 500)
sol = pd.DataFrame(solution_ER/len(erd_reny.nodes))
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
for i in range(len(erd_reny.nodes)):
    color = np.random.uniform(0,1, size = 3)
    plt.plot(t,sol[i], c = color, ls = '-')
    plt.plot(t2, occupation_df[i],c = color, ls = '--')
plt.show()
plt.tight_layout()

print(np.sum(solution_ER[-1,:]/len(erd_reny.nodes)))
print(np.array(list(Counter(list(dict(erd_reny.degree).values())).keys()))/
      (2*len(erd_reny.edges)))