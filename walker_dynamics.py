#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:56:41 2023

@author: robertbenassai
"""

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
import seaborn as sns
import csv
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
    return(G)

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
    return(G, path_tot)  

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
    # ONLY FOR PBC_LATTICE ----------------------------------------------------
    # lattice_gap = 1
    # for i in range(N):
    #     delta_t[(i, i+N*N-N)] = lattice_gap
    # #    print((i, i+N*N-N) in list(PBC_lattice.edges))
    #     if i==0:
    #         delta_t[(i,i+N-1)] = lattice_gap 
    # #        print((i,i+N-1) in list(PBC_lattice.edges))
    #     else:
    #         delta_t[(i*N, i*N+N-1)] = lattice_gap
    # #        print((i*N, i*N+N-1) in list(PBC_lattice.edges))
    # nx.set_edge_attributes(G, delta_t, name='dt')
    #--------------------------------------------------------------------------
    
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
                # print('all walkers finished, steps needed for this realisation '
                #       + str(step)+'/'+str(max_steps))
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

# %% Adjoint Node-centric random walk
def cumulative_exp(x, rate):
    return(1-np.exp(-x*rate))

def adjoint_node_centric(G_adj, G, Dt, v, n_walk, scales, ids_edge):
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
    # upper bound for steps
    min_edge_time = max(scales.values())
    max_steps = int(Dt*min_edge_time)
    overall_adj_edge_weights = dict.fromkeys(list(G_adj.edges), 0)
    print('maximum steps allowed each realisation', max_steps)
    # overall_edge_weights account for all the edge passings summed for
    # n_walk iterations
    overall_edge_weights = dict.fromkeys(G.edges, 0)
    for i in range(0, n_walk):
        path = []
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {i: i for i in G_adj.nodes}
        prev_pos = {i: i for i in G_adj.nodes}
        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        adj_edge_weights = dict.fromkeys(list(G_adj.edges), 0)
        edge_weights = dict.fromkeys(G.edges, 0)
        # first let's difine the time used for each walker
        time = dict.fromkeys(G_adj.nodes, 0)
        path.append([[ids_edge[ind] for ind in list(pos_walkers.values())], 
                     [time[ind] for ind in pos_walkers.keys()]])
        active_walkers = dict.fromkeys(G_adj.nodes, True)
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
            neighbors = {walker:list(nx.all_neighbors(G_adj, pos_walkers[walker]))
                         for walker in active_ls}
            # now randomly move the random walker to another node only if time+edge_time
            # is < Dt
            for walker, possible_nodes in neighbors.items():
                # time of each possible edge for a given node
                curr_edge = ids_edge[pos_walkers[walker]]
                interval = np.random.exponential(1/scales[curr_edge])
                if time[walker] + interval <= Dt:
                    # random selection between the final nodes
                    sel = random.choice(possible_nodes)
                    
                    # updating edge flow dict in the node of the adjoint graph
                    final_edge = ids_edge[sel]
                    i = curr_edge[0]
                    j = curr_edge[1]
                    prev_edge = ids_edge[prev_pos[walker]]
                    if step > 0:
                        final_edge = ids_edge[sel]
                        i = curr_edge[0]
                        j = curr_edge[1]
                        prev_edge = ids_edge[prev_pos[walker]]
                        
                        if i in prev_edge and j in final_edge:
                            edge_weights[curr_edge] += 1
                            if step == 1:
                                if random.choice([True, False]):
                                    if prev_edge[1] == i:
                                        edge_weights[prev_edge] += 1
                                    else:
                                        edge_weights[prev_edge] -= 1
                        elif j in prev_edge and i in final_edge:
                            edge_weights[curr_edge] -= 1
                            if step == 1:
                                if random.choice([True, False]):
                                    if prev_edge[1] == j:
                                        edge_weights[prev_edge] += 1
                                    else:
                                        edge_weights[prev_edge] -= 1

                    # updating walker position and time
                    curr_ind = pos_walkers[walker]
                    prev_pos[walker] = curr_ind
                    if (pos_walkers[walker], sel) in adj_edge_weights.keys():
                        adj_edge_weights[(pos_walkers[walker], sel)] += 1
                    else:
                        adj_edge_weights[(sel, pos_walkers[walker])] -= 1
                    
                    pos_walkers[walker] = sel
                    time[walker] += interval
                else:
                    #update the edge flow of the final edge
                    i = curr_edge[0]
                    j = curr_edge[1]
                    prev_edge = ids_edge[prev_pos[walker]]
                    if random.choice([True, False]):
                    # thrs_prob = cumulative_exp(Dt-time[walker],1/scales[curr_edge])
                    # if random.uniform(0, 1) <= thrs_prob: 
                        
                        if i in prev_edge:
                            edge_weights[curr_edge] += 1
                        else:
                            edge_weights[curr_edge] -= 1

                    # desactivate the walker that has finished
                    active_walkers[walker] = False
                    # print(Counter(active_walkers.values()))
            path.append([[ids_edge[ind] for ind in list(pos_walkers.values())], 
                         [time[ind] for ind in pos_walkers.keys()]])
        path_tot.append(path)
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
        for adj_edge, flow in adj_edge_weights.items():
            overall_adj_edge_weights[adj_edge] += flow
        
    #Retrieving the edge flow
    # for id_edge, edge in ids_edge.items():
        
    #     i = edge[0]
    #     j = edge[1]
    #     #list of edges in the adjoint graph
    #     edges_w_flow = list(overall_adj_edge_weights.keys())
    #     nbrs_e = [(id_edge, neigh) if (id_edge, neigh) in edges_w_flow else 
    #               (neigh, id_edge) for neigh in G_adj.neighbors(id_edge)]
    #     adjac_e = []
        
    #     for adj_edge in nbrs_e:
    #         edge_1 = ids_edge[adj_edge[0]]
    #         edge_2 = ids_edge[adj_edge[1]]
    #     #     if j not in edge_1 or i not in edge_2:
    #     #         adjac_e.append(overall_adj_edge_weights[adj_edge])
    #     #     else:
    #     #         adjac_e.append(-overall_adj_edge_weights[adj_edge])
    #     # overall_edge_weights[edge] += 0.5*np.sum(adjac_e)
        
        
    #     in_flow = []
    #     out_flow = []
        
    #     for adj_edge in nbrs_e:
    #         edge_1 = ids_edge[adj_edge[0]]
    #         edge_2 = ids_edge[adj_edge[1]]
    #         if j not in edge_1:
    #             in_flow.append(overall_adj_edge_weights[adj_edge])
    #         elif j not in edge_2:
    #             in_flow.append(-overall_adj_edge_weights[adj_edge])

    #         elif i not in edge_2:
    #             out_flow.append(overall_adj_edge_weights[adj_edge])
    #         elif i not in edge_1:
    #             out_flow.append(-overall_adj_edge_weights[adj_edge])

    #     total_in_fl = np.sum(in_flow)
    #     total_out_fl = np.sum(out_flow)
    #     if total_in_fl >= 0 and total_out_fl >= 0:
    #         overall_edge_weights[edge] = total_out_fl #+ total_in_fl)/2
    #     elif total_in_fl <= 0 and total_out_fl <= 0:
    #         overall_edge_weights[edge] = total_in_fl #+ total_out_fl)/2
    #     elif total_in_fl >= 0 and total_out_fl <= 0:
    #         overall_edge_weights[edge] = 0#(total_in_fl-total_out_fl)/2
            
    #     elif total_in_fl <= 0 and total_out_fl >= 0:
    #         overall_edge_weights[edge] = 0#(total_out_fl - total_in_fl)/2
        
    
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G, path_tot)
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
#%% CONTINUUM NEWMANN WITHOUT ABSORBINGS

'''first we build the transition rates matrix'''

# building the transition rate matrix
#dict of matrix ind to node index
def build_trans_rates_matrix(G, pos, v, new_edges, short_1, short_2):
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
    # print(Counter(dists.values()))
    # ONLY FOR PBC_LATTICE ----------------------------------------------------
    lattice_gap = 1
    N = 8
    for i in range(N):
        dists[(i, i+N*N-N)] = lattice_gap
        dists[(i+N*N-N, i)] = lattice_gap
        if i==0:
            dists[(i,i+N-1)] = lattice_gap 
            dists[(i+N-1, i)] = lattice_gap 
        else:
            dists[(i*N, i*N+N-1)] = lattice_gap
            dists[(i*N+N-1, i*N)] = lattice_gap
    # -------------------------------------------------------------------------
    for new_e in new_edges:
        if new_e in dists.keys():
            if dists[new_e] == 0.25:
                dists[new_e] /= short_1
            else:
                dists[new_e] /= short_2
        else:
            orig_dist = dists[(new_e[1], new_e[0])]
            if orig_dist == 0.25:
                dists[(new_e[1], new_e[0])] /= short_1
            else:
                dists[(new_e[1], new_e[0])] /= short_2

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
        for final in nx.all_neighbors(G_und, node):
            i, j = inode_to_iarr[node], inode_to_iarr[final]
            # if (node,final) in intervals.keys():
            #     rate = intervals[(node, final)]
            # else:
            #     rate = intervals[(final, node)]
            trans_rates[j][i] = 1/(k_deg*mean_dist)
    # the diagonal is -sum of the off diag elements of the col
    i,j = np.indices(trans_rates.shape)
    diagonal = np.sum(trans_rates, axis=0)
    trans_rates[i==j] = -diagonal
    return(trans_rates)#, dists)

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
    p0 = np.ones(len(G.nodes))
    
    t = np.linspace(start=0, stop=Dt,num=20000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    _, tri_dict = find_triangles(G)
    prob_product = []
    for triangle in tri_dict.keys():
        i = triangle[0]
        j = triangle[1]
        k = triangle[2]
        prob_product.append(sol[:,i]*sol[:,j]*sol[:,k]*(trans_rates[i][j]*\
        trans_rates[j][k]*trans_rates[k][i] - trans_rates[j][i]*\
        trans_rates[k][j]*trans_rates[i][k]))
    
    prob_prod_arr = np.transpose(np.array(prob_product))
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    curls = np.trapz(prob_prod_arr, t, axis=0)
    for j, p_j in enumerate(flows):
        
        for i in range(len(G.nodes)):
            net_flow[i,j] += flows[i]*trans_rates[j][i] - p_j*trans_rates[i][j]

    net_flow *= n_walk
    cont_new = {edge: net_flow[inode_to_iarr[edge[0]]][inode_to_iarr[edge[1]]] 
                for edge in G.edges}
    nx.set_edge_attributes(G, cont_new, name = 'edge_visits')
    return(G, sol, curls)

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
# %% PLOT HODGE FUNCTION

def plot_hodge(walk_graph, grad_comp, sol_comp, har_comp, pot, div, pos):
    '''
    Plot of the hodge decomposition of a Graph

    Parameters
    ----------
    walk_graph : nx.DiGraph
        Original graph with 'edge_visits' as edge_attr.
    grad_comp : dict
        Gradient component {edge: g_ij}.
    sol_comp : dict
        Solenoidal component {edge: s_ij}..
    har_comp : dict
        Harmonic component {edge: h_ij}..
    pot : dict
        Node potentials {node: pot_i}.
    div : dict
        Divergence {node: div}.
    pos : dict
        Position of nodes {node: pos_i}.

    Returns
    -------

    '''

    edge_graph = nx.get_edge_attributes(walk_graph, 'edge_visits')
    
    g_field = np.array([walk_graph[edge[0]][edge[1]]['edge_visits'] for edge
                                     in walk_graph.edges])
    percentile = np.percentile(g_field, 95)
    percentile_g = np.percentile(list(grad_comp.values()), 95)
    percentile_s = np.percentile(list(sol_comp.values()), 95)
    percentile_h = np.percentile(list(har_comp.values()), 95)
    
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
    plt.subplots(1, 1, figsize=(9, 6))
    # plt.subplot(121)
    # plt.title('Original Graph 100%')

    # # color_p = np.abs(np.array(list(edge_graph.values())))
    # color_p = np.array(list(edge_graph.values()))
    # # colors = np.linspace(0, percentile)
    # # cmap = plt.cm.Oranges
    # colors = np.linspace(-percentile, percentile)
    # cmap = plt.cm.seismic
    # vmin = min(colors)
    # vmax = max(colors)

    # # color_div = list(div.values())
    # # colors_div = range(int(min(color_div)), int(max(color_div)))
    # # cmap_div = plt.cm.RdBu_r
    # # vmin_div = min(colors_div)
    # # vmax_div = round(max(color_div))

    # nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4, 
    #                        node_color='#D3D3D3')
    #                         # node_color=color_div, cmap=cmap_div, vmin=vmin_div,
    #                         # vmax=vmax_div)
    # nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_p,
    #                        edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
    #                        arrowsize = 5, node_size = 4)

    # sm = plt.cm.ScalarMappable(cmap=cmap, 
    #                            norm=plt.Normalize(vmin=vmin, vmax=vmax))
    # sm._A = []
    # cbar = plt.colorbar(sm)
    # cbar.set_label(r'\omega$', fontsize = 18)

    # # sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
    # #                                                               vmax=vmax_div))
    # # sm2._A = []
    # # cbar2 = plt.colorbar(sm2, location='right')
    # # cbar2.set_label(r'Node div')

    plt.subplot(111)
    color_g = np.array(list(grad_comp.values()))
    # plotting edges with color gradient

    color_pot = list(pot.values())
    cmap_pot = plt.cm.PRGn
    vmax_pot = int(np.max(color_pot))
    vmin_pot = int(np.min(color_pot))


    colors = np.linspace(0, percentile_g)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)
    # plt.title('Gradient component ' + str(round(weight_g*100, 1))+'%', fontsize = 18)
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None,
                            node_size=30, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 5, node_size = 30)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\omega_g$', fontsize = 18 )
    cbar.ax.tick_params(labelsize=18)

    sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                                  vmax=vmax_pot))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Node potentials', fontsize = 18)
    cbar2.ax.tick_params(labelsize=18)

    # colors = np.linspace(0, percentile_s)
    # cmap = plt.cm.Oranges
    # vmin = min(colors)
    # vmax = max(colors)

    # color_s = np.abs(np.array(list(sol_comp.values())))
    # plt.subplot(223)
    # plt.title('Solenoidal Component ' + str(round(weight_s*100, 1))+'%')
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

    # colors = np.linspace(0, percentile_h)
    # cmap = plt.cm.Oranges
    # vmin = min(colors)
    # vmax = max(colors)

    # color_h = np.array(list(har_comp.values()))
    # plt.subplot(224)
    # plt.title('Harmonic Component ' + str(round(weight_h*100, 1))+'%')
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
    
def reverse_negative_edges(G):
    edge_ls = list(G.edges)
    for edge in edge_ls:
        flow = G[edge[0]][edge[1]]['edge_visits']
        if flow < 0:
            G.remove_edge(*edge)
            #reversing edge
            G.add_edge(edge[1], edge[0], edge_visits = -flow)
    return G
#%%
def plot_pot_corr(pot_sim, pot_theo, x_label, y_label):
    pot_list = [(pot_sim[i], pot_theo[i]) for i in pot_sim.keys()]
    x, y = zip(*pot_list)
        
    a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    plt.figure(figsize=(8,6))
    plt.xlabel(x_label, fontsize = 20)
    plt.ylabel(y_label, fontsize = 20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.scatter(x, y, s = 18)
    plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+r'$\pm$'
             +str(round(std_err,2))+'x +'+str(round(b, 2))+'\n'+r'$r^2 = $' 
             +str(round(r_value**2, 2)))
    plt.plot(x, np.array(x), c='r', label = 'y = x')
    plt.legend(fontsize = 18)
    plt.tight_layout()
# PLOT OF PREDICTED VS SIMULATED GRADIENT COMPONENT 
def plot_grad_corr(grad_sim, grad_theo, x_label, y_label):
    grad_list = [(grad_sim[i], grad_theo[i]) for i in grad_sim.keys()]
    x, y = zip(*grad_list)
        
    a, b, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    plt.figure(figsize=(8,6))
    plt.xlabel(x_label,  fontsize = 20)
    plt.ylabel(y_label,  fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.scatter(x, y, s = 18)
    plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = '+str(round(a,2))+r'$\pm$'
             +str(round(std_err,2))+'x +'+str(round(b, 2))+'\n'+r'$r^2 = $' +
             str(round(r_value**2, 2)))
    plt.plot(x, np.array(x), c='r', label = 'y = x')
    plt.legend(fontsize = 18)
    plt.tight_layout()

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
#%% ERDOS_RENYI

'''------------------Testing with an Erdös-Rényi graph----------------------'''

''' Building the graph'''

erd_reny = nx.erdos_renyi_graph(50, 0.1, seed = 1000, directed=False)
erd_reny = erd_reny.to_directed()
out_edges = [edge for edge in erd_reny.edges if edge[1]
    < edge[0]]  # removing all outward edges
erd_reny.remove_edges_from(out_edges)
pos_ER = nx.spring_layout(erd_reny, seed = 1050)
nx.draw_networkx(erd_reny, with_labels = False, node_size = 20, pos = pos_ER)

#%% Random Geometric
# Use seed when creating the graph for reproducibility
G = nx.random_geometric_graph(50, 0.2, seed=1000)
# position is stored as node attribute data for random_geometric_graph
pos_ER = nx.get_node_attributes(G, "pos")
# transforming the graph into a digraph
erd_reny = nx.DiGraph(G)
out_edges = [edge for edge in erd_reny.edges if edge[1]
    < edge[0]]  # removing all outward edges
erd_reny.remove_edges_from(out_edges)

nx.draw_networkx(erd_reny, with_labels = False, node_size = 50, pos = pos_ER)
#%%
# Set common styling parameters
common_params = {
    "node_size": 100,
    "node_color": "#2B81F2",
    "linewidths": 2,
    "edgecolors": "black",
}

# Draw nodes and edges
nx.draw_networkx_nodes(erd_reny, pos=pos_ER, **common_params)
nx.draw_networkx_edges(erd_reny, pos=pos_ER, width=2, arrowstyle="-|>", 
                       connectionstyle="arc3,rad=0.0", node_size=100)

#%% DISCRETE WALK

#SIM DISCR
steps = 100
n_walk = 200

discr_walk_ER, occupations_disc = digraph_walkers(erd_reny.copy(), steps, n_walk)

grad_discr_sim, sol_discr_sim, har_discr_sim, pot_discr_sim, div_discr_sim = \
    hodge_decomposition(discr_walk_ER, 'edge_visits')

#THEO DISCR
trans_matrix = build_trans_matrix(erd_reny.copy())
discr_ER_theo, theo_occupations = discrete_rw_edge_flow(erd_reny.copy(), 
                                                        trans_matrix, steps, n_walk)

grad_discr_th, sol_discr_th, har_discr_th, pot_discr_th, div_discr_th = \
    hodge_decomposition(discr_ER_theo, 'edge_visits')
    
#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS CTRW

plot_pot_corr(pot_discr_sim, pot_discr_th, 'Simulated Potential', 
              'Analytical Potential')
plot_grad_corr(grad_discr_sim, grad_discr_th, 'Simulated gradient component', 
               'Analytical gradient component')
#%% NODE-CENTRIC WALK
Dt = 20
n_walk = 120
v = 1
walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
grad_ER, sol_ER, har_ER, pot_ER, div_ER = hodge_decomposition(walked_ER, 'edge_visits')

# THEORETICAL 
#%%
trans_rate_ER = build_trans_rates_matrix(erd_reny.copy(), pos_ER, v)
ER_th, solution_ER, curls_ER = solve_continuous_rw_flow(erd_reny.copy(), trans_rate_ER, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = hodge_decomposition(ER_th, 'edge_visits')
print(curls_ER)
#%%
plot_hodge(ER_th, g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new,
           pos_ER)
#%%
plot_hodge(walked_ER, grad_ER, sol_ER, har_ER, pot_ER, div_ER,
           pos_ER)
#%%
plot_pot_corr(pot_ER, pot_cont_new, 'Simulated Potential', 
              'Analytical Potential')
plot_grad_corr(grad_ER, g_cont_new, 'Simulated gradient component', 
               'Analytical gradient component')
#%%
import math
n_intervals = 500
intervals = list(np.linspace(0, Dt, n_intervals))
#dataframe where we will put the node occupation at each time interval
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

#%% DISCR PROB EVO
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='black', ls = '-.'),
                Line2D([0], [0], color = 'black', ls = '-')]

time = np.linspace(0, 101, 101)
color_p = sns.color_palette(palette = 'colorblind', n_colors = len(erd_reny.nodes))
plt.figure(figsize = (8,6))
theo_evo = np.array([list(element.values()) for element in theo_occupations])
sim_evo = np.array([list(element.values()) for element in occupations_disc])

for i in range(len(erd_reny.nodes)):
    plt.plot(time, sim_evo[:,i], ls = '-.', c = color_p[i])
    plt.plot(time, theo_evo[:,i]*n_walk, ls = '-', c = color_p[i])
plt.xlabel('time (steps)', fontsize = 18)
plt.ylabel('Occupation probabilities', fontsize = 18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

plt.legend(custom_lines, ['simulation', 'theoretical'], loc='upper center', 
           ncol = 2, bbox_to_anchor=(0.5, 1), fontsize = 15)
plt.tight_layout()
#%%
#NODE CENTRIC PROBABILITY EVO
plt.figure(figsize = (8,6))
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 500)
sol = pd.DataFrame(solution_ER/len(erd_reny.nodes))
plt.xlabel('time', fontsize = 18)
plt.ylabel('Occupation probabilities', fontsize = 18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

for i in range(len(erd_reny.nodes)):
    color = color_p[i]
    plt.plot(t,sol[i], c = color, ls = '-')
    plt.plot(t2, occupation_df[i],c = color, ls = '--')
plt.legend(custom_lines, ['simulation', 'theoretical'], loc='upper center', 
           ncol = 2, bbox_to_anchor=(0.5, 1), fontsize = 15)
plt.tight_layout()

#%% EVOLUTION OF SOLENOIDAL HARMONIC AND GRADIENT COMPONENTS WITH SIM TIME
import csv
n_walk = 20
v = 1
times = np.linspace(0.1, 201.1, 300)

with open('/Users/robertbenassai/Documents/UOC/PBC_01_ratio_evolution.csv', 
          'w') as file1:
    writer = csv.writer(file1)
    with open('/Users/robertbenassai/Documents/UOC/PBC_01_abs_avg_flow.csv', 
              'w') as file2:
        writer2 = csv.writer(file2)
        for Dt in times:
            walked_ER, paths = node_centric(PBC_lattice.copy(), Dt, v, pos_PBC, n_walk)
            grad_ER, sol_ER, har_ER, pot_ER, div_ER = hodge_decomposition(walked_ER, 'edge_visits')
            
            edge_graph = nx.get_edge_attributes(walked_ER, 'edge_visits')
            w = np.array(list(edge_graph.values()))
            wg = np.array(list(grad_ER.values()))
            ws = np.array(list(sol_ER.values()))
            wh = np.array(list(har_ER.values()))
            weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
            weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
            weight_h = np.sum(np.square(wh))/np.sum(np.square(w))
            
            mean_g = np.mean(np.array(list(grad_ER.values()))**2)
            mean_s = np.mean(np.array(list(sol_ER.values()))**2)
            mean_h = np.mean(np.array(list(har_ER.values()))**2)
            mean_cycl = np.mean((np.array(list(har_ER.values()))+\
                                 np.array(list(sol_ER.values())))**2)
            writer.writerow([Dt, weight_g, weight_s, weight_h])
            writer2.writerow([Dt, mean_g, mean_s, mean_h, mean_cycl])
        
#%%read the results from a file

weight_ls_g = []
weight_ls_s = []
weight_ls_h = []

mean_abs_g = []
mean_abs_s = []
mean_abs_h = []
mean_abs_cycl = []
times = []
times2 = []

with open('/Users/robertbenassai/Documents/UOC/PBC_01_ratio_evolution.csv', 'r') \
    as file:
        with open('/Users/robertbenassai/Documents/UOC/PBC_01_abs_avg_flow.csv', 
              'r') as file2:
            reader = csv.reader(file, delimiter=',')
            reader2 = csv.reader(file2, delimiter = ',')
            for row in reader:
                times.append(float(row[0]))
                weight_ls_g.append(float(row[1]))
                weight_ls_s.append(float(row[2]))
                weight_ls_h.append(float(row[3]))
                
            for row2 in reader2:
                times2.append(float(row2[0]))
                mean_abs_g.append(float(row2[1]))
                mean_abs_s.append(float(row2[2]))
                mean_abs_h.append(float(row2[3]))
                mean_abs_cycl.append(float(row2[4]))
#%%
from scipy.optimize import curve_fit
def func(x, a, b):
    y = b*x**a
    return y

def func2(x, a):
    y = a*np.ones_like(x)
    return y
params, pcov = curve_fit(func, times[100:], weight_ls_g[100:])
#%% PLOT OF STRENGTH EVOLUTION
st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(PBC_lattice)

plt.figure()
plt.loglog(times, weight_ls_g, color = 'b' , linestyle = '-', marker = '.', 
         label = 'gradient strength ratio')
plt.loglog(times, np.array(weight_ls_s) + np.array(weight_ls_h), color = 'r', 
         linestyle = '-', marker = '.', label = 'cyclic strength ratio')
# plt.hlines(st_ratio_g, 0, 400, color = 'b' , linestyle = '--', label = 'gradient structural ratio')
# plt.hlines((st_ratio_h+st_ratio_s), 0, 400, color = 'r', linestyle = '--', 
#             label = 'cyclic structural ratio')
plt.loglog(times[40:], params[1]*times[40:]**params[0], 'b--', label = r'$y = t^{'+
            str(round(params[0], 2))+'}$')
plt.loglog(times[40:], 1-params[1]*times[40:]**params[0], 'r--', label = r'$y = 1 - t^{'+
            str(round(params[0],2))+'}$')
plt.xlabel('Simulation time')
plt.ylabel('Strength Ratio')
# plt.tight_layout()
plt.legend(loc='upper center', ncol = 2, bbox_to_anchor=(0.5, 1.18))
perr = np.sqrt(np.diag(pcov))
print(perr)

#%% PLOT MEAN SQUARED FLOW
params_cycl, pcov_cycl = curve_fit(func, times2, np.array(mean_abs_cycl))
params_tot, pcov_tot = curve_fit(func, times2, np.array(mean_abs_h)+np.array(mean_abs_s))
params_g, pcov_g = curve_fit(func2, times2[100:], mean_abs_g[100:])

perr_cycl = np.sqrt(np.diag(pcov_cycl))
perr_tot = np.sqrt(np.diag(pcov_tot))
print(perr_cycl)

plt.figure()
plt.loglog(times2, mean_abs_g, color = 'b' , linestyle = '-', marker = '.', 
         label = 'mean squared gradient component')
plt.loglog(times2, np.array(mean_abs_cycl), color = 'r', 
         linestyle = '-', marker = '.', label = 'mean squared solenoidal component')
# plt.plot(times2, np.array(mean_abs_s)+np.array(mean_abs_s), 
#           color = 'g',linestyle = '-', marker = '.', label = 'harmonic mean squared flow')

plt.loglog(times2, func2(np.array(times2), params_g[0]), 'b--', 
         label = r'$\left < \omega^2 \right > ='+str(round(params_g[0], 2))+'$')
plt.loglog(times2, func(np.array(times2), params_cycl[0], params_cycl[1]), 'r--', 
           label = r'$\left < \omega^2 \right > = \:('+str(round(params_cycl[1],2))+
           '\pm'+str(round(perr_cycl[1],2))+') t^{'+str(round(params_cycl[0],2))+
           '\pm'+str(round(perr_cycl[0],2))+'}$')

# plt.plot(times2, func(np.array(times2), params_tot[0], params_tot[1]), 'g--', 
#            label = r'$\left < \omega^2 \right > = \:('+str(round(params_tot[1],2))+
#            '\pm'+str(round(perr_tot[1],2))+') t^{'+str(round(params_tot[0],2))+
#            '\pm'+str(round(perr_tot[0],2))+'}$')
# func(np.array(times2), params_cycl[0])
plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\left < \omega^2 \right >$')
plt.legend()
plt.tight_layout()

#%%EVOLUTION OF GRADIENT AND CYCLIC FLOW

n_walk = 20
v = 1
times = np.linspace(0.1, 201.1, 300)

with open('/Users/robertbenassai/Documents/UOC/ER_50_01_grad_evolution.csv', 
          'w') as file1:
    writer = csv.writer(file1)
    with open('/Users/robertbenassai/Documents/UOC/ER_50_01_sol_evolution.csv', 
              'w') as file2:
        writer2 = csv.writer(file2)
        with open('/Users/robertbenassai/Documents/UOC/ER_50_01_har_evolution.csv', 
                  'w') as file3:
            writer3 = csv.writer(file3)
            for Dt in times:
                walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
                grad_ER, sol_ER, har_ER, pot_ER, div_ER = hodge_decomposition(walked_ER, 'edge_visits')
                
                grad_list = [grad_ER[key] for key in sorted(grad_ER.keys())]
                sol_list = [sol_ER[key] for key in sorted(sol_ER.keys())]
                har_list = [har_ER[key] for key in sorted(har_ER.keys())]
                
                writer.writerow([Dt] + grad_list)
                writer2.writerow([Dt] + sol_list)
                writer3.writerow([Dt] + har_list)
#%%READING AND PLOTTING

grad_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/ER_50_01_grad_evolution.csv', 
                       sep = ',', header=None)
sol_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/ER_50_01_sol_evolution.csv', 
                       sep = ',', header=None)
har_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/ER_50_01_har_evolution.csv', 
                      sep = ',', header=None)

plt.figure()
plt.xlabel('total simulaiton time')
plt.ylabel(r'$\omega^2$')
for i in range(1, len(grad_evo.T)):
    plt.plot(times,(har_evo.iloc[:, i] + sol_evo.iloc[:, i])**2)
#%% ADJOINT MATRIX

def build_adjoint_graph(G:nx.DiGraph()):
    '''
    Builds the adjacent graph of G.

    Parameters
    ----------
    G : nx.DiGraph()
        Input graph.

    Returns
    -------
    adjoint_graph : nx.Graph
        Adjoint graph of G.

    '''
    #redefine edges as nodes with a dictionary: edge_id: edge
    edge_ids = {edge: i for i, edge in enumerate(list(G.edges))}
    ids_edge = {i: edge for i, edge in enumerate(list(G.edges))}
    #now we need a dictionary that has node: [(id_edge1, id_edge_2), ...]
    new_edges = {}
    adjoint_graph = nx.Graph()
    adjoint_graph.add_nodes_from(list(edge_ids.values()))
    
    for node in G.nodes:
        adj_edge_id_ls = []
        if len(list(nx.all_neighbors(G, node))) > 1:
            adj_edges = list(G.in_edges(node))+list(G.out_edges(node))
        for i in range(len(adj_edges)):
            edge_1 = adj_edges.pop(0)
            for edge_2 in adj_edges:
                adj_edge_id_ls.append((edge_ids[edge_1], edge_ids[edge_2]))
        adjoint_graph.add_edges_from(adj_edge_id_ls)

        new_edges[node] = adj_edge_id_ls
    return (adjoint_graph, ids_edge)

#%% 

# building the transition rate matrix
#dict of matrix ind to node index
def build_trans_adjoint_rates_matrix(G_adj:nx.Graph(), rates:dict, ids_edge:dict):
    '''
    Builds the transition rates matrix or generator matrix for an unweighted 
    random walk in a topological graph for the adjoint graph

    Parameters
    ----------
    G : nx.Graph
        Graph from which the transition rates will be computed.
    rates : Dict
        dictionary of node:rate. (rate = v/d)
    Returns
    -------
    trans_rates : npumpy ndarray.
                  transition rates matrix

    '''    
    #mapping of each node to its index in the transition rate array
    inode_to_iarr = {node:i for i, node in enumerate(G_adj.nodes)}
    #building the array
    N_e = len(G_adj.nodes)
    trans_rates = np.zeros((N_e, N_e))
    for node in G_adj.nodes:
        #Degree of the departure node
        k_deg = len(list(G_adj.neighbors(node)))
                             
        #for each neighbour we calculate the transition rate probability of going 
        for final in G_adj.neighbors(node):
            i, j = inode_to_iarr[node], inode_to_iarr[final]
            trans_rates[j][i] = rates[ids_edge[node]]/(k_deg)
    # the diagonal is -sum of the off diag elements of the col
    i,j = np.indices(trans_rates.shape)
    diagonal = np.sum(trans_rates, axis=0)
    trans_rates[i==j] = -diagonal
    return(trans_rates)


#%%
def solve_adjoint_ctrw_flow(G: nx.DiGraph(), G_adj:nx.Graph, trans_rates:np.array, 
                            ids_edge: dict, rates:dict, Dt:float, n_walk:int):
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
    N = len(list(G_adj.nodes))
    #mapping of each node to its index in the transition rate array
    inode_to_iarr = {node:i for i, node in enumerate(G_adj.nodes)}
    iarr_to_inode = {i:node for i, node in enumerate(G_adj.nodes)}
    #array of flows through the edges
    net_flow = np.zeros((N,N))
    #initial probabilities
    p0 = np.ones(len(G_adj.nodes))
    t = np.linspace(start=0, stop=Dt,num=20000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    
    #edge flow for adjacent graph
    for j, p_j in enumerate(flows):
        
        for i in range(len(G_adj.nodes)):
            net_flow[i,j] += flows[i]*trans_rates[j][i] - p_j*trans_rates[i][j]

    net_flow *= n_walk
    cont_new_adj = {adj_edge: net_flow[inode_to_iarr[adj_edge[0]]][inode_to_iarr[adj_edge[1]]] 
                for adj_edge in G_adj.edges}
    final_flow = {}
    for id_edge, edge in ids_edge.items():
        
        i = edge[0]
        j = edge[1]
        #list of edges in the adjoint graph
        edges_w_flow = list(cont_new_adj.keys())
        nbrs_e = [(id_edge, neigh) if (id_edge, neigh) in edges_w_flow else 
                  (neigh, id_edge) for neigh in G_adj.neighbors(id_edge)]
        
        #Computing the final edge flow. The edge flow is defined as the ougoing
        #flow of the edge if both the in and out flow go in the same direction.
        #If in and out flows go in opposite directions the walkers die at the edge
        #and the nodes of that link are sinks or sources.
        
        in_flow = []
        out_flow = []
        adjac_e = []
        for adj_edge in nbrs_e:
            edge_1 = ids_edge[adj_edge[0]]
            edge_2 = ids_edge[adj_edge[1]]
        #     if j not in edge_1 or i not in edge_2:
        #         adjac_e.append(cont_new_adj[adj_edge])
        #     else:
        #         adjac_e.append(-cont_new_adj[adj_edge])
        # final_flow[edge] = 0.5*np.sum(adjac_e)
        
            if j not in edge_1:
                in_flow.append(cont_new_adj[adj_edge])
            elif j not in edge_2:
                in_flow.append(-cont_new_adj[adj_edge])

            elif i not in edge_2:
                out_flow.append(cont_new_adj[adj_edge])
            elif i not in edge_1:
                out_flow.append(-cont_new_adj[adj_edge])

        total_in_fl = np.sum(in_flow)
        total_out_fl = np.sum(out_flow)
        if total_in_fl >= 0 and total_out_fl >= 0:
            final_flow[edge] = (total_out_fl + total_in_fl)/2
            # final_flow[edge] = total_out_fl
        elif total_in_fl <= 0 and total_out_fl <= 0:
            final_flow[edge] = (total_in_fl + total_out_fl)/2
            # final_flow[edge] = total_in_fl
        elif total_in_fl >= 0 and total_out_fl <= 0:
            final_flow[edge] = (total_in_fl+total_out_fl)/2
            # final_flow[edge] = 0
            
        elif total_in_fl <= 0 and total_out_fl >= 0:
            final_flow[edge] = (total_out_fl + total_in_fl)/2
            # final_flow[edge] = 0
            
    nx.set_edge_attributes(G, final_flow, name = 'edge_visits')
    return(G, sol)

#%%
erd_reny = nx.erdos_renyi_graph(5, 0.5, seed = 1000, directed=False)
erd_reny = erd_reny.to_directed()
out_edges = [edge for edge in erd_reny.edges if edge[1]
    < edge[0]]  # removing all outward edges
erd_reny.remove_edges_from(out_edges)
pos_ER = nx.spring_layout(erd_reny, seed = 1050)
nx.draw_networkx(erd_reny, with_labels = True, node_size = 500, pos = pos_ER)
#%%

adj_ER, ids_edge_ER = build_adjoint_graph(erd_reny)
pos_adj = nx.spring_layout(adj_ER, seed = 1050)
#%% PLOT OF ORIGINAL VS ADJOINT
plt.subplots(1,2, figsize = (8, 4))
plt.subplot(121)
plt.title('original')
nx.draw_networkx(erd_reny, pos = pos_ER)
plt.subplot(122)
plt.title('adjoint')
nx.draw_networkx(adj_ER)
plt.tight_layout()
#%%
node_size = 500
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
plt.subplots_adjust(hspace=0.4)

# Set common styling parameters
common_params = {
    "node_size": node_size,
    "node_color": "#AED0EE",
    "linewidths": 2,
    "edgecolors": "black",
}

# Set styling parameters for edge labels
edge_label_params = {
    "font_size": 16,
    "bbox": {"facecolor": "white", "edgecolor": "none", "pad": 0.5},
}

# Iterate over subplots

ax[0].set_title('Original', fontsize=25)

# Draw nodes and edges
nx.draw_networkx_nodes(erd_reny, pos=pos_ER, ax=ax[0], **common_params)
nx.draw_networkx_edges(erd_reny, pos=pos_ER, width=2, arrowstyle="-|>", 
                       connectionstyle="arc3,rad=0.0", ax=ax[0], node_size=node_size)

# Draw labels and edge labels
nx.draw_networkx_labels(erd_reny, pos_ER, font_size=18, font_color="black", 
                        ax=ax[0], verticalalignment='center_baseline')

ax[0].tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax[0].set_facecolor("#F7F7F7")

ax[1].set_title('Adjoint', fontsize=25)

# Draw nodes and edges
nx.draw_networkx_nodes(adj_ER, pos=pos_adj, ax=ax[1], **common_params)
nx.draw_networkx_edges(adj_ER, pos=pos_adj, width=2, arrowstyle="-|>", 
                       connectionstyle="arc3,rad=0.0", ax=ax[1], node_size=node_size)

# Draw labels and edge labels
nx.draw_networkx_labels(adj_ER, pos = pos_adj, font_size=18, font_color="black", 
                        ax=ax[1], verticalalignment='center_baseline')


# Set axis properties
ax[1].tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax[1].set_facecolor("#F7F7F7")

plt.tight_layout()
plt.show()

#%% ADJOINT THEORETICAL
Dt = 20
n_walk = 120
v = 1

rates_ER = {edge: v/np.linalg.norm(np.array(pos_ER[edge[0]])-np.array(pos_ER[edge[1]]))
        for edge in erd_reny.edges}


adj_trans_rates = build_trans_adjoint_rates_matrix(adj_ER.copy(), rates_ER, ids_edge_ER)

ER_theo_adj, sol_ER_adj = solve_adjoint_ctrw_flow(erd_reny.copy(), adj_ER.copy(), 
                                                  adj_trans_rates, ids_edge_ER,
                                                  rates_ER, Dt, n_walk)

total_fl_theo = {edge:ER_theo_adj[edge[0]][edge[1]]['edge_visits'] for edge in 
                ER_theo_adj.edges}

grad_adj, sol_adj, har_adj, pot_adj, div_adj = hodge_decomposition(ER_theo_adj,
                                                                   'edge_visits')
plot_hodge(ER_theo_adj, grad_adj, sol_adj, har_adj, pot_adj, div_adj, pos_ER)

#%% ADJOINT SIMULATION

walked_ER_edge, path_edges_ER = adjoint_node_centric(adj_ER.copy(), erd_reny.copy()
                                                     , Dt, v, n_walk, rates_ER, 
                                                     ids_edge_ER)

total_fl_sim = {edge:walked_ER_edge[edge[0]][edge[1]]['edge_visits'] for edge in 
                walked_ER_edge.edges}

grad_sim_adj, sol_sim_adj, har_sim_adj, pot_sim_adj, div_sim_adj = \
    hodge_decomposition(walked_ER_edge,'edge_visits')

plot_hodge(walked_ER_edge, grad_sim_adj, sol_sim_adj, har_sim_adj, pot_sim_adj, 
           div_sim_adj, pos_ER)
#%%

n_intervals = 500
intervals = list(np.linspace(0, Dt, n_intervals))
#dataframe where we will put the noed occupation at each time interval
occupation_df = pd.DataFrame({str(edge):np.full_like(intervals, 0) for edge 
                              in erd_reny.edges})

#we take every walk from n_walk realizations
for walk in path_edges_ER:
    # initialize a dataframe with the node where the walker jumps to and the 
    #time at which they happen of each walker
    occupation_evo_df = pd.DataFrame({str(edge):np.full_like(intervals, None) for edge 
                                  in erd_reny.edges})
    #for every jump we find the interval correponding to the jump time
    for i, step in enumerate(walk):
        index =  0
        edges = step[0]
        times = step[1]
        for edge, t in zip(edges, times):
            # node, t = item[0], item[1]
            if i == 0:
                j = 0
                occupation_evo_df.iloc[j, index] = str(edge)

            else:
                intervals.append(t)
                intervals = sorted(intervals)
                j = intervals.index(t)
                intervals.remove(t)
                occupation_evo_df.iloc[j, index] = str(edge)
            index += 1
    #filling the df with the current node at all times
    for column in occupation_evo_df.columns:
        last_edge = occupation_evo_df.loc[0, column]
        for n in range(0,n_intervals):
            # print(type(occupation_evo_df.loc[n, column]))
            if type(occupation_evo_df.loc[n, column]) != str:
                occupation_evo_df.loc[n, column] = last_edge
            else:
                last_edge = occupation_evo_df.loc[n, column]
    #now we count how many walkers are at each node at each time step
    for i_row in range(0, len(occupation_df)):
        row = occupation_evo_df.iloc[i_row]
        occupations = dict(Counter(row))
        for key, value in occupations.items():
            occupation_df.loc[i_row, key] += value
            # new_value = cell + value
            # occupation_df.loc[i_row, key] = new_value

occupation_df = occupation_df.divide(n_walk*len(list(erd_reny.edges)))

#%% PROBABILITY PLOT
custom_lines = [Line2D([0], [0], color='black', ls = '-.'),
                Line2D([0], [0], color = 'black', ls = '-')]
plt.figure()
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 500)
sol = pd.DataFrame(sol_ER_adj/len(erd_reny.edges))
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
color_p = sns.color_palette(palette = 'colorblind', n_colors = len(erd_reny.edges))
for i, edge in enumerate(erd_reny.edges):
    if i< 10:
        plt.plot(t,sol[i], c = color_p[i], ls = '-')
        plt.plot(t2, occupation_df[str(edge)],c = color_p[i], ls = '--')
    
plt.legend(custom_lines, ['simulation', 'theoretical'], loc='upper center', 
           ncol = 2, bbox_to_anchor=(0.5, 1))
plt.tight_layout()

print(np.sum(sol_ER_adj[-1,:]/len(erd_reny.edges)))
print(np.array(list(Counter(list(dict(adj_ER.degree).values())).keys()))/
      (2*len(adj_ER.edges)))

#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS ADJOINT CTRW

plot_pot_corr(pot_sim_adj, pot_adj, 'Simulated Potential', 'Analytical Potential')
plot_grad_corr(grad_sim_adj, grad_adj, 'Simulated Gradient component',
               'Analytical Gradient Component' )

#%% CONSTANT VELOCITY WALK
Dt = 20
n_walk = 120
v = 1
walked_v_ER = node_walkers(erd_reny.copy(), Dt, v, pos_ER, n_walk)
grad_cnst_v, sol_cnst_v, har_cnst_v, pot_cnst_v, div_cnst_v = \
    hodge_decomposition(walked_v_ER,'edge_visits')

#%% correlations between dynamics: PBC LATTICE
#edge centric
plot_pot_corr(pot_adj, pot_cnst_v, 'Potential Adjoint Model', 
              'Potential constant velocity')
plot_grad_corr(grad_adj, grad_cnst_v, 'Grad. comp. Adjoint Model', 
              'Grad. comp. constant velocity')

#node centric
plot_pot_corr(pot_cont_new, pot_cnst_v,'Potential Node Centric', 
              'Potential constant velocity')

plot_grad_corr(g_cont_new, grad_cnst_v,'Grad. comp. Node Centric', 
              'Grad. comp. constant velocity')
#%%BUILDING A LATTICE GRAPH WITH PERIODIC BOUNDARY CONDITIONS AND ADDING NODES
#AT THE CENTRE

'''-----------------TESTING WITH THE PBC MODIFIED LATTICE-------------------'''


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
    # new_edges += [((3,3+i/4), (3, 3+(i+1)/4))]
    new_edges += [((3+i/4,4), (3+(i+1)/4, 4))]
    for j in range(0, 4):
        new_edges += [((3+i/4, 3+j/4), (3+(i+1)/4, 3+j/4)),((3+i/4, 3+j/4), (3+i/4, 3+(j+1)/4)) 
                      ,((3+i/4, 3+j/4), (3+(i+1)/4, 3+(j+1)/4)), 
                      (((3+i/4, 3+(j+1)/4)), (3+(i+1)/4, 3+j/4))
                      ]
#higher degree lattice
# for i,j in PBC_lattice.nodes:
#     if (i != 3 or i!=4) and (j != 3 or j!=4):
#         if i != N-1 and j != N-1:
#             new_edges += [((i,j), (i+1,j+1))]
#         if i != N-1 and j != N and j != 0:
#             new_edges +=[((i,j), (i+1,j-1))]
PBC_lattice.remove_edges_from([((3,3),(3,4)), ((3, 3), (4,3)), ((3,4), (4,4)),
                               ((4,3), (4,4))])
PBC_lattice.add_nodes_from(new_vertexes)
PBC_lattice.add_edges_from(new_edges)
#pos = {i * N + j:(i, j) for i, j in PBC_lattice.nodes()}
pos_PBC = {i: j for i,j in enumerate(PBC_lattice.nodes)}
#labels = dict( ((i, j), i * N + j) for i, j in PBC_lattice.nodes() )
labels = {i: k for k, i in enumerate(PBC_lattice.nodes)}
new_edges_rel = [(labels[edge[0]], labels[edge[1]]) for edge in new_edges]
nx.relabel_nodes(PBC_lattice, labels, copy=False)

PBC_lattice = nx.DiGraph(PBC_lattice)
out_edges = [edge for edge in PBC_lattice.edges if edge[0]
    > edge[1]]  # removing all outward edges
PBC_lattice.remove_edges_from(out_edges)
#%%
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from matplotlib.patches import ConnectionPatch
import networkx as nx

# Create the figure and main axes
fig, ax = plt.subplots()

# Original plot
nx.draw_networkx(PBC_lattice, pos=pos_PBC, with_labels=False, node_size=30, ax=ax)

# Create an inset axes for the zoomed-in plot
zoom_ax = zoomed_inset_axes(ax, 2, loc=1)

# Zoomed-in inner plot
zoom_range = ((2.8, 4.2), (2.8, 4.2))  # Range for x and y coordinates respectively
nx.draw_networkx(PBC_lattice, pos=pos_PBC, with_labels=False, node_size=40, ax=zoom_ax)
zoom_ax.set_xlim(zoom_range[0])
zoom_ax.set_ylim(zoom_range[1])

# Mark the zoomed-in region on the original plot
mark_inset(ax, zoom_ax, loc1=2, loc2=4, fc="none", ec="0.5")

# Show the plot
plt.tight_layout()
plt.show()

#%% curl prediction
Dt = 100
n_walk = 100
v = 1

#theoretical
trans_rate_PBC_nc = build_trans_rates_matrix(PBC_lattice.copy(), pos_PBC, v)

PBC_th, solution_PBC, curl_theo = solve_continuous_rw_flow(PBC_lattice.copy(), trans_rate_PBC_nc, Dt, n_walk)

rev_PBC_th = reverse_negative_edges(PBC_th)

grad_PBC_nc_th, sol_PBC_nc_th, har_PBC_nc_th, pot_PBC_nc_th, div_PBC_nc_th = \
    hodge_decomposition(rev_PBC_th, 'edge_visits')

plot_hodge(rev_PBC_th, grad_PBC_nc_th, sol_PBC_nc_th, har_PBC_nc_th, pot_PBC_nc_th,
           div_PBC_nc_th, pos_PBC)
#%%
#simulation
walked_nc_PBC, paths_PBC_nc = node_centric(PBC_lattice.copy(), Dt, v, pos_PBC, n_walk)
grad_PBC_nc, sol_PBC_nc, har_PBC_nc, pot_PBC_nc, div_PBC_nc = \
    hodge_decomposition(walked_nc_PBC, 'edge_visits')

print('curl theo', curl_theo)
#%% NODE CENTRIC RW ON PBC MOD LATTICE
Dt = 70
n_walk = 20
v = 1

# SIMULATION
grad_comp_list = []
sol_comp_list = []
har_comp_list = []
for i in range(200):
    walked_nc_PBC, paths_PBC_nc = node_centric(PBC_lattice.copy(), Dt, v, pos_PBC, n_walk)
    grad_PBC_nc, sol_PBC_nc, har_PBC_nc, pot_PBC_nc, div_PBC_nc = \
        hodge_decomposition(walked_nc_PBC, 'edge_visits')
    grad_comp_list.append(grad_PBC_nc)
    sol_comp_list.append(sol_PBC_nc)
    har_comp_list.append(har_PBC_nc)
#%%   HISTOGRAMS
from scipy.stats import norm, poisson

trans_rate_PBC_nc = build_trans_rates_matrix(PBC_lattice.copy(), pos_PBC, v)#, new_edges_rel, 1, 1)
PBC_th, solution_PBC, _ = solve_continuous_rw_flow(PBC_lattice.copy(), trans_rate_PBC_nc, Dt, n_walk)
g_PBC_th, s_PBC_th, h_PBC_th, pot_PBC_th, div_PBC_th = \
    hodge_decomposition(PBC_th, 'edge_visits')

def fit_function(k, lamb):
    # The parameter lamb will be used as the fit parameter
    return poisson.pmf(k, lamb)

n_bins = 32-19
edge = (28,29)
hist_ls_g = [comp[edge] for comp in grad_comp_list]
hist_ls_s = [comp[edge] for comp in sol_comp_list]
hist_ls_h = [comp[edge] for comp in har_comp_list]
plt.figure()
n, bins_h, _ = plt.hist(hist_ls_h, n_bins, facecolor='violet', alpha=0.75,density=True)
(mu_h, sigma_h) = norm.fit(hist_ls_h)
y_h = norm.pdf(bins_h, mu_h, sigma_h)

n, bins_s, _ = plt.hist(hist_ls_s, n_bins, facecolor='violet', alpha=0.75,density=True)
(mu_s, sigma_s) = norm.fit(hist_ls_s)
y_s = norm.pdf(bins_s, mu_s, sigma_s)

n, bins_g, _ = plt.hist(hist_ls_g, n_bins, facecolor='violet', alpha=0.75,density=True)
(mu_g, sigma_g) = norm.fit(hist_ls_g)
y_g = norm.pdf(bins_g, mu_g, sigma_g)

x_poiss = np.arange(19,32)

# calculate bin centers
middles_bins = (bins_g[1:] + bins_g[:-1]) * 0.5
parameters, cov_matrix = curve_fit(fit_function, middles_bins, n)

fig, ax = plt.subplots(1,3, figsize = (10,4))
ax[0].tick_params(axis='both', which='major', labelsize=12)
sns.histplot(hist_ls_g, bins = n_bins, color = (176/255, 5/255, 44/255, 0.3), ax=ax[0], 
             stat = 'density')
ax[0].axvline(g_PBC_th[edge], 0, 25, linestyle = '--', color = 'black', 
              label = 'Analytical : '+str(round(g_PBC_th[edge], 3)))
ax[0].set_xlabel(r'$\omega_g$', fontsize = 12)
ax[0].set_ylabel('frequency', fontsize = 12)
ax[0].plot(x_poiss, fit_function(x_poiss, *parameters), 'r-', label = r'$\mu=%.3f,\ \sigma=%.3f$' %(mu_g, sigma_g))
ax[0].legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.25))


sns.histplot(hist_ls_s, bins = n_bins, color = (35/255, 89/255, 49/255, 0.3), ax=ax[1],
             stat = 'density')
ax[1].set_xlabel(r'$\omega_s$', fontsize = 12)
ax[1].set_ylabel('frequency', fontsize = 12)
ax[1].plot(bins_s, y_s, 'r-', label = r'$\mu=%.3f,\ \sigma=%.3f$' %(mu_s, sigma_s))
ax[1].legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.2))
ax[1].tick_params(axis='both', which='major', labelsize=12)


sns.histplot(hist_ls_h, bins = n_bins, color = (14/255, 12/255, 100/255, 0.3), ax=ax[2],
             stat = 'density')
ax[2].set_xlabel(r'$\omega_h$', fontsize = 12)
ax[2].set_ylabel('frequency', fontsize = 12)
ax[2].plot(bins_h, y_h, 'r-', label = r'$\mu=%.3f,\ \sigma=%.3f$' %(mu_h, sigma_h))
ax[2].legend(loc = 'upper center', bbox_to_anchor = (0.5, 1.2))
ax[2].tick_params(axis='both', which='major', labelsize=12)

plt.tight_layout()

#%% distribution of the means
means_g = []
means_s = []
means_h = []

for edge in nx.to_undirected(PBC_lattice).edges:
    if edge not in PBC_lattice.edges():
        edge_1 = (edge[1], edge[0])
    else:
        edge_1 = edge
    hist_ls_g = [comp[edge_1] for comp in grad_comp_list]
    hist_ls_s = [comp[edge_1] for comp in sol_comp_list]
    hist_ls_h = [comp[edge_1] for comp in har_comp_list]
    
    (mu_g, sigma_g) = norm.fit(hist_ls_g)
    means_g.append([dists[edge], mu_g])
    (mu_h, sigma_h) = norm.fit(hist_ls_h)
    means_s.append([dists[edge], mu_s])
    (mu_s, sigma_s) = norm.fit(hist_ls_s)
    means_h.append([dists[edge], mu_h])

fig, ax = plt.subplots(1,3, figsize = (10,4))
ax[0].scatter(*zip(*means_g))
ax[0].set_xlabel('Edge Lengths')
ax[0].set_ylabel('Mean of each edge (gradient)')

ax[1].scatter(*zip(*means_s))
ax[1].set_xlabel('Edge Lengths')
ax[1].set_ylabel('Mean of each edge (solenoidal)')

ax[2].scatter(*zip(*means_h))
ax[2].set_xlabel('Edge Lengths')
ax[2].set_ylabel('Mean of each edge (harmonic)')
#%%
short_1_lsp = np.linspace(0.25, 1, num = 1000)
short_2_lsp = np.linspace(0.3535533905932738, 1, num = 1000)

#finding the nodes that are inside the cluster
inner_nodes = []
outer_nodes = []
for node, pos in pos_PBC.items():
    if 3 <= pos[0] <=4 and 3 <= pos[1] <= 4:
        inner_nodes.append(node)
    else:
        outer_nodes.append(node)
#%%
# THEORETICAL 
Dt = 100
n_walk = 80
v = 1
results_avg_pot = []
for short_1, short_2 in zip(short_1_lsp, short_2_lsp):
    theta = (0.25/short_1-0.25)/(1-0.25)
    trans_rate_PBC_nc = build_trans_rates_matrix(PBC_lattice.copy(), pos_PBC, v, new_edges_rel, short_1, short_2)
    PBC_th, solution_PBC, _ = solve_continuous_rw_flow(PBC_lattice.copy(), trans_rate_PBC_nc, Dt, n_walk)
    g_PBC_th, s_PBC_th, h_PBC_th, pot_PBC_th, div_PBC_th = \
        hodge_decomposition(PBC_th, 'edge_visits')
    
    if short_1 == short_1_lsp[0]:
        PBC_ini = PBC_th
        g_PBC_ini = g_PBC_th
        pot_PBC_ini = pot_PBC_th
        div_PBC_ini = div_PBC_th
        
    elif short_1 == short_1_lsp[-1]:
        PBC_fin = PBC_th
        g_PBC_fin = g_PBC_th
        pot_PBC_fin = pot_PBC_th
        div_PBC_fin = div_PBC_th
        
    elif theta > 0.446 and theta < 0.448:
        PBC_zero = PBC_th
        g_PBC_zero = g_PBC_th
        pot_PBC_zero = pot_PBC_th
        div_PBC_zero = div_PBC_th
    avg_inner_pot = np.mean([pot_PBC_th[in_node] for in_node in inner_nodes])
    avg_outer_pot = np.mean([pot_PBC_th[out_node] for out_node in outer_nodes])
    results_avg_pot.append([theta, avg_outer_pot - avg_inner_pot])
#%% PLOT OF INNER-OUTER POT. DIFFERENCE
fig,ax =plt.subplots(figsize = (8,6))
ax.plot(*zip(*results_avg_pot), color = 'black')
ax.set_ylabel(r'$\Delta V_{io}$', fontsize = 18)
ax.set_xlabel(r'$\theta = \frac{d_{mod}-d_{real}}{d_{lattice}-d_{real}}$',fontsize = 18)
ax.tick_params('both', labelsize = 18)
ax.fill_between(np.linspace(0,1,50), 0, 250, color = (5/255, 51/255, 177/255, 0.3), label = 'The centre is repulsive')
ax.fill_between(np.linspace(0,1,50), 0, -200, color = (176/255, 5/255, 44/255, 0.3), label = 'The centre is attractive')
ax.axhline(y = 0, color = 'gray', linestyle = '--')
plt.xlim([0, 1])
plt.ylim([-200, 250])
ax.legend(fontsize = 15)
plt.tight_layout()
#%%
print((0.25/short_1_lsp-0.25)/(1-0.25))

#%%
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import networkx as nx

plt.figure(figsize=(16, 14))
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1], wspace=0.2)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(212)

ax1.axes.set_aspect('equal')
ax2.axes.set_aspect('equal')
ax3.axes.set_aspect('equal')

percentile_g = np.percentile(list(g_PBC_ini.values()), 95)
color_g = np.array(list(g_PBC_ini.values()))
# plotting edges with color gradient

color_pot = list(pot_PBC_ini.values())
cmap_pot = plt.cm.PRGn
vmax_pot = int(np.max(color_pot))
vmin_pot = int(np.min(color_pot))

colors = np.linspace(0, percentile_g)
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

nx.draw_networkx_nodes(PBC_ini, pos=pos_PBC, label=None,
                       node_size=50, node_color=color_pot, cmap=cmap_pot,
                       vmin=vmin_pot, vmax=vmax_pot, ax=ax1)

nx.draw_networkx_edges(PBC_ini, pos=pos_PBC, label=None, edge_color=color_g,
                      edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                      arrowsize=5, node_size=50, ax=ax1)

# sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
# sm._A = []
# cbar_ax1 = plt.colorbar(sm, ax=ax1, shrink=0.4)
# cbar_ax1.set_label(r'$\omega_g$', fontsize=15)
# cbar_ax1.ax.tick_params(labelsize=18)

# sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot, vmax=vmax_pot))
# sm2._A = []
# cbar_ax2 = plt.colorbar(sm2, ax=ax1, location='right', shrink=0.4)
# cbar_ax2.set_label(r'Node potentials', fontsize=12)
# cbar_ax2.ax.tick_params(labelsize=12)

#-------- middle plot ------
# percentile_g = np.percentile(list(g_PBC_zero.values()), 95)
color_g = np.array(list(g_PBC_zero.values()))
# plotting edges with color gradient

color_pot = list(pot_PBC_zero.values())
# cmap_pot = plt.cm.PRGn
# vmax_pot = int(np.max(color_pot))
# vmin_pot = int(np.min(color_pot))

# colors = np.linspace(0, percentile_g)
# cmap = plt.cm.Oranges
# vmin = min(colors)
# vmax = max(colors)

nx.draw_networkx_nodes(PBC_zero, pos=pos_PBC, label=None,
                       node_size=50, node_color=color_pot, cmap=cmap_pot,
                       vmin=vmin_pot, vmax=vmax_pot, ax=ax2)

nx.draw_networkx_edges(PBC_zero, pos=pos_PBC, label=None, edge_color=color_g,
                      edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                      arrowsize=5, node_size=50, ax=ax2)

# sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
# sm._A = []
# cbar_ax1 = plt.colorbar(sm, ax=ax2, shrink=0.4)
# cbar_ax1.set_label(r'$\omega_g$', fontsize=15)
# cbar_ax1.ax.tick_params(labelsize=18)

# sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot, vmax=vmax_pot))
# sm2._A = []
# cbar_ax2 = plt.colorbar(sm2, ax=ax2, location='right', shrink=0.4)
# cbar_ax2.set_label(r'Node potentials', fontsize=12)
# cbar_ax2.ax.tick_params(labelsize=12)

#-------- right plot ------
percentile_g = np.percentile(list(g_PBC_fin.values()), 95)
color_g = np.array(list(g_PBC_fin.values()))
# plotting edges with color gradient

color_pot = list(pot_PBC_fin.values())
cmap_pot = plt.cm.PRGn
vmax_pot = int(np.max(color_pot))
vmin_pot = int(np.min(color_pot))

colors = np.linspace(0, percentile_g)
cmap = plt.cm.Oranges
vmin = min(colors)
vmax = max(colors)

nx.draw_networkx_nodes(PBC_fin, pos=pos_PBC, label=None,
                       node_size=50, node_color=color_pot, cmap=cmap_pot,
                       vmin=vmin_pot, vmax=vmax_pot, ax=ax3)

nx.draw_networkx_edges(PBC_fin, pos=pos_PBC, label=None, edge_color=color_g,
                      edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                      arrowsize=5, node_size=50, ax=ax3)

plt.subplots_adjust(right=0.75)



cbar_ax1 = plt.gcf().add_axes([0.9, 0.5, 0.02, 0.4])  # Adjust the position and size as needed
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm._A = []
col_b_ax = plt.colorbar(sm, cax=cbar_ax1, orientation='vertical')
col_b_ax.set_label(r'$\omega_g$', fontsize = 15 )
cbar_ax1.tick_params(labelsize=15)

cbar_ax2 = plt.gcf().add_axes([0.77, 0.5, 0.02, 0.4])  # Adjust the position and size as needed
sm = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot, vmax=vmax_pot))
sm._A = []
col_b_ax = plt.colorbar(sm, cax=cbar_ax2, orientation='vertical')
col_b_ax.set_label('Node potentials', fontsize = 15 )
cbar_ax2.tick_params(labelsize=15)


ax4.plot(*zip(*results_avg_pot), color='black')
ax4.set_ylabel(r'$\Delta V_{io}$', fontsize=18)
ax4.set_xlabel(r'$\theta = \frac{d_{mod}-d_{real}}{d_{lattice}-d_{real}}$', fontsize=18)
ax4.tick_params('both', labelsize=18)
ax4.fill_between(np.linspace(0, 1, 50), 0, 250, color=(5 / 255, 51 / 255, 177 / 255, 0.3),
                 label='The centre is repulsive')
ax4.fill_between(np.linspace(0, 1, 50), 0, -200, color=(176 / 255, 5 / 255, 44 / 255, 0.3),
                 label='The centre is attractive')
ax4.axhline(y=0, color='gray', linestyle='--')
ax4.set_xlim([0, 1])
ax4.set_ylim([-200, 250])
ax4.legend(fontsize=15)

# plt.tight_layout()
plt.show()
#%%


#%% PLOT PBC NC
#sim
plot_hodge(walked_nc_PBC, grad_PBC_nc, sol_PBC_nc, har_PBC_nc, pot_PBC_nc, 
           div_PBC_nc, pos_PBC)
#%%
#theo
plot_hodge(PBC_th, g_PBC_th, s_PBC_th, h_PBC_th, pot_PBC_th, div_PBC_th, 
           pos_PBC)

#%% DISCRETE PBC

#SIM DISCR

steps = 100
n_walk = 100

discr_PBC, _  = digraph_walkers(PBC_lattice.copy(), steps, n_walk)

grad_discr_sim, sol_discr_sim, har_discr_sim, pot_discr_sim, div_discr_sim = \
    hodge_decomposition(discr_PBC, 'edge_visits')
#%%
#THEO DISCR
steps = 100
n_walk = 100
trans_matrix_PBC = build_trans_matrix(PBC_lattice.copy())
discr_PBC_theo, _ = discrete_rw_edge_flow(PBC_lattice.copy(), trans_matrix_PBC, 
                                         steps, n_walk)

rev_discr_PBC_theo = reverse_negative_edges(discr_PBC_theo)
grad_discr_th, sol_discr_th, har_discr_th, pot_discr_th, div_discr_th = \
    hodge_decomposition(rev_discr_PBC_theo, 'edge_visits')
#%% plot discr pbc

plot_hodge(discr_PBC, grad_discr_sim, sol_discr_sim, har_discr_sim, 
           pot_discr_sim, div_discr_sim, pos_PBC)

plot_hodge(rev_discr_PBC_theo, grad_discr_th, sol_discr_th, har_discr_th, pot_discr_th,
           div_discr_th, pos_PBC)

#%% constant velocity RW
Dt = 100
n_walk = 100
v = 1
constant_v_pbc, _ = periodic_walkers(PBC_lattice.copy(), Dt, v, pos_PBC, n_walk
                                     , N, 1)

rev_constant_v_pbc = reverse_negative_edges(constant_v_pbc)
g_PBC_cont, s_PBC_cont, h_PBC_cont, pot_PBC_cont, div_PBC_cont =\
hodge_decomposition(rev_constant_v_pbc, 'edge_visits')
#%%
plot_hodge(rev_constant_v_pbc, g_PBC_cont, s_PBC_cont, h_PBC_cont, pot_PBC_cont, 
           div_PBC_cont, pos_PBC)
#%% edge centric
Dt = 20
n_walk = 100
v = 1
rates_PBC = {edge: v/np.linalg.norm(np.array(pos_PBC[edge[0]])-np.array(pos_PBC[edge[1]]))
        for edge in PBC_lattice.edges}

# ONLY FOR PBC_LATTICE ----------------------------------------------------
lattice_gap = 1
N = 8
for i in range(N):
    rates_PBC[(i, i+N*N-N)] = v/lattice_gap
    rates_PBC[(i+N*N-N, i)] = v/lattice_gap
    if i==0:
        rates_PBC[(i,i+N-1)] = v/lattice_gap 
        rates_PBC[(i+N-1, i)] = v/lattice_gap 
    else:
        rates_PBC[(i*N, i*N+N-1)] = v/lattice_gap
        rates_PBC[(i*N+N-1, i*N)] = v/lattice_gap
#-----------------------------------------------------------------------
#THEORETICAL ADJ PBC
adj_PBC, ids_edge_PBC = build_adjoint_graph(PBC_lattice)
adj_trans_rates_PBC = build_trans_adjoint_rates_matrix(adj_PBC.copy(), 
                                                       rates_PBC, ids_edge_PBC)

PBC_theo_adj, sol_PBC_adj = solve_adjoint_ctrw_flow(PBC_lattice.copy(), adj_PBC.copy(), 
                                                  adj_trans_rates_PBC, ids_edge_PBC,
                                                  rates_PBC, Dt, n_walk)

total_fl_theo = {edge:PBC_theo_adj[edge[0]][edge[1]]['edge_visits'] for edge in 
                PBC_theo_adj.edges}

rev_PBC_theo_adj = reverse_negative_edges(PBC_theo_adj)

grad_adj_PBC, sol_adj_PBC, har_adj_PBC, pot_adj_PBC, div_adj_PBC = \
    hodge_decomposition(rev_PBC_theo_adj,'edge_visits')
plot_hodge(rev_PBC_theo_adj, grad_adj_PBC, sol_adj_PBC, har_adj_PBC, 
           pot_adj_PBC, div_adj_PBC, pos_PBC)
#%%

walked_PBC_edge, _ = adjoint_node_centric(adj_PBC.copy(), PBC_lattice.copy()
                                                     , Dt, v, n_walk, rates_PBC, 
                                                     ids_edge_PBC)

total_fl_sim = {edge:walked_PBC_edge[edge[0]][edge[1]]['edge_visits'] for edge in 
                walked_PBC_edge.edges}

grad_sim_adj_PBC, sol_sim_adj_PBC, har_sim_adj_PBC, pot_sim_adj_PBC, div_sim_adj_PBC = \
    hodge_decomposition(walked_PBC_edge,'edge_visits')

plot_hodge(walked_PBC_edge, grad_sim_adj_PBC, sol_sim_adj_PBC, har_sim_adj_PBC, 
           pot_sim_adj_PBC, div_sim_adj_PBC, pos_PBC)
#%% correlations between dynamics: PBC LATTICE
#edge centric
plot_pot_corr(pot_adj_PBC, pot_PBC_cont, 'Potential Adjoint Model', 
              'Potential constant velocity')
plot_grad_corr(grad_adj_PBC, g_PBC_cont, 'Grad. comp. Adjoint Model', 
              'Grad. comp. constant velocity')

#node centric
plot_pot_corr(pot_PBC_th, pot_PBC_cont,'Potential Node Centric', 
              'Potential constant velocity')

plot_grad_corr(g_PBC_th, g_PBC_cont,'Grad. comp. Node Centric', 
              'Grad. comp. constant velocity')

#%% Test distance vs degree: LOW DEGREE HIGH CLUSTERING

node_numb = 40
np.random.seed(1000)
clust_1 = np.random.normal(0, 0.2, size=(node_numb, 2))
clust_2 = np.random.normal(0, 2, size=(node_numb, 2))
nodes_test = np.concatenate((clust_1, clust_2))
# plt.scatter(nodes_test[:, 0], nodes_test[:, 1])
pos_dict = {i: pos for i, pos in enumerate(nodes_test)}
def compl_waxman_prob(d, L):
    '''
    

    Parameters
    ----------
    d : float
        Edge distance.
    L : float
        Max distance.

    Returns
    -------
    1-exp(-d/L).

    '''
    alpha = 2
    norm_factor = -1/(1-np.exp(-1/alpha))
    p = 0.5*norm_factor*(np.exp(-d/(alpha*L)) - 1)
    return(p)

def connect_compl_waxman(pos_nodes):
    np.random.seed(1000)
    dist_dict = {}
    edge_list = []
    for i, node_1 in enumerate(pos_nodes):
        for j, node_2 in enumerate(pos_nodes[i+1:]):
            dist = np.linalg.norm(node_1-node_2)
            dist_dict[(i,i+j+1)] = dist
    max_dist = np.max(list(dist_dict.values()))
    for edge, dist_e in dist_dict.items():
        # if np.random.choice([True, False, False]):
            # if dist_e <= 0.5*max_dist:
                if np.random.uniform(0, 1) <= compl_waxman_prob(dist_e, max_dist):
                    edge_list.append(edge)
    return(edge_list)

low_deg_high_clust = nx.DiGraph()
low_deg_high_clust.add_nodes_from(pos_dict.keys())
edge_list = connect_compl_waxman(nodes_test)
low_deg_high_clust.add_edges_from(edge_list)

plt.figure()
color_div = [nx.degree(nx.to_undirected(low_deg_high_clust), node) for node in low_deg_high_clust.nodes]
colors_div = range(int(min(color_div)), int(max(color_div)))
cmap_div = plt.cm.Oranges
vmin_div = min(colors_div)
vmax_div = round(max(color_div))

nx.draw_networkx(low_deg_high_clust, pos = pos_dict, with_labels = False,
                 node_size = 30, node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                 vmax=vmax_div)
sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                              vmax=vmax_div))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node degree')

#%% NODE-CENTRIC WALK
Dt = 20
n_walk = 120
v = 5
walked_test_nc, _ = node_centric(low_deg_high_clust.copy(), Dt, v, pos_dict, n_walk)
grad_test_nc, sol_test_nc, har_test_nc, pot_test_nc, div_test_nc = hodge_decomposition(walked_test_nc, 'edge_visits')

# THEORETICAL 

trans_rate_test_nc = build_trans_rates_matrix(low_deg_high_clust.copy(), 
                                              pos_dict, v)

test_nc_th, _ = solve_continuous_rw_flow(low_deg_high_clust, 
                                              trans_rate_test_nc, Dt, n_walk)

g_cont_test_nc, s_cont_test_nc, h_cont_test_nc, pot_cont_test_nc, div_cont_test_nc \
    = hodge_decomposition(test_nc_th, 'edge_visits')

#%% PLOT

plot_hodge(walked_test_nc, grad_test_nc, sol_test_nc, har_test_nc, pot_test_nc,
           div_test_nc, pos_dict)
#%%
fig, ax = plt.subplots()
x, y = zip(*[[nx.degree(nx.to_undirected(low_deg_high_clust), node), pot] for node, pot 
 in pot_test_nc.items()])

ax.scatter(x,y)
#%%Configuration model + realocation according to degree

deg_list = [PBC_lattice.degree()[i] for i in PBC_lattice.nodes]
CM = nx.DiGraph(nx.configuration_model(deg_list), seed = 1001)
out_edges = [edge for edge in CM.edges if edge[1]
    < edge[0]]  # removing all outward edges
CM.remove_edges_from(out_edges)
CM.remove_edges_from(nx.selfloop_edges(CM))
#Reallocate the nodes according to degree using the positions of the two gaussians.
#Lowest degree nodes at the centre higher degree outside.

node_numb = 40
np.random.seed(1001)
clust_1 = np.random.normal(0, 0.2, size=(node_numb, 2))
clust_2 = np.random.normal(0, 2, size=(node_numb, 2))
nodes_test = np.concatenate((clust_1, clust_2))
#%%
nodes_test = list(pos_PBC.values())
sorted_coordinates = sorted(nodes_test, key=lambda point: np.linalg.norm(point-np.array([3.5,3.5])))
print(sorted_coordinates, len(list(PBC_lattice.nodes())), len(deg_list))
CM_degree = sorted(CM.degree(), key = lambda x:x[1])
pos_nodes = {}

for i, coord in zip(CM_degree, sorted_coordinates):
    pos_PBC[i[0]] = coord


plt.figure()
color_div = [nx.degree(nx.to_undirected(CM), node) for node in CM.nodes]
colors_div = range(int(min(color_div)), int(max(color_div)))
cmap_div = plt.cm.Oranges
vmin_div = min(colors_div)
vmax_div = round(max(color_div))

nx.draw_networkx(CM, pos = pos_PBC, with_labels = False,
                 node_size = 30, node_color=color_div, cmap=cmap_div, vmin=vmin_div,
                 vmax=vmax_div)
sm2 = plt.cm.ScalarMappable(cmap=cmap_div, norm=plt.Normalize(vmin=vmin_div,
                                                              vmax=vmax_div))
sm2._A = []
cbar2 = plt.colorbar(sm2, location='right')
cbar2.set_label(r'Node degree')

#%% NODE-CENTRIC WALK
Dt = 20
n_walk = 120
v = 1
walked_test_nc, _ = node_centric(CM.copy(), Dt, v, pos_nodes, n_walk)
grad_test_nc, sol_test_nc, har_test_nc, pot_test_nc, div_test_nc = \
    hodge_decomposition(walked_test_nc, 'edge_visits')
    
#%% CONSTANT SPEED
Dt = 20
n_walk = 120
v = 1
walked_v_CM = node_walkers(CM.copy(), Dt, v, pos_nodes, n_walk)
grad_cnst_v, sol_cnst_v, har_cnst_v, pot_cnst_v, div_cnst_v = \
    hodge_decomposition(walked_v_CM,'edge_visits')
    
#%%
# THEORETICAL 
Dt = 20
n_walk = 120
v = 1
trans_rate_test_nc = build_trans_rates_matrix(CM.copy(), 
                                              pos_PBC, v)

test_nc_th, _ = solve_continuous_rw_flow(CM.copy(), trans_rate_test_nc, Dt, n_walk)

g_cont_test_nc, s_cont_test_nc, h_cont_test_nc, pot_cont_test_nc, div_cont_test_nc \
    = hodge_decomposition(test_nc_th, 'edge_visits')

#%% PLOT

plot_hodge(walked_v_CM, grad_cnst_v, sol_cnst_v, har_cnst_v, pot_cnst_v, div_cnst_v, pos_nodes)
#%% THEO NC
plot_hodge(test_nc_th, g_cont_test_nc, s_cont_test_nc, h_cont_test_nc, pot_cont_test_nc, div_cont_test_nc , pos_PBC)
#%%
plot_hodge(walked_test_nc, grad_test_nc, sol_test_nc, har_test_nc, pot_test_nc, div_test_nc , pos_nodes)
#%%
fig, ax = plt.subplots()
x, y = zip(*[[nx.degree(nx.to_undirected(CM), node), pot] for node, pot 
 in pot_cont_test_nc.items()])

ax.scatter(x,y)
#%% comparison with neigh distances
CM_und = PBC_lattice.to_undirected()
neigh_mean_d = {}
#distances of the edges in the graph
dists = {edge: np.linalg.norm(np.array(pos_PBC[edge[0]])-np.array(pos_PBC[edge[1]]))/v
        for edge in CM_und.edges}
for node in CM.nodes():
    neigh_dists = np.array([dists[(node,final)] if (node,final) in dists.keys()
                            else dists[(final,node)] for final 
                            in CM_und.neighbors(node)])
    neigh_mean_d[node] = np.mean(neigh_dists)





