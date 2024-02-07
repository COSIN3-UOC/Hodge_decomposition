#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  1 11:56:41 2023

@author: robertbenassai
"""


import pandas as pd
import scipy
import os
os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from sympy import LeviCivita
from itertools import combinations, permutations
from scipy.spatial import Delaunay
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import seaborn as sns
import csv
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
    path_tot = []
    #overall edge flow
    overall_edge_weights = {edge: 0 for edge in list(G.edges)}
    for i in range(0, n_walk):
        print('realisation '+str(i+1)+'/'+str(n_walk))
        path = []
        # let's add a list of n_walk random walkers
        initial_nodes = list(G.nodes)
        # initial_nodes += initial_nodes
        # path stores the position of each walker in each time step
        path.append([initial_nodes, np.zeros_like(initial_nodes)])
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

            path.append([final_nodes, np.ones_like(final_nodes)*step])

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
            
        path_tot.append(path)
        
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
        # for node in G.nodes:
        #     overall_node_weights[node] += weights[node]
        # set node value as the number of visits of each value
        # print(weights)
    # nx.set_node_attributes(G, overall_node_weights, name='weights')
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    return(G, path_tot)
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
    
    path_tot = []
    
    # calculating the time intrinsec to each edge
    
    delta_t = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]
              ]))/v for edge in G.edges}
    nx.set_edge_attributes(G, delta_t, name='dt')
    
    # upper bound for steps
    
    min_edge_time = min(delta_t.values())
    max_steps = int(Dt/min_edge_time)
    print('maximum steps allowed each realisation', max_steps)
    
    # overall_edge_weights account for all the edge passings summed for
    # n_walk iterations
    
    overall_edge_weights = dict.fromkeys(G.edges, 0)
    
    for i in range(0, n_walk):
        
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        pos_walkers = {}
        walk_iter = 0
        for node in G.nodes:
            degree_n = G.degree[node]
            for k in range(degree_n):
                pos_walkers[walk_iter] = node
                walk_iter += 1
        path = []
        # edge_weights counts the amount of times a walker has visited each 
        # edge in 1 realisation
        edge_weights = dict.fromkeys(G.edges, 0)
        
        # first let's difine the time used for each walker
        time = dict.fromkeys(pos_walkers.keys(), 0)
        path.append([list(pos_walkers.values()), list(time.values())])
        active_walkers = dict.fromkeys(pos_walkers.keys(), True)
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
            neighbors = [{n: (G[pos_walkers[walker]][n]['dt'] if (pos_walkers[
                        walker], n) in list(G.edges) else G[n][pos_walkers[
                        walker]]['dt']) for n in nx.all_neighbors(G, 
                        pos_walkers[walker])} for walker in active_ls]
                                                                  
            # now randomly move the random walker to another node only if time+edge_time
            # is < Dt
            
            for walker, goto_nodes in zip(active_ls, neighbors):
                # time of each possible edge for a given node
                
                if np.any(time[walker] + np.array(list(goto_nodes.values()))
                          <= Dt):
                    
                    # list of the possible final nodes
                    
                    possible_nodes = list(goto_nodes.keys())
                    
                    # random selection between the final nodes
                    
                    sel = random.choice(possible_nodes)
                    
                    #if the choice surpasses the allowed time we deactivate the walker
                    
                    if time[walker] + goto_nodes[sel]>Dt:
                        
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
                    
                else:
                    # desactivate the walker that has finished
                    
                    active_walkers[walker] = False
                    
            path.append([list(pos_walkers.values()), list(time.values())])
        path_tot.append(path)
        
        for edge in G.edges:
            overall_edge_weights[edge] += edge_weights[edge]
            
    nx.set_edge_attributes(G, overall_edge_weights, name='edge_visits')
    
    return(G, path_tot)

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
        # dictionary of walker index:current_node the walker is at
        
        pos_walkers = {i: i for i in G.nodes}

        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        
        edge_weights = dict.fromkeys(G.edges, 0)
        
        # first let's difine the time used for each walker
        
        time = dict.fromkeys(pos_walkers.keys(), 0)
        path.append([list(pos_walkers.values()), list(time.values())])
        
        # first let's difine the time used for each walker
        
        active_walkers = dict.fromkeys(pos_walkers.keys(), True)
        
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
            path.append([list(pos_walkers.values()), list(time.values())])
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
    # print('maximum steps allowed each realisation', max_steps)
    # overall_edge_weights account for all the edge passings summed for
    # n_walk iterations
    overall_edge_weights = dict.fromkeys(G.edges, 0)
    for i in range(0, n_walk):
        path = []
        # print('realisation '+str(i+1)+'/'+str(n_walk))
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
    Function that performs an edge centric RW. The function puts n = 
    len(list(G.nodes)) random walkers in a Digraph transitable in both 
    directions and moves them during Dt time. 

    Parameters
    ----------
    G_adj : nx.DiGraph
        Adjoint or line graph of G.
    G : nx.DiGraph
        Input graph to walk on.
    Dt : Float
        Maximum moving time.
    v : Float
        velocity of the walkers.
    n_walk : Integer
        number of realisations of the random walk during Dt max time'.
    scales : dict
        Rates (\lambda) of each edge {edge:rate}
    ids_edge : dict
        Dict with the adjoint node label and corresponding edge 
        {adjoint node: edge}

    Returns
    -------
    G: nx.DiGraph with the edge flows in an attribute called 'edge_visits'.
    path_tot: List with the path of the walkers
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
    # ini_pos = {i: i for i in np.random.choice(G_adj.nodes, 1,#len(G.nodes), 
    #                                               replace = False)}
    for i in range(0, n_walk):
        path = []
        print('realisation '+str(i+1)+'/'+str(n_walk))
        # dictionary of walker index:current_node the walker is in
        
        # pos_walkers = {i: i for i in np.random.choice(G_adj.nodes, len(G.nodes), 
        #                                               replace = False)}
        pos_walkers = {}
        
        walk_ind = 0
        for edge in G_adj.nodes:
            
            pos_walkers[walk_ind] = edge
            pos_walkers[walk_ind+1] = edge
            
            walk_ind += 2
            
        prev_pos = pos_walkers.copy()
        # edge_weights counts the amount of times a walker has visited
        # each edge in 1 realisation
        adj_edge_weights = dict.fromkeys(list(G_adj.edges), 0)
        edge_weights = dict.fromkeys(G.edges, 0)
        # first let's difine the time used for each walker
        time = dict.fromkeys(pos_walkers.keys(), 0)
        path.append([[ids_edge[ind] for ind in list(pos_walkers.values())], 
                      [time[ind] for ind in pos_walkers.keys()]])
        active_walkers = dict.fromkeys(pos_walkers.keys(), True)
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
                # interval = np.sqrt(np.random.normal(1/scales[curr_edge], 0.2)**2)
                # interval = 1/scales[curr_edge]
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
                    # if random.choice([True, False]):
                    thrs_prob = cumulative_exp(Dt-time[walker],1/scales[curr_edge])
                    if random.uniform(0, 1) <= thrs_prob: 
                        
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
        # for adj_edge, flow in adj_edge_weights.items():
        #     overall_adj_edge_weights[adj_edge] += flow
    
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
def build_trans_rates_matrix(G, pos, v):#, new_edges, short_1, short_2):
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
    # # ONLY FOR PBC_LATTICE ----------------------------------------------------
    # lattice_gap = 1
    # N = 8
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
    # for new_e in new_edges:
    #     if new_e in dists.keys():
    #         if dists[new_e] == 0.25:
    #             dists[new_e] /= short_1
    #         else:
    #             dists[new_e] /= short_2
    #     else:
    #         orig_dist = dists[(new_e[1], new_e[0])]
    #         if orig_dist == 0.25:
    #             dists[(new_e[1], new_e[0])] /= short_1
    #         else:
    #             dists[(new_e[1], new_e[0])] /= short_2

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
    stat_flow = np.zeros(np.shape(trans_rates))
    #initial probabilities
    p0 = np.ones(len(G.nodes))
    
    t = np.linspace(start=0, stop=Dt,num=20000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    # _, tri_dict = find_triangles(G)
    # prob_product = []
    # for triangle in tri_dict.keys():
    #     i = triangle[0]
    #     j = triangle[1]
    #     k = triangle[2]
    #     prob_product.append(sol[:,i]*sol[:,j]*sol[:,k]*(trans_rates[i][j]*\
    #     trans_rates[j][k]*trans_rates[k][i] - trans_rates[j][i]*\
    #     trans_rates[k][j]*trans_rates[i][k]))
    
    # prob_prod_arr = np.transpose(np.array(prob_product))
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    # curls = np.trapz(prob_prod_arr, t, axis=0)
    stationary_probs = sol[-1][:]
    for j, p_j in enumerate(flows):
        
        for i in range(len(G.nodes)):
            net_flow[i,j] += flows[i]*trans_rates[j][i] - p_j*trans_rates[i][j]
            stat_flow[i,j] = stationary_probs[i]*trans_rates[j][i] - \
            stationary_probs[j]*trans_rates[i][j]
            

    net_flow *= n_walk
    cont_new = {edge: net_flow[inode_to_iarr[edge[0]]][inode_to_iarr[edge[1]]] 
                for edge in G.edges}
    nx.set_edge_attributes(G, cont_new, name = 'edge_visits')
    return(G, sol, stat_flow)#, curls)

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
    vmax_pot = np.max(color_pot)
    vmin_pot = np.min(color_pot)


    colors = np.linspace(0, max(color_g))
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)
    # plt.title('Gradient component ' + str(round(weight_g*100, 1))+'%', fontsize = 18)
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None,
                            node_size=30, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 9, node_size = 30, width = 1.5)

    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    # cbar.set_label(r'$\omega_g$', fontsize = 18 )
    cbar.set_label(r'edge probability', fontsize = 18 )
    cbar.ax.tick_params(labelsize=18)

    sm2 = plt.cm.ScalarMappable(cmap=cmap_pot, norm=plt.Normalize(vmin=vmin_pot,
                                                                  vmax=vmax_pot))
    sm2._A = []
    cbar2 = plt.colorbar(sm2, location='right')
    cbar2.set_label(r'Node probabilities', fontsize = 18)
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
    plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = 'f'{a:.2f}'+r'$\pm$'
             +''f'{std_err:.2f}'+'x + 'f'{b:0.1f}'+'\n'+r'$r^2 = $' 
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
    plt.plot(x, a*np.array(x)+b, c = 'black', label = 'y = 'f'{a:.2f}'+r'$\pm$'
             +''f'{std_err:.2f}'+'x + 'f'{b:0.1f}'+'\n'+r'$r^2 = $' 
             +str(round(r_value**2, 2)))
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

nx.draw_networkx(erd_reny, with_labels = True, node_size = 50, pos = pos_ER)

#%%
dists = np.array([np.linalg.norm(np.array(pos_ER[edge[1]])-np.array(pos_ER[edge[0]]))
         for edge in erd_reny.edges])

n, bins, _ = plt.hist(dists, bins = 20,range=(0,1))
bin_width = bins[1] - bins[0]
# sum over number in each bin and mult by bin width, which can be factored out
integral = bin_width * sum(n[0:1])
print(integral)
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
steps = 10
n_walk = 20

discr_walk_ER, occupations_disc = digraph_walkers(erd_reny.copy(), steps, n_walk)

grad_discr_sim, sol_discr_sim, har_discr_sim, pot_discr_sim, div_discr_sim = \
    hd.hodge_decomposition(discr_walk_ER, 'edge_visits')
    
hd.plot_hodge(discr_walk_ER, grad_discr_sim, sol_discr_sim, har_discr_sim, 
              pot_discr_sim, div_discr_sim, pos_ER)
#%%
#THEO DISCR
trans_matrix = build_trans_matrix(erd_reny.copy())
discr_ER_theo, theo_occupations = discrete_rw_edge_flow(erd_reny.copy(), 
                                                        trans_matrix, steps, n_walk)

grad_discr_th, sol_discr_th, har_discr_th, pot_discr_th, div_discr_th = \
    hd.hodge_decomposition(discr_ER_theo, 'edge_visits')
    
#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS DTRW

plot_pot_corr(pot_discr_sim, pot_discr_th, 'Simulated Potential', 
              'Analytical Potential')
plot_grad_corr(grad_discr_sim, grad_discr_th, 'Simulated gradient component', 
               'Analytical gradient component')
#%% NODE-CENTRIC WALK

Dt = 200
n_walk = 200
v = 1

# THEORETICAL 
trans_rate_ER = build_trans_rates_matrix(erd_reny.copy(), pos_ER, v)
ER_th, solution_ER, curls_ER = solve_continuous_rw_flow(erd_reny.copy(), trans_rate_ER, Dt, n_walk)
g_cont_new, s_cont_new, h_cont_new, pot_cont_new, div_cont_new = hd.hodge_decomposition(ER_th, 'edge_visits')
print(np.max(curls_ER))
#%% SIMULATION

walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
grad_ER, sol_ER, har_ER, pot_ER, div_ER = hd.hodge_decomposition(walked_ER, 'edge_visits')
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
#%% Probability evolution for the node centric
import math

def get_node_probabilities_nc(G, Dt, n_intervals, paths, numb_walkers, n_walk):
    
    '''
    Get the evolution of the node probabilities for a random walk performed on the nodes

    Parameters
    ----------
    G : nx.DiDraph or Graph
        Graph on which the dynamics have been performed.
    Dt : float
        Total walking time.
    n_intervals : Int
        number of time windows in which to calculate the occupation of each node.
    paths : list
        list of lists comming from the random walk functions.
    numb_walkers : int
        Total number of walkers that move.

    Returns
    -------
    occupation_df : data-frame with the occupation probability of each node at each interval

    '''
    
    intervals = list(np.linspace(0, Dt, n_intervals))
    #dataframe where we will put the node occupation at each time interval
    occupation_df = pd.DataFrame({node:np.full_like(intervals, 0) for node 
                                  in G.nodes})
    
    #we take every walk from n_walk realizations
    for k, walk in enumerate(paths):
        print(k)
        # initialize a dataframe with the node where the walker jumps to and the 
        #time at which they happen of each walker
        occupation_evo_df = pd.DataFrame({node:np.full_like(intervals, None) for node 
                                      in range(numb_walkers)})
        #for every jump we find the interval correponding to the jump time
        for i, step in enumerate(walk):
            index =  0
            nodes = step[0]
            times = step[1]
            for node, t in zip(nodes, times):

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
                    # print(last_node)
                    occupation_evo_df.loc[n, column] = last_node
                else:
                    last_node = occupation_evo_df.loc[n, column]
                    
        #now we count how many walkers are at each node at each time step
        for i_row in range(0, len(occupation_df)):
            row = occupation_evo_df.iloc[i_row]
            occupations = dict(Counter(row))
            for key, value in occupations.items():
                occupation_df.loc[i_row, key] += value
    
    occupation_df = occupation_df.divide(n_walk*numb_walkers)
    
    return(occupation_df)
#%% DISCR PROB EVO
n_intervals = 500
numb_walkers = len(erd_reny.nodes)
occupation_df = get_node_probabilities_nc(erd_reny, Dt, n_intervals, paths, numb_walkers)

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='black', ls = '-.'),
                Line2D([0], [0], color = 'black', ls = '-')]

time = np.linspace(0, 101, 101)
color_p = sns.color_palette(palette = 'colorblind', n_colors = len(erd_reny.nodes))
plt.figure(figsize = (8,6))
theo_evo = np.array([list(element.values()) for element in theo_occupations])
sim_evo = np.array([list(element.values()) for element in occupations_disc])
print(nx.degree(erd_reny))
for i in range(len(erd_reny.nodes)):
    plt.plot(time, sim_evo[:,i]/(10*n_walk), ls = '-.', c = color_p[i])
    plt.plot(time, theo_evo[:,i]/10, ls = '-', c = color_p[i])
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
with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/RG/RG_ratio_evolution.csv', 
          'w') as file1:
    writer = csv.writer(file1)
    with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/RG/RG_abs_avg_flow.csv', 
              'w') as file2:
        writer2 = csv.writer(file2)
        for Dt in times:
            print(Dt)
            walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
            grad_ER, sol_ER, har_ER, pot_ER, div_ER = hd.hodge_decomposition(walked_ER, 'edge_visits')
            
            edge_graph = nx.get_edge_attributes(walked_ER, 'edge_visits')
            w = np.array(list(edge_graph.values()))
            wg = np.array(list(grad_ER.values()))
            ws = np.array(list(sol_ER.values()))
            wh = np.array(list(har_ER.values()))
            weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
            weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
            weight_h = np.sum(np.square(wh))/np.sum(np.square(w))
            
            mean_tot = np.mean(np.array(list(edge_graph.values()))**2)
            mean_g = np.mean(np.array(list(grad_ER.values()))**2)
            mean_s = np.mean(np.array(list(sol_ER.values()))**2)
            mean_h = np.mean(np.array(list(har_ER.values()))**2)
            mean_cycl = np.mean((np.array(list(har_ER.values()))+\
                                 np.array(list(sol_ER.values())))**2)
            writer.writerow([Dt, weight_g, weight_s, weight_h])
            writer2.writerow([Dt, mean_g, mean_s, mean_h, mean_cycl, mean_tot])
        
#%%read the results from a file

weight_ls_g = []
weight_ls_s = []
weight_ls_h = []

mean_abs_g = []
mean_abs_s = []
mean_abs_h = []
mean_abs_cycl = []
mean_abs_tot = []
times = []
times2 = []

with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/ER/ER_50_01_abs_avg_flow.csv', 'r') \
    as file:
        with open('/Users/robertbenassai/Documents/UOC/ER_50_01_abs_avg_flow.csv', 
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
                mean_abs_tot.append(float(row2[5]))
#%%
from scipy.optimize import curve_fit
def func(x, b):
    # y = b*x**a
    y = b*x
    return y

def func2(x, a):
    y = a*np.ones_like(x)
    return y
params, pcov = curve_fit(func, times[100:], weight_ls_g[100:])
#%% PLOT OF STRENGTH EVOLUTION
# st_ratio_g, st_ratio_s, st_ratio_h = structural_ratios(PBC_lattice)

perr_weights = np.sqrt(np.diag(pcov))

plt.figure(figsize= (10,8))
plt.loglog(times, weight_ls_g, color = 'b' , linestyle = '-', marker = '.', 
         label = 'Gradient strength ratio')
plt.loglog(times, np.array(weight_ls_s) + np.array(weight_ls_h), color = 'r', 
         linestyle = '-', marker = '.', label = 'Cyclic strength ratio')

# plt.hlines(st_ratio_g, 0, 400, color = 'b' , linestyle = '--', label = 'gradient structural ratio')
# plt.hlines((st_ratio_h+st_ratio_s), 0, 400, color = 'r', linestyle = '--', 
#             label = 'cyclic structural ratio')
plt.loglog(times[5:], params[1]*times[5:]**params[0], 'b--', label = r'$y = t_1^{'+
            str(round(params[0], 2))+'\pm'+str(round(perr_weights[0], 2))+'}$')
plt.loglog(times[5:], 1-params[1]*times[5:]**params[0], 'r--', label = r'$y = 1 - t_1^{'+
            str(round(params[0],2))+'\pm'+str(round(perr_weights[0], 2))+'}$')
plt.xlabel(r'Final simulation time $(t_1)$', fontsize = 22)
plt.ylabel('Strength Ratio', fontsize = 22)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.legend(loc='upper center', ncol = 2, bbox_to_anchor=(0.5, 1.15), 
           fontsize = 20, framealpha=1)
plt.tight_layout()

perr = np.sqrt(np.diag(pcov))
print(perr)

#%% PLOT MEAN SQUARED FLOW
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

params_cycl, pcov_cycl = curve_fit(func, times2, np.array(mean_abs_cycl))
params_tot, pcov_tot = curve_fit(func, times2, np.array(mean_abs_tot))
params_g, pcov_g = curve_fit(func2, times2[100:], mean_abs_g[100:])

perr_cycl = np.sqrt(np.diag(pcov_cycl))
perr_tot = np.sqrt(np.diag(pcov_tot))
print(perr_cycl)

fig, ax = plt.subplots(figsize = (10,8))
ax.plot(times2, mean_abs_g, color = 'b' , linestyle = '-', marker = '.', 
         label = r'$\left< \omega_g^2 \right >$')
ax.plot(times2, np.array(mean_abs_cycl), color = 'r', 
         linestyle = '-', marker = '.', label = r'$\left< \omega_{cycl}^2 \right >$')
ax.plot(times2, np.array(mean_abs_tot), color = 'g', 
          linestyle = '-', marker = '.', label =  r'$\left< \omega_{tot}^2 \right >$')

ax.plot(times2, func2(np.array(times2), params_g[0]), 'b--', 
         label = r'$\left < \omega^2_g \right > = 'f'{params_g[0]:.2f}'+ 
         '\pm 'f'{perr_tot[0]:.1}'+'$')
ax.plot(times2, func(np.array(times2), params_cycl[0], params_cycl[1]), 'r--', 
           label = r'$\left < \omega^2_{cycl} \right > = \:('f'{params_cycl[1]:.0f}'+
           '\pm 'f'{perr_cycl[1]:.0f})' + 't_1^{'+str(round(params_cycl[0],2))+
           '\pm'+str(round(perr_cycl[0],2))+'}$')

ax.plot(times2, func(np.array(times2), params_tot[0], params_tot[1]), 'g--', 
            label =r'$\left < \omega^2_{cycl} \right > = \:('f'{params_tot[1]:.0f}'+
            '\pm 'f'{perr_tot[1]:.0f})' + 't_1^{'+str(round(params_tot[0],2))+
            '\pm'+str(round(perr_tot[0],2))+'}$')

# func(np.array(times2), params_cycl[0])
ax.set_xlabel(r'Final simulation time $(t_1)$', fontsize = 22)
ax.set_ylabel(r'$\left < \omega^2 \right >$', fontsize = 22)
# plt.xticks(fontsize = 20)
# plt.yticks(fontsize = 20)
ax.tick_params(labelsize = 20)
plt.legend(loc = 'upper center', bbox_to_anchor=(0.5, 1.15), fontsize = 18, 
           ncols = 2, framealpha=1)
plt.tight_layout()


#%%EVOLUTION OF GRADIENT AND CYCLIC SQUARED FLOWS

n_walk = 200
v = 1
times = np.linspace(0.1, 101.1, 100)
# times = np.array(range(101))
with open('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_grad_evolution_w2.csv', 
          'w') as file1:
    writer = csv.writer(file1)
    with open('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_sol_evolution_w2.csv', 
              'w') as file2:
        writer2 = csv.writer(file2)
        with open('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_har_evolution_w2.csv', 
                  'w') as file3:
            writer3 = csv.writer(file3)
            for Dt in times:
                print(Dt)
                grad_list = np.zeros(len(erd_reny.edges))
                sol_list = np.zeros(len(erd_reny.edges))
                har_list = np.zeros(len(erd_reny.edges))
                # for i in range(20):
                walked_ER, paths = node_centric(erd_reny.copy(), Dt, v, 
                                                pos_ER, n_walk)
                # walked_ER, paths = digraph_walkers(erd_reny.copy(), Dt, 
                #                                    n_walk)
                grad_ER, sol_ER, har_ER, pot_ER, div_ER = \
                    hd.hodge_decomposition(walked_ER, 'edge_visits')
                
                grad_list += np.array([grad_ER[key]**2 for key in 
                                       sorted(grad_ER.keys())])/n_walk
                
                sol_list += np.array([sol_ER[key]**2 for key in 
                                      sorted(sol_ER.keys())])/n_walk
                
                har_list += np.array([har_ER[key]**2 for key in 
                                      sorted(har_ER.keys())])/n_walk
                
                writer.writerow([Dt] + list(grad_list))
                writer2.writerow([Dt] + list(sol_list))
                writer3.writerow([Dt] + list(har_list))
#%%
n_walk = 100
v = 1
times = np.linspace(0.1, 101.1, 300)
trans_rate_RG = build_trans_rates_matrix(erd_reny.copy(), pos_ER, v)

with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/RG/RG_50_02_grad_evolution_th_w2.csv', 
          'w') as file1:
    writer = csv.writer(file1)
    with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/RG/RG_50_02_sol_evolution_th_w2.csv', 
              'w') as file2:
        writer2 = csv.writer(file2)
        with open('/Users/robertbenassai/Documents/UOC/Evolution_MSF/RG/RG_50_02_har_evolution_th_w2.csv', 
                  'w') as file3:
            writer3 = csv.writer(file3)
            for Dt in times:
                print(Dt)
                grad_list = np.zeros(len(erd_reny.edges))
                sol_list = np.zeros(len(erd_reny.edges))
                har_list = np.zeros(len(erd_reny.edges))
                # for i in range(60):
                RG_th, solution_RG, st_RG = solve_continuous_rw_flow(
                    erd_reny.copy(), trans_rate_RG, Dt, n_walk)
                
                grad_ER, sol_ER, har_ER, pot_ER, div_ER = \
                    hd.hodge_decomposition(RG_th, 'edge_visits')
                
                grad_list += np.array([grad_ER[key]**2 for key in 
                                       sorted(grad_ER.keys())])/n_walk
                
                sol_list += np.array([sol_ER[key]**2 for key in 
                                      sorted(sol_ER.keys())])/n_walk
                
                har_list += np.array([har_ER[key]**2 for key in 
                                      sorted(har_ER.keys())])/n_walk
                
                writer.writerow([Dt] + list(grad_list))
                writer2.writerow([Dt] + list(sol_list))
                writer3.writerow([Dt] + list(har_list))
#%%READING AND PLOTTING

key_dict = {key: i+1 for i, key in enumerate(sorted(erd_reny.edges))}
print(key_dict)
#simulation

#RG
grad_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_grad_evolution_w2.csv', 
                        sep = ',', header=None)
sol_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_sol_evolution_w2.csv', 
                        sep = ',', header=None)
har_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_har_evolution_w2.csv', 
                      sep = ',', header=None)

# ER

# grad_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/ER/ER_50_01_grad_evolution.csv', 
#                         sep = ',', header=None)
# sol_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/ER/ER_50_01_sol_evolution.csv', 
#                         sep = ',', header=None)
# har_evo = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/ER/ER_50_01_har_evolution.csv', 
#                       sep = ',', header=None)

#theortical

grad_evo_th = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_grad_evolution_th_w2.csv', 
                        sep = ',', header=None)
sol_evo_th = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_sol_evolution_th_w2.csv', 
                        sep = ',', header=None)
har_evo_th = pd.read_csv('/Users/robertbenassai/Documents/UOC/project_HHD/Evolution_MSF/RG/RG_50_02_sol_evolution_th_w2.csv', 
                      sep = ',', header=None)


# plt.figure()
# plt.xlabel('total simulaiton time')
# plt.ylabel(r'$\omega^2$')
times = np.linspace(0.1, 201.1, 100)

diff_ls = []
err_ls = []
# total_comp = har_evo+sol_evo+grad_evo
# mean_tot = np.mean(np.array(total_comp)[:, 1:])
# theo_MSD = np.zeros_like(times)
# for track_edge in erd_reny.edges:
track_edge = (29,36)
total_comp = har_evo.iloc[:, key_dict[track_edge]]+sol_evo.iloc[:, key_dict[track_edge]]\
                                            +grad_evo.iloc[:, key_dict[track_edge]]
   
slope = solution_ER[-1][track_edge[0]]*trans_rate_ER[track_edge[1]]\
    [track_edge[0]]+solution_ER[-1][track_edge[1]]*trans_rate_ER\
    [track_edge[0]][track_edge[1]]

slope /= len(erd_reny.nodes)

params, pcov = curve_fit(func, times, total_comp/(len(erd_reny.nodes)))
perr_weights = np.sqrt(np.diag(pcov))
    
    # diff_ls.append([round(slope, 5), params[0]])    
    # err_ls.append(perr_weights[0])
    # theo_MSD += slope*times/erd_reny.number_of_edges()
    
plt.plot(times, har_evo.iloc[:, key_dict[track_edge]]/(len(erd_reny.nodes)), 'b', label = 'Harmonic')
plt.plot(times, grad_evo.iloc[:, key_dict[track_edge]]/(len(erd_reny.nodes)), 'g', label = 'Gradient')
plt.plot(times, sol_evo.iloc[:, key_dict[track_edge]]/(len(erd_reny.nodes)), 'r', label = 'Solenoidal')
plt.plot(times, total_comp/(len(erd_reny.nodes)), 'orange', label = 'Total')

# plt.plot(times, np.mean(np.array(har_evo)[:, 1:], axis = 1)/(len(erd_reny.nodes)), 'b', label = 'Harmonic')
# plt.plot(times, np.mean(np.array(grad_evo)[:, 1:], axis = 1)/(len(erd_reny.nodes)), 'g', label = 'Gradient')
# plt.plot(times, np.mean(np.array(sol_evo)[:, 1:], axis = 1)/(len(erd_reny.nodes)), 'r', label = 'Solenoidal')
# plt.plot(times, np.mean(np.array(total_comp)[:, 1:], axis = 1)/(len(erd_reny.nodes)), 'orange', label = 'Total')

# times = np.linspace(0.1, 101.1, 300)
# plt.plot(times,har_evo_th.iloc[:, key_dict[track_edge]], 'b--', label = 'Harmonic')
# plt.plot(times,grad_evo_th.iloc[:, key_dict[track_edge]], 'g--', label = 'Gradient')
# plt.plot(times,sol_evo_th.iloc[:, key_dict[track_edge]], 'r--', label = 'Solenoidal')

plt.plot(times, slope*times, color = 'gray', label = 'prediction: '
          +str(round(slope,2))+'t', ls = '-.')

# plt.plot(times, params[1]*times**params[0], '--', color = 'black', label = r'$y = ('\
#           +str(round(params[1], 2))+'\pm'+str(round(perr_weights[1], 2))+') t_1^{'+\
#           str(round(params[0], 2))+'\pm'+str(round(perr_weights[0], 2))+'}$')

plt.plot(times, params[0]*times, '--', color = 'black', label = r'$y = ('\
          +str(round(params[0], 2))+'\pm'+str(round(perr_weights[0], 3))+') t_1$')
plt.legend()

#%%
plt.figure()
plt.xlabel('birth/death')
plt.ylabel('fit')
plt.errorbar(*zip(*diff_ls), yerr = err_ls, ls = "None")
plt.scatter(*zip(*diff_ls), s = 10)

solution_ER[-1]
#%% ADJOINT MATRIX

'''ADJOINT RANDOM WALKS'''

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
    # iarr_to_inode = {i:node for i, node in enumerate(G_adj.nodes)}
    #array of flows through the edges
    net_flow = np.zeros((N,N))
    stat_flow = np.zeros(np.shape(trans_rates))
    #initial probabilities
    p0 = np.ones(N)*2
    # p0[np.random.choice(G_adj.nodes, len(G.nodes), replace = False)] = 1
    t = np.linspace(start=0, stop=Dt,num=20000)
    sol = odeint(func=dp_dt, y0=p0, t=t, args = (trans_rates, None))
    
    #INTEGRATION
    #The transition rates are constant, thus the integral (p_i*R_{ij} - p_j*R_{ji}) can
    #be computed by first integrating the p's and then doing the calculations
    flows = np.trapz(sol, t, axis=0)
    stationary_probs = sol[-1][:]
    #edge flow for adjacent graph
    for j, p_j in enumerate(flows):
        
        for i in range(len(G_adj.nodes)):
            net_flow[i,j] += flows[i]*trans_rates[j][i] - p_j*trans_rates[i][j]
            stat_flow[i,j] = stationary_probs[i]*trans_rates[j][i] - \
            stationary_probs[j]*trans_rates[i][j]
            
    net_flow *= n_walk
    cont_new_adj = {adj_edge: net_flow[inode_to_iarr[adj_edge[0]]][inode_to_iarr[adj_edge[1]]] 
                for adj_edge in G_adj.edges}
    final_flow = {}
    for id_edge, edge in ids_edge.items():
        
        i = edge[0]
        j = edge[1]
        #list of edges in the adjoint graph
        edges_w_flow = list(cont_new_adj.keys())
        #list of adjacent edges to "edge"
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
            # final_flow[edge] = max(total_out_fl, total_in_fl)
            # final_flow[edge] = total_out_fl
        elif total_in_fl <= 0 and total_out_fl <= 0:
            final_flow[edge] = (total_in_fl + total_out_fl)/2
            # final_flow[edge] = min(total_out_fl, total_in_fl)
            # final_flow[edge] = total_in_fl
        elif total_in_fl >= 0 and total_out_fl <= 0:
            final_flow[edge] = (total_in_fl+total_out_fl)/2
            # final_flow[edge] = 0
            
        elif total_in_fl <= 0 and total_out_fl >= 0:
            final_flow[edge] = (total_out_fl + total_in_fl)/2
            # final_flow[edge] = 0
            
    nx.set_edge_attributes(G, final_flow, name = 'edge_visits')
    return(G, sol, stat_flow)

#%%

''' EDGE CENTRIC CTRWs'''

# building the test-graph 

#ERDÖS-RÉNYI
erd_reny = nx.erdos_renyi_graph(50, 0.1, seed = 1000, directed=False)
erd_reny = erd_reny.to_directed()
out_edges = [edge for edge in erd_reny.edges if edge[1]
    < edge[0]]  # removing all outward edges
erd_reny.remove_edges_from(out_edges)
pos_ER = nx.spring_layout(erd_reny, seed = 1050)
nx.draw_networkx(erd_reny, with_labels = True, node_size = 500, pos = pos_ER)

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

#%% NxN lattice

N = 3
lattice = nx.grid_2d_graph(N,N,periodic=False)

pos_lat = {i: j for i,j in enumerate(lattice.nodes)}
#labels = dict( ((i, j), i * N + j) for i, j in PBC_lattice.nodes() )
labels = {i: k for k, i in enumerate(lattice.nodes)}
nx.relabel_nodes(lattice, labels, copy=False)

lattice = nx.DiGraph(lattice)
out_edges = [edge for edge in lattice.edges if edge[0]
    > edge[1]]  # removing all outward edges
lattice.remove_edges_from(out_edges)

nx.set_node_attributes(lattice, pos_lat, name='pos')
#%% CONVERT THE GRAPH TO THE LINE GRAPH (ADJOINT)

adj_ER, ids_edge_ER = build_adjoint_graph(lattice)
pos_adj = nx.spring_layout(adj_ER, seed = 1050)
#%% PLOT OF ORIGINAL VS ADJOINT

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
nx.draw_networkx_nodes(lattice, pos=pos_lat, ax=ax[0], **common_params)
nx.draw_networkx_edges(lattice, pos=pos_lat, width=2, arrowsize = 20, arrowstyle="-|>", 
                       connectionstyle="arc3,rad=0.0", ax=ax[0], node_size=node_size)

# Draw labels and edge labels
nx.draw_networkx_labels(lattice, pos_lat, font_size=18, font_color="black", 
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
n_walk = 200
v = 0.5

rates_ER = {edge: v/np.linalg.norm(np.array(pos_lat[edge[0]])-np.array(pos_lat[edge[1]]))
        for edge in lattice.edges}

adj_trans_rates = build_trans_adjoint_rates_matrix(adj_ER.copy(), rates_ER, ids_edge_ER)

ER_theo_adj, sol_ER_adj, stat_adj_flows = solve_adjoint_ctrw_flow(lattice.copy(), adj_ER.copy(), 
                                                  adj_trans_rates, ids_edge_ER,
                                                  rates_ER, Dt, n_walk)
total_fl_theo = {edge:ER_theo_adj[edge[0]][edge[1]]['edge_visits'] for edge in 
                ER_theo_adj.edges}

grad_adj, sol_adj, har_adj, pot_adj, div_adj = hd.hodge_decomposition(ER_theo_adj,
                                                                   'edge_visits')
hd.plot_hodge(ER_theo_adj, grad_adj, sol_adj, har_adj, pot_adj, div_adj, pos_lat)

print(np.max(stat_adj_flows), np.min(stat_adj_flows))
#%% ADJOINT SIMULATION

walked_ER_edge, path_edges_ER = adjoint_node_centric(adj_ER.copy(), lattice.copy()
                                                     , Dt, v, n_walk, rates_ER, 
                                                     ids_edge_ER)

total_fl_sim = {edge:walked_ER_edge[edge[0]][edge[1]]['edge_visits'] for edge in 
                walked_ER_edge.edges}

grad_sim_adj, sol_sim_adj, har_sim_adj, pot_sim_adj, div_sim_adj = \
    hd.hodge_decomposition(walked_ER_edge,'edge_visits')

plot_hodge(walked_ER_edge, grad_sim_adj, sol_sim_adj, har_sim_adj, pot_sim_adj, 
           div_sim_adj, pos_lat)

#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS ADJOINT CTRW

plot_pot_corr(pot_sim_adj, pot_adj, 'Simulated Potential', 'Analytical Potential')
plot_grad_corr(grad_sim_adj, grad_adj, 'Simulated Gradient component',
               'Analytical Gradient Component' )

#%% GETTING EDGE PROBABILITIES FOR THE EDGE CENTRIC
def get_edge_probabilities_ec(G, Dt, n_intervals, path, numb_walkers, n_walk):
    
    '''
    Get the evolution of the edge probabilities for a random walk performed on the edges

    Parameters
    ----------
    G : nx.DiDraph or Graph
        Graph on which the dynamics have been performed.
    Dt : float
        Total walking time.
    n_intervals : Int
        number of time windows in which to calculate the occupation of each node.
    paths : list
        list of lists comming from the random walk functions.
    numb_walkers : int
        Total number of walkers that move.

    Returns
    -------
    occupation_df : data-frame with the occupation probability of each edge at each interval

    '''
    
    intervals = list(np.linspace(0, Dt, n_intervals))
    #dataframe where we will put the noed occupation at each time interval
    occupation_df_ec = pd.DataFrame({str(edge):np.full_like(intervals, 0) for edge 
                                  in G.edges})
    
    #we take every walk from n_walk realizations
    for k, walk in enumerate(path):
        print(k)
        # initialize a dataframe with the node where the walker jumps to and the 
        #time at which they happen of each walker
        occupation_evo_df = pd.DataFrame({i:np.full_like(intervals, None) for i 
                                      in range(numb_walkers)})
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
        for i_row in range(0, len(occupation_df_ec)):
            row = occupation_evo_df.iloc[i_row]
            occupations = dict(Counter(row))
            for key, value in occupations.items():
                occupation_df_ec.loc[i_row, key] += value
    
    occupation_df_ec = occupation_df_ec.divide(n_walk*numb_walkers)
    
    return(occupation_df_ec)
#%%
n_intervals = 200
numb_walkers = 2*len(lattice.edges)
occupation_df_ec = get_edge_probabilities_ec(lattice, Dt, n_intervals, 
                                             path_edges_ER, numb_walkers, n_walk)
#%% PROBABILITY PLOT
custom_lines = [Line2D([0], [0], color='black', ls = '-.'),
                Line2D([0], [0], color = 'black', ls = '-')]
plt.figure()
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 200)
sol = pd.DataFrame(sol_ER_adj/(2*len(lattice.edges)))
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
color_p = sns.color_palette(palette = 'colorblind', n_colors = len(lattice.edges))
for i, edge in enumerate(lattice.edges):
    if i< 10:
        plt.plot(t,sol[i], c = color_p[i], ls = '-')
        plt.plot(t2, occupation_df_ec[str(edge)],c = color_p[i], ls = '--')
    
plt.legend(custom_lines, ['simulation', 'theoretical'], loc='upper center', 
           ncol = 2, bbox_to_anchor=(0.5, 1))
plt.tight_layout()

print(np.sum(sol_ER_adj[-1,:]/(2*len(lattice.edges))))

#%% GETTING NODE PROBABILITIES FOR THE EDGE CENTRIC

n_intervals = 200

def get_node_probabilities_ec(G, Dt, n_intervals, path, numb_walkers, n_walk):
    
    '''
    Get the evolution of the node probabilities for a random walk performed on the edges

    Parameters
    ----------
    G : nx.DiDraph or Graph
        Graph on which the dynamics have been performed.
    Dt : float
        Total walking time.
    n_intervals : Int
        number of time windows in which to calculate the occupation of each node.
    paths : list
        list of lists comming from the random walk functions.
    numb_walkers : int
        Total number of walkers that move.

    Returns
    -------
    occupation_df : data-frame with the occupation probability of each node at each interval
    '''
    
    intervals = list(np.linspace(0, Dt, n_intervals))
    #dataframe where we will put the node occupation at each time interval
    occupation_df_ec_n = pd.DataFrame({node:np.full_like(intervals, 0) for node 
                                  in G.nodes})
    
    #we take every walk from n_walk realizations
    for iterat, walk in enumerate(path):
        print(iterat)
        # initialize a dataframe with the node where the walker jumps to and the 
        #time at which they happen of each walker
        occupation_evo_df = pd.DataFrame({node:np.full_like(intervals, None) for node 
                                      in range(numb_walkers)})
        
        last_j = np.zeros(numb_walkers)
        last_node = list(np.zeros(numb_walkers))
        #for every jump we find the interval correponding to the jump time
        for i, step in enumerate(walk[:-1]):
            #index of the walker
            index =  0
    #       current edge
            edges = step[0]
            times = step[1]
            #next edge
            edges_1 = walk[i+1][0]
            
            for edge, t in zip(edges, times):
                #next node
                edge_1 = edges_1[index]
                    
                if i == 0:
                    j = 0
                    
                    vertex_ls = list(np.intersect1d(edge, edge_1))
                    
                    if len(vertex_ls) == 1:
                        occupation_evo_df.iloc[j, index] = vertex_ls[0]
                        last_node[index] = vertex_ls[0]
                    else:
                        print(vertex_ls, edge, edge_1)
                    
                else:
    #               finding the corresponding time interval
                    intervals.append(t)
                    intervals = sorted(intervals)
                    last_j_i = int(last_j[index])
                    j = intervals.index(t)
                    
                    if j == n_intervals:
                        j -= 1
            
                    intervals.remove(t)
                    
    #               filling the df between the two jumping times with the last edge
                    if last_j_i != j:
                        for t_it in range(last_j_i, j):
                            # print(last_j_i, j)
                            occupation_evo_df.iloc[t_it, index] = last_node[index]
                    
                    
                    vertex_ls = list(np.intersect1d(edge, edge_1))
                    
                    if len(vertex_ls) == 1:
                        occupation_evo_df.iloc[j, index] = vertex_ls[0]
                        last_node[index] = vertex_ls[0]
                    else:
                        
                        if edge_1 == edge:
                        # print('walker ended '+str((node_1, node)))
                            for t_it in range(j, n_intervals):
                                # print(t_it, n_intervals)
                                occupation_evo_df.iloc[t_it, index] = last_node[index]
    
                if i == len(walk[:-1])-1:
                    for t_it in range(j, n_intervals):
                        occupation_evo_df.iloc[t_it, index] = last_node[index]
                last_j[index] = j
                index += 1      
        
        #now we count how many walkers are at each node at each time step
        
        for i_row in range(0, len(occupation_df_ec_n)):
            row = occupation_evo_df.iloc[i_row]
            occupations = dict(Counter(row))
            for key, value in occupations.items():
                # print(key)
                occupation_df_ec_n.loc[i_row, key] += value
    
    occupation_df_ec_n = occupation_df_ec_n.divide(n_walk*numb_walkers)
    
    return(occupation_df_ec_n)

#%%
numb_walkers = 2*len(lattice.edges)

occupation_df_ec_n = get_node_probabilities_ec(G, Dt, n_intervals, path_edges_ER, 
                                               numb_walkers, n_walk)
#%% CONSTANT VELOCITY WALK
Dt = 20
n_walk = 200
v = 0.5
walked_v_ER, path_tot_cnst_v = node_walkers(lattice.copy(), Dt, v, pos_lat, n_walk)
grad_cnst_v, sol_cnst_v, har_cnst_v, pot_cnst_v, div_cnst_v = \
    hd.hodge_decomposition(walked_v_ER,'edge_visits')
#%%

plot_hodge(walked_v_ER, grad_cnst_v, sol_cnst_v, har_cnst_v, pot_cnst_v, 
           div_cnst_v, pos_lat)
total_fl_cnst_v = {edge:walked_v_ER[edge[0]][edge[1]]['edge_visits'] for edge in 
                walked_v_ER.edges}
#%% retreiving edge occupation probabilities from the constant speed RW

def get_edge_probs_nc(G, Dt, n_intervals, path, numb_walkers, n_walk):
    
    '''
    Get the evolution of the edge probabilities for a random walk performed on the nodes

    Parameters
    ----------
    G : nx.DiDraph or Graph
        Graph on which the dynamics have been performed.
    Dt : float
        Total walking time.
    n_intervals : Int
        number of time windows in which to calculate the occupation of each node.
    paths : list
        list of lists comming from the random walk functions.
    numb_walkers : int
        Total number of walkers that move.

    Returns
    -------
    occupation_df : data-frame with the occupation probability of each edge at each interval
    '''
    intervals = list(np.linspace(0, Dt, n_intervals))
    #dataframe where we will put the node occupation at each time interval
    occupation_df = pd.DataFrame({edge:np.full_like(intervals, 0) for edge 
                                  in G.edges})
    
    #we take every walk from n_walk realizations
    for iterat, walk in enumerate(path):
        print(iterat)
        # initialize a dataframe with the node where the walker jumps to and the 
        #time at which they happen of each walker
        occupation_evo_df = pd.DataFrame({node:np.full_like(intervals, None) for node 
                                      in range(numb_walkers)})
        
        last_j = np.zeros(numb_walkers)
        last_edge = list(np.zeros(numb_walkers))
        #for every jump we find the interval correponding to the jump time
        for i, step in enumerate(walk[:-1]):
            #index of the walker
            index =  0
    #       current edge
            nodes = step[0]
            times = step[1]
            #next edge
            nodes_1 = walk[i+1][0]
            times_1 = walk[i+1][1]
            
            for node, t in zip(nodes, times):
                #next node
                node_1 = nodes_1[index]
                    
                if i == 0:
                    j = 0
                    if (node, node_1) in G.edges:
                        occupation_evo_df.iloc[j, index] = str((node, node_1))
                        last_edge[index] = str((node, node_1))
                    elif (node_1, node) in G.edges:
                        occupation_evo_df.iloc[j, index] = str((node_1, node))
                        last_edge[index] = str((node_1, node))
                    else:
                        occupation_evo_df.iloc[j, index] = last_edge[index]
                        print('problem '+str((node_1, node)))
                else:
    #               finding the corresponding time interval
                    intervals.append(t)
                    intervals = sorted(intervals)
                    last_j_i = int(last_j[index])
                    j = intervals.index(t)
                    
                    if j == n_intervals:
                        j -= 1
            
                    intervals.remove(t)
                    
    #               filling the df between the two jumping times with the last edge
                    if last_j_i != j:
                        for t_it in range(last_j_i, j):
                            # print(last_j_i, j)
                            occupation_evo_df.iloc[t_it, index] = last_edge[index]
                    
                    if (node, node_1) in G.edges:
                        # print('j = ', j)
                        occupation_evo_df.iloc[j, index] = str((node, node_1))
                        last_edge[index] = str((node, node_1))
                    elif (node_1, node) in G.edges:
                        # print('j = ', j)
                        occupation_evo_df.iloc[j, index] = str((node_1, node))
                        last_edge[index] = str((node_1, node))
                    elif node_1 == node:
                        # print('walker ended '+str((node_1, node)))
                        for t_it in range(j, n_intervals):
                            # print(t_it, n_intervals)
                            occupation_evo_df.iloc[t_it, index] = last_edge[index]

                        
                    else:
                        print('edge not in graph', (node, node_1))
                        break

                if i == len(walk[:-1])-1:
                    for t_it in range(j, n_intervals):
                        occupation_evo_df.iloc[t_it, index] = last_edge[index]
                last_j[index] = j
                index += 1      

        #now we count how many walkers are at each node at each time step
        for i_row in range(0, len(occupation_df)):
            row = occupation_evo_df.iloc[i_row]
            occupations = dict(Counter(row))
            for key, value in occupations.items():
                # print(key)
                occupation_df.loc[i_row, eval(key)] += value
    
    occupation_df = occupation_df.divide(n_walk*numb_walkers)
    
    return(occupation_df)

#%%
n_intervals = 200
numb_walkers = 2*len(lattice.edges)
occupation_df = get_edge_probs_nc(lattice, Dt, n_intervals, path_tot_cnst_v, 
                                  numb_walkers, n_walk)
#%% Edge probability plot edge-centric vs constant velocity

custom_lines = [Line2D([0], [0], color='black', ls = '--'),
                Line2D([0], [0], color = 'black', ls = '-'),
                Line2D([0], [0], color = 'black', ls = '-.')]
dists = {edge: np.linalg.norm(np.array(pos_ER[edge[0]])-np.array(pos_ER[edge[1]]))/v
        for edge in lattice.edges}
plt.figure()
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 200)
sol = pd.DataFrame(sol_ER_adj/(2*len(lattice.edges)))
plt.xlabel('Time')
plt.ylabel('Edge Occupation Probabilities')
color_p = sns.color_palette(palette = 'colorblind', n_colors = len(lattice.edges))
for i, edge in enumerate(lattice.edges):
    # if i < 4:
        print(dists[edge], color_p[i])
        print(edge)
        plt.plot(t,sol[i], c = color_p[i], ls = '-')
        plt.plot(t2, np.ones_like(t2)*1/len(lattice.edges), 'black', ls = '-.')
        plt.plot(t2, np.ones_like(t2)*3/44, 'black', ls = '-.')
        plt.plot(t2, np.ones_like(t2)*5/44, 'black', ls = '-.')
        # plt.plot(t2,occupation_df_ec[str(edge)]/2, c = color_p[i], ls = '-', 
        #          label = str(edge))
        plt.plot(t2, occupation_df[edge],c = color_p[i], ls = '--', 
                 label = str(edge))
    
plt.legend(custom_lines, ['simulation cnst. vel.', 'theoretical edge-centric', 
                          'stationary probabilities'],
           loc='upper center', ncol = 1, bbox_to_anchor=(0.8, 0.9))
plt.tight_layout()

print(np.sum(sol_ER_adj[-1,:]))
print(np.sum(occupation_df.iloc[-1,:]))

#%% NODE PROBABILITIES FOR CONSTANT VELOCITY

n_intervals = 200
numb_walkers = 2*len(lattice.edges)
occupation_df_node = get_node_probabilities_nc(lattice, Dt, n_intervals, 
                                               path_tot_cnst_v, numb_walkers, 
                                               n_walk)
#%% PLOT OF NODE PROBABILITIES IN THE CONSTANT VELOCITY

plt.figure(figsize = (8,6))
t2 = np.linspace(0, Dt, 200)
plt.xlabel('time', fontsize = 18)
plt.ylabel('Node Occupation Probabilities', fontsize = 18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

custom_lines = [Line2D([0], [0], color='black', ls = '--'),
                Line2D([0], [0], color = 'black', ls = '-'),
                Line2D([0], [0], color = 'black', ls = '-.')]
deg_2 = 0
for node in lattice.nodes:
    deg_2 += nx.degree(lattice, node)**2
    
A = deg_2-2*len(list(lattice.edges))

for i, node in enumerate(lattice.nodes):
    color = color_p[i]
    deg = nx.degree(lattice, node)
    
    plt.plot(t2, occupation_df_ec_n[i],c = color, ls = '-')
    plt.plot(t2, np.ones_like(t2)*deg*(deg-1)/A, c = 'black', ls = '-.', 
             linewidth=2)
    plt.plot(t2, occupation_df_node[i],c = color, ls = '--')
plt.legend(custom_lines, ['cnst. v', 'edge-centric', 'predicted stationary edge c.'], 
           loc='upper center', 
           ncol = 3, bbox_to_anchor=(0.5, 1.1), fontsize = 15)
plt.tight_layout()

# checking if the probabilities sum up to 1

tot = 0
for node in lattice.nodes:
    deg = nx.degree(lattice, node)
    tot += deg*(deg-1)/A
print(tot)
print(np.sum(occupation_df_node.iloc[199]))

#%%

pred_stat_flows = np.zeros((len(lattice.nodes), len(lattice.nodes)))

for edge in lattice.edges:
    pred_stat_flows[edge[0]][edge[1]] = (nx.degree(lattice, edge[0])-
                                        nx.degree(lattice, edge[1]))/A

print(np.max(pred_stat_flows))
#%% correlations between dynamics
#edge centric
plot_pot_corr(pot_sim_adj, pot_cnst_v, 'Potential edge-centric Model', 
              'Potential constant velocity')
plot_grad_corr(grad_sim_adj, grad_cnst_v, 'Grad. comp. edge-centric Model', 
              'Grad. comp. constant velocity')
print(len(erd_reny.nodes), len(erd_reny.edges))
#%%
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

nx.draw_networkx(PBC_lattice, pos=pos_PBC, with_labels=True, node_size=30)

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

stat_probs = {i: prob/len(PBC_lattice.nodes) for i, prob in enumerate(solution_PBC[-1])}

rev_PBC_th = reverse_negative_edges(PBC_th)

grad_PBC_nc_th, sol_PBC_nc_th, har_PBC_nc_th, pot_PBC_nc_th, div_PBC_nc_th = \
    hd.hodge_decomposition(rev_PBC_th, 'edge_visits')

plot_hodge(rev_PBC_th, grad_PBC_nc_th, sol_PBC_nc_th, har_PBC_nc_th, stat_probs,
           div_PBC_nc_th, pos_PBC)
#%%
#simulation
walked_nc_PBC, paths_PBC_nc = node_centric(PBC_lattice.copy(), Dt, v, pos_PBC, n_walk)
grad_PBC_nc, sol_PBC_nc, har_PBC_nc, pot_PBC_nc, div_PBC_nc = \
    hd.hodge_decomposition(walked_nc_PBC, 'edge_visits')

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
        hd.hodge_decomposition(walked_nc_PBC, 'edge_visits')
    grad_comp_list.append(grad_PBC_nc)
    sol_comp_list.append(sol_PBC_nc)
    har_comp_list.append(har_PBC_nc)
#%% HISTOGRAMS
from scipy.stats import norm, poisson

trans_rate_PBC_nc = build_trans_rates_matrix(PBC_lattice.copy(), pos_PBC, v, new_edges_rel, 1, 1)
PBC_th, solution_PBC, _ = solve_continuous_rw_flow(PBC_lattice.copy(), trans_rate_PBC_nc, Dt, n_walk)
g_PBC_th, s_PBC_th, h_PBC_th, pot_PBC_th, div_PBC_th = \
    hd.hodge_decomposition(PBC_th, 'edge_visits')

def fit_function(k, lamb):
    # The parameter lamb will be used as the fit parameter
    return poisson.pmf(k, lamb)

n_bins = 32-19
edge = (65,69)
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
ax[0].plot(bins_g, y_g, 'r-', label = r'$\mu=%.3f,\ \sigma=%.3f$' %(mu_g, sigma_g))
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
        hd.hodge_decomposition(PBC_th, 'edge_visits')
    
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

discr_PBC, path_discr  = digraph_walkers(PBC_lattice.copy(), steps, n_walk)

occ_df_discr = get_node_probabilities_nc(PBC_lattice, steps, 100, path_discr, 
                                         len(PBC_lattice.nodes), n_walk)

grad_discr_sim, sol_discr_sim, har_discr_sim, pot_discr_sim, div_discr_sim = \
    hd.hodge_decomposition(discr_PBC, 'edge_visits')
#%%
#THEO DISCR
steps = 100
n_walk = 100
trans_matrix_PBC = build_trans_matrix(PBC_lattice.copy())
discr_PBC_theo, prob_evo_discr_pbc = discrete_rw_edge_flow(PBC_lattice.copy(), trans_matrix_PBC, 
                                         steps, n_walk)
stat_probs = {item[0]: item[1]/len(PBC_lattice.nodes) for item in prob_evo_discr_pbc[-1].items()}

rev_discr_PBC_theo = reverse_negative_edges(discr_PBC_theo)
grad_discr_th, sol_discr_th, har_discr_th, pot_discr_th, div_discr_th = \
    hd.hodge_decomposition(rev_discr_PBC_theo, 'edge_visits')
#%% plot discr pbc

plot_hodge(discr_PBC, grad_discr_sim, sol_discr_sim, har_discr_sim, 
            dict(occ_df_discr.iloc[-1]), div_discr_sim, pos_PBC)

print(sum(list(stat_probs.values())))

plot_hodge(rev_discr_PBC_theo, grad_discr_th, sol_discr_th, har_discr_th, stat_probs,
           div_discr_th, pos_PBC)

#%% constant velocity RW
Dt = 70
n_walk = 100
v = 1
constant_v_pbc, paths_cnst_v_pbc = periodic_walkers(PBC_lattice.copy(), Dt, v, 
                                                    pos_PBC, n_walk, N, 1)

rev_constant_v_pbc = reverse_negative_edges(constant_v_pbc)

#%%
numb_walkers = len(PBC_lattice.nodes)

occ_df_nc_pbc = get_node_probabilities_nc(PBC_lattice, Dt, 200, paths_cnst_v_pbc, 
                                          numb_walkers, n_walk)

g_PBC_cont, s_PBC_cont, h_PBC_cont, pot_PBC_cont, div_PBC_cont =\
hd.hodge_decomposition(rev_constant_v_pbc, 'edge_visits')
#%%
stat_probs_cnst_v = dict(occ_df_nc_pbc.iloc[-1])
plot_hodge(rev_constant_v_pbc, g_PBC_cont, s_PBC_cont, h_PBC_cont, stat_probs_cnst_v, 
           div_PBC_cont, pos_PBC)
#%% edge centric
Dt = 70
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

PBC_theo_adj, sol_PBC_adj, _ = solve_adjoint_ctrw_flow(PBC_lattice.copy(), adj_PBC.copy(), 
                                                  adj_trans_rates_PBC, ids_edge_PBC,
                                                  rates_PBC, Dt, n_walk)

total_fl_theo = {edge:PBC_theo_adj[edge[0]][edge[1]]['edge_visits'] for edge in 
                PBC_theo_adj.edges}

# -- getting theoretical stationary node probabilities for edge-centr.

iarr_to_edge = {i:ids_edge_PBC[node] for i, node in enumerate(adj_PBC.nodes)}

stat_arr = np.zeros((len(PBC_lattice.nodes), len(PBC_lattice.nodes)))

for i, prob in enumerate(sol_PBC_adj[-1]):
    
    stat_arr[iarr_to_edge[i][0]][iarr_to_edge[i][1]] = prob/(4*len(PBC_lattice.edges))
    
stat_node_arr = (np.sum(stat_arr, axis = 0) + np.sum(stat_arr, axis = 1))

print(np.sum(stat_node_arr))
deg_ls = [nx.degree(nx.to_undirected(PBC_lattice), node) for node in PBC_lattice.nodes]

stat_node_probs_ec = {node: prob for node, prob in enumerate(stat_node_arr)}

# ---

# stationary edge probabilities

stat_edge_probs_arr = sol_PBC_adj[-1]/(2*len(PBC_lattice.edges))

stat_edge_probs = {iarr_to_edge[i]:prob for i, prob in 
                   enumerate(stat_edge_probs_arr)}

print(np.sum(stat_edge_probs_arr))
# --

rev_PBC_theo_adj = reverse_negative_edges(PBC_theo_adj)

grad_adj_PBC, sol_adj_PBC, har_adj_PBC, pot_adj_PBC, div_adj_PBC = \
    hd.hodge_decomposition(rev_PBC_theo_adj,'edge_visits')
    
plot_hodge(rev_PBC_theo_adj, stat_edge_probs, sol_adj_PBC, har_adj_PBC, 
           stat_node_probs_ec, div_adj_PBC, pos_PBC)
#%%
plt.figure()
plt.scatter(deg_ls, stat_node_arr)
plt.xlabel('degree')
plt.ylabel('probability')

#%% SIMULATION EDGE CENTRIC ON PBC LATTICE

walked_PBC_edge, path_pbc_ec = adjoint_node_centric(adj_PBC.copy(), 
                                                    PBC_lattice.copy(), Dt, v, 
                                                    n_walk, rates_PBC, 
                                                    ids_edge_PBC)

total_fl_sim = {edge:walked_PBC_edge[edge[0]][edge[1]]['edge_visits'] for edge in 
                walked_PBC_edge.edges}

#getting node probabillities

n_intervals = 200
numb_walkers = 2*len(PBC_lattice.edges)

occ_probs_ec_pbc = get_node_probabilities_ec(PBC_lattice, Dt, n_intervals, 
                                              path_pbc_ec, numb_walkers, n_walk)
#%% PLOTTING GRAPH WITH STATIONARY PROBABILITIES

rev_PBC_sim_adj = reverse_negative_edges(walked_PBC_edge)

stat_probs_ec = dict(occ_probs_ec_pbc.iloc[-1])

grad_sim_adj_PBC, sol_sim_adj_PBC, har_sim_adj_PBC, pot_sim_adj_PBC, div_sim_adj_PBC = \
    hd.hodge_decomposition(rev_PBC_sim_adj,'edge_visits')

plot_hodge(rev_PBC_sim_adj, grad_sim_adj_PBC, sol_sim_adj_PBC, har_sim_adj_PBC, 
           stat_probs_ec, div_sim_adj_PBC, pos_PBC)

#%% EVOLUTION OF THE NODE PROBABILITIES EC VS CNST. V

plt.figure(figsize = (8,6))
t2 = np.linspace(0, Dt, 200)
plt.xlabel('time', fontsize = 18)
plt.ylabel('Node Occupation Probabilities', fontsize = 18)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)

custom_lines = [Line2D([0], [0], color='black', ls = '--'),
                Line2D([0], [0], color = 'black', ls = '-'),
                Line2D([0], [0], color = 'black', ls = '-.')]

for i, node in enumerate(lattice.nodes):
    if i<5:
        color = color_p[i]
        deg = nx.degree(lattice, node)
        
        plt.plot(t2, occ_probs_ec_pbc[i],c = color, ls = '-')
        plt.plot(t2, np.ones_like(t2)*stat_node_probs_ec[i], c = 'black', ls = '-.', 
                 linewidth=2)
        plt.plot(t2, occ_df_nc_pbc[i],c = color, ls = '--')
    
plt.legend(custom_lines, ['cnst. v', 'edge-centric', 'predicted stationary edge c.'], 
           loc='upper center', 
           ncol = 3, bbox_to_anchor=(0.5, 1.1), fontsize = 15)

plt.tight_layout()

print(np.sum(occ_df_nc_pbc.iloc[-1]))
#%% correlations between dynamics: PBC LATTICE
#edge centric
plot_pot_corr(pot_adj_PBC, pot_sim_adj_PBC, 'Potential Adjoint Model', 
              'Potential constant velocity')
plot_grad_corr(grad_adj_PBC, grad_sim_adj_PBC, 'Grad. comp. Adjoint Model', 
              'Grad. comp. constant velocity')
#%%
#node centric
plot_pot_corr(pot_PBC_th, pot_PBC_cont,'Potential Node Centric', 
              'Potential constant velocity')

plot_grad_corr(g_PBC_th, g_PBC_cont,'Grad. comp. Node Centric', 
              'Grad. comp. constant velocity')

#%% Test distance vs degree: LOW DEGREE HIGH CLUSTERING

node_numb = 40
np.random.seed(1000)
clust_1 = np.random.normal(3.5, 0.2, size=(node_numb, 2))

#NxN lattice ----
N = 8
clust_lattice = nx.grid_2d_graph(N,N,periodic=True)

pos_lat = {i: j for i,j in enumerate(clust_lattice.nodes)}
#labels = dict( ((i, j), i * N + j) for i, j in PBC_lattice.nodes() )
labels = {i: k for k, i in enumerate(clust_lattice.nodes)}
nx.relabel_nodes(clust_lattice, labels, copy=False)

clust_lattice = nx.DiGraph(clust_lattice)
out_edges = [edge for edge in clust_lattice.edges if edge[0]
    > edge[1]]  # removing all outward edges
clust_lattice.remove_edges_from(out_edges)

nx.set_node_attributes(clust_lattice, pos_lat, name='pos')
# ---
# x,y coords of points

tri = Delaunay(clust_1)
# plt.figure()
# plt.triplot(clust_1[:, 0], clust_1[:, 1], tri.simplices)
# plt.legend()
# plt.show()

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
graph = nx.DiGraph(list(edges))
# plt.figure()
# plot graph
# dictionary of node:position
pointIDXY = dict(zip(range(len(clust_1)), clust_1))
nx.set_node_attributes(graph, pointIDXY, name='pos')

clust_lattice = nx.disjoint_union(clust_lattice, graph)

clust_lattice.add_edges_from([(28, 102), (28, 96), (27, 65), (27, 75), (36, 74), 
                              (36, 91), (35, 67), (35, 71)])
pos_lat = nx.get_node_attributes(clust_lattice, 'pos')
nx.draw_networkx(clust_lattice, pos = pos_lat)

