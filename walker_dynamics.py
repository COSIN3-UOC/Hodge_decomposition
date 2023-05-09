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

# %% Adjoint Node-centric random walk

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
                    # final_edge = ids_edge[sel]
                    # i = curr_edge[0]
                    # j = curr_edge[1]
                    # prev_edge = ids_edge[prev_pos[walker]]
                    # if step > 0:
                    #     final_edge = ids_edge[sel]
                    #     i = curr_edge[0]
                    #     j = curr_edge[1]
                    #     prev_edge = ids_edge[prev_pos[walker]]
                        
                    #     if i in prev_edge and j in final_edge:
                    #         edge_weights[curr_edge] += 1
                    #         if step == 1:
                    #             if prev_edge[1] == i:
                    #                 edge_weights[prev_edge] += 1
                    #             else:
                    #                 edge_weights[prev_edge] -= 1
                    #     elif j in prev_edge and i in final_edge:
                    #         edge_weights[curr_edge] -= 1
                    #         if step == 1:
                    #             if prev_edge[1] == j:
                    #                 edge_weights[prev_edge] += 1
                    #             else:
                    #                 edge_weights[prev_edge] -= 1

                    # # updating walker position and time
                    # curr_ind = pos_walkers[walker]
                    # prev_pos[walker] = curr_ind
                    if (pos_walkers[walker], sel) in adj_edge_weights.keys():
                        adj_edge_weights[(pos_walkers[walker], sel)] += 1
                    else:
                        adj_edge_weights[(sel, pos_walkers[walker])] -= 1
                    
                    pos_walkers[walker] = sel
                    time[walker] += interval
                else:
                    #update the edge flow of the final edge
                    # i = curr_edge[0]
                    # j = curr_edge[1]
                    # prev_edge = ids_edge[prev_pos[walker]]
                    
                    # if i in prev_edge:
                    #     edge_weights[curr_edge] += 1
                    # else:
                    #     edge_weights[curr_edge] -= 1

                    # desactivate the walker that has finished
                    active_walkers[walker] = False
                    # print(Counter(active_walkers.values()))
            path.append([[ids_edge[ind] for ind in list(pos_walkers.values())], 
                         [time[ind] for ind in pos_walkers.keys()]])
        path_tot.append(path)
        # for edge in G.edges:
        #     overall_edge_weights[edge] += edge_weights[edge]
        for adj_edge, flow in adj_edge_weights.items():
            overall_adj_edge_weights[adj_edge] += flow
        
    #Retrieving the edge flow
    for id_edge, edge in ids_edge.items():
        
        i = edge[0]
        j = edge[1]
        #list of edges in the adjoint graph
        edges_w_flow = list(overall_adj_edge_weights.keys())
        nbrs_e = [(id_edge, neigh) if (id_edge, neigh) in edges_w_flow else 
                  (neigh, id_edge) for neigh in G_adj.neighbors(id_edge)]
        adjac_e = []
        
        for adj_edge in nbrs_e:
            edge_1 = ids_edge[adj_edge[0]]
            edge_2 = ids_edge[adj_edge[1]]
            if j not in edge_1 or i not in edge_2:
                adjac_e.append(overall_adj_edge_weights[adj_edge])
            else:
                adjac_e.append(-overall_adj_edge_weights[adj_edge])
        overall_edge_weights[edge] += 0.5*np.sum(adjac_e)
    
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
    intervals = {edge: np.linalg.norm(np.array(pos[edge[0]])-np.array(pos[edge[1]]))/v
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
        # neigh_dists = np.array([dists[(node,final)] if (node,final) in dists.keys()
        #                         else dists[(final,node)] for final 
        #                         in G_und.neighbors(node)])
        # mean_dist = np.mean(neigh_dists)
        for final in G_und.neighbors(node):
            i, j = inode_to_iarr[node], inode_to_iarr[final]
            if (node,final) in intervals.keys():
                rate = intervals[(node, final)]
            else:
                rate = intervals[(final, node)]
            trans_rates[j][i] = 1/(k_deg*rate)
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
    
    t = np.linspace(start=0, stop=Dt,num=20000)
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
# %% PLOT HODGE FUNCTION

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
    cmap_pot = plt.cm.PRGn
    vmax_pot = max([int(np.max(color_pot)), int(np.abs(np.min(color_pot)))])
    vmin_pot = -vmax_pot


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
walked_ER, paths = edge_centric(erd_reny.copy(), Dt, v, pos_ER, n_walk)
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
    print(new_edges)
    return adjoint_graph, ids_edge

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
    # for initial in range(len(list(G.nodes))):
    # p0 = np.ones(len(G.nodes))
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
        adjac_e = []
        
        for adj_edge in nbrs_e:
            edge_1 = ids_edge[adj_edge[0]]
            edge_2 = ids_edge[adj_edge[1]]
            if j not in edge_1 or i not in edge_2:
                adjac_e.append(cont_new_adj[adj_edge])
            else:
                adjac_e.append(-cont_new_adj[adj_edge])
        final_flow[edge] = 0.5*np.sum(adjac_e)
        
    nx.set_edge_attributes(G, final_flow, name = 'edge_visits')
    return(G, sol)

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
#%%
print(structural_ratios(erd_reny))
adj_ER, ids_edge_ER = build_adjoint_graph(erd_reny)

Dt = 4
n_walk = 200

plt.subplots(1,2, figsize = (8, 4))
plt.subplot(121)
plt.title('original')
nx.draw_networkx(erd_reny, pos = pos_ER)
plt.subplot(122)
plt.title('adjoint')
nx.draw_networkx(adj_ER)
plt.tight_layout()
#%%
rates_ER = {edge: v/np.linalg.norm(np.array(pos_ER[edge[0]])-np.array(pos_ER[edge[1]]))
        for edge in erd_reny.edges}


adj_trans_rates = build_trans_adjoint_rates_matrix(adj_ER.copy(), rates_ER, ids_edge_ER)

ER_theo_adj, sol_ER_adj = solve_adjoint_ctrw_flow(erd_reny, adj_ER.copy(), 
                                                  adj_trans_rates, ids_edge_ER,
                                                  rates_ER, Dt, n_walk)

grad_adj, sol_adj, har_adj, pot_adj, div_adj = hodge_decomposition(ER_theo_adj,
                                                                   'edge_visits')
plot_hodge(ER_theo_adj, grad_adj, sol_adj, har_adj, pot_adj, div_adj, pos_ER)

#%%

walked_ER_edge, path_edges_ER = adjoint_node_centric(adj_ER.copy(), erd_reny.copy()
                                                     , Dt, v, n_walk, rates_ER, 
                                                     ids_edge_ER)
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

plt.figure()
t = np.linspace(0, Dt, 20000)
t2 = np.linspace(0, Dt, 500)
sol = pd.DataFrame(sol_ER_adj/len(erd_reny.edges))
plt.xlabel('time')
plt.ylabel('Occupation probabilities')
for i, edge in enumerate(erd_reny.edges):
    color = np.random.uniform(0,1, size = 3)
    plt.plot(t,sol[i], c = color, ls = '-')
    plt.plot(t2, occupation_df[str(edge)],c = color, ls = '--')
plt.show()
plt.tight_layout()

print(np.sum(sol_ER_adj[-1,:]/len(erd_reny.edges)))
print(np.array(list(Counter(list(dict(adj_ER.degree).values())).keys()))/
      (2*len(adj_ER.edges)))

#%% PLOT OF PREDICTED VS SIMULATED POTENTIALS ADJOINT CTRW
pot_list = [(pot_sim_adj[i], pot_adj[i]) for i in pot_sim_adj.keys()]
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
#%% PLOT OF PREDICTED VS SIMULATED GRADIENT COMPONENT ADJOINT CTRW
grad_list = [(grad_sim_adj[i], grad_adj[i]) for i in grad_sim_adj.keys()]
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

