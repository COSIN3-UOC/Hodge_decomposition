#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:18:27 2024

@author: robertbenassai
"""

import os
os.environ["SCIPY_USE_PROPACK"] = "1"
import time
import csv
import momepy
import pandas as pd
import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import hodge_dec.hodge_decomposition as hd
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix
from scipy.sparse.linalg import svds

# WALKERS ON A DIRECTED GRAPH (walk during amount of steps)
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
# ABSORBING RANDOM WALK ON A DIRECTED GRAPH (walk during amount of steps)
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

#  OPTIMIZED NODE_WALKERS (time/n_wak = 11s for clusters of 100 points)

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


#  Node-centric random walk

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


def plot_hodge(walk_graph, grad_comp, pot, pos):

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
    # ws = np.array(list(sol_comp.values()))
    # wh = np.array(list(har_comp.values()))
    weight_g = np.sum(np.square(wg))/np.sum(np.square(w))
    # weight_s = np.sum(np.square(ws))/np.sum(np.square(w))
    # weight_h = np.sum(np.square(wh))/np.sum(np.square(w))

    # print(weight_g, weight_s, weight_h, weight_g+weight_s+weight_h)
    
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
    # vmax_pot = np.max(color_pot)
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
                            vmin=vmin_pot, vmax=vmax_pot, nodelist = list(pot.keys()))
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 5, node_size = 6, edgelist = list(grad_comp.keys()))

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
    
    # colors = np.linspace(0, percentile)
    # cmap = plt.cm.Oranges
    # vmin = min(colors)
    # vmax = max(colors)

    #color_s = np.abs(np.array(list(sol_comp.values())))
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


# MAPS, GEOPANDAS

def distr_to_nx(distr_ind:int, path_edges: str, path_distr: str, edge_attr 
                = np.empty(shape = (1)), name :str = 'length', 
                path_nodes: str = '/Users/robertbenassai/Documents/UOC/'+
                'project_HHD/xarxaneta/data/paris/shp/nodes/nodes.shp'):
    
    # edges and nodes of the whole graph
    bcn_edges = gpd.read_file(path_edges, crs="EPSG:25831")


    bcn_nodes = gpd.read_file(path_nodes,
                          crs="EPSG:25831")
    
    if distr_ind == 0: # graph of the entire city
    
        # initialize graph of the district
        distr = nx.DiGraph()
        distr.add_nodes_from(pd.to_numeric(bcn_nodes['uid'], downcast='integer'))
        
        ind_to_pos = {}
        
        for i in range(len(bcn_nodes)):
            
            pos = (bcn_nodes.loc[i,'geometry'].x, bcn_nodes.loc[i,'geometry'].y) 
            
            if np.all(pos not in ind_to_pos.values()):
                
                ind_to_pos[int(bcn_nodes.loc[i,'uid'])] = pos
                
            else:
                distr.remove_node(int(bcn_nodes.loc[i,'uid']))
        
        # ind_to_pos = {int(bcn_nodes.loc[i,'uid']):(bcn_nodes.loc[i,'geometry'].x,
        #                                       bcn_nodes.loc[i,'geometry'].y) 
        #               for i in range(len(bcn_nodes))}
    
        # for i, j in zip(bcn_edges['i'], bcn_edges['j']):
            
        #     distr.add_edge(int(i), int(j), length =\
        #                     np.linalg.norm(np.array(ind_to_pos[int(i)]) - 
        #                                   np.array(ind_to_pos[int(j)])))

                
        for e in bcn_edges['to_from']:
            if e == None:
                continue
            ind = ''
            i = ''
            for a in e:
                
                if a.isdigit():
                    ind += a
                if a == ',':
                    i = int(ind)
                    ind = ''
            j = int(ind)
            
            # print(i, j)
            
            if i == '' or j == '' or i == j: #sometimes edges only have 1 node or are self-edges
                continue
            
            if i not in distr.nodes or j not in distr.nodes:
                continue
            
            distr.add_edge(int(i), int(j), length =\
                            np.linalg.norm(np.array(ind_to_pos[int(i)]) - 
                                          np.array(ind_to_pos[int(j)])))

            i = j = ''
            
        
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
            
    
    print(list(nx.selfloop_edges(distr)))
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
    
    #checking for self loops
    self_loops = []
    print('checking for self loops and multiedges...')
    for edge in distr.edges:
        if edge[0] == edge[1]:
            print('self loop found', edge)
            self_loops.append(edge)
        #remove multiedges
        if (edge[1], edge[0]) in distr.edges:
            if edge[0]>edge[1]:# remove only one of the edges
                self_loops.append(edge)
    distr.remove_edges_from(self_loops)
    
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


    
# DISCRETE NON ABSORBING NEWMANN: ANALYTIC DISCRETE RW
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

# raising the array to the amount of steps
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
    visit_prob = trans_probabilities @ np.ones_like(deg_arr)
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

# CTRW

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


# SYSTEM OF ODES
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
    dy = R @ p
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
    
    t = np.linspace(start=0, stop=Dt,num=1000)
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


def find_triangles(G:nx.classes.DiGraph):
    """
    Returns a list of all triangles in the graph G.
    Parameters
    ----------
    G : nx.DiGraph
        Digraph.
    
    Returns
    -------
    count_Triangle : int
                     number of triangles in the graph
    triangles : dict
                dictionary with the triangle label (i1,i2,i3): list of triangle edges
    """
    triangles = {}
    count_Triangle = 0
    for node in G.nodes():
        # Get the neighbors of the node
        neighbors = set(list(G.neighbors(node)) +list(G.predecessors(node)))
        for neighbor in neighbors:
            # Get the common neighbors of the node and its neighbor
            common_neighbors = neighbors.intersection(set(list(G.neighbors(neighbor)) 
                                                          +list(G.predecessors(neighbor))))
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
                    count_Triangle+= 1
    return count_Triangle, triangles

def pseudo_inverse_svd(matrix, tol=1e-5):
    # Compute the SVD of the matrix
    K = np.min(matrix.shape)-1
    U, s, Vt = np.linalg.svd(matrix, full_matrices = False)

    # Invert the singular values, with a threshold for numerical stability
    s_inv = np.array([1/si if si > tol else 0 for si in s])
    # Compute the pseudo-inverse
    return (Vt.T * s_inv).dot(U.T)

def pseudo_inverse_svd_sparse(matrix, tol=1e-5, k=None):
    
    # Determine the number of singular values to compute
    if k is None:
        k = min(matrix.shape) - 1
    
    # Compute the SVD of the sparse matrix
    u, s, vt = svds(matrix, k=k)
    
    # Invert the singular values, with a threshold for numerical stability
    s_inv = np.sort(np.array([1/si if si > tol else 0 for si in s]))
    
    # Construct the pseudo-inverse from the SVD components
    # Note: Need to ensure multiplication with sparse/dense arrays works correctly
    pseudo_inv = (vt.T @ (np.diag(s_inv))) @(u.T)
    
    return pseudo_inv

def create_grad_arr(G):
    # Prepare lists to construct COO format matrix
    rows = []  # Edge index
    cols = []  # Node index
    data = []  # Element value (1, -1)

    node_to_idx = {node: idx for idx, node in enumerate(G.nodes())}
    edge_to_idx = {edge: idx for idx, edge in enumerate(G.edges())}

    for edge in G.edges():
        edge_idx = edge_to_idx[edge]
        
        # For edge[1] (destination node), value is 1
        rows.append(edge_idx)
        cols.append(node_to_idx[edge[1]])
        data.append(1)
        
        # For edge[0] (source node), value is -1
        rows.append(edge_idx)
        cols.append(node_to_idx[edge[0]])
        data.append(-1)

    # Number of edges and nodes
    n_edges = len(G.edges())
    n_nodes = len(G.nodes())
    
    # Create a COO-format sparse matrix
    grad_arr_sparse = coo_matrix((data, (rows, cols)), shape=(n_edges, n_nodes), 
                                 dtype=float)

    return grad_arr_sparse


def hodge_decomposition(G, attr_name):
    '''
    Performs the Hodge decomposition on a directed and weighted graph

    Parameters
    ----------
    G : nx.DiGraph
        Networkx directed graph.
    attr_name : str
        Name of the edge attribute that represent the flow to be decomposed.

    Returns
    -------
    grad_comp: dict {edge: gradient comp}
    sol_comp: dict {edge: solenoidal comp}
    har_comp: dict {edge: harmonic comp}
    pot_nodes: dict {node: potential}
    div: dict {node: divergence}
    '''

# vector of the edge attribute considered where each component is the value of the
# attr of a given edge
    g_field = np.transpose(np.array([G[edge[0]][edge[1]][attr_name] for edge
                                     in G.edges]))
    
    field_arr = nx.adjacency_matrix(G, nodelist=None, dtype=None, weight=attr_name)
    
    div_arr = np.sum(field_arr, axis = 1) - np.sum(field_arr, axis = 0)
# Computing divergence
    div = {node: div_arr[i] for i, node in enumerate(G.nodes)}
    # div = {node: 0 for node in G.nodes}


    # for node in G.nodes:
    #     for in_edge in G.in_edges(node):
    #         div[node] -= G[in_edge[0]][in_edge[1]][attr_name]
    #     for out_edge in G.out_edges(node):
    #         div[node] += G[out_edge[0]][out_edge[1]][attr_name]
    # nx.set_node_attributes(G, div, name='div')
    
    # for node in G.nodes:
    #     if np.abs(div[node] - div2[node]) > 0.000001:
    #         print(node, div[node], div2[node])
# GRADIENT OPERATOR
    # grad_arr_np = np.array([[1 if node == edge[1] else -1 if node == edge[0] else 0 for node in G.nodes] for edge in G.edges])
    # grad_arr = create_grad_arr(G)
    inc_matr = nx.incidence_matrix(G, oriented = True)
    grad_arr = inc_matr.T
    # print('gradient op defined', inc_matr == grad_arr.T)
# LAPLACIAN MATRIX
    # lap = (grad_arr.T @ grad_arr).A
    lap = nx.laplacian_matrix(G.to_undirected()).todense()

    print('lap')
    # print('checking solvability',np.linalg.cond(lap))
    # OBTAINING THE GRADIENT COMPONENT
    # compute the pseudo inverse of the laplacian to solve the syst of equations
    # computationally intensive for large problems
    # Apinv = pseudo_inverse_svd(lap)
    Apinv = np.linalg.pinv(lap)
    # apply the pseudo-inverse to the rhs vector to obtain the `solution'
    # div_arr = np.transpose(np.array(list(div.values())))
    # print('divergence\n',div_arr)

    pot_field = np.squeeze(np.asarray(Apinv @ (-div_arr)))
    # print('error in node potentials', lap.dot(pot_field)+div_arr)
    # gradient of the potential field
    # pot_field.insert(remove_index, 0)
    pot_nodes = {n: pot for n, pot in zip(G.nodes, pot_field)}
    # print('node potentials', pot_nodes)
    grad_comp_arr = grad_arr @ pot_field
    grad_comp = {edge: grad_comp_arr[i] for i, edge in enumerate(G.edges)}
    
    print('finished computing gradient component')

# SOLENOIDAL COMPONENT
# Calculating the edge-2 matrix or oriented-face incidence matrix which
# corresponds to the curl operator
    n_tri, tri_dict = find_triangles(G)
    # print('number of triangles', n_tri)
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
        curl_op = np.array(curl_op, dtype = float)
        curl_op_sparse = csc_matrix(curl_op)
        print('finished computing curl op')

        rot_arr_inv = (np.linalg.pinv(curl_op @
            curl_op.T)) @ curl_op

        # computing the curl of the graph
        # rot = curl_op.dot(g_field)
        # print('big rot arr\n',curl_op.dot(np.transpose(curl_op)))
        # print('curl_op',curl_op)
        # print('field',g_field)
        # print('solvability of curl comp', np.linalg.cond(rot_arr))
        # rot_pinv = np.linalg.pinv(rot_arr)

        # solving the system of equations
        print('computing triangle potentials')

        tri_pot = rot_arr_inv @ g_field
        
        print('finished computing triangle potentials')

        # tri_pot = np.squeeze(np.asarray(rot_pinv.dot(rot)))
        # print('error curl component',
        #       curl_op.dot(np.transpose(curl_op)).dot(tri_pot)-rot)
        # solenoidal component is delta1* (transpose of curl op) of the potential:
        g_s = np.transpose(curl_op).dot(tri_pot)
        # print('curl of graph', rot)
        sol_comp = {edge: comp for edge, comp in zip(G.edges, g_s)}

        # for edge, comp in zip(G.edges, g_s):
        #     sol_comp[edge] = comp
    else:
        g_s = np.transpose(np.zeros(np.shape(grad_comp_arr)))
        sol_comp = {edge: 0 for edge in G.edges}
        # for edge in G.edges:
        #     sol_comp[edge] = 0
    
    print('finished computing solenoidal component')

# HARMONIC COMPONENT

    if n_tri != 0:
        big_arr = grad_arr @ inc_matr + curl_op_sparse.T @ curl_op_sparse
    else:
        big_arr = grad_arr @ inc_matr
#        print(big_arr)

    g_har_vec = g_field - grad_comp_arr - g_s
    thrs = 10**(-10)
    if np.all(np.abs(big_arr @ g_har_vec) < thrs):

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

# if __name__ == '__main__':
    
#     'PARIS'
    
#     distr_ind = 0
#     path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/data/paris/shp/edges/edges.shp'
#     path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'
    
#     bcn_graph, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr)


if __name__ == '__main__':
    
    'BOSTON'
    
    distr_ind = 0
    path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/data/boston/shp/edges/edges.shp'
    path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'
    path_nodes = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/data/boston/shp/nodes/nodes.shp'
    
    bcn_graph, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr, path_nodes = path_nodes)

    print('starting walk')
    v = 1.42
    Dt = 15*60
    n_walk = 20
    
    #THEO
    trans_rates_bcn = build_trans_rates_matrix(bcn_graph, ind_to_pos, v)
    
    print('transition rates defined. Solving diff eqs...')
    
    walk_bcn, _ = solve_continuous_rw_flow(bcn_graph.copy(), trans_rates_bcn,
                                            Dt, n_walk)

    with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/boston/tot.csv', 'w') as tot_f:
        writer = csv.writer(tot_f)
        e_flow = nx.get_edge_attributes(walk_bcn, 'edge_visits')
        for edge, w in e_flow.items():
            writer.writerow([edge[0], edge[1], w])
                
    # print('walk finished. Adjusting link direction to flow...')
    

#     rev_bcn = hd.reverse_negative_edges(walk_bcn)

#     print('reverse finished')
    

#     grad_bcn, sol_bcn, har_bcn, pot_bcn, div_bcn = hd.hodge_decomposition(
#         rev_bcn, 'edge_visits')

#     with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/grad.csv', 'w') as grad_f:
#         writer = csv.writer(grad_f)
#         for edge, wg in grad_bcn.items():
#             writer.writerow([edge[0], edge[1], wg])
            
#     with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/pot.csv', 'w') as pot_f:
#         writer = csv.writer(pot_f)
#         for node, p in pot_bcn.items():
#             writer.writerow([node, p])
    
#     with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/div.csv', 'w') as div_f:
#         writer = csv.writer(div_f)
#         for node, d in div_bcn.items():
#             writer.writerow([node, d])
    
#     print('HHD finished')

#     plot_hodge(rev_bcn, grad_bcn, sol_bcn, har_bcn, pot_bcn, div_bcn, ind_to_pos)


# #%%

    # print('starting walk')
    # v = 1.42
    # Dt = 15*60
    # n_walk = 20
    
    
    # tot_flow_dict = {}

    # with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/tot.csv', 'r') as tot_f:
    #     reader = csv.reader(tot_f)
    #     for row in reader:
    #         tot_flow_dict[(int(row[0]), int(row[1]))] = float(row[2])
                
    # print('walk finished. Adjusting link direction to flow...')
    
    # nx.set_edge_attributes(bcn_graph, tot_flow_dict,'edge_visits')
    # rev_bcn = hd.reverse_negative_edges(bcn_graph)

    # print('reverse finished')
    

    # grad_bcn, sol_bcn, har_bcn, pot_bcn, div_bcn = hodge_decomposition(
    #     rev_bcn, 'edge_visits')

    # with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/grad.csv', 'w') as grad_f:
    #     writer = csv.writer(grad_f)
    #     for edge, wg in grad_bcn.items():
    #         writer.writerow([edge[0], edge[1], wg])
            
    # with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/pot.csv', 'w') as pot_f:
    #     writer = csv.writer(pot_f)
    #     for node, p in pot_bcn.items():
    #         writer.writerow([node, p])
    
    # with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/paris/div.csv', 'w') as div_f:
    #     writer = csv.writer(div_f)
    #     for node, d in div_bcn.items():
    #         writer.writerow([node, d])
    
    # print('HHD finished')

    # plot_hodge(rev_bcn, grad_bcn, sol_bcn, har_bcn, pot_bcn, div_bcn, ind_to_pos)

##%% BARCELONA

# 'ALL THE CITY'

# distr_ind = 0
# path_bcn = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/edges_clean_net_willum.shp'
# path_distr = '/Users/robertbenassai/Documents/UOC/project_HHD/BCN_UNITATS_ADM/0301040100_Districtes_UNITATS_ADM.shp'
# p_nodes = '/Users/robertbenassai/Documents/UOC/project_HHD/xarxaneta/nodes_clean_net_willum.shp'

# bcn_graph, ind_to_pos = distr_to_nx(distr_ind, path_bcn, path_distr, path_nodes = p_nodes)

# tot_flow_dict = {}
# #%%
# with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/bcn/tot_fl.csv', 'r') as tot_f:
#     reader = csv.reader(tot_f)
#     for row in reader:
#         if (int(row[0]), int(row[1])) in bcn_graph.edges:
            
#             tot_flow_dict[(int(row[0]), int(row[1]))] = float(row[2])
            
#         elif (int(row[1]), int(row[0]))  in bcn_graph.edges:
            
#             tot_flow_dict[(int(row[1]), int(row[0]))] = -float(row[2])
#         else:
#             print(int(row[0]), int(row[1]))
                        
        
# nx.set_edge_attributes(bcn_graph, tot_flow_dict,'edge_visits')
# rev_bcn = hd.reverse_negative_edges(bcn_graph)    

# grad_bcn = {}

# with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/bcn/w_g.csv', 'r') as tot_f:
#     reader = csv.reader(tot_f)
#     for row in reader:
#         if (int(row[0]), int(row[1])) in rev_bcn.edges:
            
#             grad_bcn[(int(row[0]), int(row[1]))] = float(row[2])
#         else:
#             grad_bcn[(int(row[1]), int(row[0]))] = -float(row[2])

# pot_bcn = {}

# with open('/Users/robertbenassai/Documents/UOC/project_HHD/hodge_decomposition_lib/hodge_dec/tests/decompositions/bcn/pot.csv', 'r') as tot_f:
#     reader = csv.reader(tot_f)
#     for row in reader:
#         pot_bcn[int(row[0])] = -float(row[1])


# count = 0
# for edge in rev_bcn.edges:
#     if abs(grad_bcn[edge]-rev_bcn[edge[0]][edge[1]]['edge_visits']) > 0.1:
#         count += 1
#         print(edge)
# print(count, len(rev_bcn.edges))

# #%%       
# plot_hodge(rev_bcn, grad_bcn, pot_bcn, ind_to_pos)

