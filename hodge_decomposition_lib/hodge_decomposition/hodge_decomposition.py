
from sympy import LeviCivita
import scipy
import networkx as nx
import numpy as np

def find_triangles(G:nx.classes.digraph.DiGraph):
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

def hodge_decomposition(G:nx.classes.digraph.DiGraph, attr_name: str):
    """
    Decomposes a grah according to hodge decomposition.
    Parameters
    ----------
    G : nx.DiGraph
        Digraph with edge attributes to be decomposed.
    attr_name : string
        Name of the edge attribute to be decomposed.

    Returns
    -------
    grad_comp : dict
                gradient component.
    sol_comp : dict
                solenoidal component.
    har_comp : dict
               Harmonic component.
    pot_nodes : dict
                Node potential of the gradient component.
    div : dict
          Divergence of each node

    """
    if type(attr_name) != str():
        print(TypeError)
        exit()
#vector of the edge attribute considered where each component is the value of the
#attr of a given edge       
    g_field = np.transpose(np.array([G[edge[0]][edge[1]][attr_name] for edge 
                                     in G.edges]))
       
# Computing divergence
    div = {node:0 for node in G.nodes}
    
    for node in G.nodes:
        for in_edge in G.in_edges(node):
            div[node] -= G[in_edge[0]][in_edge[1]][attr_name]
        for out_edge in G.out_edges(node):
            div[node] += G[out_edge[0]][out_edge[1]][attr_name]
    nx.set_node_attributes(G, div, name = 'div')
    
# GRADIENT OPERATOR
    grad_arr = []
    for edge in G.edges:
        row = []
        for node in G.nodes:
            if edge[1]==node:
                row.append(1)
            elif edge[0]==node:
                row.append(-1)
            else:
                row.append(0)
        grad_arr.append(row)
        row = []
    grad_arr = np.array(grad_arr)
    print('gradient op\n',grad_arr)
#LAPLACIAN MATRIX  
    lap = np.transpose(grad_arr).dot(grad_arr)
    # print('checking solvability',np.linalg.cond(lap))
    #OBTAINING THE GRADIENT COMPONENT
    # compute the pseudo inverse of the laplacian to solve the syst of equations
    # computationally intensive for large problems
    Apinv = np.linalg.pinv(lap)
    # apply the pseudo-inverse to the rhs vector to obtain the `solution'
    div_arr = np.transpose(np.array(list(div.values())))
    #print('divergence\n',div_arr)
    
    pot_field = np.squeeze(np.asarray(Apinv.dot(-div_arr)))
    print('error in node potentials', lap.dot(pot_field)+div_arr)
    #gradient of the potential field
    # pot_field.insert(remove_index, 0)
    pot_nodes = {n:pot for n, pot in enumerate(pot_field)}
    # print('node potentials', pot_nodes)
    grad_comp_arr = np.transpose(grad_arr.dot(pot_field))
    grad_comp = {edge:grad_comp_arr[i] for i, edge in enumerate(G.edges)}
# SOLENOIDAL COMPONENT
#Calculating the edge-2 matrix or oriented-face incidence matrix which
#corresponds to the curl operator
    n_tri, tri_dict = find_triangles(G)
    print(tri_dict)
    print('number of triangles', n_tri)
    if n_tri != 0:
        curl_op = []
        #creating the curl operator:
        for triangle in tri_dict.keys():
            row = []
            for edge in G.edges:
#                print(triangle)
                if (triangle[1] == min(edge) and triangle[2] == max(edge)) or\
                (triangle[0] == min(edge) and triangle[1] == max(edge)):
                    if edge[0]< edge[1]:
                        row.append(1)
                    else:
                        row.append(-1)
                elif triangle[0] == min(edge) and triangle[2] == max(edge):
                    if edge[0]< edge[1]:
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
        rot_arr_inv = np.linalg.pinv(curl_op.dot(np.transpose(curl_op))).dot(curl_op)

#        print(curl_op)
        #computing the curl of the graph
        rot = curl_op.dot(g_field)
        print('big rot arr\n',curl_op.dot(np.transpose(curl_op)))
        # print('curl_op',curl_op)
        # print('field',g_field)
        # print('solvability of curl comp', np.linalg.cond(rot_arr))
        # rot_pinv = np.linalg.pinv(rot_arr)
        
        #solving the system of equations
        tri_pot = rot_arr_inv.dot(g_field)
        # tri_pot = np.squeeze(np.asarray(rot_pinv.dot(rot)))
        print('error curl component',
              curl_op.dot(np.transpose(curl_op)).dot(tri_pot)-rot)
        #solenoidal component is delta1* (transpose of curl op) of the potential:
        g_s = np.transpose(curl_op).dot(tri_pot)
        print('curl of graph', 
              rot)
        sol_comp = {}
        for edge, comp in zip(G.edges, g_s):
            sol_comp[edge] = comp
    else:
        g_s = np.transpose(np.zeros(np.shape(grad_comp_arr)))
        sol_comp= {}
        for edge in G.edges:
            sol_comp[edge] = 0
            
#HARMONIC COMPONENT

    if n_tri != 0:
        big_arr = grad_arr.dot(np.transpose(grad_arr)) +\
            np.transpose(curl_op).dot(curl_op)
    else:
        big_arr = grad_arr.dot(np.transpose(grad_arr))
#        print(big_arr)
    
    g_har_vec = g_field - grad_comp_arr - g_s
    print(g_har_vec)
    thrs = 10**(-10)
    if np.all(np.abs(big_arr.dot(g_har_vec)) < thrs):
        
        har_comp = {}
        for edge, har in zip(G.edges, g_har_vec):
            har_comp[edge] = har
    else:
        print('problem in the harmonic component')
        print('error of the harmonic component', big_arr.dot(g_har_vec))
        # print('divergence of the harmonic component', 
        #       grad_arr.dot(np.transpose(grad_arr)).dot(g_har_vec))
        # print('curl of the harmonic component', 
        #       curl_op.dot(g_har_vec))
        
        har_comp = {}
        for edge, har in zip(G.edges, g_har_vec):
            har_comp[edge] = har
    return grad_comp, sol_comp, har_comp, pot_nodes, div