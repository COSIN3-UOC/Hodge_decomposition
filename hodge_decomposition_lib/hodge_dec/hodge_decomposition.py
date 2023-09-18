
import scipy
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

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
        # print('curl of graph', rot)
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

#PLOT HODGE FUNCTION

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
    percentile_pot = np.percentile(list(pot.values()), 95)
    

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
    plt.title('Original Graph 100%', fontsize = 20)

    # color_p = np.abs(np.array(list(edge_graph.values())))
    color_p = np.array(list(edge_graph.values()))
    # colors = np.linspace(0, percentile)
    # cmap = plt.cm.Oranges
    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
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
    cbar.set_label(r'$\omega$', fontsize = 18)
    cbar.ax.tick_params(labelsize=18)


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
                            node_size=8, node_color=color_pot, cmap=cmap_pot,
                            vmin=vmin_pot, vmax=vmax_pot)
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_g,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax, 
                           arrowsize = 5, node_size = 8)

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
    plt.subplot(223)
    plt.title('Solenoidal Component ' + str(round(weight_s*100, 1))+'%', fontsize = 20)
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
                           node_color='#D3D3D3')
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_s,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 4)


    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_s\right|$', fontsize = 18)
    cbar.ax.tick_params(labelsize=18)
    
    
    colors = np.linspace(0, percentile)
    cmap = plt.cm.Oranges
    vmin = min(colors)
    vmax = max(colors)

    color_h = np.array(list(har_comp.values()))
    plt.subplot(224)
    plt.title('Harmonic Component ' + str(round(weight_h*100, 1))+'%', fontsize = 20)
    nx.draw_networkx_nodes(walk_graph, pos=pos, label=None, node_size=4,
                           node_color='#D3D3D3')
    nx.draw_networkx_edges(walk_graph, pos=pos, label=None, edge_color=color_h,
                           edge_cmap=cmap, edge_vmin=vmin, edge_vmax=vmax,
                           arrowsize = 5, node_size = 4)
    sm = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'$\left|\omega_h\right|$', fontsize = 18)
    cbar.ax.tick_params(labelsize=18)
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

# Structural_ratios
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