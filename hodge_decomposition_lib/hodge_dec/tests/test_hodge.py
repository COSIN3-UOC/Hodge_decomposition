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



def test_hodge():
    
    g = nx.DiGraph()

    # SQUARE
    # g.add_nodes_from([0, 1, 2, 3])
    # g.add_edges_from([(0,1), (1,2), (2,3), (3,0)])
    # attr = {(0,1): 2, (1,2): 2, (2,3): 2, (3,0):-2}

    # 2 triangles
    g.add_nodes_from([0, 1, 2, 3, 4])
    g.add_edges_from([(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 4), (0, 4)])
    attr = {(0, 1): 3, (0, 2): 3, (1, 2): 2, (1, 3): 1, (2, 3): 1, (3, 4): -2, (0, 4): 1}
    # g.add_edges_from([(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])
    # attr = {(0, 1): 3, (0, 2): 3, (1, 2): 2, (1, 3): 1, (2, 3): 1}

    # line
    # g.add_nodes_from([0, 1, 2, 3])
    # g.add_edges_from([(0,1), (1,2), (2,3)])
    # attr = {(0,1): 2,(1,2): 2, (2,3):2}


    nx.set_edge_attributes(g, attr, name='edge_visits')
    g_rev = hd.reverse_negative_edges(g)
    grad_comp, sol_comp, har_comp, pot, div = hd.hodge_decomposition(g_rev, 'edge_visits')
    
    pos = nx.planar_layout(g)
    
    #plot
    hd.plot_hodge(g_rev, grad_comp, sol_comp, har_comp, pot, div, pos)
    
#%% Test run
test_hodge()