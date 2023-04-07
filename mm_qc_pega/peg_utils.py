"""
Contains graph algorithms for use in qc-ldpc optimisation
"""
import numpy as np 
import galois
from itertools import combinations
import matplotlib.pyplot as plt

def rk_edge_local_girth_layer(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, cn_degrees = None):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for current_cn_index in range(int(enumerated_cn_max[t])):
        if t == 0:
            max_girth = np.array([-np.inf])

        if not G.has_cyclical_edge_set((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index #Keep track of path of this branch of the DFS-search tree 
            enumerated_cn_max[t+1] = current_cn_index   #In next step in DFS-search limit search to enumerated cn_s

            G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
        
            if not cn_degrees is None:
                #Limit maximum cn_degree if given by assigning minimum girth
                if len(G.get_adjecent_cn(current_cn_index)) == np.max(cn_degrees) + 1:
                    new_girth = -np.inf
                else:
                    new_girth = shortest_cycles(G, current_cn_index, current_vn_index)
            else:
                #Else calculate metric (girth) for this branch of the DFS-search tree
                new_girth = shortest_cycles(G, current_cn_index, current_vn_index)

            #Update if this edge decreased the girth of the node
            girths[t+1] = min(girths[t], new_girth)

            if max_girth[0] <= girths[t+1]:
                if t == rk-1:
                    #We have found a complete branch of the DFS search tree with currently maximal girth, meaning all nodes in the path has this potential rk-girth.
                    cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
                    max_girth[0] = girths[t+1]

                else:
                    #Branch not complete yet but still potentially better than current girth, keep searching     
                    rk_edge_local_girth_layer(G, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, cn_degrees)      

            G.remove_cyclical_edge_set(current_cn_index, current_vn_index)

def rk_edge_local_girth(G, current_vn_index, rk, cn_degrees = None):
    """
    Calculate the maximum girth possible when adding an edge from current_vn_index to each check node, with a look-ahead depth of rk. 
    """
    t = 0
    enumerated_cn_indexes = np.zeros(rk+1, dtype=int) #s in article
    enumerated_cn_max = np.zeros(rk+1, dtype=int) #u in article
    girths = np.zeros(rk+1) 
    max_girth = np.array([-np.inf])
    cn_girths = np.full(G.n_cn, -np.inf)

    enumerated_cn_max[0] = G.n_cn
    girths[0] = np.inf 

    #Calculate local girths from all variable nodes and from the circulant for i = 0, N, m-N
    rk_edge_local_girth_layer(G, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, cn_degrees)

    return max_girth, cn_girths


# Non finished approximative versions
# def rk_edge_local_girth_layer_gcd(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, vn_distances, cn_distances, cn_degrees = None):
#     for current_cn_index in range(int(enumerated_cn_max[t])):
        
#         if not G.has_cyclical_edge_set((current_cn_index, current_vn_index)):
#             enumerated_cn_indexes[t] = current_cn_index #Keep track of path of this branch of the DFS-search tree 
#             enumerated_cn_max[t+1] = current_cn_index   #In next step in DFS-search limit search to 

#             G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
        
#             if not cn_degrees is None:
#                 if len(G.get_adjecent_cn(current_cn_index)) == np.max(cn_degrees) + 1:
#                     new_girth == -np.inf
#             else:
#                 new_girth = shortest_cycles_gcd(G, current_cn_index, current_vn_index, vn_distances, cn_distances)

#             girths[t+1] = min(girths[t], new_girth)

#             if max_girth[0] <= girths[t+1]:
#                 if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
#                     cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
#                     max_girth[0] = girths[t+1]

#                 else: 
#                     vn_indexes = list(range(G.n_cn, G.n_nodes))
#                     vn_distances = shortest_distances(G, current_vn_index, vn_indexes)
                    
#                     cn_distances = {}
#                     T = np.arange(G.N)
#                     for ci in range(0, G.n_cn, G.N):
#                         cn_indexes = list(G.proto_index(ci)*G.N + G.proto_value(ci + T))
#                         cn_distances[ci] = shortest_distances(G, ci, cn_indexes)
                
#                     rk_edge_local_girth_layer(G, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, vn_distances, cn_distances)
#             else:
#                 pass        

#             G.remove_cyclical_edge_set(current_cn_index, current_vn_index)

# def rk_edge_local_girth_gcd(G, current_vn_index, rk, cn_degrees = None):
#     """
#     Calculate the maximum girth possible when adding an edge from current_vn_index to each check node, with a look-ahead depth of rk. 
#     """
#     t = 0
#     enumerated_cn_indexes = np.zeros(rk+1, dtype=int) #s in article
#     enumerated_cn_max = np.zeros(rk+1, dtype=int) #u in article
#     girths = np.zeros(rk+1) 
#     max_girth = np.array([-np.inf])
#     cn_girths = np.full(G.n_cn, -np.inf)

#     enumerated_cn_max[0] = G.n_cn
#     girths[0] = np.inf 

#     #Calculate local girths from all variable nodes and from the circulant for i = 0, N, m-N
#     vn_indexes = list(range(G.n_cn, G.n_nodes))
#     vn_distances = shortest_distances(G, current_vn_index, vn_indexes)
    
#     cn_distances = {}
#     T = np.arange(G.N)
#     for ci in range(0, G.n_cn, G.N):
#         cn_indexes = list(G.proto_index(ci)*G.N + G.proto_value(ci + T))
#         cn_distances[ci] = shortest_distances(G, ci, cn_indexes)

#     rk_edge_local_girth_layer_gcd(G, current_vn_index, rk, t, 
#                         enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, vn_distances, cn_distances)

#     return max_girth, cn_girths

def shortest_distances(G, node,  stop_nodes = []):
    """
    Shortest path to all listed stop nodes starting from node. Implemented with a BFS algorithm.
    """
    Q = [node]
    explored = set(Q)
    distance = 0
    distances = {}

    while Q and stop_nodes:
        distance += 1
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if not adjecent_node in explored:
                    explored.add(adjecent_node)
                    adjecent_nodes.append(adjecent_node)
                    
                    distances[adjecent_node] = distance

                    if adjecent_node in stop_nodes:
                        distances[adjecent_node] = distance
                        stop_nodes.remove(adjecent_node)

        Q = adjecent_nodes

    return distances

def shortest_cycles(G, node, stop_node = None):
    """
    Shortest cycles starting in node and passing through each adjecent node.

    Performs a BFS search, saving the path taken up to that path. This path is necessarily minimal.
    If a stop node is provided, the algorithm stops when the cycle containing this node is found and this cycle distance is returned.
    """
    Q = [node]
    explored = {node : []}
    distances = {}

    while Q:
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if not adjecent_node in explored:
                    explored[adjecent_node] = explored[node] + [adjecent_node]
                    adjecent_nodes.append(adjecent_node)
                
                else:
                    overlap = [n1==n2 for n1, n2 in zip(explored[adjecent_node], explored[node])]
                    if not any(overlap) and len(explored[node]) > 1 :
                        n1 = explored[node][0]
                        n2 = explored[adjecent_node][0]

                        distance = len(explored[node]) + len(explored[adjecent_node]) + 1
                        distances[n1] = distance
                        distances[n2] = distance

                        if stop_node in explored[node]:
                            return distance

        Q = adjecent_nodes
    
    if stop_node is None:
        return distances
    else:
        return np.inf

def shortest_cycles_gcd(G, cn, vn, vn_distances, cn_distances):
    """Approximation of local edge girth for short cycles"""    
    T = np.arange(G.N)
    cn_indexes = list(G.proto_index(cn)*G.N + G.proto_value(cn + T))
    vn_indexes = list(G.proto_index(vn)*G.N + G.proto_value(vn + T) + G.n_cn)
    distances = shortest_distances(G, vn, cn_indexes.copy())
    
    delta = [distances.get(key, np.inf) + 1 for key in cn_indexes]
    result = delta[0]

    for t in range(1, G.N-1):
        t1 = np.mod(G.N-t, G.N)
        cond1 = delta[t] + delta[t1]
        
        t1 = np.mod(t-cn, G.N)
        t2 = np.mod(-cn, G.N)
        temp = cn_distances.get(cn_indexes[t1], np.inf)
        if temp == np.inf:
            cond2 = np.inf
        else:   
            cond2 = vn_distances.get(vn_indexes[t], np.inf) + temp.get(vn_indexes[t2], np.inf) + 2
            
        cond3 = delta[t]*G.N / np.gcd(G.N, t)

        result = min([cond1, cond2, cond3, result])

    return result

def graph_stats(G,file_name=None):
        """Plot H and girth distribution"""
        H = G.get_H()

        girths = []
        for vn in range(G.n_cn, G.n_nodes):
            cycles = shortest_cycles(G, vn)
            if cycles:
                girths.append(min(cycles.values()))
            else:
                girths.append(-10)

        print(girths)
        
        parity_eqs = np.zeros((G.m*G.N,G.n*G.N))
        parity_eqs[:,-G.m*G.N:] = 0.2

        # plt.subplot(2,1,1)
        # plt.imshow(1 - np.maximum(H, parity_eqs), cmap="gray")
        # plt.subplot(2,1,2)
        # plt.hist(girths)
        # plt.show()
        

def make_invertable(G):
    """
    Reorder the QC_tanner_graph G such that the last m columns of parity equations in the protograph form a matrix invertible in GF(N+1) 
    <==> last n_cn equations in H invertible in GF(2)
    Solution found by deleting rows while checking invertibility
    """
    H = galois.GF2(G.get_H().astype(int))
    assert np.linalg.matrix_rank(H) == G.n_cn, "Bad matrix rank"
    
    non_invertible_H = H[:,-9*2:]

    H_indexes = np.arange(H.shape[1])

    for i in np.arange(0, H.shape[1], G.N):
        indexes = np.arange(i,(i+G.N))
        H_new = H.copy()
        H_new[:, indexes] = 0
    
        if np.linalg.matrix_rank(H_new) == H.shape[0]:
            H = H_new 
            H_indexes[indexes] = -1

    H_indexes = H_indexes[H_indexes >= 0]

    invertible_H = np.take(H, H_indexes, axis = 1)
    np.linalg.inv(invertible_H)

    reminding_columns = [i for i in range(G.n_vn) if i not in H_indexes]
    new_order = reminding_columns + list(H_indexes)
    G_invertible = G.reordered(new_order)

    H = galois.GF2(G_invertible.get_H().astype(int))
    np.linalg.inv(H[:,-G.n_cn:])
    assert np.linalg.matrix_rank(H) == G.n_cn, "Bad matrix rank"

    return G_invertible

def discretize_polynomial(polynomial,n_nodes,D):
    diff = np.sum(D)-n_nodes
    discrete_prob = 1/n_nodes
    if diff > 0:
        while diff != 0:
            tmp_arr = np.abs(polynomial/(discrete_prob) - \
                (0.5+ np.floor(polynomial/(discrete_prob))))
            idx = np.amin(tmp_arr)
            D[idx] -= 1
            diff -= 1
            polynomial[idx] -= discrete_prob
    elif diff < 0:
        while diff != 0:
            tmp_arr = np.abs(polynomial/(discrete_prob) - 
                   (0.5 + np.floor(polynomial/(discrete_prob))))
            idx = np.amin(tmp_arr)
            D[idx] += 1
            diff += 1
            polynomial[idx] += discrete_prob

    return D


def to_degree_distribution(polynomial, n_nodes):
    assert np.sum(polynomial) - 1 <1e-7
    D = np.around(polynomial*n_nodes).astype(int)
    x = np.sum(D)
    
    if x != n_nodes:
        D = discretize_polynomial(polynomial,n_nodes,D)

    C = np.zeros(n_nodes,dtype=int)
    current_idx = 0
    for degree, coeff in enumerate(D):
        if coeff:   
            C[current_idx:current_idx+coeff] = degree + 1 
            current_idx = current_idx+coeff
    return C

def vn_polynomial_repr(vn_polynomial):
    s = ""
    for degree, coeff in enumerate(vn_polynomial):
        if coeff:
            s += f"{coeff}*x^{degree+1} + "

    return s[:-2]

def strategy1(max_girth, cn_girths, G, vn_index, cn_degrees):
    #Enforce single weight circulant matrices
    max_girth = np.max(cn_girths)
    for cn in range(G.n_cn):
        if G.has_cyclical_edge_set((cn, vn_index)):
            cn_girths[cn] == -np.inf
            
    #Keep edges of maximal r-girth
    girth_survivors = np.argwhere(cn_girths == np.max(cn_girths))
    
    #Keep edges of maximal shortest cycle after adding an edge
    # cn_local_girths = np.zeros(survivors.size)
    # for i in range(survivors.size):
    #     cn_index = int(survivors[i])
    #     if G.add_cyclical_edge_set(cn_index, vn_index):
    #         cn_local_girths[i] = shortest_cycles(G, cn_index, vn_index)
    #         G.remove_cyclical_edge_set(cn_index, vn_index)
    # survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
 
    #Keep edges with maximal free degree
    current_degrees = G.get_check_degrees()
    free_degrees = cn_degrees - current_degrees
    free_degrees_survivors = np.take(free_degrees, girth_survivors)

    survivors = girth_survivors[free_degrees_survivors == np.max(free_degrees_survivors)]
    survivor = int(np.random.choice(survivors))
    
    if free_degrees[survivor] == 0:
        cond1 = (cn_degrees > cn_degrees[survivor]) 
        cond2 = (current_degrees[survivor] >= current_degrees)
        candidate_nodes = np.argwhere(np.logical_and(cond1, cond2))
        
        if len(candidate_nodes) > 0:
            swap_node = np.random.choice(candidate_nodes.flatten())
            x = G.proto_index(swap_node)
            new_val = cn_degrees[swap_node]
            cn_degrees[int(x*G.N):int((x+1)*G.N)] = cn_degrees[survivor]
            cn_degrees[survivor] = new_val
        else:
            cn_girths[survivor] = -np.inf
            max_girth = np.max(cn_girths)
            survivor = strategy1(max_girth, cn_girths, G, vn_index, cn_degrees)

    return survivor