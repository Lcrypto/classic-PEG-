"""
Contains classes Tanner_graph and QC_tanner_graph.
"""

import matplotlib.pyplot as plt 
import numpy as np

class Tanner_graph:
    """
    Sparse implementation of a tanner graph with dimensions m x n
    
    Check nodes are stored in self.nodes[0:n_cn] and variable nodes are stored
    in self.nodes[n_cn:n_nodes], but external methods should not care about this, meaning
    variable nodes are indexed from 0,..., n_vn-1 externally.
    """

    def __init__(self, m, n):
        """Creates a graph with m*N check nodes, n*N variable nodes and no edges. If swc is set to true the Single Weight Constraint is enforced, 
        meaning each circulant matrix has weight 1 or 0"""
        assert m > 0 and n > 0, "m and n must be positive integers."
    
        self.n_nodes = int((m+n))
        self.n_cn = m
        self.n_vn = n 
        self.m = m
        self.n = n
        self.nodes = [set() for i in range(self.n_nodes)]

    def __repr__(self):
        return f"{self.n_cn} CNs, {self.n_vn} VNs, {sum([len(node) for node in self.nodes])//2} edges"

    def assert_edge(self, edge):
        assert 0 <= edge[0] < self.n_cn, f"Edge non existent, check node index {edge[0]} out of bounds."
        assert 0 <= edge[1] < self.n_nodes, f"Edge non existent, variable node index {edge[1]} out of bounds."

    def assert_vn_node(self, vn):
        assert 0 <= vn < self.n_vn, "Variable node index out of bounds."

    def assert_cn_node(self, cn):
        assert 0 <= cn < self.n_cn, "Check node index out of bounds."

    def add_edges(self, edges):
        """Add edges defined as pairs of check nodes and variable nodes to the graph"""
        for edge in edges:
            self.assert_edge(edge)
            self.nodes[edge[0]].add(edge[1] + self.n_cn) 
            self.nodes[edge[1] + self.n_cn].add(edge[0])

    def remove_edges(self, edges):
        """Remove edges defined as pairs of check nodes and variable nodes from the graph"""
        for edge in edges:
            assert self.has_edge(edge), "Cannot remove non existent edge"
            self.nodes[edge[0]].remove(edge[1] + self.n_cn) 
            self.nodes[edge[1] + self.n_cn].remove(edge[0])

    def get_adjecent(self, node):
        """Returns all adjecent nodes of node index node"""
        return self.nodes[int(node)]

    def get_adjecent_cn(self, cn):
        """Returns adjecent nodes of check node index cn"""
        self.assert_cn_node(cn)
        return self.get_adjecent(cn)

    def get_adjecent_vn(self, vn):
        """Returns adjecent nodes of variable node index vn"""
        self.assert_vn_node(vn)
        return self.get_adjecent(vn + self.n_cn)

    def has_edge(self, edge):
        """Returns true if the graph contains the edge (ci, vi)"""
        self.assert_edge(edge)
        return edge[1] + self.n_cn in self.nodes[edge[0]]

    def get_check_degrees(self) -> list:
        """Returns the degree of all check nodes of the graph"""
        return np.array([len(self.nodes[i]) for i in range(self.n_cn)])

    def get_var_degrees(self) -> list:
        """Returns the degree of all check nodes of the graph"""
        return np.array([len(self.nodes[i]) for i in range(self.n_cn, self.n_nodes)])

    def get_H(self):
        """Generates a dense representation of the graph"""
        H = np.zeros((self.n_cn, self.n_vn))
        
        for i, nodes in enumerate(self.nodes[0:self.n_cn]):
            for j in nodes:
                H[i,j-self.n_cn] = 1

        return H

    def plot(self):
        """Graphical representation of the tanner graph"""
        width = 500
        height = 250
        border = 100
        vn_coords = np.stack((np.linspace(0, width, self.n_vn), np.full(self.n_vn, height)))
        cn_coords = np.stack((np.linspace(0, width, self.n_cn), np.full(self.n_cn, 0)))
        
        #Init figure
        plt.figure()
        plt.xlim(-border, width+border)
        plt.ylim(-border, height + border)
        
        #Plot nodes
        plt.scatter(vn_coords[0,:],vn_coords[1,:], s = 40, c = "black", marker = "o")
        plt.scatter(cn_coords[0,:],cn_coords[1,:], s = 40, c = "black", marker = "s")
        
        #Plot edges
        for cn, vns in enumerate(self.nodes[:self.n_cn]):
            for vn in vns:
                vn_x = vn_coords[0,vn-self.n_cn]
                vn_y = vn_coords[1,vn-self.n_cn]
                cn_x = cn_coords[0, cn]    
                cn_y = cn_coords[1, cn]    
                plt.plot([vn_x, cn_x], [vn_y, cn_y], c = "black")

        plt.show()

class QC_tanner_graph(Tanner_graph):
    """
    Extends the tanner graph to work with QC-codes.
    """

    def __init__(self, m, n, N):
        """
        Sparse implementation of a qc tanner graph with protograph of dimensions m x n, 
        a scaling factor N, and single weight constraint on circulant matrixes.
        
        Check nodes are stored in self.nodes[0:n_cn] and variable nodes are stored
        in self.nodes[n_cn:n_nodes], but external methods should not care about this, meaning
        variable nodes are indexed from 0,..., n_vn-1 externally.
        """
        assert m > 0 and n > 0 and N > 0, "m, n and N must be positive integers."
        

        self.n_nodes = int((m+n)*N)
        self.n_cn = int(m*N)
        self.n_vn = int(n*N) 
        self.N = int(N)  
        self.m = m
        self.n = n

        self.nodes = [set() for i in range(self.n_nodes)]
        self.proto = np.full((m, n), -1)

    def save(self, filename):
        """
        Saves the protograph as a .qc file, as defined here: https://aff3ct.readthedocs.io/en/latest/user/simulation/parameters/codec/ldpc/decoder.html
        Appends metadata at end
        """
        data = f"{self.n} {self.m} {self.N}\n\n"
        for row in self.proto:
            for char in row:
                data = data + str(char) + " "
            data += "\n"
            f = open(filename, "w")

        f.write(data)
        f.close() 

    @staticmethod
    def read(filename):
        """
        Creates a new graph from a .qc file
        """
        meta_data = np.loadtxt(filename, max_rows=1)
        n = int(meta_data[0])
        m = int(meta_data[1])
        N = int(meta_data[2])

        proto = np.loadtxt(filename, skiprows=2)

        G = QC_tanner_graph(m, n, N)

        for i in range(m):
            for j in range(n):
                if not proto[i,j] == -1:
                    cn = i*N + proto[i,j]
                    vn = j*N

                    G.add_cyclical_edge_set(cn, vn)
        G.proto = proto

        return G

    def proto_index(self, node):
        """Returns the index of node in the protograph"""
        return np.floor(node/self.N)

    def proto_value(self, node):
        """Returns the shift of the 1-weighted circulant matrix containing node"""
        return np.mod(node, self.N)

    def has_cyclical_edge_set(self, edge):
        """Checks if the graph has a cyclical edge set for the given edge"""
        i = int(self.proto_index(edge[0]))
        j = int(self.proto_index(edge[1]))
        return not self.proto[i,j] == -1

    def add_cyclical_edge_set(self, cn_index, vn_index):
        """Adds a cyclical edge set pi(ci, vi, N) to the graph, returning true on success and false otherwise"""
        self.assert_edge((cn_index, vn_index))
      
        i = int(self.proto_index(cn_index))
        j = int(self.proto_index(vn_index))
        if self.proto[i,j] == -1:
            self.proto[i,j] = np.mod(self.proto_value(vn_index) - self.proto_value(cn_index), self.N)

            t = np.arange(self.N)
            check_nodes = self.proto_index(cn_index)*self.N + self.proto_value(cn_index + t)
            variable_nodes = self.proto_index(vn_index)*self.N + self.proto_value(vn_index + t)
            self.add_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))
            return True 
        return False

    def remove_cyclical_edge_set(self, cn_index, vn_index):
        """Removes a cyclical edge set pi(ci, vi, N) from the graph, returning true on success and false otherwise"""
        self.assert_edge((cn_index, vn_index))

        i = int(self.proto_index(cn_index))
        j = int(self.proto_index(vn_index))
        if not self.proto[i,j] == -1:
            self.proto[i,j] = -1

            t = np.arange(self.N)
            check_nodes = i*self.N + self.proto_value(cn_index + t)
            variable_nodes = j*self.N + self.proto_value(vn_index + t)
            self.remove_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))
            return True 

        return False

    def reordered(self, index_list):
        """Returns a new graph with an equivalent tanner graph and reordered variable nodes according to the index list"""
        G = QC_tanner_graph(self.m, self.n, self.N)
        index_list = np.array(index_list)

        for i, row in enumerate(self.proto):
            for j, shift in enumerate(row):
                if not shift == -1:
                    vn_index = int(np.argwhere(index_list == j*self.N) + shift)
                    cn_index = int(i * self.N)
                    G.add_cyclical_edge_set(cn_index, vn_index)
                
        return G

