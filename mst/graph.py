import numpy as np
import heapq
from typing import Union




class Graph:

    def __init__(self, adjacency_mat: Union[np.ndarray, str]):
        """
    
        Unlike the BFS assignment, this Graph class takes an adjacency matrix as input. `adjacency_mat` 
        can either be a 2D numpy array of floats or a path to a CSV file containing a 2D numpy array of floats.

        In this project, we will assume `adjacency_mat` corresponds to the adjacency matrix of an undirected graph.
    
        """
        if type(adjacency_mat) == str:
            self.adj_mat = self._load_adjacency_matrix_from_csv(adjacency_mat)
        elif type(adjacency_mat) == np.ndarray:
            self.adj_mat = adjacency_mat
        else: 
            raise TypeError('Input must be a valid path or an adjacency matrix')
        self.mst = None

    def _load_adjacency_matrix_from_csv(self, path: str) -> np.ndarray:
        with open(path) as f:
            return np.loadtxt(f, delimiter=',')

    def construct_mst(self):
        """
    
        TODO: Given `self.adj_mat`, the adjacency matrix of a connected undirected graph, implement Prim's 
        algorithm to construct an adjacency matrix encoding the minimum spanning tree of `self.adj_mat`. 
            
        `self.adj_mat` is a 2D numpy array of floats. Note that because we assume our input graph is
        undirected, `self.adj_mat` is symmetric. Row i and column j represents the edge weight between
        vertex i and vertex j. An edge weight of zero indicates that no edge exists. 
        
        This function does not return anything. Instead, store the adjacency matrix representation
        of the minimum spanning tree of `self.adj_mat` in `self.mst`. We highly encourage the
        use of priority queues in your implementation. Refer to the heapq module, particularly the 
        `heapify`, `heappop`, and `heappush` functions.

        """
        self.mst = None

        # Initializing visited node list
        visited = set()
        # Initializing list for holding edge weights and connecting nodes
        neighbor_list = []

        
        # Setting adj matrix to a variable
        mat = self.adj_mat

        # Checking for numpy matrix
        if not isinstance(mat, np.ndarray):
            raise Exception('Adjacency matrix must be a numpy array!')
        
        # Checking for 2D matrix
        if mat.ndim != 2:
            raise Exception('Matrix is not 2D')

        # Checking if adj mat is empty
        if np.all(mat == 0):
            raise ValueError('Adjacency matrix has no edges')
        
        # Initializing self.mstß
        self.mst = np.zeros_like(mat)
        
        # Extracting number of rows OR number of nodes in the graph
        row_len, _ = mat.shape

        # Set a random start node
        start_node = 0
        
        # Adding start node to the visited set
        visited.add(start_node)

        
        # While our loop has not traversed all nodes in the graph
        for start_node in range(row_len):
            if start_node not in visited:
                visited.add(start_node)
            
            # Loop to find neighbors of start 
            for i in range(row_len):
                edge = mat[start_node, i] 
                if edge != 0:
                    heapq.heappush(neighbor_list, (edge, (start_node, i)))

            while neighbor_list:
                weight, (parent, child) = heapq.heappop(neighbor_list)

                if child not in visited:
                    visited.add(child)
                    start_node = child
                    self.mst[parent, child] = self.mst[child, parent] = weight
                    

                    for i in range(row_len):
                        if mat[child, i] != 0 and i not in visited:
                            heapq.heappush(neighbor_list, (mat[child, i], (child, i)))
        
            # Check if all nodes have been visited
            if len(visited) == row_len:
                break  # All nodes visited, exit the loop
        












        
        


