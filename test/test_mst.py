import pytest
import numpy as np
from mst import Graph
from sklearn.metrics import pairwise_distances


def check_mst(adj_mat: np.ndarray, 
              mst: np.ndarray, 
              expected_weight: int, 
              allowed_error: float = 0.0001):
    """
    
    Helper function to check the correctness of the adjacency matrix encoding an MST.
    Note that because the MST of a graph is not guaranteed to be unique, we cannot 
    simply check for equality against a known MST of a graph. 

    Arguments:
        adj_mat: adjacency matrix of full graph
        mst: adjacency matrix of proposed minimum spanning tree
        expected_weight: weight of the minimum spanning tree of the full graph
        allowed_error: allowed difference between proposed MST weight and `expected_weight`

    TODO: Add additional assertions to ensure the correctness of your MST implementation. For
    example, how many edges should a minimum spanning tree have? Are minimum spanning trees
    always connected? What else can you think of?

    """

    def approx_equal(a, b):
        return abs(a - b) < allowed_error

    total = 0
    for i in range(mst.shape[0]):
        for j in range(i+1):
            total += mst[i, j]
    
    assert approx_equal(total, expected_weight), 'Proposed MST has incorrect expected weight'
    
    
    # Additional assertion to check for correct number of edges
    num_edges = (np.count_nonzero(mst)) / 2

    num_nodes = mst.shape[0]

    # Assert that number of edges is the expected number
    assert num_edges == num_nodes - 1
    
    # Assert that number of nodes in initial matrix matches the MST matrix
    assert num_nodes == adj_mat.shape[0]
    
    # Assert number of edges in MST is less than or equal to number of edges in adj_mat
    assert num_edges <= np.count_nonzero(adj_mat)
    
    

def test_mst_small():
    """
    
    Unit test for the construction of a minimum spanning tree on a small graph.
    
    """
    file_path = './data/small.csv'
    g = Graph(file_path)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 8)

    
def test_mst_single_cell_data():
    """
    
    Unit test for the construction of a minimum spanning tree using single cell
    data, taken from the Slingshot R package.

    https://bioconductor.org/packages/release/bioc/html/slingshot.html

    """
    file_path = './data/slingshot_example.txt'
    coords = np.loadtxt(file_path) # load coordinates of single cells in low-dimensional subspace
    dist_mat = pairwise_distances(coords) # compute pairwise distances to form graph
    g = Graph(dist_mat)
    g.construct_mst()
    check_mst(g.adj_mat, g.mst, 57.263561605571695)
    
    
def test_mst_student():
    """
    
    TODO: Write at least one unit test for MST construction.
    
    """

    # Disconnected input adj matrix
    mat = np.array([[0, 1, 1, 1, 0],
       [1, 0, 1, 1, 0],
       [1, 1, 0, 1, 0],
       [1, 1, 1, 0, 0],
       [0, 0, 0, 0, 0]])

    g = Graph(mat)
    g.construct_mst()
    num_edges = (np.count_nonzero(g.mst)) / 2
    num_nodes = g.mst.shape[0]

    # Assert that number of edges WILL NOT be 1 less than number of nodes, because input graph is disconnected!
    assert num_edges != num_nodes - 1

    #Check if weight of MST is lower that that of adj mat
    assert np.sum(g.adj_mat) > np.sum(g.mst)

    
    # Asserting ValueError for empty matrix
    empty_graph = np.zeros((5,5))
    g = Graph(empty_graph)
    with pytest.raises(ValueError, match='Adjacency matrix has no edges'):
        g.construct_mst()
