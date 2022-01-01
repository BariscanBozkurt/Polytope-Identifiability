"""
Title: poly_iden.py

Identifiable Polytope Characterization

Code Writer: Bariscan Bozkurt (Koç University - EEE & Mathematics)

Date: 08.03.2021
"""


# import sys
from sage.all import *
import numpy as np
from tqdm import tqdm
import pypoman 
import itertools
from sympy.utilities.iterables import multiset_permutations
from IPython.display import display, Latex, Math

class PolytopeIdentifiability():
    """
    
    A python class written based on SageMath software to check a given polytope is identifiable or not.
    The identifiability definition is given in the following references. Let P be a polytope with vertex
    representation, s.t., v_1, ..., v_m, be the vertices of the polytope. Let V = [v_1 ... v_m] in Re^{m \\times n}
    A polytope is defined to be identifiable [1,2] if and only if for any permutation matrix R in Re^{m \\times m}
    the solution to the equation GV = VR holds only for a G in Re^{n \\times n} which is signed permutation 
    (multiplication of a diagonal matrix with +-1 with a permutation). For further reference to signed permutation,
    check [6]. In this algorithm, we will use the edge colored graph automorphism group to find the linear automorphism
    group [3,4] of the polytope which will enable us to check whether a given polytope is identifiable or not in a faster way.
    The code is written based on the SageMath built in functions (graph-automorphism & restricted automorphism, etc.). In order
    to use this code, you need the SageMath installed, check [7,8]. It is recommended to use a virtual environment or a conda 
    environment with SageMath in it. Check the requirements.txt file come with this project.

    Parameters:
    ===============================================================================
    V         :                                     An (m by n) matrix whose colums are the vertices of a polytope

    Methods   :
    ================================================================================
    check_sign_perm(self,R)

    lin_group_generators(self)

    check_identifiability(self)

    identifiability_spoiler_matrices(self)


    REFERENCES

    [1]Gokcan Tatli and Alper T. Erdogan,  “Polytopic matrix factorization:  Determinant maximization based criterionand identifiability,”IEEE Transactions on Signal Processing,2020

    [2]Gokcan Tatli and Alper T. Erdogan, ”Generalized polytopic matrix factorization”,IEEE ICASSP, 2021

    [3]D Bremner,  MD Sikiric,  A Sch ̈urmann,  ”Polyhedralrepresentation conversion up to symmetries”,Polyhedral Com-putation, CRM Proceedings & Lecture Notes, 45-72, 2009

    [4]D Bremner,  MD Sikiri ́c,  DV Pasechnik,  T Rehn,  ASch ̈urmann, ”Computing symmetry groups of polyhedra”,LMSJournal  of  computation  and  mathematics  17  (1),  565-581,2014
    
    [5]B.D.  McKay,  A.  Piperno,  ”Practical  Graph  Isomor-phism, II”, J. Symbolic Computation, 2013, Preprint versionat https://arxiv.org/pdf/1301.1493.pdf

    [6]https://en.wikipedia.org/wiki/Generalized_permutation_matrix

    [7]https://www.sagemath.org/

    [8]https://doc.sagemath.org/html/en/installation/conda.html
    """
    def __init__(self,V):
        self.V = V
        self.m = V.T.shape[0]
        self.n = V.T.shape[1]
        self.lin_group = None
        self.G_spoiler = None
        self.Perm_spoiler = None
        
    def check_sign_perm(self,perm_mat):
        """Check whether a given matrix is sgined permutation matrix or not."""
        n = perm_mat.shape[0]
        if (np.linalg.norm((np.abs(perm_mat)>1e-6).sum(axis = 0) - np.ones((1,n))) < 1e-6):
            return True
        else: 
            return False

    def poly_linear_automorphism_group(self, output = "matrix"):
            V_ = self.V
            Poly = Polyhedron(vertices = V_.T)

            index0 = 0
            V = [v.homogeneous_vector()[:-1]  for v in Poly.Vrepresentation()]
            Q = V_ @ V_.T
            Qplus = np.linalg.inv(Q)
            C = np.array(V_.T @ Qplus @ V_, dtype = np.float32)
            C = C * (np.abs(C) > 1e-7)

            G = Graph()
            for i in range(len(V)):
                for j in range(i+1, len(V)):
                    G.add_edge(index0+i, index0+j, C[i,j])
            permgroup = G.automorphism_group(edge_labels=True)

            if output == "permutation":
                return permgroup
            elif output == "matrix":
                permgroup = permgroup.gens()

            # Compute V+ = Vt Q+ as list of row vectors
            # Vplus = list(matrix(V) * Qplus)  # matrix(V) is Vt
            Vplus = matrix(V_.T @ Qplus)
            # Vplus = V_.T @ Qplus
            # Compute W = 1 - V V+
            # W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

            # Convert the permutation group to a matrix group.
            # If P is a permutation, then we return the matrix
            # B = (V P V+) + W.
            #
            # If output == "matrix", we loop over the generators of the group.
            # Otherwise, we loop over all elements.
            matrices = []
            for perm in permgroup:
                A = sum(V[perm(i)].column() * Vplus[i].row() for i in range(len(V)))
                matrices.append(A)

            for mat in matrices:
                mat.set_immutable()
                
            return [np.array(element) for element in matrices]

    def lin_group_generators_old(self):
        lin_group = []
        
        V = self.V
        m = self.m
        n = self.n

        p = Polyhedron(vertices = V.T)
        
        if np.allclose(V.mean(axis = 1), np.zeros(n)):
            affine_grp = p.restricted_automorphism_group(output = 'matrix')
            lin_gens = affine_grp.gens()
            for i in range(len(lin_gens)):
                lin_group.append(np.array(matrix(lin_gens[i]))[0:n,0:n])
        else:
            affine_grp = p.restricted_automorphism_group(output = 'matrixlist')
            order = len(affine_grp)
            for i in range(order):
                G = np.array(matrix(affine_grp[i]))
                if np.allclose(G[:-1,n], np.zeros(n)):
                    lin_group.append(G[0:n,0:n])

        lin_group = np.array(lin_group)
        self.lin_group = lin_group
        
        return lin_group

    def check_identifiability_old(self, verbose = True, printing = False):
        lin_group = self.lin_group_generators_old()
        order = lin_group.shape[0]

        check_val = 0
        if verbose:
            for i in tqdm(range(order)):
                G = lin_group[i]
                if not self.check_sign_perm(G):
                    check_val += 1
                    self.G_spoiler = G
                    break
        else:
            for i in (range(order)):
                G = lin_group[i]
                if not self.check_sign_perm(G):
                    check_val += 1
                    self.G_spoiler = G
                    break        

        identifiability = (not (check_val > 0))
        self.identifiability = identifiability
        if printing:
            if identifiability:
                print('Identifiable')
            else:
                print('Not identifiable')
        return identifiability

    def lin_group_generators(self):
        lin_group = self.poly_linear_automorphism_group()
        self.lin_group = lin_group
        return lin_group

    def check_identifiability(self, verbose = True, printing = False, return_generator_order = False):
        lin_group = self.lin_group_generators()
        order = len(lin_group)

        check_val = 0
        if verbose:
            for i in tqdm(range(order)):

                G = lin_group[i]
                if np.linalg.matrix_rank(G) == G.shape[0]:
                    if not self.check_sign_perm(G):
                        check_val += 1
                        self.G_spoiler = G
                        print(G)
                        break
        else:
            for i in (range(order)):
                G = lin_group[i]
                if np.linalg.matrix_rank(G) == G.shape[0]:
                    if not self.check_sign_perm(G):
                        check_val += 1
                        self.G_spoiler = G
                        break        

        identifiability = (not (check_val > 0))
        self.identifiability = identifiability
        if printing:
            if identifiability:
                print('Identifiable')
            else:
                print('Not identifiable')
        if return_generator_order:
            return order, identifiability
        else:
            return identifiability
    
    def identifiability_spoiler_matrices(self):
        if self.G_spoiler is not None:
            G_spoiler = self.G_spoiler
            V = self.V
            V_perm = G_spoiler @ V
            m = self.m
            perm = []
            for i in range(m):
    #             perm.append(np.where(np.all((V_perm - V[:,i][np.newaxis].T) == 0, axis = 0))[0][0])
                perm.append(np.where(np.all(np.abs(V_perm - V[:,i][np.newaxis].T) < 1e-6, axis = 0))[0][0])
            
            P_spoiler = np.eye(m)[perm]
            return G_spoiler, P_spoiler
        else:
            if self.identifiability is not None:
                print('Polytope is Identifiable. There exist no identifiability spoiler matrix')
                return None, None
            else:
                print('You must run check_identifiability function in this class before using this function!!!')
                return None, None


def check_sign_perm(perm_mat):
    """Check whether a given matrix is sgined permutation matrix or not."""
    n = perm_mat.shape[0]
    if (np.linalg.norm((np.abs(perm_mat)>1e-6).sum(axis = 0) - np.ones((1,n))) < 1e-6):
        return True
    else: 
        return False


def CheckPolytopeBruteForce(V):
    """
    Brute force algorithm to check whether a given polytope is identifiable or not.
    (Check the above references for identifiability definition). In this function,
    all possible permutation matrices (m! possibilities) are tried to check the identifiability
    of a polytope.
    """
    n = V.shape[0]
    m = V.shape[1]
    checkval = 0
    indxs = np.arange(m)
#     print('{:d} Vertices:'.format(m))
    Vp = np.linalg.pinv(V)
    for i,indxp in tqdm(enumerate(multiset_permutations(indxs))):
        Vindxp = V[:, indxp]
        G = Vindxp@Vp
        if (np.linalg.norm(G@V-Vindxp) < 1e-6):
            GG = G*(abs(G) > 1e-3)
            # if (np.linalg.norm(np.abs(np.sign(GG)).sum(axis = 0) - np.ones((1,n))) > 1e-6):
            #     # print('Not identifiable, the index is {}'.format(i))
            #     checkval += 1
            if not check_sign_perm(GG):
                checkval +=1
#                 print('Not Identifiable!')
                break
    
    if (checkval > 0):
#         print('Not Identifiable!')
        return False
    else:
#         print('Identifiable!')
        return True
    

def generate_random_practical_polytope(dim = 3, number_of_antisparse = None, number_of_nonnegative = None, number_of_relative_sparsity_constraint = 2, return_vertices = True):
    try:
        A = []
        b = []
        if number_of_antisparse is None:
            number_of_antisparse = np.random.randint(dim)
        if number_of_nonnegative is None:
            number_of_nonnegative = dim - number_of_antisparse

        if (number_of_nonnegative + number_of_antisparse) != dim:
            print('Number of antisparse and nonnegative dimensions should add up to the Polytope Dimension (which is input dim).')
            print('Adjust your inputs accordinly')
            return (None, None), None

        dim_indices = [i for i in range(dim)]
        antisparse_dims = np.random.choice(dim_indices, number_of_antisparse, replace = False)
        nonnegative_dims = np.random.choice(dim_indices, number_of_nonnegative, replace = False)

        for j in antisparse_dims:
            row1 = [0 for _ in range(dim)]
            row2 = row1.copy()
            row1[j] = 1
            A.append(row1)
            b.append(1)
            row2[j] = -1
            A.append(row2)
            b.append(1)

        for j in nonnegative_dims:
            row1 = [0 for _ in range(dim)]
            row2 = row1.copy()
            row1[j] = 1
            A.append(row1)
            b.append(1)
            row2[j] = -1
            A.append(row2)
            b.append(0)

        for _ in range(number_of_relative_sparsity_constraint):

            relative_sparse_dims = np.random.choice(dim_indices, np.random.randint(2,dim), replace = False)
            row = np.zeros(dim)
            pm_one = [[1,-1] for _ in range(relative_sparse_dims.shape[0])]
            for i in itertools.product(*pm_one):
                row_copy = row.copy()
                row_copy[relative_sparse_dims] = i
                A.append(list(row_copy))
                b.append(1)

        A = np.array(A)
        b = np.array(b)
        if return_vertices:
            vertices = pypoman.compute_polytope_vertices(A, b)
            V = np.array([list(v) for v in vertices]).T

            return (A,b), V
        else:
            return (A,b)
    except:
        print('Random polytope generation failed. Try the function again')
        if return_vertices:
            return (None, None), None
        else:
            return (None, None)
        
def polytope_linear_automorphism_group_v2(V_, output = "matrix"):
    P = Polyhedron(vertices = V_.T)
    index0 = 0
    V = [v.homogeneous_vector()[:-1]  for v in P.Vrepresentation()]
    # Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()
    Q = V_ @ V_.T
    Qplus = np.linalg.inv(Q)
    C = np.array(V_.T @ Qplus @ V_, dtype = np.float32)
    C = C * (np.abs(C) > 1e-7)

    G = Graph()
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            G.add_edge(index0+i, index0+j, C[i,j])
    permgroup = G.automorphism_group(edge_labels=True)

    if output == "permutation":
        return permgroup
    elif output == "matrix":
        permgroup = permgroup.gens()

    # Compute V+ = Vt Q+ as list of row vectors
    # Vplus = list(matrix(V) * Qplus)  # matrix(V) is Vt
    Vplus = matrix(V_.T @ Qplus)

    # Compute W = 1 - V V+
    W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

    # Convert the permutation group to a matrix group.
    # If P is a permutation, then we return the matrix
    # B = (V P V+) + W.
    #
    # If output == "matrix", we loop over the generators of the group.
    # Otherwise, we loop over all elements.
    matrices = []
    for perm in permgroup:
        A = sum(V[perm(i)].column() * Vplus[i].row() for i in range(len(V)))
        matrices.append(A + W)

    for mat in matrices:
        mat.set_immutable()
        
    return MatrixGroup(matrices)

def polytope_linear_automorphism_group_v3(V_, output = "matrix"):
    P = Polyhedron(vertices = V_.T)
    index0 = 0
    V = [v.homogeneous_vector()[:-1]  for v in P.Vrepresentation()]
    # Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()
    Q = V_ @ V_.T
    Qplus = np.linalg.inv(Q)
    C = np.array(V_.T @ Qplus @ V_, dtype = np.float32)
    C = C * (np.abs(C) > 1e-7)

    G = Graph()
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            G.add_edge(index0+i, index0+j, C[i,j])
    permgroup = G.automorphism_group(edge_labels=True)

    if output == "permutation":
        return permgroup
    elif output == "matrix":
        permgroup = permgroup.gens()

    # Compute V+ = Vt Q+ as list of row vectors
    # Vplus = list(matrix(V) * Qplus)  # matrix(V) is Vt
    Vplus = matrix(V_.T @ Qplus)

    # Compute W = 1 - V V+
    W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

    # Convert the permutation group to a matrix group.
    # If P is a permutation, then we return the matrix
    # B = (V P V+) + W.
    #
    # If output == "matrix", we loop over the generators of the group.
    # Otherwise, we loop over all elements.
    matrices = []
    for perm in permgroup:
        A = sum(V[perm(i)].column() * Vplus[i].row() for i in range(len(V)))
        #matrices.append(A + W)
        matrices.append(A + W)

    for mat in matrices:
        mat.set_immutable()
        
    return [np.array(element) for element in matrices]

def display_matrix(array):
    data = ''
    for line in array:
        if len(line) == 1:
            data += ' %.3f &' % line + r' \\\n'
            continue
        for element in line:
            data += ' %.3f &' % element
        data += r' \\' + '\n'
    display(Math('\\begin{bmatrix} \n%s\end{bmatrix}' % data))