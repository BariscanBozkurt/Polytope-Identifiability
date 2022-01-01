import sys
from sage.all import *
import numpy as np
from tqdm import tqdm
import pandas as pd
import pypoman 
import os
from poly_iden import *

from time import time 

import warnings
warnings.filterwarnings("ignore")

np.random.seed(2022)

if not os.path.exists('poly_dfs'):
    os.mkdir('poly_dfs')

polytopes_df_dim3_6 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
polytopes_df_dim7 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
polytopes_df_dim8 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
polytopes_df_dim9 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
polytopes_df_dim10 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
# polytopes_df_dim11 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])
# polytopes_df_dim12 = pd.DataFrame(columns = ['Dimension', 'NumberOfVertices', 'Vertices','H-repr(A-b)','Number_of_Antisparse','Number_of_Nonnegative','Number_of_Sparsity_Constraint'])

dim = 3
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim3_6 = polytopes_df_dim3_6.append(df_append, ignore_index = True)
            i += 1

dim = 4
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim3_6 = polytopes_df_dim3_6.append(df_append, ignore_index = True)
            i += 1


dim = 5
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim3_6 = polytopes_df_dim3_6.append(df_append, ignore_index = True)
            i += 1


dim = 6
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim3_6 = polytopes_df_dim3_6.append(df_append, ignore_index = True)
            i += 1


dim = 7
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim7 = polytopes_df_dim7.append(df_append, ignore_index = True)
            i += 1


dim = 8
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim8 = polytopes_df_dim8.append(df_append, ignore_index = True)
            i += 1


dim = 9
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim9 = polytopes_df_dim9.append(df_append, ignore_index = True)
            i += 1


dim = 10
print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
n_polytopes = 750
i = 0
while i < n_polytopes:
    number_of_relative_sparsity_constraint = np.random.randint(2,dim)
    number_of_antisparse = np.random.randint(dim)
    number_of_nonnegative = dim - number_of_antisparse
    try:
        (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
    except:
        pass

    if A is not None:
        NumberOfVertices = V.shape[1]
        if NumberOfVertices < 250:
            H_repr = (A,b)
            Vertices = V
            df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
                        'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
                        'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
            polytopes_df_dim10 = polytopes_df_dim10.append(df_append, ignore_index = True)
            i += 1


# dim = 11
# print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
# n_polytopes = 100
# i = 0
# while i < n_polytopes:
#     number_of_relative_sparsity_constraint = np.random.randint(1,6)
#     number_of_antisparse = np.random.randint(dim)
#     number_of_nonnegative = dim - number_of_antisparse
#     try:
#         (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
#     except:
#         pass

#     if A is not None:
#         NumberOfVertices = V.shape[1]
#         H_repr = (A,b)
#         Vertices = V
#         df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
#                      'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
#                      'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
#         polytopes_df_dim11 = polytopes_df_dim11.append(df_append, ignore_index = True)
#         i += 1


# dim = 12
# print('\n\n{} dimensional Polytope generation started!\n\n'.format(dim))
# n_polytopes = 100
# i = 0
# while i < n_polytopes:
#     number_of_relative_sparsity_constraint = np.random.randint(1,6)
#     number_of_antisparse = np.random.randint(dim)
#     number_of_nonnegative = dim - number_of_antisparse
#     try:
#         (A,b), V = generate_random_practical_polytope(dim, number_of_antisparse, number_of_nonnegative, number_of_relative_sparsity_constraint)
#     except:
#         pass

#     if A is not None:
#         NumberOfVertices = V.shape[1]
#         H_repr = (A,b)
#         Vertices = V
#         df_append = {'Dimension' : dim, 'NumberOfVertices': NumberOfVertices, 'Vertices' : V, 
#                      'H-repr(A-b)': H_repr, 'Number_of_Antisparse': number_of_antisparse, 
#                      'Number_of_Nonnegative': number_of_nonnegative, 'Number_of_Sparsity_Constraint': number_of_relative_sparsity_constraint}
#         polytopes_df_dim12 = polytopes_df_dim12.append(df_append, ignore_index = True)
#         i += 1


polytopes_df_dim3_6.to_pickle('poly_dfs/practical_polytopes_dim3_6.pkl')
polytopes_df_dim7.to_pickle('poly_dfs/practical_polytopes_dim7.pkl')
polytopes_df_dim8.to_pickle('poly_dfs/practical_polytopes_dim8.pkl')
polytopes_df_dim9.to_pickle('poly_dfs/practical_polytopes_dim9.pkl')
polytopes_df_dim10.to_pickle('poly_dfs/practical_polytopes_dim10.pkl')
# polytopes_df_dim11.to_pickle('poly_dfs/practical_polytopes_dim11.pkl')
# polytopes_df_dim12.to_pickle('poly_dfs/practical_polytopes_dim12.pkl')