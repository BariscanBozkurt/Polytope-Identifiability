import sys
import os
import numpy as np
from tqdm import tqdm
import pandas as pd
from time import time 
import matplotlib.pyplot as plt
from sympy.utilities.iterables import multiset_permutations
import logging
import os
import warnings
warnings.filterwarnings("ignore")


## Set the Logger
code_name = 'practical_polytopes_identifiability_decision_brute_force'

####################################################
experiment_name = 'practical_brute_force'
#################################################

if not os.path.exists('Results_and_Logs'):
    os.mkdir('Results_and_Logs')

if not os.path.exists('Results_and_Logs/Results'):
    os.mkdir('Results_and_Logs/Results')

if not os.path.exists(os.path.join('Results_and_Logs', experiment_name)):
    os.mkdir(os.path.join('Results_and_Logs', experiment_name))
    
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s | | %(levelname)s | | %(message)s')

logger_file_name = os.path.join(experiment_name, code_name)
# logger_file_name = input('Enter the Logger File Name: ')
# logger_file_name = 'Trial'
logger_file_name = os.path.join('Results_and_Logs', logger_file_name)
file_handler = logging.FileHandler(logger_file_name,'w')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

logger.info('Code started. The Experiment Name is {}.'.format(experiment_name))



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
    for i,indxp in (enumerate(multiset_permutations(indxs))):
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

results_path = './poly_dfs'

df_paths = os.listdir(results_path)

poly_df = pd.DataFrame()
for i,df_path in enumerate(df_paths):
    df_append = pd.read_pickle(os.path.join(results_path, df_path))
    poly_df = poly_df.append(df_append, ignore_index = False)

poly_df_for_bf = poly_df.loc[poly_df['NumberOfVertices'] <= 8].copy()
poly_df_for_bf = poly_df_for_bf.append(poly_df.loc[(poly_df['NumberOfVertices'] == 9)].sample(150), ignore_index = True)
poly_df_for_bf = poly_df_for_bf.append(poly_df.loc[(poly_df['NumberOfVertices'] == 10)].sample(50), ignore_index = True)
poly_df_for_bf = poly_df_for_bf.append(poly_df.loc[(poly_df['NumberOfVertices'] == 11)].sample(40), ignore_index = True)

del poly_df

iden_time = []

for i in (range(poly_df_for_bf.shape[0])):
    logger.info('Current Iteration Number : {}'.format(i))
    try:
        V = poly_df_for_bf.iloc[i]['Vertices']
        logger.info('Number of Vertices : {}'.format(V.shape[1]))
        n_vertices = V.shape[1]
        t0 = time()
        iden = CheckPolytopeBruteForce(V)
        t1 = time()
        iden_time.append(t1 - t0)
    except Exception as e:
        logger.info('Exception at {}'.format(i))
        logger.info('Exception is {}'.format(e))

poly_df_for_bf['Iden_Decision_Time_BF'] = iden_time

poly_df_for_bf.to_pickle('Results_and_Logs/Results/practical_polytopes_results_brute_force.pkl')