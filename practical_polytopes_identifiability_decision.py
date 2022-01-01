import sys
from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import pypoman 
from time import time 
import itertools
from poly_iden import *
import logging
import os
import warnings
warnings.filterwarnings("ignore")

## Set the Logger
code_name = 'practical_polytopes_identifiability_decision'

####################################################
experiment_name = 'practical_dim3_6'
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

##################################################################################
poly_df1 = pd.read_pickle('poly_dfs/practical_polytopes_dim3_6.pkl')
##################################################################################

iden_list = []

for i in (range(poly_df1.shape[0])):
    try:
        V = poly_df1.iloc[i]['Vertices']
        n_vertices = V.shape[1]
        t0 = time()
        p = PolytopeIdentifiability(V)
        n_aut_order, iden = p.check_identifiability(verbose = False, return_generator_order = True)
        t1 = time()
        # iden_list.append((i, iden,n_aut_order, t1 - t0))
        if not iden:
            Gspoiler, Pspoiler = p.identifiability_spoiler_matrices()
            iden_list.append((i, iden,n_aut_order, t1 - t0, Gspoiler, Pspoiler))
        else:
            iden_list.append((i, iden, n_aut_order, t1 - t0, None, None))
            
        del V
        del p
        del n_aut_order
        del iden

    except Exception as e:
        iden_list.append((None, None, None, None, None, None))
        logger.info('The encountered error is: {}'.format(e))

    if (i % 1) == 0:
        logger.info("The number of vertices is {}".format(n_vertices))
        logger.info("{} number of iterations done".format(i))
        


poly_df1['Identifiable'] = [e[1] for e in iden_list]
poly_df1['N_Aut_Gen']  =[e[2] for e in iden_list]
poly_df1['Iden_Decision_Time'] = [e[3] for e in iden_list]
poly_df1['Gspoiler'] = [e[4] for e in iden_list]
poly_df1['Pspoiler'] = [e[5] for e in iden_list]

##############################################################################################
poly_df1.to_pickle('Results_and_Logs/Results/practical_polytopes_results_dim3_6.pkl')
##############################################################################################


logger.info('The code is done.')