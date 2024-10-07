import pandas as pd
import numpy as np
import json
import statistics
import cvxpy
import re
import sys
import os

# be sure to run run.py before running this file
data_organism = "Escherichia_coli"
prokode_dir = "C:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE"
network_loc = prokode_dir + "/src/network.json"
data_dir = prokode_dir + "/beta_training/GEO_expression_data/" + data_organism

# organism specs
genome_len = 4500000
init_cell_volume = 1.8 # μm^3 (Bionumbers)

sys.path.append('..')

# global jazz
num_data = 10
sensor_normal_dist = 10
basal_rate = 3
decay_rate = np.log(2)/300
Np = 6000
Kd_p = 0.1
Nns = 4600000

def function_to_test(data_t0, protein_data_t0, dt):
     # rate_of_protein_creation = some_function_dependent_on_num_ribo()
     # overall_protein_half_life = some_constant_func_independent_of_other_stuff()

     # overall_rate_of_protein_change = rate_of_protein_creation * data_t0  - overall_protein_half_life * protein_data_t0

     # rate_of_mRNA_creation =
     # overall_mRNA_half_life = some_constant_func_independent_of_other_stuff()

     # overall_rate_of_mRNA_change = rate_of_mRNA_creation - overall_mRNA_half_life * data_t0

     # protein_data_t1 = protein_data_t0 + overall_rate_of_protein_change * dt
     # data_t1 = data_t0 + overall_rate_of_mRNA_change * dt
     return data_t1

def protein_amnts_from_mRNA_amnts(mRNA_amnts):
     protein_amnts = mRNA_amnts
     return protein_amnts

from scipy.optimize import curve_fit
from src.maths import *

def function_for_timepoint(data_t0, data_t1, dt, gene_key, tf_key):
     # MOST CRITICAL PART TO DETERMINE!!!
     protein_data_t0 = protein_amnts_from_mRNA_amnts(data_t0)

     for i in range(len(gene_key)):
          gene = gene_key[i]
          gene_mRNA_t0 = data_t0[i]

          gene_info_dict = search_network_json(network_loc, gene)
          regulators_dict = gene_info_dict["regulators"]
          Kd_rnap_gene = regulators_dict["polymerase"]["Kd"]

          overall_mRNA_change_rate = (data_t1[i] - data_t0[i]) / dt
          coefficient_arr, beta_all = beta_from_overall_mRNA(gene, gene_mRNA_t0, overall_mRNA_change_rate, protein_data_t0, regulators_dict, gene_key, tf_key, genome_length, cell_volume, Kd_rnap_gene)

     return coefficient_arr, beta_all

def fit_function(coefficient_matrix, beta_all_arr):
     # clean matrices
     empty_columns = []
     for i in zip(*coefficient_matrix):
          column = zip(*coefficient_matrix)[i]
          if column == [0] * len(coefficient_matrix):
               empty_columns.append(i)
          else:
               continue

     # curve fitting parameters
     func_to_fit = beta_function
     n_features = len(coefficient_matrix[0])
     p = [1] * n_features
     beta_arr, covar_uncertainty = curve_fit(func_to_fit, coefficient_matrix.T, beta_all_arr, p0 = p)

     for index in empty_columns:
          beta_arr[index] = None

     return beta_arr, covar_uncertainty



for data_file in ["C:/Users/cryst/LOFScreening/archive/PROKODE/ProkODE/beta_training/GSE10159_results.csv"]: # os.listdir(data_dir):
     # get data
     data_df = pd.read_csv(data_file, index_col=0)
     gene_key = list(data_df.index) # may have to subtracted header
     tf_key = [ tf for val in json.load(open(network_loc, 'r')).values() for tf in val["regulators"].keys() ]
     tf_key = list(set(tf_key)) # filter repeats

     # prep stuffs
     # beta_all (1, n_samples) coefficient_matrix (n_samples, n_coefficients)
     beta_all_arr, coefficient_matrix = ([], [])
     results_matrix = []

     # lil preprocessing
     groups = [ [] for i in range(10) ]
     for i in range(len(data_df.axes[1])):
          colname = data_df.columns[i]
          group_n = re.search(" \\| (\\d)", colname).group(1)
          data_df = data_df.rename(columns = {colname: colname[:colname.find(" |")]})
          groups[int(group_n)].append(i)
     groups = [x for x in groups if x != []]

     # do the stuff
     for i in range(len(groups)):
          group_data_indecies = sorted(groups[i])

          col_names = []

          for n in group_data_indecies[:-1]:
               # input: pair of time points
               data_t0 = list(data_df.iloc[:,n])
               data_t1 =  list(data_df.iloc[:,n + 1])
               dt = int(data_df.columns[n + 1]) - int(data_df.columns[i])

               cell_volume = get_cell_volume(init_cell_volume, int(data_df.columns[i]))

               coefficient_arr, beta_all = function_for_timepoint(data_t0, data_t1, dt, gene_key, tf_key)
               coefficient_matrix.append(coefficient_arr)
               beta_all_arr.append(beta_all)

               col_names.append(f"{data_df.columns[n]} to {data_df.columns[n + 1]}")

     # overall function fitting for all data points in an organisms
     beta_arr, covar = fit_function()

     # export to csv
     results_matrix = np.array(results_matrix)
     beta_names = tf_key
     results_df = pd.DataFrame(results_matrix.T, index = beta_names, columns = col_names)
     results_df.to_csv(f"{prokode_dir}/beta_training/results/{data_file[data_file.find("beta_training\\")+14:data_file.find(".py")]}_group{i}.csv")



