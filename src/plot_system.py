### UNDER DEVELOPMENT ###
import json
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from maths import * 

# network.json formatting:
# {
#     "tg1":[
#         "tgdecay":float,
#         "regulators":[
#             "tf1": [
#                 "beta":float,
#                 "kd_tf":float,
#             ]
#             "tf2":[

#             ]
#             ...
#         ]
#     ],
#     "tg2":[

#     ],
#     ...
# } 

# constants
temperature = 200 # kelvin
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 0.1 # WHAT IS THE BINDING AFFINITY OF RIBO TO DALGARNO SEQ????
len_taken_by_rnap = 30 # nt
len_taken_by_ribo = 30 # aa

def score_to_K(score):
    delta_G = score
    return np.exp(delta_G / (1.98722 * temperature))

def get_mRNA_decay_info(network_dict):
    #search for mrna degrading genes: PNPase, Ribonuclease III, Ribonuclease Y, Ribonuclease G, Ribonuclease E, RNase J, RNase II, RNase R
    degrading_prot_ref = {
        "RNase E":{"names":["rne"], "score":0.1},
        "RNase G":{"names":[pnp], "score":0.1},
        "RNase Y":{"names":[rne; smbB; ams; hmp1; RNase E], "score":0.1},
        "RNase III":{"names":["rnc"], "score":0.1},
        "RNase J":{"names":[rne; smbB; ams; hmp1; RNase E], "score":0.1},
        "RNase R":{"names":["rnr", "vacB", "yjeC"], "score":0.1},
        "RNase II":{"names":["rnb"], "score":0.1},
        "PNPase":{"names":["pnp"], "score":0.1},
    }
    try: 
        network_dict[]
    return 

def transcription_rate(gene, protein_amnts, gene_key, N_rnap, Kd_rnap, gene_info_dict, genome_len):
    genome_len = 4.5e6

    # calculate beta_all
    beta_all = 0
    for tf, tf_info in gene_info_dict["regulators"].items():
        N_tf = protein_amnts[gene_key.index(tf)]
        P = N_tf / (N_tf + score_to_K(tf_info["score"]))
        beta = tf_info["beta"]
        beta_all += beta * P

    # calculate max transcription rate
    transcript_len = gene_info_dict["transcript length"]
    max_txn_rate = 1 / (len_taken_by_rnap / elongation_rate) # transcripts/s

    # calculate basal rnap binding prob
    P_rnap_basal = N_rnap / (genome_len * (N_rnap + Kd_rnap))

    return beta_all * P_rnap_basal * max_txn_rate

def translation_rate(gene, protein_amnts, N_ribo, gene_info_dict):
    # calculate max translation rate
    mRNA_len = gene_info_dict["mRNA length"]
    max_translation_rate = 1 / (len_taken_by_ribo / peptide_rate)

    # calculate ribo binding prob
    P_ribo_bound = N_ribo / (N_ribo + Kd_ribo_mrna)

    return P_ribo_bound * max_translation_rate

def RNA_decay_rate(prev_total_mRNA_amnt, protein_amnts, decay_dict, gene_key):
    # calculation of natural decay component
    natural_mRNA_decay_rate = 0

    # influence of degrading proteins
    # decay_dict contains protein geneids and kds
    rate_of_mRNA_cleavage = 0
    for degrad_prot, K in decay_dict.items():
        N_dp = protein_amnts[gene_key.index(degrad_prot)]
        # K_1 is the reaction rate from [Enzyme] + [mRNA] -> [Enzyme-mRNA] (ie, K_1[Enzyme][mRNA] = [Enzyme-mRNA] create / time ) ... k_1 is in 1 / (mols * time)
        rate_of_mRNA_cleavage += K * prev_total_mRNA_amnt * N_dp / (1 + K * prev_total_mRNA_amnt)
        # what is the relationship of multiple degrading proteins? it should be additive, right?
        # how can you identify a protien as degrading?
        # what is the relationship of degrading prots component and decay (ie half life) component

    return rate_of_mRNA_cleavage + natural_mRNA_decay_rate

def protein_decay_rate(gene, protein_amnts, decay_dict, gene_key):
    return 0.01

def update_amnts(gene_key:list, prev_mRNA_amnts:list, prev_protein_amnts:list, N_rnap, N_ribo, mRNA_decay_dict, network_dict, dt:float):
    mRNA_amnts = []
    protein_amnts = []

    mRNA_decay_rate = RNA_decay_rate(sum(prev_mRNA_amnts), prev_protein_amnts, mRNA_decay_dict, gene_key)

    for n in range(len(gene_key)):
        gene = gene_key[n]
        gene_info_dict = network_dict[gene]

        # mRNA calculation
        mRNA_creation_rate = transcription_rate(gene, prev_protein_amnts, gene_key, N_rnap, 0.1, gene_info_dict, genome_len=0)
        total_mRNA_rate = mRNA_creation_rate - mRNA_decay_rate
        mRNA_amnt = prev_mRNA_amnts[n] + total_mRNA_rate * dt

        # protein calculation
        prot_creation_rate = translation_rate(gene, prev_protein_amnts)
        prot_decay_rate = protein_decay_rate("protein", gene, prev_protein_amnts)
        total_prot_rate = prot_creation_rate - prot_decay_rate * prev_protein_amnts[n]
        protein_amnt = prev_protein_amnts[n] + total_prot_rate * dt

        mRNA_amnts.append(mRNA_amnt)
        protein_amnts.append(protein_amnt)

    # update rnap and ribosome populations

    return mRNA_amnts, protein_amnts, N_rnap, N_ribo

def get_plotting_data(network_dict, gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt):
    x_coords = [0]
    y_coords = [init_protein_amnts]

    # initialize variables
    time = 0
    mRNA_amnts = init_mRNA_amnts
    protein_amnts = init_protein_amnts
    N_rnap = 10000
    N_ribo = 15000

    mRNA_decay_dict = get_mRNA_decay_info()

    for i in range(math.floor(max_time / dt)):
        time += dt
        x_coords.append(time)
        
        mRNA_amnts, protein_amnts, N_rnap, N_ribo = update_amnts(gene_key, mRNA_amnts, protein_amnts, N_rnap, N_ribo, mRNA_decay_dict, network_dict, dt)

        y_coords.append(protein_amnts)

    return x_coords, y_coords

def plot_system(network_loc,  dt, max_time):
    # read inputs files
    network_dict:dict = json.load(open(network_loc, 'r'))
    gene_key = list(network_dict.keys())

    # default to no initial mRNA
    # try: read file except:
    init_mRNA_amnts = [0] * len(gene_key)

    # read input init protein file
    init_protein_amnts = [np.random.randint(0,1000) for i in range(len(gene_key))]

    x_coords, y_coords_tuple = get_plotting_data(network_dict, gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt)

    for i in range(2):
        plt.plot(x_coords,[gene_exps[i] for gene_exps in y_coords_tuple], random.choice(['r-','g-','b-']), linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()

plot_system('C:/Users/cryst/LOFScreening/archive/PROKODE/src/network.json',1,10)