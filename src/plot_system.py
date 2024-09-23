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
elongation_rate = 60 # nt/s
peptide_rate = 20 # aa/s
Kd_ribo_mrna = 

def transcription_rate(gene, protein_amnts, gene_key, N_rnap, Kd_rnap, gene_info_dict, genome_len):
    genome_len = 4.5e6

    # calculate beta_all
    beta_all = 0
    for tf, tf_info in gene_info_dict["regulators"].items():
        N_tf = protein_amnts[gene_key.index(tf)]
        P = N_tf / (N_tf + tf_info["kd_tf"])
        beta = tf_info["beta"]
        beta_all += beta * P

    # calculate max transcription rate
    transcript_len = gene_info_dict["transcript length"]
    max_txn_rate = 1 / (transcript_len / elongation_rate) # transcripts/s

    # calculate basal rnap binding prob
    P_rnap_basal = N_rnap / (genome_len * (N_rnap + Kd_rnap))

    return beta_all * P_rnap_basal * max_txn_rate

def translation_rate(gene, protein_amnts, N_ribo, gene_info_dict):
    # calculate max translation rate
    mRNA_len = gene_info_dict["mRNA length"]
    max_translation_rate = 1 / (mRNA_len * peptide_rate)

    # calculate ribo binding prob
    P_ribo_bound = N_ribo / (N_ribo + Kd_ribo_mrna)

    return P_ribo_bound * max_translation_rate

def decay_rate(gene, protein_amnts, decay_dict, gene_key):
    # calculation of natural decay component


    # influence of degrading proteins
    # decay_dict contains protein geneids and kds
    for degrad_prot, kd in decay_dict.items():
        N_dp = protein_amnts[gene_key.index(degrad_prot)]
        # what is the relationship of multiple degrading proteins? it should be additive, right?
        # how can you identify a protien as degrading?
        # what is the relationship of degrading prots component and decay (ie half life) component

def update_amnts(gene_key:list, prev_mRNA_amnts:list, prev_protein_amnts:list, N_rnap, N_ribo, network_dict, dt:float):
    mRNA_amnts = []
    protein_amnts = []

    for n in range(len(gene_key)):
        gene = gene_key[n]

        # mRNA calculation
        mRNA_creation_rate = transcription_rate(gene, prev_protein_amnts)
        mRNA_decay_rate = decay_rate("mRNA", gene, prev_protein_amnts)
        total_mRNA_rate = mRNA_creation_rate - mRNA_decay_rate * prev_mRNA_amnts[n]
        mRNA_amnt = prev_mRNA_amnts[n] + total_mRNA_rate * dt

        # protein calculation
        prot_creation_rate = translation_rate(gene, prev_protein_amnts)
        prot_decay_rate = decay_rate("protein", gene, prev_protein_amnts)
        total_prot_rate = prot_creation_rate - prot_decay_rate * prev_protein_amnts[n]
        protein_amnt = prev_protein_amnts[n] + total_prot_rate * dt

        mRNA_amnts.append(mRNA_amnt)
        protein_amnts.append(protein_amnt)

    # update rnap and ribosome populations
    rate_rnap_decay = decay_rate("polymerase",)

    return mRNA_amnts, protein_amnts, N_rnap, N_ribo

def get_plotting_data(gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt):
    x_coords = [0]
    y_coords = [init_protein_amnts]

    # initialize variables
    time = 0
    mRNA_amnts = init_mRNA_amnts
    protein_amnts = init_protein_amnts
    N_rnap = 10000
    N_ribo = 15000

    for i in math.floor(max_time / dt):
        time += dt
        x_coords.append(time)
        
        mRNA_amnts, protein_amnts, N_rnap, N_ribo = update_amnts(gene_key, mRNA_amnts, protein_amnts, N_rnap, N_ribo, Kd_rnap, dt)

        y_coords.append(protein_amnts)

    return x_coords, y_coords

def plot_system(gene_key, init_mRNA_amnts, init_protein_amnts, dt, max_time):
    x_coords, y_coords_tuple = get_plotting_data(gene_key, init_mRNA_amnts, init_protein_amnts, max_time, dt)

    for i in range(10):
        plt.plot(t,s[:,i],random.choice(['r-','g-','b-']), linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()
