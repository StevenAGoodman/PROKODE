import json
import random
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd

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

def config_network_json(gene_arr, decay_file, tfbs_file):
    # needed inputs: lists of tgs, gene - decay rates ref file,

    # implicit file dependencies: 
    decay_df = pd.read_csv(decay_file) # names=['gene', 'prot decay rate']
    topology_df = pd.read_csv(tfbs_file) # names=['tf','tg','kd', 'beta]
    topology_df = topology_df.reset_index()

    output = {}

    for tg in gene_arr:
        tgdecay = decay_df.loc[decay_df['gene']==tg, 'decay']# query decay rates file

        arr = {"tg_decay":tgdecay.to_list()[0]}
        reg_arr = {}
        for _, tf_row in topology_df[topology_df['tg']==tg].iterrows():
            beta = tf_row['beta']
            kdtf = tf_row['kd']
            reg_arr[tf_row['tf']] = {"beta":beta,"kd_tf":kdtf}
        arr["regulators"] = reg_arr
        
        output[tg] = arr

    json.dump(output, open('network.json', 'w'), indent=4)

# def get_basal(Nns):
    # Np = 
    # Kdp =
    return 0. # Np / (Nns * Kdp)

def plot_system(params, s0, tg_arr):
    # check if params len mathces nodes

    def f(s,t):
        print(s)
        diffeqs = []

        for tg in params:
            param = params[tg]
            tg_decay = param['tg_decay']
            regulators = param['regulators']
            
            probspace_taken = 0 ## NOT ACCURATE
            R_trans = 0
            for tf in regulators:
                N_tf = s[tg_arr.index(tf)]
                print(tg_arr.index(tf))
                beta = regulators[tf]["beta"]
                kd_tf = regulators[tf]["kd_tf"]
                bindProb = N_tf / (N_tf + kd_tf)
                assert N_tf >= 0, (N_tf, f'{s}')

                R_trans += bindProb * (beta * E_basal)
                probspace_taken = bindProb + probspace_taken - (probspace_taken * bindProb)
            print(probspace_taken)
            assert probspace_taken < 1
            probspace_remaining = 1 - probspace_taken
            R_trans += (probspace_remaining) * E_basal
            d = R_trans - tg_decay * s[tg_arr.index(tg)]
            diffeqs.append(d)
        return diffeqs

    t = np.arange(0,1000,0.1)

    s = odeint(f,s0,t)
    print(t)

    for i in range(10):
        plt.plot(t,s[:,i],random.choice(['r-','g-','b-']), linewidth=2.0)
    # plt.plot(t,s[:,2], 'g-',linewidth=2.0)
    plt.xlabel("t")
    plt.ylabel("S[N,C]")
    plt.legend(["N","C",'d'])
    plt.show()

genome_df = pd.read_csv('./inputs/annotation.csv')
gene_arr = list(set(genome_df['geneid'].to_list()))
# config_network_json(gene_arr,'decay_rates.csv', 'tfbs.csv')
params = json.load(open('network.json', 'r'))
s0 = [random.randint(0,100) for i in range(len(gene_arr))] 
E_basal = 0.00434782608696
print(len(params))
plot_system(params, s0, gene_arr)