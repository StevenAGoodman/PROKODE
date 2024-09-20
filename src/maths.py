### UNDER DEVELOPMENT ###
import numpy as np

def transcription_rate(gene, protein_amnts):

    return None
def translation_rate(gene, protein_amnts):
    return None
def decay_rate(molecule, gene_key, gene, molecule_amnts):
    if molecule == "mRNA":
        gene_index = gene_key.find(gene)
        gene_molecule_amnt = molecule_amnts[gene_index]
        degrading_proteins = # screening for decaying protiens
        for dp in degrading_proteins:
            N_dp = protein_amnts[dp]
            Kd_dp = 
            # some function of ndp kddp and nmrna

    elif molecule == "protein"

    else:
        raise Exeption("cannot find decay rate of specified molecule")
