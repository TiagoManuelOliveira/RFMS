#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__version__ = "0.1b"

## Imports
import random
import numpy as np
import multiprocessing as mp
import time


def generate_seed_pool(sequence_sizes, seed):
    '''
    Generates a pool of seeds to be used in multiprocessing of sequence generation.
    It uses numpy randint to create an array of random numbers to be used as seeds to be fed into the
    sequence generation def to insure that each each thread has a reproductable result. The function also insures
    that each int in the array is unique so that you dont get equal sequences chunks
    :param sequence_sizes: list of numpy arrays with the sizes of each sequence chuck for each org
    :param seed: seed used in random to be fed into numpy internal random generator
    :return: list of arrays with a seed for each sequence chunk for each org
    '''
    np.random.seed(seed)
    orgs_seeds = []

    for sequence_size in sequence_sizes:
        print(sequence_size.size)
        not_unique = True
        not_verified = True
        while not_unique and not_verified:
            org_seed = np.random.randint(1, 10 ** 6, sequence_size.size)
            not_unique = org_seed.size == np.unique(org_seed).size
            i = 0
            for seeds in orgs_seeds:
                for i in seeds:
                    if i in org_seed:
                        i += 1
            if i == 0:
                not_verified = False

        orgs_seeds.append(org_seed)
    return orgs_seeds



def generate_sequence_size(orgs_use, orgs):
    '''
    Generates a list of arrays with the sequence length of each sequence chunk with a maximum size of 10**6, for each
    org to be simulated.
    The full sequence size range is taken from iniparser object containing the seq_min and seq_max for all of the orgs
    to be simulated and then that size is divided in sequence chunks with a maximum size of 10**6 to be fed into
    each thread of sequence generation.
    :param orgs_use: list of orgs to be used
    :param orgs: iniparse orgs database
    :return: list of arrays with sequence size of each chunk for each org to be simulated
    '''
    size = 10**6
    sequence_sizes = []
    for org in orgs_use:
        pieces = []
        seq_size = random.randint(int(orgs[org]["seq_min"]), int(orgs[org]["seq_max"]))
        for i in range(0, seq_size // size):
            pieces.append(size)
        extra_piece = seq_size % size
        if extra_piece != 0:
            pieces.append(extra_piece)
        sequence_sizes.append(np.asarray(pieces))
    return sequence_sizes

def generate_header(org, probs, sequence_size, seed):
    '''
    Generates the header of each sequence in the FASTA file
    :param org: string of org simulated
    :param probs: array of probabilities used in the sequence generation
    :param sequence_size: string of size of the generated sequence
    :param seed: string of seed used in the generation of the sequence
    :return: string with the header of the sequence to be written in the fasta file
    '''
    probs = str(probs).replace("[","").replace("]","").replace(",",":")
    header = ">"
    header += str(org + "|")
    header += str(probs.replace(" ", ":") + "|")
    header += str(str(sequence_size) + "|")
    header += str(seed)
    return header

def write_sequence(fasta, seq, org, probs, sequence_size, seed):
    '''
    writter function of the sequence
    :param fasta: string with the name of the fasta file
    :param seq: string of the full sequence to be written
    :param org: string with the org name
    :param probs: array of probabilites used in the generation of seq string
    :param sequence_size: interger with the size of the simualted sequence
    :param seed: interger with the seed used in the generation of the sequence
    :return:
    '''
    header = generate_header(org, probs, sequence_size, seed)
    try:
        with open(fasta, "x") as file:
            file.write(header + "\n")
            file.write(seq)
    except:
        with open(fasta, "a") as file:
            file.write("\n\n")
            file.write(header + "\n")
            file.write(seq)



def simulate_reps ():
    '''
    simulates sequence reps
    :return:
    '''
    pass

def simulate_sequence(prob_seq, seq_size, seed):
    '''
    Simulates a sequence from a given seed using a prob array and returns a string with a DNA sequence.
    :param prob_seq: array with probabilites for each base
    :param seq_size: int with the sequence size to be generated
    :param seed: int with the seed to be used by the generator
    :return: string with the generated DNA sequence
    '''
    random.seed(seed)
    prob_seq = prob_seq.cumsum()
    seq = ""
    for i in range(seq_size):
        prob = random.random()
        if prob < prob_seq[0]:
            seq += "A"
        elif prob >= prob_seq[0] and prob < prob_seq[1]:
            seq += "T"
        elif prob >= prob_seq[1] and prob < prob_seq[2]:
            seq += "C"
        elif prob >= prob_seq[2] and prob < prob_seq[3]:
            seq += "G"
    return seq


def simulate_sample(orgs_use, orgs, seed, threads, fasta):
    '''
    Simulates the DNA sequence of each organism in a given metagenomic sample. It handles multithreading and ensures
    reproductability.
    It saves the end sequences in a Fasta File
    :param orgs_use: List of strings of the orgs to be used
    :param orgs: iniparser object containing information about every org the database
    :param seed: interger with the seed to be used in the simulation
    :param threads: number of threads to be used in the simulation
    :param fasta: string with the name of the FASTA file
    :return: ####
    '''
    sequence_sizes = generate_sequence_size(orgs_use, orgs)
    orgs_seeds = generate_seed_pool(sequence_sizes, seed)
    pool = mp.Pool(threads)
    for i_org in range(0, len(orgs_use)):
        seq = ""
        org = orgs_use[i_org]
        probs = np.array([float(orgs[org]["a"]), float(orgs[org]["t"]), float(orgs[org]["c"]), float(orgs[org]["g"])])
        if probs.sum() != 1.0:
            print("The sum of base pair probabilities of ", org, "doesn't equal to 1.")
            probs_scaled = probs / probs.min()
            probs = probs_scaled / probs_scaled.sum()
        seq = pool.starmap(simulate_sequence,
                           [(probs, sequence_sizes[i_org][i_seq], orgs_seeds[i_org][i_seq])
                            for i_seq in range(0, sequence_sizes[i_org].size)])

        seq = "".join(seq)
        write_sequence(fasta, seq, org, probs, sequence_sizes[i_org].sum(), seed)
    pool.close()


def main ():
    import configparser
    orgs_file = "../INIs/orgs.ini"
    orgs = configparser.ConfigParser()
    orgs.read(orgs_file)
    seed = random.randint(1, 10 ** 6)
    random.seed(seed)
    simulate_sample(["HS", "MM"], orgs, seed, 6, "teste.fasta")

########################################################################################################################
if __name__ == '__main__':
    main()
