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
import random, time, tqdm, sys
import numpy as np
import multiprocessing as mp



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
        not_unique = True
        while not_unique:
            org_seed = np.random.randint(1, 10 ** 6, sequence_size.size)
            if org_seed.size == np.unique(org_seed).size:
                not_unique=False

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
    i = 80
    while i > len(seq):
        seq = seq[:i] + "\n" + seq[i:]
    try:
        with open(fasta, "x") as file:
            file.write(header + "\n")
            #file.writelines(seq)
            file.write(seq)
    except:
        with open(fasta, "a") as file:
            file.write("\n")
            file.write(header + "\n")
            file.writelines(seq)
            file.write(seq)

def reverse_complement(sequence):
    '''

    :param sequence:
    :return:
    '''
    complement = {65:84,
                  84:65,
                  67:71,
                  71:67}
    new_sequence = bytearray()
    for base in sequence:
        new_sequence.extend(complement[base])
    return new_sequence[::-1]

def rep_sequences(sequence, reps):
    elem, elem_reps = reps
    elem_reps = [i for i in range(0, elem_reps)]
    rep_seqs = {}
    while len(elem_reps) > 0:
        for chunk in elem:
            for rep in chunk:
                if rep[2] in elem_reps:
                    if rep[4] == 0:
                        rep_seqs[rep[2]] = sequence[rep[0]:rep[1]]
                        rep_seqs[-int(rep[2])] = reverse_complement(rep_seqs[rep[2]])
                        elem_reps.remove(rep[2])
    return rep_seqs



def substitute_rep(sequence, list_pos):
    '''
    :param sequence:
    :param list_pos:
    :return:
    '''
    rep_seq = {}
    i = 0
    print(len(list_pos))
    for pos in list_pos:
        print(i)
        i+=1
        if pos[0] in rep_seq.keys():
            sub_rep = False
            if pos[3] == "-":
                sub_rep = reverse_complement(rep_seq[pos[2]])
            if sub_rep:
                sequence = sequence[:pos[0]] + sub_rep + sequence[(pos[1]):]
            else:
                sequence = sequence[:pos[0]] + sequence[pos[0]:pos[1]] + sequence[(pos[1]):]

        else:
            rep_seq[pos[0]] = sequence[pos[0]:(pos[1])]
    return sequence

def simulate_reps_pos(sequence_sizes, orgs_use, orgs):
    '''
    simulates sequence reps
    :return:
    '''
    reps = {}
    for i_org in range(0, len(orgs_use)):
        org = orgs_use[i_org]
        sequence_size = sequence_sizes[i_org].sum()
        chunks = len(sequence_sizes)
        elem_reps = random.randint(int(orgs[org]["min_rep_elem"]), int(orgs[org]["max_rep_elem"]))
        pos = [[] for i in range(0, chunks)]
        for elem in range(0, elem_reps):
            rep_size = random.randint(int(orgs[org]["min_rep_size"]), int(orgs[org]["max_rep_size"]))
            nr_reps = random.randint(int(orgs[org]["min_rep_nr"]), int(orgs[org]["max_rep_nr"]))
            direct = False
            for nr in range(0, nr_reps):
                if float(orgs[org]["prob_rep_direction"]) <= random.random():
                    direction = "+"
                else:
                    direction = "-"

                if direct:
                    pivot = end
                else:
                    pivot = random.randint(0, sequence_size-rep_size-1)
                end = pivot + rep_size

                chunk_p = pivot // 10**6
                chunk_e = end // 10**6
                pivot = pivot % 10 ** 6
                end = end % 10 ** 6
                if chunk_p == chunk_e:
                    pos[chunk_p].append((pivot, end + 1, elem, direction, 0))
                else:
                    pos[chunk_p].append((pivot, 10**6 + 1, elem, direction, 1))
                    pos[chunk_e].append(0, end + 1, elem, direction, 2)


                if float(orgs[org]["prob_direct_rep"]) <= random.random():
                    direct = True
                else:
                    direct = False
        for i in pos:
            random.shuffle(i)
        reps[org] = (pos, elem_reps)
    return reps


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
    seq = bytearray()
    for i in range(seq_size):
        prob = random.random()
        if prob < prob_seq[0]:
            seq.append(65)
        elif prob >= prob_seq[0] and prob < prob_seq[1]:
            seq.append(84)
        elif prob >= prob_seq[1] and prob < prob_seq[2]:
            seq.append(67)
        elif prob >= prob_seq[2] and prob < prob_seq[3]:
            seq.append(71)
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
    print("Sequence Sizes Generated")
    orgs_seeds = generate_seed_pool(sequence_sizes, seed)
    print("Seeds Generated")
    #reps = simulate_reps_pos(sequence_sizes, orgs_use, orgs)
    print("Repetitions generated")

    print("Using ", threads, " threads to generate sequences")
    for i_org in tqdm.trange(0, len(orgs_use), desc="Simulating Full Sample Sequence"):
        pool = mp.Pool(threads)
        org = orgs_use[i_org]
        probs = np.array([float(orgs[org]["a"]), float(orgs[org]["t"]), float(orgs[org]["c"]), float(orgs[org]["g"])])
        if probs.sum() != 1.0:
            print("The sum of base pair probabilities of ", org, "doesn't equal to 1.")
            probs_scaled = probs / probs.min()
            probs = probs_scaled / probs_scaled.sum()
        seq = pool.starmap(simulate_sequence,
                           [(probs, sequence_sizes[i_org][i_seq], orgs_seeds[i_org][i_seq])
                            for i_seq in range(0, sequence_sizes[i_org].size)])
        print(type(seq[0]))

        #seq = b''.join(seq).decode()
        #rep_seqs = rep_sequences(seq, reps[org])
        #seq = substitute_rep(seq, reps[org])
        #print("seq_subbed")
        pool.close()
        write_sequence(fasta, seq, org, probs, sequence_sizes[i_org].sum(), seed)


    print("end simulation")


def test ():
    import configparser
    orgs_file = "../INIs/orgs.ini"
    orgs = configparser.ConfigParser()
    orgs.read(orgs_file)
    seed = random.randint(1, 10 ** 6)
    random.seed(seed)
    simulate_sample(["HS", "MM"], orgs, seed, 6, "teste.fasta")

########################################################################################################################
if __name__ == '__main__':
    test()