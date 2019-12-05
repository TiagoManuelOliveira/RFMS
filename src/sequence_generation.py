#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__version__ = "0.2b"

## Imports
import random, tqdm
import numpy as np
import multiprocessing as mp
import bitarray

#TODO: Substitutir bytearrays por bitarrays

def generate_seed_pool(sequence_sizes, seed):
    '''
    Generates a pool of seeds to be used in multiprocessing of sequence generation.
    It uses numpy randint to create an array of random numbers to be used as seeds to be fed into the
    sequence generation def to insure that each each thread has a reproducible result. The function also insures
    that each int in the array is unique so that you don't get equal sequences chunks
    :param sequence_sizes: list of numpy arrays with the sizes of each sequence chuck for each org
    :param seed: seed used in random to be fed into numpy internal random generator
    :return: list of arrays with a seed for each sequence chunk for each org
    '''
    np.random.seed(seed)
    orgs_seeds = []
    for sequence_size in sequence_sizes:
        not_unique = True
        while not_unique: # Creates random array of numbers and checks if every seed is unique
            org_seed = np.random.randint(1, 10 ** 6, sequence_size.size)
            if org_seed.size == np.unique(org_seed).size:
                not_unique=False

        orgs_seeds.append(org_seed)

    return orgs_seeds



def generate_sequence_size(orgs_use, orgs, size=10**6):
    '''
    Generates a list of arrays with the sequence length of each sequence chunk with a maximum size of 10**6, for each
    org to be simulated.
    The full sequence size range is taken from iniparser object containing the seq_min and seq_max for all of the orgs
    to be simulated and then that size is divided in sequence chunks with a maximum size of 10**6 to be fed into
    each thread of sequence generation.
    :param orgs_use: list of orgs to be used
    :param orgs: iniparse orgs database
    :param size: Int size of each chunk of the DNA sequence
    :return: list of arrays with sequence size of each chunk for each org to be simulated
    '''
    sequence_sizes = []
    for org in orgs_use:
        pieces = []
        seq_size = random.randint(int(orgs[org]["seq_min"]), int(orgs[org]["seq_max"]))
        for i in range(0, seq_size // size): # All chunks with size = size
            pieces.append(size)
        extra_piece = seq_size % size
        if extra_piece != 0: # Last chunk size
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

def write_sequence(fasta, seq, org, probs, sequence_size, seed, linebreak=80):
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
        with open(fasta, "x") as file: # If FASTA file doesn't exist. "X" mode throws error if file exists
            file.write(header + "\n")
            i = 0
            j = linebreak
            while i < len(seq):
                file.write(seq[i:j] + "\n")
                i+=linebreak
                j+=linebreak
    except:
        with open(fasta, "a") as file: # If Fasta file exists, appends new sequence to end of file.
            i = 0
            j = linebreak
            file.write("\n")
            file.write(header + "\n")
            while i < len(seq):
                file.write(seq[i:j] + "\n")
                i += linebreak
                j += linebreak

def reverse_complement(sequence):
    '''
    Function that creates the reverse complement of a DNA sequence in bytearray format
    :param sequence: DNA sequence in Bytearray format
    :return: Reverse complement DNA sequence in bytearray format
    '''
    # complement dic has the int encoding the DNA letters in Byte format with 255 encoding
    complement = {65:84,
                  84:65,
                  67:71,
                  71:67}
    new_sequence = bytearray()
    for base in sequence:
        new_sequence.append(complement[base])
    return new_sequence[::-1]

def rep_sequences(sequence, reps):
    '''
    Stores all repetition sequences for an org
    :param sequence: DNA sequence used to obtain the rep sequences
    :param reps:Tuple with list of list of tuples containing the positions of each repetition. The 2nd element of the
    tuple is an int with the number of repetition elements
    :return: dictionary with key being the rep element nr and value being the DNA sequence. negative keys mean
    the complementary reverse sequence of the rep element
    '''
    elem, elem_reps = reps
    elem_reps = [i for i in range(0, elem_reps)] # List that just stores the elements to be obtained from main sequence
    rep_seqs = {}
    while len(elem_reps) > 0: # While there still exists elements to be obtained
        i_chunk = 0
        for chunk in elem:
            for rep in chunk:
                if rep[2] in elem_reps:
                    if rep[4] == 0:
                        rep_seqs[rep[2]] = sequence[i_chunk][rep[0]:rep[1]]
                        rep_seqs[-int(rep[2])] = reverse_complement(sequence[i_chunk][rep[0]:rep[1]])
                        elem_reps.remove(rep[2])
            i_chunk += 1
    return rep_seqs



def substitute_rep(arguments):
    '''
    Function introduces repetition sequences into a DNA sequence
    :param sequence: DNA sequence, intended to be be in Bytearray format but it can also be used with any iterable type
    :param pos: list of tuples containing the positions and repetition information (pivot, end, repetition element
    number, direction and if [full rep, 1st part of rep or last part of rep])
    :param elem: int with number of repetitions
    :param rep_seqs: dictionary where keys are rep elem numbers and values are the rep DNA sequences in iterable
    format. Negative keys represent the reverse complementary of the repetition DNA sequence.
    :return: DNA sequence with repetitions inserted into respective positions
    '''
    sequence, pos, elem, rep_seqs = arguments
    for rep in pos:
        if rep[3] == "+":
            sub = rep_seqs[rep[2]]
        else:
            sub = rep_seqs[-int(rep[2])]
        if rep[4] > 0:
            sub_part_size = rep[1] - rep[0]
            if rep[4] == 1:
                sub = sub[0:sub_part_size]
            elif rep[4] == 2:
                size_sub = len(sub)
                sub = sub[size_sub-sub_part_size:size_sub]

        sequence[rep[0]:rep[1]] = sub

    return sequence


def getKey(item):
    '''
    helper function to order iterable inside another iterable by its first position (ex: order matrix by the first
    column value)
    :param item: iterable
    :return: iterable first element
    '''
    return item[0]

def simulate_reps_pos(sequence_sizes, orgs_use, orgs):
    '''
    Simulates repetition positions for all organisms in sample
    :param sequence_sizes: list of arrays with size of each chunk of DNA for each Organism
    :param orgs_use: Organisms in sample
    :param orgs: INIparser object with information of every organism
    :return:
    '''
    #TODO: Optimize this function. It takes way too much time. Maybe some paralelization
    #TODO: Function with large numbers ocupies too much space. Try to substitute some flags by boleans or bytes
    reps = {}
    for i_org in tqdm.trange(0, len(orgs_use), desc="Simulating sequences for organisms", position=0):
        org = orgs_use[i_org]
        sequence_size = sequence_sizes[i_org].sum()
        chunks = sequence_sizes[i_org].size
        elem_reps = random.randint(int(orgs[org]["min_rep_elem"]), int(orgs[org]["max_rep_elem"]))
        pos = [[] for i in range(0, chunks)]
        for elem in tqdm.trange(0, elem_reps, desc="Generating repetitions information for organism", position=1):
            rep_size = random.randint(int(orgs[org]["min_rep_size"]), int(orgs[org]["max_rep_size"]))
            nr_reps = random.randint(int(orgs[org]["min_rep_nr"]), int(orgs[org]["max_rep_nr"]))
            direct = False # Tandem repeat flag
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
                    pos[chunk_e].append((0, end + 1, elem, direction, 2))

                if float(orgs[org]["prob_direct_rep"]) <= random.random():
                    direct = True
                else:
                    direct = False
        for i in range(0,len(pos)):
            random.shuffle(pos[i])
            pos[i] = sorted(pos[i], key=getKey)

        reps[org] = (pos, elem_reps)
    return reps


def simulate_sequence(arguments):
    '''
    Simulates a sequence from a given seed using a prob array and returns a string with a DNA sequence.
    :param prob_seq: array with probabilites for each base
    :param seq_size: int with the sequence size to be generated
    :param seed: int with the seed to be used by the generator
    :return: Byearray with the generated DNA sequence
    '''
    prob_seq, seq_size, seed = arguments
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


def simulate_sample(orgs_use, orgs, seed, threads, fasta, linebreak=80):
    '''
    Simulates the DNA sequence of each organism in a given metagenomic sample. It handles multithreading and ensures
    reproducibility.
    It saves the end sequences in a Fasta File with a default line break of 80 characters in the DNA sequence
    :param orgs_use: List of strings of the orgs to be used
    :param orgs: iniparser object containing information about every org the database
    :param seed: interger with the seed to be used in the simulation
    :param threads: number of threads to be used in the simulation
    :param fasta: string with the name of the FASTA file
    :return: Writes in a FASTA file the simulated DNA sequences
    '''
    sequence_sizes = generate_sequence_size(orgs_use, orgs)
    print("Sequence Sizes Generated")
    orgs_seeds = generate_seed_pool(sequence_sizes, seed)
    print("Seeds Generated")
    reps = simulate_reps_pos(sequence_sizes, orgs_use, orgs)
    print("Repetitions generated")

    print("Using ", threads, " threads to generate sequences")
    for i_org in tqdm.trange(0, len(orgs_use), desc="Simulating Organism Sequence", position=0):
        org = orgs_use[i_org]
        pos, elem = reps[org]
        probs = np.array([float(orgs[org]["a"]), float(orgs[org]["t"]), float(orgs[org]["c"]), float(orgs[org]["g"])])
        if probs.sum() != 1.0:
            print("The sum of base pair probabilities of ", org, "doesn't equal to 1.")
            probs_scaled = probs / probs.min()
            probs = probs_scaled / probs_scaled.sum()

        generation_arguments = []
        for i_seq in range(0,sequence_sizes[i_org].size):
            generation_arguments.append((probs, sequence_sizes[i_org][i_seq], orgs_seeds[i_org][i_seq]))
        with mp.Pool(threads) as pool:
            seq = list(tqdm.tqdm(pool.imap(simulate_sequence, generation_arguments),
                                 total=sequence_sizes[i_org].size, desc="Chunks of sequence simulated", position=1 ))

        rep_seqs = rep_sequences(seq, reps[org])
        rep_arguments = []
        for i_rep in range(0,len(seq)):
            rep_arguments.append((seq[i_rep], pos[i_rep], elem, rep_seqs))
        with mp.Pool(threads) as pool:
            new_seq = list(tqdm.tqdm(pool.imap(substitute_rep, rep_arguments),
                                     total=len(seq), desc="Inserting repetitions", position=2))
        seq = b''.join(new_seq).decode()
        write_sequence(fasta, seq, org, probs, sequence_sizes[i_org].sum(), seed, linebreak)

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