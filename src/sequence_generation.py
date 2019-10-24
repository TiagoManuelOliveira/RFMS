#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"

## Imports
import random

def write_sequence(seq):
    pass

def simulate_reps ():
    '''
    simulates sequence reps
    :return:
    '''
    pass

def simulate_sequence(a, t, c, g, seq_size):
    '''
    simulate a sequence with n parametrs
    :return:
    '''
    for i in range(seq_size):
        seq = ""
        prob = random.random()
        if prob < a:
            seq += "A"
        elif prob >= a and prob < (t+a):
            seq += "T"
        elif prob >= (t+a) and prob < (c+t+a):
            seq += "C"
        elif prob >= (c+t+a) and prob < (c+t+a+g):
            seq += "G"
    return seq


def simulate_sample(orgs_use, orgs, fasta):
    '''
    simulate sample of n sequences
    :return:
    '''
    for org in orgs_use:
        a = float(orgs[org][a])
        t = float(orgs[org][t])
        c = float(orgs[org][c])
        g = float(orgs[org][g])
        seq_size = random.randint(orgs[org]["seq_min"], orgs[org]["seq_max"])
        seq = simulate_sequence(a, t, c, g, seq_size)

def main ():
    pass

########################################################################################################################
if __name__ == '__main__':
    random.seed(111)
    main()