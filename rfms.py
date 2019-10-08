#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__version__ = "RFMS - 0.01a"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"

## Imports
import argparse
import os
import configparser # use to store profiles in INI format
import random


## Argument Requirements

ngstec_l = ["myseq", "hiseq", "novaseq", "iontorrent", "pacbio"]
profiles_l = ["profile1", "profile2"]
sample_l = []


## Arguments Input

parser = argparse.ArgumentParser(description="Reference-free metagenomic simulator")
parser.add_argument("-s", "--seed", action="store",help="Set simulation seed", default=False)
parser.add_argument("-n", "--ngs", action="store", help="Type of NGS technology used", default=False, required=False, choices=ngstec_l)
parser.add_argument("-f", "--fasta", action="store", help="Output file name for FASTA file", default=False, required=False)
parser.add_argument("-fq", "--fastq", action="store", help="Output file name for FASTQ file", default=False, required=False)
parser.add_argument("-o", "--org", action="append", help="Organism info -> <org>;<nr de seqs>;<SeqMin:SeqMax>;<sampleType>;<FreqsA:T:G:C>;<freq Reps>", default=False, required=False)
parser.add_argument("-p", "--profile", choices=profiles_l, help="Use a predefined profile", default=False, action="append")
parser.add_argument("-v", "--version", action="version", version=__version__)

args = parser.parse_args()

## Files

working_dir = os.path.join(os.environ.get("HOME"), ".RFMS")
profiles_file = os.path.join(os.environ.get("HOME"),".RFMS", "profiles.ini")
orgs_file = os.path.join(os.environ.get("HOME"),".RFMS", "orgs.ini")

profiles = configparser.ConfigParser()
profiles.read(profiles_file)
orgs = configparser.ConfigParser()
orgs.read(orgs_file)


def orgs_arg_parser(orgs_l):
    '''

    :param orgs_l:
    :return:
    '''
    org = {}
    if len(orgs_l) == 1:
        pass
    if len(orgs_l) < 3:
        org["seq_nr"] = orgs_l[1]
    if len(orgs_l) < 4:
        org["seq_min"] = int(orgs_l[2].split(":")[0])
        org["seq_max"] = int(orgs_l[2].split(":")[1])
    if len(orgs_l) < 5:
        if orgs_l[3] in sample_l:
            org["sample_type"] = orgs_l[3]
        else:
            f"Invalid sample type"
            return False
    if len(orgs_l) < 6:
        if len(orgs_l[4].split(":")) == 4:
            org["A"] = orgs_l[4].split(":")[0]
            org["T"] = orgs_l[4].split(":")[1]
            org["G"] = orgs_l[4].split(":")[2]
            org["C"] = orgs_l[4].split(":")[3]
    if len(orgs_l) < 7:
        freq_reps = orgs_l[5]
    if len(orgs_l) >= 7:
        return False
    return org

## main
def arg_parser():
    '''
    Argmumentes Parser and Reader of INI files
    '''
    #TODO Add verification to the existance of output files or if they are valid
    #TODO Add verification to the existance of specific orgs in INI file
    if args.seed:
        random.seed(args.seed)
    if args.ngs:
        ngs = args.ngs
    if args.profile:
        profile = args.profile
    if args.org:
        orgs_dic = {}
        for entry in args.org:
            org_l = entry.split(";")
            if org_l[0] in orgs.sections():
                orgs_arg_parser_result = orgs_arg_parser(org_l)
                if orgs_arg_parser_result:
                    orgs_dic[org_l[0]] = orgs_arg_parser_result
                else:
                    return f"Organism entry not valid"
            else:
                return f"Organism not valid"

########################################################################################################################
if __name__ == '__main__':
    arg_parser()
    print(args)