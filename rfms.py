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

ngstec_l = ["Illumina-Iseq","Illumina-MiniSeq", "Illumina-MiSeq", "Illumina-HiSeq", "Illumina-NovaSeq","IT-Chip",
            "BGI-MGISEQ-2000", "BGI-BGISEQ", "BGI-MGISEQ-T7", "PACBIO-Sequel"]
profiles_l = ["profile1", "profile2"]
sample_l = []
orgs_sample = []

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

profiles_file = "INIs/profiles.ini"
orgs_file = "INIs/orgs.ini"
profiles = configparser.ConfigParser()
profiles.read(profiles_file)
orgs = configparser.ConfigParser()
orgs.read(orgs_file)

def profile_arg_parser(profile):
    pass

def org_arg_new(org_l):
    pass


def orgs_arg_parser(orgs_l):
    '''
    :param orgs_l:
    :return:
    '''
    if len(orgs_l) >= 7:
        raise ValueError("Invalid Organism Entry")
    if len(orgs_l) < 3:
        orgs.set(orgs_l[0], "seq_nr", "orgs_l[1]")
    if len(orgs_l) < 4:
        orgs.set(orgs_l[0], "seq_min"), int(orgs_l[2].split(":")[0])
        orgs.set(orgs_l[0], "seq_max"), int(orgs_l[2].split(":")[1])
    if len(orgs_l) < 5:
        if orgs_l[3] in sample_l:
            orgs.set(orgs_l[0], "sample_type", orgs_l[3])
        else:
            raise ValueError("Invalid Sample Type")
    if len(orgs_l) < 6:
        if len(orgs_l[4].split(":")) == 4:
            orgs.set(orgs_l[0], "a", orgs_l[4].split(":")[0])
            orgs.set(orgs_l[0], "t", orgs_l[4].split(":")[1])
            orgs.set(orgs_l[0], "g", orgs_l[4].split(":")[2])
            orgs.set(orgs_l[0], "c", orgs_l[4].split(":")[3])
    if len(orgs_l) < 7:
        orgs.set(orgs_l[0], "freq_reps", orgs_l[5])

## main
def arg_parser():
    '''
    Argmumentes Parser and Reader of INI files
    '''
    #TODO Add verification to the existance of output files or if they are valid
    ngs = args.ngs
    if args.seed:
        random.seed(args.seed)
    if args.profie == False and args.org == False:
        raise ValueError("Insert a Valid Profile and/or ORG type")
    else:
        if args.profile:
            if args.profile in profiles.sections():
                profile_arg_parser(args.profile)
            else:
                raise ValueError("Insert a Valid Profile")
        if args.org:
            for entry in args.org:
                org_l = entry.split(";")
                if org_l[0] == "new":
                    org_arg_new(org_l)
                elif org_l[0] in orgs.sections():
                    orgs_arg_parser(org_l)
                else:
                    raise ValueError("Organism Not Valid.")
    if args.fasta:
        fasta = args.fasta
    else:
        fasta = ""
    if args.fastq:
        fastq = args.fastq
    else:
        fastq = ""


def main():
    arg_parser()
    print(args)

########################################################################################################################
if __name__ == '__main__':
    main()