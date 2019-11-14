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
__description__ = "Reference-free metagenomic simulator"

## Imports
import argparse
from datetime import datetime
import configparser # use to store profiles in INI format
import random
import src.sequence_generation as sequence_generation
import multiprocessing as mp

## Argument Requirements
ngstec_l = ["Illumina-Iseq","Illumina-MiniSeq", "Illumina-MiSeq", "Illumina-HiSeq", "Illumina-NovaSeq","IT-Chip",
            "PACBIO-Sequel"]
profiles_l = ["profile1", "profile2"]
sample_l = []
orgs_sample = []

## Needed Objects
orgs_use = []

## Arguments Input

def arg_parser():
    parser = argparse.ArgumentParser(prog="rfms",
                                     description=__version__ + " - " + __description__)
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument("-s", "--seed", action="store", help="Set simulation seed", default=False)
    parser.add_argument("-t", "--threads", action="store", help="Number of threads to use. If none is inputed, program "
                                                                "will use 1 thread", default= 1)
    generation = parser.add_argument_group("Generation", "Arguments related to Sequence Generation")
    generation.add_argument("-f", "--fasta", action="store", help="Output file name for FASTA file", default=False,
                        required=False)
    generation.add_argument("-o", "--org", action="append",
                            help="Organism info -> <org>;<nr de seqs>;<SeqMin:SeqMax>;<sampleType>;"
                                 "<FreqsA:T:G:C>;<freq Reps>", default=False, required=False)
    generation.add_argument("-p", "--profile", choices=profiles_l, help="Use a predefined profile", default=False,
                            action="append")
    sequencing = parser.add_argument_group("Sequencing", "Arguments related to Sequencing Simulation")
    sequencing.add_argument("-fq", "--fastq", action="store", help="Output file name for FASTQ file", default=False,
                            required=False)
    sequencing.add_argument("-n", "--ngs", action="store", help="Type of NGS technology used", default=False,
                            required=False, choices=ngstec_l)
    return parser.parse_args()

## Files
profiles_file = "INIs/profiles.ini"
orgs_file = "INIs/orgs.ini"
profiles = configparser.ConfigParser()
profiles.read(profiles_file)
orgs = configparser.ConfigParser()
orgs.read(orgs_file)

def profile_arg_parser(profile):
    '''
    Function checks "profile" section from profiles ini file and changes parameters from
    "orgs" config object. If org is not present in "orgs" object sections it adds them.
    :param profile: str of name section present in "profiles.ini"
    '''
    profile_orgs_nr = len(profiles[profile]["org"].split(";"))
    items = ["a", "t", "g", "c", "freq_reps", "seq_min", "seq_max", "seq_nr", "sample_type"]
    for i in range(0, profile_orgs_nr):
        org_used = profiles[profile]["org"].split(";")[i]
        orgs_use.append(org_used)
        if org_used in orgs.sections():
            for item in items:
                if profiles[profile][item].split(";")[i].strip(" ") != "":
                    orgs.set(org_used, item, profiles[profile][item].split(";")[i].strip(" "))
        else:
            orgs.add_section(org_used)
            for item in items:
                if profiles[profile][item].split(";")[i].strip(" ") != "":
                    orgs.set(org_used, item, profiles[profile][item].split(";")[i].strip(" "))
                else:
                    raise ValueError("Invalid entry in profile.ini. Please check", profile, "values.")

def org_arg_new(orgs_l):
    '''
    Function adds new a new section to "orgs" config object and sets its items.
    If the org information is incomplete it raises a ValueError.
    :param orgs_l: list of cli args with the information about the org.
    '''
    if len(orgs_l) == 6:
        orgs.add_section(orgs_l[0])
        orgs.set(orgs_l[0], "seq_nr", int(orgs_l[1]))
        orgs.set(orgs_l[0], "seq_min"), int(orgs_l[2].split(":")[0])
        orgs.set(orgs_l[0], "seq_max"), int(orgs_l[2].split(":")[1])
        if orgs_l[3] in sample_l:
            orgs.set(orgs_l[0], "sample_type", orgs_l[3])
        else:
            raise ValueError("Invalid Sample Type")
        if len(orgs_l[4].split(":")) == 4:
            orgs.set(orgs_l[0], "a", float(orgs_l[4].split(":")[0]))
            orgs.set(orgs_l[0], "t", float(orgs_l[4].split(":")[1]))
            orgs.set(orgs_l[0], "g", float(orgs_l[4].split(":")[2]))
            orgs.set(orgs_l[0], "c", float(orgs_l[4].split(":")[3]))
        else:
            raise ValueError("Invalid Basepair Frequencies")
        orgs.set(orgs_l[0], "freq_reps", int(orgs_l[5]))
    else:
        raise ValueError("Invalid Organism Entry")

def orgs_arg_parser(orgs_l):
    '''
    Function that changes "orgs" config object values if section exists.
    The fucntion changes the values if values are given
    :param orgs_l: list of cli args with the information about the org.
    '''
    if len(orgs_l) >= 7:
        raise ValueError("Invalid Organism Entry")
    if len(orgs_l) < 3:
        if len(orgs_l[1].strip(" ")) > 0:
            orgs.set(orgs_l[0], "seq_nr", int(orgs_l[1]))
    if len(orgs_l) < 4:
        if len(orgs_l[2].strip(" ")) > 0:
            orgs.set(orgs_l[0], "seq_min"), int(orgs_l[2].split(":")[0])
            orgs.set(orgs_l[0], "seq_max"), int(orgs_l[2].split(":")[1])
    if len(orgs_l) < 5:
        if len(orgs_l[3].strip(" ")) > 0:
            if orgs_l[3] in sample_l:
                orgs.set(orgs_l[0], "sample_type", orgs_l[3])
            else:
                raise ValueError("Invalid Sample Type")
    if len(orgs_l) < 6:
        if len(orgs_l[4].strip(" ")) > 0:
            if len(orgs_l[4].split(":")) == 4:
                orgs.set(orgs_l[0], "a", float(orgs_l[4].split(":")[0]))
                orgs.set(orgs_l[0], "t", float(orgs_l[4].split(":")[1]))
                orgs.set(orgs_l[0], "g", float(orgs_l[4].split(":")[2]))
                orgs.set(orgs_l[0], "c", float(orgs_l[4].split(":")[3]))
            else:
                raise ValueError("Invalid Basepair Frequencies")
    if len(orgs_l) < 7:
        if len(orgs_l[5].strip(" ")) > 0:
            orgs.set(orgs_l[0], "freq_reps", float(orgs_l[5]))
## main
def arg_handler():
    '''
    Argmumentes Parser and Reader of INI files
    :return: :ngs: str with ngs technology to be simulated
    :return: :seed: int seed value to be used in the simulation. If no seed is given in arg, seed = False
    :return :fasta: str with name of fasta file
    :return :fastq str with name of fastq file
    '''
    args = arg_parser()
    ngs = args.ngs
    if args.seed:
        random.seed(args.seed)
        seed = args.seed
    else:
        seed = random.randint(1,10**6)
        random.seed(seed)
    if args.threads > mp.cpu_count():
        threads = mp.cpu_count()
        print("Number of threads indicated exceeds number of available threads, Program will use ",
              mp.cpu_count(), " instead")
    else:
        threads = args.threads
    if args.profile == False and args.org == False:
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
                if org_l[0] in orgs.sections():
                    orgs_arg_parser(org_l)
                else:
                    org_arg_new(org_l)
                orgs_use.append(org_l[0])
    if args.fasta:
        fasta = args.fasta
    else:
        fasta = "rfms_" + str(datetime.datetime.now()).replace(" ", "_")[:-10]+".fasta"
    if args.fastq:
        fastq = args.fastq
    else:
        fastq = "rfms_" + str(datetime.datetime.now()).replace(" ", "_")[:-10]+".fastq"
    return (ngs, fasta, fastq, threads)


def main():
    ngs, seed, fasta, fastq, threads = arg_handler()
    sequence_generation.simulate_sample(orgs_use, orgs, seed, threads, fasta)

########################################################################################################################
if __name__ == '__main__':
    main()