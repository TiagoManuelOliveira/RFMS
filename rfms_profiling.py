#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__version__ = "RFMS_profiling - 0.1"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"
__description__ = "Reference-free metagenomic simulator profiling tool"

## Imports
import argparse
from datetime import datetime
import configparser # use to store profiles in INI format
import art
from src.generate_profiles import profile_report


parser = argparse.ArgumentParser(prog="rfms_profiling",description=__version__ + " - " + __description__)
parser.add_argument("-v", "--version", action="version", version=__version__)
profiling = parser.add_argument_group("Generation", "Profiling of genomes")
profiling.add_argument("-f", "--fasta", action="store", help="Input file name for FASTA file", default=False,
                        required=True)
profiling.add_argument("-l", "--level", action="store", help="Level for profiling", default=False,
                        required=True)
profiling.add_argument("-t", "--threshold", action="store", help="complexity log2 threshold", default=1)
profiling.add_argument("-w", "--window", action="store", help="window length for complexity profile", default=1000)
profiling.add_argument("-d", "--drop", action="store", help="drop of bps in complexity profile", default=0)
profiling.add_argument("-op", "--output_profile", action="store", help="Output file name for profile", default=False,
                       required=False)
profiling.add_argument("-s", "--store", action="store", help="Bolean used to indicate of profile should be stored "
                                                             "in ORGs INI", default=False, required=False)
profiling.add_argument("-n", "--org_name", help="Organism name", default="Org", required=False)
args = parser.parse_args()

## Files
profiles_file = "INIs/profiles.ini"
orgs_file = "INIs/orgs.ini"
profiles = configparser.ConfigParser()
profiles.read(profiles_file)
orgs = configparser.ConfigParser()
orgs.read(orgs_file)

## main
def arg_handler():
    fasta = str(args.fasta)
    try:
        level = int(args.level)
    except:
        raise ValueError("Level must be an integral")
    if args.output_profile:
        output_profile = str(args.output_profile) + ".profile"
    else:
        output_profile = args.output_profile
    if args.store:
        store = str(args.store)
    else:
        store = args.store
    org_name = (str(args.org_name))
    window = args.window
    threshold = args.threshold
    drop = args.drop
    return fasta, level, output_profile, store, org_name, window, threshold, drop

if __name__ == '__main__':
    fasta, level, output_profile, store, org_name, window, threshold, drop= arg_handler()
    rep_sizes= profile_report(fasta, level, window, threshold, drop, multiple=True, show=True)