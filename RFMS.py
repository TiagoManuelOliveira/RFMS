#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__version__ = "0.01"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"

## Imports
import argparse
import os

## Arguments Input

parser = argparse.ArgumentParser(description="Reference-free metagenomic simulator")
parser.add_argument("-s", "--seed", action="store",help="Set simulation seed")
parser.add_argument("-n", "--ngs", action="store", help="Type of NGS technology used")
parser.add_argument("-f", "--fasta", action="store", help="Output file name for FASTA file")
parser.add_argument("-fq", "--fastq", action="store", help="Output file name for FASTQ file")
parser.add_argument("-o", "--org", action="append", help="Organism info -> <org>;<nr de seqs>;<SeqMin:SeqMax>;<sampleType>;<FreqsA:T:G:C>;<freq Reps>")
parser.add_argument("-p", "--profile", action="store", help="Use a predefined profile")

args = parser.parse_args()


def main():
    pass



if __name__ == '__main__':
    main()
    print(args)