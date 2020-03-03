#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Tiago Oliveira"
__copyright__ = ""
__credits__ = ["Tiago Oliveira"]
__license__ = "GPL-3.0"
__version__ = "RFMS - 0.1b"
__maintainer__ = "Tiago Oliveira"
__email__ = "tiagomanuel28@gmail.com"
__status__ = "Production"
__description__ = "Reference-free metagenomic simulator"


import os, subprocess
import numpy as np
import matplotlib.pyplot as plt
import src.file_cleanup as cleanup


def profile_fasta(fasta, level, window, threshold, drop):
    cleanup.create_folder("graphs")
    fasta_name = str(fasta.split(".")[0] +
                     "_l" + str(level) +
                     "_w" + str(window) +
                     "_t" + str(threshold) +
                     "_d" + str(drop))

    profile_graph = str("graphs/" + fasta_name)
    profile_histogram = str("graphs/"+ fasta_name) + "_hist.pdf"
    command = ["sh", "src/gto_rfms_profile_regions.sh", fasta, str(level), str(window), str(threshold), str(drop)]
    subprocess.call(command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    command = ["sh", "src/gto_rfms_profile_graph.sh", profile_graph]
    subprocess.call(command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    rep_sizes = []
    with open(str(fasta+".positions"), "r") as positions_file:
        for line in positions_file:
            if len(line) > 0:
                rep_sizes.append(int(line.split(":")[1])-int(line.split(":")[0]))

    cleanup.clean_profiling()

    rep_sizes = np.array(rep_sizes)

    print(rep_sizes)
    hist, bin_edges = np.histogram(rep_sizes, bins=int(len(rep_sizes)))
    print(hist, bin_edges)
    plt.hist(rep_sizes, bins=len(rep_sizes))
    plt.savefig(profile_histogram)
    plt.show()


    return rep_sizes
