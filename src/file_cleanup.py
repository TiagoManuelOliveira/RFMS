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

import os

def create_folder(folder):
    if os.path.exists(folder) == False:
        print("Folder ", folder, " doesnt exist...\n>", folder, " folder created")
        os.makedirs(folder)

def clean_profiling():
    # files created during profiling process
    files = ["A_D",
             "A_min",
             "A_R",
             "FIL_UB.x",
             "FIL_UB_N.x",
             "FIL_UB_R.x",
             "IDXES",
             "PROFILE_N",
             "SEQ",
             "SEQ.co",
             "SEQ.iae",
             "SEQ_R",
             "SEQ_R.co",
             "SEQ_r.iae",
             "SEQ_R.iae",
             "SEQ_R_UB",
             "SEQ_UB"]
    # find .positions files
    for file in os.listdir():
        if file[-10:] == ".positions":
            files.append(file)
    for file in files:
        if file in os.listdir():
            os.remove(file)

if __name__ == '__main__':
    clean_profiling()