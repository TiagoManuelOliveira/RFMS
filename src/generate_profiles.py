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
from matplotlib import ticker, gridspec
import tqdm

def get_reps(fasta):
    rep_sizes = []
    with open(str(fasta + ".positions"), "r") as positions_file:
        for line in positions_file:
            if len(line) > 0:
                size = int(line.split(":")[1]) - int(line.split(":")[0])
                if size >= 20:
                    rep_sizes.append(size)

    return  rep_sizes

def get_profile():
    coordinates = list()
    compression = list()
    with open("PROFILE_N", "r") as compression_result:
        for line in compression_result:
            if len(line) > 0:
                coordinates.append(int(line.split("\t")[0]))
                compression.append(float(line.split("\t")[1]))

    return coordinates, compression

def profile_draw_graphs(ax1, ax2, rep_sizes, hist_max, name, threshold, show=False, nested=False, bar_space=None, minor_size=4):
    # Get compression profile information
    coordinates, compression = get_profile()
    x_scale = 10 ** (len(str(len(coordinates))) - 1)
    x_limit = (int(len(coordinates) / x_scale) + 1) * x_scale

    # Compression profile
    ax1.plot(coordinates, compression, linewidth=0.5, color="b")
    ax1.plot([0, x_limit], [threshold, threshold], 'k-', lw=0.5, label="_not in legend", color="g")
    ax1.set_xlim([0,x_limit])
    ax1.set_ylim([0, 2])
    ax1.grid(True, alpha=0.5)
    ax1.set_ylabel("Bps")
    ax1.set_xlabel=("Length (Bps)")
    ax1.set_title("Low Complexity Regions")

    # Histogram
    ax2.hist(rep_sizes, bins=len(rep_sizes), align="mid", rwidth=bar_space, edgecolor='k')
    ax2.set_ylabel("Frequency of Regions Length")
    ax2.set_xlabel("Length (Bps)")
    ax2.set_title("Low Complexity Regions Sizes")
    ax2.grid(axis='y', alpha=0.5)
    ax2.set_ylim([0,hist_max+1])
    ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.tick_params(which='minor', length=1, labelsize=minor_size)
    ax2.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%d'))

    if not nested:
        plt.savefig(name)
        if show:
            plt.show()



def profile_fasta(fasta, level, window, threshold, drop):
    command = ["sh", "src/gto_rfms_profile_regions.sh", fasta, str(level), str(window), str(threshold), str(drop)]
    subprocess.call(command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

def profile_report(fasta, level, window, threshold, drop, multiple=False, show=False):
    cleanup.create_folder("graphs")
    fasta_name = str(fasta.split(".")[0] +
                     "_l" + str(level) +
                     "_w" + str(window) +
                     "_t" + str(threshold) +
                     "_d" + str(drop))

    profile_graph = str("graphs/" + fasta_name) + ".pdf"
    profile_report = str("graphs/" + fasta_name) + "_report.pdf"

    if multiple:
        bar = tqdm.tqdm(total=5)
        f = plt.figure(figsize=(8.27, 11.69))
        plt.rc('axes', labelsize=6)
        plt.rc('axes', titlesize=6)  # fontsize of the axes title
        plt.rc('xtick', labelsize=4)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=4)  # fontsize of the tick labels
        main_grid = gridspec.GridSpec(15,7, wspace=0.4, hspace=1)
        #########
        grid1 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4,:3])
        profile_fasta(fasta, level, int(int(window)/4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        ax1 = plt.subplot(grid1[:2, 0:])
        ax2 = plt.subplot(grid1[2:, 0:])
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_report, float(threshold), False, True, 0.5, 2)
        bar.update()
        #####
        grid2 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4, 4:])
        profile_fasta(fasta, level, int(int(window) / 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        ax3 = plt.subplot(grid2[:2, 0:])
        ax4 = plt.subplot(grid2[2:, 0:])
        profile_draw_graphs(ax3, ax4, rep_sizes, hist.max(), profile_report, float(threshold), False, True, 0.5, 2)
        bar.update()
        #####
        grid3 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[5:9, 2:5])
        profile_fasta(fasta, level, int(window), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        ax5 = plt.subplot(grid3[:2, 0:])
        ax6 = plt.subplot(grid3[2:, 0:])
        profile_draw_graphs(ax5, ax6, rep_sizes, hist.max(), profile_report, float(threshold), False, True, 0.5, 2)
        bar.update()
        #####
        grid4 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, :3])
        profile_fasta(fasta, level, int(int(window) * 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        ax7 = plt.subplot(grid4[:2, 0:])
        ax8 = plt.subplot(grid4[2:, 0:])
        profile_draw_graphs(ax7, ax8, rep_sizes, hist.max(), profile_report, float(threshold), False, True, 0.5, 2)
        bar.update()
        #####
        grid5 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, 4:])
        profile_fasta(fasta, level, int(int(window) * 4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        ax9 = plt.subplot(grid5[:2, 0:])
        ax10 = plt.subplot(grid5[2:, 0:])
        profile_draw_graphs(ax9, ax10, rep_sizes, hist.max(), profile_report, float(threshold), multiple, show, 0.5, 2)
        bar.update()
    else:
        profile_fasta(fasta, level, window, threshold, drop) # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)

        grid = plt.GridSpec(4, 3, wspace=0.4, hspace=1)
        ax1 = plt.subplot(grid[:2, 0:])
        ax2 = plt.subplot(grid[2:, 0:])
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_graph, float(threshold), show, multiple, show)

    cleanup.clean_profiling() # Cleanup