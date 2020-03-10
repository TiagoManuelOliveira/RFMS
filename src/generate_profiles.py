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


import subprocess
import numpy as np
import matplotlib.pyplot as plt
import src.file_cleanup as cleanup
from matplotlib import ticker, gridspec
import tqdm

def get_reps(fasta):
    '''
    Gets low complexity regions from .positions file, calculates de length and returns a list with the values
    :param fasta: name of the FASTA file profiled for complex regions
    :return: list with length values
    '''
    rep_sizes = []
    with open(str(fasta + ".positions"), "r") as positions_file:
        for line in positions_file:
            if len(line) > 0:
                size = int(line.split(":")[1]) - int(line.split(":")[0])
                if size >= 20:
                    rep_sizes.append(size)

    return rep_sizes

def get_profile():
    '''
    Gets the pair of coordinates and compression value of profiling low complexity regions
    :return: tuple with lists, with coordinates and compression value
    '''
    coordinates = list()
    compression = list()
    with open("PROFILE_N", "r") as compression_result:
        for line in compression_result:
            if len(line) > 0:
                coordinates.append(int(line.split("\t")[0]))
                compression.append(float(line.split("\t")[1]))

    return coordinates, compression

def profile_draw_graphs(ax1, ax2, rep_sizes, hist_max, name, threshold, nested=False, bar_space=None, minor_size=4):
    '''
    Draws subplots for profile reporting, given 2 subplots (ax1 and ax2)
    :param ax1: Subplot from a GridSpec Object. Used to draw compression profile plot
    :param ax2: Subplot from a GridSpec Object. Used to draw histogram with the length frequency of low
    complexity regions
    :param rep_sizes: list with lengths low complexity regions
    :param hist_max: int with maximum value of y axis to be used in ax2
    :param name: str with name of the file to be saved in .pdf format
    :param threshold: float with threshold used during profiling, to be drawn as a horizontal line ax1
    :param nested: Bol with information if the GridSpec object is nested
    :param bar_space: float with space to be used for separation of bars in histogram (ax2)
    :param minor_size: int with font size of minor label in x axis of histogram (ax2)
    '''
    # Get compression profile information
    coordinates, compression = get_profile()
    x_scale = 10 ** (len(str(len(coordinates))) - 1)
    x_limit = (int(len(coordinates) / x_scale) + 1) * x_scale

    # Compression profile (plot)
    ax1.plot(coordinates, compression, linewidth=0.5, color="b")
    ax1.plot([0, x_limit], [threshold, threshold], 'k-', lw=0.5, label="_not in legend", color="g")
    ax1.set_xlim([0,x_limit])
    ax1.set_ylim([0, 2])
    ax1.grid(True, alpha=0.5)
    ax1.set_ylabel("Bps")
    ax1.set_xlabel=("Length (Bps)")
    ax1.set_title("Low Complexity Regions")

    # Histogram
    if len(rep_sizes) == 0:
        ax2.hist(rep_sizes, bins=1, align="mid", rwidth=bar_space, edgecolor='k')
    else:
        ax2.hist(rep_sizes, bins=len(rep_sizes), align="mid", rwidth=bar_space, edgecolor='k')
    ax2.set_ylabel("Frequency of Regions Length")
    ax2.set_xlabel("Length (Bps)")
    ax2.set_title("Low Complexity Regions Sizes")
    ax2.grid(axis='y', alpha=0.5)
    ax2.set_ylim([0,hist_max+1])
    ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax2.tick_params(which='minor', length=1, labelsize=minor_size)
    ax2.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%d'))

    # Save figure in .pdf format
    if not nested:
        plt.savefig(name)



def profile_fasta(fasta, level, window, threshold, drop):
    '''
    Runs modified gto_profile_regions.sh for rfms
    :param fasta: str with name of the fasta file to profiled
    :param level: int with level to be used by gto in profiling
    :param window: int with length of the window to be used in profiling
    :param threshold: float with threshold to be used to acess compression (value between 0 and 2)
    :param drop: int with drop (sensibility of compression)
    :return:
    '''
    #TODO: NEED TO ADD TRY METHOD IN CASE PROFILING FAILS DURING .SH RUN
    command = ["sh", "src/gto_rfms_profile_regions.sh", fasta, str(level), str(window), str(threshold), str(drop)]
    subprocess.call(command, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

def profile_report(fasta, level, window, threshold, drop, multiple=False, show=False):
    '''
    Functions generates a report for the profiling of low complexity regions in a genomic sequence using gto
    compression tools
    :param fasta: str with name of the .FASTA file to profiled
    :param level: int with level to be used internally in GTO for compression profiling
    :param window: int with length of window to be used for compression profiling
    :param threshold: float with threshold for compression profiling (between 0.0 and 2.0)
    :param drop: int with drop of compression profiling (sensibility of compression)
    :param multiple: bol with information if a multiple window report should be done (default=False)
    :param show: bol with information if the graph should be displayed in a window after completion (default=False)
    :return:
    '''
    cleanup.create_folder("graphs")
    fasta_name = str(fasta.split(".")[0] +
                     "_l" + str(level) +
                     "_w" + str(window) +
                     "_t" + str(threshold) +
                     "_d" + str(drop))

    profile_graph = str("graphs/" + fasta_name) + ".pdf"
    profile_report = str("graphs/" + fasta_name) + "_report.pdf"

    if multiple: # Multiple window report
        # Setting up tqdm progress bar
        bar = tqdm.tqdm(total=5)
        # Initialization of figure with size of an A4 page
        f = plt.figure(figsize=(8.27, 11.69))
        plt.rc('axes', labelsize=6)
        plt.rc('axes', titlesize=6)  # fontsize of the axes title
        plt.rc('xtick', labelsize=4)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=4)  # fontsize of the tick labels
        # Creation of a main GridSpec object
        main_grid = gridspec.GridSpec(15,7, wspace=0.4, hspace=1.5)

        ## Report for window / 4
        # Nesting of GridSpec object
        grid1 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4,:3])
        #Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) / 4))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window)/4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax1 = plt.subplot(grid1[:2, 0:])
        ax2 = plt.subplot(grid1[2:, 0:])
        ax1.text(0.5, 1.3, str("Window: " + str(int(window)/4)), size=10, ha="center", transform=ax1.transAxes)
        #Drawing of subplots in figure
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_report, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window / 2
        # Nesting of GridSpec object
        grid2 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[:4, 4:])
        # Handling of progress bars
        description = str("Generating report for Window: "+ str(int(window) / 2))
        bar.set_description(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) / 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax3 = plt.subplot(grid2[:2, 0:])
        ax4 = plt.subplot(grid2[2:, 0:])
        ax3.text(0.5, 1.3, str("Window: " + str(int(window)/2)), size=10, ha="center", transform=ax3.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax3, ax4, rep_sizes, hist.max(), profile_report, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window
        # Nesting of GridSpec object
        grid3 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[5:9, 2:5])
        # Handling of progress bars
        description = str("Generating report for Window: "+ str(int(window)))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(window), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax5 = plt.subplot(grid3[:2, 0:])
        ax6 = plt.subplot(grid3[2:, 0:])
        ax5.text(0.5, 1.3, str("Window: " + str(int(window))), size=10, ha="center", transform=ax5.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax5, ax6, rep_sizes, hist.max(), profile_report, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window * 2
        # Nesting of GridSpec object
        grid4 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, :3])
        # Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) * 2))
        bar.set_description_str(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) * 2), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax7 = plt.subplot(grid4[:2, 0:])
        ax8 = plt.subplot(grid4[2:, 0:])
        ax7.text(0.5, 1.3, str("Window: " + str(int(window)*2)), size=10, ha="center", transform=ax7.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax7, ax8, rep_sizes, hist.max(), profile_report, float(threshold), True, 0.5, 2)
        # Progress update
        bar.update()

        ## Report for window * 4
        # Nesting of GridSpec object
        grid5 = gridspec.GridSpecFromSubplotSpec(4, 3, wspace=0.4, hspace=1, subplot_spec=main_grid[10:, 4:])
        # Handling of progress bars
        description = str("Generating report for Window: " + str(int(window) * 4))
        bar.set_description(desc=description)
        # Profiling
        profile_fasta(fasta, level, int(int(window) * 4), threshold, drop)  # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Creation of Subplots in nested GridSpec
        ax9 = plt.subplot(grid5[:2, 0:])
        ax10 = plt.subplot(grid5[2:, 0:])
        ax9.text(0.5, 1.3, str("Window: " + str(int(window)*4)), size=10, ha="center", transform=ax9.transAxes)
        # Drawing of subplots in figure
        profile_draw_graphs(ax9, ax10, rep_sizes, hist.max(), profile_report, float(threshold), False, 0.5, 2)
        # Handling of progress bars
        bar.update()
        description = "Finished compiling profile report"
        bar.set_description(desc=description)

    else: #If not nested
        # Setting up tqdm progress bar
        bar = tqdm.tqdm(total=2, desc="Generating complexity profile")
        # Profiling
        profile_fasta(fasta, level, window, threshold, drop) # Run profile pipeline
        rep_sizes = get_reps(fasta)
        np_reps = np.array(rep_sizes)
        bin_size = int(len(rep_sizes))
        if bin_size<=0:
            bin_size = 1
        hist, bin_edges = np.histogram(np_reps, bins=bin_size)
        # Handling of progress bars
        bar.set_description_str(desc="Generating profile report")
        bar.update()
        # Creation of GridSpec object
        grid = plt.GridSpec(4, 3, wspace=0.4, hspace=2)
        #title="Complexity profile of " + str(fasta) + \
        #      "\nlevel: " + str(level) + \
        #      ", window: " + str(window) + \
        #      ", threshold" + str(threshold) + \
        #      ", drop: " + str(drop)
        #plt.suptitle(title, fontsize=12)
        # Creation of Subplots
        ax1 = plt.subplot(grid[:2, 0:])
        ax2 = plt.subplot(grid[2:, 0:])
        # Drawing of subplots in figure
        profile_draw_graphs(ax1, ax2, rep_sizes, hist.max(), profile_graph, float(threshold), multiple, show)
        # Handling of progress bars
        bar.set_description_str(desc="Finished compiling profile report")
        bar.update()

    cleanup.clean_profiling() # Cleanup of files produced during profiling
    if show: # Show figure in a seperate window
        plt.show()