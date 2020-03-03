#!/bin/bash
#
# ==============================================================================
# |                                                                            |
# |             THIS PROGRAM DRAWS COMPLEXITY PROFILES GRAPHS                  |
# |             ==================================================             |
# |                                                                            |
# |                      ./gto_rfms_profile_graph.sh name.pdf                  |
# |                                                                            |
# |                     FILE NEEDED TO THE COMPUTATION:                        |
# |                                                                            |
# |                    $1: NAME OF OUTPUT PDF FILE                             |
# |                                                                            |
# ==============================================================================
#
# ==============================================================================
# =================================== GRAPH ====================================
#
PDF="${1}.pdf"
gnuplot << EOF
    reset
    set terminal pdfcairo enhanced color font 'Verdana,12'
    set output "$PDF"
    set style line 101 lc rgb '#000000' lt 1 lw 4
    set border 3 front ls 101
    set tics nomirror out scale 0.75
    set format '%g'
    set size ratio 0.1
    unset key
    set yrange [:2]
    set xrange [:]
    set xtics auto
    set ytics 0.5
    set grid
    set ylabel "Bps"
    set xlabel "Length"
    set border linewidth 1.5
    set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5 ps 0.4 # --- blue
    set style line 2 lc rgb '#0060ad' lt 1 lw 4 pt 6 ps 0.4 # --- green
    plot "PROFILE_N" using 1:2 with lines ls 1
EOF