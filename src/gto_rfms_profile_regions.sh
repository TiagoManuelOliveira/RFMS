#!/bin/bash
#
# ==============================================================================
# |                                                                            |
# |             THIS PROGRAM COMPUTES COMPLEXITY PROFILES WITH GTO             |
# |             ==================================================             |
# |                                                                            |
# |   ./gto_rfms_profile_regions.sh sequence.fa level window threshold drop    |
# |                                                                            |
# |                     FILE NEEDED TO THE COMPUTATION:                        |
# |                                                                            |
# |                    $1: FILE WITH THE GENOME IN FASTA                       |
# |                    $2: INTERGER WITH LEVEL                                 |
# |                    $3: INTERGER WITH THE WINDOW LENGTH                     |
# |                    $4: FLOAT WITH LOG THRESHOLD LEVEL                      |
# |                    $5: INTERGER WITH DROP LENGTH                           |
# |                                                                            |
# ==============================================================================
#
# ==============================================================================
# ================================ DEFINITIONS =================================
#
LEVEL="$2";
#
GET_GTO=0;
RUN_COMPARISON=1;
#
# ==============================================================================
# ================================== GET GTO ===================================
#
if [[ "$GET_GTO" -eq "1" ]];
  then
  conda install -c cobilab gto --yes
  fi
#
# ==============================================================================
# =============================== RUN COMPARISON ===============================
#
if [[ "$RUN_COMPARISON" -eq "1" ]];
  then
  gto_fasta_rand_extra_chars < $1 | gto_fasta_to_seq > SEQ;
  gto_reverse < SEQ > SEQ_R;
  #
  # GET WINDOW SIZE BY SEQUENCE SIZE
  LENGTH=`gto_info < SEQ | grep "Number of sym" | awk '{ print $5}'`;
  WINDOW_SIZE="$3"
  echo "WINDOW SIZE: $WINDOW_SIZE";
  #
  gto_geco -v -l $LEVEL -e SEQ
  gto_geco -v -l $LEVEL -e SEQ_R
  #
  gto_upper_bound -u 2 < SEQ.iae   > SEQ_UB
  gto_upper_bound -u 2 < SEQ_R.iae > SEQ_R_UB
  #
  gto_filter -d "$5" -w $WINDOW_SIZE -c < SEQ_UB   > FIL_UB.x
  gto_filter -d "$5" -w $WINDOW_SIZE -c < SEQ_R_UB > FIL_UB_R.x
  #
  tac FIL_UB_R.x > FIL_UB_N.x
  awk '{print $1;}' FIL_UB.x   > IDXES
  awk '{print $2;}' FIL_UB.x   > A_D
  awk '{print $2;}' FIL_UB_N.x > A_R
  #
  gto_min -f A_D -s A_R > A_min
  #
  paste -d '\t' IDXES A_min > PROFILE_N
  #
  gto_segment -t "$4" < PROFILE_N > $1.positions
  #

  fi
#
#
# ==============================================================================
# ==============================================================================
