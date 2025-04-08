#!/bin/bash

source /scratch/groups/noahr/tami/codis_panel/config.sh

ml R/4.2

for n_ind_rep in {1..100}; do
  for fraction_refs in 0.75; do

    out_dir="$DIR_DATA_PARTITIONS/${fraction_refs}/"
    Rscript 2a_make_ref_and_test_ids.R "$fraction_refs" "$n_ind_rep" "$out_dir"

  done
done
