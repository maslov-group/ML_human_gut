#!/bin/sh
for i in {1..100}; do
    qsub -v "var1=$i" qsub_MCMC_add_and_delete_links.pbs
done
