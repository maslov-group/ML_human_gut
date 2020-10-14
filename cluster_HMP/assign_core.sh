#!/bin/sh
for i in {1..100}; do
    sbatch "--export=var1=$i" qsub_MCMC_add_links.sbatch
done
