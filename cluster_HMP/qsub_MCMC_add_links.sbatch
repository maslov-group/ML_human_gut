#!/bin/bash
#
#SBATCH -t 04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J myjob
#SBATCH --partition=secondary
#SBATCH --get-user-env
#
# End of embedded SBATCH options
# set directory to current one
#new_folder=$PBS_O_WORKDIR
new_folder=$SLURM_SUBMIT_DIR
cd $new_folder

. /etc/profile.d/modules.sh # add this line to make sure that the slurm script aware of the modules system yet
module load anaconda/2019-Oct/3
# Run the executable (a.out)
python ./MCMC_add_links.py $var1
