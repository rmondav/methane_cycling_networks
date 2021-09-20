#! /bin/bash

#SBATCH -A snic2020-5-529
#SBATCH -J CH4cyc_netwks
#SBATCH -t 20:30:00
#SBATCH -p core -n 4

module load bioinfo-tools R/4.0.4 R_packages/4.0.4
maindir="/proj/snic2020-16-196/nobackup/CHtrophy_network_nobackup/methanotrophy_networks"
cd $maindir/
Rscript --no-restore --no-save "$maindir/scripts/CH4cycling_networks.R"

