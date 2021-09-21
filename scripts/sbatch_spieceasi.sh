#! /bin/bash

#SBATCH -A snic2020-5-529
#SBATCH -J spiec.easi
#SBATCH -t 03:30:00
#SBATCH -p core -n 4

module load bioinfo-tools R/4.0.4 R_packages/4.0.4
maindir="/proj/snic2020-16-196/nobackup/CHtrophy_network_nobackup/methanotrophy_networks"
cd $maindir/processed_tables/
Rscript --no-restore --no-save "$maindir/scripts/SPIEC-EASI.R"
