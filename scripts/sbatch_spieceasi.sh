#! /bin/bash

#SBATCH -A snic2020-15-261
#SBATCH -J spiec.easi
#SBATCH -t 4:00:00
#SBATCH -p core -n 4

module load bioinfo-tools R/4.0.4 R_packages/4.0.4
datadir="/proj/snic2020-16-196/nobackup/CHtrophy_network_nobackup/methanotrophy_networks/processed_tables/"
cd $datadir
Rscript --no-restore --no-save SPIEC-EASI.R

