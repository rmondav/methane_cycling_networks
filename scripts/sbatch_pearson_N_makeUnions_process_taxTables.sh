#! /bin/bash

#SBATCH -A snic2020-5-529
#SBATCH -J pearsNunion
#SBATCH -t 00:40:00
#SBATCH -p core -n 1

module load bioinfo-tools R/4.0.4 R_packages/4.0.4
maindir="/proj/snic2020-16-196/nobackup/CHtrophy_network_nobackup/methanotrophy_networks"
cd $maindir/
Rscript --no-restore --no-save "scripts/Pearson_correlation.R" &&
Rscript --no-restore --no-save "scripts/make_unioned_tables.R"
Rscript --no-restore --no-save "scripts/process_OTU_tax_tables_for_plotting.R"


