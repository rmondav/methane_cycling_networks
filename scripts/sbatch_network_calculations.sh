#! /bin/bash

#SBATCH -A snic2020-5-529
#SBATCH -J netwrk_calcs
#SBATCH -t 06:00:00
#SBATCH -p core -n 4

module load bioinfo-tools R/4.0.4 R_packages/4.0.4
maindir="/proj/snic2020-16-196/nobackup/CHtrophy_network_nobackup/methanotrophy_networks"
cd $maindir/processed_tables/
Rscript --no-restore --no-save "$maindir/scripts/SPIEC-EASI.R" &&
cd ..
Rscript --no-restore --no-save "$maindir/scripts/SPARCC.R" &&
Rscript --no-restore --no-save "scripts/Pearson_correlation.R" &&
Rscript --no-restore --no-save "scripts/make_unioned_tables.R" &&
Rscript --no-restore --no-save "scripts/process_OTU_tax_tables_for_plotting.R"

