# Author: Jennifer Nguyen

tissue=$1
qtl_list=$2
date=$3
Rscript /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/scripts/get_qtl_stats_for_randvars.R --tissue $tissue --qtl_list $qtl_list --date $date