# Author: Jennifer Nguyen

export PATH=/frazer01/home/jennifer/software/R-4.3.3/bin/:$PATH

date=$1

Rscript /projects/CARDIPS/analysis/epigenome_resource/analyses/jennifer/ipscore_unique_qtls/mashr/scripts/mashr_fit.R --date $date
