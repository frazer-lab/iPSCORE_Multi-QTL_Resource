###########################
This directory contains the scripts and notebooks used to map QTLs across the eight molecular datasets. 

The same linear mixed model (limix) was used to map eQTLs, caQTLs, and haQTLs. The runs differed by:
    - the window size (variants within 1MB of the gene body were tested in eQTLs, and variants within 100kb of each peak were tested for caQTLs and haQTLs)
    - the number of PEER factors (they were optimized for each dataset). 

Fill the config file with the paths and input files required for QTL mapping.

The functions.R file is workhorse for QTL mapping. It contains custom scripts to prepare and run different steps of the pipeline.

The jupyter notebooks were used to write and execute scripts on the cluster.
