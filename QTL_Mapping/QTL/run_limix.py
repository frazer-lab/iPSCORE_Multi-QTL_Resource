from numpy        import ones, stack
from numpy.random import RandomState
from pandas       import DataFrame
from limix.qtl    import scan
from limix.qtl    import iscan

import limix
import numpy  as np
import pandas as pd
import sys
import os

tmp_folder   = sys.argv[1]
element_id   = sys.argv[2]
model        = sys.argv[3]
distribution = sys.argv[4]

infile_data = tmp_folder + "/TMP." + element_id + ".gt.h5"
infile_covs = tmp_folder + "/TMP." + element_id + ".cov.csv"
infile_kin  = tmp_folder + "/TMP." + element_id + ".kin.csv"

data = limix.io.hdf5.read_limix(infile_data)
M    = pd.read_csv(infile_covs, index_col = 0)
K    = pd.read_csv(infile_kin , index_col = 0)

print(M.columns)

Y = data['phenotype']
G = data['genotype' ]

if model == "scan":
    out = scan(G, Y, distribution, K, M=M, verbose=False)
elif model != "scan":
    nrow, ncol = M.shape
    E0         = np.ones((nrow,1))
    E1 = M.loc[:,model]
    out        = iscan(G, Y, distribution, K, M=M, E0=E0, E1=E1, verbose=False)

effsize = pd.DataFrame(out.effsizes['h2'])
pvals   = pd.DataFrame(out.stats)
effsize = effsize[effsize['effect_type']=="candidate"]

if model != "scan":
    effsize = effsize[effsize['env']=="env1_0"]

effsize = effsize[["test", "effsize", "effsize_se"]]
effsize = effsize.set_index('test')

if model == "scan":
    outdf   = effsize.join(pvals[["pv20"]])
elif model != "scan":
    outdf   = effsize.join(pvals[["pv21"]])


outfile = tmp_folder + "/TMP." + element_id + ".qtl.csv"
outdf.to_csv(outfile)

print("Saved: {}".format(outfile))
