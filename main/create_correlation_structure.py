import os
import numpy as np
import sys
sys.path.append('../scripts/')
sys.path.append('utils/')
from readgx import readgx

gxfile = '/scratch/sbanerj/trans-eqtl/input/gtex/expression/gtex_protein_coding_normalized_lmcorrected_withage_aa.txt'
gx, donors, genes = readgx(gxfile)

C = np.cov(gx.T)
C.shape

from scipy.linalg import eigh
W, Q = eigh(C)
Wsparse = W.copy()
Wsparse[np.where(W < 0)] = 0
correct = np.allclose(C, Q @ np.diag(Wsparse) @ Q.T)
print(correct)

Apart = Q @ np.diag(np.sqrt(Wsparse))
outfile = '/scratch/sbanerj/trans-eqtl/simulation/gtex_correlation_lmcorrected_withage_aa_8core.npy'
np.save(outfile, Apart)
