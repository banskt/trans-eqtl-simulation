import numpy as np

def readgx(gxfile):
    exprs = list()
    donors = list()
    genes = list()
    with open(gxfile) as mfile:
        donors = mfile.readline().strip().split("\t")[1:]
        for line in mfile:
            linesplit = line.strip().split("\t")
            gene = linesplit[0].strip()
            genes.append(gene)
            expr = np.array([float(x) for x in linesplit[1:]])
            exprs.append(expr)
    exprs = np.transpose(np.array(exprs))
    return exprs, donors, genes
