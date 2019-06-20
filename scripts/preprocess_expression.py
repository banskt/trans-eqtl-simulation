import os
import argparse
import collections
import numpy as np
import scipy.stats as stats
from sklearn import linear_model

def parse_args():

    parser = argparse.ArgumentParser(description = 'Preprocess gene expression')

    parser.add_argument('-g',
                        type = str, 
                        required = True,
                        dest = 'gxfile',
                        metavar = 'FILE',
                        help = 'input gene expression file')

    parser.add_argument('-o',
                        type = str, 
                        required = True,
                        dest = 'outfile',
                        metavar = 'FILE',
                        help = 'output gene expression file')

    parser.add_argument('-c',
                        type = str, 
                        required = True,
                        dest = 'covfile',
                        metavar = 'FILE',
                        help = 'input covariates file')


    parser.add_argument('-m',
                        nargs = '*',
                        default = ['norm'],
                        dest = 'methods',
                        help = 'which method to apply: norm, qn, cclm or any combination')
                        

    opts = parser.parse_args()
    return opts


def normalize_center(Y):
    # requires  G x N matrix, where G is the number of genes, and Y is the number of samples
    # returns   G x N matrix
    newY = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return newY


def read_expression(filename):
    # returns G x N gene expression
    #         list of N donor ids
    #         list of G gene names
    gx = list()
    genes = list()
    with open(filename, 'r') as infile:
        header = next(infile)
        for line in infile:
            linesplit = line.strip().split()
            thisgx = np.array([float(x) for x in linesplit[1:]])
            gx.append(thisgx)
            genes.append(linesplit[0])
    gx = np.array(gx)
    donors = header.strip().split()[1:]
    return gx, donors, genes


def write_expression(GX, donors, genes, filename):
    with open(filename, 'w') as fout:
        header = 'gene_ids\t' + '\t'.join(donors) + '\n'
        fout.write(header)
        for i in range(GX.shape[0]):
            line = '{:s}\t'.format(genes[i]) + '\t'.join(['{:g}'.format(x) for x in GX[i, :]]) + '\n'
            fout.write(line)


def normalize_quantiles(M):
    """
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")  

    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return M


def inverse_quantile_normalize(M):
    """
    After quantile normalization of samples, standardize expression of each gene
    """
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


def gtex_v6_preproc(GX):
    M = normalize_quantiles(GX)
    R = inverse_quantile_normalize(M)
    return R


def covcorrlm(X, W):
    # inputs: X of shape G x N (G = number of genes, N = number of samples)
    #         W of shape C x N (C = number of covariates)
    Wnorm = (W - np.mean(W, axis = 1).reshape(-1, 1)) / np.std(W, axis = 1).reshape(-1, 1)
    linreg = linear_model.LinearRegression()
    linreg.fit(Wnorm.T, X.T)
    Xcorr = X - linreg.predict(Wnorm.T).T
    return Xcorr

def covcorrlasso(X, W, alpha = 0.05):
    # inputs: X of shape G x N (G = number of genes, N = number of samples)
    #         W of shape C x N (C = number of covariates)
    Wnorm = (W - np.mean(W, axis = 1).reshape(-1, 1)) / np.std(W, axis = 1).reshape(-1, 1)
    Xcorr = np.zeros(X.shape)
    #coefs = np.zeros((X.shape[0], W.shape[0]))
    for i in range(X.shape[0]):
        lassoreg = linear_model.Lasso(alpha = alpha)
        lassoreg.fit(Wnorm.T, X[i, :])
        Xcorr[i] = X[i, :] - lassoreg.predict(Wnorm.T).T
    return Xcorr


def pp_options(method, gx, cov = None):
    # a pythonic dictionary runs through all functions and creates the full dictionary
    # hence if-else is used.
    if method == 'norm':
        print('   > Center and scale')
        res = normalize_center(gx)
    elif method == 'qn':
        print('   > QN')
        res = gtex_v6_preproc(gx)
    elif method == 'cclm':
        print('   > Covariate correction with linear model')
        res = covcorrlm(gx, cov)
    elif method == 'cclasso':
        print('   > Adaptive covariate correction using LASSO')
        res = covcorrlasso(gx, cov)
    return res


if __name__ == '__main__':

    opts = parse_args()
    
    # Read gene expression
    print ("Reading gene expression.")
    gx, gxdonors, genes = read_expression(opts.gxfile)

    # Read covariates
    if any(['cc' in x for x in opts.methods]):
        print ("Reading covariates.")
        cov, covdonors, covnames = read_expression(opts.covfile)
        covidx = np.array([covdonors.index(x) for x in gxdonors])
        covsort = cov[:, covidx]
    else:
        cov = None

    # Perform the preprocessing
    gxpp = collections.defaultdict(lambda: None)
    gxpp['raw'] = gx

    for method in opts.methods:
        print(f'Applying {method}')
        msteps = method.split('_')
        prevstep = 'raw'
        curstep = ''

        for mstep in msteps:
            curstep = f'{mstep}' if curstep == '' else f'{curstep}_{mstep}'
            print(f'  Step: {curstep}')
            if gxpp[curstep] is None:
                gxpp[curstep] = pp_options(mstep, gxpp[prevstep], cov = cov)
            prevstep = curstep


    # Write output
    print ("Writing all outputs.")
    oprefix = os.path.splitext(opts.outfile)[0]
    osuffix = os.path.splitext(opts.outfile)[1]
    for method in opts.methods:
        ofilename = f'{oprefix}_{method}{osuffix}'
        write_expression(gxpp[method], gxdonors, genes, ofilename)
