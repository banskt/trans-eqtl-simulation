import numpy as np
import collections
from scipy import stats

def plist2dict(rsids, pvals):
    pvals = np.array(pvals)
    min_nonzero = np.min(pvals[np.nonzero(pvals)])
    pvals[pvals == 0] = min_nonzero
    res = dict()
    for i, key in enumerate(rsids):
        res[key] = -np.log10(pvals[i])
    return res

def scorelist2dict(rsids, scores):
    scores = np.array(scores)
    minscore = np.min(scores)
    res = dict()
    for i, key in enumerate(rsids):
        res[key] = scores[i] - minscore
    return res

def tejaas(filepath):
    rsids = list()
    qscales = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr   = line.strip().split("\t")
            rsid  = arr[0]
            p     = float(arr[5])
            q     = float(arr[2])
            mu    = float(arr[3])
            sigma = float(arr[4])
            qscale = (q - mu) / sigma
            rsids.append(rsid)
            qscales.append(qscale)
    pvalarr = 2.0 * (1 - stats.norm.cdf(np.abs(qscales)))
    pvals = list(pvalarr)
    res = plist2dict(rsids, pvals)
    return res

def jpa(filepath):
    pvals = list()
    rsids = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr = line.strip().split("\t")
            rsids.append(arr[0])
            pvals.append(float(arr[2]))
    res = plist2dict(rsids, pvals)
    return res


def jpa_scores(filepath):
    qscores = list()
    rsids = list()
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr = line.strip().split("\t")
            rsids.append(arr[0])
            qscores.append(float(arr[1]))
    res = scorelist2dict(rsids, qscores)
    return res


def matrixeqtl_with_fdr(filepath):
    rescoll = collections.defaultdict(list)
    fdrcoll = collections.defaultdict(list)
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            pval = float(arr[4])
            fdr = float(arr[5])
            rescoll[rsid].append(pval)
            fdrcoll[rsid].append(fdr)
    pvals = list()
    rsids = list()
    fdrs  = list()
    for key, val in rescoll.items():
        rsids.append(key)
        pvals.append(np.min(val))
        fdrs.append(np.min(fdrcoll[key]))
    res = plist2dict(rsids, pvals)
    return res, fdrs


def matrixeqtl(filepath):
    res, fdrs = matrixeqtl_with_fdr(filepath)
    return res


def matrixeqtl_fdr(filepath):
    rescoll = collections.defaultdict(list)
    with open(filepath, 'r') as mfile:
        next(mfile)
        for line in mfile:
            arr  = line.strip().split("\t")
            rsid = arr[0]
            pval = float(arr[5])
            rescoll[rsid].append(pval)
    pvals = list()
    rsids = list()
    for key, val in rescoll.items():
        rsids.append(key)
        pvals.append(np.min(val))
    res = plist2dict(rsids, pvals)
    return res


##def matrixeqtl(filepath):
##    pvals = list()
##    rsids = list()
##    with open(filepath, 'r') as mfile:
##        next(mfile)
##        for line in mfile:
##            arr  = line.strip().split("\t")
##            rsid = arr[0]
##            pval = float(arr[4])
##            if rsid not in rsids:
##                rsids.append(rsid)
##                pvals.append(pval)
##    res = plist2dict(rsids, pvals)
##    return res

##def matrixeqtl_fdr(filepath):
##    pvals = list()
##    rsids = list()
##    with open(filepath, 'r') as mfile:
##        next(mfile)
##        for line in mfile:
##            arr  = line.strip().split("\t")
##            rsid = arr[0]
##            fdr  = float(arr[5])
##            if rsid not in rsids:
##                rsids.append(rsid)
##                pvals.append(fdr)
##    res = plist2dict(rsids, pvals)
##    return res


##def matrixeqtl_old(filepath):
##    res = dict()
##    with open(filepath, 'r') as mfile:
##        next(mfile)
##        for line in mfile:
##            arr  = line.strip().split("\t")
##            rsid = arr[0]
##            pval = float(arr[4])
##            if rsid not in res:
##                res[rsid] = -np.log10(pval)
##    
##    return res

##def matrixeqtl_fdr_old(filepath):
##    res = dict()
##    with open(filepath, 'r') as mfile:
##        next(mfile)
##        for line in mfile:
##            arr  = line.strip().split("\t")
##            rsid = arr[0]
##            fdr  = float(arr[5])
##            if rsid not in res:
##                res[rsid] = -np.log10(fdr)
##    return res
