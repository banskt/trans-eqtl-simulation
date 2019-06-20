import argparse
import numpy as np
import os
import collections
import precision_recall_scores as roc
import load_results
from statsmodels.stats import multitest


INFO_FIELDS = ['rsid', 'stat', 'causality']
class ValidateResult(collections.namedtuple('_ValidateResult', INFO_FIELDS)):
    __slots__ = ()


def parse_args():

    parser = argparse.ArgumentParser(description="Calculate ROC plot values for a set of simulations")

    parser.add_argument('--startsim',
                        type = int,
                        dest = 'startsim',
                        help = 'Index of the starting simulation')

    parser.add_argument('--endsim',
                        type = int,
                        dest = 'endsim',
                        help = 'Index of the last simulation')

    parser.add_argument('--method',
                        type = str,
                        dest = 'method',
                        help = 'Name of the method to analyze')

    parser.add_argument('--which',
                        nargs = '*',
                        type = str,
                        dest = 'whichplots',
                        help = 'which statistics to plot: x and y separated with underscore, e.g. fpr_tpr')

    parser.add_argument('--ntop',
                        nargs = '*',
                        type = int,
                        dest = 'ntop',
                        help = 'number of SNPs to consider for discovering trans-eQTLs, accepts list of integers')

    parser.add_argument('--srcdir',
                        type = str,
                        dest = 'srcdir',
                        help = 'where can you find the simulation directories')

    parser.add_argument('--mdir',
                        type = str,
                        dest = 'method_dir',
                        help = 'where can you find this method\'s output')

    parser.add_argument('--outdir',
                        type = str,
                        dest = 'outdir',
                        help = 'where to keep the analysis results')

    parser.add_argument('--outprefix',
                        type = str,
                        dest = 'outprefix',
                        help = 'prefix of the output file')

    opts = parser.parse_args()
    return opts


def read_data(key, filename, truth):
    valres = list()
    fdr = list()
    if key == 'jpa':
        for key, value in load_results.jpa_scores(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'matrixeqtl':
        tmpres, fdr = load_results.matrixeqtl_with_fdr(filename)
        for key, value in tmpres.items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'rr':
        for key, value in load_results.tejaas(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    return valres, fdr


def get_power(valres, simfdr, method, fdr):
    if method == 'jpa' or method == 'rr':
        keepres = [x for x in valres if not np.isnan(x.stat)]
        nitems = len(keepres)
        pvals = np.array([np.power(10.0, -x.stat) for x in keepres])
        true  = np.array([x.causality for x in keepres])
        pred, pcorr = multitest.fdrcorrection(pvals, alpha = fdr)
        truepred = np.array([True if x and true[i] else False for i, x in enumerate(pred)])
        power = sum(truepred) / sum(true)
    if method == 'matrixeqtl':
        pcorr = np.array(simfdr)
        true = np.array([x.causality for x in valres])
        truepred = np.array([True if x <= fdr and true[i] else False for i, x in enumerate(pcorr)])
        power = sum(truepred) / sum(true)
    return power


def get_validpred(valres, simfdr, method, fdr, ntop):
    nselect = ntop
    if method == 'jpa' or method == 'rr':
        keepres = [x for x in valres if not np.isnan(x.stat)]
        pvals = np.array([np.power(10.0, -x.stat) for x in keepres])
        _, pcorr = multitest.fdrcorrection(pvals, alpha = fdr)
    if method == 'matrixeqtl':
        keepres = valres
        pcorr = np.array(simfdr)

    isort = np.argsort(pcorr)
    if pcorr[isort[ntop - 1]] > fdr:
        if pcorr[isort][0] >= fdr:
            nselect = 0
        else:
            nselect = np.where(pcorr[isort] < fdr)[0][-1] + 1
    iselect = isort[:nselect]
    predselect = [keepres[i].causality for i in iselect]
    truepred = sum(predselect)
    falsepred = nselect - truepred
    return truepred, falsepred

opts = parse_args()

simlist = list(range(opts.startsim, opts.endsim))
nsel    = list()
tpr     = list()
ppv     = list()
fpr     = list()
valids  = list()
thres   = list()

rocdata = {'nsel': nsel, 'tpr': tpr, 'ppv': ppv, 'fpr': fpr, 'valids': valids, 'thres': thres}
fdrlist = [0.002, 0.05, 0.10, 0.20, 0.50]
power = dict()
truepos = dict()
falsepos = dict()

for fdr in fdrlist:
    power[fdr] = list()
    truepos[fdr] = dict()
    falsepos[fdr] = dict()
    for ntop in opts.ntop:
        truepos[fdr][ntop] = list()
        falsepos[fdr][ntop] = list()


if not os.path.exists(opts.outdir):
    os.makedirs(opts.outdir)

outfileprefix = os.path.join(opts.outdir, opts.outprefix)

nsim = 0
for i, simidx in enumerate(simlist):
    simdir = "sim{:03d}".format(simidx)
    sim_trans_file  = os.path.join(opts.srcdir, simdir, "input/expression.trans")
    sim_cis_file    = os.path.join(opts.srcdir, simdir, "input/expression.cis")

    simoutmdir = os.path.join(opts.srcdir, simdir, opts.method_dir)
    if opts.method == 'jpa':        simoutfile = os.path.join(simoutmdir, "all_jpa_pvals.txt")
    if opts.method == 'rr':         simoutfile = os.path.join(simoutmdir, "rr.txt")
    if opts.method == 'matrixeqtl': simoutfile = os.path.join(simoutmdir, "trans_eqtl.txt")

    if os.path.exists(simoutfile):
        nsim += 1
        true_trans = list()
        with open(sim_trans_file, 'r') as fin:
            for line in (fin):
                linesplit = line.strip().split()
                true_trans.append(linesplit[1])
    
        simdata, simfdr = read_data(opts.method, simoutfile, true_trans)
        _nsel, _tpr, _ppv, _valids, _thres, _fpr = roc.confusion_matrix(simdata)
        rocdata['nsel'].append(_nsel)
        rocdata['tpr'].append(_tpr)
        rocdata['ppv'].append(_ppv)
        rocdata['fpr'].append(_fpr)
        rocdata['valids'].append(_valids)
        rocdata['thres'].append(_thres)
        for fdr in fdrlist:
            power[fdr].append(get_power(simdata, simfdr, opts.method, fdr))
            for ntop in opts.ntop:
                thistp, thisfp = get_validpred(simdata, simfdr, opts.method, fdr, ntop)
                truepos[fdr][ntop].append(thistp)
                falsepos[fdr][ntop].append(thisfp)
    #else:
    #    print ("{:s} does not exist".format(simoutfile))

print ("{:s} --> {:d} simulations".format(opts.outprefix, nsim))
if nsim > 0:
    ninterp = int(2 * np.max([len(x) for x in rocdata['thres']]))
    for desc in opts.whichplots:
        xstr, ystr = tuple(desc.split('_'))
        xx = rocdata[xstr]
        yy = rocdata[ystr]
        xmax = int(np.max([max(x) for x in xx]) + 0.5)
        xvals = np.linspace(0, xmax, ninterp)
        reverse = False
        if xx[0][-1] < xx[0][0]: reverse = True
        if not reverse:
            yest_arr = np.array([np.interp(xvals, xx[i], yy[i]) for i in range(nsim)])
        else:
            yest_arr = np.array([np.interp(xvals, xx[i][::-1], yy[i][::-1]) for i in range(nsim)])
        yest = np.mean(yest_arr, axis = 0)
        yerr = np.std (yest_arr, axis = 0)
        outfile = "{:s}_{:s}_{:s}.txt".format(outfileprefix, xstr, ystr)
        with open(outfile, 'w') as fout:
            fout.write("{:s} {:s} {:s}_std\n".format(xstr, ystr, ystr))
            for x, y, e in zip(xvals, yest, yerr):
                fout.write("{:g} {:g} {:g}\n".format(x, y, e))
    # Power at different FDR
    outfile = "{:s}_power.txt".format(outfileprefix)
    with open(outfile, 'w') as fout:
        fout.write("FDR Power Error\n")
        for fdr in fdrlist:
            powermean = np.mean(power[fdr])
            powerstd  = np.std(power[fdr])
            fout.write("{:g} {:g} {:g}\n".format(fdr, powermean, powerstd))

    #TP, FP at different FDR
    outfile = "{:s}_tpfp.txt".format(outfileprefix)
    with open(outfile, 'w') as fout:
        fout.write("Ntop FDR TPmean FPmean TPstd FPstd\n")
        for fdr in fdrlist:
            for ntop in opts.ntop:
                tpmean = np.mean(truepos[fdr][ntop])
                tpstd  = np.std(truepos[fdr][ntop])
                fpmean = np.mean(falsepos[fdr][ntop])
                fpstd  = np.std(falsepos[fdr][ntop])
                fout.write("{:g} {:d} {:g} {:g} {:g} {:g}\n".format(fdr, ntop, tpmean, fpmean, tpstd, fpstd))

else:
    print ("No result files. Could not create output.")
