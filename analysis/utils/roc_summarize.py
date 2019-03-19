import argparse
import numpy as np
import os
import collections
import precision_recall_scores as roc
import load_results


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

    parser.add_argument('--shuffle',
                        action='store_true',
                        dest = 'shuffle',
                        help = 'whether to load results from shuffled analysis')

    parser.add_argument('--which',
                        nargs = '*',
                        type = str,
                        dest = 'whichplots',
                        help = 'which statistics to plot: x and y separated with underscore, e.g. fpr_tpr')

    parser.add_argument('--srcdir',
                        type = str,
                        dest = 'srcdir',
                        help = 'where can you find the simulation directories')

    parser.add_argument('--outdir',
                        type = str,
                        dest = 'outdir',
                        help = 'where to keep the analysis results')

    opts = parser.parse_args()
    return opts


def read_data(key, filename, truth):
    valres = list()
    if key == 'jpa':
        for key, value in load_results.jpa_scores(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'matrixeqtl':
        for key, value in load_results.matrixeqtl(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'rr':
        for key, value in load_results.tejaas(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    return valres


perm_sbeta = '0.01'

opts = parse_args()

simlist = list(range(opts.startsim, opts.endsim))
nsel    = list()
tpr     = list()
ppv     = list()
fpr     = list()
valids  = list()
thres   = list()

rocdata = {'nsel': nsel, 'tpr': tpr, 'ppv': ppv, 'fpr': fpr, 'valids': valids, 'thres': thres}

if not os.path.exists(opts.outdir):
    os.makedirs(opts.outdir)

nsim = 0
for i, simidx in enumerate(simlist):
    simdir = "sim{:03d}".format(simidx)
    sim_trans_file  = os.path.join(opts.srcdir, simdir, "input/expression.trans")
    sim_cis_file    = os.path.join(opts.srcdir, simdir, "input/expression.cis")

    if opts.method == 'jpa':        simoutfile = os.path.join(opts.srcdir, simdir, "tejaas/jpa/all_jpa_pvals.txt")
    if opts.method == 'matrixeqtl': simoutfile = os.path.join(opts.srcdir, simdir, "matrixeqtl/trans_eqtl.txt")
    if opts.method == 'rr':         simoutfile = os.path.join(opts.srcdir, simdir, "tejaas/permnull_sb{:s}/rr.txt".format(perm_sbeta))
    if opts.shuffle:
        if opts.method == 'jpa':        simoutfile = os.path.join(opts.srcdir, simdir, "matrixeqtl_rand/trans_eqtl.txt")
        if opts.method == 'matrixeqtl': simoutfile = os.path.join(opts.srcdir, simdir, "matrixeqtl_rand/trans_eqtl.txt")
        if opts.method == 'rr':         simoutfile = os.path.join(opts.srcdir, simdir, "tejaas_rand/permnull_sb{:s}/rr.txt".format(perm_sbeta))

    if os.path.exists(simoutfile):
        nsim += 1
        true_trans = list()
        with open(sim_trans_file, 'r') as fin:
            for line in (fin):
                linesplit = line.strip().split()
                true_trans.append(linesplit[1])
    
        simdata = read_data(opts.method, simoutfile, true_trans)
        _nsel, _tpr, _ppv, _valids, _thres, _fpr = roc.confusion_matrix(simdata)
        rocdata['nsel'].append(_nsel)
        rocdata['tpr'].append(_tpr)
        rocdata['ppv'].append(_ppv)
        rocdata['fpr'].append(_fpr)
        rocdata['valids'].append(_valids)
        rocdata['thres'].append(_thres)

    else:
        print ("{:s} does not exist".format(simoutfile))

ninterp = int(2 * np.max([len(x) for x in rocdata['thres']]))
print ("Using {:d} simulations for ROC analysis.".format(nsim))
if nsim > 0:
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
        filename = "{:s}_{:s}_{:s}.txt".format(opts.method, xstr, ystr)
        if opts.shuffle: filename = "shuffled_{:s}_{:s}_{:s}.txt".format(opts.method, xstr, ystr)
        outfile = os.path.join(opts.outdir, filename)
        print (outfile)
        with open(outfile, 'w') as fout:
            fout.write("{:s} {:s} std_col2\n".format(xstr, ystr))
            for x, y, e in zip(xvals, yest, yerr):
                fout.write("{:g} {:g} {:g}\n".format(x, y, e))
else:
    print ("No result files. Could not create output.")
