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
                        nargs = *,
                        type = str,
                        help = 'which statistics to plot: x and y separated with underscore, e.g. fpr_tpr')

    parser.add_argument('--srcdir',
                        type = str,
                        dest = 'srcdir',
                        help = 'where can you find the simulation directories')

    parser.add_argument('--outdir',
                        type = str,
                        detst = 'outdir',
                        help = 'where to keep the analysis results')

    opts = parser.parse_args()
    return opts


def read_data(key, filename, truth):
    valres = list()
    if key == 'jpa':
        for key, value in load_results.jpa_scores(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'meqtl':
        for key, value in load_results.matrixeqtl_fdr(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    if key == 'rr':
        for key, value in load_results.tejaas(filename).items():
            caus = 1 if key in truth else 0
            valres.append(ValidateResult(rsid = key, stat = value, causality = caus))
    return valres


simlist = list(range(opts.startsim, opts.endsim))
simdata = [None for x in simlist]
nsel    = [None for x in simlist]
tpr     = [None for x in simlist]
ppv     = [None for x in simlist]
fpr     = [None for x in simlist]
valids  = [None for x in simlist]
thres   = [None for x in simlist]

rocdata = {'nsel': nsel, 'tpr': tpr, 'ppv': ppv, 'fpr': fpr, 'valids': valids, 'thres': thres}

if not os.path.exists(opts.outdir):
    os.makedirs(opts.outdir)

for i, simidx in enumerate(simlist):
    simdir = "sim{:03d}".format(simidx)
    sim_trans_file  = os.path.join(opts.srcdir, simdir, "input/expression.trans")
    sim_cis_file    = os.path.join(opts.srcdir, simdir, "input/expression.cis")
    rr_outfile      = os.path.join(opts.srcdir, simdir, "tejaas/permnull_sb{:s}/rr.txt".format(perm_sbeta))
    meqtl_outfile   = os.path.join(opts.srcdir, simdir, "matrixeqtl/trans_eqtl.txt")
    jpa_outfile     = os.path.join(opts.srcdir, simdir, "tejaas/jpa/all_jpa_pvals.txt")
    if opts.shuffle:
        meqtl_outfile  = os.path.join(opts.srcdir, simdir, "matrixeqtl_rand/trans_eqtl.txt")
        jpa_outfile    = os.path.join(opts.srcdir, simdir, "tejaas_rand/jpa/all_jpa_pvals.txt")
        rr_outfile     = os.path.join(opts.srcdir, simdir, "tejaas_rand/permnull_sb{:s}/rr.txt".format(perm_sbeta))
    
    true_trans = list()
    with open(sim_trans_file, 'r') as fin:
        for line in (fin):
            linesplit = line.strip().split()
            true_trans.append(linesplit[1])

    if opts.method == 'jpa':        simdata[i] = read_data('jpa', jpa_outfile, true_trans)
    if opts.method == 'matrixeqtl': simdata[i] = read_data('meqtl', meqtl_outfile, true_trans)
    if opts.method == 'rr':         simdata[i] = read_data('rr', rr_outfile, true_trans)

    _nsel, _tpr, _ppv, _valids, _thres, _fpr = roc.confusion_matrix(res[i])
    rocdata['nsel'][i] = _nsel
    rocdata['tpr'][i] = _tpr
    rocdata['ppv'][i] = _ppv
    rocdata['fpr'][i] = _valids
    rocdata['valids'][i] = _thres
    rocdata['thres'][i] = _fpr

nsim = len(simlist)
ninterp = int(2 * np.max([len(x) for x in rocdata['thres']]))

for desc in opts.whichplots:
    xstr, ystr = tuple(desc.split('_'))
    xx = rocdata[xstr]
    yy = rocdata[ystr]
    xmax = int(np.max([max(x) for x in xx]) + 0.5)
    xvals = np.linspace(0, xmax, ninterp)
    reverse = False
    if yy[0][-1] < yy[0][0]: reverse = True
    if not reverse:
        yest_arr = np.array([np.interp(xvals, xx[i], yy[i]) for i in range(nsim)])
    else:
        yest_arr = np.array([np.interp(xvals, xx[i][::-1], yy[i][::-1]) for i in range(nsim)])
    yest = np.mean(yest_arr, axis = 0)
    yerr = np.std (yest_arr, axis = 0)
    filename = "{:s}_{:s}_{:s}.txt".format(opts.method, xstr, ystr)
    if opts.shuffle: filename = "shuffled_{:s}_{:s}_{:s}.txt".format(opts.method, xstr, ystr)
    outfile = os.path.join(opts.outdir, filename)
    with open(outfile, 'w') as fout:
        fout.write("{:s} {:s} std_col2\n".format(xstr, ystr))
        for x, y, e in zip(xvals, yvals, yerr):
            fout.write("{:g} {:g} {:g}\n".format(x, y, e))
