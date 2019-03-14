import os
import gzip
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from readvcf import ReadVCF

from scipy.interpolate import interp1d
import mpl_stylesheet
mpl_stylesheet.banskt_presentation(fontfamily = 'system', colors = 'kelly')

def parse_args():

    parser = argparse.ArgumentParser(description = 'Simulate gene expression using input parameters')

    parser.add_argument('-n',
                        type = int,
                        required = True,
                        dest = 'ngene',
                        metavar = 'INT',
                        help = 'number of genes')


    parser.add_argument('-s',
                        type = int,
                        required = True,
                        dest = 'nsample',
                        metavar = 'INT',
                        help = 'number of samples')


    parser.add_argument('-f',
                        type = int,
                        required = True,
                        dest = 'ncf',
                        metavar = 'INT',
                        help = 'number of confounding factors (CF)')


    parser.add_argument('-c',
                        type = int,
                        required = True,
                        dest = 'ncis',
                        metavar = 'INT',
                        help = 'number of cis genes')


    parser.add_argument('-t',
                        type = int,
                        required = True,
                        dest = 'ntf',
                        metavar = 'INT',
                        help = 'number of transciption factors (TF)')


    parser.add_argument('-g',
                        type = int,
                        required = True,
                        dest = 'ntrans',
                        metavar = 'INT',
                        help = 'number of target genes for each TF')


    parser.add_argument('--tfcis',
                        nargs = 2,
                        type = float,
                        default = [4.0, 0.1],
                        dest = 'tfcis_params',
                        metavar = 'FLOAT',
                        help = 'Shape (k) and scale (theta) parameters for the Gamma distribution from which cis-effects of TF are sampled')


    parser.add_argument('--tftrans',
                        nargs = 2,
                        type = float,
                        default = [20, 0.02],
                        dest = 'tftrans_params',
                        metavar = 'FLOAT',
                        help = 'Shape (k) and scale (theta) parameters for the Gamma distribution from which trans-effects are sampled')


    parser.add_argument('--cf',
                        nargs = 3,
                        type = float,
                        default = [0.5, 0, 0.1],
                        dest = 'cf_params',
                        metavar = 'FLOAT',
                        help = 'Sparsity (pi), mean (mu) and standard deviation (sigma) for the spike-and-slab distribution from which effects of CF are sampled ')


    parser.add_argument('--gxcorr',
                        type = str,
                        dest = 'gxcorr_file',
                        metavar = 'FILE',
                        help = 'Numpy file (.npy) for Q @ sqrt(W) matrix for generating background expression with covariance S, where S = Q @ W @ Q*')


    parser.add_argument('--noise',
                        nargs = 2,
                        type = float,
                        default = [100, 0.01],
                        dest = 'noise_params',
                        metavar = 'FLOAT',
                        help = 'Shape (k) and scale (theta) parameters for the Gamma distribution from which standard deviation of noise is sampled')


    parser.add_argument('--cis',
                        nargs = 1,
                        type = float,
                        default = [0.6],
                        dest = 'cis_params',
                        metavar = 'FLOAT',
                        help = 'Fixed value of the cis-effects')


    parser.add_argument('--geno',
                        type = str,
                        required = True,
                        dest = 'gtfile',
                        metavar = 'FILE',
                        help = 'gunzipped VCF file of the genotype')


    parser.add_argument('--out',
                        type = str,
                        required = True,
                        dest = 'outfile',
                        metavar = 'FILE',
                        help = 'filename of the simulated gene expression')


    opts = parser.parse_args()
    return opts


def spike_and_slab(pi, mu, sigma, size=None):
    n = size if size is not None else 1
    res = np.zeros(n)
    for i in range(n):
        mrand = np.random.uniform(0, 1)
        if mrand <= pi:
            res[i] = np.random.normal(mu, sigma)
        else:
            res[i] = 0
    if size is None:
        return res[0]
    else:
        return res
    
def sample_sign(n):
    return np.random.choice([-1, 1], size=n, p=[0.5, 0.5])


# Simulate noise for n genes and s samples.
def simulate_noise(n, s, params):
    X = np.zeros((n, s))
    sigma = np.random.gamma(params[0], params[1], size = n)
    for i in range(n):
        X[i, :] = np.random.normal(0, sigma[i], size = s)
    return X

# Simulate f counfounding factors for n genes and s samples
def simulate_confounders(n, s, f, params):
    CF = np.zeros((f, s))
    beta = np.zeros((n, f))
    for i in range(f):
        CF[i, :] = np.random.normal(0, 1, size = s)
        beta[:, i] = spike_and_slab(params[0], params[1], params[2], size = n)
    X = np.einsum('ij, jk', beta, CF)
    return X, CF, beta


def sample_correlation(corrfile, n, s):
    Xrand = np.random.normal(0, 1, size = s * n).reshape(n, s)
    QsqrtW = np.load(corrfile)
    X = np.dot(QsqrtW, Xrand)
    return Xrand, X 


def simulate_cis(GT, icis, itf, cisparams, tfparams):
    n = GT.shape[0]
    s = GT.shape[1]
    X = np.zeros((n, s))
    itf_idx = np.array([np.where(icis == i)[0][0] for i in itf])
    bcis = np.repeat(cisparams, icis.shape[0])
    btf  = np.random.gamma(tfparams[0], tfparams[1], size = itf.shape[0])
    bvec = np.multiply(bcis, sample_sign(bcis.shape[0]))
    bvec[itf_idx] = np.multiply(btf, sample_sign(btf.shape[0]))
    for i, idx in enumerate(icis):
        thisgt = GT[idx, :]
        X[idx, :] = bvec[i] * thisgt
    return X, bvec


def simulate_trans(GX, itf, itrans, params):
    btrans = np.zeros((GX.shape[0], itf.shape[0]))
    for i, targetidx in enumerate(itrans):
        btrans[targetidx, i] = np.random.gamma(params[0], params[1], size = targetidx.shape[0])
    X = np.einsum('ij, jk', btrans, GX[itf, :])
    return X, btrans


def normalize_expr(Y):
    # requires a G x N matrix, where G is the number of genes, and Y is the number of samples
    newY = (Y - np.mean(Y, axis = 1).reshape(-1, 1)) / np.std(Y, axis = 1).reshape(-1, 1)
    return newY


def write_expression(GX, donors, filename):
    with open(filename, 'w') as fout:
        header = 'gene_ids\t' + '\t'.join(donors) + '\n'
        fout.write(header)
        for i in range(GX.shape[0]):
            line = 'ENSG{:06d}\t'.format(i) + '\t'.join(['{:g}'.format(x) for x in GX[i, :]]) + '\n'
            fout.write(line)


def write_background(GX, donors, filename):
    bg_filename = os.path.splitext(filename)[0] + '.bground'
    with open(bg_filename, 'w') as fout:
        header = 'gene_ids\t' + '\t'.join(donors) + '\n'
        fout.write(header)
        for i in range(GX.shape[0]):
            line = 'ENSG{:06d}\t'.format(i) + '\t'.join(['{:g}'.format(x) for x in GX[i, :]]) + '\n'
            fout.write(line)


def write_gtf(snpinfo, filename):
    gtf_filename = os.path.splitext(filename)[0] + '.gtf.gz'
    genepos_filename = os.path.splitext(filename)[0] + '.genepos'
    with gzip.open(gtf_filename, 'w') as fout, open(genepos_filename, 'w') as fposout:
        header = "##description: placeholder gtf file generated by simulation\n"
        fout.write(header.encode('utf-8'))
        header = "geneid\tchr\tleft\tright\n"
        fposout.write(header)
        for i in range(len(snpinfo)):
            chrm = snpinfo[i].chrom
            start = snpinfo[i].bp_pos
            end = start
            attr = "gene_id \"ENSG{0:06d}.1\"; transcript_id \"ENSG{0:06d}.1\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"GEN{0:06d}\";".format(i)
            line = "chr{:d}\tSIMULATION\tgene\t{:d}\t{:d}\t.\t+\t.\t{:s}\n".format(chrm, start, end, attr)
            fout.write(line.encode('utf-8'))
            line = "ENSG{:06d}\tchr{:d}\t{:d}\t{:d}\n".format(i, chrm, start, end)
            fposout.write(line)


def write_eQTLs(snpinfo, cis_genes, tf_genes, trans_genes, filename):
    ceqtl_filename = os.path.splitext(filename)[0] + '.cis'
    teqtl_filename = os.path.splitext(filename)[0] + '.trans'
    with open(ceqtl_filename, 'w') as fout:
        for i, cis in enumerate(cis_genes):
            if cis not in tf_genes:
                snp = snpinfo[cis]
                line = "chr{:d}\t{:s}\t{:d}\tENSG{:06d}\n".format(snp.chrom, snp.varid, snp.bp_pos, cis)
                fout.write(line)
    with open(teqtl_filename, 'w') as fout:
        for i, tf in enumerate(tf_genes):
            snp = snpinfo[tf]
            gnet = ','.join(['{:d}'.format(x) for x in trans_genes[i]])
            line = "chr{:d}\t{:s}\t{:d}\tENSG{:06d}\t{:s}\n".format(snp.chrom, snp.varid, snp.bp_pos, tf, gnet)
            fout.write(line)


def write_confounders(CF, betacf, donors, filename):
    cf_filename = os.path.splitext(filename)[0] + '.cf'
    with open(cf_filename, 'w') as fout:
        header = 'cf_ids\t' + '\t'.join(donors) + '\n'
        fout.write(header)
        for i in range(CF.shape[0]):
            line = 'CF{:04d}\t'.format(i) + '\t'.join(['{:g}'.format(x) for x in CF[i, :]]) + '\n'
            fout.write(line)
    betacf_filename = os.path.splitext(filename)[0] + '.betacf'
    with open(betacf_filename, 'w') as fout:
        header = 'gene_ids\t' + '\t'.join(['CF{:04d}'.format(i) for i in range(betacf.shape[1])]) + '\n'
        fout.write(header)
        for i in range(betacf.shape[0]):
            line = 'ENSG{:06d}\t'.format(i) + '\t'.join(['{:g}'.format(x) for x in betacf[i, :]]) + '\n'
            fout.write(line)
    #np.savetxt(betacf_filename, betacf)


def write_ciseffects(CT, cis_genes, filename):
    cis_filename = os.path.splitext(filename)[0] + '.betacis'
    with open(cis_filename, 'w') as fout:
        header = 'gene_ids\tcis_effect\n'
        fout.write(header)
        for i in range(CT.shape[0]):
            line = f'ENSG{cis_genes[i]:06d}\t{CT[i]}\n'
            fout.write(line)


def write_transeffects(BT, itf, filename):
    trans_filename = os.path.splitext(filename)[0] + '.betatrans'
    with open(trans_filename, 'w') as fout:
        header = 'gene_ids\t' + '\t'.join(['{:d}'.format(x) for x in itf]) + '\n'
        fout.write(header)
        for i in range(BT.shape[0]):
            line = 'ENSG{:06d}\t'.format(i) + '\t'.join(['{:g}'.format(x) for x in BT[i, :]]) + '\n'
            fout.write(line)


def plot_components(ax, data, xmax, xmin, nbin, ninterp, legend, interp = True, ymax = 0.0):
    bins = np.linspace(xmin, xmax, nbin)
    xbin = [(bins[i] + bins[i+1]) / 2 for i in range(bins.shape[0] - 1)] # centers of the bins
    x = np.linspace(xbin[0], xbin[-1], ninterp)
    nrow = data.shape[0]
    ncol = data.shape[1]
    ymaxcal = 0
    for i in range(nrow):
        yhist, _ = np.histogram(data[i, :], bins = bins)
        yfreq = yhist / ncol
        if interp:
            fy = interp1d(xbin, yfreq, kind='cubic')
            y = fy(x)
            ax.plot(x, y, alpha = 0.3)
            ymaxnew = np.max(y)
        else:
            ax.plot(xbin, yfreq, alpha = 0.3)
            ymaxnew = np.max(yfreq)
        if ymaxnew > ymaxcal:
            ymaxcal = ymaxnew
    if ymax == 0:
        ymaxset = round(ymaxcal + 0.05, 1)
    else:
        ymaxset = ymax
    ax.set_ylim([0, ymaxset])
    ax.text(0.1, 0.85, legend, horizontalalignment='left', transform=ax.transAxes)
    return ymaxset

    
if __name__ == '__main__':

    opts = parse_args()
    
    # Read genotype
    readvcf = ReadVCF(opts.gtfile, 0, opts.ngene)
    GT = readvcf.dosage

    # Select cis-genes from all available genes
    cis_genes = np.sort(np.random.choice(opts.ngene, opts.ncis, replace=False))

    # Select transciption factors (TF) from the cis-gens
    tf_idx = np.sort(np.random.choice(opts.ncis, opts.ntf, replace=False))
    tf_genes = np.array([cis_genes[i] for i in tf_idx])

    # Select trans-genes from all non-TF genes
    nontf_genes = np.array([i for i in range(opts.ngene) if i not in tf_genes])
    trans_genes = [np.sort(np.random.choice(nontf_genes, opts.ntrans, replace=False)) for i in tf_genes]

    # Simulate noise and CF
    if opts.gxcorr_file is None:
        print ("Generating noise and confounder from prior distributions using supplied options.")
        GX_noise = simulate_noise(opts.ngene, opts.nsample, opts.noise_params)
        GX_cf, confounders, cf_effectsize = simulate_confounders (opts.ngene, opts.nsample, opts.ncf, opts.cf_params)
    else:
        print ("Generating background gene expression by sampling GTEx data")
        GX_noise, GX_cf = sample_correlation(opts.gxcorr_file, opts.ngene, opts.nsample)

    # Simulate cis-effects
    GX_cis, cis_effectsize = simulate_cis (GT, cis_genes, tf_genes, opts.cis_params, opts.tfcis_params)

    # Sum expression components
    if opts.gxcorr_file is None:
        print ("Summing noise, confounders and cis-effects")
        GX_tmp = GX_noise + GX_cf + GX_cis
    else:
        print ("Summing background and cis-effects")
        GX_tmp = GX_cf + GX_cis

    # Simulate trans-effects
    GX_trans, trans_effectsize = simulate_trans (GX_tmp, tf_genes, trans_genes, opts.tftrans_params)

    # Total normalized gene expression
    GX_raw = GX_tmp + GX_trans
    GX = normalize_expr(GX_raw)

    # Output
    write_expression(GX, readvcf.donor_ids, opts.outfile)
    write_gtf(readvcf.snpinfo, opts.outfile)
    write_eQTLs(readvcf.snpinfo, cis_genes, tf_genes, trans_genes, opts.outfile)
    if opts.gxcorr_file is None:
        write_confounders(confounders, cf_effectsize, readvcf.donor_ids, opts.outfile)
    else:
        write_background(GX_cf, readvcf.donor_ids, opts.outfile)
    write_ciseffects(cis_effectsize, cis_genes, opts.outfile)
    write_transeffects(trans_effectsize, tf_genes, opts.outfile)

    # Plot the components
    fig = plt.figure(figsize = (16, 18))
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    ax6 = fig.add_subplot(326)
    
    nbin = 30
    ninterp = 100
    ngene_plot = 100
    xmax = int(max(abs(np.min(GX[:ngene_plot, :])), abs(np.max(GX[:ngene_plot, :]))) + 0.5)
    xmin = -xmax
    
    ymax = plot_components(ax4, GX_cf   [:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Confounding')
    _    = plot_components(ax1, GX_raw  [:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Total - with trans', ymax = ymax)
    _    = plot_components(ax2, GX_tmp  [:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Total - without trans', ymax = ymax)
    _    = plot_components(ax3, GX_noise[:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Noise', ymax = ymax)
    _    = plot_components(ax5, GX_cis  [:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Cis', ymax = 1.0)
    _    = plot_components(ax6, GX_trans[:ngene_plot, :], xmax, xmin, nbin, ninterp, 'Trans', ymax = 1.0)
    
    plt.tight_layout()
    plotfile = os.path.splitext(opts.outfile)[0] + '_components.pdf'
    plt.savefig(plotfile, bbox_inches='tight')
