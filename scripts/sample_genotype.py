import argparse
import gzip
import random
import os

def parse_args():

    parser = argparse.ArgumentParser(description="Select random SNPs from genotype gunzipped VCF files")

    parser.add_argument('--input',
                        nargs='*',
                        type=str,
                        dest='vcfgz_files',
                        help='list of input files')

    parser.add_argument('--out',
                        type=str,
                        dest='outfile',
                        help='output gunzipped VCF file')

    parser.add_argument('--ntarget', '-n',
                        type=int,
                        dest='ntarget',
                        help='number of SNPs to be selected')

    opts = parser.parse_args()
    return opts


def meta_info(filename):
    count = 0
    rsidlist = list()
    with gzip.open(filename, 'r') as vcfstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[0] != '#':
                count += 1
                rsidlist.append(linestrip.split("\t")[2])
    return count, rsidlist


def get_headers(filename):
    lines = list()
    with gzip.open(filename, 'r') as vcfstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[0] == '#':
                lines.append(linestrip + "\n")
            else:
                break
    return lines


def get_snps(filename, n, rsidlist):
    choose = random.sample(rsidlist, n)
    lines = list()
    mafs = list()
    rsids = list()
    chrms = list()
    bpposs = list()
    alleles = list()
    with gzip.open(filename, 'r') as vcfstream:
        for line in vcfstream:
            linestrip = line.decode().strip()
            if linestrip[0] != '#':
                linesplit = linestrip.split("\t")
                rsid = linesplit[2]
                if rsid in choose:
                    rsids.append(rsid)
                    lines.append(linestrip + "\n")
                    dsidx = linesplit[8].split(':').index("DS")
                    ds = [x.split(":")[dsidx] for x in linesplit[9:]]
                    ds_notna = [float(x) for x in ds if x != "."]
                    maf = sum(ds_notna) / 2.0 / len(ds_notna)
                    mafs.append(maf)
                    chrms.append(linesplit[0])
                    bpposs.append(linesplit[1])
                    alleles.append(f"{linesplit[2]} {linesplit[3]}")
    return lines, mafs, rsids, chrms, bpposs, alleles


if __name__ == '__main__':

    opts = parse_args()
    
    nsnps = list()
    allrsid = list()
    rsidlist = list()

    maflist = list()
    chrlist = list()
    bplist = list()
    allelelist = list()

    for vcfgz_file in opts.vcfgz_files:
        n, rsids = meta_info(vcfgz_file)
        nsnps.append(n)
        allrsid.append(rsids)
        print (f'{vcfgz_file}: {n} SNPs')

    nselect = [int(x * opts.ntarget / sum(nsnps) + 0.5) for x in nsnps] # split in proper ratio
    nselect[-1] = opts.ntarget - sum(nselect[:-1]) # adjust, algebraically it is impossible to have opts.ntarget < sum(nselect[:-1])

    headers = get_headers(opts.vcfgz_files[0])
    print (f'Read {len(headers)} header lines')
    print(f'Selecting {nselect} SNPs from the VCF files')

    with gzip.open(opts.outfile, 'w') as outfile:
        for line in headers:
            outfile.write(line.encode('utf-8'))
        for i, vcfgz_file in enumerate(opts.vcfgz_files):
            snplines, mafs, rsids, chrms, bpposs, alleles = get_snps(vcfgz_file, nselect[i], allrsid[i])
            for i, line in enumerate(snplines):
                outfile.write(line.encode('utf-8'))
                maflist.append(mafs[i])
                rsidlist.append(rsids[i])
                chrlist.append(chrms[i])
                bplist.append(bpposs[i])
                allelelist.append(alleles[i])

    maffile = os.path.splitext(os.path.splitext(opts.outfile)[0])[0] + '.altfreq' # split .gz and then .vcf
    with open(maffile, 'w') as outfile:
        for i, rsid in enumerate(rsidlist):
            outfile.write(f"{chrlist[i]} {bplist[i]} {rsid} {allelelist[i]} {maflist[i]:5.3f}\n")

    shuffile = os.path.splitext(maffile)[0] + '.shuf'
    donors = headers[-1].strip().split('\t')[9:]
    random.shuffle(donors)
    with open(shuffile, 'w') as outfile:
        for donor in donors:
            outfile.write(f"{donor}\n")
