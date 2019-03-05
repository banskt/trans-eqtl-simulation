# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Script for benchmarking by F. Simonetti
# Minor modifications by Saikat Banerjee to include file in bash pipeline, remove input of chromosome

library("optparse")
library(MatrixEQTL)

read_vcf <- function (vcf_file_name, maf_filter=T, maf_thres=0.1) {
    message("Reading Genotype ...")
    message(vcf_file_name)

    # read two times the vcf file, first for the columns names, second for the data
    # Read the vcf file to a string
    tmp_vcf_strings = readLines(vcf_file_name)
    header_colnum = grep("#CHROM", tmp_vcf_strings)
    total_colnum  = length(tmp_vcf_strings)
    data_colnum   = total_colnum - header_colnum

    header_string = tmp_vcf_strings[-(header_colnum + 1):-total_colnum]
    colnames<-unlist(strsplit(header_string,"\t"))

    snps_data = read.table(text="", col.names = colnames)
    colnames(snps_data) = colnames
    data_strings  = tmp_vcf_strings[ (header_colnum + 1): total_colnum]
    for (i in 1 : data_colnum) {
        tmpstrvals = unlist(strsplit(data_strings[i], "\t"))
        for (j in 1:9) {
            snps_data[i, j] = tmpstrvals[j]
        }
        for (j in 10:length(colnames)) {
            snps_data[i, j] = as.numeric(unlist(strsplit(tmpstrvals[j], ":"))[3])
        }
    }

    # get SNPs positions for cis and trans analysis (before cropping the snp matrix)
    row_names = snps_data[,3]
    snpspos = snps_data[,c(3,1,2)]
    snpspos[,2] = paste("chr", snpspos[,2], sep="")
    colnames(snpspos) = c("snpid","chr","pos")
    snpspos[,3] = as.numeric(snpspos[,3])

    snps_data = snps_data[, 10:ncol(snps_data)]
    rownames(snps_data) = row_names

    snps_mat = as.matrix(snps_data)

    return (list(snps_mat, snpspos))
}

option_list = list(
    make_option(c("-g", "--genotype"), type="character", default=NULL,
              help="path to genotype file name", metavar="character"),
    make_option(c("-u", "--selectdonors"), type="character", default=NULL,
              help="path to selected donor ids file", metavar="character"),
    make_option(c("-i", "--geneinfo"), type="character", default=NULL,
              help="path to gene info file name", metavar="character"),
    make_option(c("-s", "--datatype"), type="character", default=NULL,
              help="type of dataset [gtex, cardiogenics]", metavar="character"),
    make_option(c("-e", "--expression"), type="character", default=NULL,
              help="path to gene expression file name", metavar="character"),
    make_option(c("-p", "--pvalcis"), type="numeric", default=1e-3,
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-t", "--pvaltrans"), type="numeric", default=1e-3,
              help="pvalue threshold for cis eQTLs [default= %default]", metavar="number"),
    make_option(c("-o", "--outfilecis"), type="character", default=NULL,
              help="path to output file for cis-eQTLs", metavar="character"),
    make_option(c("-q", "--outfiletrans"), type="character", default=NULL,
              help="path to output file for trans-eQTLs", metavar="character"),
    make_option(c("-m", "--model"), type="character", default="modelLINEAR",
              help="Model to use from modelANOVA, modelLINEAR, or modelLINEAR_CROSS [default \"%default\"]", metavar="character"),
    make_option(c("-r", "--randomize"), action="store_true", default=FALSE,
              help="Randomize the gene expression"),
    make_option(c("-R", "--shufflewith"), type="character", default=NULL,
              help="file with shuffled donor ids for genotype shuffling", metavar="character"),
    make_option(c("-w", "--window"), type="numeric", default=1e6,
              help="length of cis window in base pairs [default= %default]", metavar="number")
);

opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

if(is.null(opt$geneinfo)) {
    stop("Gene info file is missing");
} else {
    genepos = read.table(opt$geneinfo, header = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"));
}

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Output file name
output_file_name_cis = opt$outfilecis;
output_file_name_tra = opt$outfiletrans;

SNP_file_name = opt$genotype;
res = read_vcf(SNP_file_name);
snps_mat = res[[1]] #genotype matrix
snpspos  = res[[2]] #SNP position info

if (!is.null(opt$shufflewith)) {
    message("Shuffling genotype using supplied donor IDs");
    shuffled_colnames = scan(opt$shufflewith, what="", sep="\n")
    colnames(snps_mat) = shuffled_colnames;
}

snps = SlicedData$new();
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$CreateFromMatrix( snps_mat );

## Load gene expression data
expression_file_name = opt$expression;
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

if (opt$randomize) {
    message("Randomizing gene expression values by permuting samples")
    expr_mat = as.matrix(gene);
    shuffled = t(apply(expr_mat, 1, sample));
    colnames(shuffled) = colnames(gene);
    gene$CreateFromMatrix(shuffled);
}

# match columns of samples from genotype and expression
if (!is.null(opt$selectdonors)) {
    message("Selecting donors from user list")
    selectdonors = read.csv(file=opt$selectdonors, sep=" ", stringsAsFactors=F, header=F, comment.char="#")
    col_index1 = match(selectdonors[[1]], colnames(gene))
    col_index2 = match(selectdonors[[1]], colnames(snps))
} else {
    common_donors = intersect(colnames(gene), colnames(snps))
    col_index1 = match(common_donors, colnames(gene))
    col_index2 = match(common_donors, colnames(snps))
}
if (all(colnames(snps)[col_index2] == colnames(gene)[col_index1])) {
    gene$ColumnSubsample(col_index1)
    snps$ColumnSubsample(col_index2)
}

# Covariates file name. Set to character() for no covariates
if (is.null(opt$covariates)) {
    cvrt = SlicedData$new();
} else {

    ## Load covariates
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(opt$covariates);
    }
}

if (opt$datatype == "cardiogenics") {
    #trim genepos ensembl ids
    ensembls = sapply(genepos$geneid, function(x){strsplit(x, '.', fixed=T)})
    ensembls = as.character(sapply(ensembls, "[[", 1))
    genepos$geneid = ensembls
}


# Only associations significant at this level will be saved
pvOutputThreshold_cis = opt$pvalcis;
pvOutputThreshold_tra = opt$pvaltrans;

# Distance for local gene-SNP pairs
cisDist = opt$window;
message( cisDist );

# Error covariance matrix. Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

## Run the analysis
me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name_tra,
        pvOutputThreshold = pvOutputThreshold_tra,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        cisDist = cisDist,
        snpspos = snpspos,
        genepos = genepos,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);

## Results:
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat(paste('Detected eQTLs:', me$all$neqtls, '\n', sep=""));
