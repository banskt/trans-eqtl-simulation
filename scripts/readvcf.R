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
