#' get allele count matrix for all single cells in bam file
#' @description this is a wrapper function for the scAlleleCount.py python script
#' @param snps data.frame of four colums: chromosome, position, reference allele, alternate allele
#' @param bamFile name of bamfile, must be sorted and indexed
#' @param cellBarcodes character vector of barcodes to scan for
#' @param verbose logical verbosity, default: TRUE
#' @param scAlleleCountExec name of scAlleleCount.py python script
#' @return list of three sparse matrices: refmat, altmat, and covmat, note thatn covmat will not be the sum of refmat and altmat because of genomic locations where the snp has a base othern than ref or alt
getFastCellAlleleCount <- function(snps,bamFile,cellBarcodes, verbose=TRUE, scAlleleCountExec = './scAlleleCount.py'){
    ## check the input
    if(class(snps) != "data.frame") stop("snps is not a data.frame")
    if(ncol(snps) != 4) stop("snps must have exactly 4 columns");
    if(class(snps[[2]]) != 'integer') stop('snps column 2 must be integer chromosomal positions')
    if(!all(sort(unique(snps[[3]])) == c("A","C","G","T"))) stop('snps column 3 must by A,T,G, or C')
    if(!all(sort(unique(snps[[4]])) == c("A","C","G","T"))) stop('snps column 4 must by A,T,G, or C')

    ## convert integers to chars if they have been read as numeric
    snps[[1]] <- as.character(snps[[1]])

    ## filter out snps that are actually indels and warn
    kv <- (sapply(snps[[3]],nchar) == 1 & sapply(snps[[4]], nchar) == 1)
    if(!all(kv)) warning('Some variants were indels, removing');
    snps <- snps[kv,]
    
    ## check barcodes
    if(class("barcodes") != "character") stop('barcodes must be a character vector')

    ## Get a working tmp directory
    w <- tempdir()
    
    ## todo improve names so that there are no collisions if thie
    ## is run in parallel in a single session
    
    ## make barcodes files
    tmp.barcodes.file <- file.path(w,'barcodes.txt')
    tmp.snps.file <- file.path(w, 'snps.txt')
    tmp.output.prefix <- file.path(w, 'out')
    tmp.output.prefix

    ## write the input for the python command
    write.table(cellBarcodes,tmp.barcodes.file, quote=F, col.names=F, row.names=F)
    write.table(snps, tmp.snps.file, quote=F, col.names=F, row.names=F)

    ## Build command
    verbose.arg = ifelse(verbose,' -v ','')

    cmd <- paste0(scAlleleCountExec, verbose.arg ,' --snps ', tmp.snps.file, ' --barcodes ', tmp.barcodes.file,
                  ' --output-format mm --max-depth 9999999 --output-prefix ', tmp.output.prefix, ' ',
                  bamFile)

    ## Call python program
    cmd.ret <- system(cmd)
    if (cmd.ret != 0) {
        stop(paste0('An error occured while running ',scAlleleCountExec))
    }

    ## read the output matrices
    refmatname <- paste0(tmp.output.prefix, 'refmat.mtx')
    altmatname <- paste0(tmp.output.prefix, 'altmat.mtx')
    covmatname <- paste0(tmp.output.prefix, 'covmat.mtx')

    ## read the output matrices
    refmat <- Matrix::readMM(refmatname)
    altmat <- Matrix::readMM(altmatname)
    covmat <- Matrix::readMM(covmatname)

    ## Fix column and row names
    snpnames <- paste0(snps$V1, ':', snps$V2, '-', snps$V2)
    rownames(refmat) <- rownames(altmat) <- rownames(covmat) <- cellBarcodes
    colnames(refmat) <- colnames(altmat) <- colnames(covmat) <- snpnames

    ## transpose for compatibility
    refmat <- Matrix::t(refmat);
    altmat <- Matrix::t(altmat);
    covmat <- Matrix::t(covmat);

    ## Cleanup
    file.remove(tmp.barcodes.file)
    file.remove(tmp.snps.file)
    file.remove(refmatname,altmatname,covmatname)
    
    ## return
    list(refmat=refmat, altmat=altmat, covmat=covmat)
}
