export <- function (x, dest_dir=getwd(), dest_dir_create=TRUE, delim=c("tsv","csv"), input.data=TRUE)
{
    
    if ((!dir.exists(dest_dir))&&(dest_dir_create)) {
        dir.create(dest_dir)
    }
    
    delim = match.arg(delim)
    if (delim=="tsv") {
        sep="\t"
        ext="tsv"
    } else if (delim=="csv") {
        sep=","
        ext="csv"
    }
    
    if (x$method=="deconstructSigs") {
        args.desc = c("method","sig.bkg.adj","signature.cutoff","noise",
        "neg.binom.size","n_rdm")
        args.value = list(x$method,x$sig.bkg.adj,x$signature.cutoff,x$noise,
        x$neg.binom.size,x$n_rdm)
    } else if (x$method=="MutationalPatterns") {
        args.desc = c("method","sig.bkg.adj","signature.cutoff","noise","neg.binom.size","n_rdm")
        args.value = list(x$method,x$sig.bkg.adj,x$signature.cutoff,x$noise,
        x$neg.binom.size,x$n_rdm)
    }
    args.value[sapply(args.value, is.null)] <- NA
    output = rbind(c("argument","value"),cbind(args.desc,unlist(args.value)))
    write(t(output),ncolumns=ncol(output),
    file=file.path(dest_dir,paste0("arguments.",ext)),sep=sep)
    
    sig_label = names(x$weights[[1]])
    n_sig = length(sig_label)
    sample_label = names(x$weights)
    n_sample = length(sample_label)
    
    if (input.data) {
        output = rbind(c("channel",sig_label),cbind(rownames(x$sig),x$sig))
        outfile = file.path(dest_dir,paste0("sig.",ext))
        write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)
        
        output = rbind(c("channel",sample_label),cbind(rownames(x$spectra),x$spectra))
        outfile = file.path(dest_dir,paste0("spectra.",ext))
        write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)
    }
    
    tmp = t(matrix(unlist(x$weights),byrow=TRUE,ncol=length(x$weights[[1]])))
    output = rbind(c("signatures",sample_label),cbind(sig_label,tmp))
    outfile = file.path(dest_dir,paste0("weights.",ext))
    write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)

    for (i_sample in 1:n_sample) {
        output = rbind(c("rdm_sample",sig_label),cbind(rownames(x$weights_rdm[[i_sample]]),x$weights_rdm[[i_sample]]))
        outfile = file.path(dest_dir,paste0("weights_rdm_",sample_label[i_sample],".",ext))
        write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)
    }
    
    output = rbind(c("signatures",sample_label),cbind(sig_label,x$map_pval))
    outfile = file.path(dest_dir,paste0("map_pval.",ext))
    write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)
    
    if (x$method=="MutationalPatterns") {
        tmp = t(matrix(unlist(x$counts),byrow=TRUE,ncol=length(x$counts[[1]])))
        output = rbind(c("signatures",sample_label),cbind(sig_label,tmp))
        outfile = file.path(dest_dir,paste0("counts.",ext))
        write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)

        for (i_sample in 1:n_sample) {
            output = rbind(c("rdm_sample",sig_label),cbind(rownames(x$counts_rdm[[i_sample]]),x$counts_rdm[[i_sample]]))
            outfile = file.path(dest_dir,paste0("counts_rdm_",sample_label[i_sample],".",ext))
            write(t(output),ncolumns=ncol(output),file=outfile,sep=sep)
        }
    }
}
