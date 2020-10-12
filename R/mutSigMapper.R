mutSigMapper = function(spectra, sig=NULL, ref=c("cosmic_v2","cosmic_v3","cosmic_v3_exome","cosmic_v3.1","mutagen53"), method=c("MutationalPatterns","deconstructSigs"),
    sig.bkg.adj=c("none","1/exome","1/genome","exome/genome","genome/exome","custom"), sig.bkg.adj.custom=NULL,
    signature.cutoff=0.06,
    noise=c("poisson","neg.binom"), neg.binom.size=NULL, n_rdm=200,
    save_obj=FALSE, dest_dir=getwd(), dest_dir_create=TRUE, dest_dir_create_recur=FALSE,
    dest_obj="mutSigMap.Robj")
{
    if (save_obj) {
        if ((!dir.exists(dest_dir))&&(dest_dir_create)) {
            dir.create(dest_dir, showWarnings=FALSE, recursive=dest_dir_create_recur)
        }
    }
    if (!((min(spectra)>=0)&(max(spectra)>1))) {
        stop("Spectra must contain mutation counts")
    }
    if (is.null(sig)) {
        ref = match.arg(ref)
        sig = read.table(system.file("extdata",paste0(match.arg(ref),".txt"),package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
    }
    method = match.arg(method)
    noise = match.arg(noise)
    if (noise=="neg.binom") {
        neg.binom.size = as.numeric(neg.binom.size)
        if (!(neg.binom.size>0)) {
            stop("neg.binom.size must be a positive (integer or real) number")
        }
    }
    
    spectra = spectra[match(rownames(sig),rownames(spectra)),,drop=FALSE]
    if (!identical(rownames(sig),rownames(spectra))) {
        stop("Sig vs spectra substitution labeling mismatch")
    }
    
    sig.bkg.adj = match.arg(sig.bkg.adj)
    if (sig.bkg.adj=="1/exome") {
        ref_bkg_exome = read.table(system.file("extdata","bkg_exome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_exome = ref_bkg_exome[match(rownames(sig),rownames(ref_bkg_exome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_exome))) {
            stop("Sig vs ref_bkg_exome substitution labeling mismatch")
        }
        for (i_sig in 1:ncol(sig)) {
            sig[,i_sig] = sig[,i_sig]/ref_bkg_exome[,"Counts"]
            sig[,i_sig] = sig[,i_sig]/sum(sig[,i_sig])
        }
    } else if (sig.bkg.adj=="1/genome") {
        ref_bkg_genome = read.table(system.file("extdata","bkg_genome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_genome = ref_bkg_genome[match(rownames(sig),rownames(ref_bkg_genome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_genome))) {
            stop("Sig vs ref_bkg_genome substitution labeling mismatch")
        }
        for (i_sig in 1:ncol(sig)) {
            sig[,i_sig] = sig[,i_sig]/ref_bkg_genome[,"Counts"]
            sig[,i_sig] = sig[,i_sig]/sum(sig[,i_sig])
        }
    } else if (sig.bkg.adj=="exome/genome") {
        ref_bkg_exome = read.table(system.file("extdata","bkg_exome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_exome = ref_bkg_exome[match(rownames(sig),rownames(ref_bkg_exome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_exome))) {
            stop("Sig vs ref_bkg_exome substitution labeling mismatch")
        }
        ref_bkg_genome = read.table(system.file("extdata","bkg_genome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_genome = ref_bkg_genome[match(rownames(sig),rownames(ref_bkg_genome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_genome))) {
            stop("Sig vs ref_bkg_genome substitution labeling mismatch")
        }
        for (i_sig in 1:ncol(sig)) {
            sig[,i_sig] = sig[,i_sig]*ref_bkg_exome[,"Counts"]/ref_bkg_genome[,"Counts"]
            sig[,i_sig] = sig[,i_sig]/sum(sig[,i_sig])
        }
    } else if (sig.bkg.adj=="genome/exome") {
        ref_bkg_exome = read.table(system.file("extdata","bkg_exome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_exome = ref_bkg_exome[match(rownames(sig),rownames(ref_bkg_exome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_exome))) {
            stop("Sig vs ref_bkg_exome substitution labeling mismatch")
        }
        ref_bkg_genome = read.table(system.file("extdata","bkg_genome.txt",package="mutSigMapper"),header=TRUE,row.names=1,stringsAsFactors=FALSE,sep="\t")
        ref_bkg_genome = ref_bkg_genome[match(rownames(sig),rownames(ref_bkg_genome)),]
        if (!identical(rownames(sig),rownames(ref_bkg_genome))) {
            stop("Sig vs ref_bkg_genome substitution labeling mismatch")
        }
        for (i_sig in 1:ncol(sig)) {
            sig[,i_sig] = sig[,i_sig]*ref_bkg_genome[,"Counts"]/ref_bkg_exome[,"Counts"]
            sig[,i_sig] = sig[,i_sig]/sum(sig[,i_sig])
        }
    } else if (sig.bkg.adj=="custom") {
        sig.bkg.adj.custom = sig.bkg.adj.custom[match(rownames(sig),rownames(sig.bkg.adj.custom)),]
        if (!identical(rownames(sig),rownames(sig.bkg.adj.custom))) {
            stop("Sig vs sig.bkg.adj.custom substitution labeling mismatch")
        }
        for (i_sig in 1:ncol(sig)) {
            sig[,i_sig] = sig[,i_sig]*sig.bkg.adj.custom[,"Counts"]
            sig[,i_sig] = sig[,i_sig]/sum(sig[,i_sig])
        }
    }
    
    if (method=="MutationalPatterns") {
        mutSigMap = mutSigMapper_MutationalPatterns(spectra=spectra, sig=sig, n_rdm=n_rdm, noise=noise, neg.binom.size=neg.binom.size, signature.cutoff=signature.cutoff, sig.bkg.adj=sig.bkg.adj)
    } else if (method=="deconstructSigs") {
        mutSigMap = mutSigMapper_deconstructSigs(spectra=spectra, sig=sig, n_rdm=n_rdm, noise=noise, neg.binom.size=neg.binom.size, signature.cutoff=signature.cutoff, sig.bkg.adj=sig.bkg.adj)
    }
    mutSigMap$call <- match.call()
    class(mutSigMap) <- "mutSigMapper"
    if (save_obj) {
        save(mutSigMap,file=file.path(dest_dir,dest_obj))
    }
    mutSigMap
}
