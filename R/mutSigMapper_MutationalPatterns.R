mutSigMapper_MutationalPatterns = function(spectra, sig, n_rdm, noise, neg.binom.size,  signature.cutoff, sig.bkg.adj)
{
    if (noise=="poisson") {
        add_noise = function(n_rdm,mu,size) {
            rpois(n=n_rdm,lambda=mu)
        }
    } else if (noise=="neg.binom") {
        add_noise = function(n_rdm,mu,size) {
            rnbinom(n=n_rdm,mu=mu,size=size)
        }
    }
    spectra_input = spectra
    sig_input = sig
    spectra = as.matrix(spectra)
    channel = rownames(spectra)
    n_channel = length(channel)
    sample_label = colnames(spectra)
    n_sample = length(sample_label)
    sig = as.matrix(sig)
    sig_label = colnames(sig)
    n_sig = length(sig_label)
    sample_rdm_template = matrix(rep(0,n_channel*n_rdm),ncol=n_rdm)
    weights_rdm_template = matrix(rep(0,n_rdm*n_sig),ncol=n_sig)
    rownames(weights_rdm_template) = paste0("rdm_",1:n_rdm)
    colnames(weights_rdm_template) = sig_label
    weights = vector("list",n_sample)
    names(weights) = sample_label
    weights_rdm = vector("list",n_sample)
    names(weights_rdm) = sample_label
    map_pval = matrix(rep(NA,n_sig*n_sample),ncol=n_sample)
    rownames(map_pval) = sig_label
    colnames(map_pval) = sample_label
    for (i_sample in 1:n_sample) {
        sample = spectra[,i_sample,drop=FALSE]
        index = which(sample>0)
        n_index = length(index)
        if (n_index<1) {
            weights[[i_sample]] = rep(0,n_sig)
            names(weights[[i_sample]]) = sig_label
            weights_rdm[[i_sample]] = weights_rdm_template
        } else {
            weights[[i_sample]] = as.numeric(fit_to_signatures(sample,sig)$contribution)
            names(weights[[i_sample]]) = sig_label
            sample_rdm = sample_rdm_template
            for (i_index in 1:n_index) {
                sample_rdm[index[i_index],] = add_noise(n_rdm,as.numeric(sample[index[i_index]]),neg.binom.size)
            }
            weights_rdm[[i_sample]] = weights_rdm_template
            for (i_rdm in 1:n_rdm) {
                if (sum(sample_rdm[,i_rdm])>0) {
                    weights_rdm[[i_sample]][i_rdm,] = as.numeric(fit_to_signatures(sample_rdm[,i_rdm,drop=FALSE],sig)$contribution)
                }
            }
        }
    }
    # we normalize the output contributions
    counts = weights
    counts_rdm = weights_rdm
    for (i_sample in 1:n_sample) {
        weights[[i_sample]] = weights[[i_sample]]/max(sum(weights[[i_sample]]),1)
        weights_rdm[[i_sample]] = weights_rdm[[i_sample]]/pmax(rowSums(weights_rdm[[i_sample]]),1)
        map_pval[,i_sample] = (apply(weights_rdm[[i_sample]]<=signature.cutoff,2,sum)+1)/(n_rdm+1)
    }
    list(spectra=spectra_input, sig=sig_input, method="MutationalPatterns", sig.bkg.adj=sig.bkg.adj, 
    signature.cutoff=signature.cutoff, noise=noise, neg.binom.size=neg.binom.size,
    n_rdm=n_rdm, weights=weights, weights_rdm=weights_rdm, map_pval=map_pval,
    counts=counts, counts_rdm=counts_rdm)
}
