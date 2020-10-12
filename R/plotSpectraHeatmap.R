plotSpectraHeatmap = function(x, spectra.set=NULL, sig.set=NULL, type=c("weights","cosine","entropy","correlation"), cor.method = c("pearson","kendall","spearman"), signif=TRUE, pval=0.05, signif.sig.only=FALSE, signif.sample.only=FALSE, xlab=NULL, ylab=NULL, color=NULL, breaks=NULL, breaks.scale.adj=TRUE, breaks.scale.adj.max=c("q99","q95","max"), scale="none", Rowv=FALSE, Colv=TRUE, na.color=NULL, cexRow = NULL, srtRow=0, cexCol=NULL, srtCol=45, margins=c(5,5), key=TRUE, key.title=NA, dendrogram="none", trace="none", notecol = "red", cor.notecol = "white", notecex = 1.5, main.title=NULL, col.title="black", line.title=1, cex.title=0.95, subtitle=NULL, col.subtitle="red", line.subtitle=-1, cex.subtitle=0.65, ...)
{
    weights = x$weights
    weights_rdm = x$weights_rdm
    map_pval = x$map_pval
    sample_label = names(weights)
    n_sample = length(sample_label)
    sig_label = colnames(weights_rdm[[1]])
    n_sig = length(sig_label)
    type = match.arg(type)
    if (type=="weights") {
        heatmap_mat = t(matrix(unlist(weights),byrow=TRUE,ncol=n_sig))
    } else if (type=="cosine") {
        ### from https://stats.stackexchange.com/questions/31565/is-there-an-r-function-that-will-compute-the-cosine-dissimilarity-matrix
        cos.sim=function(ma, mb){
            mat=tcrossprod(ma, mb)
            t1=sqrt(apply(ma, 1, crossprod))
            t2=sqrt(apply(mb, 1, crossprod))
            mat / outer(t1,t2)
        }
        heatmap_mat = cos.sim(t(x$sig),t(x$spectra))
    } else if (type=="correlation") {
        cor.method = match.arg(cor.method)
        heatmap_mat = cor(x$sig,x$spectra,method=cor.method)
    } else if (type=="entropy") {
        epsilon = 1.e-6
        sig_norm = t(t(x$sig)/colSums(x$sig))
        sig_norm[which(sig_norm<epsilon,arr.ind=TRUE)] = epsilon
        spectra_norm = t(t(x$spectra)/colSums(x$spectra))
        spectra_norm[which(spectra_norm<epsilon,arr.ind=TRUE)] = epsilon
        KLD = function(x,y) sum(x*log2(x/y)) # Kullback-Leibler divergence
        JSD = function(x,y) { # Jensen-Shannon divergence
            z = (x+y)/2
            (KLD(x,z)+KLD(y,z))/2
        }
        heatmap_mat = matrix(rep(NA,n_sig*n_sample),ncol=n_sample)
        for (i_sig in 1:n_sig) {
            for (i_sample in 1:n_sample) {
                heatmap_mat[i_sig,i_sample] = JSD(sig_norm[,i_sig],spectra_norm[,i_sample])
            }
        }
    }
    rownames(heatmap_mat) = sig_label
    colnames(heatmap_mat) = sample_label
    if (!is.null(sig.set)) {
        sig.index = match(sig.set,rownames(heatmap_mat))
        sig.index = sig.index[!is.na(sig.index)]
        heatmap_mat = heatmap_mat[sig.index,]
        map_pval = map_pval[sig.index,]
        n_sig = nrow(heatmap_mat)
    }
    if (!is.null(spectra.set)) {
        sample.index = match(spectra.set,colnames(heatmap_mat))
        sample.index = sample.index[!is.na(sample.index)]
        heatmap_mat = heatmap_mat[,sample.index]
        map_pval = map_pval[,sample.index]
        n_sample = ncol(heatmap_mat)
    }
    breaks.scale.adj.max = match.arg(breaks.scale.adj.max)
    if (type%in%c("weights","cosine","entropy")) {
        if (is.null(breaks)) {
            if (breaks.scale.adj) {
                if (breaks.scale.adj.max=="q99") {
                    scale_max = quantile(heatmap_mat,probs=0.99)
                } else if (breaks.scale.adj.max=="q95") {
                    scale_max = quantile(heatmap_mat,probs=0.95)
                } else if (breaks.scale.adj.max=="max") {
                    scale_max = max(heatmap_mat)
                }
                breaks = seq(0,scale_max,by=scale_max/10)
            } else {
                breaks = seq(0,1,by=0.1)
            }
        }
        if (is.null(color)) {
            if (type%in%c("weights","cosine")) {
                color = suppressWarnings(colorRampPalette(brewer.pal(length(breaks)-1,"YlGnBu")))
            } else {
                color = suppressWarnings(colorRampPalette(rev(brewer.pal(length(breaks)-1,"YlGnBu"))))
            }
        }
        if (is.null(na.color)) {
            na.color = "black"
        }
    } else if (type=="correlation") {
        if (is.null(breaks)) {
            if (breaks.scale.adj) {
                if (breaks.scale.adj.max=="q99") {
                    scale_max = quantile(abs(heatmap_mat),probs=0.99)
                } else if (breaks.scale.adj.max=="q95") {
                    scale_max = quantile(abs(heatmap_mat),probs=0.95)
                } else if (breaks.scale.adj.max=="max") {
                    scale_max = max(abs(heatmap_mat))
                }
                breaks = seq(0,scale_max,by=scale_max/5)
                breaks = unique(c(-breaks,breaks))
                breaks = breaks[order(breaks)]
            } else {
                breaks = seq(-1,1,by=0.2)
            }
        }
        if (is.null(color)) {
            color = redgreen
        }
        if (is.null(na.color)) {
            na.color = "white"
        }
        notecol = cor.notecol
    }
    if (is.null(cexRow)) {
        cexRow = min(1.6/log(n_sig),0.5)
    }
    if (is.null(cexCol)) {
        cexCol = min(1.6/log(n_sample),0.5)
    }
    if (is.null(main.title)) {
        if (type=="weights") {
            main.title = "Signature vs spectra: weights matrix"
        } else if (type=="cosine") {
            main.title = "Signature vs spectra: cosine similarity matrix"
        } else if (type=="correlation") {
            main.title = paste0("Signature vs spectra: correlation (",cor.method,") matrix")
        } else if (type=="entropy") {
            main.title = "Signature vs spectra: Jensen-Shannon divergence matrix"
        }
    }
    if (signif) {
        annot = matrix(rep("",n_sig*n_sample),ncol=n_sample)
        annot[which(map_pval<pval,arr.ind=TRUE)] = "*"
        if (signif.sig.only) {
            row_index = sort(unique(which(annot=="*",arr.ind=TRUE)[,1]))
            heatmap_mat = heatmap_mat[row_index,,drop=FALSE]
            annot = annot[row_index,,drop=FALSE]
        }
        if (signif.sample.only) {
            col_index = sort(unique(which(annot=="*",arr.ind=TRUE)[,2]))
            heatmap_mat = heatmap_mat[,col_index,drop=FALSE]
            annot = annot[,col_index,drop=FALSE]
        }
        if (min(dim(heatmap_mat))>1) {
            heatmap.2(heatmap_mat, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=color,
            breaks = breaks, symkey=FALSE, dendrogram = dendrogram, margins=margins, cexRow=cexRow,
            srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace, key.title=key.title,
            xlab= xlab, ylab = ylab,
            cellnote = annot, notecol = notecol, notecex = notecex, ...)
            title(main=main.title, col.main=col.title, line=line.title, cex.main=cex.title)
            if (is.null(subtitle)) {
                subtitle = paste0("* p-value < ",pval)
            }
            title(main=subtitle,col.main=col.subtitle, line=line.subtitle, cex.main=cex.subtitle)
        } else {
            warning("Number of signatures and/or samples selected is < 2")
        }
    } else {
        if (min(dim(heatmap_mat))>1) {
            heatmap.2(heatmap_mat, scale=scale, Rowv=Rowv, Colv=Colv, na.color=na.color, col=color,
            breaks = breaks, symkey=FALSE, dendrogram = dendrogram, margins=margins, cexRow=cexRow,
            srtRow=srtRow, cexCol=cexCol, srtCol=srtCol, key=key, trace=trace, key.title=key.title,
            xlab= xlab, ylab = ylab, ...)
            title(main=main.title, col.main=col.title, line=line.title, cex.main=cex.title)
            if (!is.null(subtitle)) {
                title(main=subtitle,col.main=col.subtitle, line=line.subtitle, cex.main=cex.subtitle)
            }
        } else {
            warning("Number of signatures and/or samples selected is < 2")
        }
    }
    list(sig_spectra=as.data.frame(heatmap_mat))
}
