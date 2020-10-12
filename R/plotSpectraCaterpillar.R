plotSpectraCaterpillar = function(x, spectra.set=NULL, sig.set=NULL, order_sig_by_pval=TRUE, signif=TRUE, pval=0.05, CI=0.95, xlab=NULL, ylab=NULL, cexRow=NULL, cex.lab=0.95, main.title=NULL, col.title="black", line.title=1.5, cex.title=0.85, subtitle=NULL, col.subtitle="darkgray", line.subtitle=0.5, cex.subtitle=0.55, legend=TRUE)
{

    weights = x$weights
    weights_rdm = x$weights_rdm
    if ((CI<=0.5)|(CI>=1)) {
        stop("CI must be 0.5<CI<1")
    }
    n_sig = ncol(weights_rdm[[1]])
    weights_rdm_CI_template = matrix(rep(NA,3*n_sig),ncol=3)
    rownames(weights_rdm_CI_template) = colnames(weights_rdm[[1]])
    colnames(weights_rdm_CI_template) = c("lower_CI","median","upper_CI")

    if (is.null(spectra.set)) {
        sample.index = 1:length(weights)
    } else {
        sample.index = match(spectra.set,names(weights))
        sample.index = sample.index[!is.na(sample.index)]
    }
    for (i_sample in sample.index) {
        weights_rdm_CI = weights_rdm_CI_template
        for (i_sig in 1:n_sig) {
            weights_rdm_CI[i_sig,] = quantile(weights_rdm[[i_sample]][,i_sig],probs=c((1-CI)/2,0.5,(1+CI)/2))
        }
        labels = names(weights[[i_sample]])
        map_pval_sample = x$map_pval[,i_sample]
        if (signif) {
            annot = rep(" ",n_sig)
            annot[map_pval_sample<pval] = " *"
            labels = paste0(labels,annot)
            if (is.null(subtitle)) {
                subtitle = paste0("* p-value < ",pval)
            }
        }
        if (!is.null(sig.set)) {
            o = match(sig.set,names(weights[[i_sample]]))
            o = o[!is.na(o)]
            weights[[i_sample]] = weights[[i_sample]][o]
            weights_rdm_CI = weights_rdm_CI[o,]
            labels = labels[o]
            map_pval_sample = map_pval_sample[o]
        }
        if (order_sig_by_pval) {
            o = order(-map_pval_sample,weights[[i_sample]]) 
        } else {
            o = rev(1:length(labels))
        }
        weights[[i_sample]] = weights[[i_sample]][o]
        weights_rdm_CI = weights_rdm_CI[o,]
        labels = labels[o]
        if (is.null(xlab)) {
            if (x$method=="deconstructSigs") {
                xlab = "signature weight"
            } else if (x$method=="MutationalPatterns") {
                xlab = "signature relative contribution"
            }
        }
        if (is.null(ylab)) {
            ylab = ""
        }
        if (is.null(cexRow)) {
            cexRow = min(1.6/log(length(labels)),0.5)
        }
        if (is.null(main.title)) {
            if (x$method=="deconstructSigs") {
                main_title = paste0("Signature weights for sample ",names(weights)[i_sample])
            } else if (x$method=="MutationalPatterns") {
                main_title = paste0("Signature relative contributions for sample ",names(weights)[i_sample])
            }
        } else {
            main_title = main.title
        }

        x_range = range(weights[[i_sample]],weights_rdm_CI)
        y_range = c(1,length(labels))
        plot(x_range,y_range,type="n",xlab=xlab,ylab=ylab,yaxt="n",cex.lab=cex.lab)

        axis(side=2,at=1:length(labels),labels=labels,las=2,cex.axis=cexRow,tck=-0.005)
        abline(h=(1:length(labels)),col="lightgray",lty="dotted")
        title(main_title,col.main=col.title,line=line.title,cex.main=cex.title)
        title(subtitle,col.main=col.subtitle,line=line.subtitle,cex.main=cex.subtitle)
        points(weights_rdm_CI[,"median"],1:length(labels),pch=16,cex=1,col="blue",lty=3,type="p")
        suppressWarnings(
        arrows(weights_rdm_CI[,"lower_CI"],1:length(labels),weights_rdm_CI[,"upper_CI"],1:length(labels),
        col="blue",length=0.025,angle=90,code=3,lty=2)
        )
        points(weights[[i_sample]],1:length(labels),pch=4,cex=1,col="red",lty=3,type="p")
        if (legend) {
            legend("topright",inset=c(0,-0.165),legend=c("observed",paste0(round(CI*100),"% CI")),
            pch=c(4,16),col=c("red","blue"),xpd=TRUE,cex=0.75)
        }
    }
}
