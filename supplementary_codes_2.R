N_GE <- function(Kernel,         # some genomic, enviromic or enviromic x genomic kernel
                 fraction =.98,  # expected fraction explained (Misztal, 2016)
                 plot=TRUE,      # plot  svd?
                 svd.print=FALSE, # save svd?
                 title = 'Effective N size',
                 ylab = 'Fraction of Kernel Variance',
                 xlab= 'Individuals'
                 ){
  svd.KERNEL <- svd(Kernel, nu = nrow(Kernel), nv = ncol(Kernel))
  SVdD<-cumsum((svd.KERNEL$d[1:ncol(Kernel)])^2/sum(svd.KERNEL$d^2))
  N <- length(SVdD[which(SVdD < fraction)])
  if(isTRUE(plot)){
    plot(SVdD,lwd=2,pch=15,xlab=xlab,bty='l', cex.lab=1.5, cex.axis=1.5, cex.main=2.2,
         main=title,ylab=ylab)
    lines(SVdD,lwd=2)
    abline(v = N,lwd=3,col='red')
  }
  
  cat(paste0('--------------------------------------','\n'))
  cat(paste0('Explained variance ',100*fraction,'% \n'))
  cat(paste0('Effective Ne:',N,'\n'))
  cat(paste0('--------------------------------------','\n'))
  if(isFALSE(svd.print)) return(list(Ne=N))
  if(isTRUE(svd.print)) return(list(Ne=N,svd.d = svd.KERNEL$d,svd.u=svd.KERNEL$u,svd.v=svd.KERNEL$v))
}
