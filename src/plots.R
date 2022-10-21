###########################################################################
# Plots
###########################################################################

plot_path <- function(gamma,beta,names=NULL){

  d <- data.frame(beta)
  if(!is.null(names)){
    colnames(d) <- names
  }
  d$gamma <- gamma
  p <- ggplot(melt(d,id.vars = c('gamma')),
              aes(x=gamma,y=value,color=variable)) +
    geom_line() +
    scale_color_discrete(name = "") +
    scale_x_continuous(breaks=seq(0,1,0.25),
                       labels = expression(hat(beta)[CD],
                                           '0.25',
                                           '0.5',
                                           '0.75',
                                           hat(beta)[OLS])) +
    ylab('Coefficient') +
    xlab(expression('Regularization path' ~ hat(beta)[gamma])) +
    ggtitle('Regularization path for causal regularization')

  return(p)

}

# out-of-sample simulation plot
out_plot_for<-function(target,B,n,init_shift,gamma,beta,
                       absdi=NULL,
                       sdi=NULL,
                       iri=NULL,
                       infdi=NULL,
                       name=NULL){
  shift<-seq(100,500,100)
  on <- 10000 # out-of-sample number of examples
  ors<-list()
  ranges<-NULL
  j<-1
  for (s in shift){
    
    odata<-gen(target,B,on,0,s) # out-of-sample data
    
    or <- NULL
    for (i in seq(1,dim(beta)[1],1)){
      b<-t(t(beta[i,]))
      or<-c(or,risk(odata$X,odata$y, b)/s)
    }
    
    ranges<-c(ranges,range(or))
    ors[[j]]<-or
    
    j<-j+1
  }
  
  r = c(min(ranges),max(ranges))
  
  out_plot <- function(){
    plot(gamma,ors[[1]],type = "l",
         ylab="normalized out-of-sample risk",ylim=r,
         main=paste(name,", n =",n,",","shift = ",init_shift,sep= " "))
    
    for(j in seq_along(shift)[-1]){
      lines(gamma,ors[[j]],lty=1)  
    }
    
    if(!is.null(sdi)){
      lines(c(gamma[absdi],gamma[absdi]),r,col="green")
      lines(c(gamma[sdi],gamma[sdi]),r,col="red")
      lines(c(gamma[iri],gamma[iri]),r,col="blue")  
      lines(c(gamma[infdi],gamma[infdi]),r,col="magenta")
    }
    
  }
  
  return(out_plot)
}

metric_plot <- function(gamma,diffs,metric,color,prefix,suffix){
  
  m <- diffs[[metric]]
  m_se <- diffs[[paste(metric,'se',sep='_')]]
  rg <- range(m)
  
  if(!is.null(m_se)){
    m_upper <- m + m_se
    m_lower <- m - m_se
    
    min_y <- min(m_lower)
    max_y <- max(m_upper)
    rg <- c(min_y,max_y)
    
    plot(gamma,m,type = "b",lty=1,ylim=rg,ylab="")
    mtext(paste(prefix,"in-sample",suffix), side=2, line=2.2, cex=1.5)
    
    lines(gamma,m_upper,type = "l",lty=2)
    lines(gamma,m_lower,type = "l",lty=2)  
    
    folds <- diffs[['folds']]
    for (i in seq(1,length(folds)) ) {
      lines(gamma,folds[[i]][[metric]],type = "b",lty=1,
            col=scales::alpha("black", 0.1))
    }
    
  }else{
    
    plot(gamma,m,type = "l",
         ylab=paste(prefix,"in-sample absolute difference"))
    
  }
  
  m_index <- diffs[[paste(metric,'i',sep='_')]]
  
  if(!is.null(m_index)){
    lines(c(gamma[m_index],gamma[m_index]),rg,col=color)  
  }
  
}

trade_off_plot<- function(gamma,diffs,prefix,name){
  
  if(is.null(prefix)){
    prefix<-""
  }
  
  par(mfrow=c(1,4))
  metric_plot(gamma = gamma, diffs = diffs,metric = 'absd',
              color = 'green',prefix=prefix,suffix='absolute difference')
  metric_plot(gamma = gamma, diffs = diffs,metric = 'sd',
              color = 'red',prefix=prefix,suffix='squared difference')
  metric_plot(gamma = gamma, diffs = diffs,metric = 'infd',
              color = 'magenta',prefix=prefix,suffix='infinity difference')
  metric_plot(gamma = gamma, diffs = diffs,metric = 'ir',
              color = 'blue',prefix=prefix,suffix='pool risk')
  
  mtext(name, outer=TRUE,  cex=1.5, line=-2)
  
}

plots <- function(target,B,n,init_shift,diffs,gamma,beta,prefix,name,simulation=FALSE){
  
  # path plot
  par(mfrow=c(1,1))
  plot_path(gamma,beta,name,
            absdi=diffs[['absd_i']],
            sdi=diffs[['sd_i']],
            iri=diffs[['ir_i']],
            infdi=diffs[['infd_i']])
  
  # trade-off plot
  trade_off_plot(gamma = gamma,diffs=diffs,name=name,prefix = prefix)
  
  # out-of-sample risk plot
  if(simulation){
    out_plot <- out_plot_for(target,B,n,init_shift,gamma,beta,
                             absdi=diffs[['absd_i']],
                             sdi=diffs[['sd_i']],
                             iri=diffs[['ir_i']],
                             infdi=diffs[['infd_i']],
                             name = name)  
    par(mfrow=c(1,1))
    out_plot()
  }
  
}