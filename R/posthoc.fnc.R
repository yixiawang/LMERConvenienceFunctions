posthoc.fnc<-function(
    model,   #unquoted name of model
    factor, # quoted name of factor
    ndigits = 4,
    MCMC=FALSE,
    nsim=10000,
    addPlot=FALSE
){
  if(MCMC){
    pdf(file="posthoc.fnc_MCMC_plots.pdf")
  }

  data<-model@frame

  levs=sort(levels(data[,factor]))
  posthoc=list()
  for(lev in levs){
    cat(paste("processing level",paste("''",lev,"''",sep=""),"(",grep(lev,levs),"of",length(levs),") ...\n"),sep=" ")
    data[,factor]=relevel(data[,factor],lev)
    model=update(model,.~.,data=data)
    if(class(model)=="lm"){
      posthoc[[lev]]=round(summary(model)$coefficients,ndigits)
    } else {
      posthoc[[lev]]=round(summary(model)@coefs,ndigits)
      if(MCMC){
        model.mcmc=pvals.fnc(model,nsim,ndigits=ndigits,withMCMC=TRUE,addPlot=addPlot)
        posthoc[[paste(lev,"MCMC",sep="")]]=model.mcmc$fixed
        posthoc[[paste(lev,"mcmcMat",sep="")]]=model.mcmc$mcmc
      }
    }
  }
  cat("\n")
  if(MCMC){
    dev.off()
  }
  return(posthoc)
}


