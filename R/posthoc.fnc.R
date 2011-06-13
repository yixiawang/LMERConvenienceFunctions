posthoc.fnc<-function(
    model,   #unquoted name of model
    factor, # quoted name of factor
    two.tailed=TRUE, # or FALSE for 1-tailed
    num.comp=NULL, # an integer for the number of comparisons
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
      temp<-round(summary(model)$coefficients,ndigits)
      posthoc[[lev]]<-temp
    } else {
      temp<-as.data.frame(round(summary(model)@coefs,ndigits))

      dims<-NULL
      anova.table=anova(model)
      rank.X=qr(model@X)$rank

      # get upper-bound df
      udf<-nrow(model@frame)-rank.X

      # get lower-bound df
      model.ranef<-ranef(model)
      lower.bound<-0
      for(i in 1:length(names(model.ranef))){
	      dims<-dim(model.ranef[[i]])
	      lower.bound<-lower.bound+dims[1]*dims[2]
      }
      ldf<-nrow(model@frame)-rank.X-lower.bound
     
      multip<-ifelse(two.tailed,2,1)
      temp[,"udf"]<-udf
      temp[,"ldf"]<-ldf
      temp[,"up.p.unadj."]<-round(multip*(1-pt(abs(temp[,"t value"]),udf)),ndigits)
      temp[,"low.p.unadj."]<-round(multip*(1-pt(abs(temp[,"t value"]),ldf)),ndigits)
      
      if(!is.null(num.comp)){
           temp[,"up.p.adj."]<-round(num.comp*(multip*(1-pt(abs(temp[,"t value"]),udf))),ndigits)
           temp[,"low.p.adj."]<-round(num.comp*(multip*(1-pt(abs(temp[,"t value"]),ldf))),ndigits)
      }

      posthoc[[lev]]<-temp

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


