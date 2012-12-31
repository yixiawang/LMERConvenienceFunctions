fitLMER.fnc <- function(model=as.character(),
    backfit.on="F", #can also be "t"
    item=FALSE, # can be an item identifier such as "Item" or "Word"
    ran.effects=list(ran.intercepts=as.character(),
                     slopes=as.character(),
		     corr=as.character(),
                     by.vars=as.character()),
    alpha=0.05,
    if.warn.not.add=TRUE,
    llrt=FALSE,
    prune.ranefs=TRUE, # remove or not random effects for which the variable is not in the fixed effects structure. 
    p.value="upper", # or "lower"
    t.threshold=2,
    set.REML.FALSE=TRUE,
    reset.REML.TRUE=TRUE,
    log.file.name=file.path(tempdir(),paste("fitLMER_log_",gsub(":","-",
	gsub(" ","_",date())),".txt",sep="")) # or other path and file name or FALSE
    ){

  current.dir=getwd()
  temp.dir=tempdir()
  tempdir()

  if(log.file.name!=FALSE)sink(file=log.file.name,split=TRUE)  
  cat("======================================================\n")
  cat("===              backfitting fixed effects         ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=model,item=item,alpha=alpha,llrt=llrt,
	prune.ranefs=prune.ranefs,p.value=p.value,set.REML.FALSE=set.REML.FALSE,
	reset.REML.TRUE=FALSE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=model,item=item,alpha=alpha,llrt=llrt,
	prune.ranefs=prune.ranefs,t.threshold=t.threshold,set.REML.FALSE=set.REML.FALSE,
	reset.REML.TRUE=FALSE,log.file=FALSE)
  }


  cat("======================================================\n")
  cat("===            forwardfitting random effects       ===\n")
  cat("======================================================\n")
  mod=ffRanefLMER.fnc(model=mod,ran.effects=ran.effects,alpha=alpha,if.warn.not.add=if.warn.not.add,log.file=FALSE)


  cat("======================================================\n")
  cat("===            re-backfitting fixed effects        ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=mod,item=FALSE,alpha=alpha,llrt=llrt,
	prune.ranefs=prune.ranefs,p.value=p.value,set.REML.FALSE=FALSE,
	reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=mod,item=FALSE,alpha=alpha,llrt=llrt,
	prune.ranefs=prune.ranefs,t.threshold=t.threshold,set.REML.FALSE=FALSE,
	reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }

  options(warn=0)
  sink(file=NULL,type="message")        

  if(log.file.name!=FALSE){
    cat("log file is",log.file.name,"\n")
    sink(file=NULL)
  }

  return(model=mod)
}

