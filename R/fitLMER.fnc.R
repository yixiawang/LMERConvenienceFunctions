fitLMER.fnc <-
function(
    model=as.character(),
    data=as.character(),
    backfit.on="F", #can also be "t"
    item=FALSE, # can be an item identifier such as "Item" or "Word"
    ran.effects=list(ran.intercepts=as.character(),
                     slopes=as.character(),
                     by.vars=as.character()),
    alpha=0.05,
    llrt=FALSE,
    p.value="Upper", # or "Lower"
    t.threshold=2,
    set.REML.FALSE=TRUE,
    reset.REML.TRUE=TRUE,
    file.name=paste("full_backfitLMER_log_",gsub(":","_",gsub(" ","_",date())),".txt",sep="") # or other path and file name or FALSE
    ){

  current.dir=getwd()
  temp.dir=tempdir()
  tempdir()
  setwd(temp.dir)

  if(file.name!=FALSE)sink(file=file.name,split=TRUE)  
  cat("======================================================\n")
  cat("===              backfitting fixed effects         ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=model,data=data,item=item,alpha=alpha,llrt=llrt,p.value=p.value,set.REML.FALSE=set.REML.FALSE,reset.REML.TRUE=FALSE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=model,data=data,item=item,alpha=alpha,llrt=llrt,t.threshold=t.threshold,set.REML.FALSE=set.REML.FALSE,reset.REML.TRUE=FALSE,log.file=FALSE)
  }


  cat("======================================================\n")
  cat("===            forwardfitting random effects       ===\n")
  cat("======================================================\n")
  mod=ffRanefLMER.fnc(model=mod,data=data,ran.effects=ran.effects,alpha=alpha)


  cat("======================================================\n")
  cat("===            re-backfitting fixed effects        ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=mod,data=data,item=FALSE,alpha=alpha,llrt=llrt,p.value=p.value,set.REML.FALSE=FALSE,reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=mod,data=data,item=FALSE,alpha=alpha,llrt=llrt,t.threshold=t.threshold,set.REML.FALSE=FALSE,reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }

  if(file.name!=FALSE){
    sink(file=NULL)
    cat("Log file saved in directory",temp.dir,"\n")
  }
  setwd(current.dir)

  return(model=mod)
}

