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
    t.threshold=2,
    set.REML.FALSE=TRUE,
    reset.REML.TRUE=TRUE,
    file.name=paste("./full_backfitLMER_log_",gsub(" ","_",date()),".txt",sep="") # can also be FALSE to disable
    ){
  if(file.name!=FALSE)sink(file=file.name,split=TRUE)  
  cat("======================================================\n")
  cat("===              backfitting fixed effects         ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=model,data=data,item=item,alpha=alpha,set.REML.FALSE=set.REML.FALSE,reset.REML.TRUE=FALSE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=model,data=data,item=item,alpha=alpha,t.threshold=t.threshold,set.REML.FALSE=set.REML.FALSE,reset.REML.TRUE=FALSE,log.file=FALSE)
  }
  cat("saving model to file",paste(getwd(),"/","temp_bfFixefLMER.RData",sep=""),"\n")
  save(mod,file="temp_bfFixefLMER.RData")
  cat("======================================================\n")
  cat("===            forwardfitting random effects       ===\n")
  cat("======================================================\n")
  mod=ffRanefLMER.fnc(model=mod,data=data,ran.effects=ran.effects,alpha=alpha)
  cat("saving model to file",paste(getwd(),"/","temp_ffRanefLMER.RData",sep=""),"\n")
  save(mod,file="temp_ffRanefLMER.RData")
  cat("======================================================\n")
  cat("===            re-backfitting fixed effects        ===\n")
  cat("======================================================\n")
  if(backfit.on=="F"){
    mod=bfFixefLMER_F.fnc(model=mod,data=data,item=FALSE,alpha=alpha,set.REML.FALSE=FALSE,reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }else{
    mod=bfFixefLMER_t.fnc(model=mod,data=data,item=FALSE,alpha=alpha,t.threshold=t.threshold,set.REML.FALSE=FALSE,reset.REML.TRUE=reset.REML.TRUE,log.file=FALSE)
  }
  cat("saving model to file",paste(getwd(),"/","temp_re_bfFixefLMER.RData",sep=""),"\n")
  save(mod,file="temp_re_bfFixefLMER.RData")
  if(file.name!=FALSE)sink(file=NULL)
  return(model=mod)
}

