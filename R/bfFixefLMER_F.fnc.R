bfFixefLMER_F.fnc <-
function(
        model=as.character(),
        data=as.character(),
        item=FALSE, # can be an item identifier such as "Item" or "Word"
        alpha=0.05,
        set.REML.FALSE=TRUE,
        reset.REML.TRUE=TRUE,
        log.file=paste("./fixef_backfit_log_",gsub(" ","_",date()),".txt",sep="") # or other path and file name or FALSE
                ){
  if(length(model)==0){
    stop("please supply a value to the ''model'' argument")
  }

  if(length(data)==0){
    stop("please supply a value to the ''data'' argument")
  }

  if(length(item)==0){
    stop("please supply a value to the ''item'' argument")
  }

  if(length(alpha)==0){
    stop("please supply a value to the ''alpha'' argument")
  }

  if(length(set.REML.FALSE)==0){
    stop("please supply a value to the ''set.REML.FALSE'' argument")
  }

  if(length(reset.REML.TRUE)==0){
    stop("please supply a value to the ''reset.REML.TRUE'' argument")
  }

  if(as.vector(model@call[1])=="glmer()"){
    stop("glmer models not yet supported")
  }

  if(set.REML.FALSE){
    cat("setting REML to FALSE\n")
    model=update(model,.~.,REML=FALSE)
  }

  options(warn=1)

  unlink("temp.txt")
  sink(file=NULL,type="message")   

  if(log.file!=FALSE)sink(file=log.file,split=TRUE)

  if(item!=FALSE){
    cat("checking",paste("by-",item,sep=""),"random intercepts\n")
    # Determine whether to include Items as a random effect
    eval(parse(text=paste("model.updated=update(model,.~.+(1|",item,"))",sep="")))
    if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
      cat("\tadding",paste("by-",item,sep=""),"random intercepts to model\n")
      model=model.updated
    } else {
      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
      cat("\tnot adding",paste("by-",item,sep=""),"random intercepts to model\n")
    }
  }

  # Get initial model coefficients and summary
  coefficients=row.names(pamer.fnc(model))
  smry=pamer.fnc(model)
  smry.temp=pamer.fnc(model)

  # Determine the interaction orders for each term
  temp=strsplit(coefficients,":")
  names(temp)=coefficients
  intr.order=list()
  for(i in coefficients){
    intr.order[[i]]=length(temp[[i]])
  }
  intr.order=as.data.frame(unlist(intr.order))
  colnames(intr.order)="Order"
  intr.order$Coef=row.names(intr.order)
  row.names(intr.order)=1:nrow(intr.order)

  orders=sort(as.numeric(unique(intr.order$Order)),decreasing=TRUE) # decreasing

  smry.temp$Order=intr.order$Order
  count=1

  for(order in orders){
    cat("processing model terms of interaction level",order,"\n")
    # Keep only terms with interaction order i
    keepers=as.character(row.names(smry.temp[smry.temp$Order==order,]))
    smry.temp2=smry.temp[keepers,]
    
    if(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),"Upper-bound p-value"]<=alpha){
      cat("\tall terms of interaction level",order,"significant\n")
    }else{

      while(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),"Upper-bound p-value"][1]>alpha){
        cat("\titeration",count,"\n")
	evaluated.term=row.names(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),])
        cat("\t\tp-value for term",paste("\"",evaluated.term,"\"",sep=""),"=",smry.temp2[evaluated.term,"Upper-bound p-value"],">",alpha,"\n")
  
        # Is the fixed-effect term part of a higher-order interaction? If it is, skip it and proceed with backfitting, otherwise evaluate it.
        if(length(grep("FALSE",unique(unlist(strsplit(evaluated.term,":"))%in%unlist(strsplit(intr.order[intr.order$Order>order,"Coef"],":")))))<1){
          cat("\t\tpart of higher-order interaction\n")
          cat("\t\tskipping term\n")
          keepers=row.names(smry.temp2)
          keepers=keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),])),keepers)]
          smry.temp2=smry.temp2[keepers,]
          smry.temp2=na.omit(smry.temp2)
        } else {
          cat("\t\tnot part of higher-order interaction\n")
  
          # Fit less complex nested model
          eval(parse(text=paste("m.temp=update(model,.~.-",row.names(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),]),")",sep="")))
  
          # Which model should be kept? The more or the less complex one?
          if(as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"])<=alpha){
            cat("\t\tlog-likelihood ratio test p-value =",as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]),"<=",alpha,"\n")
            cat("\t\tskipping term\n")
            keepers=row.names(smry.temp2)
            keepers=keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[,"Upper-bound p-value"]==max(smry.temp2[,"Upper-bound p-value"]),])),keepers)]
            smry.temp2=smry.temp2[keepers,]
          } else {
            cat("\t\tlog-likelihood ratio test p.value =",as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]),">",alpha,"\n")
            cat("\t\tremoving term\n")
            model=m.temp
    
            # Get initial model coefficients and summary
            coefficients=row.names(pamer.fnc(model))
            smry=pamer.fnc(model)
            smry.temp=pamer.fnc(model)
          
            # Determine the interaction orders for each term
            temp=strsplit(coefficients,":")
            names(temp)=coefficients
            intr.order=list()
            for(i in coefficients){
              intr.order[[i]]=length(temp[[i]])
            }
            intr.order=as.data.frame(unlist(intr.order))
            colnames(intr.order)="Order"
            intr.order$Coef=row.names(intr.order)
            row.names(intr.order)=1:nrow(intr.order)
      
            smry.temp$Order=intr.order$Order
  
            # Keep only terms with interaction order i
            keepers=as.character(row.names(smry.temp[smry.temp$Order==order,]))
            smry.temp2=smry.temp2[keepers,]
            smry.temp2=na.omit(smry.temp2)
          }
        }
        count=count+1
        if(nrow(smry.temp2)==0){
          break
        }
      }
    }
  }
  if(reset.REML.TRUE){
    cat("resetting REML to TRUE\n")
    model=update(model,.~.,REML=TRUE)
  }

  if(log.file!=FALSE)sink(file=NULL)

  return(model=model)
}

