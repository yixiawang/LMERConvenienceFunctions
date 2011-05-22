bfFixefLMER_t.fnc <-
function(
        model,
        data,
        item=FALSE, # otherwise, put between quotes an item identifier such as "Item" or "Word"
        alpha=0.05,
	llrt=FALSE, # or TRUE to have an extra step of log-likelihood ratio testing
        t.threshold=2,
        set.REML.FALSE=TRUE,
        reset.REML.TRUE=TRUE,
        log.file=file.path(tempdir(),paste("bfFixefLMER_t_log_",gsub(":","-",gsub(" ","_",date())),".txt",sep="")) # or other path and file name or FALSE
                ){
  #if(length(model)==0){
    #stop("please supply a value to the ''model'' argument")
  #}

  #if(length(data)==0){
    #stop("please supply a value to the ''data'' argument")
  #}

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

  if(llrt && set.REML.FALSE){
    cat("setting REML to FALSE\n")
    model=update(model,.~.,REML=FALSE)
  }

  current.dir=getwd()
  temp.dir=tempdir()
  tempdir()
  #setwd(temp.dir)

  options(warn=1)

  unlink(file.path(temp.dir,"temp.txt"))
  sink(file=NULL,type="message")   

  if(log.file!=FALSE)sink(file=log.file,split=TRUE)

  if(item!=FALSE){
    cat("checking",paste("by-",item,sep=""),"random intercepts\n")
    # Determine whether to include Items as a random effect
    model.updated<-NULL
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
  if(as.vector(model@call[1])=="glmer()"){
    odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
    data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
    temp.lmer=update(model,.~.,family="gaussian",data=data)
    coefficients=row.names(anova(temp.lmer))
    data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
  } else {
    coefficients=row.names(anova(model))
  }

  if(is(model,"mer")){
    model<-asS4(model)
    smry=as.data.frame(summary(model)@coefs)
  } else {
      stop("the input model is not a mer object\n")
  }

  # Keep only the level of an interaction term that has the greatest absolute t-value
  smry.temp=smry[-c(1:nrow(smry)),]
  for(coef in coefficients){
    intr.order=length(unlist(strsplit(coef,":")))
    orig.coef=coef
    coef=gsub(":","\\.\\*:",coef)
    coef=gsub("(^.*$)","\\1\\.\\*",coef)
    coef=gsub("\\(","\\\\(",coef)
    coef=gsub("\\)","\\\\)",coef)
    smry.temp2=smry[grep(paste("^",coef,sep=""),row.names(smry)),]
    smry.temp3=smry.temp2
    smry.temp3=smry.temp3[-c(1:nrow(smry.temp3)),]
    for(j in 1:nrow(smry.temp2)){
      if(length(unlist(strsplit(row.names(smry.temp2)[j],":")))==intr.order){
        smry.temp3=rbind(smry.temp3,smry.temp2[j,]) 
      }
    }
    smry.temp2=smry.temp3
    smry.temp2[,3]=abs(smry.temp2[,3])
    smry.temp2=smry.temp2[smry.temp2[,3]==max(smry.temp2[,3]),]
    row.names(smry.temp2)=orig.coef
    smry.temp=rbind(smry.temp,smry.temp2)
  }

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
    
    if(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),3]>=t.threshold){
      cat("\tall terms of interaction level",order,"significant\n")
    }

    while(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),3]<t.threshold){
      cat("\titeration",count,"\n")
      cat("\t\tt-value for term",paste("\"",row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),]),"\"",sep=""),"=",smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),3],"<",t.threshold,"\n")

      # Is the fixed-effect term part of a higher-order interaction? If it is, skip it and proceed with backfitting, otherwise evaluate it.
      if(length(grep(gsub("\\)","",gsub("\\(","",as.character(row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),])))),gsub("\\)","",gsub("\\(","",intr.order[intr.order$Order>order,"Coef"]))))>0){
        cat("\t\tpart of higher-order interaction\n")
        cat("\t\tskipping term\n")
        keepers=row.names(smry.temp2)
        keepers=keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),])),keepers)]
        smry.temp2=smry.temp2[keepers,]
      }else{
	cat("\t\tnot part of higher-order interaction\n")

        # Fit less complex nested model
	m.temp<-NULL
        eval(parse(text=paste("m.temp=update(model,.~.-",row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),]),")",sep="")))

	if(llrt){
        	# Which model should be kept? The more or the less complex one?
        	if(as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"])<=alpha){
          		cat("\t\tlog-likelihood ratio test p-value =",as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]),"<=",alpha,"\n")
          		cat("\t\tskipping term\n")
          		keepers=row.names(smry.temp2)
          		keepers=keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),])),keepers)]
          		smry.temp2=smry.temp2[keepers,]
        	} else {
          		cat("\t\tlog-likelihood ratio test p.value =",as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]),">",alpha,"\n")
          		cat("\t\tremoving term\n")
          		model=m.temp
  	
          		if(as.vector(model@call[1])=="glmer()"){
            			odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
            			data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
            			temp.lmer=update(model,.~.,family="gaussian",data=data)
            			coefficients=row.names(anova(temp.lmer))
            			data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
          		} else {
            			coefficients=row.names(anova(model))
          		}

  			smry2=as.data.frame(summary(asS4(model))@coefs)
    	
          		# Keep only the level of an interaction term that has the greatest absolute t-value
          		smry.temp2=smry2[-c(1:nrow(smry2)),]
          		for(coef in coefficients){
            			intr.order=length(unlist(strsplit(coef,":")))
            			orig.coef=coef
            			coef=gsub(":","\\.\\*:",coef)
            			coef=gsub("(^.*$)","\\1\\.\\*",coef)
            			coef=gsub("\\(","\\\\(",coef)
            			coef=gsub("\\)","\\\\)",coef)
            			smry.temp3=smry2[grep(paste("^",coef,sep=""),row.names(smry2)),]
            			smry.temp4=smry.temp3
            			smry.temp4=smry.temp4[-c(1:nrow(smry.temp4)),]
            			for(j in 1:nrow(smry.temp3)){
              				if(length(unlist(strsplit(row.names(smry.temp3)[j],":")))==intr.order){
                				smry.temp4=rbind(smry.temp4,smry.temp3[j,]) 
              				}
            			}
            			smry.temp3=smry.temp4
            			smry.temp3[,3]=abs(smry.temp3[,3])
            			smry.temp3=smry.temp3[smry.temp3[,3]==max(smry.temp3[,3]),]
            			row.names(smry.temp3)=orig.coef
            			smry.temp2=rbind(smry.temp2,smry.temp3)
          		}
        	
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
    	
          		smry.temp2$Order=intr.order$Order
    	
          		# Keep only terms with interaction order i
          		keepers=as.character(row.names(smry.temp2[smry.temp2$Order==order,]))
          		smry.temp2=smry.temp2[keepers,]
          	}
	  }else{
          	cat("\t\tremoving term\n")
          	model=m.temp
  	
          	if(as.vector(model@call[1])=="glmer()"){
            		odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
            		data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
            		temp.lmer=update(model,.~.,family="gaussian",data=data)
            		coefficients=row.names(anova(temp.lmer))
            		data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
          	} else {
            		coefficients=row.names(anova(model))
          	}

  		smry2=as.data.frame(summary(asS4(model))@coefs)
    	
          	# Keep only the level of an interaction term that has the greatest absolute t-value
          	smry.temp2=smry2[-c(1:nrow(smry2)),]
          	for(coef in coefficients){
            		intr.order=length(unlist(strsplit(coef,":")))
            		orig.coef=coef
            		coef=gsub(":","\\.\\*:",coef)
            		coef=gsub("(^.*$)","\\1\\.\\*",coef)
            		coef=gsub("\\(","\\\\(",coef)
            		coef=gsub("\\)","\\\\)",coef)
            		smry.temp3=smry2[grep(paste("^",coef,sep=""),row.names(smry2)),]
            		smry.temp4=smry.temp3
            		smry.temp4=smry.temp4[-c(1:nrow(smry.temp4)),]
            		for(j in 1:nrow(smry.temp3)){
              			if(length(unlist(strsplit(row.names(smry.temp3)[j],":")))==intr.order){
                			smry.temp4=rbind(smry.temp4,smry.temp3[j,]) 
              			}
            		}
            		smry.temp3=smry.temp4
            		smry.temp3[,3]=abs(smry.temp3[,3])
            		smry.temp3=smry.temp3[smry.temp3[,3]==max(smry.temp3[,3]),]
            		row.names(smry.temp3)=orig.coef
            		smry.temp2=rbind(smry.temp2,smry.temp3)
          	}
        
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
    
          	smry.temp2$Order=intr.order$Order
    	
          	# Keep only terms with interaction order i
          	keepers=as.character(row.names(smry.temp2[smry.temp2$Order==order,]))
          	smry.temp2=smry.temp2[keepers,]
            }
	  }
        count=count+1
        if(nrow(smry.temp2)==0){
          break
        }
      }
  }

  if(llrt && reset.REML.TRUE){
    cat("resetting REML to TRUE\n")
    model=update(model,.~.,REML=TRUE)
  }

  if(log.file!=FALSE){
    sink(file=NULL)
    cat("Log file saved in directory",temp.dir,"\n")
  }
  #setwd(current.dir)

  return(model=model)
}

