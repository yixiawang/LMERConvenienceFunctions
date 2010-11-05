ffRanefLMER.fnc <-
function(
          model=as.character(),
          data=as.character(),
          ran.effects=list(ran.intercepts=as.character(), # or can specify a vector
                           slopes=as.character(), # of random effects to consider, e.g.,
                           by.vars=as.character()), # c("(0+Length|Subject)","(1+Frequency|Subject)")
          alpha=0.05,
          log.file=paste("./ranef_forwardfit_log_",gsub(" ","_",date()),".txt",sep="") # or other path and file name or FALSE
          ){
  unlink("temp.txt")
  sink(file=NULL,type="message")   
  
  if(length(model)==0){
    stop("please supply a value to the ''model'' argument")
  }

  if(length(data)==0){
    stop("please supply a value to the ''data'' argument")
  }

  if(length(alpha)==0){
    stop("please supply a value to the ''alpha'' argument")
  }
  
  options(warn=1)

  unlink("temp.txt")
  sink(file=NULL,type="message")   
  
  if(log.file!=FALSE)sink(file=log.file,split=TRUE)

  if(is.list(ran.effects)){

    if(is.null(ran.effects$ran.intercepts)){
      ran.effects$ran.intercepts=as.character()
    }

    if(is.null(ran.effects$slopes)){
      ran.effects$slopes=as.character()
    }

    if(is.null(ran.effects$by.vars)){
      ran.effects$by.vars=as.character()
    }

    # Get model coefficients and summary
    if(as.vector(model@call[1])=="glmer()"){
      odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
      data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
      temp.lmer=update(model,.~.,family="gaussian",data=data)
      coefficients=row.names(anova(temp.lmer))
      data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
    } else {
      coefficients=row.names(anova(model))
    }
    coefficients=unique(unlist(strsplit(coefficients,":")))
    coefficients=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",coefficients)
    coefficients=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",coefficients)
    coefficients=c(names(model@flist),coefficients)
    
    if(length(ran.effects$ran.intercepts)>0){
      intercepts=ran.effects$ran.intercepts
      ### Random intercepts ###
      cat("\t===     random intercepts     ===\n")
      # Determine which ranefs to include as a random effects
      for(intercept in intercepts){
        if(!intercept%in%coefficients){
          cat("Warning:",intercept,"not part of model coefficients\n")
          cat("\tskipping\n")
        } else {
          cat("evaluating addition of", paste("by-",intercept,sep=""),"random intercepts to model\n")
          
            wngs=file("temp.txt",open="w+",blocking=TRUE)
            sink(wngs,type="message")
      
          # Fit more complex model
          eval(parse(text=paste("model.updated=update(model,.~.+(1|",intercept,"))",sep="")))
      
        # Should the model with the more complex random-effects structure be kept?
          warn=readLines(wngs)
          unlink("temp.txt")
          sink(file=NULL,type="message")   
      
          if(length(warn)>0){
            if(length(grep("Warning",warn,value=TRUE))>0){
              cat(paste("\t",warn,"\n",sep=""))
              cat("\tnot adding",paste("by-",intercept,sep=""),"random intercepts to model\n")
            }
          } else {
              if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
                cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                cat("\tadding",paste("by-",intercept,sep=""),"random intercepts to model\n")
                model=model.updated
              } else {
                cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                cat("\tnot adding",paste("by-",intercept,sep=""),"random intercepts to model\n")
              }
          }
        }
      }
    }
  
    ### Random slopes ###
    cat("\t===       random slopes       ===\n")
#    if(length(ran.effects$slopes)==0){
#      if(as.vector(model@call[1])=="glmer()"){
#        odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
#        data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
#        temp.lmer=update(model,.~.,family="gaussian",data=data)
#        coefficients=row.names(anova(temp.lmer))
#        data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
#        slopes=unique(unlist(strsplit(row.names(anova(temp.lmer)),":")))
#        slopes=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",slopes)
#      } else {
#        slopes=unique(unlist(strsplit(row.names(anova(model)),":")))
#        slopes=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",slopes)
#      }
#    } else { 
#      slopes=ran.effects$slopes
#    }

    if(length(ran.effects$slopes)>0){
      slopes=ran.effects$slopes
      by.vars=ran.effects$by.vars
  
      for(slope in slopes){
        if(!slope%in%coefficients){
          cat("Warning:",slope,"not part of model coefficients\n")
          cat("\tskipping\n")
        } else {
          for(var in by.vars){
            cat("evaluating addition of",slope,"|",var,"random slopes to model\n")
            if(slope!=var){
                wngs=file("temp.txt",open="w+")
                sink(wngs,type="message")
                if(!var%in%coefficients){
                  cat("\tWarning:",var,"not part of model coefficients\n")
                  cat("\tskipping\n")
                } else {
                  # Fit more complex model
                  eval(parse(text=paste("model.updated=update(model,.~.+(0+",slope,"|",var,"))",sep="")))
            
                  # Should the model with the more complex random-effects structure be kept?
                  warn=readLines(wngs)
                  unlink("temp.txt")
                  sink(file=NULL,type="message")        
            
                  if(length(warn)>0){
                    if(length(grep("Warning",warn,value=TRUE))>0){
                      cat(paste("\t",warn,"\n",sep=""))
                      cat("\tnot adding",slope,"|",var,"random slopes to model\n")
                    }
                  } else {      
                    if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
                      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                      cat("\tadding",slope,"|",var,"random slopes to model\n")
                      model=model.updated
                    } else {
                      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                      cat("\tnot adding",slope,"|",var,"random slopes to model\n")
                    }
                 }
               }
             } # close if(slope!=var) loop
           } # close for(var in by.vars) loop
         }
       } # close for(slope in slopes) loop
     } # close if(length(ran.effects$slopes)>0)
    }else{
    for(ranef in ran.effects){
          wngs=file("temp.txt",open="w+",blocking=TRUE)
          sink(wngs,type="message")

          cat("evaluating addition of",ranef,"to model\n")
          # Fit more complex model
          eval(parse(text=paste("model.updated=update(model,.~.+",ranef,")",sep="")))
    
          # Should the model with the more complex random-effects structure be kept?
          warn=readLines(wngs)
          unlink("temp.txt")
          sink(file=NULL,type="message")        
          
          if(length(warn)>0){
            if(length(grep("Warning",warn,value=TRUE))>0){
              cat(paste("\t",warn,"\n",sep=""))
              cat("\tnot adding",ranef,"to model\n")
            }
          } else {      
            if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
              cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              cat("\tadding",ranef,"to model\n")
              model=model.updated
            } else {
              cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              cat("\tnot adding",ranef,"to model\n")
            }
          }
      }
   }
   if(log.file!=FALSE)sink(file=NULL)

   return(model=model)
}

