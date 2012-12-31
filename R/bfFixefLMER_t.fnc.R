bfFixefLMER_t.fnc<-function (model, item = FALSE, alpha = 0.05, llrt = FALSE, 
                            prune.ranefs = TRUE, t.threshold = 2, set.REML.FALSE = TRUE, 
                            reset.REML.TRUE = TRUE, log.file = file.path(tempdir(), 
                            paste("bfFixefLMER_t_log_",gsub(":", "-", gsub(" ", "_", date())),
                            ".txt", sep = ""))) 
{
  # perform some initial checks
  if(length(item)==0){
    stop("please supply a value to argument \"item\".\n")
  }
  if(length(alpha)==0){
    stop("please supply a value to argument \"alpha\".\n")
  }
  if(length(set.REML.FALSE)==0){
    stop("please supply a value to argumnet \"set.REML.FALSE\".\n")
  }
  if(length(reset.REML.TRUE)==0){
    stop("please supply a value to argument \"reset.REML.TRUE\".\n")
  }
  if("factor"%in%attributes(attributes(model@frame)$terms)$dataClasses){
    stop("factor variable in model terms, please use function \"bfFixefLMER_F.fnc\" instead.\n")
  }

  # set REML to FALSE if llrt and set.REML.FALSE are TRUE
  if(llrt&&set.REML.FALSE){
    cat("setting REML to FALSE\n")
    model<-update(model,.~.,REML=FALSE)
  }

  # get data frame from model
  data<-model@frame

  # set 'statistic' to 't-value'
  statistic<-"t-value"

  # create a temporary file to which messages will be diverted
  temp.dir<-tempdir()
  tempdir()
  options(warn=1)
  unlink(file.path(temp.dir,"temp.txt"))
  sink(file=NULL,type="message")

  # create log file
  if(log.file!=FALSE) 
    sink(file=log.file,split=TRUE)

  # determine whether (1|Item) is warranted
  if(item!=FALSE){
    cat("checking",paste("by-",item,sep=""),"random intercepts\n")
    model.updated<-NULL
    # update model with by-item random intercepts        
    eval(parse(text=paste("model.updated=update(model,.~.+(1|",item,"))",sep="")))
    # log-likelihood ratio test
    if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model, 
      model.updated)[2,"Pr(>Chisq)"]),"\n")
      cat("\tadding",paste("by-",item,sep=""),"random intercepts to model\n")
      model<-model.updated
    }else{
      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model, 
      model.updated)[2,"Pr(>Chisq)"]),"\n")
      cat("\tnot adding",paste("by-",item,sep=""),"random intercepts to model\n")
    }
  }

   # for GLMER, get the following coefficients table
   if(as.vector(model@call[1])=="glmer()"){
     statistic<-"z-value"
     odv<-data[,as.character(unlist(as.list(model@call))$formula[2])]
     data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
     temp.lmer<-update(model,.~.,family="gaussian",data=data)
     coefs<-row.names(anova(temp.lmer))
     data[,as.character(unlist(as.list(model@call))$formula[2])]<-odv
   }else{   
     coefs<-row.names(anova(model))
   } 

   # for LMER, get this coefficients table
   if(is(model,"mer")){
     model<-asS4(model)
     smry<-as.data.frame(summary(model)@coefs)
   }else{
     stop("the input model is not a mer object\n")
   }

   # format coefficients table
   smry.temp<-smry[-c(1:nrow(smry)),]
   smry<-as.data.frame(summary(model)@coefs)    
   for(coef in coefs){
      # get interaction order        
      intr.order<-length(unlist(strsplit(coef, ":")))
      orig.coef<-coef
      coef<-gsub(":","\\.\\*:",coef)
      coef<-gsub("(^.*$)","\\1\\.\\*",coef)
      coef<-gsub("\\(","\\\\(",coef)
      coef<-gsub("\\)","\\\\)",coef)
      smry.temp2<-smry[grep(paste("^",coef,sep=""),row.names(smry)),]
      smry.temp3<-smry.temp2
      smry.temp3<-smry.temp3[-c(1:nrow(smry.temp3)),]
    for(j in 1:nrow(smry.temp2)){
     if(length(unlist(strsplit(row.names(smry.temp2)[j],":")))==intr.order){
       smry.temp3<-rbind(smry.temp3,smry.temp2[j,])
     }
    }
    smry.temp2<-smry.temp3
    smry.temp2[,3]=abs(smry.temp2[,3])
    smry.temp2<-smry.temp2[smry.temp2[,3]==max(smry.temp2[,3]),]
    row.names(smry.temp2)<-orig.coef
    smry.temp<-rbind(smry.temp,smry.temp2)
  }

  temp<-strsplit(coefs,":")
  names(temp)<-coefs
  intr.order<-list()
  for(i in coefs){
    intr.order[[i]]<-length(temp[[i]])
  }
  intr.order<-as.data.frame(unlist(intr.order))
  colnames(intr.order)<-"Order"
  intr.order$Coef<-row.names(intr.order)
  row.names(intr.order)<-1:nrow(intr.order)
  orders<-sort(as.numeric(unique(intr.order$Order)),decreasing=TRUE)
  smry.temp$Order<-intr.order$Order
  # set iteration counter to 1    
  count<-1
  # begin backfitting    
  for(order in orders){
      cat("processing model terms of interaction level",order,"\n")
      # which rows of the current summary should be kept      
      keepers<-as.character(row.names(smry.temp[smry.temp$Order==order,]))
      smry.temp2<-smry.temp[keepers,]
      # check if any term of the current order is non-significant      
      if(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),3][1]>=t.threshold){
        cat("\tall terms of interaction level",order,"significant\n")
      }
      # if at least one term is not significant at the current 
      # order, backfit at that order      
      while(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),3][1]<t.threshold){
        cat("\titeration", count, "\n")
        cat("\t\t",statistic,"for term",paste("\"",row.names(smry.temp2[smry.temp2[, 
          3]==min(smry.temp2[,3]),])[1],"\"",sep=""), 
          "=",smry.temp2[smry.temp2[,3]==min(smry.temp2[, 
          3]),3][1],"<",t.threshold,"\n")
        # changed: check for occurence of the exact evaluated terms in 
        # each higher order interaction instead of checking occurence 
        # of the pattern -- Johannes Rasijn
        fac<-unlist(strsplit(gsub("\\)","",gsub("\\(","", 
          as.character(row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[,3]),
          ]))[1])),":"))
        hoi<-gsub("\\)","",gsub("\\(","",intr.order[intr.order$Order>order,"Coef"]))
        tt<-c()
        for(i in 1:length(hoi)){
          tt[i]<-length(grep("FALSE",unique(fac %in% unlist(strsplit(hoi[i],":")))))<1
	}
        if(length(grep("TRUE",tt))!=0){
          # changed until here -- Johannes Rasijn
          cat("\t\tpart of higher-order interaction\n")
          cat("\t\tskipping term\n")
          # if part of higher-order term, keep it, i.e., 
          # remove it from variable keepers          
          keepers<-row.names(smry.temp2)
          keepers<-keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[, 
            3]==min(smry.temp2[,3]),]))[1],keepers)]
          smry.temp2<-smry.temp2[keepers,]
        }else{
          cat("\t\tnot part of higher-order interaction\n")
          m.temp <- NULL
          # if not part of higher-order interaction,
          # fit simpler model without the term                
          eval(parse(text=paste("m.temp=update(model,.~.-", 
            row.names(smry.temp2[smry.temp2[,3]==min(smry.temp2[, 
            3]),])[1],")",sep="")))
          if(llrt){
          # if llrt is TRUE, perform log-likelihood ratio test
            if(as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"])<=alpha){
            cat("\t\tlog-likelihood ratio test p-value =", 
              as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]), 
              "<=",alpha,"\n")
            cat("\t\tskipping term\n")
            # if llrt is not significant, keep term in model, i.e., 
            # remove it from the keepers frame            
            keepers<-row.names(smry.temp2)
            cat("length =",length(keepers),"\n")
            keepers<-keepers[-which(keepers==as.character(row.names(smry.temp2[smry.temp2[, 
              3]==min(smry.temp2[,3]),]))[1])]
            smry.temp2<-smry.temp2[keepers, ]
            # a reduction term makes the script shorter (and easier to look 
            # through), compared to repeating the same commands -- Johannes Rasijn        
            reduction<-FALSE
          }else{
            reduction<-TRUE
            cat("\t\tlog-likelihood ratio test p.value =", 
              as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]), 
              ">",alpha,"\n")
          }
        }else{
          reduction<-TRUE
        }
        if(reduction){
        # remove term from model and get new summary                    
          cat("\t\tremoving term\n")
          model<-m.temp
          if (as.vector(model@call[1])=="glmer()"){
            odv<-data[,as.character(unlist(as.list(model@call))$formula[2])]
            data[,as.character(unlist(as.list(model@call))$formula[2])]<-rnorm(nrow(data),0,1)
            temp.lmer<-update(model,.~.,family="gaussian",data=data)
            coefs<-row.names(anova(temp.lmer))
            data[,as.character(unlist(as.list(model@call))$formula[2])]<-odv
          }else{
            coefs<-row.names(anova(model))
            # if no terms in model, break from current loop                    
            if(length(coefs)==0){
              break
            }
          }
          smry2<-as.data.frame(summary(asS4(model))@coefs)
          smry.temp2<-smry2[-c(1:nrow(smry2)),]
          # get new interaction orders                  
          for(coef in coefs){
            intr.order<-length(unlist(strsplit(coef,":")))
            orig.coef<-coef
            coef<-gsub(":","\\.\\*:",coef)
            coef<-gsub("(^.*$)","\\1\\.\\*",coef)
            coef<-gsub("\\(","\\\\(",coef)
            coef<-gsub("\\)","\\\\)",coef)
            smry.temp3<-smry2[grep(paste("^",coef,sep=""),row.names(smry2)),]
            smry.temp4<-smry.temp3
            smry.temp4<-smry.temp4[-c(1:nrow(smry.temp4)),]
            for(j in 1:nrow(smry.temp3)){
              if (length(unlist(strsplit(row.names(smry.temp3)[j],":")))==intr.order){
                smry.temp4<-rbind(smry.temp4,smry.temp3[j,])
              }
            }
            smry.temp3<-smry.temp4
            smry.temp3[,3]<-abs(smry.temp3[,3])
            smry.temp3<-smry.temp3[smry.temp3[,3]==max(smry.temp3[,3]),]
            row.names(smry.temp3)<-orig.coef
            smry.temp2<-rbind(smry.temp2,smry.temp3)
          }
          temp<-strsplit(coefs,":")
          names(temp)<-coefs
          intr.order<-list()
          for(i in coefs){
            intr.order[[i]]<-length(temp[[i]])
          }
          intr.order<-as.data.frame(unlist(intr.order))
          colnames(intr.order)<-"Order"
          intr.order$Coef<-row.names(intr.order)
          row.names(intr.order)<-1:nrow(intr.order)
          smry.temp2$Order<-intr.order$Order
          keepers<-as.character(row.names(smry.temp2[smry.temp2$Order==order,]))
          # smry.temp needs to be changed when the model is updated -- Johannes Ransijn
          smry.temp=smry.temp2
          # -- Johannes Ransijn
          smry.temp2<-smry.temp2[keepers,]
        }
      }
      # increase iteration counter by 1      
      count<-count+1
      # break from loop if there are no more terms at current order      
      if(nrow(smry.temp2)==0){
        break
      }
    }
  }
  
  # reset REML to TRUE
  if(llrt&&reset.REML.TRUE) {
    cat("resetting REML to TRUE\n")
    model<-update(model,.~.,REML=TRUE)
  }
  # prune random effects, i.e., remove those random
  # effects for which there are no counterparts in
  # fixed-effect structure  
  if(prune.ranefs){
    cat("pruning random effects structure ...\n")
    # get random effects in model    
    split1<-gsub(" ","",model@call)[2]
    split2<-unlist(strsplit(split1,"\\~"))[2]
    split2<-gsub("\\)\\+\\(","\\)_____\\(",split2)
    split2<-gsub("\\(0\\+","\\(0\\&\\&\\&",split2)
    split2<-gsub("\\(1\\+","\\(0\\&\\&\\&",split2)
    split3<-unlist(strsplit(split2,"\\+"))
    split4<-grep("\\|",split3,value=TRUE)
    split4<-gsub("\\&\\&\\&","\\+",split4)
    split5<-unlist(strsplit(split4,"_____"))
    m.ranefs<-split5
    coefs<-row.names(anova(model))
    # determine which ranefs to remove    
    m.ranef.variables<-gsub("\\((.*)\\|.*","\\1",m.ranefs)
    m.ranef.variables<-gsub(".\\+(.*)","\\1",m.ranef.variables)
    m.ranef.variables<-m.ranef.variables[m.ranef.variables!="1"]
    m.ranef.variables<-m.ranef.variables[m.ranef.variables!="0"]
    # create vector of random effects    
    ranef.to.remove<-vector("character")
    for(m.ranef.variable in m.ranef.variables){
      if (!m.ranef.variable%in%coefs){
        ranef.to.remove<-c(ranef.to.remove,m.ranef.variable)
      }
    }
    # if there are ranefs for which there are no fixefs,
    # remove them from model    
    if(length(ranef.to.remove)>0){
      rtr<-vector("character")
      for(iii in ranef.to.remove){
        rtr<-c(rtr,grep(iii,m.ranefs,value=TRUE))
        cat(paste("\t",iii," in random effects structure, \n\tbut not in fixed effects structure\n",sep=""))
        cat("\t\tremoving",iii,"from model ...\n")
       }
       eval(parse(text=paste("model<-update(model,.~.-",paste(rtr,collapse="-"),",data=data)",sep="")))
    }else{
      cat("\tnothing to prune\n")
    }
  }

  # close error message capturing connection  
  options(warn=0)
  sink(file=NULL,type="message")
  if(log.file!=FALSE){
    cat("log file is",log.file,"\n")
    sink(file=NULL)
  }

  # close log file connection  
  return(model=model)
}
