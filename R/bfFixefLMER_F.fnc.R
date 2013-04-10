bfFixefLMER_F.fnc<-function(model,
	item = FALSE,
	method=c("F","llrt","AIC", "BIC","relLik.AIC","relLik.BIC"),
	threshold=NULL,
	alpha=NULL,
	alphaitem=NULL,
	prune.ranefs=TRUE,
	p.value="upper",
	set.REML.FALSE=TRUE,
	keep.single.factors=FALSE,
	reset.REML.TRUE=TRUE,
	log.file=NULL
){
	###################################
	### perform some initial checks ###
	###################################
	if(length(item)==0){
		stop("please supply a value to argument \"item\".\n")
	}
	if(length(set.REML.FALSE)==0){
		stop("please supply a value to argumnet \"set.REML.FALSE\".\n")
	}
	if(length(reset.REML.TRUE)==0){
		stop("please supply a value to argument \"reset.REML.TRUE\".\n")
	}
	if(!method[1]%in%c("F","llrt","AIC", "BIC","relLik.AIC","relLik.BIC")){
		stop("please supply a proper method name (F llrt AIC BIC, relLik.AIC, or relLik.BIC).\n")}
	if(as.vector(model@call[1])=="glmer()"){
		stop("sorry, glmer models not yet supported")
	}

	##########################
	### some preliminaries ###
	##########################
	# determine threshold
  	if(is.null(threshold)){
		threshold<-c(0.05,0.05,5,5,4,4)[match(method[1],
			c("F","llrt","AIC", "BIC","relLik.AIC","relLik.BIC"))]
	}

	# determine alpha
	if(is.null(alpha)){
		if(method[1]=="F"){
			alpha<-threshold
		}else{
			alpha<-0
		}
	}

	# determine alpha item
	if(is.null(alphaitem)){
		if(method[1]=="F"||method[1]=="llrt"){
			alphaitem<-threshold
		}else{
			alphaitem<-0.05
		}
	}

	# determine log.file name	
	if(is.null(log.file)){
		log.file<-file.path(tempdir(),paste("bfFixefLMER_F_log_",
			gsub(":","-",gsub(" ","_",date())),".txt",sep=""))
	}

	# get data frame from model
	data<-model@frame

	# set REML to FALSE if method != F and set.REML.FALSE = TRUE
		if(method[1]!="F"&set.REML.FALSE){
  		cat("setting REML to FALSE\n")
  		model<-update(model,.~.,REML=FALSE)
	}
  
	# set up things to capture error messages
	options(warn=1)
	temp.dir<-tempdir()
	tempdir()
	unlink(file.path(temp.dir,"temp.txt"))
	sink(file=NULL,type="message")

	# open connection to log file
	if (log.file!=FALSE) 
		sink(file=log.file,split=TRUE)

	############################################################
	### test whether by-item random intercepts are warranted ###
	############################################################
	if(item!=FALSE){
		cat("checking",paste("by-",item,sep=""),"random intercepts\n")
		model.updated<-NULL
		# update model with by-item random intercepts
		eval(parse(text=paste("model.updated<-update(model,.~.+(1|",
			item,"))",sep="")))
		# log-likelihood ratio test
		if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alphaitem){
			cat("  log-likelihood ratio test p-value =",as.vector(anova(model, 
				model.updated)[2,"Pr(>Chisq)"]),"\n")
			cat("  adding",paste("by-",item,sep=""),
				"random intercepts to model\n")
			model<-model.updated
		}else{
			cat("  log-likelihood ratio test p-value =",as.vector(anova(model, 
				model.updated)[2,"Pr(>Chisq)"]),"\n")
			cat("  not adding",paste("by-",item,sep=""),
				"random intercepts to model\n")
		}
	}

	#########################
	### begin backfitting ###
	#########################
	# get model terms
	coefs<-row.names(anova(model))

	# get model summary and format it
	smry<-pamer.fnc(model)
	smry.temp<-pamer.fnc(model)
	temp<-strsplit(coefs, ":")
	names(temp)<-coefs

	# get interaction order
	intr.order<-list()
	for (i in coefs){
		intr.order[[i]]<-length(temp[[i]])
	}
	intr.order<-as.data.frame(unlist(intr.order))
	colnames(intr.order)<-"Order"
	intr.order$Coef<-row.names(intr.order)
	row.names(intr.order)<-1:nrow(intr.order)
	orders<-sort(as.numeric(unique(intr.order$Order)),decreasing=TRUE)
 
	# Do not evaluate single factors if 
	# keep.single.factors = TRUE -- Johannes Ransijn
  	if(keep.single.factors){
		orders<-orders[orders!=1]
	}
	smry.temp$Order<-intr.order$Order

	# set iteration counter to 1
	count<-1
	order<-orders[1]
  
  	# begin backfitting
  	for(order in orders){
		cat("processing model terms of interaction level",order,"\n")
		# which rows of the current summary should be kept
		keepers<-as.character(row.names(smry.temp[smry.temp$Order==order,]))
		smry.temp2<-smry.temp[keepers,]

		# check if any term of the current order is non-significant
		if(smry.temp2[smry.temp2[,paste(p.value,".p.val", 
		sep="")]==max(smry.temp2[,paste(p.value,".p.val", 
		sep="")]),paste(p.value,".p.val",sep="")]<=alpha){
	      cat("  all terms of interaction level",order,"significant\n")
	    }else{
			# if at least one term is not significant at the current 
			# order, backfit at that order
			while(smry.temp2[smry.temp2[,paste(p.value,".p.val", 
	        sep="")]==max(smry.temp2[,paste(p.value, 
	        ".p.val",sep="")]),paste(p.value,".p.val", 
	        sep="")][1]>=alpha){
				cat("  iteration",count,"\n")
				evaluated.term <-row.names(smry.temp2[smry.temp2[, 
					paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, 
					paste(p.value, ".p.val", sep = "")]), ])
				cat("    p-value for term",paste("\"",evaluated.term, 
	 				"\"",sep=""),"=",smry.temp2[evaluated.term, 
	 				paste(p.value,".p.val",sep="")],">",alpha,"\n")
				# changed: check for occurrence of the evaluated term 
				# in each higher order interaction instead of checking for 
				# occurrence in a combination of all terms occurring in higher 
				# order interaction -- Johannes Ransijn
				hoi<-gsub("\\)","",gsub("\\(","",intr.order[intr.order$Order>
					order,"Coef"]))
				tt<-data.frame()
				for(i in 1:length(hoi)){
					tt[i,1]<-length(grep("FALSE",
						unique(unlist(strsplit(evaluated.term, 
						":"))%in%unlist(strsplit(hoi[i],":")))))<1
				}
	          	if (length(grep("TRUE",tt[,1]))!=0){
	          		#changed until here -- Johannes Ransijn
	            	cat("    part of higher-order interaction\n")
	            	cat("    skipping term\n")
	            	# if part of higher-order term, keep it, i.e., 
	            	# remove it from variable keepers
	            	keepers<-row.names(smry.temp2)
	            	keepers<-keepers[-grep(as.character(row.names(
						smry.temp2[smry.temp2[,paste(p.value,".p.val",
						sep="")]==max(smry.temp2[,paste(p.value,".p.val",
						sep = "")]), ])),keepers)]
	            	smry.temp2<-smry.temp2[keepers,]
	            	smry.temp2<-na.omit(smry.temp2)
	          	}else{
	            	cat("    not part of higher-order interaction\n")
	            	m.temp<-NULL
	            	# if not part of higher-order interaction,
	            	# fit simpler model without the term
	            	eval(parse(text=paste("m.temp<-update(model,.~.-", 
	             		row.names(smry.temp2[smry.temp2[,paste(p.value, 
	             		".p.val",sep="")]==max(smry.temp2[,paste(p.value,
	             		".p.val",sep="")]),]),")",sep="")))
	            	if(method[1]=="llrt"){
	            		# if llrt is TRUE, perform log-likelihood ratio test
	            		if(as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"])
						<=threshold){
	               			cat("    log-likelihood ratio test p-value =", 
	                    		as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]), 
	                    		"<=",threshold,"\n")
	                		cat("    skipping term\n")
	                		# if llrt is not significant, keep term in 
							# model, i.e., remove it from the keepers frame
	                		keepers<-row.names(smry.temp2)
	                		cat("length =",length(keepers),"\n")
	                		keepers<-keepers[-which(keepers==as.character(
								row.names(smry.temp2[smry.temp2[,3]==
								min(smry.temp2[,3]),]))[1])]
	                		smry.temp2<-smry.temp2[keepers,]
	                		# a reduction term makes the script shorter 
							# (and easier to look through), compared to 
							# repeating the same commands -- Johannes Ransijn
	                		reduction<-FALSE
	              		}else{
	                		reduction<-TRUE
	                		cat("    log-likelihood ratio test p.value =", 
	                    		as.vector(anova(model,m.temp)[2,"Pr(>Chisq)"]), 
	                    		">",threshold,"\n")
	              		}
					}
	            	
	            	if(method[1]%in%c("AIC","BIC")){
	            		# if criterion is used, compare criteria
						m.temp.ic<-attributes(summary(m.temp))$AICtab[1,
							method[1]]
						model.ic<-attributes(summary(model))$AICtab[1,method[1]]
						ic.diff<-m.temp.ic-model.ic

						# if the difference is equal to or larger than
						# threshold, keep the more complex model
						if(ic.diff>=threshold){
							reduction<-FALSE
						}else{
							reduction<-TRUE
						}

						if(!reduction){
	                		cat(paste("    ",method[1]," simple = ",
								round(m.temp.ic),"; ",method[1]," complex = ",
								round(model.ic),"; decrease = ",round(ic.diff),
								" >= ",threshold,"\n",sep=""))
	                		cat("    skipping term\n")
	                		# if the tested parameter reduces the criterion 
							# value, keep term in model, i.e., remove it 
							# from the keepers frame
	                		keepers<-row.names(smry.temp2)
	                		cat("length =",length(keepers),"\n")
	                		keepers<-keepers[-which(keepers==as.character(
								row.names(smry.temp2[smry.temp2[,3]==min(
								smry.temp2[,3]),]))[1])]
	                		smry.temp2<-smry.temp2[keepers, ]
	                		# a reduction term makes the script shorter (and 
							# easier to look through), compared to repeating 
							# the same commands -- Johannes Ransijn        
	              		}else{
	                		cat(paste("    ",method[1]," simple = ",
								round(m.temp.ic),"; ",method[1]," complex = ",
								round(model.ic),"; decrease = ",round(ic.diff),
								" < ",threshold,"\n",sep=""))
	              		}
	            	}
	            
	            	if(method[1]%in%c("relLik.AIC","relLik.BIC")){
						# get criterion
						ic<-gsub("relLik.(.*)","\\1",method[1])
						m.temp.ic<-eval(parse(text=paste(ic,"(m.temp)",
							sep="")))
						model.ic<-eval(parse(text=paste(ic,"(model)",
							sep="")))
						my.relLik<-relLik(m.temp,model,ic)["relLik"]
						# set reduction to FALSE for now. Will be
						# set to the appropriate value below.
						reduction<-FALSE

	                	cat(paste("    ",ic," simple = ",
							round(m.temp.ic),"; ",ic," complex = ",
							round(model.ic),"; relLik = ",round(my.relLik),
							"\n",sep=""))

						# if AIC value of m.temp is smaller or equal to
						# AIC of model, keep simpler model. Otherwise
						# evaluate whether relLik of complex model is
						# at or above threshold.
						if(m.temp.ic<=model.ic){
							reduction<-TRUE
						}

	            		if(!reduction){
							if(my.relLik>=threshold){
	                			cat(paste("    ","relLik of ",
									round(my.relLik)," >= ",threshold,
									"\n",sep=""))
	                			cat("    skipping term\n")
	                			# if the tested parameter reduces the criterion 
								# value, keep term in model, i.e., 
	                			# remove it from the keepers frame            
	                			keepers<-row.names(smry.temp2)
	                			cat("length =",length(keepers),"\n")
	                			keepers<-keepers[-which(keepers==
									as.character(row.names(smry.temp2[smry.temp2[,
									3]==min(smry.temp2[,3]),]))[1])]
	                			smry.temp2<-smry.temp2[keepers,]
	                			# a reduction term makes the script shorter (and 
								# easier to look through), compared to repeating 
								# the same commands -- Johannes Ransijn        
	              			}else{
	                			cat(paste("    ","relLik of ",
									round(my.relLik)," < ",threshold,
									"\n",sep=""))
								reduction<-TRUE
							}
						}else{
							cat(paste("    simple model more likely than",
									 "complex model\n"))
						}
	              	}

	            	if(method[1]=="F"){
						reduction<-TRUE
	            	}  
	            
	            	if(reduction){
	            		# remove term from model and get new summary
	            		cat("    removing term\n")
	            		model<-m.temp
	            		coefs<-row.names(anova(model))
	            		# if no terms in model, break from current loop
	            		if(length(coefs)==0){
	              			break
	            		}
	              		smry<-pamer.fnc(model)
	              		smry.temp<-pamer.fnc(model)
	              		temp<-strsplit(coefs, ":")
	              		names(temp)<-coefs
	              		# get new interaction orders
	              		intr.order<-list()
	              		for (i in coefs){
	                		intr.order[[i]]<-length(temp[[i]])
	              		}
	              		intr.order<-as.data.frame(unlist(intr.order))
	              		colnames(intr.order)<-"Order"
	              		intr.order$Coef<-row.names(intr.order)
	              		row.names(intr.order)<-1:nrow(intr.order)
	              		smry.temp$Order<-intr.order$Order
	              		keepers<-as.character(row.names(smry.temp[
							smry.temp$Order==order,]))
	              		#new -- added by Johannes Ransijn            
	              		smry.temp2<-smry.temp[keepers,]
	              		smry.temp2<-na.omit(smry.temp2)
	            	}
	          	}
	          	# increase iteration counter by 1
	          	count<-count + 1
	          	# break from loop if there are no more terms at current order
	          	if(nrow(smry.temp2)==0){
	          		break
				}
			}
		}
	}

  	# reset REML to TRUE
  	if (reset.REML.TRUE){
    	cat("resetting REML to TRUE\n")
    	model<-update(model,.~.,REML=TRUE)
  	}
  
	############################
  	### prune random effects ###
	############################
	# i.e., remove those random
  	# effects for which there are no counterparts in
  	# fixed-effect structure
  	if(prune.ranefs){
    	cat("pruning random effects structure ...\n")
    	# get random effects in model
    	split1<-gsub(" ","",model@call)[2]
    	split2<-unlist(strsplit(split1,"\\~"))[2]
    	split2<-gsub("\\)\\+\\(", "\\)_____\\(",split2)
    	split2<-gsub("\\(0\\+","\\(0\\&\\&\\&",split2)
    	split2<-gsub("\\(1\\+","\\(0\\&\\&\\&",split2)
    	split3<-unlist(strsplit(split2,"\\+"))
    	split4<-grep("\\|",split3,value=TRUE)
    	split4<-gsub("\\&\\&\\&","\\+",split4)
    	split5<-unlist(strsplit(split4,"_____"))
    	# create vector of random effects
    	m.ranefs<-vector("character")
    	for(vb in 1:length(split5)){
        	if(length(which(unlist(strsplit(split5[vb],""))=="|"))>0){
                	m.ranefs<-c(m.ranefs,split5[vb])
        	}
    	}
    	# get model terms
    	coefs<-row.names(anova(model))
    	# determine which ranefs to remove
    	m.ranef.variables<-gsub("\\((.*)\\|.*","\\1",m.ranefs)
    	m.ranef.variables<-gsub(".\\+(.*)","\\1",m.ranef.variables)
    	m.ranef.variables<-m.ranef.variables[m.ranef.variables!="1"]
    	m.ranef.variables<-m.ranef.variables[m.ranef.variables!="0"]
    	ranef.to.remove<-vector("character")
    	for(m.ranef.variable in m.ranef.variables){
      		if(!m.ranef.variable%in%coefs){
        		ranef.to.remove<-c(ranef.to.remove,m.ranef.variable)
      		}
    	}
    	# if there are ranefs for which there are no fixefs,
    	# remove them from model
    	if(length(ranef.to.remove)>0){
      	rtr<-vector("character")
      	for(iii in ranef.to.remove){
        	rtr<-c(rtr,grep(iii,m.ranefs,value=TRUE))
        	cat(paste("  ",iii," in random effects structure,
				  but not in fixed effects structure\n",sep=""))
        	cat("    removing",iii,"from model ...\n")
      	}
      	eval(parse(text=paste("model<-update(model,.~.-", 
        	paste(rtr,collapse="-"),",data=data)",sep="")))
    	}else{
      	cat("  nothing to prune\n")
    	}
  	}

	##########################
	### wrap everything up ###
	##########################
  	# close error message capturing connection
  	options(warn=0)
  	sink(file=NULL,type="message")

  	# close log file connection
  	if(log.file!=FALSE){
    	cat("log file is",log.file,"\n")
    	sink(file=NULL)
  	}

  	# return backfitted model
  	return(model=model)
}
