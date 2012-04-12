mcposthoc.fnc<-function(model,var,two.tailed=TRUE,
    mcmc=FALSE,nsim=10000,ndigits=4,mc.cores=1,...){
	# check if multicore is available
	if(!"multicore" %in% .packages(all.available = TRUE)){
		stop("package multicore not available on this machine\n")
	}

        if(mcmc){
	        # check if languageR is available
	        if(!"languageR" %in% .packages(all.available = TRUE)){
		        stop("package languageR not available on this machine\n")
	        }
	        require(languageR,quietly=TRUE)
        }

	# load library
	require(multicore,quietly=TRUE)

	if(!is.list(var)){
		stop("argument \"var\" must be a list\n")
	}
	
	if(is.null(names(var))){
		stop("\"var\" must be a named list\n")
	}

        data<-model@frame
	
	# create list of posthoc tests to do
	for(k in 1:length(names(var))){
		x<-t(unique(eval(parse(text=paste("data[,var[[k]]]",collapse=",")))))
		colnames(x)<-paste(names(var)[k],1:ncol(x),sep="__")
                if(length(var[[k]])==1)rownames(x)<-var[k]
		for(kk in 1:length(x[1,])){
			if(kk==1){
                                if(k==1){
				        ph.list<-list()
                                }
				tmp<-x[,kk]
				names(tmp)<-rownames(x)
				ph.list[[colnames(x)[kk]]]<-tmp
			}else{
				tmp<-x[,kk]
				names(tmp)<-rownames(x)
				ph.list[[colnames(x)[kk]]]<-tmp
			}

		}
	}

        if(!mcmc){
		my.update<-function(model,ph.list.element){
			# do posthoc for an element of ph.list --> relevel and update model
				for(j in 1:(length(names(ph.list.element)))){
					predictor<-names(ph.list.element)[j]
					data[,predictor]<-relevel(data[,predictor],
						as.character(ph.list.element[j]))
				}
				return(invisible(summary(update(model,.~.,data=data))@coefs))
		}
	
		ph.smrys<-mclapply(ph.list,FUN=function(x)my.update(model=model,ph.list.element=x),
				mc.cores=mc.cores,...)
	
		# create ph.names
		ph.names<-gsub("(.*)__.*","\\1",names(ph.smrys))
		for(ii in 1:length(ph.names)){
			ph.names[ii]<-paste(ph.names[ii],paste(ph.list[[ii]],collapse="_"),sep="_")
		}
		names(ph.smrys)<-ph.names
	
		for(mm in 1:length(ph.names)){
			# store summary data in object temp
			temp<-as.data.frame(ph.smrys[[mm]])
	
			# get rank of X matrix
			rank.X=qr(model@X)$rank
	
			# get upper-bound df
			udf<-nrow(model@frame)-rank.X
	
			# get lower-bound df
			model.ranef<-ranef(model)
			lower.bound<-0
			for(nn in 1:length(names(model.ranef))){
				dims<-dim(model.ranef[[nn]])
				lower.bound<-lower.bound+dims[1]*dims[2]
			}
			ldf<-nrow(model@frame)-rank.X-lower.bound
	     	
			multip<-ifelse(two.tailed,2,1)
			temp[,"udf"]<-udf
			temp[,"ldf"]<-ldf
			temp[,"up.p.unadj."]<-round(multip*(1-pt(abs(temp[,"t value"]),udf)),ndigits)
			temp[,"low.p.unadj."]<-round(multip*(1-pt(abs(temp[,"t value"]),ldf)),ndigits)
	      
			# put temp back in summary object
			ph.smrys[[mm]]<-temp
		}
        }


        if(mcmc){
		my.mc.update<-function(model,ph.list.element,nsim=nsim,
                                    ndigits=ndigits,...){
			# do posthoc for an element of ph.list --> relevel and update model
				for(j in 1:(length(names(ph.list.element)))){
					predictor<-names(ph.list.element)[j]
					data[,predictor]<-relevel(data[,predictor],
						as.character(ph.list.element[j]))
				}

                                updated.model<-update(model,.~.,data=data)
                                model.mcmc<-pvals.fnc(object=updated.model,nsim=nsim,
                                    ndigits=ndigits,withMCMC=FALSE,addPlot=FALSE,...)
				return(invisible(model.mcmc$fixed))
		}
	
		ph.smrys<-mclapply(ph.list,FUN=function(x)my.mc.update(model=model,
                        ph.list.element=x,nsim=nsim,ndigits=ndigits),mc.cores=mc.cores,...)
	
		# create ph.names
		ph.names<-gsub("(.*)__.*","\\1",names(ph.smrys))
		for(ii in 1:length(ph.names)){
			ph.names[ii]<-paste(ph.names[ii],paste(ph.list[[ii]],collapse="_"),sep="_")
		}
		names(ph.smrys)<-ph.names
        }

        res<-list(n=nrow(model@frame),var=var,summaries=ph.smrys)
        class(res)<-"mcposthoc"
	# return results
	return(res)
}

