mcposthoc.fnc<-function(model,data,var,two.tailed=TRUE,num.comp=NULL,ndigits=4,verbose=FALSE){
	# check if multicore is available
	if(!"multicore" %in% .packages(all.available = TRUE)){
		stop("package multicore not available on this machine\n")
	}

	# load library
	require(multicore,quietly=TRUE)

	if(!is.list(var)){
		stop("argument ''var'' must be a list\n")
	}
	
	if(is.null(names(var))){
		stop("''var'' must be a named list\n")
	}
	
	orig.var<-var
	to.collect<-vector("character")	
	new.var<-NULL
	ph.smrys<-NULL

	for(k in 1:length(names(orig.var))){
		var<-orig.var[[k]]
		# create new variable
		eval(parse(text=paste("new.var<-as.character(unique(paste(",gsub(" ","",paste("data[,'",var,"']",collapse=",")),",sep='___')))",sep="")))
	
		# do posthoc -- for each level of new.var, update model
		for(i in 1:length(new.var)){
			temp<-unlist(strsplit(new.var[i],"___"))
			for(j in 1:(length(temp))){
				eval(parse(text=paste("data[,'",var[j],"']<-relevel(data[,'",var[j],"'],'",temp[j],"')",sep="")))
			}
			nm.tmp<-paste(names(orig.var)[k],"_",gsub("___","_",new.var[i]),sep="")
			if(verbose){
				cat("sending-out job",nm.tmp,"\n")
			}
			eval(parse(text=paste(nm.tmp,"<-mcparallel(summary(update(model,.~.,data=data))@coefs,name=nm.tmp)",sep="")))
			to.collect<-c(to.collect,nm.tmp)
		}
	}
	
	# collect jobs
	if(verbose)cat("waiting to collect jobs ...\n")
	eval(parse(text=paste("ph.smrys<-collect(jobs=list(",paste(to.collect,collapse=","),"),wait=TRUE)",sep="")))

	ph.names<-names(ph.smrys)
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
      
		if(!is.null(num.comp)){
			temp[,"up.p.adj."]<-round(num.comp*(multip*(1-pt(abs(temp[,"t value"]),udf))),ndigits)
			temp[,"low.p.adj."]<-round(num.comp*(multip*(1-pt(abs(temp[,"t value"]),ldf))),ndigits)
		}

		# put temp back in summary object
		ph.smrys[[mm]]<-temp
	}

	# return results
	return(list(n=nrow(model@frame),var=orig.var,summaries=ph.smrys))
}

