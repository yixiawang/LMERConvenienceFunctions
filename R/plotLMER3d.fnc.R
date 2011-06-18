plotLMER3d.fnc<-function(model,
			pred,
			intr,
			plot.type="contour", # or "persp" or "persp3d"
			xlim = range(x, na.rm = TRUE),
			ylim = range(y, na.rm = TRUE),
			zlim=range(z, na.rm = TRUE), 
			xlab="", 
			ylab="", 
			zlab="", 
			cex=1,
			fun=NA, 
			n=30,
			color="topo",
			lit=TRUE,
			theta=0,
			phi=0,
			contourstepsize=0.2,
			play3d=FALSE,   # or TRUE or a list with first element axis, e.g., c(0,0,1), second element rpm, 
					# e.g., 4, and third element duration, e.g., 20.			
			ref.surf=FALSE,
			xaxt="s",
			yaxt="s"
){

	temp.dir<-tempdir()
	model.name<-as.character(model@call)
	model.name<-gsub(" ","",model.name)
	model.name<-paste(model.name[1],"___",model.name[2],"___",model.name[3],"___",pred,"_",intr,sep="")
	model.name<-gsub("\\+","__",model.name)
	model.name<-gsub("\\:","_",model.name)
	model.name<-paste(model.name,".rda",sep="")

	if(!model.name%in%list.files(path=temp.dir,pattern="lmer___.*\\.rda$")){
		# create LMER plot and store values
		list1<-plotLMER.fnc(model=model,fun=fun,pred=pred,intr=list(intr,
			quantile(model@frame[,intr],seq(0,1,1/n)),"end",list(seq(1,length(seq(0,1,1/n)),1),
			seq(1,length(seq(0,1,1/n)),1))),n=n,withList=TRUE,verbose=TRUE)
	
		# get info from stored plotting list
		length.intr<-eval(parse(text=paste("length(list1$",pred,")",sep="")))
		x<-eval(parse(text=paste("list1$",pred,"[[1]]$X",sep="")))
		y=quantile(model@frame[,intr],seq(0,1,1/n))
	
		# create plotting matrix
		for(i in 1:length.intr){
			if(i==1){
				z<-eval(parse(text=paste("list1$",pred,"[[",i,"]]$Y",sep="")))
			}else{
				z<-cbind(z,eval(parse(text=paste("list1$",pred,"[[",i,"]]$Y",sep=""))))
			}	
		}
	
		# add row and column names to matrix
		rownames(z)<-x
		colnames(z)<-y
	
		# remove identical columns
		rem<-vector("numeric")
		for(i in 2:ncol(z)){
			if(unique(z[,i-1]==z[,i])[1]){
				rem<-c(rem,i)
			}
		}
		z<-z[,-rem]
		
		save(z,file=file.path(temp.dir,model.name))
	}else{
		load(file.path(temp.dir,model.name))
	}

	if(is.null(zlim)){
		zlim=range(z)
	}

#	trying.function<-function(z,col,zlim){
#		dev.new()
#		image(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),
#			z=z,col=col,zlim=zlim)
#	}
#

	if(plot.type=="contour"){
		contourlevels = seq(zlim[1], zlim[2], by=contourstepsize)
			
		# Determine color.
        	if(color=="heat"){
            		pal=heat.colors(50)
            		con.col=3
        	}else if(color=="topo"){
            		pal=topo.colors(50)
            		con.col=2
        	}else if(color=="cm"){
            		pal=cm.colors(50)
            		con.col=1
        	}else if(color=="terrain"){
            		pal=terrain.colors(50)
            		con.col=2
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		pal=gray(seq(0.1,0.9,length=50))
            		con.col=1
        	}else{
			stop("color scheme not recognised")
		}

	
#		trying<-try(trying.function(z=z,col=pal,zlim=zlim),silent=FALSE)
#		dev.off()
#
#		if(is.null(trying)){
#			image(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,col=pal,zlim=zlim,xlab=xlab,
#				ylab=ylab,main=zlab,cex.main=cex,cex.lab=cex,cex.axis=cex)
#			contour(z=z,zlim=zlim,add=TRUE,levels=round(contourlevels,2))
#			box()
#		}else{
#			cat("plotting, but will not use supplied x- and y-values ...\n")
			image(z=z,col=pal,zlim=zlim,main=zlab,cex.main=cex,cex.lab=cex,cex.axis=cex,
				xlab=xlab,ylab=ylab,axes=FALSE)
				#xlab=paste(xlab,"-- Random Units",sep=" "),ylab=paste(ylab,
				#"-- Random Units",sep=" "))
			if(xaxt=="s")axis(1,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(rownames(z))),2))
			if(yaxt=="s")axis(2,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(colnames(z))),2))
			contour(z=z,zlim=zlim,add=TRUE,levels=round(contourlevels,2),axes=FALSE)
			box()
		#}
	}else if (plot.type=="persp"){
		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(z)
		ncz<-ncol(z)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(50))
        	}else if(color=="topo"){
			#jet.colors <- colorRampPalette( c("purple","blue", "green","yellow","red","white") ) 
			jet.colors <- colorRampPalette(topo.colors(50)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(50))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(heat.colors(50))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
        	}else{
			stop("color scheme not recognised")
		}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)


		#trying<-try(persp(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,ticktype="detailed",col=color[facetcol],
			#phi=phi,theta=theta,xlab=xlab,ylab=ylab,zlab=zlab,zlim=zlim))

		#if(is.null(trying)){
			#persp(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,ticktype="detailed",col=color[facetcol],
				#phi=phi,theta=theta,xlab=xlab,ylab=ylab,zlab=zlab,zlim=zlim)
		#}else{
			#cat("plotting, but will not use supplied x- and y-values ...\n")
			persp(z=z,ticktype="simple",col=color[facetcol],phi=phi,theta=theta,
				zlab=zlab,zlim=zlim,xlab=xlab,ylab=ylab,axes=TRUE)
				#xlab=paste(xlab,"-- Random Units",sep=" "),
				#ylab=paste(ylab,"-- Random Units",sep=" "))
		#}
	}else{
		if(!"rgl"%in%.packages(all.available=TRUE)){
			stop("Package ''rgl'' not available.\n Please set ''plot.type'' to ''contour'' or ''persp''.\n")
		}	
		require(rgl,quietly=TRUE) 

		# the color portion of this code is adapted from the persp() help page
		#par(bg="white")
		nrz<-nrow(z)
		ncz<-ncol(z)

		# Create a function interpolating colors in the range of specified colors
        	if(color=="heat"){
            		jet.colors<-colorRampPalette(heat.colors(100))
        	}else if(color=="topo"){
			jet.colors <- colorRampPalette(topo.colors(100)) 
        	}else if(color=="cm"){
            		jet.colors<-colorRampPalette(cm.colors(100))
        	}else if(color=="terrain"){
            		jet.colors<-colorRampPalette(heat.colors(100))
        	}else if(color=="gray"||color=="bw"||color=="grey"){
            		jet.colors<-colorRampPalette(gray(seq(0.1,0.9,length=7)))
        	}else{
			stop("color scheme not recognised")
		}

		# Generate the desired number of colors from this palette
		nbcol<-100
		color<-jet.colors(nbcol)

		# Compute the z-value at the facet centres
		zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]

		# Recode facet z-values into color indices
		facetcol<-cut(zfacet,nbcol)
		facetcol=color[facetcol]

		# this portion is from the persp3d() help page
		nx=length(rownames(z))
		ny=length(colnames(z))
		options(warn=-1)
		col <- rbind(0, cbind(matrix(facetcol, nx-1, ny-1), 0))
		options(warn=0)

		# create persp3d plot
		#trying<-try(persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,col=col,zlim=zlim,
			#xlab=xlab,ylab=ylab,zlab=zlab,smooth=FALSE,lit=lit))

		#if(is.null(trying)){
			#persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,col=col,zlim=zlim,
				#xlab=xlab,ylab=ylab,zlab=zlab,smooth=FALSE,lit=lit)
		#}else{
			#cat("plotting, but will not use supplied x- and y-values ...\n")
			persp3d(z=z,col=col,zlim=zlim,zlab=zlab,smooth=FALSE,lit=lit,xlab=xlab,ylab=ylab)
				#xlab=paste(xlab,"-- Random Units",sep=" "),ylab=paste(ylab,
				#"-- Random Units",sep=" "))
			#axis(1,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(rownames(z))),2))
			#axis(2,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(colnames(z))),2))


		#}

		if(ref.surf){
			persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),matrix(0,ncol=ncol(z),nrow=nrow(z),
				byrow=TRUE),col=grey(0.9),zlim=zlim,add=TRUE,lit=lit)
		}

		if(play3d || is.list(play3d)){
			if(is.list(play3d)){
				play3d(spin3d(axis=play3d[[1]], rpm=play3d[[2]]), duration=play3d[[3]])
			}else{
				play3d(spin3d(axis=c(0,0,1),rpm=4),duration=20)
			}
		}


	}
}

