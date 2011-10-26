plotLMER3d.fnc<-function(model=NULL,
			pred,
			intr,
			plot.type="contour", # or "persp" or "persp3d"
			xlim = range(x, na.rm = TRUE),
			ylim = range(y, na.rm = TRUE),
			zlim=range(z, na.rm = TRUE), 
			xlab=NULL, 
			ylab=NULL, 
			zlab=NULL, 
			main=NULL, 
			cex=1,
			fun=NA, 
			n=30,
			color="topo",
			alpha=1,
			alpha.rs=0.65,
			alpha.u=1,
			lit=TRUE,
			theta=0,
			phi=0,
			contourstepsize=0.2,
			play3d=FALSE,   # or TRUE or a list with first element axis, e.g., c(0,0,1), second element rpm, 
					# e.g., 4, and third element duration, e.g., 20.			
			ref.surf=FALSE,
			underneath=FALSE,
			rug=FALSE,
			rug.u=FALSE,
			plot.dat="default",
			path="default",
			ret=FALSE
){
	if(is.null(model) && plot.dat=="default"){
		stop("either provide a value to argument ''model'' (an object of class ''mer'') \n 
			or provide a file name containing plotting information \n")
	}

	# set labels if NULL
	if(is.null(xlab)){
		xlab=pred
	}

	if(is.null(ylab)){
		ylab=intr
	}


     	options(warn=-1)
	if(is.null(zlab)){
		if(try(is.null(model),silent=TRUE)){
			zlab="Response"
		}else{
			zlab<-as.character(model@call)[2]
			zlab<-gsub(" ","",unlist(strsplit(zlab,"~"))[1])
		}
	}
     	options(warn=0)

	if(is.null(main)){
		if(plot.type=="contour"){
			main=zlab
		}else{
			main=""
		}
	}

	# create file name for saving in temp dir
	if(plot.dat!=FALSE){
		if(path=="default"){
			temp.dir<-tempdir()
		}else{
			temp.dir=path
		}
		if(plot.dat=="default"){
			model.name<-as.character(model@call)
			model.name<-gsub(" ","",model.name)
			model.name<-paste(model.name[1],"___",model.name[2],"___",model.name[3],"___",pred,"_",intr,sep="")
			model.name<-gsub("\\+","__",model.name)
			model.name<-gsub("\\:","_",model.name)
			model.name<-gsub("\\*","_",model.name)
			model.name<-gsub("\\^","_",model.name)
			model.name<-gsub("\\|","_",model.name)
			model.name<-gsub("\\~","_",model.name)
			model.name<-gsub("\\(","WWW",model.name)
			model.name<-gsub("\\)","WWW",model.name)
			model.name<-paste(model.name,".rda",sep="")
		}else{
			model.name=paste("lmer___",plot.dat,".rda",sep="")
		}
	}


	# get previously generated plotting info if it exists
	if(plot.dat!=FALSE){
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
			z<-as.matrix(z)
		
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
			if(length(rem)!=0){
				z<-z[,-rem]
			}
			
			save(z,file=file.path(temp.dir,model.name))
		}else{
			load(file.path(temp.dir,model.name))
		}
	}

	if(is.null(zlim)){
		zlim=range(z,na.rm=TRUE)
	}


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

	
		err<-try(image(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,col=pal,zlim=zlim,main=main,
			cex.main=cex,cex.lab=cex,cex.axis=cex,xlab=xlab,ylab=ylab,axes=TRUE),silent=FALSE)

		if(length(err)>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			image(z=z,col=pal,zlim=zlim,main=main,cex.main=cex,cex.lab=cex,cex.axis=cex,
				axes=TRUE,xlab=paste(xlab,"-- Random Units",sep=" "),
				ylab=paste(ylab,"-- Random Units",sep=" "))
			#if(xaxt=="s")axis(1,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(rownames(z))),2))
			#if(yaxt=="s")axis(2,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(colnames(z))),2))
		}
		rm(err)

		err<-try(contour(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,zlim=zlim,
			add=TRUE,levels=round(contourlevels,2),axes=FALSE),silent=TRUE)
		if(length(err)>0){
			contour(z=z,zlim=zlim,add=TRUE,levels=round(contourlevels,2),axes=FALSE)
		}
		rm(err)

		if(rug){
			xy<-expand.grid(as.numeric(rownames(z)),as.numeric(colnames(z)))
			points(xy[,1],xy[,2],pch=19,cex=0.05)
		}

		box()
		if(ret){
			return(list(z=z,col=pal))
		}
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
            		jet.colors<-colorRampPalette(terrain.colors(50))
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


		err<-try(persp(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,
			ticktype="detailed",col=color[facetcol],phi=phi,theta=theta,
			zlab=zlab,zlim=zlim,xlab=xlab,ylab=ylab,main=main,axes=TRUE)->res,silent=FALSE)

		
		if(length(err)>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			persp(z=z,ticktype="detailed",col=color[facetcol],phi=phi,theta=theta,
				zlab=zlab,zlim=zlim,xlab=paste(xlab,"-- Random Units",sep=" "),
				ylab=paste(ylab,"-- Random Units",sep=" "),main=main,axes=TRUE)->res
		}
		rm(err)


		if(rug){
			xy<-expand.grid(as.numeric(rownames(z)),as.numeric(colnames(z)))
			temp<-vector("numeric")
			for(i in 1:nrow(xy)){
				temp<-c(temp,z[as.character(xy$Var1[i]),as.character(xy$Var2[i])])
			}
			points(trans3d(xy[,1],xy[,2],temp,pmat=res),pch=19,cex=0.5)
		}

		if(ret){
			return(list(z=z,col=color[facetcol]))
		}
	}else{
		if(!"rgl"%in%.packages(all.available=TRUE)){
			stop("Package ''rgl'' not available.\n Please set ''plot.type'' to ''contour'' or ''persp''.\n")
		}	
		require(rgl,quietly=TRUE) 

		dev.new()

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
            		jet.colors<-colorRampPalette(terrain.colors(100))
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
		col <- rbind(0, cbind(matrix(facetcol, nx-1, ny-1), 0))

		# create persp3d plot
		err<-try(persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=z,col=col,zlim=zlim,
				zlab=zlab,main=main,alpha=alpha,smooth=FALSE,lit=lit,xlab=xlab,ylab=ylab),silent=FALSE)

		if(length(err)>0){
			cat("\tplotting anyways, but will not use supplied x- and y-values ...\n")
			persp3d(z=z,col=col,zlim=zlim,zlab=zlab,main=main,alpha=alpha,smooth=FALSE,
				lit=lit,xlab=paste(xlab,"-- Random Units",sep=" "),ylab=paste(ylab,
				"-- Random Units",sep=" "))
			#axis(1,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(rownames(z))),2))
			#axis(2,at=quantile(seq(0,1,.1)),labels=round(quantile(as.numeric(colnames(z))),2))
		}
		rm(err)


		if(ref.surf){
			persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),matrix(mean(z),ncol=ncol(z),
				nrow=nrow(z),byrow=TRUE),col="grey",alpha=alpha.rs,zlim=zlim,add=TRUE,lit=lit)
		}

		if(underneath){
			persp3d(x=as.numeric(rownames(z)),y=as.numeric(colnames(z)),z=matrix(min(z),nrow(z),ncol(z)),
				col=col,alpha=alpha.u,add=TRUE,smooth=FALSE,lit=lit,zlim=zlim)	
			if(rug.u){
				for(i in 1:nrow(z)){
					plot3d(x=as.numeric(rownames(z))[i],y=as.numeric(colnames(z)),
						z=matrix(min(z),nrow(z),ncol(z))[i,],add=TRUE,col="black",size=3)
				}
			}
		}

		if(rug){
			for(i in 1:nrow(z)){
				plot3d(x=as.numeric(rownames(z))[i],y=as.numeric(colnames(z)),z=z[i,],
					add=TRUE,col="black",size=3)
			}
		}


		if(play3d || is.list(play3d)){
			if(is.list(play3d)){
				play3d(spin3d(axis=play3d[[1]], rpm=play3d[[2]]), duration=play3d[[3]])
			}else{
				play3d(spin3d(axis=c(0,0,1),rpm=4),duration=20)
			}
		}

		if(ret){
			return(list(z=z,col=col))
		}

		dev.off()
	}
}
