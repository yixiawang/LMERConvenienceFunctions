pamer.fnc <-
function(model,ndigits=4){ 
	dims<-NULL
	rank.X=qr(model@X)$rank
	anova.table=anova(model)
	anova.table=cbind(anova.table,"Upper Den. DF"=nrow(model@frame)-rank.X)
	p.values.upper=as.numeric()
	p.values.lower=as.numeric()
	for(term in row.names(anova.table)){
		# upper-bound p-values
		p.values.upper=c(p.values.upper,round(1-pf(anova.table[term,"F value"],anova.table[term,"Df"],nrow(model@frame)-rank.X),ndigits))

		# loer-bound p-values
		model.ranef<-ranef(model)
		lower.bound<-0
		for(i in names(model.ranef)){
			eval(parse(text=paste("dims<-dim(model.ranef$",i,")",sep="")))
			lower.bound<-lower.bound+dims[1]*dims[2]
		}
		p.values.lower=c(p.values.lower,1-pf(anova.table[term,"F value"],anova.table[term,"Df"],nrow(model@frame)-rank.X-lower.bound))


	}
	anova.table=round(cbind(anova.table,"Upper p value"=p.values.upper,
		"Lower Den. DF"=nrow(model@frame)-rank.X-lower.bound,"Lower p value"=p.values.lower),
		ndigits)

	return(anova.table)
}

