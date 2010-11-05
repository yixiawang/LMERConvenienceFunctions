pamer.fnc <-
function(model,ndigits=4){ 
	rank.X=qr(model@X)$rank
	anova.table=anova(model)
	anova.table=cbind(anova.table,"Den. DF"=nrow(model@frame)-rank.X)
	p.values=as.numeric()
	for(term in row.names(anova.table)){
		p.values=c(p.values,round(1-pf(anova.table[term,"F value"],anova.table[term,"Df"],nrow(model@frame)-rank.X),ndigits))
	}
	anova.table=cbind(anova.table,"Upper-bound p-value"=p.values)

	return(anova.table)
}

