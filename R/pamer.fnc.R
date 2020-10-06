pamer.fnc<-function (model, ndigits = 4) 
{
  if (length(rownames(anova(model))) == 0) {
    cat("nothing to evaluate: model has only an intercept.\n\n")
    cat("printing model fixed effects:\n")
    fixef(model)
  }
  else {
    dims <- NULL
    rank.X = qr(model@pp$X)$rank
    anova.table = anova(model)
    anova.table = cbind(anova.table, upper.den.df = nrow(model@frame) - rank.X)
    p.values.upper = as.numeric()
    p.values.lower = as.numeric()
    for (term in row.names(anova.table)) {
      p.values.upper = c(p.values.upper, round(1 - pf(anova.table[term,"F value"], anova.table[term, "npar"], nrow(model@frame) - rank.X), ndigits))
      model.ranef <- ranef(model)
      lower.bound <- 0
      for (i in 1:length(names(model.ranef))) {
        dims <- dim(model.ranef[[i]])
        lower.bound <- lower.bound + dims[1] * dims[2]
      }
      p.values.lower = c(p.values.lower, 1 - pf(anova.table[term, "F value"], anova.table[term, "npar"], nrow(model@frame) - rank.X - lower.bound))
    }
    aov.table <- as.data.frame(anova(model))
	if(!as.vector(model@call[1]) == "glmer()"){
    	dv <- gsub(" ", "", gsub("(.*)~.*", "\\1", as.character(model@call)[2]))
    	ss.tot <- sum((model@frame[, dv] - mean(model@frame[, dv]))^2)
    	expl.dev <- vector("numeric")
    	for (i in rownames(aov.table)) {
      	expl.dev <- c(expl.dev, aov.table[i, 2]/ss.tot)
    	}
    	names(expl.dev) <- rownames(aov.table)
	}else{
		expl.dev<-rep(NA,nrow(aov.table))
	}
    anova.table = round(cbind(anova.table, upper.p.val = p.values.upper, 
		lower.den.df = nrow(model@frame) - rank.X - lower.bound, 
		lower.p.val = p.values.lower, `expl.dev.(%)` = expl.dev * 100), ndigits)
    return(anova.table)
  }
}
