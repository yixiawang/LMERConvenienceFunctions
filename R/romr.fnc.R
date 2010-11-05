romr.fnc <-
function(model,data){
	data$rstand = as.vector(scale(resid(model))) 
	row.names(data)=1:nrow(data)
	outliers=as.numeric(row.names(data[abs(data$rstand)>2.5,]))
	data0=data
	data=data[-outliers,,drop=TRUE]
	cat("n.removed =",(nrow(data0)-nrow(data)),"\n")
	cat("percent.removed =",(nrow(data0)-nrow(data))/nrow(data0)*100,"\n")
	return(list(data=data,data0=data0,n.removed=nrow(data0)-nrow(data),percent.removed=(nrow(data0)-nrow(data))/nrow(data0)*100))
}
