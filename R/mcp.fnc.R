mcp.fnc <-
function(model,data){
	par(mfrow=c(2,2))
	data$rstand = as.vector(scale(resid(model))) 
	plot(density(data$rstand),main="") 
	qqnorm(data$rstand,pch = ".",main ="") 
	qqline(data$rstand,col="red")
	plot(data$rstand~fitted(model),pch=".")
	abline(h = c(-2.5, 2.5),col="red")
	dffits = abs(resid(model, "dffits")) 
	plot(dffits, type="h")
	par(mfrow=c(1,1))
}

