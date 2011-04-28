mcp.fnc <-
function(model,trim=2.5){
	data<-model@frame
	par(mfrow=c(2,2))
	data$rstand = as.vector(scale(resid(model))) 
	plot(density(data$rstand),main="") 
	qqnorm(data$rstand,pch = ".",main ="") 
	qqline(data$rstand,col="red")
	plot(data$rstand~fitted(model),pch=".")
	abline(h = c(-trim, trim),col="red")
	dffits = abs(resid(model, "dffits")) 
	plot(dffits, type="h")
	par(mfrow=c(1,1))
}

