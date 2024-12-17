log10_axes <- function(axe,min_ax,max_ax,numbers,linewd,orient) {
axis(axe,at=seq(min_ax,max_ax,by=(max_ax-min_ax)),tcl=0,labels=NA,lwd.ticks=NA, lwd=linewd, las=orient)
ticks <- length(numbers)
#if ((numbers[2]/numbers[1])>=increment)	{
for (i in 1:ticks)	{
	l <- numbers[i]
	axis(2,at=log10(l),tcl=-.3,labels=numbers[i],lwd=0,lwd.ticks=linewd,las=orient)
	}
for (j in 2:9)	{
	l <- j*numbers[1]
	while (log10(l)<=max_ax)	{
#		l <- numbers[i]+j*numbers[i]
		axis(2,at=log10(l),tcl=-.15,labels=NA,lwd=0,lwd.ticks=linewd)
		l <- 10*l
		}
	}
}

geometric_distribution_min_freq <- function(lambda, min){

g <- geodist <- 1-exp(-lambda)
ct <- 1
while (g>min)	{
	ct <- ct+1
	g <- (1-exp(-lambda))*(exp(-lambda))^(ct-1)
	geodist <- c(geodist,g)
	}
return (geodist)
}

geometric_distribution_N_entities <- function(lambda, N){

geodist <- vector(length=N)
for (ct in 1:N)
	geodist[ct] <- (1-exp(-lambda))*(exp(-lambda))^(ct-1)
geodist <- geodist/sum(geodist)
return (geodist)
}

zipf_distribution <- function(mandel, N)	{
zipfdist <- vector(length=N)
for (i in 1:N)
	zipfdist[i] <- i^-(1*mandel)
zipfdist <- zipfdist/sum(zipfdist)
return (zipfdist)
}

lognormal_distribution <- function(sigma, N)	{
logndist <- vector(length=N)
aa <- vector(length=N)
for (i in 1:N)	aa[i] <- (N+1-i)/(N+1)
for (i in 1:N)
	logndist[i] = sigma^qnorm(aa[i])

logndist <- logndist/sum(logndist)
return (logndist)
}

Pielous_J <- function(distribution, S)	{
J <- H <- 0
for (i in 1:S)	H <- H+(distribution[i]*log(distribution[i]))
J=-1*H/log(S)
return (J)
}

random_crap <- function()	{
set_J <- 0.8
lambda <- 0.05		# standard exponential expectation
gg <- geometric_distribution_N_entities(lambda,N=100)
gg_J <- Pielous_J(gg,100)
while (gg_J>set_J)	{
	lambda <- lambda*1.01
	gg <- geometric_distribution_N_entities(lambda,N=100)
	gg_J <- Pielous_J(gg,100)
	}

mandel <- 0.75		# decay parameter of Zip
zz <- zipf_distribution(mandel,N=100)
zz_J <- Pielous_J(zz,100)
while(zz_J>set_J)	{
	mandel <- mandel*1.01
	zz <- zipf_distribution(mandel,N=100)
	zz_J <- Pielous_J(zz,100)
	}

sigma <- 2.5			# magnitude of change ever standard deviation
ll <- lognormal_distribution(sigma,N=100)
ll_J <- Pielous_J(ll,100)
while(ll_J>set_J)	{
	sigma <- 1.01*sigma
	ll <- lognormal_distribution(sigma,N=100)
	ll_J <- Pielous_J(ll,100)
	}

theta <- 0.7
aaa <- zsm(99,theta,0.01)
hh <- sort(aaa,decreasing=TRUE)
hh_J <- Pielous_J(hh,99)
while (hh_J>set_J)	{
	theta <- theta*1.005
	aaa <- zsm(99,theta,0.01)
	hh <- sort(aaa,decreasing=TRUE)
	hh_J <- Pielous_J(hh,99)
	}

mxy <- ceiling(10*max(ll,zz,gg,hh))/10
mny <- 10^ceiling(log10(min(ll,zz,gg,hh)))
numbers <- c(10^-4,10^-3,10^-2,10^-1)
plot(NA,type='n',axes=FALSE,main="",xlab="Rank Abundance",ylab="Prop. Specimens",xlim=c(0,100),ylim=c(log10(mny),log10(mxy)),cex.lab=1)
log10_axes(axe=2,min_ax=log10(mny),max_ax=log10(mxy),numbers,linewd=1.1,orient=2)
axis(1,at=seq(0,100,by=1),tcl=-0.1,labels=FALSE,lwd=0,lwd.ticks=1.1)
axis(1,at=seq(0,100,by=5),tcl=-0.2,labels=FALSE,lwd=0,lwd.ticks=1.1)
axis(1,at=seq(0,100,by=10),tcl=-0.3,labels=TRUE,lwd=1.1,lwd.ticks=1.1)

lines(xx,log10(gg),col="black",lwd=6)
lines(xx,log10(gg),col="yellow",lwd=4)

lines(xx,log10(hh),col="black",lwd=6)
lines(xx,log10(hh),col="orange",lwd=4)

lines(xx,log10(zz),col="black",lwd=6)
lines(xx,log10(zz),col="purple",lwd=4)

lines(xx,log10(ll),col="black",lwd=6)
lines(xx,log10(ll),col="blue",lwd=4)


}