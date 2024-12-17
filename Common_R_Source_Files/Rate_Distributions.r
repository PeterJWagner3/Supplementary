library(combinat)
library(gtools)
library(sads)
library(untb)

replace_zeros <- function(v,repl=ZERO)	{
for (i in 1:length(v))	if (v[i]==0)	v[i] <- repl
return(v)
}

# Sampling ####
prob_missing_given_lognormal <- function(p,median_sampling,stdev_mag,ncoll)	{
#r <- median_sampling*stdev_mag^qnorm(p)
return((1-median_sampling*stdev_mag^qnorm(p))^ncoll)
#return((1-(median_sampling*stdev_mag^qnorm(p)))^ncoll)
}

# Diversification give Gamma or Beta ####
#rates <- as.numeric(colnames(lnl_rates)); lnl_rates <- lnl_orig;
accersi_best_gamma_rates <- function(lnl_rates,rates)	{
lnl_rates_all <- colSums(lnl_rates);
ml_rate <- rates[match(max(lnl_rates_all),lnl_rates_all)];
cr <- vector(length=nrow(lnl_rates));
for (nn in 1:nrow(lnl_rates))	cr[nn] <- rates[match(max(lnl_rates[nn,]),lnl_rates[nn,])];
lpr <- log(dgamma(rates,shape=1)/sum(dgamma(rates,shape=1)))
#plot(rates,(lpr+lnl_rates_all))
#plot(rates,dgamma(rates,shape,rate))
best_gamma <- data.frame(alpha=as.numeric(0),beta=as.numeric(0));
shape <- 0.05;
newmaxlgl <- maxlgl <- -MAXNO;
while(newmaxlgl>=maxlgl)	{
	new_shape_lnl <- alflgl <- -MAXNO;
	drate <- -(shape/ml_rate)/10;
#	rate <- 2*shape/ml_rate;
	rate <- shape/ml_rate;
	while (new_shape_lnl>=alflgl & rate>0 & (shape/rate)<5 & (shape/rate)>0.01)	{
		lpr <- log(dgamma(rates[rates>0],shape,rate)/sum(dgamma(rates[rates>0],shape,rate)));
		lpr <- c(-MAXNO,lpr);
		new_lnl_rates <- rbind(lnl_rates,lpr);
		overall_lnl_rates <- colSums(new_lnl_rates);
		overall_l_rates <- exp(overall_lnl_rates);
		new_shape_lnl <- log(sum(overall_l_rates));
		if (alflgl<new_shape_lnl)	{
			alflgl <- max(lpr);
			best_rate <- rate;
			rate <- rate+drate;
			} else if ((rate-drate)==(shape/ml_rate))	{
			drate <- -drate;
			new_shape_lnl <- alflgl;
#			rate <- 2*shape/ml_rate;
			rate <- shape/ml_rate;
			rate <- rate+drate;
			}
		}
	newmaxlgl <- alflgl;
	if (maxlgl < newmaxlgl)	{
		best_gamma$alpha <- shape;
		best_gamma$beta <- best_rate;
		maxlgl <- newmaxlgl;
		shape <- shape+0.05;
		}
	}
plot(rates,dgamma(rates,shape=best_gamma$alpha,rate=best_gamma$beta))
}

prob_gamma_turnover <- function(rate,shape,invscale,pmiss,SS1,two_timers,gap_fillers)	{
# calculates P[data | rate] x P[rate | gamma]
# rate: particular rate of loss between bins
# shape: shape parameter of gamma
# invshape: inverse of shape parameter of gamma
# pmiss: probability of missing a taxon in the second bin
# SS1: vector giving numbers of taxa in first bin for different groups
# two_timers: vector giving number of taxa in first bin seen in second bin
# gap_fillers: number of taxa in bin 1 not seen in other bin, but seen right before/after
continuous=TRUE
groups <- length(SS1)

shift <- 1
denom <- last <- 0
x <- 0.01
while (shift>0.000001)	{
	densf <- dgamma(x,shape,invscale)
	denom <- denom+densf
	if (densf<last)	shift <- densf
	last <- densf
	x <- x+0.01
	}
#denom <- integrate(dgamma,lower=0,upper=100,shape=shape,rate=invscale)$value
lnl_rate <- vector(length=groups)
for (g in 1:groups)	{
	hS2 <- seq(two_timers[g]+gap_fillers[g],SS1[g],by=1)
	lnl_rate[g] <- (log(prob_turnover_rate_given_sampling(rate,pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=TRUE))
		+
		log(dgamma(rate,shape,invscale)/denom))
	}
return(exp(sum(lnl_rate)))
}

prob_gamma_turnover_two_parameters <- function(shape,invscale,pmiss,SS1,two_timers,gap_fillers,tiles)	{
# calculates P[data | rate] x P[rate | gamma]
# rate: particular rate of loss between bins
# shape: shape parameter of gamma
# invshape: inverse of shape parameter of gamma
# pmiss: probability of missing a taxon in the second bin
# SS1: vector giving numbers of taxa in first bin for different groups
# two_timers: vector giving number of taxa in first bin seen in second bin
# gap_fillers: number of taxa in bin 1 not seen in other bin, but seen right before/after
continuous=TRUE
groups <- length(SS1)

shift <- 1
denom <- last <- 0
x <- 0.01
while (shift>0.000001)	{
	densf <- dgamma(x,shape,invscale)
	denom <- denom+densf
	if (densf<last)	shift <- densf
	last <- densf
	x <- x+0.01
	}
#denom <- integrate(dgamma,lower=0,upper=100,shape=shape,rate=invscale)$value
lnl_rate <- vector(length=groups)
for (g in 1:groups)	{
	hS2 <- seq(two_timers[g]+gap_fillers[g],SS1[g],by=1)
	lnl_rate[g] <- (log(prob_turnover_rate_given_sampling(rate,pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=TRUE))
		+
		log(dgamma(rate,shape,invscale)/denom))
	}
return(exp(sum(lnl_rate)))
}

log_prob_turnover_conditional_gamma <- function(p0,pmiss,SS1,two_timers,gap_fillers,tiles)	{
shape <- p0[1]
invscale <- p0[2]
#print(c(shape,invscale))
pcond_rates <- riemann_gamma(shape,invscale,tiles)
grates <- qgamma(seq(1/(tiles+1),(tiles)/(tiles+1),by=1/(tiles+1)),shape,invscale)
groups <- length(SS1)
l_rate <- vector(length=groups)
for (g in 1:groups)	{
#	hS2 <- seq(two_timers[g]+gap_fillers[g],SS1[g],by=1)	# get range of possibly shared taxa
	l_rate[g] <- 0
	for (r in 1:tiles)	{
		pturn <- prob_turnover_rate_given_sampling(grates[r],pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=TRUE)
		if (pturn<=0)	pturn <- MINNO
		p_data <- exp(log(pturn)+log(pcond_rates[r]))
		l_rate[g] <- l_rate[g]+p_data
		}
	if (l_rate[g]<=0)	l_rate[g] <- MINNO
	}

return(sum(log(l_rate)))
}

get_initial_gamma_from_mean_and_variance <- function(average,variance)	{
invscale <- average/variance
shape <- average*invscale
return(c(shape,invscale))
}

get_initial_beta_from_mean_and_variance <- function(average,variance)	{
# mean = scale_a/(scale_a+scale_b)
# variance = (scale_a*scale_b)/([scale_a+scale_b]^2 x [scale_a+scale_b+1])
prec <- 0.001
scale_a <- seq(prec,1-prec,by=prec)
scale_b <- ((scale_a/average)-scale_a)
beta_var <- (scale_a*scale_b)/(((scale_a+scale_b)^2)*(scale_a+scale_b+1))
beta_var_miss <- abs(beta_var-variance)
best_beta <- match(min(beta_var_miss),beta_var_miss)
return(c(prec*best_beta,beta_var[best_beta]))
}

log_prob_turnover_conditional_beta <- function(p0,pmiss,SS1,two_timers,gap_fillers,tiles)	{
shape_a <- p0[1]	# alpha
shape_b <- p0[2]	# beta
#print(c(shape_a,shape_b))
brates <- get_beta_percentiles(prob_vector=seq(1/(tiles+1),(tiles)/(tiles+1),by=1/(tiles+1)),shape_a,shape_b)
#dbeta(brates,shape_a,shape_b)/sum(dbeta(brates,shape_a,shape_b))
pcond_rates <- riemann_beta(shape_a,shape_b,tiles)
groups <- length(SS1)
l_rate <- vector(length=groups)
for (g in 1:groups)	{
	l_rate[g] <- 0
#	debug <- vector(length=tiles)
	for (r in 1:tiles)	{
		pturn <- prob_turnover_rate_given_sampling(rate=brates[r],pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=FALSE)
		if (pturn<=0)	pturn <- MINNO
		p_data <- exp(log(pturn)+log(pcond_rates[r]))
#		debug[r] <- exp(log(pturn)+log(pcond_rates[r]))
		l_rate[g] <- l_rate[g]+p_data
		}
	if (l_rate[g]<=0)	l_rate[g] <- MINNO
	}
return(sum(log(l_rate)))
}

prob_turnover_given_beta <- function(rate,shape_a,shape_b,pmiss,S1,two_timers,gap_fillers)	{
if (rate==0)	{
	pturn <- dbinom(two_timers,S1,pmiss)
	}	else {
	pturn <- dbeta(rate,shape_a,shape_b)*prob_turnover_rate_given_sampling(rate,pmiss,S1=S1,two_timers=two_timers,gap_fillers=gap_fillers,continuous=FALSE)
	}
return(pturn)
}

smooth_beta_turnover_by_groups <- function(SS1,twotimers,gaptfillers,pfind)	{
groups <- length(SS1)
gg <- 1+(10*groups)
trates <- seq(0.1/gg,(gg-1)/gg,by=1/gg)
tr <- length(trates)
ptb <- ptd <- vector(length=tr)
for (g in 1:groups)	{
	if (g==1)	{
		for (r in 1:tr)	ptb[r] <- prob_turnover_rate_given_sampling(trates[r],pmiss=1-pfind,S1=SS1[g],two_timer=twotimers[g],gap_filler=gapfillers[g],continuous=FALSE)
		ptb <- SS1[g]*ptb/sum(ptb)
		}	else	{
		for (r in 1:tr)	ptd[r] <- prob_turnover_rate_given_sampling(trates[r],pmiss=1-pfind,S1=SS1[g],two_timer=twotimers[g],gap_filler=gapfillers[g],continuous=FALSE)
		ptd <- SS1[g]*ptd/sum(ptd)
		ptb <- rbind(ptb,ptd)
		}
	}
smooth <- colSums(ptb)
smoothie <- vector(length=groups)
for (g in 1:groups)	{
	a <- 1+(10*(g-1))
	b <- 10*g
	smoothie[g] <- sum(smooth[a:b])
	}
}

log_likelihood_gamma_turnover_given_different_groups <- function(g0,weighted_distribution,groups,max_rate)	{
rateline <- seq(max_rate/groups,max_rate,by=max_rate/groups)-max_rate/(2*groups)
gamma_dist <- dgamma(rateline,g0[1],g0[2])/sum(dgamma(rateline,g0[1],g0[2]))
return(sum(weighted_distribution*log(gamma_dist)))
}

log_likelihood_beta_turnover_given_different_groups <- function(b0,weighted_distribution,groups)	{
beta_dist <- riemann_beta(shape_a=b0[1],shape_b=b0[2],tiles=groups)
return(sum(weighted_distribution*log(beta_dist)))
}

saturated_turnover_given_different_groups <- function(weighted_distribution)	{
wd <- weighted_distribution[weighted_distribution>0]
xyx <- wd/sum(wd)
return(sum(weighted_distribution*log(xyx)))
}
#p0 <- bg1_result_o$par
#p0 <- c(A,0.6405495)
#A <- 3
#dobby <- 0.6407948*qgamma(xxx,A,A)
#plot(dobby,dgamma(dobby/0.6407948,A,A),xlim=c(0,max(spec_origs)))
log_prob_turnover_conditional_one_param_scaled_gamma <- function(p0,pmiss,SS1,two_timers,gap_fillers,tiles)	{
shape <- invscale <- p0[1]
rescale <- p0[2]
#print(c(shape,invscale))
pcond_rates <- riemann_gamma(shape,invscale,tiles)
grates <- rescale*qgamma(seq(1/(tiles+1),(tiles)/(tiles+1),by=1/(tiles+1)),shape,invscale)
groups <- length(SS1)
l_rate <- vector(length=groups)
for (g in 1:groups)	{
#	hS2 <- seq(two_timers[g]+gap_fillers[g],SS1[g],by=1)	# get range of possibly shared taxa
	l_rate[g] <- 0
	for (r in 1:tiles)	{
		pturn <- prob_turnover_rate_given_sampling(grates[r],pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=TRUE)
		if (pturn<=0)	pturn <- MINNO
		p_data <- exp(log(pturn)+log(pcond_rates[r]))
		l_rate[g] <- l_rate[g]+p_data
		}
	if (l_rate[g]<=0)	l_rate[g] <- MINNO
	}

return(sum(log(l_rate)))
}

#p0 <- c(g0b[1],g0b[2],best_uni_for_families/(g0b[1]/g0b[2]))
#p0 <- bg2_result_o$par
log_prob_turnover_conditional_rescaled_gamma <- function(p0,pmiss,SS1,two_timers,gap_fillers,tiles)	{
shape <- p0[1]
invscale <- p0[2]
rescale <- p0[3]
#print(c(shape,invscale))
pcond_rates <- riemann_gamma(shape,invscale,tiles)
grates <- rescale*qgamma(seq(1/(tiles+1),(tiles)/(tiles+1),by=1/(tiles+1)),shape,invscale)
groups <- length(SS1)
l_rate <- vector(length=groups)
for (g in 1:groups)	{
#	hS2 <- seq(two_timers[g]+gap_fillers[g],SS1[g],by=1)	# get range of possibly shared taxa
	l_rate[g] <- 0
	for (r in 1:tiles)	{
		pturn <- prob_turnover_rate_given_sampling(grates[r],pmiss,S1=SS1[g],two_timers[g],gap_fillers[g],continuous=TRUE)
		if (pturn<=0)	pturn <- MINNO
		p_data <- exp(log(pturn)+log(pcond_rates[r]))
		l_rate[g] <- l_rate[g]+p_data
		}
	if (l_rate[g]<=0)	l_rate[g] <- MINNO
	}

return(sum(log(l_rate)))
}

{}