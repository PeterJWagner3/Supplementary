# some basic stats & probability routines that I wrote for myself because I didn't like the way that R did them...
# Rossetta Stone:
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal
# scribere: write

#### Non parametrics ####
mann_whitney <- function(data_v,category)	{
ranked <- rank(data_v)
cats <- unique(category)
n1 <- sum(category==cats[1])
n2 <- sum(category==cats[2])
nt <- length(category)
R1 <- sum(ranked[(1:nt)[category==cats[1]]])
U1 <- R1 -((n1*(n1+1))/2)
R2<- sum(ranked[(1:nt)[category==cats[2]]])
U2 <- R2 -((n2*(n2+1))/2)
mU <- (n1*n2)/2
unique_ranks <- sort(unique(ranked))
numr <- 0
for (i in 1:length(unique_ranks))	{
	ti <- sum(ranked %in% unique_ranks[i])
	numr <- numr + ((ti^3)-ti)
	}
numr <- numr/(nt*(nt-1))
sigU <- ((n1*n2/12)*(nt+1-numr))^0.5	# expected standard deviation in ranks corrected for ties
#(n1*n2*(n1+n2+1)/12)^0.5
if (abs(U2-mU)>0)	{
	pval <- 2*integrate(dnorm,lower=abs(U2-mU)/sigU,upper=10)$value;
	} else	{
	pval <- 1.0;  # if we have only the same numbers in both categories, then it's a perfect match
	}
output <- data.frame(U=as.numeric(max(U2,U1)),pval=as.numeric(pval),stringsAsFactors = F);
return(output)
}

#### Probability 101 ####
Poisson_rate_to_probability <- function(expectation)	{
return(1-exp(-expectation))
}

probability_to_Poisson_rate <- function(pn)	{
return(-1*log(1-pn))
}

binomial_support_bars_ext <- function (n,N,support=1.0)	{
ml_Rate <- n/N;
rates <- seq(0.001,0.999,by=0.001);
rate_loglikelihood <- log(dbinom(n,N,rates));
max_lnL <- max(rate_loglikelihood);
rate_support <- rate_loglikelihood-max_lnL;
ok_rates <- rates[rate_support>=-support];
return(c(min(ok_rates),max(ok_rates)));
}

#### Row & Column Statistics ####
colMax <- function(data_matrix)	{
col_maxes <- c();
for (i in 1:ncol(data_matrix))	col_maxes <- c(col_maxes,max(data_matrix[,i]));
return(col_maxes)
}
	
colMin <- function(data_matrix)	{
col_mins <- c();
for (i in 1:ncol(data_matrix))	col_mins <- c(col_mins,min(data_matrix[,i]));
return(col_mins)
}
	
# this works like rowSums
rowProds <- function(dmat,use_logs=T)	{
if (ncol(dmat)==1)	return(dmat);
rps <- c();
if (use_logs)	{
	lndmat <- log(dmat);
	rps <- exp(rowSums(lndmat));
	} else	{
	for (i in 1:nrow(dmat))
		rps <- c(rps,prod(dmat[i,]));
	}
return(rps);
}

colProds <- function(dmat,use_logs=T)	{
if (nrow(dmat)==1) return(dmat);
rps <- c();
if (use_logs)	{
	lndmat <- log(dmat);
	rps <- exp(colSums(lndmat));
	} else	{
	for (i in 1:nrow(dmat))
		rps <- c(rps,prod(dmat[,i]));
	}
return(rps);
}

colMedians <- function(data_matrix)	{
col_medians <- c();
for (i in 1:ncol(data_matrix))	col_medians <- c(col_medians,median(data_matrix[,i]));
return(col_medians)
}

rowMedians <- function(data_matrix)	{
row_medians <- c();
for (i in 1:nrow(data_matrix))	row_medians <- c(row_medians,median(data_matrix[i,]));
if (!is.null(rownames(data_matrix)))	names(row_medians) <- rownames(data_matrix);
return(row_medians)
}

rowMins <- function(m,exclude="")	{
rmn <- c();
for (i in 1:nrow(m))
	rmn <- c(rmn,min(m[i,!m[i,] %in% exclude]));
return(rmn);
}

rowMaxs <- function(m,exclude="")	{
rmx <- c();
for (i in 1:nrow(m))
	rmx <- c(rmx,max(m[i,!m[i,] %in% exclude]));
return(rmx);
}

colMins <- function(m)	{
m <- clear_na_from_matrix(m,-MAXNO);
cmn <- c();
for (i in 1:ncol(m))
	cmn <- c(cmn,min(m[!m[i,] %in% exclude,i]));
return(cmn);
}

colMaxs <- function(m)	{
m <- clear_na_from_matrix(m,-MAXNO);
cmx <- c();
for (i in 1:ncol(m))
	cmx <- c(cmx,max(m[,i]));
return(cmx);
}

rowVars <- function(m)	{
row_variances <- vector(length=nrow(m));
for (i in 1:nrow(m))	row_variances[i] <- var(simplify2array(m[i,]));
return(row_variances);
}

colVars <- function(m)	{
col_variances <- vector(length=ncol(m));
for (i in 1:ncol(m))	col_variances[i] <- var(m[,i]);
ncol(col_variances);
}

colcumsum <- function(m,descending=T)	{
cumulative_m <- m;
mm <- nrow(m);
for (i in 1:ncol(m))	{
	if (descending)	{
		cumulative_m[,i] <- cumsum(m[,i]);
		} else	{
		cumulative_m[,i] <- cumsum(m[mm:1,i]);
		}
	}
return(cumulative_m);
}

rowcumsum <- function(m,descending=T)	{
cumulative_m <- m;
mm <- ncol(m);
for (i in 1:nrow(m))	{
	if (descending)	{
		cumulative_m[i,] <- cumsum(m[i,]);
		} else	{
		cumulative_m[i,] <- cumsum(m[i,mm:1])[mm:1];
		}
	}
return(cumulative_m);
}

#### Model Fitting ####
accersi_best_polynomial <- function (xx,yy)	{
# evaluated using residuals after Royall 1997; see Finarelli (2006)
n <- length(xx);

model_b <- lm(yy ~ xx);
best_model_summary <- model_summary <- summary(model_b);
ESS <- sum(model_summary$residuals^2);
loglikelihood <- best_loglikelihood <- (n*log(ESS/n)/-2);

pl <- 1;
improvement <- 1000;
while (improvement > 2)	{
	pl <- pl+1;
	model_t <- lm(yy ~ poly(xx,pl,raw=T));
	model_summary <- summary(model_t);
	ESS <- sum(model_summary$residuals^2);
	new_loglikelihood <- (n*log(ESS/n)/-2);
	loglikelihood <- c(loglikelihood,new_loglikelihood);
	improvement <- new_loglikelihood-best_loglikelihood;
	if (improvement>2)	{
		best_loglikelihood <- new_loglikelihood;
		model_b <- model_t;
		best_model_summary <- model_summary;
		}
	}
sxx <- sort(xx);
expectations <- best_model_summary$coefficients[1];
for (i in 2:pl)
	expectations <- expectations+(best_model_summary$coefficients[i]*(sxx^(i-1)));
output <- list(best_model_summary,loglikelihood,expectations);
names(output) <- c("model_summary","loglikelihood","expectations");
return(output);
}

two_sided_credible_intervals <- function(posteriors,cred_int=0.95)	{
peak <- match(max(posteriors),posteriors);
best_posteriors <- sort(posteriors,decreasing = T);
interval_credibility <- array(0,dim=length(best_posteriors));
for (n in 1:length(posteriors))	{
	h <- match(best_posteriors[n],posteriors);
	interval_credibility[h] <- max(interval_credibility)+posteriors[h];
	}

credible_bounds <- (1:length(interval_credibility))[interval_credibility < 1-(1-cred_int)/2];
return(credible_bounds);
}


#### Regression ####
#x <- -abs(branch_durations$mid[2:nrow(branch_durations)])
#y <- log(best_rates_per_branch[2:nrow(branch_durations)]);
reduced_major_axis_regression <- function(x,y)	{
rma_decay <- sqrt(var(y))/sqrt(var(x));
dd <- lm(y~x);
if (dd$coefficients[2]<0)	rma_decay <- -1*rma_decay;

mean_x <- mean(x);
mean_y <- mean(y);

y_str <- mean_y+(min(x)-mean_x)*rma_decay;
y_end <- mean_y+(max(x)-mean_x)*rma_decay;
intercept <- mean_y - (mean_x*rma_decay);

#rma <- data.frame(slope=as.numeric(rma_decay),intercept=as.numeric(intercept));
rma <- data.frame(slope=as.numeric(rma_decay),mean_x=as.numeric(mean_x),mean_y=as.numeric(mean_y),stringsAsFactors = F);
return(rma);
}
#### Information Theory 101 ####
# first-order AIC
AIC <- function(lnL,k)	{
# L: log-likelihood; # k: parameters; # n: data points
return((-2*lnL)+(2*k))
}

# second-order (modified) AIC
modified_AIC <- function(lnL,k,n)	{
# L: log-likelihood; # k: parameters; # n: data points
#if (n==k)	n <- n+2
if (is.na(n / (n - k - 1)) || (n / (n - k - 1)<0) || (n / (n - k - 1))==Inf)	{
	aic_c <- AIC(lnL,k);
	}	else	aic_c <- (-2*lnL) + (2*k)*(n / (n - k - 1));
return(aic_c);
}

#aics <- aics[1,]
akaikes_weights <- function(aics)	{
delta_aics <- aics-min(aics);
return(exp(-0.5 *delta_aics)/sum(exp(-0.5 *delta_aics)));
}
#### Odds & Ends
geomean <- function(x)	{
lnx <- log(x);
return(exp(mean(lnx)));
}