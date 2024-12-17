# generate geometric (= exponential) distribution
geometric_distribution <- function(decay)	{
if (decay>1)	decay <- (1/decay)
S <- round(1+((log(10^-9)-log(1-decay))/log(decay)))
ranks <- (1:S)
prop <- (1-decay)*(decay^(ranks-1))
return(prop)
}

# generate geometric (= exponential) distribution
exponential_distribution <- function(decay)	{
#print(decay);
d <- min(decay,1/decay)
S <- accersi_richness_with_prop_greater_than_x_given_exponential(decay,MINEXPN)
ranks <- (1:S)
prop <- (1-d)*(d^(ranks-1))
return(prop)
}

# estimate log-series distribution given alpha and N
logseries_distribution <- function(alpha,J)	{
# J is the true (incompletely sampled) population size
# alpha = Fisher's Biodiversity Index
prob_n_finds <- dls((1:J),J,alpha)	# dls from sads
N <- 0
S <- 2
while (N<J)	{
	S <- S+1
	N <- N_given_S_and_Prob_N(prob_n_finds,S)
	}
taxa_w_n_finds <- S*prob_n_finds;
sad <- convert_expected_to_relative_abundance(taxa_w_n_finds);
return (sad);
}

# calcualate lognormal for S entities that increase in magnitude by mag every standard deviation
lognormal_distribution <- function(mag, S)	{
fi <- (S:1)/(S+1);
prop <- mag^(qnorm(fi,0,1))/sum(mag^(qnorm(fi,0,1)));
return(prop);
}

# generate Zipf distribution for S individuals with log-decay = zg
zipf_distribution <- function(zg, S)	{
return(((1:S)^-zg)/sum((1:S)^-zg));
#zpf <- ranks^-zg
#prop <- zpf/sum(zpf)
#return(prop)
}

# generate Zipf-Mandelbrot distribution for S individuals with log-decay = zg
zipf_mandelbrot_distribution <- function(zg,zb,S)	{
return((1/((1:S)+zb)^zg)/sum(1/((1:S)+zb)^zg))
}

# get gamma distribution
gamma_distribution <- function(a, b, S)	{
p0 <- vector(length=S)
p0[1] <- S/(S+1)
for(i in 2:S) p0[i] <- p0[i-1]-(1/(S+1))
prop <- qgamma(p0,a,b)/sum(qgamma(p0,a,b))
return(prop)
}

# get beta distribution
beta_distribution <- function(shape1, shape2, S)	{
# NOTE:I get some wonky assed results using dbeta & qbeta
ranks <- (1:S)/(S+1)
numer <- (ranks^(shape1-1))*((1-ranks)^(shape2-1))
#numer[numer==0] <- ZERO
# use log gammas, as gamma does not like numbers much past 100
Bab <- exp((lgamma(shape1)+lgamma(shape2))-lgamma(shape1+shape2))
cdf <- numer/Bab
pdf <- cdf/sum(cdf)
return(pdf)
}
