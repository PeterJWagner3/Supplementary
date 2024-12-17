Jaccard_Similarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (0.0);
	} else		return (shared/(shared+S1+S2));
}

Sorensen_Similarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (0.0);
	} else	{
	return ((2*shared)/((2*shared)+S1+S2));
	}
}

Simpson_Similarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (0.0);
#	} else	return (shared/(shared+min(S1,S2)));
	} else	return (shared/min(S1,S2));
}

Nested_Resultant_Similarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
x=abs(S1-S2)/((2*shared)+S1+S2);
y=shared/(shared+min(S1,S2));
if (shared==0)	{
	return (0.0);
	} else	{
	return(1-(x*y));
	}
#return(1-((abs(S1-S2)/((2*shared)+S1+S2))*(shared/(shared+min(S1,S2)))));
#return (1-(((Smax-Smin)/((2*Sshared))) * (Sshared/(Sshared-Smin))));
}

Lennon_Similarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (0.0);
	} else	{
	return (1-(2*abs(S1-S2)/((2*shared)+S1+S2)));
	}
#return (1-(2*(Smax-Smin)/((2*Sshared))));
}

Jaccard_Disimilarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (1.0);
	} else	{
	return ((S1+S2)/(shared+S1+S2));
	}
}

Sorensen_Disimilarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (1.0);
	} else	{
	return ((S1+S2)/((2*shared)+S1+S2));
	}
}

Simpson_Disimilarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (1.0);
	} else	{
	return (min(S1,S2)/min(S1,S2));
	}
}

NestedResultant_Disimilarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);

x=abs(S1-S2)/((2*shared)+S1+S2);
y=shared/(shared+min(S1,S2));
if (shared==0)	{
	return (1.0);
	} else	{
	return(x*y);
	}
#return((abs(S1-S2)/((2*shared)+S1+S2))*(shared/(shared+min(S1,S2))));
#return (((Smax-Smin)/((2*Sshared))) * (Sshared/(Sshared-Smin)));
}

Lennon_Disimilarity <- function(shared, S1, S2)	{
shared=max(shared,0);
S1=max(S1,0);
S2=max(S2,0);
if (shared==0)	{
	return (1.0);
	} else	{
	return (2*abs(S1-S2)/((2*shared)+S1+S2));
	}
#return (2*(Smax-Smin)/((2*Sshared)));
}

chao_shared_richness <- function(tad1, tad2) {
#####################################################################
#	tad1: taxon occurrences in Bin X
#	tad2: taxon occurrences in Bin Y
#	spid1: taxon numbers for taxa in Bin X
#	spid2: taxon numbers for taxa in Bin Y
#	S1: observed taxa in Bin X
#	S2: observed taxa in Bin Y
#	N1: for Bin X, total finds if abundance data or total collections if locality data
#	N2: for Bin Y, total finds if abundance data or total collections if locality data
#	variables follow Chao et al. 2005 insofar as possible
#########################################################
spid1 <- names(tad1);
spid2 <- names(tad2);
S1 <- length(tad1);
S2 <- length(tad2);
N1 <- sum(tad1);
N2 <- sum(tad2);
shared <- X <- Y <- vector(length=S1);		# finds in list X for shared taxa	#
shared_taxa <- spid1[spid1 %in% spid2];
fj1 <- sum(tad2[names(tad2) %in% shared_taxa]==1);			# taxa from list 1 found once in list 2	#
fj2 <- sum(tad2[names(tad2) %in% shared_taxa]==2);			# taxa from list 1 found twice in list 2	#
f1k <- sum(tad1[names(tad1) %in% shared_taxa]==1);			# taxa from list 2 found once in list 2	#
f2k <- sum(tad1[names(tad1) %in% shared_taxa]==2);			# taxa from list 2 found twice in list 1	#
f11 <- n;			# taxa found once each in both lists 1 & 2	#
f22 <- n;			# taxa found twice each in both lists 1 & 2	#
D12 <- length(shared_taxa); # observed shared
for (a in 1:S1)	{
	for (b in 1:S2)	{
		if (spid1[a]==spid2[b])	{
			shared[D12]=spid1[a];
			X[D12]=tad1[a];		# sum all of the finds for shared taxa from list X	#
			Y[D12]=tad2[b];		# sum all of the finds for shared taxa from list Y	#
			# tally shared species that are singletons in one or both	#
			if (Y[D12]==1)	{
				fj1 <- fj1+1;			# taxa from list X found once in list Y	#
				if (X[D12]==1)	f11 <- f11+1;
				} else if (Y[D12]==2)	{
				fj2 <- fj2+1;			# taxa from list X found twice in list Y	#
				if (X[D12]==2)	f22 <- f22+1;
				}
			if (X[D12]==1)	{
				f1k <- 1+f1k;			# taxa from list Y found once in list X
				} else if (X[D12]==2) {
				f2k <- 1+f2k;			# taxa from list Y found twice in list X
				}
			b  <- S2;	# terminate search: it's been found!	#
			D12 <= 1+D12;
			}	# end search for match	#
		}	# end search to see if bin B has a match with bin A	#
	}

#if (fj1==0)	fj1=1;			# make sure these numbers all are at least 1	#
if (fj2==0)	fj2 <- 1;
#if (f1k==0)	f1k=1;
if (f2k==0)	f2k <- 1;
#if (f11==0)	f11=1;
if (f22==0)	f22 <- 1;

K1 <- (N1-1)/N1;
K2 <- (N2-1)/N2;

S12 <- D12 + ((K1*f1k*f1k)/(2*f2k))+((K2*fj1*fj1)/(2*fj2))+(((K1*K2*f11*(f11-1)))/(4*(f22+1)));

return (S12);
# find shared taxa	#
}


chao_shared_UV <- function(list1, list2)	{
#int a, b;
#int D12=0;
#double fx1=0.0f, fx2=0.0, f1y=0.0f, f2y=0.0f;
#double sumXi=0, sumYi=0, sumXiOne=0, sumYiOne=0;
#double nx, ny;
#//double U, V;
#double *X, *Y;				/* arrays of finds for shared species in X and Y			*/
#double *UV;
#double ii, ij, ik, im, ji, jj, jk, jm;

#/* U = Prop. of List X occurrences that are co-occurrences + (([lists in Y] - 1)/[lists in Y]) x Prop. of occurrences in List X that are singletons in Y	*/

UV <- vector(length=2);

#UV=dvector(2);
#//shared=ivector(S1);
X=dvector(S1);		/* finds in list X for shared taxa	*/
Y=dvector(S1);		/* finds in list Y for shared taxa	*/
for (a=0; a<S1; ++a)	{
	for (b=0; b<S2; ++b)	{
		if (spid1[a]==spid2[b])	{
//			shared[D12]=spid1[a];
			X[D12]=tad1[a];					/* sum all of the finds for shared taxa from list X	*/
			Y[D12]=tad2[b];					/* sum all of the finds for shared taxa from list Y	*/
			if (Y[D12]==1)	++fx1;			/* taxa from list X found once in list Y	*/
			else if (Y[D12]==2)	++fx2;		/* taxa from list X found twice in list Y	*/
			if (X[D12]==1)	++f1y;			/* taxa from list Y found once in list X	*/
			else if (X[D12]==2)	++f2y;		/* taxa from list Y found twice in list X	*/

			++D12;
			b=S2;	/* terminate search: it's been found!	*/
			}	/* end search for match	*/
		}
	}

nx=sumulvector(tad1,S1);	/* total finds in X	*/
ny=sumulvector(tad2,S2);	/* total finds in Y	*/

if (f2y==0)	f2y=1.0f;
if (fx2==0)	fx2=1.0f;

/* U = [proportion of finds in X from species also found in Y] + [(Ycoll - 1)/Ycoll] x [shared Y singletons / (2 x shared Y doubletons)] x [prop. of finds in X from Y singletons]	*/
ii=sumdvector(X,D12);
ii/=nx;										/* proportion of occurrences in X from species shared with Y		*/
ij=((double) (z-1))/((double) z);			/* (max poss. occurrences in Y - 1) / max possible occurrences in Y	*/
ik=fx1/(2*fx2);								/* Y singletons / (2 x Y doubletons) assuming shared with X			*/
im=0.0f;									/* prop. of occurrences in X from singletons in Y					*/
for (a=0; a<D12; ++a)	if (Y[a]==1)	im+=X[a]/((double) nx);
UV[0]= ii + (ij*ik*im);

/* V = [proportion of finds in Y from species also found in X] + [(Xcoll - 1)/Xcoll] x [shared X singletons / (2 x shared X doubletons)] x [prop. of finds in Y from X singletons]	*/
ji=sumdvector(Y,D12);
ji/=ny;
jj=((double) (w-1))/((double) w);
jk=f1y/(2*f2y);
jm=0.0f;
for (a=0; a<D12; ++a)	if (X[a]==1)	jm+=Y[a]/((double) ny);
UV[1]= ji + (jj*jk*jm);

return (UV);
/* find shared taxa	*/
}


ord_example_function <- function()	{
ord_example <- read.csv(file.choose(),header=T,stringsAsFactors = F);
if (min(ord_example$formation)<1)	ord_example$formation <- ord_example$formation+1-min(ord_example$formation);
list_nos <- unique(ord_example$formation);
for (ln in 1:length(list_nos))	{
	this_list <- subset(ord_example,ord_example$formation==list_nos[ln]);
	these_taxa <- unique(this_list$taxon);
	these_finds <- c();
	for (tt in 1:length(these_taxa))	these_finds <- c(these_finds,sum(this_list$taxon %in% these_taxa[tt]));
	dummy <- data.frame(site=rep(ln,length(these_taxa)),taxon=as.numeric(these_taxa),finds=as.numeric(these_finds));
	dummy <- dummy[order(-dummy$finds,dummy$taxon),];
	if (ln==1)	{
		ord_lists <- dummy;
		} else	{
		ord_lists <- rbind(ord_lists,dummy);
		}
	}
list1 <- subset(ord_lists,ord_lists$site==6)
list2 <- subset(ord_lists,ord_lists$site==7)
#chao_shared_richness <- function(tad1, tad2, spid1, spid2, S1, S2, n1, n2)	{
chao_shared_richness <- function(list1, list2)	{
# list$taxon: taxon names or numbers
# list$finds: taxon finds;
if (nrow(list2) > nrow(list1) || (nrow(list2)==nrow(list1) && sum(list2$finds)>sum(list1$finds)))	{
	dlist <- list1;
	list1 <- list2;
	list2 <- dlist;
	}
obs_shared <- list1$taxon[list1$taxon %in% list2$taxon];
Ss_obs <- length(obs_shared);
list1_sh <- list1[list1$taxon %in% obs_shared,];
list2_sh <- list2[list2$taxon %in% obs_shared,];
list2_sh <- list2_sh[order(list2_sh$taxon,list1_sh$taxon),];
fk1 <- sum(list2_sh$finds==1);	# taxa on 1st list that appear once in 2nd list
fk2 <- sum(list2_sh$finds==2);	# taxa on 1st list that appear twice in 2nd list
fj1 <- sum(list1_sh$finds==1);	# taxa on 2nd list that appear once in 1st list
fj2 <- sum(list1_sh$finds==2);	# taxa on 2nd list that appear twice in 1st list
f11 <- sum((list1_sh$finds==1)*(list2_sh$finds==1));
f22 <- sum((list1_sh$finds==2)*(list2_sh$finds==2));

#if (fj1==0)	fj1=1;			#make sure these numbers all are at least 1	*/
if (fj2==0)	fj2 <- 1;
#if (f1k==0)	f1k=1;
if (f2k==0)	f2k <- 1;
#if (f11==0)	f11=1;
if (f22==0)	f22 <- 1;

K1 <- (w-1)/w;
K2 <- (z-1)/z;

S12 <- ((K1*f1k*f1k)/(2*f2k))+((K2*fj1*fj1)/(2*fj2))+(((K1*K2*f11*(f11-1)))/(4*(f22+1)));

return (Ss_obs+S12);
#find shared taxa	*/
}
}

