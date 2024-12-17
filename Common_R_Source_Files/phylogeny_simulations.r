T <- TRUE;

# Simulate Tree Evolution ####
simulate_cladogram <- function(otu,tree_type=c("vector","matrix"))	{
htu <- otu+1
branch_durations <- vector(length=(2*otu)-1)
branch_durations[1:2] <- 1
if (tree_type=="matrix")	{
	cladogram <- c(1,2)
	for (s in 3:otu)	{
		htu <- htu+1
		split <- ceiling(runif(1)/(1/(s-1)))
		new_cl <- c(split,s)
		cladogram <- rbind(cladogram,new_cl)
		# find where the tree needs to be altered
		alter <- which(cladogram==split,arr.ind=TRUE)
		alter <- alter[order(alter[,1],decreasing=FALSE),]
		# replace split species with htu number
		cladogram[alter[1,1],alter[1,2]] <- htu
		cladogram[alter[1,1],] <- sort(cladogram[alter[1,1],])
		branch_durations[htu] <- branch_durations[split]
		branch_durations[split] <- 0
		branch_durations[array(cladogram)[array(cladogram)<=otu]] <- branch_durations[array(cladogram)[array(cladogram)<=otu]]+1
		}
	} else	{
	cladogram <- rep(-1,(2*otu)-1)
	cladogram[1:2] <- htu
	branch_durations <- vector(length=(2*otu)-1)
	branch_durations[1:2] <- 1
	for (s in 3:otu)	{
		htu <- htu+1
		split <- ceiling(runif(1)/(1/(s-1)))
		cladogram[htu] <- cladogram[split]
		cladogram[split] <- cladogram[s] <- htu
		branch_durations[htu] <- branch_durations[split]
		branch_durations[split] <- 0
		branch_durations[1:s] <- branch_durations[1:s]+1
#		s <- s+1
		}
	}
output <- list(cladogram,branch_durations)
names(output) <- c("Cladogram","Rel_Branch_Lengths")
return (output)
}

# routine to evolve S contemporaneous (e.g., extant) taxa
evolve_to_standing_richness_S <- function(S,lambda,mu,bifurcation=F,temp_prec=0.1)	{
# cumulative and standing diversity
spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
dbr <- 1
while (max(standing) < S)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	if (daughters > 0)	{
		if (bifurcation)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			remaining <- c(remaining,cumulative+(1:daughters))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		}
	}

sim_species_present <- accersi_species_present_in_bins(ranges,temp_prec=0.1)
if (!is.na(match(S,standing)))	{
	present <- match(S,standing)
	}	else	{
	present <- match(max(standing),standing)
	}
sampled <- sim_species_present[present,]

if (length(sampled)>S)	{
	xx <- sampled[order(birth[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral)
venn_tree <- venn_tree_info$Venn_Tree
branchings <- venn_tree_info$Branchings
branchings[S+1] <- 0
vtree <- transform_venn_tree_to_vector_tree(venn_tree)
divergence_dates <- birth[venn_tree_info$First_Taxon]
divergence_dates[S+1] <- max(divergence_dates)
divergence_dates <- divergence_dates-(min(divergence_dates))
divergence_dates[S+1] <- 0
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates)
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates")
return(output)
}

# routine to evolve S taxa sampled from throughout the clade's history
evolve_to_sampled_richness_S_modified <- function(S,p,q,psi=0,freqrat=-1,anagenesis=0,pulsed_turnover=F,bifurcation=F,temp_prec=0.05)	{
# S: sampled richness to evolve
# r: birth rate
# q: death rate
# psi: sampling rate
# freqrat: proportion of taxa sampled
# anagenesis: rate of anagenesis relative to death rate.
#	pulsed_turnover: if true, then cladogenesis & extinction happen at integers immitating stages
#	bifurcation: if true, then cladogenesis is bifurctating; otherwise, it is budding
#	temp_prec: the precision of originations & extinctions

# if freqrat is provided rather than sampling rate:
if (freqrat>0 && psi==0)	psi <- q/((1/freqrat)-1);

lambda <- p/(1+anagenesis);	# cladogenesis rate
mu <- q/(1+anagenesis);			# extinction rate
darwin <- anagenesis*mu;		# anagenesis rate;
maxfinds <- dbr <- spc <- cumulative <- standing <- 1;
birth <- 0.0;
#life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
if (pulsed_turnover)	duration <- ceiling(duration);
death <- c(birth[spc]+duration);
ancestral <- c(0);
remaining <- c(1);
found <- 0;	# this will keep track of how many species we have sampled
sampled <- c();	# list of sampled species
fa <- la <- occurrences <- vector(length=1)	# occurrences per species
all_finds <- array(0,dim=c(1,1));
#all_finds <- array(finds,dim=c(1,length(finds)));
keep_going <- TRUE;
while (keep_going)	{
	spc <- remaining[1];
	if (nrow(all_finds)<spc)	all_finds <- rbind(all_finds,array(0,dim=c(1,ncol(all_finds))));
#	print(c(dbr,spc)); # for debugging.
	duration <- abs((death[spc]-birth[spc]));
	mini_mes <- rpois(1,darwin*duration);	# number of anagenetic pseudospeciations/pseudoextinctions
	if (mini_mes>0)	{
		regens <- birth[spc]+sort(runif(mini_mes))*duration;  # place anagenetic events;
		new_ends <- c(regens,death[spc]);
		chuckies <- (cumulative+1):(cumulative+mini_mes);	# number anagenetics
		added <- rep(0,daughters)
		occurrences <- c(occurrences,added)
		ancestral <- c(ancestral,rep(spc,mini_mes));
		birth <- c(birth,regens);
		death[spc] <- regens[1];
		death <- c(death,new_ends[2:length(new_ends)]);
		duration <- death[spc]-birth[spc];		# duration of morphospecies
		remaining <- c(remaining,cumulative+(1:mini_mes));	# add anagenetics to the queue
		cumulative <- cumulative+mini_mes;
		added <- rep(0,mini_mes)
		occurrences <- c(occurrences,added);
		fa <- c(fa,added)
		la <- c(la,added)
		all_finds <- rbind(all_finds,array(0,dim=c(mini_mes,maxfinds)));
#		rownames(all_finds) <- 1:cumulative;
		rownames(all_finds) <- 1:nrow(all_finds);
		}
	daughters <- rpois(1,(lambda*duration));	# cladogenetic daugthers;
	if (daughters > 0)	{
		if (bifurcation)	{
			splits_dates <- birth[spc]+sort(duration*runif(daughters));
			if (pulsed_turnover)
				split_dates <- floor(splits_dates);
			splits <- temp_prec*ceiling(splits_dates/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc;
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
				death <- c(death,(splits[i] + duration));
				if (pulsed_turnover)
					death[length(death)] <- ceiling(death[length(death)]);
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (fossils > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,array(0,dim=c(daughters*2,maxfinds)));
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						added <- max(occurrences)-maxfinds;
						n <- spc-1;
						all_finds <- cbind(all_finds,matrix(0,n,added));
#						if (nrow(all_finds)>n)
#							all_finds <- all_finds[1:n,];
#						all_finds <- cbind(all_finds,rep(0,length(occurrences)));
						maxfinds <- max(occurrences);
#						all_finds <- cbind(all_finds,matrix(0,length(occurrences),added));
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap);
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f];
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f];
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
#			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			birth <- c(birth,birth[spc]+sort(temp_prec*ceiling((duration*runif(daughters))/temp_prec)))
			ancestral <- c(ancestral,rep(spc,daughters));
			for (d in 1:daughters)	{
				duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+duration
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,array(0,dim=c(daughters,maxfinds)));
			rownames(all_finds) <- 1:nrow(all_finds);
			}
		cumulative <- max(remaining);
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
#	if ((daughters == 0 || !bifurcation) && fossils > 0)	{
	fossils <- rpois(1,psi*duration);					# finds for morphospecies;
	if (fossils > 0)	{
#	if (fossils > 0)	
		finds <-  sort(birth[spc]+runif(fossils)*duration);
		# we found the beastie
		occurrences[spc] <- fossils;
		fa[spc] <- min(finds);
		la[spc] <- max(finds);
		sampled <- c(sampled,spc);
		if (maxfinds < occurrences[spc])	{
			added <- occurrences[spc] - maxfinds;
			if (spc>1)	{
				all_finds <- cbind(all_finds,array(0,dim=c(nrow(all_finds),added)));
#				all_finds <- all_finds[1:(spc-1),];
#				all_finds <- rbind(all_finds,finds);
				} else	{
				all_finds <- array(finds,dim=c(1,length(finds)));
				}
#			if (length(finds)>0)
			all_finds[spc,] <- finds;
			maxfinds <- max(occurrences);
			} else	{
			all_finds[spc,1:length(finds)] <- finds;
			}
		}# else	{
#		occurrences[spc] <- 0;
#		all_finds[spc,] <- rep(0,ncol(all_finds));
#		}
	found <- length(sampled);
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]; # remove species from queue
		}	else if (found < S)	{
		dbr <- dbr+1;
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- 0.0;
		duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
		death <- c(birth[spc]+duration)
		if (pulsed_turnover)	{
			death <- ceiling(death);
			duration <- death-birth;
			} # case for pulsed turnovers 
		ancestral <- c(0);
		remaining <- c(1);
		# clear sampling vectors
		occurrences <- sampled <- c()	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1;
		all_finds <- array(0,dim=c(1,1));
		}
	found <- length(sampled);
	relv_fas <- sort(fa[sampled])
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE);
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral);
venn_tree <- venn_tree_info$Venn_Tree;
branchings <- venn_tree_info$Branchings;
branchings[S+1] <- 0;
vtree <- transform_venn_tree_to_vector_tree(venn_tree);
divergence_dates <- birth[venn_tree_info$First_Taxon];
if (sampled[1]==1)	{
	divergence_dates <- divergence_dates-divergence_dates[1];
	}	else	{
	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1);
	}
divergence_dates[S+1] <- 0;
durations <- cbind(birth[sampled],death[sampled]);
strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec));
#print(branchings[1:S]);
#print(divergence_dates);
otu_finds <- all_finds[sampled,colSums(all_finds[sampled,])>0]
original <- list(ancestral,birth,death,occurrences,all_finds);
names(original) <- c("ancestor_raw","birth","death","nfinds","finds");
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled],otu_finds,original);
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","No_Finds","Finds","Full_Simulation");
return(output);
}

# routine to evolve S taxa sampled from throughout the clade's history
evolve_to_sampled_richness_S <- function(S,p,q,psi=0,freqrat=-1,anagenesis=0,pulsed_turnover=F,bifurcation=F,temp_prec=0.1)	{
# S: sampled richness to evolve
# r: birth rate
# q: birth rate
# psi: sampling rate	
# freqrat: proportion of taxa sampled
# anagenesis: anagenesis rate
if (freqrat>0 && psi==0)	psi <- q/((1/freqrat)-1);

lambda <- p/(1+anagenesis);	# cladogenesis rate
mu <- q/(1+anagenesis);			# extinction rate
darwin <- anagenesis*mu;		# anagenesis rate;
maxfinds <- dbr <- spc <- cumulative <- standing <- 1;
birth <- 0.0;
#life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
if (pulsed_turnover)	duration <- ceiling(duration);
death <- c(birth[spc]+duration);
ancestral <- c(0);
remaining <- c(1);
found <- 0;	# this will keep track of how many species we have sampled
sampled <- c();	# list of sampled species
fa <- la <- occurrences <- vector(length=1)	# occurrences per species
all_finds <- array(0,dim=c(1,1));
#all_finds <- array(finds,dim=c(1,length(finds)));
keep_going <- T;
while (keep_going)	{
	spc <- remaining[1];
	if (nrow(all_finds)<spc)	all_finds <- rbind(all_finds,array(0,dim=c(1,ncol(all_finds))));
#	print(c(dbr,spc));
	duration <- abs((death[spc]-birth[spc]));
	mini_mes <- rpois(1,darwin*duration);
	if (mini_mes>0)	{
		regens <- birth[spc]+sort(runif(mini_mes))*duration;  # place anagenetic events;
		new_ends <- c(regens,death[spc]);
		chuckies <- (cumulative+1):(cumulative+mini_mes);	# number anagenetics
		added <- rep(0,daughters)
		occurrences <- c(occurrences,added)
		ancestral <- c(ancestral,rep(spc,mini_mes));
		birth <- c(birth,regens);
		death[spc] <- regens[1];
		death <- c(death,new_ends[2:length(new_ends)]);
		duration <- death[spc]-birth[spc];		# duration of morphospecies
		remaining <- c(remaining,cumulative+(1:mini_mes));	# add anagenetics to the queue
		cumulative <- cumulative+mini_mes;
		added <- rep(0,mini_mes)
		occurrences <- c(occurrences,added);
		fa <- c(fa,added)
		la <- c(la,added)
		all_finds <- rbind(all_finds,array(0,dim=c(mini_mes,maxfinds)));
#		rownames(all_finds) <- 1:cumulative;
		rownames(all_finds) <- 1:nrow(all_finds);
		}
	daughters <- rpois(1,(lambda*duration));	# cladogenetic daugthers;
	if (daughters > 0)	{
		if (bifurcation)	{
			splits_dates <- birth[spc]+sort(duration*runif(daughters));
			if (pulsed_turnover)
				split_dates <- floor(splits_dates);
			splits <- temp_prec*ceiling(splits_dates/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc;
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
				death <- c(death,(splits[i] + duration));
				if (pulsed_turnover)
					death[length(death)] <- ceiling(death[length(death)]);
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (fossils > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,array(0,dim=c(daughters*2,maxfinds)));
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						added <- max(occurrences)-maxfinds;
						n <- spc-1;
						all_finds <- cbind(all_finds,matrix(0,n,added));
#						if (nrow(all_finds)>n)
#							all_finds <- all_finds[1:n,];
#						all_finds <- cbind(all_finds,rep(0,length(occurrences)));
						maxfinds <- max(occurrences);
#						all_finds <- cbind(all_finds,matrix(0,length(occurrences),added));
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap);
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f];
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f];
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
#			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			birth <- c(birth,birth[spc]+sort(temp_prec*ceiling((duration*runif(daughters))/temp_prec)))
			ancestral <- c(ancestral,rep(spc,daughters));
			for (d in 1:daughters)	{
				duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+duration
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,array(0,dim=c(daughters,maxfinds)));
			rownames(all_finds) <- 1:nrow(all_finds);
			}
		cumulative <- max(remaining);
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
#	if ((daughters == 0 || !bifurcation) && fossils > 0)	{
	fossils <- rpois(1,psi*duration);					# finds for morphospecies;
	if (fossils > 0)	{
#	if (fossils > 0)	
		finds <-  temp_prec*ceiling(sort(birth[spc]+runif(fossils)*duration)/temp_prec);
		# we found the beastie
		occurrences[spc] <- fossils;
		fa[spc] <- min(finds);
		la[spc] <- max(finds);
		sampled <- c(sampled,spc);
		if (maxfinds < occurrences[spc])	{
			added <- occurrences[spc] - maxfinds;
			if (spc>1)	{
				all_finds <- cbind(all_finds,array(0,dim=c(nrow(all_finds),added)));
#				all_finds <- all_finds[1:(spc-1),];
#				all_finds <- rbind(all_finds,finds);
				} else	{
				all_finds <- array(finds,dim=c(1,length(finds)));
				}
			if (nrow(all_finds))
			all_finds[spc,] <- finds;
			maxfinds <- max(occurrences);
			}
		}# else	{
#		occurrences[spc] <- 0;
#		all_finds[spc,] <- rep(0,ncol(all_finds));
#		}
	found <- length(sampled);
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]; # remove species from queue
		}	else if (found < S)	{
		dbr <- dbr+1;
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- 0.0;
		duration <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
		death <- c(birth[spc]+duration)
		if (pulsed_turnover)	{
			death <- ceiling(death);
			duration <- death-birth;
			} # case for pulsed turnovers 
		ancestral <- c(0);
		remaining <- c(1);
		# clear sampling vectors
		occurrences <- sampled <- c()	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1;
		all_finds <- array(0,dim=c(1,1));
		}
	found <- length(sampled);
	relv_fas <- sort(fa[sampled])
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral);
venn_tree <- venn_tree_info$Venn_Tree;
branchings <- venn_tree_info$Branchings;
branchings[S+1] <- 0;
vtree <- transform_venn_tree_to_vector_tree(venn_tree);
#divergence_dates <- birth[venn_tree_info$First_Taxon];
#if (sampled[1]==1)	{
#	divergence_dates <- divergence_dates-divergence_dates[1];
#	}	else	{
#	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1);
#	}
#divergence_dates[S+1] <- 0;
divergence_dates <- accersi_simulated_divergence_times(sampled,ancestral,birth);
#divergence_dates[divergence_dates$onset!=divergence_dates$end,]
durations <- data.frame(birth=as.numeric(birth[sampled]),death=as.numeric(death[sampled]));
#strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec));
strat_ranges <- data.frame(fa=as.numeric(fa),la=as.numeric(la));
#print(branchings[1:S]);
#print(divergence_dates);
original <- list(ancestral,birth,death,occurrences);
names(original) <- c("ancestor_raw","birth","death","finds");
strat_ranges <- data.frame(fa=as.numeric(strat_ranges[sampled,1]),la=as.numeric(strat_ranges[sampled,2]));
#strat_ranges$la[occurrences[sampled]==1] <- strat_ranges$fa[occurrences[sampled]==1];
n <- cbind(strat_ranges,birth[sampled],death[sampled])
#n[strat_ranges$fa<birth[sampled],]
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled],original);
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","Finds","Full_Simulation");
return(output);
}

evolve_to_sampled_richness_S_simple <- function(S,lambda,mu,freqrat,bifurcation=F,temp_prec=0.1)	{
# cumulative and standing diversity
maxfinds <- dbr <- spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
found <- 0	# this will keep track of how many species we have sampled
sampled <- c()	# list of sampled species
fa <- la <- occurrences <- vector(length=1)	# occurrences per species
all_finds <- matrix(0,1,maxfinds)
keep_going <- TRUE
while (keep_going)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	fossils <- rpois(1,freqrat*life)
	if (fossils > 0)	finds <-  sort(birth[spc]+runif(fossils)*life)
	if (daughters > 0)	{
		if (bifurcation)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (fossils > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,matrix(0,daughters*2,maxfinds))
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						maxfinds <- occurrences[ap]
						all_finds <- cbind(all_finds,rep(0,length(occurrences)))
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap)
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f]
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f]
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
#			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			birth <- c(birth,birth[spc]+sort(temp_prec*ceiling((life*runif(daughters))/temp_prec)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,matrix(0,daughters,maxfinds))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if ((daughters == 0 || !bifurcation) && fossils > 0)	{
		occurrences[spc] <- fossils
		fa[spc] <- min(finds)
		la[spc] <- max(finds)
		sampled <- c(sampled,spc)
		if (maxfinds < occurrences[spc])	{
			added <- occurrences[spc] - maxfinds
			maxfinds <- occurrences[spc]
			all_finds <- cbind(all_finds,matrix(0,length(occurrences),added))
			all_finds[spc,] <- finds
			}	else	{
			all_finds[spc,1:fossils] <- finds	
			}
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		# clear sampling vectors
		sampled <- c()	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1
		all_finds <- matrix(0,1,maxfinds)
		}
	
	found <- length(sampled)
	relv_fas <- sort(fa[sampled])
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral)
venn_tree <- venn_tree_info$Venn_Tree
branchings <- venn_tree_info$Branchings
branchings[S+1] <- 0
vtree <- transform_venn_tree_to_vector_tree(venn_tree)
divergence_dates <- birth[venn_tree_info$First_Taxon]
if (sampled[1]==1)	{
	divergence_dates <- divergence_dates-divergence_dates[1]
	}	else	{
	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1)
	}
divergence_dates[S+1] <- 0
durations <- cbind(birth[sampled],death[sampled])
strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec))
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled])
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","Finds")
return(output)
}

simulate_extinction_date_given_variable_extinction_rates <- function(bday,mus,timescale,tbins,temp_prec)	{
b_bin <- sum(bday>=timescale);
nbins <- length(timescale);
d_bin <- b_bin-1;
ndy <- TRUE;
survived <- F;
while (ndy)	{
	d_bin <- d_bin+1;
	if (b_bin>=nbins)	{
		bint <- abs(timescale[d_bin])+tbins[d_bin]-bday;
		} else if (d_bin==nbins)	{
		bint <- tbins[d_bin];
		} else if (d_bin==b_bin)	{
		bint <- abs(timescale[d_bin+1]-bday);
		} else	{
		bint <- tbins[d_bin];
		}
#	print(paste(mus[d_bin],bint,ndy))
	if (runif(1)<=(1-dpois(0,mus[d_bin]*bint)))	ndy <- F;
	if (d_bin==nbins)	{
		ndy <- F;
		survived <- TRUE;
		}
	}
if (survived)	{
	dday <- timescale[d_bin]+tbins[d_bin];
	} else if (d_bin>b_bin)	{
	dday <- timescale[d_bin]+runif(1)*tbins[d_bin];
	} else	{
	dday <- bday+runif(1)*abs(timescale[d_bin+1]-bday);
	}
if (temp_prec>0)	dday <- round(dday/temp_prec,0)*temp_prec;
return(dday);
}

#tscale <- timescale;
simulate_events_given_varying_event_rates_over_time <- function(bday,dday,rates_over_time,tscale,tbins,temp_prec)	{
# bday: simulated origination time
# dday: simulated extinction time
#	rates_over_time: rates per chronostratigraphic unit
#	tscale: timescale providing onsets of chronostratigraphic units
# tbins: durations of chronostratigraphic units
# temp_prec: temporal precision of dates to be simulated (default: 0.1)
nbins <- length(tscale);
b_bin <- sum(bday>=tscale);
d_bin <- sum(dday>tscale);
if (tscale[nbins]!=0)	tscale[nbins+1] <- 0;
evdates <- c();
for (bn in b_bin:d_bin)	{
	bnn <- b_bin+bn-1;
	onset <- max(bday,tscale[bn]);
	if (bn==birth_bin[spc] && birth_bin[spc]==death_bin[spc])	{
		aspan <- abs(dday-bday);
		} else if (bn==b_bin)	{
		aspan <- abs(tscale[bn+1]-bday);
		} else if (bn==d_bin)	{
		aspan <- abs(dday-tscale[bn]);
		} else	{
		aspan <- tbins[bn];
		}
	bdaughters <- rpois(1,duration*rates_over_time[bn]);
#	bfossils <- rpois(1,duration*psis[bn]);
	if (bdaughters>0)	evdates <- c(evdates,sort(onset+runif(bdaughters)*aspan));
	}
if (length(evdates)>0 && temp_prec>0)	{
	evdates <- ceiling(evdates/temp_prec)*temp_prec;
	evdates[evdates<=bday] <- bday+temp_prec;
	}
return(evdates);
}

# editted 2024-08-22 for punc-eek project.
evolve_to_sampled_richness_S_shifting_rates <- function(S,ps,qs,psis,timescale,anagenesis=0,bifurcation=FALSE,temp_prec=0.01,debug=FALSE)	{
# cumulative and standing diversity
lambdas <- ps/(1+anagenesis);	# cladogenesis rate
mus <- qs/(1+anagenesis);			# extinction rate
darwin <- anagenesis*mean(mus);		# anagenesis rate;
nbins <- length(timescale);
if (timescale[1]>timescale[2])	timescale <- -timescale;
if (timescale[nbins]==0)	{
	tbins <- timescale[2:nbins]-timescale[1:(nbins-1)];
	last_bin <- tbins[nbins];
	nbins <- nbins-1;
	timescale <- timescale[1:nbins];
	} else	{
	tbins <- c(timescale[2:nbins],0)-timescale;
	last_bin <- abs(0-timescale[nbins]);
	names(tbins) <- names(timescale);
	}
timescale <- timescale-timescale[1];
bin_ends <- timescale+tbins;
the_end <- bin_ends[nbins];
birth <- timescale[1];
death[1] <- simulate_extinction_date_given_variable_extinction_rates(bday=timescale[1],mus,timescale,tbins,temp_prec=temp_prec);
while (death[1]<=0)
	death[1] <- simulate_extinction_date_given_variable_extinction_rates(bday=timescale[1],mus,timescale,tbins,temp_prec=temp_prec);
#print(death-birth)
duration <- death-birth;
birth_bin <- maxfinds <- dbr <- spc <- cumulative <- standing <- 1;
death_bin <- sum(death[1]>timescale);
ancestral <- c(0);
remaining <- c(1);
tries <- found <- 0	# this will keep track of how many species we have sampled
sampled <- c();	# list of sampled species
fa <- la <- occurrences <- vector(length=1);	# occurrences per species
all_finds <- matrix(0,1,maxfinds);
keep_going <- TRUE;
while (keep_going)	{
	spc <- remaining[1];
	duration <- abs((death[spc]-birth[spc]));
	mini_mes <- 0;
	if (anagenesis>0)	mini_mes <- rpois(1,darwin*duration);	# number of anagenetic pseudospeciations/pseudoextinctions
	# if anagenetic speciation happened, then breakup lineages into anagenetic events
	if (mini_mes>0)	{
		regens <- birth[spc]+sort(runif(mini_mes))*duration;  # place anagenetic events;
		new_ends <- c(regens,death[spc]);
		chuckies <- (cumulative+1):(cumulative+mini_mes);	# number anagenetics
		added <- rep(0,mini_mes);
		occurrences <- c(occurrences,added)
		ancestral <- c(ancestral,rep(spc,mini_mes));
		birth <- c(birth,regens);
		death[spc] <- regens[1];
		death_bin[spc] <- sum(death[spc]>timescale)
		death <- c(death,new_ends[2:length(new_ends)]);
		duration <- death[spc]-birth[spc];		# duration of morphospecies
		for (ch in 1:chuckies)	{
			birth_bin[chuckies[ch]] <- sum(birth[chuckies[ch]]>timescale);
			death_bin[chuckies[ch]] <- sum(death[chuckies[ch]]>timescale);
			}
		remaining <- c(remaining,cumulative+(1:mini_mes));	# add anagenetics to the queue
		cumulative <- cumulative+mini_mes;
		added <- rep(0,mini_mes)
		occurrences <- c(occurrences,added);
		fa <- c(fa,added)
		la <- c(la,added)
		all_finds <- rbind(all_finds,array(0,dim=c(mini_mes,maxfinds)));
#		rownames(all_finds) <- 1:cumulative;
		rownames(all_finds) <- 1:nrow(all_finds);
		}
	# divvy this up to allow for different origination rates at different points in species lifetime
	
	# first simulate sampling;
	fossils <- simulate_events_given_varying_event_rates_over_time(bday=birth[spc],dday=death[spc],rates_over_time=psis,tscale=timescale,tbins,temp_prec = temp_prec);
	occurrences[spc] <- nfinds <- length(fossils);
	if (nfinds>0)	{
		found <- found+1;
		if (maxfinds<nfinds)	{
			all_finds <- cbind(all_finds,array(0,dim=c(cumulative,nfinds-maxfinds)));
			maxfinds <- nfinds;
			}
		all_finds[spc,1:nfinds] <- fossils;
		fa[spc] <- min(fossils);
		la[spc] <- max(fossils);
		}
	
	# second, simulate cladogenesis;
	bdates <- simulate_events_given_varying_event_rates_over_time(bday=birth[spc],dday=death[spc],rates_over_time=lambdas,tscale=timescale,tbins,temp_prec = temp_prec);
	bdates <- bdates[bdates<the_end];
	daughters <- length(bdates);
	if (sum(bdates<birth[spc])>0)	{
		print("Call The Doctor!!!");
		break;
		}
	if (daughters > 0)	{
		if (bifurcation)	{
#			splits <- birth[spc]+sort(duration*runif(daughters))
#			splits <- temp_prec*ceiling(splits/temp_prec)
			# trick here: extract finds after splits & assign them to new ancestor;
			# trick here: new ancestor does not join "remaining"
			splits <- bdates;
			new_lines <- resid_lines <- c();
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1));
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)));
			prior_anc <- spc;
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1])
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1];
					}
				duration <- temp_prec*ceiling((-log(runif(1))/mus)/temp_prec)
				death <- c(death,(splits[i] + duration))
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (nfinds > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,matrix(0,daughters*2,maxfinds))
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						maxfinds <- occurrences[ap]
						all_finds <- cbind(all_finds,rep(0,length(occurrences)))
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap)
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f]
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f]
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
			birth <- c(birth,bdates);
			ancestral <- c(ancestral,rep(spc,daughters));
			ddates <- c();
			for (d in 1:daughters)	{
				ddates[d] <- simulate_extinction_date_given_variable_extinction_rates(bday=bdates[d],mus,timescale,tbins,temp_prec=temp_prec);
				death <- c(death,ddates[d]);
				death_bin <- c(death_bin,sum(ddates[d]>timescale));
				birth_bin <- c(birth_bin,sum(bdates[d]>timescale));
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,matrix(0,daughters,maxfinds));
			}
		cumulative <- max(remaining)
		rownames(all_finds) <- 1:cumulative;
		durations <- cbind(birth,death)
#		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	
	if (sum(all_finds[spc,]!=0)>0)	{
		fa[spc] <- min(all_finds[spc,all_finds[spc,]!=0]);
		la[spc] <- max(all_finds[spc,all_finds[spc,]!=0]);
		sampled <- unique(c(sampled,spc));
		}

	found <- length(sampled);
	if (debug)	print(found);
	staxa <- 1:length(fa);
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)];
		remaining <- remaining[order(birth[remaining])];
		clade_so_far <- rbind(staxa,ancestral,birth,death,fa,la);
		print(clade_so_far);
		}	else if (found < S)	{
		tries <- tries+1;
		dbr <- dbr+1
		# reboot..... :-(     
		spc <- cumulative <- standing <- 1
		birth <- timescale[1];
		birth_bin <- 1;
		death <- simulate_extinction_date_given_variable_extinction_rates(bday=timescale[1],mus,timescale,tbins,temp_prec = temp_prec);
		while (death<=0)	death <- simulate_extinction_date_given_variable_extinction_rates(bday=timescale[1],mus,timescale,tbins,temp_prec = temp_prec);
		death_bin <- sum(death[spc]>timescale);
		duration <- c(death[spc]-birth[spc])
		ancestral <- c(0);
		remaining <- c(1);
		# clear sampling vectors
		sampled <- c();	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1;
		all_finds <- matrix(0,1,maxfinds)
		}
#	cbind(order(birth[remaining]),birth[remaining])

	relv_fas <- sort(fa[fa>0]);
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

#plot(birth[fa>0][1:38],fa[fa>0])
# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral);
venn_tree <- venn_tree_info$Venn_Tree;
branchings <- venn_tree_info$Branchings;
branchings[S+1] <- 0;
vtree <- transform_venn_tree_to_vector_tree(venn_tree);
divergence_dates <- birth[venn_tree_info$First_Taxon];
if (sampled[1]==1)	{
	divergence_dates <- divergence_dates-divergence_dates[1]
	}	else	{
	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1)
	}
divergence_dates[S+1] <- 0
durations <- cbind(birth[sampled],death[sampled])
strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec))
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled])
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","Finds")
return(output)
}

# routine to evolve S taxa sampled from throughout the clade's history
# dbins <- rep(3,20);
# timescale <- c(cumsum(dbins[length(dbins):1])[length(dbins):1],0);
# lambdas <- c(2,sqrt(2),rep(1,18)); mus <- rep(1,20); phis <- rep(1,20);
evolve_to_sampled_richness_S_shifting_rates_old <- function(S,lambdas,mus,phis,timescale,bifurcation=F,temp_prec=0.1)	{
# cumulative and standing diversity
if (timescale[1]>timescale[2])	timescale <- -timescale;
dbins <- abs(timescale[2:length(timescale)]-timescale[1:(length(timescale)-1)]);
nbins <- length(dbins);
maxfinds <- dbr <- spc <- cumulative <- standing <- 1;
birth <- timescale[1]+c(temp_prec);
duration <- temp_prec*ceiling((-log(runif(1))/mus[1])/temp_prec);
death <- c(birth[spc]+duration);
birth_bin <- 1;
death_bin <- sum(death>timescale);
ancestral <- c(0);
remaining <- c(1);
found <- 0	# this will keep track of how many species we have sampled
sampled <- c();	# list of sampled species
fa <- la <- occurrences <- vector(length=1);	# occurrences per species
all_finds <- matrix(0,1,maxfinds);
keep_going <- TRUE;
while (keep_going)	{
	spc <- remaining[1]
	duration <- abs((death[spc]-birth[spc]));
	daughters <- rpois(1,(lambdas*duration));
	fossils <- rpois(1,phis*duration);
	if (fossils > 0)	finds <-  sort(birth[spc]+runif(fossils)*duration)
	if (daughters > 0)	{
		if (bifurcation)	{
			splits <- birth[spc]+sort(duration*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				duration <- temp_prec*ceiling((-log(runif(1))/mus)/temp_prec)
				death <- c(death,(splits[i] + duration))
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (fossils > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,matrix(0,daughters*2,maxfinds))
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						maxfinds <- occurrences[ap]
						all_finds <- cbind(all_finds,rep(0,length(occurrences)))
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap)
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f]
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f]
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
#			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			birth <- c(birth,birth[spc]+sort(temp_prec*ceiling((duration*runif(daughters))/temp_prec)));
			ancestral <- c(ancestral,rep(spc,daughters));
			for (d in 1:daughters)	{
				duration <- temp_prec*ceiling((-log(runif(1))/mus)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+duration
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,matrix(0,daughters,maxfinds))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if ((daughters == 0 || !bifurcation) && fossils > 0)	{
		occurrences[spc] <- fossils
		fa[spc] <- min(finds)
		la[spc] <- max(finds)
		sampled <- c(sampled,spc)
		if (maxfinds < occurrences[spc])	{
			added <- occurrences[spc] - maxfinds
			maxfinds <- occurrences[spc]
			all_finds <- cbind(all_finds,matrix(0,length(occurrences),added))
			all_finds[spc,] <- finds
			}	else	{
			all_finds[spc,1:fossils] <- finds	
			}
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		duration <- temp_prec*ceiling((-log(runif(1))/mus)/temp_prec)
		death <- c(birth[spc]+duration)
		ancestral <- c(0)
		remaining <- c(1)
		# clear sampling vectors
		sampled <- c()	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1
		all_finds <- matrix(0,1,maxfinds)
		}
	
	found <- length(sampled)
	relv_fas <- sort(fa[sampled])
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral)
venn_tree <- venn_tree_info$Venn_Tree
branchings <- venn_tree_info$Branchings
branchings[S+1] <- 0
vtree <- transform_venn_tree_to_vector_tree(venn_tree)
divergence_dates <- birth[venn_tree_info$First_Taxon]
if (sampled[1]==1)	{
	divergence_dates <- divergence_dates-divergence_dates[1]
	}	else	{
	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1)
	}
divergence_dates[S+1] <- 0
durations <- cbind(birth[sampled],death[sampled])
strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec))
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled])
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","Finds")
return(output)
}

# routine to evolve S contemporaneous (e.g., extant) taxa
MBL_style_to_standing_richness_S <- function(S,lambda,mu,bifurcation=F,temp_prec=0.1)	{
# cumulative and standing diversity
spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
dbr <- 1
while (max(standing) < S)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	if (daughters > 0)	{
		if (bifurcation)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			remaining <- c(remaining,cumulative+(1:daughters))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		}
	}
output <- cbind(ancestral,birth,death)
#output <- list(vtree,venn_tree,branchings,divergence_dates)
#names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates")
return(output)
}

logistic_to_standing_richness_S_plus <- function(R=2,mu=1.0,K=25,bifurcation=F,temp_prec = 0.1,post_exponential=1,max_life=MAXNO)	{
# 2022-09-06: Written to provide post-peak history for punctuated tip-dating talk
spc <- cumulative <- standing <- 1;
birth <- 0;
life <- -log(runif(1))/mu;
while (life>max_life) life <- -log(runif(1))/mu;
death <- c(birth[spc]+life);
ancestral <- c(0);
remaining <- c(1);
dbr <- 1;
bin <- 1;	# general stage number
started <- 0;
elapsed <- temp_prec;
origination_rates <- c();
while (max(standing) < K)	{
	lambda <- immediate_expected_origination_given_logistic_constant_extinction(mu,R,K,S=standing);
	origination_rates <- c(origination_rates,lambda);
	for (s in 1:standing)	{
		spc <- remaining[s];
		base <- max(started,birth[spc]);
		top <- min(elapsed,death[spc]);
		daughters <- rpois(1,(lambda*abs(top-base)));
#		if (birth[spc]<=started & death[spc]>elapsed)	{
#			# species ranges through interval;
#		  daughters <- rpois(1,(lambda*temp_prec))
#			} else if (birth[spc]<=started & death[spc]<=elapsed)	{
#			daughters <- rpois(1,(lambda*(death[spc]-started)));
#			}	# if taxon appears partway through time slice, then reduce chance of daughters
#			} else if (birth[spc]>=started & death[spc]>elapsed)	{
#			daughters <- rpois(1,(lambda*(elapsed-birth[spc])))
#			}	# if taxon appears partway through time slice, then reduce chance of daughters
		if (daughters > 0)	{
			if (bifurcation)	{
				splits <- started+runif(1)*temp_prec
				new_lines <- resid_lines <- c()
				resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
				new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
				prior_anc <- spc
				for (i in 1:daughters)	{
					rb <- resid_lines[i]
					lb <- new_lines[i]
					ancestral <- c(ancestral,prior_anc,prior_anc)
					birth <- c(birth,splits[i],splits[i])
					if (i < daughters)  {
						death <- c(death,splits[i+1]-temp_prec)
						}	else	{
						death <- c(death,death[spc])
						death[spc] <- splits[1]-temp_prec
						}
					life <- -log(runif(1))/mu;
					while (life>max_life) life <- -log(runif(1))/mu;
					death <- c(death,(splits[i] + life))
					prior_anc <- rb	# the right branch leads to the next split
					}
				remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
				xx <- max(started,birth[spc])
				birth <- c(birth,xx+sort(temp_prec*runif(daughters)))
				ancestral <- c(ancestral,rep(spc,daughters))
				for (d in 1:daughters)	{
					life <- -log(runif(1))/mu;
					while (life>max_life) life <- -log(runif(1))/mu;
					death[d+cumulative] <- birth[d+cumulative]+life;
					}
				remaining <- c(remaining,cumulative+(1:daughters))
				}
			cumulative <- max(remaining)
			ranges <- cbind(birth,death)
			standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
			}
		}
	remaining <- remaining[death[remaining]>elapsed]
	if (length(remaining) == 0)	{
#		remaining <- remaining[2:length(remaining)]
#		}	else if (max(standing) < S)	{
		dbr <- dbr+1	# trials needed.
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- 0;
		life <- -log(runif(1))/mu;
		while (life>max_life) life <- -log(runif(1))/mu;
		death <- c(birth[spc]+life);
		ancestral <- 0;
		remaining <- 1;
		started <- 0;
		elapsed <- temp_prec;
		origination_rates <- c();
		} # end reboot
	standing <- length(remaining);
	started <- started+temp_prec;
	elapsed <- elapsed+temp_prec;
#	print(c(standing,"•",remaining));
	}
extra_time <- elapsed+(post_exponential*elapsed);
while (elapsed < extra_time)  {
  lambda <- immediate_expected_origination_given_logistic_constant_extinction(mu,R,K,S=standing);
  origination_rates <- c(origination_rates,lambda);
  for (s in 1:standing)	{
    spc <- remaining[s];
    base <- max(started,birth[spc]);
    top <- min(elapsed,death[spc]);
    daughters <- rpois(1,(lambda*abs(top-base)));
    #		if (birth[spc]<=started & death[spc]>elapsed)	{
    #			# species ranges through interval;
    #		  daughters <- rpois(1,(lambda*temp_prec))
    #			} else if (birth[spc]<=started & death[spc]<=elapsed)	{
    #			daughters <- rpois(1,(lambda*(death[spc]-started)));
    #			}	# if taxon appears partway through time slice, then reduce chance of daughters
    #			} else if (birth[spc]>=started & death[spc]>elapsed)	{
    #			daughters <- rpois(1,(lambda*(elapsed-birth[spc])))
    #			}	# if taxon appears partway through time slice, then reduce chance of daughters
    if (daughters > 0)	{
      if (bifurcation)	{
        splits <- started+runif(1)*temp_prec
        new_lines <- resid_lines <- c()
        resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
        new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
        prior_anc <- spc
        for (i in 1:daughters)	{
          rb <- resid_lines[i]
          lb <- new_lines[i]
          ancestral <- c(ancestral,prior_anc,prior_anc)
          birth <- c(birth,splits[i],splits[i])
          if (i < daughters)  {
            death <- c(death,splits[i+1]-temp_prec)
            }	else	{
            death <- c(death,death[spc])
            death[spc] <- splits[1]-temp_prec
            }
          life <- -log(runif(1))/mu;
          while (life>max_life) life <- -log(runif(1))/mu;
          death <- c(death,(splits[i] + life))
          prior_anc <- rb	# the right branch leads to the next split
        }
        remaining <- c(remaining,new_lines)
      }	else	{
        # for budding model, just get appearance and death of each daughter
        xx <- max(started,birth[spc])
        birth <- c(birth,xx+sort(temp_prec*runif(daughters)))
        ancestral <- c(ancestral,rep(spc,daughters))
        for (d in 1:daughters)	{
          life <- -log(runif(1))/mu;
          while (life>max_life) life <- -log(runif(1))/mu;
          death[d+cumulative] <- birth[d+cumulative]+life;
        }
        remaining <- c(remaining,cumulative+(1:daughters))
      }
      cumulative <- max(remaining)
      ranges <- cbind(birth,death)
      standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
    }
  }
  remaining <- remaining[death[remaining]>elapsed]
  standing <- length(remaining);
  started <- started+temp_prec;
  elapsed <- elapsed+temp_prec;
#	print(c(standing,"•",remaining));
#	print(c(standing,"•",remaining));
  }

history <- data.frame(ancestral=as.numeric(ancestral),birth=as.numeric(birth),death=as.numeric(death));
totu <- nrow(history);
ancestors <- sort(unique(history$ancestral[history$ancestral>0]));
nanc <- length(ancestors);
venn_tree <- array(0,dim=c(nanc,totu));
all_spc <- 1:totu;
for (i in 1:nanc)	{
	this_node <- c(ancestors[i],all_spc[history$ancestral==ancestors[i]]);
	grandparents <- which(venn_tree==ancestors[i],arr.ind=TRUE)[,1]
	venn_tree[i,1:length(this_node)] <- this_node;
	these_desc <- this_node[!this_node %in% ancestors[i]]
	ndesc <- length(these_desc);
	gp <- 0;
	while (gp<length(grandparents))	{
		gp <- gp+1;
		pdesc <- sum(venn_tree[grandparents[gp],]>0);
		venn_tree[grandparents[gp],(pdesc+1):(pdesc+ndesc)] <- these_desc;
		}
	}
output <- list(history,venn_tree,origination_rates);
names(output) <- c("history","venn_tree","lambdas");
return(output)
}

logistic_to_standing_richness_S <- function(R=2,mu=1.0,K=25,bifurcation=F,temp_prec = 0.1,max_life=MAXNO)	{
spc <- cumulative <- standing <- 1;
birth <- c(temp_prec);
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec);
while (life>max_life) life <- -log(runif(1))/mu;
death <- c(birth[spc]+life);
ancestral <- c(0);
remaining <- c(1);
dbr <- 1;
bin <- 1;	# general stage number
started <- temp_prec;
elapsed <- 2*temp_prec;
origination_rates <- c();
while (max(standing) < K)	{
	lambda <- immediate_expected_origination_given_logistic_constant_extinction(mu,R,K,S=standing);
	origination_rates <- c(origination_rates,lambda);
	for (s in 1:standing)	{
		spc <- remaining[s]
		if (birth[spc]<=started)	{
			daughters <- rpois(1,(lambda*temp_prec))
			} else	{
			daughters <- rpois(1,(lambda*(elapsed-birth[spc])))
			}	# if taxon appears partway through time slice, then reduce chance of daughters
		if (daughters > 0)	{
			if (bifurcation)	{
				splits <- started+runif(1)*temp_prec
				new_lines <- resid_lines <- c()
				resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
				new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
				prior_anc <- spc
				for (i in 1:daughters)	{
					rb <- resid_lines[i]
					lb <- new_lines[i]
					ancestral <- c(ancestral,prior_anc,prior_anc)
					birth <- c(birth,splits[i],splits[i])
					if (i < daughters)  {
						death <- c(death,splits[i+1]-temp_prec)
						}	else	{
						death <- c(death,death[spc])
						death[spc] <- splits[1]-temp_prec
						}
					life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
					while (life>max_life) life <- -log(runif(1))/mu;
					death <- c(death,(splits[i] + life))
					prior_anc <- rb	# the right branch leads to the next split
					}
				remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
				xx <- max(started,birth[spc])
				birth <- c(birth,xx+sort(temp_prec*runif(daughters)))
				ancestral <- c(ancestral,rep(spc,daughters))
				for (d in 1:daughters)	{
					life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
					while (life>max_life) life <- -log(runif(1))/mu;
					death[d+cumulative] <- birth[d+cumulative]+life
					}
				remaining <- c(remaining,cumulative+(1:daughters))
				}
			cumulative <- max(remaining)
			ranges <- cbind(birth,death)
			standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
			}
		}
	remaining <- remaining[death[remaining]>elapsed]
	if (length(remaining) == 0)	{
#		remaining <- remaining[2:length(remaining)]
#		}	else if (max(standing) < S)	{
		dbr <- dbr+1	# trials needed.
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		while (life>max_life) life <- -log(runif(1))/mu;
		death <- c(birth[spc]+life)
		ancestral <- 0
		remaining <- 1
		started <- 0
		elapsed <- temp_prec
		origination_rates <- c();
		} # end reboot
	standing <- length(remaining);
	started <- started+temp_prec;
	elapsed <- elapsed+temp_prec;
#	print(c(standing,"•",remaining));
	}
history <- data.frame(ancestral=as.numeric(ancestral),birth=as.numeric(birth),death=as.numeric(death));
totu <- nrow(history);
ancestors <- sort(unique(history$ancestral[history$ancestral>0]));
nanc <- length(ancestors);
venn_tree <- array(0,dim=c(nanc,totu));
all_spc <- 1:totu;
for (i in 1:nanc)	{
	this_node <- c(ancestors[i],all_spc[history$ancestral==ancestors[i]]);
	grandparents <- which(venn_tree==ancestors[i],arr.ind=TRUE)[,1]
	venn_tree[i,1:length(this_node)] <- this_node;
	these_desc <- this_node[!this_node %in% ancestors[i]]
	ndesc <- length(these_desc);
	gp <- 0;
	while (gp<length(grandparents))	{
		gp <- gp+1;
		pdesc <- sum(venn_tree[grandparents[gp],]>0);
		venn_tree[grandparents[gp],(pdesc+1):(pdesc+ndesc)] <- these_desc;
		}
	}
output <- list(history,venn_tree,origination_rates);
names(output) <- c("history","venn_tree","lambdas");
return(output)
}

accersi_simulated_divergence_times <- function(sampled,ancestral,birth)	{
# first, get the venn tree for all simulated taxa leading to sampled
raw_venn_tree <- accersi_raw_venn_tree_from_ancestor_list(sampled,ancestral=ancestral);
notu <- length(raw_venn_tree[1,]);
ttl_clades <- rnodes <- nrow(raw_venn_tree);
raw_ancestors <- as.numeric(rownames(raw_venn_tree));
ttl_prog <- vector(length=ttl_clades);
for (i in 1:ttl_clades)	ttl_prog[i] <- sum(raw_venn_tree[i,]>0);
singletons <- (1:rnodes)[ttl_prog==1];
nsingletons <- length(singletons);
i <- nsingletons;
otu_divergences <- data.frame(onset=as.numeric(birth[sampled]),end=as.numeric(birth[sampled]));
while (i>0)	{
	node_no <- singletons[i]
	spc <- raw_venn_tree[node_no,1];
	otu <- match(spc,sampled);
	all_nodes <- sort(which(raw_venn_tree==spc,arr.ind = TRUE)[,1],decreasing=TRUE);
	lo_node <- min(all_nodes[ttl_prog[all_nodes]==1]);
	lo_anc <- as.numeric(rownames(raw_venn_tree)[lo_node]);
	otu_divergences$onset[otu] <- birth[lo_anc];
	i <- i-1;
	}
rownames(otu_divergences) <- sampled;

reduction <- remove_singletons_from_venn_tree(venn_tree=raw_venn_tree);
condensed_raw_venn_tree <- reduction$Condensed_Venn_Tree;
cnodes <- nrow(condensed_raw_venn_tree);
final_raw_venn_tree <- unique(condensed_raw_venn_tree);
node_divergences <- data.frame(onset=as.numeric(birth[as.numeric(rownames(final_raw_venn_tree))]),end=as.numeric(birth[as.numeric(rownames(final_raw_venn_tree))]));
nnodes <- nrow(final_raw_venn_tree);
lost_nodes <- sort((1:cnodes)[!rownames(condensed_raw_venn_tree) %in% rownames(final_raw_venn_tree)],decreasing=TRUE);
ln <- length(lost_nodes);
i <- 0;
while (i<ln)	{
	i <- i+1;
	frvt <- row.match(condensed_raw_venn_tree[lost_nodes[i],],final_raw_venn_tree);
	node_divergences$end[frvt] <- birth[as.numeric(rownames(condensed_raw_venn_tree)[lost_nodes[i]])];
	}
rownames(node_divergences) <- rownames(final_raw_venn_tree);
divergence_dates <- rbind(otu_divergences,node_divergences);
return(divergence_dates);
}

# routine to convert simulation into phylogeny
accersi_simulated_venn_tree_info <- function(sampled,ancestral,true_births=NULL)	{
# first, get the venn tree for all simulated taxa leading to sampled
raw_venn_tree <- accersi_raw_venn_tree_from_ancestor_list(sampled,ancestral=ancestral)
#write.table(raw_venn_tree,"raw_tree.xls",sep="\t")
# now, reduce it to unique nodes, adding redundant nodes to branch lengths
#	for species, this will be venn_tree nodes in which just that specees appears
#	for nodes, this will duplicate vectors
if (is.null(true_births))	{}
otus <- length(raw_venn_tree[1,]);
reduction <- remove_singletons_from_venn_tree(venn_tree=raw_venn_tree)
#write.table(reduction$Condensed_Venn_Tree,"raw_clades_only.xls",sep="\t")
condensed_raw_venn_tree <- unique(reduction$Condensed_Venn_Tree)
#write.table(condensed_raw_venn_tree,"reduced_raw_tree.xls",sep="\t")
#condensed_raw_venn_tree <- reduction$Condensed_Venn_Tree
unique_nodes <- dim(condensed_raw_venn_tree)[1]
branchings <- reduction$Branchings
first_taxon <- raw_venn_tree[1,]
node_ancestors <- as.numeric(rownames(condensed_raw_venn_tree))
for (n in 1:otus)	{
	u <- 1
	while (u < branchings[n])	{
		first_taxon[n] <- ancestral[first_taxon[n]]
		u <- u+1
		}
	}
first_taxon <- c(first_taxon,node_ancestors)
branchings <- c(branchings,rep(0,unique_nodes))
# now, take duplicate nodes and reduce them, adding to their branch lengths
for (u in 2:unique_nodes)	{
	htu <- otus + u;
	row.is.a.match <- apply(raw_venn_tree, 1, identical, condensed_raw_venn_tree[u,]) 
	match.idx <- sort(which(row.is.a.match));
	branchings[htu] <- sum(row.is.a.match);
	}
#	first_taxon <- c(first_taxon,all_ancestors[match.idx[1]])
#	if (is.na(match(node_ancestors[htu],all_ancestors)))	{
#		first_taxon[htu] <- all_ancestors[match.idx[1]]
#		}	else	{
#		first_taxon[htu] <- condensed_raw_venn_tree[u,1]
#		branchings[condensed_raw_venn_tree[u,1]] <- 0
#		}
#	u <- u+1
#	first_taxon
#	branchings
	# add something here to track down original ancestor using all_ancestors vector
#	which(condensed_raw_venn_tree[u,1]==raw_venn_tree[,1],arr.ind=TRUE)
#	i<-1:unique_nodes
#	which(condensed_raw_venn_tree[u,],raw_venn_tree[i,])
#	}

# there might be ancestral species sampled: if so, then they have
#	branch durations of zero
all_ancestors <- accersi_all_ancestors_for_sampled_taxa(sampled,ancestral)
sampled_ancestors <- intersect(sampled,all_ancestors);
s_anc <- length(sampled_ancestors);
s <- 0;
while (s < s_anc)	{
	s <- s+1;
	htu <- otus+match(sampled_ancestors[s],condensed_raw_venn_tree[,1]);
	}
#if (!is.na(match(otus,all_ancestors)))	{
## PICK UP HERE!!!!

for (s in 1:otus)	{
	if (!is.na(match(sampled[s],all_ancestors)))	{
		sampled_ancestors <- c(sampled_ancestors,sampled[s])
		htu <- otus+match(sampled[s],condensed_raw_venn_tree[,1])
		### ancestral branches already have at least one
		###		add any additional branches
		branchings[htu] <- branchings[htu] + (branchings[s] - 1)
		branchings[s] <- 0
		}
	}

# now, change all of the sampled taxon numbers to 1…S
venn_tree <- match(condensed_raw_venn_tree[1,],raw_venn_tree[1,])
for (n in 2:unique_nodes)	{
	xx <- match(condensed_raw_venn_tree[n,],raw_venn_tree[1,])
	xx[is.na(xx)] <- 0
	venn_tree <- rbind(venn_tree,xx)
	}
# name rows after simulated founder
#rownames(venn_tree) <- first_taxon[(otus+1):length(first_taxon)]
rownames(venn_tree) <- rownames(condensed_raw_venn_tree)

output <- list(venn_tree,branchings,first_taxon)
names(output) <- c("Venn_Tree","Branchings","First_Taxon")
return(output)
}

# routine to eliminate nodes with only one sampled descendant
remove_singletons_from_venn_tree <- function(venn_tree)	{
sampled <- sort(venn_tree[1,])
init_nodes <- dim(venn_tree)[1]
otus <- length(sampled)
branchings <- rep(1,otus)
fs <- c()
for (n in 1:init_nodes)
	if (sum(venn_tree[n,]>0)==1)
		fs <- c(fs,n)
# add to branch lengths of species with singleton representations
for (sn in 1:length(fs))	{
	n <- fs[sn]
	spc <- match(venn_tree[n,1],sampled)
	branchings[spc] <- branchings[spc]+1
	}
# list retained nodes
ret_nodes <- (1:init_nodes)[!(1:init_nodes) %in% fs]
red_venn_tree <- venn_tree[ret_nodes,]
output <- list(red_venn_tree,branchings)
names(output) <- c("Condensed_Venn_Tree","Branchings")
return(output)
}

#mine <- 1:6 
#table.combos <- matrix(data = 1:12, nrow = 10, ncol = 6, byrow=TRUE) 
#row.is.a.match <- apply(table.combos, 1, identical, mine) 
#match.idx <- which(row.is.a.match) 
#total.matches <- sum(row.is.a.match) 

# routine to get basic venn tree from initial simulation
accersi_raw_venn_tree_from_ancestor_list <- function(sampled,ancestral)	{
notu <- length(sampled)
all_ancestors <- accersi_all_ancestors_for_sampled_taxa(sampled,ancestral)
all_ancestors_per_otu <- list_all_ancestors_for_all_sampled_taxa(sampled,ancestral)
raw_venn_tree <- matrix(0,length(all_ancestors),notu)
all_desc <-rep(0,length(all_ancestors))
for (i in 1:notu)	{
	spc <- sampled[i]
	nodes <- match(all_ancestors_per_otu[[i]],all_ancestors)
	if (!is.na(match(spc,all_ancestors)))	nodes <- c(nodes,match(spc,all_ancestors))
	all_desc[nodes] <- all_desc[nodes]+1
	for (n in 1:length(nodes))	{
#		print(c(n,nodes[n],all_desc[nodes[n]]))
		raw_venn_tree[nodes[n],all_desc[nodes[n]]] <- spc
		}
#	raw_venn_tree[,1:3]
#	i <- i+1
	}
rNodes <- dim(raw_venn_tree)[1]
rownames(raw_venn_tree) <- all_ancestors
while (identical(raw_venn_tree[1,],raw_venn_tree[2,]))	{
	raw_venn_tree <- raw_venn_tree[2:rNodes,]
	rNodes <- rNodes - 1
	}
return(raw_venn_tree)	# it is fine at this point: see what happens to it....
}

# routine to list all simulated ancestors for "sampled" simulated taxa
accersi_all_ancestors_for_sampled_taxa <- function(sampled,ancestral)	{
otu <- length(sampled)
added <- c()
for (i in 1:otu)	{
	spc <- sampled[i]
	added <- sort(unique(c(added,accersi_all_ancestors_for_taxon(spc,ancestral))),decreasing=FALSE)
	}
return(added)
}

# routine to list all simulated ancestors for each "sampled" simulated taxa
list_all_ancestors_for_all_sampled_taxa <- function(sampled,ancestral)	{
otu <- length(sampled)
all_anc <- list()
for (i in 1:otu)	{
	spc <- sampled[i]
	all_anc[[i]] <- sort(accersi_all_ancestors_for_taxon(spc,ancestral),decreasing=FALSE)
	}
return(all_anc )
}

# routine to list all simulated ancestors for a "sampled" simulated taxa
accersi_all_ancestors_for_taxon <- function(spc,ancestral)	{
anc <- ancestral[spc]
while (min(anc)>=1)	{
	anc <- c(anc,ancestral[min(anc)])
	}
return(anc)
}

# ROUTINES RELATED INVERSE MODELLING WITH COMPATIBILITY ####
# simulate character evolution up to N steps & tally compatibility
#		2017-02-22: use beta distribution to "smooth" P[compatibility|steps]!!!
#		2017-08-10: make sure initial compatibility is maximum or close to it.....
evolve_compatibility_over_N_changes <- function(N,init_chmatrix,venn_tree,branchings,states,types,hidden_reversals=TRUE,UNKNOWN,INAP,repl=1,print_freq=10)	{
# N: maximum number of steps
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# hidden_reversals: if true, then character can change 2+ times per sampled branch
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- ncol(venn_tree);	# number of observed taxa in tree
nodes <- nrow(venn_tree);
#simchmatrix <- matrix(0,notu,nchars)
nchars <- ncol(init_chmatrix);
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10
# get those branches where there can be change (sampled ancestors are zero, with nodal branch ≥ 1)
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times

# each branch gets one additional representative per unsampled ancestor.
# 		With perfect sampling, all taxa are entered once.
ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- sum(branchings)
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb);	# branch changes per character
st_rich <- matrix(0,nchars,max(states));	# richness of each state

mx_ch <- count_scored_otu_per_character(chmatrix=simchmatrix)

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
# for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes
for (ch in 1:nchars)  {
	# scramble order in which branches are sampled
	if (hidden_reversals==TRUE)	{
		x <- permute(subset(pabm[ch,],pabm[ch,]>0))
		}	else	{
		x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
		}
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	}
	# this routine should be unnecessary, and it crashes sometimes.
#	if (length(x) < tb)	{
#	if (pcc[ch] < tb)	{
#		y <- 1+pcc[ch]
#		pabm[ch,y:tb] <- 0	# problem here sometimes...
#		pabm[ch,(pcc[ch]+1):tb] <- 0	# problem here sometimes...
#		}
#	if (pcc[ch]<length(ab)	for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
#	}

# first, make sure that all states appear
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=N)
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- count_states(simchmatrix);
#for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# makes sure that all branches have 1+ change
# branch_changes[unique(ab)]
#  PUT THIS ROUTINE IN OTHER MODULES!
if (sum(branch_changes[unique_ab]==0)>0)	{
	needy <- unique(ab)[branch_changes[unique(ab)]==0]
	nb <- length(needy)
	for (b in 1:nb)	{
		# get first branch in need of a change
		br <- needy[b]
		# target characters where it would change earliest
		candidates <- which(pabm==br,arr.ind=TRUE)
		cn <- 1
		ch <- candidates[cn,1]
		while (states[ch]<2
			   || char_changes[ch]>=mx_ch[ch]
			   || char_changes[ch]>=candidates[cn,2]
			   || is.na(match(br,pabm[ch,])))	{
			cn <- cn+1
			ch <- candidates[cn,1]
			}	# make sure that this is an appropriate character
		# determine shifts
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		# routine for tips
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- 1
		char_changes[ch] <- char_changes[ch]+1
		# now, flip over branches
		pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
		pabm[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}
#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1;
simcompat[delta] <- sum(character_compats)/2;
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta]);
#if (repl>0 && print_status)	print(c(repl,delta,simcompat[delta]))

if (print_freq<0)	{
	print_counters <- seq(delta,N,(N-delta)*print_freq)
	} else	{
	print_counters <- seq(print_freq,N,by=print_freq)[seq(print_freq,N,by=print_freq)>delta];
	}
counter <- 0; p_c <- length(print_counters);
for (d in (delta+1):N)	{
	counter <- counter+1;
	# add the br != 0 check!!!
	br <- 0;
	while (br==0)	{
		ch <- ceiling(nchars*runif(1));
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		}
	
	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
#	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP);
	new_compat <- compatibility_of_a_character(c1=ch,chmatrix=simchmatrix,states=states,types=types);
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	dcompat <- sum(new_compat - prior_compat)
	simcompat[d] <- simcompat[d-1]+dcompat
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
#	if (repl>0 && counter%%print_freq==0)	print(c(repl,d,simcompat[d]))
	if (d %in% print_counters)	{
		pp <- match(d,print_counters);
		if (pp==1)	{
			paste(100*round(100*pp/p_c,0)/100,"% done",sep="");
			} else if (last_pp<10)	{
			cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			flush.console();
			print(paste(100*round(100*pp/p_c,0)/100,"% done",sep=""));
			} else	{
			cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			flush.console();
			print(paste(100*round(100*pp/p_c,0)/100,"% done",sep=""));
			}
		last_pp <- round(100*pp/p_c,0);
		}
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	if (sum(prior_compat==new_compat)!=nchars)	{
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}	# modify compatibility matrix of cha}racters affected by change
		}
	}

steps_v_compat <- data.frame(steps=as.numeric(steps_v_compat[,1]),compat=as.numeric(steps_v_compat[,2]));
return (steps_v_compat);
}

evolve_compatibility_over_N_changes_old <- function(N,init_chmatrix,venn_tree,branchings,states,types,hidden_reversals=TRUE,UNKNOWN,INAP,repl=1,print_status=TRUE,print_bar=TRUE,print_freq=10)	{
# N: maximum number of steps
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# hidden_reversals: if true, then character can change 2+ times per sampled branch
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- ncol(venn_tree);	# number of observed taxa in tree
nodes <- nrow(venn_tree);
#simchmatrix <- matrix(0,notu,nchars)
nchars <- ncol(init_chmatrix);
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10
# get those branches where there can be change (sampled ancestors are zero, with nodal branch ≥ 1)
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times

# each branch gets one additional representative per unsampled ancestor.
# 		With perfect sampling, all taxa are entered once.
ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
#	if (b==1)	{
#		# ab: availabe branches
#		ab <- rep(br,branchings[br])
#		}	else	{
#		ab <- c(ab,rep(br,branchings[br]))
#		}
	}
tb <- sum(branchings)
#tb <- length(ab)			# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))	# richness of each state

mx_ch <- count_scored_otu_per_character(chmatrix=simchmatrix)

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
# for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes
for (ch in 1:nchars)  {
	# scramble order in which branches are sampled
	if (hidden_reversals==TRUE)	{
		x <- permute(subset(pabm[ch,],pabm[ch,]>0))
		}	else	{
		x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
		}
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	}
	# this routine should be unnecessary, and it crashes sometimes.
#	if (length(x) < tb)	{
#	if (pcc[ch] < tb)	{
#		y <- 1+pcc[ch]
#		pabm[ch,y:tb] <- 0	# problem here sometimes...
#		pabm[ch,(pcc[ch]+1):tb] <- 0	# problem here sometimes...
#		}
#	if (pcc[ch]<length(ab)	for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
#	}

# first, make sure that all states appear
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=N)
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- count_states(simchmatrix);
#for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# makes sure that all branches have 1+ change
# branch_changes[unique(ab)]
#  PUT THIS ROUTINE IN OTHER MODULES!
if (sum(branch_changes[unique_ab]==0)>0)	{
	needy <- unique(ab)[branch_changes[unique(ab)]==0]
	nb <- length(needy)
	for (b in 1:nb)	{
		# get first branch in need of a change
		br <- needy[b]
		# target characters where it would change earliest
		candidates <- which(pabm==br,arr.ind=TRUE)
		cn <- 1
		ch <- candidates[cn,1]
		while (states[ch]<2
			   || char_changes[ch]>=mx_ch[ch]
			   || char_changes[ch]>=candidates[cn,2]
			   || is.na(match(br,pabm[ch,])))	{
			cn <- cn+1
			ch <- candidates[cn,1]
			}	# make sure that this is an appropriate character
		# determine shifts
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		# routine for tips
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- 1
		char_changes[ch] <- char_changes[ch]+1
		# now, flip over branches
		pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
		pabm[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}
#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1;
simcompat[delta] <- sum(character_compats)/2;
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta]);
if (repl>0 && print_status && !print_bar)	print(c(repl,delta,simcompat[delta]))

if (print_freq<0)	{
	print_counters <- seq(delta,N,(N-delta)*print_freq)
	} else	{
	print_counters <- seq(print_freq,N,by=print_freq)[seq(print_freq,N,by=print_freq)>delta];
	}
counter <- 0;
for (d in (delta+1):N)	{
	counter <- counter+1
	# add the br != 0 check!!!
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1));
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		}
	
	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)

	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	dcompat <- sum(new_compat - prior_compat)
	simcompat[d] <- simcompat[d-1]+dcompat
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (repl>0 && counter%%print_freq==0)	print(c(repl,d,simcompat[d]))
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	if (sum(prior_compat==new_compat)!=nchars)	{
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}	# modify compatibility matrix of cha}racters affected by change
		}
	}
#output <- list(d,simcompmat)
rownames(steps_v_compat) <- rep("",dim(steps_v_compat)[1])
colnames(steps_v_compat) <- c("Steps","Compat")
return (steps_v_compat)
}

# simulate character evolution up to a certain compatibility
# returns # steps required.  This will be a fraction if compatibility is passed
#evolve_to_particular_compatibilty <- function(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=-1)	{
evolve_to_particular_compatibilty_old <- function(init_chmatrix,ttl_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=-1)	{
# init_chmatrix: initial character matrix to be emulated
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
nchars <- ncol(init_chmatrix)
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10

unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)

ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	# this causes problems sometimes and should not be necessary
#	if (pcc[ch] < tb)	pabm[ch,(pcc[ch]+1):tb] <- 0
	}

# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# makes sure that all branches have 1+ change
#  PUT THIS ROUTINE IN OTHER MODULES!
needy <- unique(ab)[branch_changes[unique(ab)]==0]
nb <- length(needy)
b <- 1
while (b < nb)	{
	# get first branch in need of a change
	br <- needy[b]
	# target characters where it would change earliest
	candidates <- which(pabm==br,arr.ind=TRUE)
	cn <- 1
	ch <- candidates[cn,1]
	while (states[ch]<2
		   || char_changes[ch]>=mx_ch[ch]
		   || char_changes[ch]>=candidates[cn,2]
		   || is.na(match(br,pabm[ch,])))	{
		cn <- cn+1
		ch <- candidates[cn,1]
		}	# make sure that this is an appropriate character
	# determine shifts
	if (states[ch]==2)	{
		dstates <- c(1,0)
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	
		} 
	# routine for tips
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			cc <- simchmatrix[br,ch]
			simchmatrix[br,ch] <- dstates[cc+1]
			} else {
			ndg <- br-notu
			for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
				sp <- venn_tree[ndg,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	branch_changes[br] <- 1
	char_changes[ch] <- char_changes[ch]+1
	# now, flip over branches
	pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
	pabm[ch,char_changes[ch]] <- br
	b <- b+1
	}	# end case where branch needed changes

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
simcharacter_compats <- rowSums(simcompmat)-1
simcompat <- vector(length=delta)
simcompat[delta] <- sum(simcharacter_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
simcompat_test <- simcompat	# for debugging.
print_status
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	br <- 0
	br_counter <- 0
	while (br==0 && br_counter<100)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		br_counter <- br_counter+1
		}
	if (br_counter <= 99)	{

		prior_char <- simchmatrix[,ch]
		if (states[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
#		prior_compat <- simcompmat[ch,]
		prior_compat <- simcompmat
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
		#change_character_on_branch(ch,branchings[b])
			}
	# new vector of compatibilities
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	#	simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		simcompmat[,ch] <- simcompmat[ch,] <- new_compat
		simcharacter_compats <- rowSums(simcompmat)-1
		simcompat <- c(simcompat,sum(simcharacter_compats)/2)
		simcompat_test[d] <- total_compatibility(chmatrix=simchmatrix,types,UNKNOWN,INAP);
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompmat)-nchars)/2
		newmatcomp <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	#	simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		if (simcompat[d]>ttl_compat)	{
		# update character compatibilities that might have been altered
			branch_changes[br] <- branch_changes[br]+1
			char_changes[ch] <- char_changes[ch]+1
			simcompmat[ch,] <- simcompmat[,ch] <- new_compat
	#		for (c2 in 1:nchars)	{
	#			if (new_compat[c2]!=prior_compat[c2])	{
	#				if (new_compat[c2]==0)	{
	#					character_compats[c2] <- character_compats[c2]-1
	#					} else if (new_compat[c2]==1)	{
	#					character_compats[c2] <- character_compats[c2]+1
	#					}
	#				}	# case of mismatch
	#			}
			} else	if (simcompat[d]<ttl_compat)	{
				
			}
		}
	xxx <- cbind(simcompat,simcompat_test)
	print(c(xxx[d,],ttl_compat))
	}

	#simcompat[d-1]
	#simcompat[d]
	#d+(ttl_compat-simcompat[d])/(simcompat[d-1]-ttl_compat)
	### START HERE!!!!
return (d+(ttl_compat-simcompat[d])/(simcompat[d-1]-ttl_compat))
}

# simulate character evolution up to a certain compatibility
evolve_up_to_compatibility <- function(init_chmatrix,ttl_compat,venn_tree,branchings,nchars,states,types,UNKNOWN,INAP,repl=-1)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
simchmatrix <- init_chmatrix
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0;
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
#if (maxsteps==-1)	{
#	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
#	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
#	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
#	}
#simcompat <- vector(length=maxsteps)
ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
#	if (b==1)	{
#		ab <- rep(br,branchings[br])
#		}	else	{
#		ab <- c(ab,rep(br,branchings[br]))
#		}
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}
#rm(ab,list=subset(bcc[1,],bcc[1,]>0))

#delta <- sum(char_changes)
#simcompatmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
simcompatmat <- matrix(1,nchars,nchars)
prior_compat <- character_compats <- nchars-1
#simcompat[delta] <- sum(character_compats)/2
d <- 0
simcompat <- c()
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
#			if (d>1)	{
#				simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
#				}	else	simcompat <- ((nchars^2)-nchars)/2
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all branches have change
counter <- 0
if (min(branch_changes[unique_ab])==0)	{
	needy_brs <- unique_ab[branch_changes[unique_ab]==0]
	for (b in 1:length(needy_brs))	{
		br <- needy_brs[b]	### make sure that is in all such routines!!!!
#		if (branch_changes[br]==0)	{
		counter <- counter+1
		d <- d+1
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		prior_compat <- simcompatmat[ch,]
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	#		bcc[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))
if (simcompat[d]<ttl_compat)	{
	### find where we overshot
	
	}	else	{
	counter <- 1
	while (simcompat[d]>ttl_compat)	{
	### while loop added 2017-06-03
		br <- 0
		while (br==0)	{
			ch <- ceiling(nchars*runif(1))
			while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
			c <- char_changes[ch]+1
			br <- pabm[ch,c]
			}
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
	
		if (states[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
	#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompatmat)-nchars)/2
		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		simcompatmat[ch,] <- new_compat
		simcompatmat[,ch] <- new_compat
		simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		counter <- counter+1
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}
		}
	#output <- list(d,simcompatmat)
	}

return ((d-1)+abs(simcompat[d-1]-ttl_compat)/abs(simcompat[d-1]-simcompat[d]))
}
	
# simulate character evolution up to a certain compatibility
emulate_observed_compatibility <- function(init_chmatrix,venn_tree,branchings,chtypes,ttl_compat=NULL,nstates=NULL,maxsteps=NULL,UNKNOWN=-11,INAP=-22,repl=-1)  {
# init_chmatrix: real character matrix
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# ttl_compat: compatibility to which to evolve (this can be left blank)
# nstates: number of states (this can be left blank)
# chtypes: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do (this can be left blank)
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- nrow(init_chmatrix);	# number of observed taxa in tree
nchar <- ncol(init_chmatrix);
nodes <- nrow(venn_tree);
if (is.null(nstates)) nstates <- count_states(init_chmatrix);
if (is.null(maxsteps)) {
  maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(chmatrix=mod_matrix,n_states=nstates);
  maxsteps <- round((sum(maxsteps_per_char)+sum(taxa_scored))/2,0);
  }
for (sp in 1:notu)  init_chmatrix[sp,init_chmatrix[sp,]<0 & !init_chmatrix[sp,] %in% c(UNKNOWN,INAP)] <- UNKNOWN;
if (is.null(ttl_compat)) ttl_compat <- total_compatibility(chmatrix=init_chmatrix,types=chtypes,UNKNOWN=UNKNOWN,INAP=INAP);

simchmatrix <- init_chmatrix;
simchmatrix[simchmatrix>0] <- 0;	# added 2017-08-10
unique_ab <- (1:length(branchings))[branchings>0];
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes);
char_changes <- vector(length=nchars);
simcompat <- vector(length=maxsteps);

for (b in 1:length(unique_ab))	{
	br <- unique_ab[b];
	if (b==1)	{
		ab <- rep(br,branchings[br]);
		}	else	{
		ab <- c(ab,rep(br,branchings[br]));
		}
	}
tb <- length(ab);		# total branches; = branches above, so why both?
fb <- length(unique_ab);		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb);	# branch changes per character
st_rich <- matrix(0,nchars,max(nstates));
mx_ch <- vector(length=nchars);
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1;
		}
	if (nstates[ch]>=2)	mx_ch[ch] <- st_rich[ch,1];
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}
#rm(ab,list=subset(bcc[1,],bcc[1,]>0))
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (nstates[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(nstates[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(nstates[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all states appear
for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		while (nstates[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		if (nstates[ch]==2)	{
			dstates <- c(1,0)
			} else if (chtypes[ch]==0)	{
			dstates <- scramble_multistates(nstates[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
#		bcc[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,nstates,chtypes,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	ch <- ceiling(nchars*runif(1))
	while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
	c <- char_changes[ch]+1
	br <- pabm[ch,c]
	
	prior_char <- simchmatrix[,ch]
	
	if (nstates[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (chtypes[ch]==0)	{
		dstates <- scramble_multistates(nstates[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,];
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]];
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP);
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat));
	if (repl>0 && counter%%10==0)	print(c(d,simcompat[d]));
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]));
	if (simcompat[d]>ttl_compat)	{
	# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1;
		char_changes[ch] <- char_changes[ch]+1;
		simcompmat[ch,] <- new_compat;
		simcompmat[,ch] <- new_compat;
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1;
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1;
					}
				}	# case of mismatch
			}
		}	else if (simcompat[d]<ttl_compat) {
		simchmatrix[,ch] <- prior_char;
		d <- d-1;
		}
	}

output <- list(simchmatrix,steps_v_compat);
names(output) <- c("simulated_matrix","steps_v_compatibility");

return (output);
}

emulate_observed_compatibility_old <- function(init_chmatrix,venn_tree,branchings,chtypes,ttl_compat=NULL,nstates=NULL,maxsteps=NULL,UNKNOWN=-11,INAP=-22,repl=-1)	{
# init_chmatrix: real character matrix
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# ttl_compat: compatibility to which to evolve (this can be left blank)
# nstates: number of states (this can be left blank)
# chtypes: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do (this can be left blank)
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- nrow(init_chmatrix);	# number of observed taxa in tree
nchar <- ncol(init_chmatrix);
nodes <- nrow(venn_tree);
if (is.null(nstates)) nstates <- count_states(init_chmatrix);
if (is.null(maxsteps)) {
  maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(chmatrix=mod_matrix,n_states=nstates);
  maxsteps <- round((sum(maxsteps_per_char)+sum(taxa_scored))/2,0);
  }
for (sp in 1:notu)  init_chmatrix[sp,init_chmatrix[sp,]<0 & !init_chmatrix[sp,] %in% c(UNKNOWN,INAP)] <- UNKNOWN;
if (is.null(ttl_compat)) ttl_compat <- total_compatibility(chmatrix=init_chmatrix,types=chtypes,UNKNOWN=UNKNOWN,INAP=INAP);

simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes);
char_changes <- vector(length=nchars);
simcompat <- vector(length=maxsteps);

for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	if (b==1)	{
		ab <- rep(br,branchings[br])
		}	else	{
		ab <- c(ab,rep(br,branchings[br]))
		}
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(nstates))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (nstates[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}
#rm(ab,list=subset(bcc[1,],bcc[1,]>0))
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (nstates[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(nstates[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(nstates[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all states appear
for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		while (nstates[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		if (nstates[ch]==2)	{
			dstates <- c(1,0)
			} else if (chtypes[ch]==0)	{
			dstates <- scramble_multistates(nstates[ch]);
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
#		bcc[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,nstates,chtypes,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1;
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	ch <- ceiling(nchars*runif(1))
	while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
	c <- char_changes[ch]+1
	br <- pabm[ch,c]
	
	prior_char <- simchmatrix[,ch]
	
	if (nstates[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (chtypes[ch]==0)	{
		dstates <- scramble_multistates(nstates[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
	# update compatibility =NULLmatrix
	# tally=-22 new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (repl>0 && counter%%10==0)	print(c(d,simcompat[d]))
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	if (simcompat[d]>ttl_compat)	{
	# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}
		}	else if (simcompat[d]<ttl_compat) {
		simchmatrix[,ch] <- prior_char
		d <- d-1
		}
	}

#output <- list(d,simcompmat)
return (steps_v_compat)
}

# simulate character evolution replicating (or nearly so!) observed compatibility
### modified from evolve_compatibility_over_N_changes on 2017-08-29
## two steps: 1) race up to compatibility;
##			  2) if overshat, backup and try 100 times till it's found
emuli_observed_compatibility <- function(init_chmatrix,ttl_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=-1)	{
# initial_character matrix
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
nchars <- ncol(init_chmatrix)
notu <- nrow(init_chmatrix)	# number of observed taxa in tree
nodes <- nrow(venn_tree)
simchmatrix <- init_chmatrix
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0
unique_ab <- (1:length(branchings))[branchings>0]

branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)

ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
mx_ch <- vector(length=nchars);
for (ch in 1:nchars)	mx_ch[ch] <- sum(pabm[ch,]>0)
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}

simcompatmat <- matrix(1,nchars,nchars);
prior_compat <- character_compats <- nchars-1;

d <- 0
simcompat <- c()

# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0);

# second, make sure that all branches have change
counter <- 0
if (min(branch_changes[unique_ab])==0)	{
	needy_brs <- unique_ab[branch_changes[unique_ab]==0]
	for (b in 1:length(needy_brs))	{
		br <- needy_brs[b]	### make sure that is in all such routines!!!!
		counter <- counter+1
		d <- d+1
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		prior_compat <- simcompatmat[ch,]
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}
				}
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
#steps_v_compat <- c(d,simcompat[d])
d <- length(simcompat)
if (repl>0)	print(c(repl,d,simcompat[d]))
counter <- 1
while (simcompat[d]>ttl_compat)	{
### while loop added 2017-06-03
#		simchmatrix_last <- simchmatrix
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		}
	prior_char <- simchmatrix[,ch]
	prior_compat <- simcompatmat[ch,]

	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompatmat)-nchars)/2
	d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
	simcompatmat[ch,] <- new_compat
	simcompatmat[,ch] <- new_compat
	simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
	if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
	counter <- counter+1
#	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	# update character compatibilities that might have been altered
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	for (c2 in 1:nchars)	{
		if (new_compat[c2]!=prior_compat[c2])	{
			if (new_compat[c2]==0)	{
				character_compats[c2] <- character_compats[c2]-1
				} else if (new_compat[c2]==1)	{
				character_compats[c2] <- character_compats[c2]+1
				}
			}	# case of mismatch
		}
	}

if (simcompat[d]<ttl_compat)	{
	simchmatrix_closest <- simchmatrix
	simchmatrix_closest[,ch] <- prior_char;
	closest <- simcompat[d-1]-ttl_compat
	max_br <- ncol(pabm)
	}	# set up closest compatibility with too much: we'll use that in a pinch

counter <- 0
while (simcompat[d]!=ttl_compat && counter<100)	{
	# Undo prior run (if we improved, then we are keeping it the same here)
	simcompatmat[,ch] <- simcompatmat[ch,] <- prior_compat
	simchmatrix[,ch] <- prior_char
	branch_changes[br] <- branch_changes[br]-1
	char_changes[ch] <- char_changes[ch]-1
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		brbr <- pabm[ch,c:mx_ch[ch]]
		rem_poss_br <- pabm[ch,c:mx_ch[ch]][pabm[ch,c:mx_ch[ch]]>0]
		if (length(rem_poss_br)>0)	{
			cc <- 0
			while (cc==0)	cc <- ceiling(runif(1)*length(rem_poss_br))
			br <- rem_poss_br[cc]
			}
		}
	prior_char <- simchmatrix[,ch]
	prior_compat <- simcompatmat[ch,]

	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	if (br<=notu  && br>0)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else if (br>notu) {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompatmat)-nchars)/2
#		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
	simcompatmat[ch,] <- new_compat
	simcompatmat[,ch] <- new_compat
	simcompat[d] <- sum(simcompatmat[lower.tri(simcompatmat)])
	if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
	counter <- counter+1
#	steps_v_compat[dd,2] <- simcompat[d]
	# update character compatibilities that might have been altered
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	for (c2 in 1:nchars)	{
		if (new_compat[c2]!=prior_compat[c2])	{
			if (new_compat[c2]==0)	{
				character_compats[c2] <- character_compats[c2]-1
				} else if (new_compat[c2]==1)	{
				character_compats[c2] <- character_compats[c2]+1
				}
			}	# case of mismatch
		}
#	print(c(simcompat[d],ttl_compat))
	
	# if improvement without getting it exactly right
	if (closest > (simcompat[d]-ttl_compat) && (simcompat[d]-ttl_compat)>0)	{
		closest <- simcompat[d]-ttl_compat
		simchmatrix_closest <- simchmatrix
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
		d <- d+1
		simcompat <- c(simcompat,ttl_compat-1)
#		dd <- dd+1
		branch_changes[br] <- branch_changes[br]+1	# these will both be decremented above
		char_changes[ch] <- char_changes[ch]+1		# these will both be decremented above
		}
	}

output <- list(simchmatrix,simcompatmat,d)
names(output) <- c("Simulated_Char_Matrix","Simulated_Compatibility_Matrix","Steps")
return (output)
}

evolve_to_particular_compatibilty <- function(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN=-11,INAP=-22,repl=1)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
simchmatrix <- matrix(0,notu,nchars);
unique_ab <- (1:length(branchings))[branchings>0];
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=maxsteps)

for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	if (b==1)	{
		ab <- rep(br,branchings[br])
		}	else	{
		ab <- c(ab,rep(br,branchings[br]))
		}
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(chmatrix,nchars,notu,branchings,rich,ab) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	pabm[ch,(pcc[ch]+1):tb] <- 0
	}

# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- length(unique(simchmatrix[,ch]))

# makes sure that all branches have 1+ change
#  PUT THIS ROUTINE IN OTHER MODULES!
needy <- unique(ab)[branch_changes[unique(ab)]==0]
nb <- length(needy)
for (b in 1:nb)	{
	# get first branch in need of a change
	br <- needy[b]
	# target characters where it would change earliest
	candidates <- which(pabm==br,arr.ind=TRUE)
	cn <- 1
	ch <- candidates[cn,1]
	while (states[ch]<2
		   || char_changes[ch]>=mx_ch[ch]
		   || char_changes[ch]>=candidates[cn,2]
		   || is.na(match(br,pabm[ch,])))	{
		cn <- cn+1
		ch <- candidates[cn,1]
		}	# make sure that this is an appropriate character
	# determine shifts
	if (states[ch]==2)	{
		dstates <- c(1,0)
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	
		} 
	# routine for tips
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			cc <- simchmatrix[br,ch]
			simchmatrix[br,ch] <- dstates[cc+1]
			} else {
			ndg <- br-notu
			for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
				sp <- venn_tree[ndg,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	branch_changes[br] <- 1
	char_changes[ch] <- char_changes[ch]+1
	# now, flip over branches
	pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
	pabm[ch,char_changes[ch]] <- br
	}	# end case where branch needed changes

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	ch <- ceiling(nchars*runif(1))
	while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
	c <- char_changes[ch]+1
	br <- pabm[ch,c]
	
	prior_char <- simchmatrix[,ch]
	
	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (repl>0 && counter%%10==0)	print(c(d,simcompat[d]))
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	if (simcompat[d]>ttl_compat)	{
	# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}
		}	else if (simcompat[d]<ttl_compat) {
		simchmatrix[,ch] <- prior_char
		d <- d-1
		}
	}

### START HERE!!!!
max_states <- 1+max(simchmatrix)
live_states <- (0:max(simchmatrix))
state_dists <- matrix(nchars,max_states)
for (ch in 1:nchars)
	for (st in 1:states[ch])
		state_dists[ch,st] <- sum(simchmatrix[,ch]==live_states[st])

output <- list(char_changes,character_compats)

return (steps_v_compat)
}

# get a matrix with the same compatibility as an observed matrix
evolve_discrete_characters_to_observed_compatibility <- function(init_chmatrix,venn_tree,branchings,chtypes,ttl_compat=NULL,nstates=NULL,UNKNOWN=-11,INAP=-22,repl=-1,full_output=F)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nstates: number of states for each character
# chtypes: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
nchars <- ncol(init_chmatrix);
notu <- nrow(init_chmatrix);
nodes <- nrow(venn_tree);
if (is.null(ttl_compat))	ttl_compat <- total_compatibility(init_chmatrix,chtypes);
simchmatrix <- init_chmatrix;
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0;
branchings[notu+1] <- 0;
ttl_branches <- length(branchings);
#branchings[135:ttl_branches]
unique_ab <- (1:length(branchings))[branchings>0];
#unique_ab <- unique_ab[!unique_ab %in% (notu+1)]; # do not include basal branch as that shouldn't count.
branch_changes <- vector(length=notu+nodes);
char_changes <- vector(length=nchars);
available_branches <- c();
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	available_branches <- c(available_branches,rep(br,branchings[br]))
	}
tb <- length(available_branches)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb);	# branch changes per character
st_rich <- matrix(0,nchars,max(nstates))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
#		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		if (!init_chmatrix[s,ch] %in% c(UNKNOWN,INAP))	st_rich[ch,1] <- st_rich[ch,1]+1;
		}
	if (nstates[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

total_progeny <- tally_node_richness_from_venn_tree(venn_tree);
#rich <- vector(length=nodes);
#for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0);
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,available_branches,venn_tree); # possible available branches matrix
pcc <- vector(length=nchars);  # possible character changes

# set up branch change order for characters (pabm);
for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)));
	pcc[ch] <- length(x);
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c];
	if (pcc[ch]<length(available_branches))   for (c in (pcc[ch]+1):length(available_branches))   pabm[ch,c] <- 0;
	}

simcompatmat <- matrix(1,nchars,nchars);
prior_compat <- character_compats <- nchars-1;
#simcompat[delta] <- sum(character_compats)/2
d <- 0
simcompat <- c()
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (nstates[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    use_br <- sort(pabm[ch,1:(nstates[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(nstates[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
#				for (r in 1:rich[ndd])	{
				for (r in 1:total_progeny[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
#			if (d>1)	{
#				simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
#				}	else	simcompat <- ((nchars^2)-nchars)/2
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- count_states(simchmatrix);
#simstates <- vector(length=nchars);
#for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0);

# second, make sure that all branches have change
needy_brs <- unique_ab[branch_changes[unique_ab]==0];
counter <- b <- 0;
while (b < length(needy_brs))	{
	b <- b+1;
	br <- needy_brs[b]	### make sure that is in all such routines!!!!
#		if (branch_changes[br]==0)	{
	counter <- counter+1;
	mbap <- which(pabm==br,arr.ind=T);
	mbap <- mbap[order(mbap[,2]),];
	ch <- mbap[1,1];  # grab the "earliest" change that should have happened.
#	ch <- ceiling(nchars*runif(1));
	# shuffle branch forward to be the nth change for that character
	pabm[ch,] <- as.numeric(c(c(pabm[ch,1:char_changes[ch]],br),pabm[ch,][!pabm[ch,] %in% c(pabm[ch,1:char_changes[ch]],br)]));
#	rbind(pabm[ch,],c(c(pabm[ch,1:char_changes[ch]],br),pabm[ch,][!pabm[ch,] %in% c(pabm[ch,1:char_changes[ch]],br)]))
#	pabm[ch,(char_changes[ch]+2):mbap[1,2]] <- pabm[ch,(char_changes[ch]+1):(mbap[1,2]-1)];
#	pabm[ch,char_changes[ch]] <- br;
	# make sure that this is a character that we want to change;
#	while (nstates[ch]<2 || char_changes[ch]>=mx_ch[ch] || simchmatrix[br,ch] %in% c(UNKNOWN,INAP))	ch <- ceiling(nchars*runif(1));
	prior_compat <- simcompatmat[ch,];
	if (nstates[ch]==2)	{
		dstates <- c(1,0);
		} else if (chtypes[ch]==0)	{
		dstates <- scramble_multistates(nstates[ch]);
		} else if (chtypes[ch]==1)	{
		ordered_shift <- c(-1,1)[max(1,ceiling(runif(1)*2))]
		dstates <- ((1:nstates[ch])-1)+ordered_shift;
		dstates[dstates<0] <-0;
		dstates[dstates>=nstates[ch]] <- nstates[ch]-1;
		}
	names(dstates) <- 0:(nstates[ch]-1);
	if (br<=notu)	{
		cc <- simchmatrix[br,ch];
#		rc <- match(cc,names(dstates));
		if (chtypes[ch]==1 && nstates[ch]>2)	{
			if (cc==0)	{
				simchmatrix[br,ch] <- 1;
				} else if (cc==(nstates[ch]-1))	{
				simchmatrix[br,ch] <- nstates[ch]-2;
				} else	{
				simchmatrix[br,ch] <- simchmatrix[br,ch]+ordered_shift;
				}
			} else	{
			simchmatrix[br,ch] <- dstates[cc+1];	# state zero is in cell 1;
			}
		} else {
		ndg <- br-notu;	# get node number;
#		for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
		for (d in 1:total_progeny[ndg])	{
#			simchmatrix[venn_tree[ndg,1:total_progeny[ndg]],ch]
			sp <- venn_tree[ndg,d];
			if (!simchmatrix[sp,ch] %in% c(UNKNOWN,INAP))	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
			}	
		}
			#change_character_on_branch(ch,branchings[b])
	branch_changes[br] <- branch_changes[br]+1;
	char_changes[ch] <- char_changes[ch]+1;
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP);
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat));
	simcompat <- c(simcompat,total_compatibility(chmatrix=simchmatrix,chtypes=chtypes,nstates=nstates,UNKNOWN,INAP));
	#		bcc[ch,char_changes[ch]] <- br
	}	# end case where branch needed changes

#third, tally compatibility at this point
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
delta <- sum(char_changes);
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))
if (simcompat[d]<ttl_compat)	{
	### find where we overshot
	
	}	else	{
	counter <- 1;
	while (simcompat[d]>ttl_compat)	{
	### while loop added 2017-06-03
		br <- 0
		while (br==0)	{
			ch <- ceiling(nchars*runif(1))
			while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
			c <- char_changes[ch]+1
			br <- pabm[ch,c]
			}
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
	
		if (nstates[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (chtypes[ch]==0)	{
			dstates <- scramble_multistates(nstates[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
	#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompatmat)-nchars)/2
		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		simcompatmat[ch,] <- new_compat
		simcompatmat[,ch] <- new_compat
		simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		counter <- counter+1
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1;
		char_changes[ch] <- char_changes[ch]+1;
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1;
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1;
					}
				}	# case of mismatch
			}
		}
	#output <- list(d,simcompatmat)
	}
if (full_output)	{
	output <- list(simchmatrix,simcompatmat,char_changes,branch_changes,steps_v_compat);
	names(output) <- c("character_matrix","compatibility_matrix","character_changes","branch_lengths","changes_vs_compatibility");
	return(output);
	} else	{
	return (simchmatrix);
	}
}

evolve_discrete_characters_to_observed_compatibility_old2 <- function(init_chmatrix,venn_tree,branchings,chtypes,ttl_compat=NULL,nstates=NULL,UNKNOWN=-11,INAP=-22,repl=-1,full_output=F)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nstates: number of states for each character
# chtypes: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
nchars <- ncol(init_chmatrix);
notu <- nrow(init_chmatrix);
nodes <- nrow(venn_tree);
if (is.null(ttl_compat))	ttl_compat <- total_compatibility(init_chmatrix,chtypes);
simchmatrix <- init_chmatrix;
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0;
branchings[notu+1] <- 0;
ttl_branches <- length(branchings);
#branchings[135:ttl_branches]
unique_ab <- (1:length(branchings))[branchings>0];
#unique_ab <- unique_ab[!unique_ab %in% (notu+1)]; # do not include basal branch as that shouldn't count.
branch_changes <- vector(length=notu+nodes);
char_changes <- vector(length=nchars);
available_branches <- c();
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	available_branches <- c(available_branches,rep(br,branchings[br]))
	}
tb <- length(available_branches)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb);	# branch changes per character
st_rich <- matrix(0,nchars,max(nstates))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
#		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		if (!chmatrix[s,ch] %in% c(UNKNOWN,INAP))	st_rich[ch,1] <- st_rich[ch,1]+1;
		}
	if (nstates[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

total_progeny <- tally_node_richness_from_venn_tree(venn_tree);
#rich <- vector(length=nodes);
#for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0);
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,available_branches,venn_tree); # possible available branches matrix
pcc <- vector(length=nchars);  # possible character changes

# set up branch change order for characters (pabm);
for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)));
	pcc[ch] <- length(x);
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c];
	if (pcc[ch]<length(available_branches))   for (c in (pcc[ch]+1):length(available_branches))   pabm[ch,c] <- 0;
	}

simcompatmat <- matrix(1,nchars,nchars);
prior_compat <- character_compats <- nchars-1;
#simcompat[delta] <- sum(character_compats)/2
d <- 0
simcompat <- c()
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (nstates[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    use_br <- sort(pabm[ch,1:(nstates[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(nstates[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
#				for (r in 1:rich[ndd])	{
				for (r in 1:total_progeny[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
#			if (d>1)	{
#				simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
#				}	else	simcompat <- ((nchars^2)-nchars)/2
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- count_states(simchmatrix);
#simstates <- vector(length=nchars);
#for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0);

# second, make sure that all branches have change
needy_brs <- unique_ab[branch_changes[unique_ab]==0];
counter <- b <- 0;
while (b < length(needy_brs))	{
	b <- b+1;
	br <- needy_brs[b]	### make sure that is in all such routines!!!!
#		if (branch_changes[br]==0)	{
	counter <- counter+1;
	mbap <- which(pabm==br,arr.ind=T);
	mbap <- mbap[order(mbap[,2]),];
	ch <- mbap[1,1];  # grab the "earliest" change that should have happened.
#	ch <- ceiling(nchars*runif(1));
	# shuffle branch forward to be the nth change for that character
	pabm[ch,] <- as.numeric(c(c(pabm[ch,1:char_changes[ch]],br),pabm[ch,][!pabm[ch,] %in% c(pabm[ch,1:char_changes[ch]],br)]));
#	rbind(pabm[ch,],c(c(pabm[ch,1:char_changes[ch]],br),pabm[ch,][!pabm[ch,] %in% c(pabm[ch,1:char_changes[ch]],br)]))
#	pabm[ch,(char_changes[ch]+2):mbap[1,2]] <- pabm[ch,(char_changes[ch]+1):(mbap[1,2]-1)];
#	pabm[ch,char_changes[ch]] <- br;
	# make sure that this is a character that we want to change;
#	while (nstates[ch]<2 || char_changes[ch]>=mx_ch[ch] || simchmatrix[br,ch] %in% c(UNKNOWN,INAP))	ch <- ceiling(nchars*runif(1));
	prior_compat <- simcompatmat[ch,];
	if (nstates[ch]==2)	{
		dstates <- c(1,0);
		} else if (chtypes[ch]==0)	{
		dstates <- scramble_multistates(nstates[ch]);
		} else if (chtypes[ch]==1)	{
		ordered_shift <- c(-1,1)[max(1,ceiling(runif(1)*2))]
		dstates <- ((1:nstates[ch])-1)+ordered_shift;
		dstates[dstates<0] <-0;
		dstates[dstates>=nstates[ch]] <- nstates[ch]-1;
		}
	names(dstates) <- 0:(nstates[ch]-1);
	if (br<=notu)	{
		cc <- simchmatrix[br,ch];
#		rc <- match(cc,names(dstates));
		if (chtypes[ch]==1 && nstates[ch]>2)	{
			if (cc==0)	{
				simchmatrix[br,ch] <- 1;
				} else if (cc==(nstates[ch]-1))	{
				simchmatrix[br,ch] <- nstates[ch]-2;
				} else	{
				simchmatrix[br,ch] <- simchmatrix[br,ch]+ordered_shift;
				}
			} else	{
			simchmatrix[br,ch] <- dstates[cc+1];	# state zero is in cell 1;
			}
		} else {
		ndg <- br-notu;	# get node number;
#		for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
		for (d in 1:total_progeny[ndg])	{
#			simchmatrix[venn_tree[ndg,1:total_progeny[ndg]],ch]
			sp <- venn_tree[ndg,d];
			if (!simchmatrix[sp,ch] %in% c(UNKNOWN,INAP))	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
			}	
		}
			#change_character_on_branch(ch,branchings[b])
	branch_changes[br] <- branch_changes[br]+1;
	char_changes[ch] <- char_changes[ch]+1;
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP);
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat));
	simcompat <- c(simcompat,total_compatibility(chmatrix=simchmatrix,chtypes=chtypes,nstates=nstates,UNKNOWN,INAP));
	#		bcc[ch,char_changes[ch]] <- br
	}	# end case where branch needed changes

#third, tally compatibility at this point
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
delta <- sum(char_changes);
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))
if (simcompat[d]<ttl_compat)	{
	### find where we overshot
	
	}	else	{
	counter <- 1;
	while (simcompat[d]>ttl_compat)	{
	### while loop added 2017-06-03
		br <- 0
		while (br==0)	{
			ch <- ceiling(nchars*runif(1))
			while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
			c <- char_changes[ch]+1
			br <- pabm[ch,c]
			}
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
	
		if (nstates[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (chtypes[ch]==0)	{
			dstates <- scramble_multistates(nstates[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
	#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompatmat)-nchars)/2
		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		simcompatmat[ch,] <- new_compat
		simcompatmat[,ch] <- new_compat
		simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		counter <- counter+1
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1;
		char_changes[ch] <- char_changes[ch]+1;
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1;
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1;
					}
				}	# case of mismatch
			}
		}
	#output <- list(d,simcompatmat)
	}
if (full_output)	{
	output <- c(simchmatrix,simcompatmat,char_changes,branch_changes,steps_v_compat);
	names(output) <- c("character_mattrix","compatibility_matrix","character_changes","branch_lengths","changes_vs_compatibility");
	return(output);
	} else	{
	return (simchmatrix);
	}
}

evolve_discrete_characters_to_observed_compatibility_old <- function(init_chmatrix,venn_tree,branchings,chtypes,ttl_compat=NULL,nstates=NULL,UNKNOWN=-11,INAP=-22,repl=-1,keep_polymorphs=F)	{
# init_chmatrix: original data matrix to be emulated
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# chtypes: types for each character (0: unordered; 1: ordered)
# ttl_compat: compatibility to which to evolve
# nstates: number of states for each character
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
nchars <- ncol(init_chmatrix);
notu <- nrow(init_chmatrix);
nodes <- nrow(venn_tree);
if (is.null(ttl_compat))	ttl_compat <- total_compatibility(init_chmatrix,chtypes);
simchmatrix <- init_chmatrix;
if (keep_polymorphs) {
	obs_polymorphs <- unique(simchmatrix[simchmatrix<0][!simchmatrix[simchmatrix<0] %in% c(UNKNOWN,INAP)]);
	rr <- cc <- c();
	op <- 0;
	while (op < length(obs_polymorphs))	{
		op <- op+1;
		rr <- c(rr,which(simchmatrix == obs_polymorphs[op],arr.ind=TRUE)[,1]);
		cc <- c(cc,which(simchmatrix == obs_polymorphs[op],arr.ind=TRUE)[,2]);
		}
	if (length(obs_polymorphs)>0)	{
		polymorph_info <- data.frame(taxon=as.numeric(rr),char=as.numeric(cc),states=as.numeric(cc))
		for (i in 1:nrow(polymorph_info))
			polymorph_info$states[i] <- length(unravel_polymorph(simchmatrix[rr[i],cc[i]]))
		}
	}

for (ch in 1:nchars)	{
	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0;
	if (!keep_polymorphs)
		simchmatrix[((simchmatrix[,ch] < 0) + (!simchmatrix[,ch] %in% c(UNKNOWN,INAP)))==2,ch] <- UNKNOWN;
	}

unique_ab <- (1:length(branchings))[branchings>0]
branch_changes <- vector(length=notu+nodes);
char_changes <- vector(length=nchars);
ab <- c();
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b];
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb);	# branch changes per character
st_rich <- matrix(0,nchars,max(nstates));
mx_ch <- vector(length=nchars)
#accersi_maximum_parsimony_steps_per_character(init_chmatrix,nstates)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
#		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		if (!simchmatrix[s,ch] %in% c(UNKNOWN,INAP) && simchmatrix[s,ch]>=0)	{
			st_rich[ch,1] <- st_rich[ch,1]+1;
			} #else if (!simchmatrix[s,ch] %in% c(UNKNOWN,INAP))	{
#			st_rich
#			}
		}
	if (nstates[ch]>=2)	mx_ch[ch] <- st_rich[ch,1];
	}
mx_ch[mx_ch>sum(branchings>0)] <- sum(branchings>0);

total_progeny <- tally_node_richness_from_venn_tree(venn_tree);
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree); # possible available branches matrix
pcc <- vector(length=nchars);  # possible character changes

# set up branch change order for characters (pabm);
for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)));
	pcc[ch] <- length(x);
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c];
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0;
	}

simcompatmat <- matrix(1,nchars,nchars);
prior_compat <- character_compats <- nchars-1;
#simcompat[delta] <- sum(character_compats)/2
d <- 0;
simcompat <- c();
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (nstates[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    use_br <- sort(pabm[ch,1:(nstates[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(nstates[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu;
				for (r in 1:total_progeny[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
#			if (d>1)	{
#				simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
#				}	else	simcompat <- ((nchars^2)-nchars)/2
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- count_states(simchmatrix);
#simstates <- vector(length=nchars);
#for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0);

# second, make sure that all branches have change
needy_brs <- unique_ab[branch_changes[unique_ab]==0];
counter <- b <- 0;
while (b < length(needy_brs))	{
	b <- b+1;
	br <- needy_brs[b]	### make sure that is in all such routines!!!!
#		if (branch_changes[br]==0)	{
	counter <- counter+1;
	mbap <- which(pabm==br,arr.ind=TRUE);
	mbap <- mbap[order(mbap[,2]),];
	ch <- mbap[1,1];
#	ch <- ceiling(nchars*runif(1));
	pabm[ch,(char_changes[ch]+2):mbap[1,2]] <- pabm[ch,(char_changes[ch]+1):(mbap[1,2]-1)];
	pabm[ch,char_changes[ch]] <- br;
	# make sure that this is a character that we want to change;
#	while (nstates[ch]<2 || char_changes[ch]>=mx_ch[ch] || simchmatrix[br,ch] %in% c(UNKNOWN,INAP))	ch <- ceiling(nchars*runif(1));
	prior_compat <- simcompatmat[ch,];
	if (nstates[ch]==2)	{
		dstates <- c(1,0);
		} else if (chtypes[ch]==0)	{
		dstates <- scramble_multistates(nstates[ch]);
		} 
	names(dstates) <- 0:(nstates[ch]-1);
	if (br<=notu)	{
		cc <- simchmatrix[br,ch];
#		rc <- match(cc,names(dstates));
		simchmatrix[br,ch] <- dstates[cc+1];	# state zero is in cell 1;
		} else {
		ndg <- br-notu;	# get node number;
#		for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
		for (d in 1:total_progeny[ndg])	{
#			simchmatrix[venn_tree[ndg,1:total_progeny[ndg]],ch]
			sp <- venn_tree[ndg,d];
			if (!simchmatrix[sp,ch] %in% c(UNKNOWN,INAP))	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
			}	
		}
			#change_character_on_branch(ch,branchings[b])
	branch_changes[br] <- branch_changes[br]+1;
	char_changes[ch] <- char_changes[ch]+1;
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP);
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat));
	simcompat <- c(simcompat,total_compatibility(simchmatrix,chtypes,UNKNOWN,INAP));
	#		bcc[ch,char_changes[ch]] <- br
	}	# end case where branch needed changes

#third, tally compatibility at this point
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
delta <- sum(char_changes);
steps_v_compat <- data.frame(steps=as.numeric(delta),compatibility=as.numeric(simcompat[delta]));
if (repl>0)	print(c(repl,delta,simcompat[delta]))
if (simcompat[d]<ttl_compat)	{
	### find where we overshot
	
	}	else	{
	counter <- 1;
	while (simcompat[d]>ttl_compat)	{
	### while loop added 2017-06-03
		br <- 0
		while (br==0)	{
			ch <- ceiling(nchars*runif(1))
			while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
			c <- char_changes[ch]+1;
			br <- pabm[ch,c];
			}
		prior_char <- simchmatrix[,ch];
		prior_compat <- simcompatmat[ch,];
	
		if (nstates[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (chtypes[ch]==0)	{
			dstates <- scramble_multistates(nstates[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
	#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]];
				} else {
				ndh <- br-notu;
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d];
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		new_compat <- compatibility_of_a_character(ch,simchmatrix,nstates,chtypes,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompatmat)-nchars)/2
		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		simcompatmat[ch,] <- new_compat
		simcompatmat[,ch] <- new_compat
		simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		counter <- counter+1
		steps_v_compat <- rbind(steps_v_compat,data.frame(steps=as.numeric(d),compatibility=as.numeric(simcompat[d])));
		# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1;
		char_changes[ch] <- char_changes[ch]+1;
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1;
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1;
					}
				}	# case of mismatch
			}
		}
	#output <- list(d,simcompatmat)
	}

output <- list(simchmatrix,simcompatmat,character_compats,char_changes,steps_v_compat);
names(output) <- c("Simulated_Character_Matrix","Simulated_Compatibility_Matrix","Simulated_Character_Compatibility","Simulated_Character_Changes","Simulated_Compatibility_per_Steps");
return (output)
}

# routine to list all branches on which character change can be reconstructed for each character
get_possible_branches_for_all_characters <- function(simchmatrix,branchings,available_branches,venn_tree,UNKNOWN=-11,INAP=-22)	{
# find branches on which each character can change.  (Accommodates unknowns & inapplicables)
# simchmatrix: character matrix
# branchings: 
# rich: number of nodes
# available_branches: number of available branches (e.g., some scored states above it)
nchars <- ncol(simchmatrix);
notu <- nrow(simchmatrix);
poss_brs <- length(available_branches);
pabm <- matrix(0,nchars,poss_brs);	# character by available branches matrix
for (ch in 1:nchars)	{
	a <- 1;
	for (b in 1:poss_brs)	{
		if (available_branches[b]<=notu)	{
			s <- available_branches[b]
			if (simchmatrix[s,ch]!=UNKNOWN && simchmatrix[s,ch]!=INAP)	{
				pabm[ch,a]<-s
				a <- a+1
				}	# case where species is scored
			} else {
			nda <- available_branches[b]-notu;  ### PROBLEM APPEARS HERE!!!!
			accept <- 0
##			  for (r in 1:rich[nda])	{
			for (r in 1:sum(venn_tree[nda,]>0))	{
				s <- venn_tree[nda,r]
				if (simchmatrix[s,ch]!=UNKNOWN && simchmatrix[s,ch]!=INAP)	{
					accept <- 1
				  	r <- sum(venn_tree[nda,])
					}	# case where a clade member is scored
				}
			if (accept==1)	{
				pabm[ch,a]<-nda+notu
				a <- a+1
				}
			}	# end case of node
		}
	}
char_names <- vector(length=nchars)
for (c in 1:nchars)	{
	if (nchars<100)	{
		if (c<10)	{
			char_names[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_names[c] <- paste("ch_",c,sep="")	
			}
		}	else	{
		if (c<10)	{
			char_names[c] <- paste("ch_00",c,sep="")
			} else if (c<100)	{
			char_names[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_names[c] <- paste("ch_",c,sep="")
			}
		}
	}

brnch_names <- vector(length=poss_brs)
for (b in 1:poss_brs)	{
	if (poss_brs<100)	{
		if (b<10)	{
			brnch_names[b] <- paste("br_0",b,sep="")
			}	else	{
			brnch_names[b] <- paste("br_",b,sep="")	
			}
		}	else	{
		if (b<10)	{
			brnch_names[b] <- paste("br_00",b,sep="")
			} else if (c<100)	{
			brnch_names[b] <- paste("br_0",b,sep="")
			}	else	{
			brnch_names[b] <- paste("br_",b,sep="")
			}
		}
	}
rownames(pabm) <- char_names
colnames(pabm) <- brnch_names
return(pabm)
}

# character change routines ####
# routine to create unique a <-> b <-> c transitions
scramble_multistates <- function(ch_states)	{
dstates <- vector(length=ch_states)
for (i in 1:ch_states)	dstates[i]<-i-1
for (i in 1:ch_states)	{
	p <- i+ceiling((ch_states-i)*runif(1))
	q <- dstates[i]
	dstates[i]<-dstates[p]
	dstates[p]<-q
	}
return(dstates)
}

flip_binaries <- function(ch_vector)	{
new_vector <- ch_vector;
new_vector[ch_vector==0] <- 1;
new_vector[ch_vector==1] <- 0;
return(new_vector);
}

evolve_unordered_multistate <- function(ch_vector,ch_states)	{
redone <- scramble_multistates(ch_states);
nch_vector <- ch_vector;
for (st in 1:ch_states)	nch_vector[ch_vector==(st-1)] <- redone[st];
return(nch_vector);
}

# routine to create P[compatibility | steps]
# Ugh: I messed with this by mistake.  Restore it!
accersi_expected_compatibility_given_steps <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=F,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation=TRUE,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=F,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
#	simulation$Branchings[(notu+1):(notu+nNodes)]
	stc <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals,UNKNOWN,INAP,repl=0)
#	pt_done <- 2
#	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
	sim_results <- rbind(sim_results,new_results)
	}
colnames(sim_results) <- 1:maxsteps
#sim_results <- sim_results/runs
return(sim_results)
}

# routine to create P[compatibility | steps] by finding X-1 & X steps crossing observed compatibility
accersi_ave_steps_to_reach_compatibility <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=F,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
obs_comp_matrix <- compatibility_matrix(init_chmatrix,states,types,UNKNOWN,INAP)
obs_compat <- sum(obs_comp_matrix[lower.tri(obs_comp_matrix)])
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation=TRUE,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=F,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
	stc <- evolve_to_particular_compatibilty(init_chmatrix,ttl_compat=obs_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=0)
	sim_results <- c(sim_results,stc)
#	simulation$Branchings[(notu+1):(notu+nNodes)]
#	stc <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals,UNKNOWN,INAP,repl=0)
#	pt_done <- 2
#	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
#	stc <- evolve_to_particular_compatibilty(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=1)
#	sim_results <- rbind(sim_results,new_results)
	}
#colnames(sim_results) <- 1:maxsteps
#sim_results <- sim_results/runs
return(sim_results)
}

# routine to create P[compatibility | steps]
accersi_most_prob_compatibility_given_steps <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=F,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP,repl=-1)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation=TRUE,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=F,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
#	simulation$Branchings[(notu+1):(notu+nNodes)]
	stc <- evolve_up_to_compatibility(ttl_compat,init_chmatrix,venn_tree,branchings,nchars,states,types,UNKNOWN,INAP,repl)
#	pt_done <- 2
	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
#	stc <- evolve_to_particular_compatibilty(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=1)
	sim_results <- rbind(sim_results,new_results)
	}
colnames(sim_results) <- 1:maxsteps
sim_results <- sim_results/runs
return(sim_results)
}

accersi_prob_compatibility_given_steps <- function(chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=F,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP)	{
notu <- dim(chmatrix)[1]
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(chmatrix,states)
	maxsteps <- (sum(maxsteps_per_char)+sum(taxa_scored))/2
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation=TRUE,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=F,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
#	simulation$Branchings[(notu+1):(notu+nNodes)]
	stc <- evolve_compatibility_over_N_changes(N=maxsteps,chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals,UNKNOWN,INAP,repl=0)
#	pt_done <- 2
	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
#	stc <- evolve_to_particular_compatibilty(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=1)
	sim_results <- rbind(sim_results,new_results)
	}
colnames(sim_results) <- 1:maxsteps
return(sim_results)
}

# routine to find standing richness
tally_richness_from_continuous_ranges <- function(ranges,temp_prec=0.1)	{
stg_ranges <- ceiling(ranges/temp_prec)
stg_richness <- vector(length=max(stg_ranges))
etus <- dim(stg_ranges)[1]
#test <- c()
for (i in 1:etus)	{
	stg_richness[(stg_ranges[i,1]:stg_ranges[i,2])] <- stg_richness[(stg_ranges[i,1]:stg_ranges[i,2])]+1
#	if (stg_ranges[i,1]<=present && stg_ranges[i,2]>=present)
#		test <- c(test,i)
	}
return(stg_richness)
}
