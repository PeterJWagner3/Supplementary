# routines for different types of disparity analyses.

pairwise_dissimilarity <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	uchmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
uchmatrix <- unique(chmatrix)
unotu <- nrow(uchmatrix)
nchars <- ncol(uchmatrix)
udis_matrix <- matrix(0,unotu,unotu)
for (sp1 in 1:(unotu-1))	{
#	print(sp1)
	udis_matrix[sp1,sp1] <- 0.0
	for (sp2 in (sp1+1):unotu)	{
		diss <- num <- 0	# diss: # differences; num: number of comparisons
		for (ch in 1:nchars)	{
			if (((uchmatrix[sp1,ch]!=UNKNOWN && uchmatrix[sp1,ch]!=INAP) && (uchmatrix[sp2,ch]!=UNKNOWN && uchmatrix[sp2,ch]!=INAP)))	{
				num <- num+1
				if (uchmatrix[sp1,ch]!=uchmatrix[sp2,ch])	{
					## different with single scored states
					if (uchmatrix[sp1,ch]>=0 && uchmatrix[sp2,ch]>=0)	{
						if (types[ch]==0 || states[ch]==2)	{
							diss <- diss+1
							} else {
							if (weight_ordered==FALSE)	{
								x <- abs(uchmatrix[sp1,ch]-uchmatrix[sp2,ch])
								} else if (weight_ordered==TRUE)	{
								x <- abs(uchmatrix[sp1,ch]-uchmatrix[sp2,ch])/(states[ch]-1)
								}
							diss <- diss+x
							}	# end case of ordered character
						} else { # end case of two non-polymorphsics
						if (uchmatrix[sp1,ch]<0)	{
							cvec1 <- unravel_polymorph(uchmatrix[sp1,ch])
							} else {
							cvec1 <- uchmatrix[sp1,ch]
							}
						if (uchmatrix[sp2,ch]<0)	{
							cvec2 <- unravel_polymorph(uchmatrix[sp2,ch])
							} else {
							cvec2 <- uchmatrix[sp2,ch]	
							}
						if (length(cvec1)<=length(cvec2))	{
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec2 %in% cvec1)/length(cvec2)
#							x <- match(cvec1,cvec2)
#							if (length(cvec1)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec1)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1)
#							x <- match(cvec2,cvec1)		# pick up here!
#							if (length(cvec2)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec2)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}			
						} # end case of 1 or 2 polymorphsics
					} # end case of disagreement
				} # end case of two coded characters
			}	# end going through characters
		udis_matrix[sp2,sp1] <- udis_matrix[sp1,sp2] <- diss/num
		}	# end comparison between sp2 & sp1
	}	#end going throuch species

notu <- nrow(chmatrix)
dis_matrix <- matrix(0,notu,notu)
xxx <- prodlim::row.match(x=as.data.frame(chmatrix),table=as.data.frame(uchmatrix));
for (u in 1:(unotu-1))	{
	yyy <- (1:notu)[xxx %in% u]
	for (u2 in 2:unotu)	{
		zzz <- (1:notu)[xxx %in% u2]
		dis_matrix[zzz,yyy] <- dis_matrix[yyy,zzz] <- udis_matrix[u,u2]
		}
	}
return (dis_matrix)
}

pairwise_similarity_discrete <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11, INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
notu <- nrow(chmatrix)
nchars <- ncol(chmatrix)
sim_matrix <- matrix(0,notu,notu)
for (sp1 in 1:(notu-1))	{
	sim_matrix[sp1,sp1] <- 0.0
	for (sp2 in (sp1+1):notu)	{
		diss <- num <- 0	# diss: # differences; num: number of comparisons
		for (ch in 1:nchars)	{
			if (((chmatrix[sp1,ch]!=UNKNOWN && chmatrix[sp1,ch]!=INAP) && (chmatrix[sp2,ch]!=UNKNOWN && chmatrix[sp2,ch]!=INAP)))	{
				num <- num+1
				if (chmatrix[sp1,ch]!=chmatrix[sp2,ch])	{
					## different with single scored states
					if (chmatrix[sp1,ch]>=0 && chmatrix[sp2,ch]>=0)	{
						if (types[ch]==0 || states[ch]==2)	{
							diss <- diss+1
							} else {
							if (weight_ordered==FALSE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])
								} else if (weight_ordered==TRUE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])/(states[ch]-1)
								}
							diss <- diss+x
							}	# end case of ordered character
						} else { # end case of two non-polymorphsics
						if (chmatrix[sp1,ch]<0)	{
							cvec1 <- unravel_polymorph(chmatrix[sp1,ch])
							} else {
							cvec1 <- chmatrix[sp1,ch]
							}
						if (chmatrix[sp2,ch]<0)	{
							cvec2 <- unravel_polymorph(chmatrix[sp2,ch])
							} else {
							cvec2 <- chmatrix[sp2,ch]	
							}
						if (length(cvec1)<=length(cvec2))	{
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec2 %in% cvec1)/length(cvec2)
#							x <- match(cvec1,cvec2)
#							if (length(cvec1)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec1)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1)
#							x <- match(cvec2,cvec1)		# pick up here!
#							if (length(cvec2)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec2)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}			
						} # end case of 1 or 2 polymorphsics
					} # end case of disagreement
				} # end case of two coded characters
			}	# end going through characters
		sim_matrix[sp2,sp1] <- sim_matrix[sp1,sp2] <- (num-diss)/num
		}	# end comparison between sp2 & sp1
	}	#end going throuch species
return (sim_matrix)
}

pairwise_differences_discrete <- function(chmatrix,UNKNOWN=-11,INAP=-22,progress_bar=T)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#	print(UNKNOWN)
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);
diff_matrix <- matrix(0,notu,notu);
if (progress_bar)	{
	progress <- round(notu*notu*(1:100)/100,0);
	comps <- 0;
	last_update <- 0;
	print("");
	}
for (sp1 in 1:(notu-1))	{
	scored_1 <- (1:nchars)[chmatrix[sp1,]!= UNKNOWN];
	scored_1 <- scored_1[chmatrix[sp1,scored_1]!= INAP];
	if (length(scored_1)>0)	{
		for (sp2 in (sp1+1):notu)	{
			scored_2 <- (1:nchars)[chmatrix[sp2,]!= UNKNOWN];
			scored_2 <- scored_2[chmatrix[sp2,scored_2]!= INAP];
			scored_both <- scored_1[scored_1 %in% scored_2]
			if (length(scored_both)>0)	{
				diff_matrix[sp1,sp2] <- diff_matrix[sp2,sp1] <- sum(chmatrix[sp1,scored_both]!=chmatrix[sp2,scored_both]);
				}
			}
		} else	{
		for (sp2 in (sp1+1):notu)	diff_matrix[sp1,sp2] <- diff_matrix[sp2,sp1] <- 0;
		}
	if (progress_bar)	{
		comps <- comps+2*notu;
		if (sum(comps>progress)>last_update)	{
			last_update <- sum(comps>progress);
			if (last_update<10)	{
				console_update <- paste("0",last_update,"%",sep="");
				} else {
				console_update <- paste(last_update,"%",sep="");
				}
			cat('\b\b\b\b',console_update);
			}
		}
	}
return (diff_matrix);
}

pairwise_differences_and_contrasts_discrete <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#	print(UNKNOWN)
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);
diff_matrix <- matrix(0,notu,notu);
comp_matrix <- matrix(0,notu,notu);
for (sp1 in 1:(notu-1))	{
	scored_1 <- (1:nchars)[chmatrix[sp1,]!= UNKNOWN];
	scored_1 <- scored_1[chmatrix[sp1,scored_1]!= INAP];
	if (length(scored_1)>0)	{
		for (sp2 in (sp1+1):notu)	{
			scored_2 <- (1:nchars)[chmatrix[sp2,]!= UNKNOWN];
			scored_2 <- scored_2[chmatrix[sp2,scored_2]!= INAP];
			scored_both <- scored_1[scored_1 %in% scored_2]
			if (length(scored_both)>0)	{
				diff_matrix[sp1,sp2] <- diff_matrix[sp2,sp1] <- sum(chmatrix[sp1,scored_both]!=chmatrix[sp2,scored_both]);
				comp_matrix[sp1,sp2] <- comp_matrix[sp2,sp1] <- length(scored_both);
				}
			}
		} else	{
		for (sp2 in (sp1+1):notu)	{
			diff_matrix[sp1,sp2] <- diff_matrix[sp2,sp1] <- comp_matrix[sp1,sp2] <- comp_matrix[sp2,sp1] <- 0;
			}
		}
	}

output <- list(diff_matrix,comp_matrix);
names(output) <- c("differences","comparisons");
return (output);
}

gower_transform <- function(dist_matrix)	{
dmat <- -0.5*dist_matrix
ave_dist <- colSums(dmat)/nrow(dmat)
ave <- mean(ave_dist)

aii <- array(1,c(nrow(dmat),nrow(dmat)))*ave_dist
ajj <- base::t(aii)
gower <- dmat-(aii+ajj)+ave

return (gower)
}

extract_off_diagonal <- function(orig_matrix)	{
offdiagonal <- c();
for (a in 1:(nrow(orig_matrix)-1))
	for (b in (a+1):nrow(orig_matrix))
		offdiagonal <- c(offdiagonal,orig_matrix[a,b]);
return(offdiagonal);
}

cumulative_disparity <- function(pairwise_dissimilarities,first_appearances_ma)	{
notu <- nrow(pairwise_dissimilarities);
if (first_appearances_ma[1]<0)
	first_appearances_ma <- -1*first_appearances_ma;
appearance_order <- sort(unique(first_appearances_ma),decreasing = TRUE);

ttl_disparity <- sum(pairwise_dissimilarities)/(notu*(notu-1));
incr_disparity <- numeric();
for (fas in 1:length(appearance_order))	{
	cum_taxa <- (1:notu)[first_appearances_ma>=appearance_order[fas]];
	cnotu <- length(cum_taxa)
	if (cnotu>1)
		incr_disparity <- c(incr_disparity,sum(pairwise_dissimilarities[cum_taxa,cum_taxa])/(cnotu*(cnotu-1)));
	}
if (sum(first_appearances_ma==appearance_order[1])>1)	{
	output <- data.frame(ma=appearance_order,cum_disparity=incr_disparity)
	} else	{
	output <- data.frame(ma=appearance_order[2:length(appearance_order)],cum_disparity=incr_disparity)
	}
return(output);
}

pairwise_disimilarity_discrete <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11, INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
notu <- nrow(chmatrix)
nchars <- ncol(chmatrix)
dis_matrix <- matrix(0,notu,notu)
for (sp1 in 1:(notu-1))	{
	dis_matrix[sp1,sp1] <- 0.0
	for (sp2 in (sp1+1):notu)	{
		diss <- num <- 0	# diss: # differences; num: number of comparisons
		for (ch in 1:nchars)	{
			if (((chmatrix[sp1,ch]!=UNKNOWN && chmatrix[sp1,ch]!=INAP) && (chmatrix[sp2,ch]!=UNKNOWN && chmatrix[sp2,ch]!=INAP)))	{
				num <- num+1
				if (chmatrix[sp1,ch]!=chmatrix[sp2,ch])	{
					## different with single scored states
					if (chmatrix[sp1,ch]>=0 && chmatrix[sp2,ch]>=0)	{
						if (types[ch]==0 || states[ch]==2)	{
							diss <- diss+1
							} else {
							if (weight_ordered==FALSE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])
								} else if (weight_ordered==TRUE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])/(states[ch]-1)
								}
							diss <- diss+x
							}	# end case of ordered character
						} else { # end case of two non-polymorphsics
						if (chmatrix[sp1,ch]<0)	{
							cvec1 <- unravel_polymorph(chmatrix[sp1,ch])
							} else {
							cvec1 <- chmatrix[sp1,ch]
							}
						if (chmatrix[sp2,ch]<0)	{
							cvec2 <- unravel_polymorph(chmatrix[sp2,ch])
							} else {
							cvec2 <- chmatrix[sp2,ch]	
							}
						if (length(cvec1)<=length(cvec2))	{
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec2 %in% cvec1)/length(cvec2)
#							x <- match(cvec1,cvec2)
#							if (length(cvec1)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec1)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1)
#							x <- match(cvec2,cvec1)		# pick up here!
#							if (length(cvec2)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec2)>1)	{
#								x <- clear_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}			
						} # end case of 1 or 2 polymorphsics
					} # end case of disagreement
				} # end case of two coded characters
			}	# end going through characters
		dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- diss/num
		}	# end comparison between sp2 & sp1
	}	#end going throuch species
return (dis_matrix)
}

#disparity from landmark data;
calculate_distance_between_taxa_given_landmarks <- function(adjusted_landmarks,ndimn)	{
notu <- nrow(adjusted_landmarks);
nland <- ncol(adjusted_landmarks)/ndimn;

distances <- array(0,dim=c(notu,notu));	# calculate Euclidean distances
for (tx1 in 1:(notu-1))	{
	for (tx2 in (tx1+1):notu)	{
		for (nl in 1:nland)	{
			this_land <- ((ndimn*(nl-1))+1):(ndimn*nl);
			dist <- 0;
			for (nd in 1:ndimn)
				dist <- dist+(adjusted_landmarks[tx1,this_land[nd]]-adjusted_landmarks[tx2,this_land[nd]])^2;
			distances[tx2,tx1] <- distances[tx1,tx2] <- distances[tx1,tx2]+sqrt(dist);
			}
		}
	}
rownames(distances) <- colnames(distances) <- rownames(adjusted_landmarks);
return(distances);
}

bootstrap_mania <- function(data_vector)	{
bootstrapped <- data_vector[ceiling(length(data_vector)*runif(length(data_vector)))];
return(data.frame(mean=as.numeric(mean(bootstrapped)),median=as.numeric(median(bootstrapped))));
}

accersi_distance_for_one_taxon_from_every_other <- function(single_taxon,all_others,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(all_others);
ootus <- nrow(all_others);
single_applicable <- (1:nchars)[!single_taxon %in% c(UNKNOWN,INAP)];
distances <- vector(length=ootus);
for (nn in 1:ootus)	{
	other_applicable <- (1:nchars)[!all_others[nn,] %in% c(UNKNOWN,INAP)];
	applicable <- single_applicable[single_applicable %in% other_applicable];
	distances[nn] <- sum(single_taxon[applicable]!=all_others[nn,applicable]);
	}	
return(distances);
}
