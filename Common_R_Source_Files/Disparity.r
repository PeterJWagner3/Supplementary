# routines for different types of disparity analyses.

pairwise_dissimilarity <- function(chmatrix,types,states=c(),weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#uchmatrix <- unique(chmatrix);
#unotu <- nrow(uchmatrix);
if (length(states)==0)
	states <- count_states(chmatrix)
nchars <- ncol(chmatrix);
notu <- nrow(chmatrix);
dis_matrix <- array(0,dim=c(notu,notu));
colnames(dis_matrix) <- rownames(dis_matrix) <- rownames(chmatrix);
for (sp1 in 1:(notu-1))	{
#	print(sp1)
	dis_matrix[sp1,sp1] <- 0.0;
	coded_sp1 <- (1:nchars)[!chmatrix[sp1,] %in% c(UNKNOWN,INAP)];
	polymorphs_sp1 <- (1:nchars)[chmatrix[sp1,] < 0];
	polymorphs_sp1 <- polymorphs_sp1[polymorphs_sp1 %in% coded_sp1];
	for (sp2 in (sp1+1):notu)	{
		coded_sp2 <- (1:nchars)[!chmatrix[sp2,] %in% c(UNKNOWN,INAP)];
		coded_both <- coded_sp1[coded_sp1 %in% coded_sp2];
		polymorphs_sp2 <- (1:nchars)[chmatrix[sp2,] < 0];
		polymorphs_sp2 <- polymorphs_sp2[polymorphs_sp2 %in% coded_sp2];
		polymorphs_either <- sort(unique(c(polymorphs_sp1,polymorphs_sp2)));
		diss <- c <- 0;	# diss: # differences; num: number of comparisons
		num <- length(coded_both);
#		uchmatrix[c(sp1,sp2),coded_both]
		while (c<length(coded_both))	{
			c <- c+1;
			ch <- coded_both[c];
			if (!ch %in% polymorphs_either)	{
				if (chmatrix[sp1,ch]!=chmatrix[sp2,ch])	{
						if (types[ch]==0 || states[ch]<=2)	{
							diss <- diss+1;
							} else if (weight_ordered) {
							diss <- diss+(abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])/(states[ch]-1));
							} else	{
							diss <- diss+abs(chmatrix[sp1,ch]-chmatrix[sp2,ch]);
							num <- num+states[ch]-2;
							}
						}
				} else	{
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
					if (sum(cvec2 %in% cvec1)==0)	{
						if (types[ch]==0 || states[ch]==2)	{
							diss <- diss+1;
							} else if (weight_ordered)	{
							diss <- diss+(min(abs(cvec2-cvec1))/(states[ch]-1));
							} else	{
							diss <- diss+min(abs(cvec2-cvec1));
							num <- num+states[ch]-2;
							}
						} # end different polymorphism
					} # end polymorphism
				} # end comparisons
		if (num>0)	dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- diss/num;
		}	# end comparison between sp2 & sp1
	}	#end going throuch species
#which(is.na(dis_matrix),arr.ind=T)
return (dis_matrix)
}

pairwise_dissimilarity_old <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	chmatrix: matrix of characters
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

pairwise_similarity_discrete <- function(chmatrix,types,states=c(),weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11, INAP=-22)	{
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
dmat <- -0.5*dist_matrix;
ave_dist <- colSums(dmat)/nrow(dmat);
ave <- mean(ave_dist);

aii <- array(1,c(nrow(dmat),nrow(dmat)))*ave_dist;
ajj <- base::t(aii);
gower <- dmat-(aii+ajj)+ave;
return (gower)
}

center_dist_matrix <- function(dist_matrix) {
n <- nrow(dist_matrix);
One <- matrix(1, n, n);
mat <- diag(n) - One/n;
matrix_centered <- mat %*% dist_matrix %*% mat;
lowtri <- lower.tri(matrix_centered);
matrix_centered[lowtri] <- t(matrix_centered)[lowtri];
return(matrix_centered)
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

pairwise_disimilarity_discrete <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#uchmatrix <- unique(chmatrix);
#unotu <- nrow(uchmatrix);
nchars <- ncol(chmatrix);
dis_matrix <- array(0,dim=c(notu,notu));
colnames(dis_matrix) <- rownames(dis_matrix) <- rownames(chmatrix);
for (sp1 in 1:(notu-1))	{
#	print(sp1)
	dis_matrix[sp1,sp1] <- 0.0;
	coded_sp1 <- (1:nchars)[!chmatrix[sp1,] %in% c(UNKNOWN,INAP)];
	polymorphs_sp1 <- (1:nchars)[chmatrix[sp1,] < 0];
	polymorphs_sp1 <- polymorphs_sp1[polymorphs_sp1 %in% coded_sp1];
	for (sp2 in (sp1+1):notu)	{
		coded_sp2 <- (1:nchars)[!chmatrix[sp2,] %in% c(UNKNOWN,INAP)];
		coded_both <- coded_sp1[coded_sp1 %in% coded_sp2];
		polymorphs_sp2 <- (1:nchars)[chmatrix[sp2,] < 0];
		polymorphs_sp2 <- polymorphs_sp2[polymorphs_sp2 %in% coded_sp2];
		polymorphs_either <- sort(unique(c(polymorphs_sp1,polymorphs_sp2)));
		diss <- c <- 0;	# diss: # differences; num: number of comparisons
		num <- length(coded_both);
#		uchmatrix[c(sp1,sp2),coded_both]
		while (c<length(coded_both))	{
			c <- c+1;
			ch <- coded_both[c];
			if (!ch %in% polymorphs_either)	{
				if (chmatrix[sp1,ch]!=chmatrix[sp2,ch])	{
						if (types[ch]==0 || states[ch]<=2)	{
							diss <- diss+1;
							} else if (weight_ordered) {
							diss <- diss+(abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])/(states[ch]-1));
							} else	{
							diss <- diss+abs(chmatrix[sp1,ch]-chmatrix[sp2,ch]);
							num <- num+states[ch]-2;
							}
						}
				} else	{
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
					if (sum(cvec2 %in% cvec1)==0)	{
						if (types[ch]==0 || states[ch]==2)	{
							diss <- diss+1;
							} else if (weight_ordered)	{
							diss <- diss+(min(abs(cvec2-cvec1))/(states[ch]-1));
							} else	{
							diss <- diss+min(abs(cvec2-cvec1));
							num <- num+states[ch]-2;
							}
						} # end different polymorphism
					} # end polymorphism
				} # end comparisons
		if (num>0)	dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- diss/num;
		}	# end comparison between sp2 & sp1
	}	#end going throuch species
#which(is.na(dis_matrix),arr.ind=T)
return (dis_matrix)
}

pairwise_disimilarity_discrete_old <- function(chmatrix,states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11, INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
notu <- nrow(chmatrix)
states <- count_states(chmatrix);
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
		dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- diss/num;
		}	# end comparison between sp2 & sp1
	}	#end going throuch species
return (dis_matrix)
}

pairwise_dissimilarity_continuous_only <- function(cmatrix,weight_ordered=TRUE,UNKNOWN=-11,INAP=-22)	{
#	cmatrix: matrix of continuous characters
notu <- nrow(cmatrix)
nchars <- ncol(cmatrix)
dis_matrix <- matrix(0,notu,notu)
# if you used blanks for missing continuous data, then these often will be NAs.
for (sp1 in 1:notu)	cmatrix[sp1,is.na(cmatrix[sp1,])] <- UNKNOWN;

if (weight_ordered)	{
	for (ch1 in 1:nchars)	{
		mnm <- min(cmatrix[cmatrix[,ch1]!=UNKNOWN,ch1]);
		mxm <- max(cmatrix[cmatrix[,ch1]!=UNKNOWN,ch1]);
		cmatrix[cmatrix[,ch1]!=UNKNOWN,ch1] <- (cmatrix[cmatrix[,ch1]!=UNKNOWN,ch1]-mnm)/(mxm-mnm);
		}
	}

for (sp1 in 1:(notu-1))	{
#	print(sp1)
	dis_matrix[sp1,sp1] <- 0.0
	for (sp2 in (sp1+1):notu)	{
		diss <- num <- 0	# diss: # differences; num: number of comparisons
		for (ch in 1:nchars)	{
			if (((cmatrix[sp1,ch]!=UNKNOWN && cmatrix[sp1,ch]!=INAP) && (cmatrix[sp2,ch]!=UNKNOWN && cmatrix[sp2,ch]!=INAP)))	{
				num <- num+1;
				diss <- diss + abs(cmatrix[sp1,ch]-cmatrix[sp2,ch]);
				} # end case of two coded characters
			}	# end going through characters
		if (num>0)	{
			dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- diss/num;
			} else	{
			dis_matrix[sp2,sp1] <- dis_matrix[sp1,sp2] <- 0;
			}
		}	# end comparison between sp2 & sp1
	}	#end going throuch species

return (dis_matrix)
}

pairwise_disimilarity_discrete_and_continuous <- function(cmatrix,chmatrix,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11, INAP=-22)	{
#	cmatrix: matrix of continuous characters (measurements)
#	chmatrix: matrix of discrete characters
#	types: 0 for unordered, 1 for ordered
# weight_ordered: if true, then the maximum distance is rescaled from 1.0 for all contrasts.
# UNKNOWN: numerical value denoting unknown (usually -11)
# INAP: numerical value denoting inapplicable (usually -22)
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);
states <- count_states(chmatrix);
nmeas <- ncol(cmatrix);
dis_matrix <- matrix(0,notu,notu);

# if you used blanks for missing continuous data, then these often will be NAs.
for (sp1 in 1:notu)	cmatrix[sp1,is.na(cmatrix[sp1,])] <- UNKNOWN;

if (weight_ordered)	{
	for (ch1 in 1:nchars)	{
		mnm <- min(cmatrix[,cmatrix[,ch1]!=UNKNOWN]);
		mxm <- max(cmatrix[,cmatrix[,ch1]!=UNKNOWN]);
		cmatrix[,cmatrix[,ch1]!=UNKNOWN] <- (cmatrix[,cmatrix[,ch1]!=UNKNOWN]-mnm)/(mxm-mnm);
		}
	}

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
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1)
							}			
						} # end case of 1 or 2 polymorphsics
					} # end case of disagreement
				} # end case of two coded characters
			}	# end going through characters
		for (mm in 1:nmeas)	{
			if (cmatrix[sp1,mm]!=UNKNOWN && cmatrix[sp2,mm]!=UNKNOWN)	{
				diss <- diss+abs(cmatrix[sp1,mm]-cmatrix[sp2,mm])/max(cmatrix[,mm])-min(cmatrix[,mm]);
				num <- num+1;
				}
			}
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

#dissimilarity_matrix <- partition_dissimilarity
bootstrap_disparity_from_matrix <- function(dissimilarity_matrix,runs=500)	{
data_vector <- extract_off_diagonal(dissimilarity_matrix);
data_vector <- data_vector[!is.na(data_vector)];
bootstrapped <- data.frame(mean=as.numeric(1:runs),median=as.numeric(1:runs));
for (i in 1:runs)	{
	baron <- bootstrap_mania(data_vector);
	bootstrapped$mean[i] <- baron$mean;
	bootstrapped$median[i] <- baron$median;
	}
return(bootstrapped);
}

bootstrap_disparity <- function(data_vector,runs=500)	{
bootstrapped <- data.frame(mean=as.numeric(1:runs),median=as.numeric(1:runs));
for (i in 1:runs)	{
	baron <- bootstrap_mania(data_vector);
	bootstrapped$mean[i] <- baron$mean;
	bootstrapped$median[i] <- baron$median;
	}
return(bootstrapped);
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

# this is for assessing association between character states & axis
recode_multistates <- function(old_states,new_code)	{
new_states <- old_states
old_coding <- sort(unique(old_states))
otus <- length(old_states)
for (os in 1:length(old_coding))	{
	new_states[(1:otus)[old_states==old_coding[os]]] <- new_code[os]
	}
return(new_states)
}

# get Kendall's correlations between (re)-ordered states and PCO axes
kendall_correlation_between_PCO_and_states <- function(state_vector,pco_vectors,ctype)	{
options(warn=-1);
kendall_taus <- c()
if (ctype==1)	{
	for (pc in 1:ncol(pco_vectors))
		kendall_taus <- c(kendall_taus,cor.test(state_vector, pco_vectors[,pc],method = "kendall")$estimate)
	} else	{
	char_states <- sort(unique(state_vector))
	notu <- length(state_vector)
	for (pc in 1:ncol(pco_vectors))	{
		median_on_pc <- c()
		for (st in 1:length(char_states))	{
			taxa_w_st <- (1:notu)[state_vector==char_states[st]]
			median_on_pc <- c(median_on_pc,median(pco_vectors[taxa_w_st,pc]))
			}
		nc <- order(median_on_pc)
		new_code <- char_states
		for (i in 1:length(nc))	new_code[nc[i]] <- i
		old_states <- state_vector
		reordered_states <- recode_multistates(old_states,new_code)
		kendall_taus <- c(kendall_taus,cor.test(reordered_states, pco_vectors[,pc],method = "kendall")$estimate)
		}
	}
options(warn=1);
return(kendall_taus)
}

# order unordered multistates by median values along another axis
order_unordered_states_by_loadings <- function(state_vector,pco_vectors)	{
char_states <- sort(unique(state_vector))
notu <- length(state_vector)
ordered_matrix <- c()
for (pc in 1:ncol(pco_vectors))	{
	median_on_pc <- c()
	for (st in 1:length(char_states))	{
		taxa_w_st <- (1:notu)[state_vector==char_states[st]]
		median_on_pc <- c(median_on_pc,median(pco_vectors[taxa_w_st,pc]))
		}
	nc <- order(median_on_pc)
	new_code <- char_states
	for (i in 1:length(nc))	new_code[nc[i]] <- i
	old_states <- state_vector
	reordered_states <- recode_multistates(old_states,new_code)
	ordered_matrix <- cbind(ordered_matrix,reordered_states)
	}
return(ordered_matrix)
}

count_pairwise_comparisons <- function(chmatrix,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#uchmatrix <- unique(chmatrix);
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);
comparisons <- array(0,dim=c(notu,notu));
rownames(comparisons) <- colnames(comparisons) <- rownames(chmatrix);
for (sp1 in 1:(notu-1))	{
	coded_sp1 <- (1:nchars)[!chmatrix[sp1,] %in% c(UNKNOWN,INAP)];
	for (sp2 in (sp1):notu)	{
		coded_sp2 <- (1:nchars)[!chmatrix[sp2,] %in% c(UNKNOWN,INAP)];
		coded_both <- coded_sp1[coded_sp1%in%coded_sp2];
		comparisons[sp1,sp2] <- comparisons[sp2,sp1] <- length(coded_both);
		}
	}
return(comparisons);
}

