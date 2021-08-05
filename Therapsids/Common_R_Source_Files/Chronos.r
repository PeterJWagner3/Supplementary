# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal
# scribere: write
#### Functions to deal with Absolute Time & Time Scales ####

# create a time scale from the finest divisions of time given a time scale with "general and "specific" intervals.
accersi_finest_timescale <- function(chronostrat_units,time_scale,regional_scale = "International",ma_fuzz=0)	{
hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units,time_scale,regional_scale,ma_fuzz);
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[match(unique(finest_chronostrat$bin_first),finest_chronostrat$bin_first),];
return(finest_chronostrat);
}

# get hierarchically arranged timescale for a particular region
# fixed 2020-04-02
# fixed 2021-03-27: set fuzz at 10 Ma
accersi_hierarchical_timescale <- function(chronostrat_units,time_scale,regional_scale="International",ma_fuzz=10)	{
chronostrat <- time_scale[as.character(time_scale$interval) %in% chronostrat_units,];
oldest <- max(time_scale$ma_lb[as.character(time_scale$interval) %in% chronostrat_units])+ma_fuzz;
youngest <- min(time_scale$ma_ub[as.character(time_scale$interval) %in% chronostrat_units])-ma_fuzz;
relv_scales <- unique(c(chronostrat$scale,regional_scale));

xxx <- subset(time_scale,time_scale$ma_lb<=oldest);
chronostrat <- subset(xxx,xxx$ma_ub>=youngest);
chronostrat <- subset(chronostrat,chronostrat$scale==regional_scale);
chronostrat <- chronostrat[chronostrat$scale %in% relv_scales,];
chronostrat <- chronostrat[order(chronostrat$ma_lb,chronostrat$ma_lb,decreasing = T),];
chronostrat <- chronostrat[order(-chronostrat$ma_lb,chronostrat$ma_lb),];

chronostrat <- subset(chronostrat,chronostrat$ma_lb!=chronostrat$ma_ub);
bin_onsets <- sort(unique(chronostrat$ma_lb),decreasing=T);
bin_ends <- sort(unique(chronostrat$ma_ub),decreasing=T);

not_bin_ends <- bin_ends[!bin_ends %in% bin_onsets];
not_bin_onsets <- bin_onsets[!bin_onsets %in% bin_ends];

n_intr <- nrow(chronostrat);	# total intervals;
if (length(bin_onsets)==length(bin_ends))	{
	bin_start_stop <- cbind(bin_onsets,bin_ends);
	parent_interval <- array("",dim=c(n_intr));
	bin_first <- match(chronostrat$ma_lb,bin_onsets);
	bin_last <- match(chronostrat$ma_ub,bin_ends);
	bin_spans <- 1+bin_last-bin_first;
	nbins <- max(bin_last);
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")
				parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
			parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	if (bs==0 && (regional_scale=="Stage Slice") || (regional_scale=="Time Slice"))	{
		slice <- chronostrat$interval;
		parent_interval <- sapply(slice,accersi_parent_intervals_for_stage_slice);
		} else	{
		parent_interval[(1:n_intr)[parent_interval==""]] <- chronostrat$interval[(1:n_intr)[parent_interval==""]];
		}
	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(parent_interval), .after = match("interval",colnames(chronostrat)));
	} else	{
	bin_boundaries <- sort(unique(c(bin_onsets,bin_ends)),decreasing = T);
	bin_first <- match(chronostrat$ma_lb,bin_boundaries);
	bin_last <- match(chronostrat$ma_ub,bin_boundaries)-1;
	bin_spans <- 1 + bin_last - bin_first;
	parent_interval <- array("",dim=c(n_intr));
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
#	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(chronostrat$interval[chronostrat_subintervals]), .after = match("interval",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")	{
#				chronostrat$parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
				parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
				}
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
#			chronostrat$parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	if (is.null(chronostrat$parent_interval))	{
		chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(parent_interval), .after = match("interval",colnames(chronostrat)));
		} else	{
		chronostrat$parent_interval <- parent_interval;
		}
	}
return(chronostrat);
}

accersi_hierarchical_timescale_old <- function(chronostrat_units,time_scale,regional_scale="International",ma_fuzz=0)	{
chronostrat <- time_scale[as.character(time_scale$interval) %in% chronostrat_units,];
oldest <- max(time_scale$ma_lb[as.character(time_scale$interval) %in% chronostrat_units])+ma_fuzz;
youngest <- min(time_scale$ma_ub[as.character(time_scale$interval) %in% chronostrat_units])-ma_fuzz;
relv_scales <- unique(c(chronostrat$scale,regional_scale));

xxx <- subset(time_scale,time_scale$ma_lb<=oldest);
chronostrat <- subset(xxx,xxx$ma_ub>=youngest);
chronostrat <- subset(chronostrat,chronostrat$scale==regional_scale);
chronostrat <- chronostrat[chronostrat$scale %in% relv_scales,];
chronostrat <- chronostrat[order(chronostrat$ma_lb,chronostrat$ma_lb,decreasing = T),];
chronostrat <- chronostrat[order(-chronostrat$ma_lb,chronostrat$ma_lb),];

chronostrat <- subset(chronostrat,chronostrat$ma_lb!=chronostrat$ma_ub);
bin_onsets <- sort(unique(chronostrat$ma_lb),decreasing=T);
bin_ends <- sort(unique(chronostrat$ma_ub),decreasing=T);

n_intr <- nrow(chronostrat);	# total intervals;
if (length(bin_onsets)==length(bin_ends))	{
	parent_interval <- array("",dim=c(n_intr));
	bin_first <- match(chronostrat$ma_lb,bin_onsets);
	bin_last <- match(chronostrat$ma_ub,bin_ends);
	bin_spans <- 1+bin_last-bin_first;
	nbins <- max(bin_last);
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")
				parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
			parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	if (bs==0 && (regional_scale=="Stage Slice") || (regional_scale=="Time Slice"))	{
		slice <- chronostrat$interval;
		parent_interval <- sapply(slice,accersi_parent_intervals_for_stage_slice);
		} else	{
		parent_interval[(1:n_intr)[parent_interval==""]] <- chronostrat$interval[(1:n_intr)[parent_interval==""]];
		}
	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(parent_interval), .after = match("interval",colnames(chronostrat)));
	} else	{
	bin_boundaries <- sort(unique(c(bin_onsets,bin_ends)),decreasing = T);
	bin_first <- match(chronostrat$ma_lb,bin_boundaries);
	bin_last <- match(chronostrat$ma_ub,bin_boundaries)-1;
	bin_spans <- 1 + bin_last - bin_first;
	parent_interval <- array("",dim=c(n_intr));
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
#	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(chronostrat$interval[chronostrat_subintervals]), .after = match("interval",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")	{
#				chronostrat$parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
				parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
				}
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
#			chronostrat$parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(parent_interval), .after = match("interval",colnames(chronostrat)));
	}
return(chronostrat);
}

accersi_parent_intervals_for_stage_slice <- function(slice)	{
j <- strsplit(slice,split="",fixed=TRUE)[[1]];
jl <- length(j);
if (sum(c("a","b","c","d") %in% j[jl])==1)	jl <- jl-1;
jl <- jl-1;
slice <- paste(j[1:jl],collapse="");
}

accersi_coarsened_time_scale <- function(hierarchical_chronostrat,interval_thesaurus)	{
final_intervals <- unique(interval_thesaurus$new_interval);
f_i <- length(final_intervals);
coarsened_time_scale <- hierarchical_chronostrat[1:f_i,];

for (fi in 1:f_i)	{
	if (sum(interval_thesaurus$new_interval==final_intervals[fi])>1)	{
		old_ints <- match(interval_thesaurus$orig_interval[interval_thesaurus$new_interval==final_intervals[fi]],hierarchical_chronostrat$interval);
		coarsened_time_scale$ma_lb[fi] <- max(hierarchical_chronostrat$ma_lb[old_ints]);
		coarsened_time_scale$ma_ub[fi] <- min(hierarchical_chronostrat$ma_ub[old_ints]);
		coarsened_time_scale$interval[fi] <- coarsened_time_scale$interval_alt[fi] <- coarsened_time_scale$st[fi] <- final_intervals[fi];
		coarsened_time_scale$parent_interval[fi] <- hierarchical_chronostrat$parent_interval[old_ints[1]];
		coarsened_time_scale$scale[fi] <- hierarchical_chronostrat$scale[old_ints[1]];
		coarsened_time_scale$color[fi] <- hierarchical_chronostrat$color[old_ints[1]];
		coarsened_time_scale$record_no[fi] <- hierarchical_chronostrat$record_no[old_ints[1]];
		coarsened_time_scale$created[fi] <- hierarchical_chronostrat$created[old_ints[1]];
		coarsened_time_scale$modified[fi] <- hierarchical_chronostrat$modified[old_ints[1]];
		} else	{
		coarsened_time_scale[fi,] <- hierarchical_chronostrat[match(final_intervals[fi],hierarchical_chronostrat$interval),];
		}
	}
coarsened_time_scale <- subset(coarsened_time_scale,!is.na(coarsened_time_scale$interval));
unaltered_bins <- interval_thesaurus$bin_old[interval_thesaurus$bin>max(coarsened_time_scale$bin_first)];
coarsened_time_scale <- rbind(coarsened_time_scale,hierarchical_chronostrat[hierarchical_chronostrat$bin_first %in% unaltered_bins,]);

coarsened_time_scale$bin_first <- coarsened_time_scale$bin_last <- 1:nrow(coarsened_time_scale);
return(coarsened_time_scale);
}

# use ages to reassign collection to an interval on a particular times scale
reassign_intervals_to_uniform_scale <- function(ma,uniform_time_scale,onset,round_extremes=T)	{
# age (in millions of years)
# uniform_time_scale: time scale that is one time line, so that no interval belongs to another interval
# onset: if T, then this ma is the lower bound; if F, then this ma is the upper bound
# round_extremes: if T, then lower bound predating the first interval goes in interval 1;
#	an upper bound postdating the last interval goes in the last
if (onset)	{
	if (ma <= max(uniform_time_scale$ma_lb) && ma > min(uniform_time_scale$ma_ub))	{
		return(uniform_time_scale$interval[sum(ma<=uniform_time_scale$ma_lb)])
		} else	{
		if (round_extremes)	{
			return(uniform_time_scale$interval[1]);
			} else	{
			return("");
			}
		}
	} else	{
	if (ma >= min(uniform_time_scale$ma_ub) && ma < max(uniform_time_scale$ma_lb))	{
		return(uniform_time_scale$interval[sum(ma<=uniform_time_scale$ma_lb)]);
	} else	{
		if (round_extremes)	{
			return(uniform_time_scale$interval[nrow(uniform_time_scale)]);
			} else	{
			return("");
			}
		}
	}
}

#update_paleodb_collection_ages
# get temporal overlap between two spans (ranges, intervals, etc.)
do_two_ranges_overlap <- function(lb_a,ub_a,lb_b,ub_b)	{
lb1 <- max(abs(lb_a));
ub1 <- min(abs(ub_a));
lb2 <- max(abs(lb_b));
ub2 <- min(abs(ub_b));
if ((lb1<=lb2 && lb1>ub2) || (lb2<=lb1 && lb2>ub1))	{
	return(T);
	} else	{
	return(F);
	}
}

accersi_temporal_overlap <- function(lb1,ub1,lb2,ub2)	{
lb1 <- round(lb1,3);
lb2 <- round(lb2,3);
ub1 <- round(ub1,3);
ub2 <- round(ub2,3);
overlap <- data.frame(ma_lb=as.numeric(0),ma_ub=as.numeric(0));
if ((lb1 <= lb2 && lb1 > ub2) || (lb2 <= lb1 && lb2 > ub1)) {
#if (fa1 > la2 || fa2 > la1)	{
	overlap$ma_lb <- min(lb1,lb2);	# uppermost point of lower bound;
	overlap$ma_ub <- max(ub1,ub2);	# lowermost point of upper bound;
	}
return (overlap);
}

accersi_temporal_overlap_multiple_cases <- function(lb,ub)	{
overlaps <- comps <- c();
for (ll in 1:(length(lb)-1))	{
	for (uu in (ll+1):length(lb))	{
		overlaps <- rbind(overlaps,accersi_temporal_overlap(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]));
		comps <- c(comps,paste(ll,"v",uu,sep=""));
		}
	}
rownames(overlaps) <- comps;
return (overlaps);
}

accersi_temporal_overlap_matrix <- function(lb,ub)	{
overlap_matrix <- array(0,dim=c(length(lb),length(lb)));
overlaps <- comps <- c();
for (ll in 1:(length(lb)-1))	{
	for (uu in (ll+1):length(lb))	{
		overlap <- accersi_temporal_overlap(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]);
		overlap_matrix[ll,uu] <- overlap_matrix[uu,ll] <- overlap$ma_lb-overlap$ma_ub;
#		comps <- c(comps,paste(ll,"v",uu,sep=""));
		}
	}
return (overlap_matrix);
}

identify_temporal_overlap_multiple_cases <- function(lb,ub)	{
#overlaps <- c();
overlapping <- data.frame(range_1=as.numeric(),range_2=as.numeric());
if (length(lb)==1)	{
	return(overlapping)
	} else	{
	for (ll in 1:(length(lb)-1))	{
		for (uu in (ll+1):length(lb))	{
			overlaps <- accersi_temporal_overlap(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]);
			if (overlaps$ma_lb!=0)	{
				dummy_frame <- data.frame(range_1=as.numeric(ll),range_2=as.numeric(uu));
				overlapping <- rbind(overlapping,dummy_frame);
				}
			}
		}
	#rownames(overlaps) <- comps;
	return (overlapping);
	}
}

# recalculate range given disjunct ranges
accersi_minimum_extension <- function(lb1,ub1,lb2,ub2)	{
lb1 <- round(lb1,3);
lb2 <- round(lb2,3);
ub1 <- round(ub1,3);
ub2 <- round(ub2,3);
range <- data.frame(ma_lb=as.numeric(0),ma_ub=as.numeric(0));
if (lb1 < ub2)	{
	range$ma_lb <- lb2;
	range$ma_ub <- ub1;
	} else	{
	range$ma_lb <- lb1;
	range$ma_ub <- ub2;
	}
return (range);
}

# recalculate range given disjunct ranges
accersi_minimum_extension_multiple_cases <- function(lb,ub)	{
ranges <- c();
for (ll in 1:(length(lb)-1))	{
	for (uu in (ll+1):length(lb))	{
		ranges <- rbind(ranges,accersi_minimum_extension(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]));
		}
	}
return (ranges);
}

accersi_minimum_joint_ranges <- function(lb,ub)	{
joint_ranges <- data.frame(ma_lb=as.numeric(),ma_ub=as.numeric());
for (ll in 1:(length(lb)-1))	{
	for (uu in (ll+1):length(lb))	{
		joint_ranges <- rbind(joint_ranges,accersi_temporal_overlap(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]));
		if (joint_ranges[1,1]==0)	{
			redo <- nrow(joint_ranges)
			joint_ranges[redo,] <- accersi_minimum_extension(lb1=lb[ll],ub1=ub[ll],lb2=lb[uu],ub2=ub[uu]);
			}
		}
	}
return(joint_ranges);
}

# replace reported stratigraphic intervals with appropriate intervals from one time scale (usually internationally)
completely_rebin_collections_with_uniform_time_scale <- function(collections,uniform_time_scale,add_bins=F)	{
uniform_time_scale$ma_lb <- 0.001*round(uniform_time_scale$ma_lb/0.001,0);
uniform_time_scale$ma_ub <- 0.001*round(uniform_time_scale$ma_ub/0.001,0);
if (!is.null(collections$ma_lb))	{
	ma_lb_col <- match("ma_lb",colnames(collections));
	ma_ub_col <- match("ma_ub",colnames(collections));
	} else	{
	ma_lb_col <- match("max_ma",colnames(collections));
	ma_ub_col <- match("min_ma",colnames(collections));
	}
collections[,ma_lb_col] <- 0.001*round(collections[,ma_lb_col]/0.001,0);
collections[,ma_ub_col] <- 0.001*round(collections[,ma_ub_col]/0.001,0);
if (!is.null(collections$interval_lb))	{
	int_lb_col <- match("interval_lb",colnames(collections));
	int_ub_col <- match("interval_ub",colnames(collections));
	} else	{
	int_lb_col <- match("early_interval",colnames(collections));
	int_ub_col <- match("late_interval",colnames(collections));
	}
age <- collections[,ma_lb_col];
#for (a in 1:length(age))	{
#	collections[a,int_lb_col] <- rebin_collection_with_time_scale(age=age[a],onset_or_end="onset",fine_time_scale=uniform_time_scale);
#	}
collections[,int_lb_col] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=uniform_time_scale);

age <- collections[,ma_ub_col];
#for (a in 1:length(age))	{
#	collections[a,int_ub_col] <- rebin_collection_with_time_scale(age=age[a],onset_or_end="end",fine_time_scale=uniform_time_scale);
#	}
collections[,int_ub_col] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=uniform_time_scale);
if (!is.null(collections$bin_lb) || add_bins)	{
	collections$bin_lb <- match(collections[,int_lb_col],uniform_time_scale$interval);
	collections$bin_ub <- match(collections[,int_ub_col],uniform_time_scale$interval);
	}
return(collections);
}

rebin_collection_with_time_scale <- function(age,onset_or_end,fine_time_scale)	{
# age: age in millions of years
# onset_or_end: "onset" for lower bound, "end" for upper bound
# fine_time_scale: dataframe where:
#	fine_time_scale$interval gives interval name
#	fine_time_scale$ma_lb gives interval onset (given Gradstein et al. 2012)
#	fine_time_scale$ma_ub gives interval end (given Gradstein et al. 2012)
#	NOTE: fine_time_scale cannot include subintervals of other intervals in that time scale:
#	  e.g., just Cambrian, Ordovician, Silurian or Sandbian, Katian, Rhuddanian, etc.
age <- round(abs(age),3);
if (onset_or_end=="onset" || onset_or_end=="Onset")	{
	return(as.character(fine_time_scale$interval[max(1,sum(age<=round(fine_time_scale$ma_lb,3)))]));
	} else	{
	return(as.character(fine_time_scale$interval[max(1,sum(age<round(fine_time_scale$ma_lb,3)))]));
	}
}

# editted 2020-03-05
redate_collections_with_direct_dates <- function(collections,finest_chronostrat,temporal_precision=0.05)	{
ncolls <- nrow(collections);
collections$direct_ma <- as.numeric(collections$direct_ma);
collections$direct_ma_error <- as.numeric(collections$direct_ma_error);
collections$direct_ma <- expello_na_from_vector(collections$direct_ma,0);
collections$direct_ma_error <- expello_na_from_vector(collections$direct_ma_error,0);
beakerheads <- (1:ncolls)[collections$direct_ma>0];
if (length(beakerheads)>0)	{
	if (is.null(collections$interval_lb))	{
		age <- collections$max_ma;
		collections$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",fine_time_scale=finest_chronostrat);
		age <- collections$min_ma;
		collections$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",fine_time_scale=finest_chronostrat);
		} else	{
		collections$interval_lb <- as.character(collections$interval_lb);	# kluge
		collections$interval_ub <- as.character(collections$interval_ub);	# kluge
		}
	if (!is.null(collections$ma_lb))	{
		collections$direct_ma_error[beakerheads][collections$direct_ma_error[beakerheads]==0] <- temporal_precision;
		collections$ma_lb[beakerheads] <- collections$direct_ma[beakerheads]+collections$direct_ma_error[beakerheads];
		collections$ma_ub[beakerheads] <- collections$direct_ma[beakerheads]-collections$direct_ma_error[beakerheads];
		age <- collections$ma_lb[beakerheads];
		collections$interval_lb[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
		age <- collections$ma_ub[beakerheads];
	#	dummy <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
		collections$interval_ub[beakerheads] <- as.character(collections$interval_ub[beakerheads])
	#	for (dd in 1:length(dummy))	{
	#		print(as.character(collections$interval_ub[beakerheads[dd]]))
	#		}
#		collections$interval_ub[beakerheads] <- rebin_collection_with_time_scale(age=collections$ma_ub[beakerheads],onset_or_end="end",fine_time_scale=finest_chronostrat);
		} else	{
		collections$max_ma[beakerheads] <- collections$direct_ma[beakerheads]+collections$direct_ma_error[beakerheads];
		collections$min_ma[beakerheads] <- collections$direct_ma[beakerheads]-collections$direct_ma_error[beakerheads];
		age <- collections$max_ma[beakerheads];
		collections$early_interval[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
		age <- collections$min_ma[beakerheads];
		collections$late_interval[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
#		collections$late_interval[beakerheads] <- rebin_collection_with_time_scale(age=collections$ma_lb[beakerheads],onset_or_end="end",fine_time_scale=finest_chronostrat);
		}
	}
return(collections);
}

lump_zone_ages_by_higher_taxon <- function(relevant_zones)	{
higher_taxa <- unique(relevant_zones$zone_type);
higher_taxon_ranges <- data.frame(taxon=as.character(higher_taxa),
								  ma_lb=as.numeric(rep(0,length(higher_taxa))),
								  ma_ub=as.numeric(rep(0,length(higher_taxa))),
								  interval_lb=as.numeric(rep("",length(higher_taxa))),
								  interval_ub=as.numeric(rep("",length(higher_taxa))),
								  stringsAsFactors = F);
for (ht in 1:length(higher_taxa))	{
	this_zone_taxon <- subset(relevant_zones,relevant_zones$zone_type==higher_taxa[ht]);
	higher_taxon_ranges$ma_lb[ht] <- max(this_zone_taxon$ma_lb);
	higher_taxon_ranges$ma_ub[ht] <- min(this_zone_taxon$ma_ub);
	higher_taxon_ranges$interval_lb[ht] <- this_zone_taxon$interval_lb[match(higher_taxon_ranges$ma_lb[ht],this_zone_taxon$ma_lb)];
	higher_taxon_ranges$interval_ub[ht] <- this_zone_taxon$interval_ub[match(higher_taxon_ranges$ma_ub[ht],this_zone_taxon$ma_ub)];
	}
return(higher_taxon_ranges);
}

accersi_age_range_given_multiple_zones <- function(relevant_zones)	{
higher_taxon_ranges <- lump_zone_ages_by_higher_taxon(relevant_zones);
ttl_hi_taxa <- length(higher_taxon_ranges$taxon);

if (ttl_hi_taxa>1)	{
	zone_comps <- (((ttl_hi_taxa^2))-ttl_hi_taxa)/2;
	minimum_zone_overlap <- data.frame(ma_lb=rep(0,zone_comps),ma_ub=rep(0,zone_comps));
	zc <- 0;
	for (ht1 in 1:(ttl_hi_taxa-1))	{
		for (ht2 in (ht1+1):ttl_hi_taxa)	{
			zc <- zc+1;
			minimum_zone_overlap[zc,] <- accersi_temporal_overlap(lb1=as.numeric(higher_taxon_ranges$ma_lb[ht1]),ub1=as.numeric(higher_taxon_ranges$ma_ub[ht1]),lb2=as.numeric(higher_taxon_ranges$ma_lb[ht2]),ub2=as.numeric(higher_taxon_ranges$ma_ub[ht2]));
			if (minimum_zone_overlap$ma_lb[zc]==0)
				minimum_zone_overlap$ma_lb[zc] <- max(higher_taxon_ranges$ma_lb[c(ht1,ht2)]);
			if (minimum_zone_overlap$ma_ub[zc]==0)
				minimum_zone_overlap$ma_ub[zc] <- min(higher_taxon_ranges$ma_ub[c(ht1,ht2)]);
			}
		}
	} else	{
	minimum_zone_overlap <- data.frame(ma_lb=max(relevant_zones$ma_lb),ma_ub=min(relevant_zones$ma_ub));
	}

zone_ranges <- data.frame(ma_lb=min(minimum_zone_overlap$ma_lb),ma_ub=max(minimum_zone_overlap$ma_ub),
						  interval_lb=relevant_zones$interval_lb[match(min(minimum_zone_overlap$ma_lb),relevant_zones$ma_lb)],interval_ub=relevant_zones$interval_ub[match(max(minimum_zone_overlap$ma_ub),relevant_zones$ma_ub)],
						  stringsAsFactors = F);
too_old_zones <- subset(relevant_zones,relevant_zones$ma_ub>=zone_ranges$ma_lb);
while (nrow(too_old_zones)>0)	{
	new_ma_lb <- min(too_old_zones$ma_lb);
	new_interval_lb <- too_old_zones$interval_lb[match(new_ma_lb,too_old_zones$ma_lb)];
	zone_ranges$ma_lb <- new_ma_lb;
	zone_ranges$interval_lb <- new_interval_lb;
#	too_old_zones <- subset(relevant_zones,relevant_zones$ma_ub>=zone_ranges$ma_lb);
	too_old_zones <- subset(too_old_zones,too_old_zones$ma_ub>new_ma_lb);
	# remove zone(s) to prevent infinite loops on very short zones
	}
too_young_zones <- subset(relevant_zones,relevant_zones$ma_lb<=zone_ranges$ma_ub);
while (nrow(too_young_zones)>0)	{
	new_ma_ub <- max(too_young_zones$ma_ub);
	new_interval_ub <- too_young_zones$interval_ub[match(new_ma_ub,too_young_zones$ma_ub)];
	zone_ranges$ma_ub <- new_ma_ub;
	zone_ranges$interval_ub <- new_interval_ub;
#	zone_ranges$ma_ub
	too_young_zones <- subset(too_young_zones,too_young_zones$ma_lb<new_ma_ub);
	# remove zone(s) to prevent infinite loops on very short zones
	}
return(zone_ranges);
}

irradiate_collections <- function(paleodb_collections,fossilworks_collections)	{
#fossilworks_collections$direct_ma[fossilworks_collections$collection_no %in% paleodb_collections$collection_no];
ncoll1 <- nrow(paleodb_collections);
relv_fossilworks <- subset(fossilworks_collections,fossilworks_collections$direct_ma>0);
relv_fossilworks <- subset(relv_fossilworks,!is.na(as.numeric(relv_fossilworks$direct_ma)));
ncoll2 <- nrow(relv_fossilworks);
godzilla <- (1:ncoll1)[paleodb_collections$collection_no %in% relv_fossilworks$collection_no];
paleodb_collections$direct_ma <- rep(0,ncoll1);
paleodb_collections$direct_ma[godzilla] <- as.numeric(relv_fossilworks$direct_ma[match(paleodb_collections$collection_no[godzilla],relv_fossilworks$collection_no)]);
paleodb_collections$direct_ma_error <- rep(0,ncoll1);
paleodb_collections$direct_ma_error[godzilla] <- as.numeric(relv_fossilworks$direct_ma_error[match(paleodb_collections$collection_no[godzilla],relv_fossilworks$collection_no)]);
paleodb_collections$direct_ma_method <- rep("",ncoll1);
paleodb_collections$direct_ma_method[godzilla] <- as.character(relv_fossilworks$direct_ma_method[match(paleodb_collections$collection_no[godzilla],relv_fossilworks$collection_no)]);

#paleodb_collections$direct_ma[match(fossilworks_collections$collection_no[fossilworks_collections$direct_ma>0][fossilworks_collections$collection_no[fossilworks_collections$direct_ma>0] %in% paleodb_collections$collection_no],paleodb_collections$collection_no)]
#paleodb_collections$direct_ma <- as.numeric(fossilworks_collections$direct_ma[match(paleodb_collections$collection_no,fossilworks_collections$collection_no)]);
#paleodb_collections$direct_ma_error <- as.numeric(fossilworks_collections$direct_ma_error[match(paleodb_collections$collection_no,fossilworks_collections$collection_no)])/2;
#paleodb_collections$direct_ma_method <- as.character(fossilworks_collections$direct_ma_method[match(paleodb_collections$collection_no,fossilworks_collections$collection_no)]);
#paleodb_collections$direct_ma <- expello_na_from_vector(paleodb_collections$direct_ma,replacement = 0);
#paleodb_collections$direct_ma_error <- expello_na_from_vector(paleodb_collections$direct_ma_error,replacement = 0);

#relevant_fossilworks <- fossilworks_collections[fossilworks_collections$collection_no %in% paleodb_collections$collection_no,]
#ncoll_fw <- nrow(relevant_fossilworks);
#ncoll_pb <- nrow(paleodb_collections);
#fwc1 <- (1:ncoll_fw)[relevant_fossilworks$max_ma_method %in% c("Ar/Ar","14C","14C (calibrated)","fission track","K-Ar","Lu-Hf","Rb-Sr","Sr isotope","U/Pb","U/Th")];
##paste(relevant_fossilworks$collection_no[fwc1],collapse=",");
#fwc2 <- (1:ncoll_fw)[relevant_fossilworks$min_ma_method %in% c("Ar/Ar","14C","14C (calibrated)","fission track","K-Ar","Lu-Hf","Rb-Sr","Sr isotope","U/Pb","U/Th")];
#fwc <- fwc1[fwc1 %in% fwc2];
#pdwc <- match(relevant_fossilworks$collection_no[fwc],paleodb_collections$collection_no);
#paleodb_collections$direct_ma[pdwc] <- (as.numeric(relevant_fossilworks$max_ma[fwc])+as.numeric(relevant_fossilworks$min_ma[fwc]))/2;
#paleodb_collections$direct_ma_error[pdwc] <- (as.numeric(relevant_fossilworks$max_ma[fwc])-as.numeric(relevant_fossilworks$min_ma[fwc]))/2;
#paleodb_collections$direct_ma_method[pdwc] <- relevant_fossilworks$max_ma_method[fwc];
#
#fwc1 <- fwc1[!fwc1 %in% fwc];
#pdwc1 <- match(relevant_fossilworks$collection_no[fwc1],paleodb_collections$collection_no);
#paleodb_collections$direct_ma[pdwc1] <- as.numeric(relevant_fossilworks$max_ma[fwc1]);
#paleodb_collections$direct_ma_error[pdwc1] <- 0;
#paleodb_collections$direct_ma_method[pdwc1] <- relevant_fossilworks$max_ma_method[fwc1];
#fwc2 <- fwc2[!fwc2 %in% fwc];
#pdwc2 <- match(relevant_fossilworks$collection_no[fwc2],paleodb_collections$collection_no);
#paleodb_collections$direct_ma[pdwc2] <- as.numeric(relevant_fossilworks$min_ma[fwc2]);
#paleodb_collections$direct_ma_error[pdwc2] <- 0;
#paleodb_collections$direct_ma_method[pdwc2] <- relevant_fossilworks$min_ma_method[fwc2];
return(paleodb_collections);
#direct_dates <- data.frame(direct_ma=as.numeric(direct_ma),direct_ma_error=as.numeric(direct_ma_error),direct_ma_method=as.character(direct_ma_method));
#direct_dates <- expello_na_from_matrix(direct_dates,"");
}

redate_paleodb_collections_with_direct_dates <- function(paleodb_collections,finest_chronostrat)	{
if (!is.numeric(paleodb_collections$direct_ma))	{
	paleodb_collections$direct_ma <- as.numeric(paleodb_collections$direct_ma);
	paleodb_collections$direct_ma_error <- as.numeric(paleodb_collections$direct_ma_error);
	paleodb_collections$direct_ma[is.na(paleodb_collections$direct_ma)] <- 0;
	paleodb_collections$direct_ma_error[is.na(paleodb_collections$direct_ma_error)] <- 0;
	}

ncoll <- nrow(paleodb_collections);
dated_collections <- (1:ncoll)[paleodb_collections$direct_ma!=0];
if (length(dated_collections)>0)	{
	redate_these <- dated_collections[paleodb_collections$max_ma[dated_collections]>(paleodb_collections$direct_ma[dated_collections]+(paleodb_collections$direct_ma_error[dated_collections]/2))];
#cbind(paleodb_collections$max_ma[dated_collections],(paleodb_collections$direct_ma[dated_collections]+paleodb_collections$direct_ma_error[dated_collections]/2));
	paleodb_collections$max_ma[redate_these] <- paleodb_collections$direct_ma[redate_these]+(paleodb_collections$direct_ma_error[redate_these]/2);
	age <- paleodb_collections$max_ma[redate_these];
	paleodb_collections$early_interval[redate_these] <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
	if (!is.null(paleodb_collections$ma_lb))	{
		redate_these <- dated_collections[paleodb_collections$ma_lb[dated_collections]>(paleodb_collections$direct_ma[dated_collections]+(paleodb_collections$direct_ma_error[dated_collections]/2))];
		paleodb_collections$ma_lb[redate_these] <- paleodb_collections$direct_ma[dated_collections]+(paleodb_collections$direct_ma_error[dated_collections]/2);
		age <- paleodb_collections$ma_lb[redate_these];
		paleodb_collections$interval_lb[redate_these] <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
		}
	redate_these <- dated_collections[paleodb_collections$min_ma[dated_collections]<(paleodb_collections$direct_ma[dated_collections]-(paleodb_collections$direct_ma_error[dated_collections]/2))];
	redate_these <- sort(unique(c(redate_these,dated_collections[paleodb_collections$min_ma[dated_collections]>paleodb_collections$max_ma[dated_collections]])));
	paleodb_collections$min_ma[redate_these] <- paleodb_collections$direct_ma[redate_these]-(paleodb_collections$direct_ma_error[redate_these]/2);
	age <- paleodb_collections$min_ma[redate_these];
	paleodb_collections$late_interval[redate_these] <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);
	if (!is.null(paleodb_collections$ma_ub))	{
		redate_these <- dated_collections[paleodb_collections$ma_ub[dated_collections]<(paleodb_collections$direct_ma[dated_collections]-(paleodb_collections$direct_ma_error[dated_collections]/2))];
		redate_these <- sort(unique(c(redate_these,dated_collections[paleodb_collections$ma_ub[dated_collections]>paleodb_collections$ma_lb[dated_collections]])));
		paleodb_collections$ma_lb[redate_these] <- paleodb_collections$direct_ma[dated_collections]-(paleodb_collections$direct_ma_error[dated_collections]/2);
		age <- paleodb_collections$ma_ub[redate_these];
		paleodb_collections$interval_ub[redate_these] <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);
		}
	}
return(paleodb_collections);
}