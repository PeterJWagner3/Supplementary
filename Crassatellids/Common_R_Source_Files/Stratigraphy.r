# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal

library(dplyr);
library(lettercase);	#devtools::install_github('decisionpatterns/lettercase')
library(paleobioDB);    #install.packages("paleobioDB", dependencies=TRUE);
library(prodlim);		#install.packages("prodlim", dependencies=TRUE);
library(Rcpp);			#devtools::install_github("https://github.com/RcppCore/Rcpp/")
library(rvest);
library(stringr);       #install.packages("stringr", dependencies=TRUE)
library(xml2);			#install.packages("xml2", dependencies=TRUE);

hell_no <- F;
MAXNO <- 1.797693e+308;

#### Routines to get information about rock units from external database ####
summarize_relevant_rock_units <- function(all_rock_nos,rock_database)	{
rock_database_matches <- match(all_rock_nos,rock_database$rock_no);
relv_rock_summary <- data.frame(rock_unit=as.character(rock_database$full_name[rock_database_matches]),rock_no=as.numeric(rock_database$rock_no[rock_database_matches]),formation_no=as.numeric(rock_database$formation_no[rock_database_matches]),status=as.numeric(rep(0,length=length(all_rock_nos))),stringsAsFactors = hell_no);
relv_rock_summary$status[relv_rock_summary$rock_no==relv_rock_summary$formation_no] <- 1;
relv_rock_summary$status[relv_rock_summary$rock_unit==""] <- 2;
relv_rock_summary$rock_unit[relv_rock_summary$rock_unit==""] <- rock_database$group[match(relv_rock_summary$rock_no[relv_rock_summary$status==2],rock_database$rock_no)];

relv_rock_summary <- relv_rock_summary[order(relv_rock_summary$formation_no,relv_rock_summary$status,relv_rock_summary$rock_no),];
return(relv_rock_summary);
}

# extract just the rock-unit numbers from PaleoDB collections that have been matched to an external database
accersi_rock_numbers_from_paleodb_collections <- function(paleodb_collections)	{
return(sort(unique(c(paleodb_collections$rock_no_sr[paleodb_collections$rock_no_sr>0],paleodb_collections$rock2_no_sr[paleodb_collections$rock2_no_sr>0]))))
}

# this is designed to deal with PaleoDB collections and an external rock database
finds_per_rock_unit_per_taxon <- function(taxon_list,relv_rock_summary,paleodb_finds,paleodb_collections,rock_database)	{
notu <- length(taxon_list);
ttl_rocks <- nrow(relv_rock_summary);
finds_per_rock <- array(0,dim=c(notu,ttl_rocks));
colnames(finds_per_rock) <- relv_rock_summary$rock_unit;
rownames(finds_per_rock) <- taxon_list;
for (tl in 1:notu)	{
#	tl <- tl+1;
	find_nos <- which(paleodb_finds==taxon_list[tl],arr.ind = T)[,1];
	rock_find_nos <- sort(paleodb_collections$rock_no_sr[match(sort(unique(paleodb_finds$collection_no[find_nos])),paleodb_collections$collection_no)]);
	rock_find_nos <- rock_find_nos[rock_find_nos>0];
	rock_nos <- unique(rock_find_nos);
	if (length(rock_nos)>0)	{
		xxx <- hist(rock_find_nos,breaks=c(1:(max(rock_find_nos)+1)),plot=F)$counts;
		yyy <- relv_rock_summary$rock_unit[match(rock_nos,relv_rock_summary$rock_no)];
		zzz <- match(yyy,colnames(finds_per_rock));
		finds_per_rock[tl,yyy] <- xxx[xxx>0];
		}
	}
return(finds_per_rock);
}

finds_per_rock_unit_per_taxon_old <- function(taxon_list,relv_rock_summary,paleodb_finds,paleodb_collections,rock_database)	{
notu <- length(taxon_list);
ttl_rocks <- nrow(relv_rock_summary);
finds_per_rock <- array(0,dim=c(notu,ttl_rocks));
colnames(finds_per_rock) <- relv_rock_summary$rock_unit;
rownames(finds_per_rock) <- taxon_list;
for (tl in 1:notu)	{
#	tl <- tl+1;
	find_nos <- which(paleodb_finds==taxon_list[tl],arr.ind = T)[,1];
	rock_find_nos <- sort(paleodb_collections$rock_no_sr[match(sort(unique(paleodb_finds$collection_no[find_nos])),paleodb_collections$collection_no)]);
	rock_find_nos <- rock_find_nos[rock_find_nos>0];
	rock_nos <- unique(rock_find_nos);
	if (length(rock_nos)>0)	{
		xxx <- hist(rock_find_nos,breaks=c(1:(max(rock_find_nos)+1)),plot=F)$counts;
		yyy <- rock_database$full_name[match(rock_nos,rock_database$rock_no)];
		zzz <- match(yyy,colnames(finds_per_rock));
		finds_per_rock[tl,yyy] <- xxx[xxx>0];
		}
	}
return(finds_per_rock);
}

#### Basic Superposition Fun ####
constrain_dates_given_superposition <- function(ages,young_low=F)	{
# ages$lb: oldest possible age
# ages$ub: youngest possible age
if (young_low)
	ages <- -1*ages;
ttl_comps <- nrow(ages);
for (tc in 1:(ttl_comps-1))	{
	tcc <- tc+1;
	ages$lb[tcc:ttl_comps][ages$lb[tcc:ttl_comps]<ages$lb[tc]] <- ages$lb[tc];
	}
for (tc in ttl_comps:2)	{
	tcc <- tc-1;
	ages$ub[1:tcc][ages$ub[1:tcc]>ages$ub[tc]] <- ages$ub[tc];
	}
if (young_low)
	ages <- -1*ages;
return(ages);
}

# routine to find contradictory superposition statements
#	given a matrix of statements about species being found
#	above one another in different sections
find_circuit <- function(case_matrix,start=1,open=T)	{
# case_matrix$x1: step 1
# case_matrix$x2: step 2
case_matrix <- unique(case_matrix);
while (nrow(case_matrix)>nrow(case_matrix[case_matrix[,1] %in% case_matrix[,2],]) || nrow(case_matrix)>nrow(case_matrix[case_matrix[,2] %in% case_matrix[,1],]))	{
	case_matrix <- case_matrix[case_matrix[,1] %in% case_matrix[,2],];
	case_matrix <- case_matrix[case_matrix[,2] %in% case_matrix[,1],];
#	print(case_matrix)
	}
if (start>nrow(case_matrix))	start <- nrow(case_matrix);
links <- c(case_matrix$x1[start],case_matrix$x2[start]);

ll <- 2;
while (open)	{
	nxt_link <- match(links[ll],case_matrix$x1)
	if (!is.na(nxt_link))	{
		links <- c(links,case_matrix$x2[nxt_link]);
		ll <- length(links);
		if (links[ll]==links[1])	open <- F;
		} else	{
		open <- F;
		links <- c();
		}
	}
return(links);
}

#this_file <- paste(study_span,"_",study_taxon,"_Find_Order_Initial.csv",sep="");
# occurrence_order <- read.csv(this_file,header=T); rownames(occurrence_order) <- occurrence_order$X; occurrence_order$X <- NULL;
# go through occurrence_order to find all taxa above or below each taxon.
# look out for contradictions along the way
# old_occurrence_order <- occurrence_order
accersi_taxa_below_and_above_taxon <- function(occurrence_order)	{
# occurrence_order: taxon by taxon matrix.
#	occurrence_order[i,j]=1: i first appears above j;
#	occurrence_order[i,j]=1: i first appears below j;
#taxa_below <- taxa_above <- list();
nspec <- nrow(occurrence_order);
spec_below <- occurrence_order;
spec_below[spec_below!=0] <- 0;
spec_above <- spec_below;
rownames(spec_below) <- rownames(spec_above) <- colnames(spec_above) <- colnames(spec_below) <- rownames(occurrence_order);
flag <- c();
print("First Pass");
progress <- round(nspec*(1:100)/100,0);
statements <- sum(occurrence_order!=0);
for (sp in 1:nspec)	{
	if (sp %in% progress)	{
		prop_done <- 100*match(sp,progress)/length(progress);
		if (prop_done<10)	{
			console_update <- paste("0",prop_done,"% done",sep="");
			} else {
			console_update <- paste(prop_done,"% done",sep="");
			}
		cat('\b\b\b\b\b\b\b\b\b',console_update);
		flush.console()
		}
	goblin_town <- (1:nspec)[occurrence_order[sp,]==1];
	aerie <- (1:nspec)[occurrence_order[sp,]==-1];
	lokis <- c(sp,goblin_town[goblin_town %in% aerie]);
	goblin_town <- goblin_town[!goblin_town %in% lokis];
	aerie <- aerie[!aerie %in% lokis];
	if (length(lokis)>1)		flag <- rbind(flag,cbind(rep(sp,length(lokis)-1),lokis[2:length(lokis)]));
	if (length(goblin_town)>0)	spec_below[sp,1:length(goblin_town)] <- goblin_town;
	if (length(aerie)>0)		spec_above[sp,1:length(aerie)] <- aerie;
	}
old_statements <- 0;
statements <- sum(spec_below>0)+sum(spec_above>0);
ttl_runs <- 0;
test1 <- test2 <- c();
while (statements>old_statements && ttl_runs<5)	{
	ttl_runs <- ttl_runs+1;
	old_statements <- statements;
	print("");
	print(paste("Iteration",ttl_runs));
	for (sp in 1:nspec)	{
		if (sp %in% progress)	{
			prop_done <- 100*match(sp,progress)/length(progress);
			if (prop_done<10)	{
				console_update <- paste("0",prop_done,"% done",sep="");
				} else {
				console_update <- paste(prop_done,"% done",sep="");
				}
			cat('\b\b\b\b\b\b\b\b\b',console_update);
			flush.console()
			}
		overlords <- which(spec_below==sp,arr.ind = T)[,1];		# species occurring above sp.
		aerie <- spec_above[sp,spec_above[sp,]>0];
		eagles <- as.vector(spec_above[overlords,])
		eagles <- sort(unique(eagles[eagles>0]));
		aerie <- unique(c(aerie,eagles));

		underlords <- which(spec_above==sp,arr.ind = T)[,1];	# species occurring below sp.
		goblin_town <- spec_below[sp,spec_below[sp,]>0];
		orcs <- as.vector(spec_below[underlords,])
		orcs <- sort(unique(orcs[orcs>0]));
		goblin_town <- unique(c(goblin_town,orcs));

		lokis <- c(sp,goblin_town[goblin_town %in% aerie]);
		goblin_town <- goblin_town[!goblin_town %in% lokis];
		aerie <- aerie[!aerie %in% lokis];
		lokis <- lokis[lokis!=sp];
		if (length(lokis)>0)		flag <- rbind(flag,cbind(rep(sp,length(lokis)),lokis));

		if (length(goblin_town)>0)	spec_below[sp,1:length(goblin_town)] <- goblin_town;
		if (length(aerie)>0)		spec_above[sp,1:length(aerie)] <- aerie;
		if (sum(is.na(spec_below))>0)	sp <- nspec;
		if (sum(is.na(spec_above))>0)	sp <- nspec;
		}
#	test1 <- rbind(test1,spec_above[match("Olenellus (Olenellus) nevadense",rownames(spec_above)),spec_above[match("Olenellus (Olenellus) nevadense",rownames(spec_above)),]>0]);
#	test2 <- rbind(test2,spec_below[match("Olenellus (Olenellus) nevadense",rownames(spec_below)),spec_below[match("Olenellus (Olenellus) nevadense",rownames(spec_below)),]>0]);
	statements <- sum(spec_below>0)+sum(spec_above>0);
	print("");
	}
new_occurrence_order <- occurrence_order;
new_occurrence_order[new_occurrence_order!=0] <- 0;
for (sp in 1:nspec)	{
	nn <- unique(spec_below[sp,spec_below[sp,]>0]);
	new_occurrence_order[sp,nn] <- 1;
	nn <- unique(spec_above[sp,spec_above[sp,]>0]);
	new_occurrence_order[sp,nn] <- -1;
	}
return(new_occurrence_order);
}


#### Routines to put probabilities on starts & stops of rock units ####
find_relevant_sections <- function(rock_unit,rock_superposition)	{
return(unique(rock_superposition$column[unique(which(rock_superposition==rock_unit,arr.ind = T)[,1])]));
}

# section <- this_section;
order_rocks_in_section <- function(section)	{
rock_nos <- unique(c(section$rock_no_sr_up,section$rock_no_sr_lo));
rock_formations <- unique(rbind(cbind(section$formation_no_up,section$rock_no_sr_up),
								cbind(section$formation_no_lo,section$rock_no_sr_lo)));
colnames(rock_formations) <- c("formation_no","rock_no");
rock_formations <- data.frame(rock_formations,stringsAsFactors = F);
srocks <- length(rock_nos);
formation_order <- vector(length=srocks);
for (sr in 1:srocks)	{
	ordered_rocks <- accersi_rock_order_in_section(rock_no=rock_nos[sr],section,topdown = F);
	ordered_formations_all <- rock_formations$formation_no[match(ordered_rocks,rock_formations$rock_no)];
	of <- 1;
	formation_order[match(ordered_formations_all[1],rock_formations$formation_no)] <- 1;
	unique_formations <- unique(ordered_formations_all);
	rock_order <- match(rock_formations$formation_no,unique_formations);
	if (sum(ordered_rocks==ordered_formations_all)>0)	{
		member_nos <- ordered_rocks[ordered_rocks!=ordered_formations_all];
		member_forms <- unique(rock_formations$formation_no[match(member_nos,rock_formations$rock_no)]);
		for (mf in 1:length(member_forms))	{
			ttl_mems <- sum(member_forms[mf]==rock_formations$formation_no[match(ordered_rocks,rock_formations$rock_no)]);
			wts <- seq(0,1-(1/ttl_mems),1/ttl_mems);
			}
		}
	while (of < length(ordered_formations_all))	{
		of <- of+1;
		if (ordered_formations_all[of]!=ordered_formations_all[of-1])	{
			this_rock <- match(ordered_formations_all[of],rock_formations$formation_no);
			lower_rock <- match(ordered_formations_all[of-1],rock_formations$formation_no);
			formation_order[this_rock] <- floor(formation_order[lower_rock])+1;
			}
		}

	ordered_formations <- unique(ordered_formations_all);
	basic_order_1 <- match(rock_formations$rock_no,ordered_formations);
	formation_order[!is.na(basic_order_1)] <- basic_order_1[!is.na(basic_order_1)];
#	if (length(ordered_rocks)>length())
	}
}

# putting rocks in order within sections
accersi_rock_order_in_section <- function(rock_no,section,topdown=F)	{
rock_order <- as.numeric(rock_no);
# find rocks underneath unit rock_no
start <- match(rock_no,section$rock_no_sr_up);
section <- subset(section,section$rock_no_sr_up!=section$rock_no_sr_lo);
while (!is.na(start))	{
	rock_order <- c(rock_order,as.numeric(section$rock_no_sr_lo[start]));
	start <- match(section$rock_no_sr_lo[start],section$rock_no_sr_up);
	}
# find rocks underneath unit rock_no
start <- match(rock_no,section$rock_no_sr_lo);
while (!is.na(start))	{
	rock_order <- c(as.numeric(section$rock_no_sr_up[start]),rock_order);
	start <- match(section$rock_no_sr_up[start],section$rock_no_sr_lo);
	}
if (!topdown)
	rock_order <- rock_order[length(rock_order):1]
return(rock_order);
}

accersi_number_of_sections <- function(rock_no_sr,superpositions)	{
links <- subset(superpositions,superpositions$rock_no_sr_up==rock_no_sr);
links <- rbind(links,subset(superpositions,superpositions$rock_no_sr_lo==rock_no_sr));
return(sum(length(unique(links$column))));
}

accersi_named_sections <- function(rock_no_sr,superpositions)	{
links <- subset(superpositions,superpositions$rock_no_sr_up==rock_no_sr);
links <- rbind(links,subset(superpositions,superpositions$rock_no_sr_lo==rock_no_sr));
return(sort(unique(links$column)));
}

accersi_rocks_in_same_section_and_time_slice <- function(rock_no,section,formations_from_members=F,temporal_precision=0.1)	{
rock_no_lowest <- rock_no_highest <- rock_no;
if (formations_from_members)	{
	this_formation <- unique(rbind(subset(section,section$formation_no_up==rock_no),
								   subset(section,section$formation_no_lo==rock_no)));
	all_rock_nos <- unique(data.frame(formation_no=as.numeric(c(section$formation_no_up,section$formation_no_lo)),
							rock_no_sr=as.numeric(c(section$rock_no_sr_up,section$rock_no_sr_lo)),
							rock_no=as.numeric(c(section$rock_no_up,section$rock_no_lo)),stringsAsFactors = F));
	this_rock_nos <- all_rock_nos[unique(which(all_rock_nos==rock_no,arr.ind = T)[,1]),];
#	this_rock_member_nos <- subset(this_rock_nos,this_rock_nos$rock_no_sr!=rock_no);
	this_rock_member_nos <- this_rock_nos[this_rock_nos$rock_no_sr!=rock_no,];
#	sapply(rock_no=this_rock_member_nos$rock_no,accersi_rock_order_in_section,section=section)
	if (nrow(this_rock_member_nos)>0)	{
		rock_seqs <- accersi_rock_order_in_section(rock_no=this_rock_member_nos$rock_no[1],section);
		useful_members <- this_rock_member_nos$rock_no[this_rock_member_nos$rock_no %in% rock_seqs];
		useful_member_order <- match(useful_members,rock_seqs);
		rock_no_highest <- useful_members[match(min(useful_member_order),useful_member_order)];
		rock_no_lowest <- useful_members[match(max(useful_member_order),useful_member_order)];
		if (is.na(rock_no_highest) && is.na(rock_no_lowest))	{
			rock_no_lowest <- rock_no_highest <- rock_no;
			}
		}
	}

rock_order_1 <- accersi_rock_order_in_section(rock_no=rock_no_lowest,section);
rock_order_2 <- accersi_rock_order_in_section(rock_no=rock_no_highest,section);
rock_order <- data.frame(rock_no=as.numeric(rock_order_1),
						 rock_unit=as.character(section$full_name_up[rock_info]),
						 interval_lb=as.character(section$interval_lb_up[rock_info]),
						 interval_ub=as.character(section$interval_ub_up[rock_info]),
						 bin_lb=as.numeric(section$bin_lb_up[rock_info]),
						 bin_ub=as.numeric(section$bin_ub_up[rock_info]),
						 stringsAsFactors = F);
nad <- (1:nrow(rock_order))[is.na(rock_order$rock_unit)];
na_buster <- match(rock_order$rock_no[nad],section$rock_no_sr_lo);
rock_order$rock_unit[nad] <- section$full_name_lo[na_buster];
rock_order$interval_lb[nad] <- section$interval_lb_lo[na_buster];
rock_order$interval_ub[nad] <- section$interval_ub_lo[na_buster];
rock_order$bin_lb[nad] <- section$bin_lb_lo[na_buster];
rock_order$bin_ub[nad] <- section$bin_ub_lo[na_buster];

this_rock_up <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
ttl_rock <- nrow(rock_order);
for (tr in 1:(ttl_rock-1))	{
	lower_rocks <- (tr+1):ttl_rock;
	rock_order$interval_ub[lower_rocks[rock_order$bin_ub[lower_rocks]>rock_order$bin_ub[tr]]] <- rock_order$interval_ub[tr];
	rock_order$bin_ub[lower_rocks[rock_order$bin_ub[lower_rocks]>rock_order$bin_ub[tr]]] <- rock_order$bin_ub[tr];
	}
for (tr in ttl_rock:2)	{
	upper_rocks <- 1:(tr-1);
	rock_order$interval_lb[rock_order$bin_lb[upper_rocks]<rock_order$bin_lb[tr]] <- rock_order$interval_lb[tr];
	rock_order$bin_lb[rock_order$bin_lb[upper_rocks]<rock_order$bin_lb[tr]] <- rock_order$bin_lb[tr];
	}

papa_bear <- (1:ttl_rock)[rock_order$bin_lb<=rock_order$bin_ub[this_rock_up]];
mama_bear <- (1:ttl_rock)[rock_order$bin_ub>=rock_order$bin_lb[this_rock_lo]];
baby_bear <- mama_bear[mama_bear %in% papa_bear];
return(rock_order[baby_bear,]);
}

accersi_relative_position_within_time_slice_full <- function(rock_no,section,radiometric_dates,rock_to_isotope_excursion,formations_from_members=T,temporal_precision=0.1)	{
rock_no_lowest <- rock_no_highest <- rock_no;
if (formations_from_members)	{
	this_formation <- unique(rbind(subset(section,section$formation_no_up==rock_no),
								   subset(section,section$formation_no_lo==rock_no)));
	all_rock_nos <- unique(data.frame(formation_no=as.numeric(c(section$formation_no_up,section$formation_no_lo)),
									  rock_no_sr=as.numeric(c(section$rock_no_sr_up,section$rock_no_sr_lo)),
									  rock_no=as.numeric(c(section$rock_no_up,section$rock_no_lo)),stringsAsFactors = F));
	this_rock_nos <- all_rock_nos[unique(which(all_rock_nos==rock_no,arr.ind = T)[,1]),];
	#	this_rock_member_nos <- subset(this_rock_nos,this_rock_nos$rock_no_sr!=rock_no);
	this_rock_member_nos <- this_rock_nos[this_rock_nos$rock_no_sr!=rock_no,];
	#	sapply(rock_no=this_rock_member_nos$rock_no,accersi_rock_order_in_section,section=section)
	if (nrow(this_rock_member_nos)>0)	{
		rock_seqs <- accersi_rock_order_in_section(rock_no=this_rock_member_nos$rock_no[1],section);
		useful_members <- this_rock_member_nos$rock_no[this_rock_member_nos$rock_no %in% rock_seqs];
		useful_member_order <- match(useful_members,rock_seqs);
		rock_no_highest <- useful_members[match(min(useful_member_order),useful_member_order)];
		rock_no_lowest <- useful_members[match(max(useful_member_order),useful_member_order)];
		if (is.na(rock_no_highest) && is.na(rock_no_lowest))	{
			rock_no_lowest <- rock_no_highest <- rock_no;
			}
		}
	}

rock_order_1 <- accersi_rock_order_in_section(rock_no=rock_no_lowest,section);
rock_order_2 <- accersi_rock_order_in_section(rock_no=rock_no_highest,section);
rock_info <- match(rock_order_1,section$rock_no_sr_up);
rock_order <- data.frame(rock_no=as.numeric(rock_order_1),
						 rock_unit=as.character(section$full_name_up[rock_info]),
						 ma_lb=as.numeric(section$up_ma_lb[rock_info]),
						 ma_ub=as.numeric(section$up_ma_ub[rock_info]),stringsAsFactors = F);
nad <- (1:nrow(rock_order))[is.na(rock_order$rock_unit)];
na_buster <- match(rock_order$rock_no[nad],section$rock_no_sr_lo);
rock_order$rock_unit[nad] <- section$full_name_lo[na_buster];
rock_order$ma_lb[nad] <- section$lo_ma_lb[na_buster];
rock_order$ma_ub[nad] <- section$lo_ma_ub[na_buster];

this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]);
#rock_order <- rbind(rock_order,data.frame(rock_no=as.numeric(20000),rock_unit=as.character("Doushantuo (Lower)"),ma_lb=as.numeric(560),ma_ub=as.numeric(541),stringsAsFactors = F))

this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
seq_rocks <- nrow(rock_order);
rownames(rock_order) <- 1:seq_rocks;
this_rock_lb_init <- rock_order$ma_lb[this_rock_lo];
this_rock_ub_init <- rock_order$ma_ub[this_rock_hi];

# older rocks must end by the latest date of younger rocks #
#rock_no    rock_unit ma_lb ma_ub			rock_no    rock_unit ma_lb ma_ub
#1    7379 Hongjingshao 519.4 515.0		1    7379 Hongjingshao 519.4 515.0
#2    2001   Qiongzhusi 519.4 516.2		2    2001   Qiongzhusi 519.4 516.2
#3    7381   Zhujiaqing 541.0 521.0		3    7381   Zhujiaqing 541.0 521.0
#4    3494     Dengying 550.0 517.6		4    3494     Dengying 550.0 521.0
seq_rocks <- nrow(rock_order);
ro <- 1;
#for (ro in 1:(nrow(rock_order)-1))	{
while (ro < seq_rocks)	{
	older_rocks <- (ro+1):seq_rocks;
	underlapping <- older_rocks[rock_order$ma_ub_mx[older_rocks]<rock_order$ma_lb_mn[ro]];
	#	rock_order$ma_ub[underlapping[rock_order$ma_ub[underlapping]<rock_order$ma_ub[ro]]] <- rock_order$ma_ub[ro];
	rock_order$ma_ub[ro] <- max(rock_order$ma_ub[1:ro]);
	rock_order$ma_lb[ro] <- min(rock_order$ma_lb[ro],min(rock_order$ma_lb[older_rocks]));
	ro <- ro+1;
	}

ro <- this_rock_hi;
while (ro < seq_rocks)	{
	older_rocks <- (ro+1):seq_rocks;
	constrained <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_ub[ro]];
	rock_order$ma_ub[constrained] <- rock_order$ma_ub[ro];
	ro <- ro+1;
	}

# younger rocks cannot start before the onset date of older rocks #
ro <- this_rock_lo;
while (ro > 1)	{
	younger_rocks <- (ro-1):1;
	restrained <- younger_rocks[rock_order$ma_lb[younger_rocks]>rock_order$ma_lb[ro]];
	rock_order$ma_lb[restrained] <- rock_order$ma_lb[ro];
	ro <- ro-1;
	}

# redo the rocks to give the upper & lower bounds on the onset and end
cnames <- colnames(rock_order);
cnames[cnames=="ma_lb"] <- "ma_lb_mx";
cnames[cnames=="ma_ub"] <- "ma_ub_mn";
colnames(rock_order) <- cnames;
rock_order <- tibble::add_column(rock_order, ma_lb_mn=rock_order$ma_ub_mn, .after = "ma_lb_mx");
rock_order <- tibble::add_column(rock_order, ma_ub_mx=rock_order$ma_lb_mx, .after = "ma_lb_mn");
for (i in 1:nrow(rock_order))	{
	if (length(which(radiometric_dates==rock_order$rock_no[i],arr.ind = T)[,1])>0)	{
		bb <- unique(which(radiometric_dates==rock_order$rock_no[i],arr.ind = T)[,1]);
		#radiometric_dates[bb,]
		radioactive_lo <- c(rock_order$ma_lb_mn[i],radiometric_dates$ma_lb[bb][radiometric_dates$ma_lb[bb]<rock_order$ma_lb_mx[i]]);
		if (rock_order$ma_lb_mn[i] < max(radioactive_lo))
			rock_order$ma_lb_mn[i] <- max(radioactive_lo);
		radioactive_hi <- c(rock_order$ma_ub_mx[i],radiometric_dates$ma_ub[bb][radiometric_dates$ma_ub[bb]>rock_order$ma_ub_mn[i]]);
		if (rock_order$ma_ub_mx[i] > min(radioactive_hi))
			rock_order$ma_ub_mx[i] <- min(radioactive_hi);
#		if (sum(radiometric_dates$ma_ub[bb]>rock_order$ma_lb_mn[i])>0)	{
#			nd <- radiometric_dates$ma_ub[bb][radiometric_dates$ma_ub[bb]>rock_order$ma_lb_mn[i]];
#			rock_order$ma_lb_mn[i] <- min(nd);
#			rock_order$ma_ub_mx[i] <- max(radiometric_dates$ma_ub[bb]);	# earliest possible LA
#			rock_order$ma_lb_mn[i] <- min(radiometric_dates$ma_ub[bb]);	# latest possible FA
#			}
		}
	}

# lower rocks must have upper bound on start time >= lower bound of radiometric dates from higher rocks
sr <- 0;
while (sr < (seq_rocks-1))	{
	sr <- sr + 1;
	rem_rocks <- (sr+1):seq_rocks;
	rock_order$ma_lb_mn[rem_rocks][rock_order$ma_lb_mn[rem_rocks]<rock_order$ma_lb_mn[sr]] <- rock_order$ma_lb_mn[sr];
	}

# higher rocks must have lower bound on end time <= upper bound of radiometric dates from lower rocks
sr <- seq_rocks;
while (sr > 1)	{
	rem_rocks <- 1:(sr-1)
#	rock_order$ma_ub_mx[rem_rocks]
	rock_order$ma_ub_mx[rem_rocks][rock_order$ma_ub_mx[rem_rocks] > rock_order$ma_ub_mx[sr]] <- rock_order$ma_ub_mx[sr];
	sr <- sr - 1;
	}

#  rock_no           rock_unit ma_lb_mn ma_lb_mx ma_ub_mn ma_ub_mx
#1   17808 Dengying (Hamajing)      635    541.0    635.0      541
#2   12134 Doushantuo (Miaohe)      560    550.4    551.8      550
#3   12134 Doushantuo (Lower)       560    550.0    560.0      550
# Hamajing cannot start before 551.8
# Lower cannot end after 550.4
#1   17808 Dengying (Hamajing)      551.8    541.0    551.8      541.0
#2   12134 Doushantuo (Miaohe)      560.0    550.4    551.8      550.0
#3   21340 Doushantuo (Lower)       560.0    550.4    560.0      550.4
#
# new problem from Hauchabfontein section.
#  rock_no              rock_unit ma_lb_mx ma_lb_mn ma_ub_mx ma_ub_mn
#1    2313                Nomtsas      550    539.3    537.1    532.0
#2   12678     Urusis (Spitzkopf)      550    545.5    538.8    538.7
#3   12679 Urusis (Feldschuhhorn)      550    545.5    550.0    538.7
#4    2330          Urusis (Huns)      550    545.5    550.0    538.7
#5   12680         Urusis (Nasep)      550    545.5    550.0    538.7
#6   14242   Nudaus (Vingerbreek)      550    545.5    550.0    538.7
# ma_ub_mn for rocks below Spitzkopf should be 545.5 or lower
# ma_lb_mn for Nomtsas & ma_ub_mx for Spitzkopf are in conflict

this_rock_lb_mx <- max(rock_order$ma_lb_mx[this_rock_lo:this_rock_hi]);
this_rock_lb_mn <- min(rock_order$ma_lb_mn[this_rock_lo:this_rock_hi]);
this_rock_ub_mx <- max(rock_order$ma_ub_mx[this_rock_lo:this_rock_hi]);
this_rock_ub_mn <- min(rock_order$ma_ub_mn[this_rock_lo:this_rock_hi]);
rock_order <- subset(rock_order,rock_order$ma_lb_mx>this_rock_ub_mn);
rock_order <- subset(rock_order,rock_order$ma_ub_mn<this_rock_lb_mx);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);

seq_rocks <- nrow(rock_order);
#poss_higher_rocks <- sum(rock_order$ma_ub[1:this_rock_hi]<rock_order$ma_ub[this_rock_hi]);
poss_higher_rocks <- sum(rock_order$ma_ub_mn[1:this_rock_hi]<rock_order$ma_ub_mn[this_rock_hi]);
nesc_higher_rocks <- this_rock_hi-(poss_higher_rocks+1);
poss_lower_rocks <-  sum(rock_order$ma_lb_mx[this_rock_lo:seq_rocks]>rock_order$ma_lb_mx[this_rock_lo]);
nesc_lower_rocks <- (seq_rocks-this_rock_lo)-poss_lower_rocks;

fuzzy_combos <- c();
plr <- nesc_lower_rocks;
while (plr <= (nesc_lower_rocks+poss_lower_rocks))	{
	fuzzy_combos <- rbind(fuzzy_combos,cbind(rep(plr,poss_higher_rocks+1),nesc_higher_rocks:(nesc_higher_rocks+poss_higher_rocks)));
	plr <- plr + 1;
	}

time_slices <- seq(this_rock_ub_mn+(temporal_precision/2),this_rock_lb_mx-(temporal_precision/2),by=temporal_precision);
ttl_time_slices <- length(time_slices);
time_slices_rev <- sort(time_slices,decreasing=T);		# oldest times first now.
ttl_span <- abs(this_rock_ub_mn-this_rock_lb_mx);
ppa <- abs(this_rock_lb_mx-time_slices)/ttl_span;	# time slices from oldest possible date as fraction of total span
#kill_fred <- (time_slices-this_rock_ub_mn)/poss_span;

latest_poss_start <- max(rock_order$ma_lb_mn[this_rock_hi:this_rock_lo]);	# get latest possible start time;
earliest_poss_start <- max(rock_order$ma_lb_mx[this_rock_hi:this_rock_lo]);	# get latest possible start time;
time_slices_st <- seq(latest_poss_start+(temporal_precision/2),earliest_poss_start-(temporal_precision/2),by=temporal_precision);
pps <- abs(this_rock_lb_mx-time_slices_st)/ttl_span;

earliest_poss_end <- min(rock_order$ma_ub_mx[this_rock_hi:this_rock_lo]);	# get earliest possible end time;
latest_poss_end <- min(rock_order$ma_ub_mn[this_rock_hi:this_rock_lo]);	# get earliest possible end time;
time_slices_en <- seq(latest_poss_end+(temporal_precision/2),earliest_poss_end-(temporal_precision/2),by=temporal_precision);
ppe <- abs(this_rock_lb_mx-time_slices_en)/ttl_span;

p_start <- p_end <- c();
p_st <- p_en <- rep(0,ttl_time_slices);
st_cells <- match(round(time_slices_st,3),round(time_slices,3));
en_cells <- match(round(time_slices_en,3),round(time_slices,3));
for (fc in 1:nrow(fuzzy_combos))	{
	events_intervening <- abs(this_rock_hi-this_rock_lo);
	events_prior_start <- fuzzy_combos[fc,1];
	events_prior_end <- events_prior_start+events_intervening+1;
	events_after_end <- fuzzy_combos[fc,2];
	events_after_start <- events_after_end+events_intervening+1;
	ttl_events <- 1+events_after_start+events_prior_start;

	p_event <- dpois(1,ttl_events*temporal_precision/ttl_span);
	pp <- pps;	# proportion of the way through the total span that each possible starting point is.
	p_st[st_cells] <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_start,sbN=events_after_start);
#	p_st[st_cells] <- p_st[st_cells]/sum(p_st[st_cells])
#	p_st[st_cells] <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start,span=ttl_events)/sum(sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start,span=ttl_events));
#	modifier_st <- sum(p_st)/dpois(ttl_events,ttl_events);
#	p_st[time_slices<latest_poss_start] <- 0;
	p_st <- p_st/sum(p_st);
	p_start <- cbind(p_start,p_st);

	pp <- ppe;
#	p_en[en_cells] <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_end,sbN=events_after_end,span=ttl_events);
	p_en[en_cells] <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_end,sbN=events_after_end);
	pp <- ppa;	#
#	modifier_en <- sum(sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_end,sbN=events_after_end));
#	modifier_en <- sum(sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_end,sbN=events_after_end,span=ttl_events))/dpois(ttl_events,ttl_events);
#	p_en <- p_en/dpois(ttl_events,ttl_events);	# condition this on N of N expected events having happened
#	p_en[time_slices_rev>earliest_poss_end] <- 0;
#	p_en <- modifier_en*(p_en/sum(p_en));
	p_en <- (p_en/sum(p_en));
	p_end <- cbind(p_end,p_en);
#	p_end <- cbind(p_end,p_en/sum(p_en));
	}

p_start <- rowSums(p_start)/nrow(fuzzy_combos);
p_end <- rowSums(p_end)/nrow(fuzzy_combos);

#cuml_p_started <- cumsum(p_start);
#cuml_p_not_ended <- 1-cumsum(p_end);
p_present <- cumsum(p_start[ttl_time_slices:1])*(1-cumsum(p_end[ttl_time_slices:1]));
p_present <- p_present[ttl_time_slices:1]
# kluge! if extinction is limited to one unit, then problems happen
if (length((1:ttl_time_slices)[p_present %in% 0])>0)	{
	if (p_end[max((1:ttl_time_slices)[p_present %in% 0])]==1)
		p_present[max((1:ttl_time_slices)[p_present %in% 0])] <- 0.5;
	}

age_probs <- data.frame(poss_age=as.numeric(time_slices),
						prob_rock_present=as.numeric(p_present[ttl_time_slices:1]),
						prob_start=as.numeric(p_start),
						prob_end=as.numeric(p_end),stringsAsFactors = F);
#output <- list(data.frame(poss_age=as.numeric(time_slices),
#							  prob_rock_present=as.numeric(p_present),
#							  prob_start=as.numeric(p_start),
#							  prob_end=as.numeric(p_end),stringsAsFactors = F),
#				   p_end_given_start);
	#names(output) <- c("age_probs","lnp_end_given_start","p_end_given_start");
#names(output) <- c("age_probs","p_end_given_start");

return(age_probs);
}

accersi_relative_position_within_time_slice <- function(rock_no,section,formations_from_members=T,temporal_precision=0.1)	{
rock_no_lowest <- rock_no_highest <- rock_no;
if (formations_from_members)	{
	this_formation <- unique(rbind(subset(section,section$formation_no_up==rock_no),
								   subset(section,section$formation_no_lo==rock_no)));
	all_rock_nos <- unique(data.frame(formation_no=as.numeric(c(section$formation_no_up,section$formation_no_lo)),
							rock_no_sr=as.numeric(c(section$rock_no_sr_up,section$rock_no_sr_lo)),
							rock_no=as.numeric(c(section$rock_no_up,section$rock_no_lo)),stringsAsFactors = F));
	this_rock_nos <- all_rock_nos[unique(which(all_rock_nos==rock_no,arr.ind = T)[,1]),];
#	this_rock_member_nos <- subset(this_rock_nos,this_rock_nos$rock_no_sr!=rock_no);
	this_rock_member_nos <- this_rock_nos[this_rock_nos$rock_no_sr!=rock_no,];
#	sapply(rock_no=this_rock_member_nos$rock_no,accersi_rock_order_in_section,section=section)
	if (nrow(this_rock_member_nos)>0)	{
		rock_seqs <- accersi_rock_order_in_section(rock_no=this_rock_member_nos$rock_no[1],section);
		useful_members <- this_rock_member_nos$rock_no[this_rock_member_nos$rock_no %in% rock_seqs];
		useful_member_order <- match(useful_members,rock_seqs);
		rock_no_highest <- useful_members[match(min(useful_member_order),useful_member_order)];
		rock_no_lowest <- useful_members[match(max(useful_member_order),useful_member_order)];
		if (is.na(rock_no_highest) && is.na(rock_no_lowest))	{
			rock_no_lowest <- rock_no_highest <- rock_no;
			}
		}
	}

rock_order_1 <- accersi_rock_order_in_section(rock_no=rock_no_lowest,section);
rock_order_2 <- accersi_rock_order_in_section(rock_no=rock_no_highest,section);
rock_info <- match(rock_order_1,section$rock_no_sr_up);
rock_order <- data.frame(rock_no=as.numeric(rock_order_1),
						 rock_unit=as.character(section$full_name_up[rock_info]),
						 ma_lb=as.numeric(section$up_ma_lb[rock_info]),
						 ma_ub=as.numeric(section$up_ma_ub[rock_info]),stringsAsFactors = F);
nad <- (1:nrow(rock_order))[is.na(rock_order$rock_unit)];
na_buster <- match(rock_order$rock_no[nad],section$rock_no_sr_lo);
rock_order$rock_unit[nad] <- section$full_name_lo[na_buster];
rock_order$ma_lb[nad] <- section$lo_ma_lb[na_buster];
rock_order$ma_ub[nad] <- section$lo_ma_ub[na_buster];

this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
seq_rocks <- nrow(rock_order);
rownames(rock_order) <- 1:seq_rocks;
this_rock_lb_init <- rock_order$ma_lb[this_rock_lo];
this_rock_ub_init <- rock_order$ma_ub[this_rock_hi];

# older rocks must end by the latest date of younger rocks #
#rock_no    rock_unit ma_lb ma_ub			rock_no    rock_unit ma_lb ma_ub
#1    7379 Hongjingshao 519.4 515.0		1    7379 Hongjingshao 519.4 515.0
#2    2001   Qiongzhusi 519.4 516.2		2    2001   Qiongzhusi 519.4 516.2
#3    7381   Zhujiaqing 541.0 521.0		3    7381   Zhujiaqing 541.0 521.0
#4    3494     Dengying 550.0 517.6		4    3494     Dengying 550.0 521.0

ro <- 1;
#for (ro in 1:(nrow(rock_order)-1))	{
while (ro < nrow(rock_order))	{
	older_rocks <- (ro+1):seq_rocks;
	underlapping <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_lb[ro]];
#	rock_order$ma_ub[underlapping[rock_order$ma_ub[underlapping]<rock_order$ma_ub[ro]]] <- rock_order$ma_ub[ro];
	rock_order$ma_ub[ro] <- max(rock_order$ma_ub[1:ro]);
	rock_order$ma_lb[ro] <- min(rock_order$ma_lb[ro],min(rock_order$ma_lb[older_rocks]));
	ro <- ro+1;
	}

ro <- this_rock_hi;
while (ro < nrow(rock_order))	{
	older_rocks <- (ro+1):seq_rocks;
	constrained <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_ub[ro]];
	rock_order$ma_ub[constrained] <- rock_order$ma_ub[ro];
	ro <- ro+1;
	}

# younger rocks cannot start before the onset date of older rocks #
ro <- this_rock_lo;
while (ro > 1)	{
	younger_rocks <- (ro-1):1;
	restrained <- younger_rocks[rock_order$ma_lb[younger_rocks]>rock_order$ma_lb[ro]];
	rock_order$ma_lb[restrained] <- rock_order$ma_lb[ro];
	ro <- ro-1;
	}
this_rock_lb <- rock_order$ma_lb[this_rock_lo];
this_rock_ub <- rock_order$ma_ub[this_rock_hi];
rock_order <- subset(rock_order,rock_order$ma_lb>this_rock_ub);
rock_order <- subset(rock_order,rock_order$ma_ub<this_rock_lb);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);

#subset(rock_order,!as.logical((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<rock_order$ma_ub[this_rock_hi])));
# in this case, it should be done in four ways.
#	1: Blueflower (Disc) is the last rock unit before 550
#	2: Risky is the last rock unit before 550
#	3: Ingta is the last rock unit before 550
#	3: Vampire is the last rock unit before 550
#  rock_no                        rock_unit ma_lb ma_ub
#1    7265                          Vampire   560   521
#2   17804                            Ingta   560   521
#3   16466                            Risky   560   541
#4   17350                Blueflower (Disc)   560   550
#5   17351            Blueflower (Yuletide)   560   550
#6   17352 Blueflower (Peritidal Carbonate)   560   550
#7   17353               Blueflower (Basal)   560   550
#8   16467                        Gametrail   560   550
seq_rocks <- nrow(rock_order);
poss_higher_rocks <- sum(rock_order$ma_ub[1:this_rock_hi]<rock_order$ma_ub[this_rock_hi]);
nesc_higher_rocks <- this_rock_hi-(poss_higher_rocks+1);
poss_lower_rocks <-  sum(rock_order$ma_lb[this_rock_lo:seq_rocks]>rock_order$ma_lb[this_rock_lo]);
nesc_lower_rocks <- (seq_rocks-this_rock_lo)-poss_lower_rocks;

fuzzy_combos <- c();
plr <- nesc_lower_rocks;
while (plr <= (nesc_lower_rocks+poss_lower_rocks))	{
	fuzzy_combos <- rbind(fuzzy_combos,cbind(rep(plr,poss_higher_rocks+1),nesc_higher_rocks:(nesc_higher_rocks+poss_higher_rocks)));
	plr <- plr + 1;
	}

#num_overlapping_higher_rocks <- sum((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<=rock_order$ma_ub[this_rock_hi]));
#num_overlapping_lower_rocks <- sum((rock_order$ma_lb>=rock_order$ma_lb[this_rock_lo]) * (rock_order$ma_ub>rock_order$ma_ub[this_rock_hi]));

time_slices <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision);
time_slices_fine <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision/10);
ttl_time_slices <- length(time_slices);
time_slices_rev <- sort(time_slices,decreasing=T);
pp <- (time_slices-this_rock_ub)/(this_rock_lb-this_rock_ub);
ttl_span <- abs(this_rock_ub-this_rock_lb);

p_start <- p_end <- c();
for (fc in 1:nrow(fuzzy_combos))	{
	events_intervening <- abs(this_rock_hi-this_rock_lo);
	events_prior_start <- fuzzy_combos[fc,1];
	events_prior_end <- events_prior_start+events_intervening+1;
	events_after_end <- fuzzy_combos[fc,2];
	events_after_start <- events_after_end+events_intervening+1;
	ttl_events <- 1+events_after_start+events_prior_start;

	min_start <- max(rock_order$ma_ub[this_rock_hi:this_rock_lo]);	# get latest possible start time;
	max_end <- min(rock_order$ma_lb[this_rock_hi:this_rock_lo]);	# get earliest possible end time;

	p_event <- dpois(1,ttl_events*temporal_precision/ttl_span);
	p_st <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_start,sbN=events_after_start);
#	p_st <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start,span=ttl_events);
	p_st[time_slices_rev<min_start] <- 0;
#	p_start <- cbind(p_start,p_st/sum(p_st));
	p_st <- p_st/dpois(ttl_events,ttl_events);	# condition this on N of N expected events having happened
	p_start <- cbind(p_start,p_st);
	#p_end <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete,p_event,prN=events_prior_end,sbN=events_after_end);
	p_en <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_end,sbN=events_after_end);
#	p_en <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_end,sbN=events_after_end,span=ttl_events);
	p_en[time_slices_rev>max_end] <- 0;
	p_en <- p_en/dpois(ttl_events,ttl_events);	# condition this on N of N expected events having happened
	p_end <- cbind(p_end,p_en);
#	p_end <- cbind(p_end,p_en/sum(p_en));
	}

p_start <- rowSums(p_start)/nrow(fuzzy_combos);
p_end <- rowSums(p_end)/nrow(fuzzy_combos);
pt <- pp;
p_end_given_start <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
time_slices <- time_slices[ttl_time_slices:1];
colnames(p_end_given_start) <- rownames(p_end_given_start) <- time_slices;
p_end_cond <- p_end_given_start;
p_event <- ttl_events/ttl_time_slices;
for (i in 1:ttl_time_slices)	{
	p_end_given_start[i,pt>=pp[i]] <- sapply(pt[pt>=pp[i]],accersi_p_end_given_start_time,ps=pp[i],p_event,intervening_events = events_intervening,subsequent_events = events_after_end,total_events = ttl_events);
	p_end_given_start[i,] <- p_end_given_start[i,]/sum(p_end_given_start[i,]);
	p_end_cond[i,] <- p_start[i]*p_end_given_start[i,];
	}

p_started_cum <- cumsum(p_start);
p_not_ended_cum <- 1-cumsum(p_end);
p_present <- p_started_cum*p_not_ended_cum;

output <- list(data.frame(poss_age=as.numeric(time_slices),
						  prob_rock_present=as.numeric(p_present),
						  prob_start=as.numeric(p_start),
						  prob_end=as.numeric(p_end),stringsAsFactors = F),
#			   lnp_end_given_start_matrix,
			   p_end_given_start);
#names(output) <- c("age_probs","lnp_end_given_start","p_end_given_start");
names(output) <- c("age_probs","p_end_given_start");

return(output)
}

accersi_relative_position_within_time_slice_gamma <- function(rock_no,section,formations_from_members=T,temporal_precision=0.1)	{
rock_no_lowest <- rock_no_highest <- rock_no;
if (formations_from_members)	{
	this_formation <- unique(rbind(subset(section,section$formation_no_up==rock_no),
								   subset(section,section$formation_no_lo==rock_no)));
	all_rock_nos <- unique(data.frame(formation_no=as.numeric(c(section$formation_no_up,section$formation_no_lo)),
							rock_no_sr=as.numeric(c(section$rock_no_sr_up,section$rock_no_sr_lo)),
							rock_no=as.numeric(c(section$rock_no_up,section$rock_no_lo)),stringsAsFactors = F));
	this_rock_nos <- all_rock_nos[unique(which(all_rock_nos==rock_no,arr.ind = T)[,1]),];
#	this_rock_member_nos <- subset(this_rock_nos,this_rock_nos$rock_no_sr!=rock_no);
	this_rock_member_nos <- this_rock_nos[this_rock_nos$rock_no_sr!=rock_no,];
#	sapply(rock_no=this_rock_member_nos$rock_no,accersi_rock_order_in_section,section=section)
	if (nrow(this_rock_member_nos)>0)	{
		rock_seqs <- accersi_rock_order_in_section(rock_no=this_rock_member_nos$rock_no[1],section);
		useful_members <- this_rock_member_nos$rock_no[this_rock_member_nos$rock_no %in% rock_seqs];
		useful_member_order <- match(useful_members,rock_seqs);
		rock_no_highest <- useful_members[match(min(useful_member_order),useful_member_order)];
		rock_no_lowest <- useful_members[match(max(useful_member_order),useful_member_order)];
		if (is.na(rock_no_highest) && is.na(rock_no_lowest))	{
			rock_no_lowest <- rock_no_highest <- rock_no;
			}
		}
	}

rock_order_1 <- accersi_rock_order_in_section(rock_no=rock_no_lowest,section);
rock_order_2 <- accersi_rock_order_in_section(rock_no=rock_no_highest,section);
rock_info <- match(rock_order_1,section$rock_no_sr_up);
rock_order <- data.frame(rock_no=as.numeric(rock_order_1),
						 rock_unit=as.character(section$full_name_up[rock_info]),
						 ma_lb=as.numeric(section$up_ma_lb[rock_info]),
						 ma_ub=as.numeric(section$up_ma_ub[rock_info]),stringsAsFactors = F);
nad <- (1:nrow(rock_order))[is.na(rock_order$rock_unit)];
na_buster <- match(rock_order$rock_no[nad],section$rock_no_sr_lo);
rock_order$rock_unit[nad] <- section$full_name_lo[na_buster];
rock_order$ma_lb[nad] <- section$lo_ma_lb[na_buster];
rock_order$ma_ub[nad] <- section$lo_ma_ub[na_buster];

this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
seq_rocks <- nrow(rock_order);
rownames(rock_order) <- 1:seq_rocks;
this_rock_lb_init <- rock_order$ma_lb[this_rock_lo];
this_rock_ub_init <- rock_order$ma_ub[this_rock_hi];

# older rocks must end by the latest date of younger rocks #
#rock_no    rock_unit ma_lb ma_ub			rock_no    rock_unit ma_lb ma_ub
#1    7379 Hongjingshao 519.4 515.0		1    7379 Hongjingshao 519.4 515.0
#2    2001   Qiongzhusi 519.4 516.2		2    2001   Qiongzhusi 519.4 516.2
#3    7381   Zhujiaqing 541.0 521.0		3    7381   Zhujiaqing 541.0 521.0
#4    3494     Dengying 550.0 517.6		4    3494     Dengying 550.0 521.0

ro <- 1;
#for (ro in 1:(nrow(rock_order)-1))	{
while (ro < nrow(rock_order))	{
	older_rocks <- (ro+1):seq_rocks;
	underlapping <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_lb[ro]];
#	rock_order$ma_ub[underlapping[rock_order$ma_ub[underlapping]<rock_order$ma_ub[ro]]] <- rock_order$ma_ub[ro];
	rock_order$ma_ub[ro] <- max(rock_order$ma_ub[1:ro]);
	rock_order$ma_lb[ro] <- min(rock_order$ma_lb[ro],min(rock_order$ma_lb[older_rocks]));
	ro <- ro+1;
	}

ro <- this_rock_hi;
while (ro < nrow(rock_order))	{
	older_rocks <- (ro+1):seq_rocks;
	constrained <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_ub[ro]];
	rock_order$ma_ub[constrained] <- rock_order$ma_ub[ro];
	ro <- ro+1;
	}

# younger rocks cannot start before the onset date of older rocks #
ro <- this_rock_lo;
while (ro > 1)	{
	younger_rocks <- (ro-1):1;
	restrained <- younger_rocks[rock_order$ma_lb[younger_rocks]>rock_order$ma_lb[ro]];
	rock_order$ma_lb[restrained] <- rock_order$ma_lb[ro];
	ro <- ro-1;
	}
# in this case, it should be done in four ways.
#	1: Blueflower (Disc) is the last rock unit before 550
#	2: Risky is the last rock unit before 550
#	3: Ingta is the last rock unit before 550
#	3: Vampire is the last rock unit before 550
#  rock_no                        rock_unit ma_lb ma_ub
#1    7265                          Vampire   560   521
#2   17804                            Ingta   560   521
#3   16466                            Risky   560   541
#4   17350                Blueflower (Disc)   560   550
#5   17351            Blueflower (Yuletide)   560   550
#6   17352 Blueflower (Peritidal Carbonate)   560   550
#7   17353               Blueflower (Basal)   560   550
#8   16467                        Gametrail   560   550

this_rock_lb <- rock_order$ma_lb[this_rock_lo];
this_rock_ub <- rock_order$ma_ub[this_rock_hi];
rock_order <- subset(rock_order,rock_order$ma_lb>this_rock_ub);
rock_order <- subset(rock_order,rock_order$ma_ub<this_rock_lb);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);

#subset(rock_order,!as.logical((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<rock_order$ma_ub[this_rock_hi])));
seq_rocks <- nrow(rock_order);
poss_higher_rocks <- sum(rock_order$ma_ub[1:this_rock_hi]<rock_order$ma_ub[this_rock_hi]);
nesc_higher_rocks <- this_rock_hi-(poss_higher_rocks+1);
poss_lower_rocks <-  sum(rock_order$ma_lb[this_rock_lo:seq_rocks]>rock_order$ma_lb[this_rock_lo]);
nesc_lower_rocks <- this_rock_hi-(poss_higher_rocks+1);

#num_overlapping_higher_rocks <- sum((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<=rock_order$ma_ub[this_rock_hi]));
#num_overlapping_lower_rocks <- sum((rock_order$ma_lb>=rock_order$ma_lb[this_rock_lo]) * (rock_order$ma_ub>rock_order$ma_ub[this_rock_hi]));

# put in part about assessing formations with members #
ttl_events <- seq_rocks+1;
if (rock_order$ma_ub[1]<rock_order$ma_ub[this_rock_hi])
	ttl_events <- ttl_events-1;
if (rock_order$ma_lb[seq_rocks]>rock_order$ma_lb[this_rock_lo])
	ttl_events <- ttl_events-1;

events_prior_start <- sum(rock_order$ma_lb[this_rock_lo:seq_rocks]==rock_order$ma_lb[this_rock_lo])-1;
events_prior_end <- sum(rock_order$ma_lb[this_rock_hi:seq_rocks]<=rock_order$ma_lb[this_rock_lo]);
events_after_start <- sum(rock_order$ma_ub[1:this_rock_lo]>=rock_order$ma_ub[this_rock_hi]);
events_after_end <- sum(rock_order$ma_ub[1:this_rock_hi]==rock_order$ma_ub[this_rock_hi])-1;
events_intervening <- events_prior_end-events_prior_start-1;
min_start <- max(rock_order$ma_ub[this_rock_hi:this_rock_lo]);	# get latest possible start time;
max_end <- min(rock_order$ma_lb[this_rock_hi:this_rock_lo]);	# get earliest possible end time;

# fractional overlap:
# rock_no	rock_unit	ma_lb	ma_ub
# 1000		Shire		521.0	515.4
# 1001		Buckland	523.8	517.6
# set up situation where we have half an event from 521.0 onwards
# Buckland cannot stop until 521 or later
# Buckland cannot start after 517.6 or later
# consider doing this with Shire appearing before 517.6 and after, with P[before] = (521.0-517.6)/(521.0-514.4)
if (num_overlapping_higher_rocks>0)	{
	upper_breaks <- rock_order$ma_lb[rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]];
	upper_breaks <- sort(upper_breaks[upper_breaks<rock_order$ma_lb[this_rock_hi]]);
	} else	{
	upper_breaks <- rock_order$ma_ub[this_rock_hi];
	}
if (num_overlapping_lower_rocks>0)	{
	lower_breaks <- rock_order$ma_ub[rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]];
	lower_breaks <- sort(lower_breaks[lower_breaks>rock_order$ma_ub[this_rock_lo]]);
	} else	{
	lower_breaks <- rock_order$ma_lb[this_rock_lo];
	}

time_slices <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision);
time_slices_fine <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision/10);
ttl_time_slices <- length(time_slices);
time_slices_rev <- sort(time_slices,decreasing=T);
pp <- (time_slices-this_rock_ub)/(this_rock_lb-this_rock_ub);
ttl_span <- abs(this_rock_ub-this_rock_lb);
#p_event <- ttl_events/length(time_slices);
#p_start <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete,p_event,prN=events_prior_start,sbN=events_after_start);
p_event <- dpois(1,ttl_events*temporal_precision/ttl_span);
p_start <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_start,sbN=events_after_start);
#p_start <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start,span=ttl_events);
p_start[time_slices_rev<min_start] <- 0;
p_start <- p_start/sum(p_start);
#p_end <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete,p_event,prN=events_prior_end,sbN=events_after_end);
p_end <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_end,sbN=events_after_end);
p_end[time_slices_rev>max_end] <- 0;
p_end <- p_end/sum(p_end);

#p_event <- dpois(1,ttl_events/length(time_slices));
#p_start <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start);
#p_start <- p_start/sum(p_start);
#p_end <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_end,sbN=events_after_end);
#p_end_very <- dpois(1,ttl_events/(2*length(time_slices)));
#p_end <- p_end/sum(c(p_end,p_end_very));

pt <- pp;
p_end_given_start <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
time_slices <- time_slices[ttl_time_slices:1];
colnames(p_end_given_start) <- rownames(p_end_given_start) <- time_slices;
p_end_cond <- p_end_given_start;
p_event <- ttl_events/ttl_time_slices;
for (i in 1:ttl_time_slices)	{
	p_end_given_start[i,pt>=pp[i]] <- sapply(pt[pt>=pp[i]],accersi_p_end_given_start_time,ps=pp[i],p_event,intervening_events = events_intervening,subsequent_events = events_after_end,total_events = ttl_events);
	p_end_given_start[i,] <- p_end_given_start[i,]/sum(p_end_given_start[i,]);
	p_end_cond[i,] <- p_start[i]*p_end_given_start[i,];
	}

p_started_cum <- cumsum(p_start);
p_not_ended_cum <- 1-cumsum(p_end);
p_present <- p_started_cum*p_not_ended_cum;

output <- list(data.frame(poss_age=as.numeric(time_slices),
						  prob_rock_present=as.numeric(p_present),
						  prob_start=as.numeric(p_start),
						  prob_end=as.numeric(p_end),stringsAsFactors = F),
#			   lnp_end_given_start_matrix,
			   p_end_given_start);
#names(output) <- c("age_probs","lnp_end_given_start","p_end_given_start");
names(output) <- c("age_probs","p_end_given_start");

return(output)
}

accersi_relative_position_within_time_slice_beta <- function(rock_no,section,formations_from_members=T,temporal_precision=0.1)	{
rock_no_lowest <- rock_no_highest <- rock_no;
if (formations_from_members)	{
	this_formation <- unique(rbind(subset(section,section$formation_no_up==rock_no),
								   subset(section,section$formation_no_lo==rock_no)));
	all_rock_nos <- unique(data.frame(formation_no=as.numeric(c(section$formation_no_up,section$formation_no_lo)),
							rock_no_sr=as.numeric(c(section$rock_no_sr_up,section$rock_no_sr_lo)),
							rock_no=as.numeric(c(section$rock_no_up,section$rock_no_lo)),stringsAsFactors = F));
	this_rock_nos <- all_rock_nos[unique(which(all_rock_nos==rock_no,arr.ind = T)[,1]),];
#	this_rock_member_nos <- subset(this_rock_nos,this_rock_nos$rock_no_sr!=rock_no);
	this_rock_member_nos <- this_rock_nos[this_rock_nos$rock_no_sr!=rock_no,];
	if (nrow(this_rock_member_nos)>0)	{
		lo_no <- match("formation_no_lo",colnames(section));
		up_no <- match("formation_no_up",colnames(section));
		# highest rock is the last one underlying another formation
		rock_no_highest <- section$rock_no_sr_lo[(1:nrow(section))[as.logical((section[,up_no]!=rock_no)*(section[,lo_no]==rock_no))]];
		rock_no_highest <- rock_no_highest[rock_no_highest!=rock_no];
		rock_no_highest <- rock_no_highest[1];
		rock_no_highest <- rock_no_highest[!is.na(rock_no_highest)][1];
		# lowest rock is the first one overlying another formation
		rock_no_lowest <- section$rock_no_sr_up[(1:nrow(section))[as.logical((section[,up_no]==rock_no)*(section[,lo_no]!=rock_no))]];
		rock_no_lowest <- rock_no_lowest[rock_no_lowest!=rock_no];
		rock_no_lowest <- rock_no_lowest[!is.na(rock_no_lowest)][1];
		if (is.na(rock_no_highest) && is.na(rock_no_lowest))	{
			rock_no_lowest <- rock_no_highest <- rock_no;
			} else {
			if (is.na(rock_no_lowest))	{
				rock_no_lowest <- rock_no_highest;
				section <- subset(section,section$rock_no_sr_lo!=rock_no);
				} else if (is.na(rock_no_highest))	{
				rock_no_highest <- rock_no_lowest;
				section <- subset(section,section$rock_no_sr_up!=rock_no);
				}
			}
#			if (is.na(rock_no_highest) || rock_no_highest!=rock_no)	{
#				section <- subset(section,section$rock_no_sr_lo!=rock_no);
#				if (is.na(rock_no_lowest)) {
#					rock_no_highest <- rock_no_lowest;
#					} else	{
#
#					}
#				}
#			if (is.na(rock_no_lowest) || rock_no_lowest!=rock_no)	{
#				section <- subset(section,section$rock_no_sr_up!=rock_no);
#				rock_no_lowest <- rock_no_highest;
#				}
#			}
		}
	}
if (1!=1)	{
	if (nrow(this_formation)>2)	{
		# routine if the formation turns up 2+ times in the section
		member_search <- unique(data.frame(full_name=as.character(c(this_formation$full_name_up,this_formation$full_name_lo)),
										   formation_no=as.numeric(c(this_formation$formation_no_up,this_formation$formation_no_lo)),
										   rock_no=as.numeric(c(this_formation$rock_no_sr_up,this_formation$rock_no_sr_lo)),
										   stringsAsFactors = F));
#		member_search <- subset(member_search,member_search$formation_no==rock_no);
		member_search <- member_search[member_search$formation_no==rock_no,];
		if (sum(member_search$rock_no!=member_search$formation_no))	{
			formations_from_members <- T;
			} else	{
			formations_from_members <- F;
			}
		} else if (nrow(this_formation)==1)	{
		member_search <- unique(data.frame(full_name=as.character(c(this_formation$full_name_up,this_formation$full_name_lo)),
										   formation_no=as.numeric(c(this_formation$formation_no_up,this_formation$formation_no_lo)),
										   rock_no=as.numeric(c(this_formation$rock_no_sr_up,this_formation$rock_no_sr_lo)),
										   stringsAsFactors = F));
		member_search <- member_search[member_search$formation_no==rock_no,];
		formations_from_members <- F;
#		formations_from_members <- F;
#		if (this_formation$formation_no_up[1]==rock_no)	{
#			rock_no_highest <- rock_no_lowest <- this_formation$rock_no_up[1];
#			} else	{
#			rock_no_highest <- rock_no_lowest <- this_formation$rock_no_lo[1];
#			}
		}

	### now, let's just
#	if (formations_from_members)	{
#		section <- subset(section,section$rock_no_sr_up!=rock_no);
#		section <- subset(section,section$rock_no_sr_lo!=rock_no);

		# now, figure out which member is highest
#		member_layers <- subset(section,section$formation_no_lo==rock_no);
#		member_layers_lo <- subset(member_layers,member_layers$formation_no_up!=rock_no);
#		rock_no_highest <- member_layers_lo$rock_no_sr_lo[1];

		# now, figure out which member is highest
#		member_layers <- subset(section,section$formation_no_up==rock_no);
#		member_layers_up <- subset(member_layers,member_layers$formation_no_lo!=rock_no);
#		rock_no_lowest <- member_layers_up$rock_no_sr_up[1];

#		if (is.na(rock_no_lowest))	{
#			rock_no_lowest <- rock_no_highest;
#			} else	{
#			rock_no_highest <- rock_no_lowest;
#			}
#		} else	{
#		if (nrow(subset(section,section$formation_no_up==rock_no))>0)	{
#			subsection <- subset(section,section$formation_no_up==rock_no);
#			subsection$rock_no_sr_up[1]
#			}
#		subset(section,section$formation_no_lo==rock_no);
#		}
	}

rock_order_1 <- accersi_rock_order_in_section(rock_no=rock_no_lowest,section);
rock_order_2 <- accersi_rock_order_in_section(rock_no=rock_no_highest,section);
rock_info <- match(rock_order_1,section$rock_no_sr_up);
rock_order <- data.frame(rock_no=as.numeric(rock_order_1),
						 rock_unit=as.character(section$full_name_up[rock_info]),
						 ma_lb=as.numeric(section$up_ma_lb[rock_info]),
						 ma_ub=as.numeric(section$up_ma_ub[rock_info]),stringsAsFactors = F);
nad <- (1:nrow(rock_order))[is.na(rock_order$rock_unit)];
na_buster <- match(rock_order$rock_no[nad],section$rock_no_sr_lo);
rock_order$rock_unit[nad] <- section$full_name_lo[na_buster];
rock_order$ma_lb[nad] <- section$lo_ma_lb[na_buster];
rock_order$ma_ub[nad] <- section$lo_ma_ub[na_buster];

this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
rock_order <- subset(rock_order,rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
seq_rocks <- nrow(rock_order);
rownames(rock_order) <- 1:seq_rocks;
this_rock_lb_init <- rock_order$ma_lb[this_rock_lo];
this_rock_ub_init <- rock_order$ma_ub[this_rock_hi];

# older rocks must end by the latest date of younger rocks #
#rock_no    rock_unit ma_lb ma_ub			rock_no    rock_unit ma_lb ma_ub
#1    7379 Hongjingshao 519.4 515.0		1    7379 Hongjingshao 519.4 515.0
#2    2001   Qiongzhusi 519.4 516.2		2    2001   Qiongzhusi 519.4 516.2
#3    7381   Zhujiaqing 541.0 521.0		3    7381   Zhujiaqing 541.0 521.0
#4    3494     Dengying 550.0 517.6		4    3494     Dengying 550.0 521.0

for (ro in 1:(nrow(rock_order)-1))	{
	older_rocks <- (ro+1):seq_rocks;
	underlapping <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_lb[ro]];
	rock_order$ma_ub[underlapping[rock_order$ma_ub[underlapping]<rock_order$ma_ub[ro]]] <- rock_order$ma_ub[ro];
	}

ro <- this_rock_hi;
while (ro < nrow(rock_order))	{
	older_rocks <- (ro+1):seq_rocks;
	constrained <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_ub[ro]];
	rock_order$ma_ub[constrained] <- rock_order$ma_ub[ro];
	ro <- ro+1;
	}

# younger rocks cannot start before the onset date of older rocks #
ro <- this_rock_lo;
while (ro > 1)	{
	younger_rocks <- (ro-1):1;
	restrained <- younger_rocks[rock_order$ma_lb[younger_rocks]>rock_order$ma_lb[ro]];
	rock_order$ma_lb[restrained] <- rock_order$ma_lb[ro];
#	constrained <- younger_rocks[rock_order$ma_ub[younger_rocks]>rock_order$ma_ub[ro]];
	ro <- ro-1;
	}

this_rock_lb <- rock_order$ma_lb[this_rock_lo];
this_rock_ub <- rock_order$ma_ub[this_rock_hi];
rock_order <- subset(rock_order,rock_order$ma_lb>this_rock_ub);
rock_order <- subset(rock_order,rock_order$ma_ub<this_rock_lb);
this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);

#subset(rock_order,!as.logical((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<rock_order$ma_ub[this_rock_hi])));
num_overlapping_higher_rocks <- sum((rock_order$ma_lb<rock_order$ma_lb[this_rock_hi]) * (rock_order$ma_ub<=rock_order$ma_ub[this_rock_hi]));
num_overlapping_lower_rocks <- sum((rock_order$ma_lb>=rock_order$ma_lb[this_rock_lo]) * (rock_order$ma_ub>rock_order$ma_ub[this_rock_hi]));

# put in part about assessing formations with members #
seq_rocks <- nrow(rock_order);
ttl_events <- seq_rocks+1;
if (rock_order$ma_ub[1]<rock_order$ma_ub[this_rock_hi])
	ttl_events <- ttl_events-1;
if (rock_order$ma_lb[seq_rocks]>rock_order$ma_lb[this_rock_lo])
	ttl_events <- ttl_events-1;

this_rock_events <- c(1+seq_rocks-this_rock_lo,2+seq_rocks-this_rock_hi);
events_prior_start <- this_rock_events[1]-1;
events_prior_end <- this_rock_events[2]-1;
events_after_start <- (ttl_events - events_prior_start)-1;
events_after_end <- (ttl_events - events_prior_end)-1;
events_intervening <- events_prior_end-events_prior_start-1;

# fractional overlap:
# rock_no	rock_unit	ma_lb	ma_ub
# 1000		Shire		521.0	515.4
# 1001		Buckland	523.8	517.6
# set up situation where we have half an event from 521.0 onwards
# Buckland cannot stop until 521 or later
# Buckland cannot start after 517.6 or later
if (num_overlapping_higher_rocks>0)	{
	upper_breaks <- rock_order$ma_lb[rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]];
	upper_breaks <- sort(upper_breaks[upper_breaks<rock_order$ma_lb[this_rock_hi]]);
#	upper_breaks <- c(upper_breaks,rock_order$ma_ub[this_rock_hi])
	} else	{
	upper_breaks <- rock_order$ma_ub[this_rock_hi];
	}
if (num_overlapping_lower_rocks>0)	{
	lower_breaks <- rock_order$ma_ub[rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]];
	lower_breaks <- sort(lower_breaks[lower_breaks>rock_order$ma_ub[this_rock_lo]]);
#	lower_breaks <- sort(c(lower_breaks,rock_order$ma_lb[this_rock_lo]));
	} else	{
	lower_breaks <- rock_order$ma_lb[this_rock_lo];
	}

time_slices <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision);
ttl_time_slices <- length(time_slices);
if (min(upper_breaks)==rock_order$ma_ub[this_rock_hi] && max(lower_breaks)==rock_order$ma_lb[this_rock_lo])	{
	pp <- (time_slices-this_rock_ub)/(this_rock_lb-this_rock_ub);
	p_event <- dpois(1,ttl_events/length(time_slices));
	p_start <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_start,sbN=events_after_start);
#	p_start <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_start,sbN=events_after_start);
	p_start <- p_start/sum(p_start);
	p_end <- sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,prN=events_prior_end,sbN=events_after_end);
#	p_end <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,prN=events_prior_end,sbN=events_after_end);
	p_end_very <- dpois(1,ttl_events/(2*length(time_slices)));
	p_end <- p_end/sum(c(p_end,p_end_very));

	pt <- pp;
	p_end_given_start <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
	time_slices <- time_slices[ttl_time_slices:1];
	colnames(p_end_given_start) <- rownames(p_end_given_start) <- time_slices;
	p_end_cond <- p_end_given_start;
	p_event <- ttl_events/ttl_time_slices;
	for (i in 1:ttl_time_slices)	{
		p_end_given_start[i,] <- sapply(pt,accersi_p_end_given_start_time,ps=pp[i],p_event,intervening_events = events_intervening,subsequent_events = events_after_end,total_events = ttl_events);
		p_end_given_start[i,] <- p_end_given_start[i,]/sum(p_end_given_start[i,]);
		p_end_cond[i,] <- p_start[i]*p_end_given_start[i,];
		}
	} else if (min(upper_breaks)>rock_order$ma_ub[this_rock_hi] &&
			   max(lower_breaks)<rock_order$ma_lb[this_rock_lo])	{
	### case where there are overlapping possible durations for both older & younger rocks
#	time_slices <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision);
	ts_parts <- ts <- c();
	all_breaks <- unique(c(this_rock_lb,lower_breaks,upper_breaks,this_rock_ub));
	br <- 1;
	while (br < length(all_breaks))	{
		ts_parts <- c(ts_parts,length(seq(all_breaks[br+1]+(temporal_precision/2),all_breaks[br]-(temporal_precision/2),by=temporal_precision)));
		ts <- c(ts,rep(1/br,ts_parts[br]));
		br <- br+1;
		}
	ts <- ts/sum(ts);
	unique_ts <- unique(ts);
	ts_parts <- c(0,ts_parts);
	p_event_partitioned <- c();
	br <- 1;
	while (br < length(all_breaks))	{
		p_event_partitioned <- c(p_event_partitioned,dpois(1,unique_ts[br]));
		br <- br+1;
		}
	p_start <- p_end <- c();
	br <- 1;
	while (br < length(all_breaks))	{
		pp <- (this_rock_lb-seq(all_breaks[br+1]+(temporal_precision/2),all_breaks[br]-(temporal_precision/2),by=temporal_precision)) / (this_rock_lb-this_rock_ub);
		events_after_start <- sum(rock_order$ma_lb==this_rock_lb);	# note: this counts the base of the formation, but this doubles for the top of the formation
		p_start <- c(p_start,sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event=p_event_partitioned[br],prN=events_prior_start,sbN=events_after_start));
#		p_start <- c(p_start,sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,p_event=p_event_partitioned[br],prN=events_prior_start,sbN=events_after_start));
		events_prior_end <- events_after_start+events_intervening;
		events_after_end <- sum(rock_order$ma_ub==this_rock_ub)-1;	# note: this counts the base of the formation, but this doubles for the top of the formation
		p_end <- c(p_end,sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event=p_event_partitioned[br],prN=events_prior_end,sbN=events_after_end));
#		p_end <- c(p_end,sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,p_event=p_event_partitioned[br],prN=events_prior_end,sbN=events_after_end));
		br <- br+1;
		}
	p_start <- p_start/sum(p_start);
	p_end_very <- dpois(1,p_event_partitioned[ub]/2);
	p_end <- p_end/sum(c(p_end,p_end_very));

	pt <- pp;
	p_end_given_start <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
	time_slices <- time_slices[ttl_time_slices:1];
	colnames(p_end_given_start) <- rownames(p_end_given_start) <- time_slices;
	p_end_cond <- p_end_given_start;
	p_event <- ttl_events/ttl_time_slices;
	for (i in 1:ttl_time_slices)	{
		#p_end_given_start[i,] <- sapply(pt,accersi_p_end_given_start_time,ps=pp[i],p_event,intervening_events = events_intervening,subsequent_events = events_after_end,total_events = ttl_events);
		#p_end_given_start[i,] <- p_end_given_start[i,]/sum(p_end_given_start[i,]);
		p_end_given_start[i,i:ttl_time_slices] <- p_end[i:ttl_time_slices]/sum(p_end[i:ttl_time_slices]);
		p_end_cond[i,] <- p_start[i]*p_end_given_start[i,];
		}
	} else if (min(upper_breaks)>rock_order$ma_ub[this_rock_hi])	{
	### case where there are overlapping possible durations for younger rocks only
#	time_slices <- seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision);
	ts_parts <- ts <- c();
	all_breaks <- unique(c(this_rock_lb,upper_breaks,this_rock_ub));
	ub <- 1;
	while (ub < length(all_breaks))	{
		ts_prt <- length(seq(all_breaks[ub+1]+(temporal_precision/2),all_breaks[ub]-(temporal_precision/2),by=temporal_precision));
		ts_parts <- c(ts_parts,ts_prt);
		ts <- c(ts,rep(1/ub,ts_parts[ub]));
		ub <- ub+1;
		}
	ts <- ts/sum(ts);
	p_event_partitioned <- dpois(1,ts[cumsum(ts_parts)]);

	p_start <- p_end <- c();
	ub <- 1;
#	for (ub in 1:length(upper_breaks))	{
	while (ub < length(all_breaks))	{
		pp <- (this_rock_lb-seq(all_breaks[ub+1]+(temporal_precision/2),all_breaks[ub]-(temporal_precision/2),by=temporal_precision)) / (this_rock_lb-this_rock_ub)
		events_after_start <- sum(rock_order$ma_lb==this_rock_lb);	# note: this counts the base of the formation, but this doubles for the top of the formation
		p_start <- c(p_start,sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event=p_event_partitioned[ub],prN=events_prior_start,sbN=events_after_start));
#		p_start <- c(p_start,sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,p_event=p_event_partitioned[ub],prN=events_prior_start,sbN=events_after_start));
		events_prior_end <- events_after_start+events_intervening;
		events_after_end <- sum(rock_order$ma_ub==this_rock_ub)-1;	# note: this counts the base of the formation, but this doubles for the top of the formation
		p_end <- c(p_end,sapply(pp,exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event=p_event_partitioned[ub],prN=events_prior_end,sbN=events_after_end));
#		p_end <- c(p_end,sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous,p_event=p_event_partitioned[ub],prN=events_prior_end,sbN=events_after_end));
		ub <- ub+1;
		}
	p_start <- p_start/sum(p_start);
	p_end_very <- dpois(1,p_event_partitioned[length(p_event_partitioned)]/2);
	p_end <- p_end/sum(c(p_end,p_end_very));

	pt <- pp;
	p_end_given_start <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
	time_slices <- time_slices[ttl_time_slices:1];
	colnames(p_end_given_start) <- rownames(p_end_given_start) <- time_slices;
	p_end_cond <- p_end_given_start;
#	p_event <- ttl_events/ttl_time_slices;
	for (i in 1:ttl_time_slices)	{
		#p_end_given_start[i,] <- sapply(pt,accersi_p_end_given_start_time,ps=pp[i],p_event,intervening_events = events_intervening,subsequent_events = events_after_end,total_events = ttl_events);
		#p_end_given_start[i,] <- p_end_given_start[i,]/sum(p_end_given_start[i,]);
		p_end_given_start[i,i:ttl_time_slices] <- p_end[i:ttl_time_slices]/sum(p_end[i:ttl_time_slices]);
		p_end_cond[i,] <- p_start[i]*p_end_given_start[i,];
		}
	}
#pg1 <- accersi_p_fa_given_finds_p_start_and_p_end(finds=1,p_onset=p_start,p_end=p_end)
#pg2 <- accersi_p_fa_given_finds_p_start_and_p_end(finds=2,p_onset=p_start,p_end=p_end)
#pg3 <- accersi_p_fa_given_finds_p_start_and_p_end(finds=3,p_onset=p_start,p_end=p_end)
#plot(pp,pg3/sum(pg3),pch=21,bg="red")
#points(pp,pg2/sum(pg2),pch=21,bg="orange")
#points(pp,pg1/sum(pg1),pch=21,bg="yellow")

p_started_cum <- cumsum(p_start);
p_not_ended_cum <- 1-cumsum(p_end);
p_present <- p_started_cum*p_not_ended_cum;

output <- list(data.frame(poss_age=as.numeric(time_slices),
						  prob_rock_present=as.numeric(p_present),
						  prob_start=as.numeric(p_start),
						  prob_end=as.numeric(p_end),stringsAsFactors = F),
#			   lnp_end_given_start_matrix,
			   p_end_given_start);
#names(output) <- c("age_probs","lnp_end_given_start","p_end_given_start");
names(output) <- c("age_probs","p_end_given_start");

return(output)
}

#### beta distribution with alpha = 1 & beta = finds nested within beta distribution with alpha=rock_position-1 and beta=ttl_rocks-1;
prob_fa_given_finds_and_rock_position <- function(finds,rock_position,ttl_rocks,poss_onset,poss_end,slc=100)	{
# finds: number of finds in this rock;
# rock_position: rock position in sequence within this time interval, from bottom to top
# ttl_rocks: number of rocks in sequence within this time interval
# poss_onset: the oldest possible age for the rock;
# poss_end: the youngest possible age for the rock;
# slc: number of slices for whole interval (slc=100 = precision=0.01)
poss_onset <- abs(poss_onset);
poss_end <- abs(poss_end);
tspan <- poss_onset-poss_end;
dates <- seq(poss_onset-tspan/slc,poss_end+tspan/slc,by=-tspan/slc);
pp <- (poss_onset-dates)/(poss_onset-poss_end);
p_rock_start <- dbeta(pp,rock_position-1,1+ttl_rocks-rock_position)/sum(dbeta(pp,rock_position-1,1+ttl_rocks-rock_position));
remaining <- ttl_rocks-rock_position;
nn <- length(p_rock_start);
p_rock_cond <- array(0,dim=c(nn,nn));
p_fa_range <- list();
for (pr in 1:nn)	{
	parts <- pr+1;
	ppp <- seq(1/parts,(parts-1)/parts,by=1/parts);
	p_fa_range <- rlist::list.append(p_fa_range,dbeta(ppp,1,finds)/sum(dbeta(ppp,1,finds)));
	}

p_fa_cond <- array(0,dim=c(nn,nn));
for (pr in 1:nn)	{
	parts <- 2+nn-pr;
	ppp <- seq(1/parts,(parts-1)/parts,by=1/parts);
	p_rock_end_cond <- dbeta(ppp,1,remaining)/sum(dbeta(ppp,1,remaining));
	p_rock_cond[pr,pr:nn] <- p_rock_start[pr] * p_rock_end_cond;

	for (prr in 1:(1+nn-pr))
		p_fa_cond[pr,pr:(pr+prr-1)] <- p_fa_cond[pr,pr:(pr+prr-1)]+p_rock_end_cond[prr]*p_fa_range[[prr]]
	p_fa_cond[pr,] <- p_fa_cond[pr,]*p_rock_start[pr];
	}
p_fa_this_rock <- colSums(p_fa_cond);
names(p_fa_this_rock) <- dates;
return(p_fa_this_rock)
}

#### beta distribution with alpha = prn+1 & beta = sbN+1;
# prN <- 0; sbN <- 4;
prob_distribution_event_at_x_given_n_prior_and_m_later_events <- function(prN,sbN,span=1,precision=0.01)	{
# beta distribution with shape1 (alpha) = prN+1 and shape2 (beta) = sbN+1
# mode is prN/(prN+sbN); so, 2nd of 4 events expected around 0.33;
pp <- seq(precision,span,by=precision)/span;
prob_dist <- dbeta(pp,prN+1,sbN+1)/sum(dbeta(pp,prN+1,sbN+1));
names(prob_dist) <- pp;
return(prob_dist);
}

exact_prob_event_at_x_given_n_prior_and_m_later_events <- function(pp,prN,sbN,span=1,precision=0.01)	{
# prN: number of prior events;
# sbN: number of subsequent events;
# beta distribution with shape1 (alpha) = prN+1 and shape2 (beta) = sbN+1
# mode is prN/(prN+sbN); so, 2nd of 4 events expected around 0.33;
# pp: proportion of the way through the interval;
ppp <- seq(precision,span,by=precision)/span;
exact_prob <- dbeta(pp,prN+1,sbN+1)/sum(dbeta(ppp,prN+1,sbN+1));
return(exact_prob);
}

# given N events over T=1, what is the probability of the nth event happening at t=pp;
#accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete <- function(pp,n,p_event,prior,later,N)	{
accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete <- function(pp,p_event,prN,sbN)	{
## p_event: exact probability of an event at point x
## pp: proportion of time elapse at x (= x / total duration)
## prN: number of prior events
## sbN: number of subsequent events
if (pp[1]<1)	{
	return(p_event*(pp^prN)*((1-pp)^sbN));
	}	else	return(1);
}

# given N events over T=1, what is the probability of the nth event happening at t=pp;
accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete_old <- function(pp,p_event,prior,later)	{
## p_event: exact probability of an event at point x
## pp: proportion of time elapse at x (= x / total duration)
## younger: number of prior events
## number of subsequent events
if (pp[1]<1)	{
	return(p_event*(pp^prior)*((1-pp)^later));
	}	else	return(1);
}

# given N events over T=1, what is the probability of the nth event happening at t=pp;
#accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous <- function(pp,inst_rate,prN,sbN,span=1)	{
accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous <- function(pp,prN,sbN,span=1)	{
## pevent: exact probability of an event at point x
## pp: proportion of time elapse at x (= x / total duration)
## prN: number of prior events
## sbN: number of subsequent events
## span: time (or whatever) over which 1+prN+sbN events should happen;
inst_rate <- (1+prN+sbN)/span;
#if (is.null(span))
#	span <- 1+prN+sbN;
if (pp[1]<1)	{
#	return(exp(log(p_event)+log(dpois(prior,total*pp[1]))));
	# p[event] x p[earlier events over pp[1]] x p[later events over 1-pp[1]]
#	p_event*dpois(prior,total*pp)*dpois(later,total*(1-pp))
#	sum(exp(log(p_event)+log(dpois(prior,total*pp))+log(dpois(later,total*(1-pp)))))
# P[event now] x P[prN events | inst_rate*pp*span] x P[sbN | inst_rate*(1-pp)*span]
#	return(exp(log(inst_rate)+
#			   log(dpois(prN,span*pp))+
#			   log(dpois(sbN,span*(1-pp)))));
	return(exp(log(dpois(prN,inst_rate*pp*span))+log(dpois(sbN,inst_rate*(1-pp)*span))))
	}	else	return(1);
}

prob_true_fa_given_multiple_sections <- function(p_spc_onsets)	{
	sec <- nrow(p_spc_onsets)
	tslices <- ncol(p_spc_onsets)
	p_spc_ttl <- colSums(p_spc_onsets)
	col_span <- (1:tslices)[p_spc_ttl>0]
	cum_prob_fnd <- p_spc_onsets
	cum_prob_not_yet <- array(1,dim=c(sec,tslices))
	for (cl in 1:sec)	{
		cum_prob_fnd[cl,] <- cumsum(p_spc_onsets[cl,])
		cum_prob_not_yet[cl,] <- 1-cum_prob_fnd[cl,]
		cum_prob_not_yet[cl,(1:tslices)[cum_prob_not_yet[cl,]<0]] <- 0
	}
	if (max(col_span)<(tslices-1))
		for (ts in (max(col_span)+1):tslices)
			for (cl in 1:sec)	cum_prob_not_yet[cl,ts] <- 0
	#prob_not_yet <- colProd(cum_prob_not_yet)
	prob_not_yet <- exp(colSums(log(cum_prob_not_yet)));
	prob_first_appearance <- array(0,dim=tslices)
	for (ts in col_span)	{
		if (ts>1)	{
			prob_first_appearance[ts] <- prob_not_yet[ts-1]-prob_not_yet[ts]
		} else	{
			prob_first_appearance[ts] <- 1-prob_not_yet[ts]
		}
	}
	return(prob_first_appearance)
}

### this should be rewritten to get the average numbers of rocks below and average numbers of rocks above
###   that handles situations where a rock unit is the top of one section and/or the bottom of another.
prob_true_la_given_multiple_sections <- function(p_spc_termini)	{
	sec <- nrow(p_spc_termini)
	tslices <- ncol(p_spc_termini)
	p_spc_ttl <- colSums(p_spc_termini)
	col_span <- (1:tslices)[p_spc_ttl>0]
	cum_prob_fnd <- p_spc_termini
	cum_prob_not_yet <- array(1,dim=c(sec,tslices))
	for (cl in 1:sec)	{
		#	cum_prob_fnd[cl,] <- sum_vector_descending(p_spc_termini[cl,]);
		cum_prob_fnd[cl,] <- cumsum(p_spc_termini[cl,ncol(p_spc_termini):1]);
		#	cum_prob_fnd[cl,] <- sum_vector_descending(p_spc_termini[cl,]);
		cum_prob_not_yet[cl,] <- 1-cum_prob_fnd[cl,]
		cum_prob_not_yet[cl,(1:tslices)[cum_prob_not_yet[cl,]<0]] <- 0
	}
	for (ts in 1:col_span[1])	for (cl in 1:sec)	cum_prob_not_yet[cl,ts] <- 0
	#cum_prob_fnd[,col_span]
	#cum_prob_not_yet[,col_span]
	#prob_not_yet <- colProd(cum_prob_not_yet)
	prob_not_yet <- exp(colSums(log(cum_prob_not_yet)))
	prob_last_appearance <- array(0,dim=tslices)
	col_span <- sort(col_span,decreasing = TRUE)
	for (ts in col_span)	prob_last_appearance[ts] <- prob_not_yet[ts+1]-prob_not_yet[ts]
	return(prob_last_appearance)
}

### use this to get P[formation present at time t]
#events_prior_start <- this_rock_events[1]-1;
#events_prior_end <- this_rock_events[2]-1;
#intervening_events <- (events_prior_end-events_prior_start)-1;
#events_after_start <- (ttl_events - events_prior_start)-1;
#subsequent_events <- events_after_end <- (ttl_events - events_prior_end)-1;
#total_events <- ttl_events;
#ps <- 0.5
#pt <- 0.75
accersi_p_end_given_start_time <- function(pt,ps,p_event,intervening_events=0,subsequent_events,total_events)	{
if (pt>ps)	{
	return(exp(log(p_event)+log(dpois(intervening_events,total_events*(pt-ps)))+log(dpois(subsequent_events,total_events*(1-pt)))));
#	p_event*dpois(intervening_events,total_events*(pt-ps))*dpois(subsequent_events,total_events*(1-pt))
#	return(exp(log(dpois(intervening_events,total_events*(pt-ps)))+log(p_event)+log(dpois(subsequent_events,total_events*(1-pt)))));
#	return(p_event*((pt-ps)^intervening_events)*((1-pt)^subsequent_events));
#	p_event*(pp^prior)*((1-pp)^later)
#	p_event*dpois(intervening_events,total_events*(pt-ps))
	} else if (pt==ps)	{
	pt <- 1-((1-pt)/2)
	return(exp(log(p_event)+log(dpois(intervening_events,total_events*(pt-ps)))+log(dpois(subsequent_events,total_events*(1-pt)))));
	} else	{
	return(exp(log(p_event)+log(dpois(intervening_events,total_events*(pt-ps)))+log(dpois(subsequent_events,total_events*(1-pt)))));
	}
}

### use this to get P[formation present at time t]
accersi_p_present_given_p_start_and_p_end <- function(p_onset,p_end)	{
cum_p_onset <- cumsum(p_onset);
cum_p_end <- cumsum(p_end);
cpe <- c(0,cum_p_end[1:(length(cum_p_end)-1)]);
cum_p_still <- 1-cpe;
return(cum_p_onset-cpe);
}

### use this to get P[fa at time t]
accersi_p_fa_given_finds_p_start_and_p_end <- function(finds,p_onset,p_end)	{
prior <- 0;
later <- finds-1;
prob_fa <- array(0,dim=length(p_onset));
for (s in 1:length(p_onset))	{
	#	for (t in 1:length(p_end))	{
	t <- s;
	while (t<=length(p_end))	{
		p_event <- 1/length(s:t);
		end <- 1+t-s;
		pp <- seq(1/(2*end),1,by=1/end);
#		pfa <- accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete(pp,p_event,prior,later);
		pfa <- exact_prob_event_at_x_given_n_prior_and_m_later_events(pp,prN=prior,sbN=later);
#		pfa <- accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous(pp,p_event,prior,later);
		prob_fa[s:t] <- prob_fa[s:t] + (p_onset[s]*p_end[t]*pfa);
		t <- t+1;
		}
	}
prob_fa <- prob_fa/sum(prob_fa);
return(prob_fa)
}

### use this to get P[la at time t]
accersi_p_la_given_finds_p_start_and_p_end <- function(finds,p_onset,p_end)	{
prior <- finds-1;
later <- 0;
prob_la <- array(0,dim=length(p_onset));
for (s in 1:length(p_onset))	{
		#	for (t in 1:length(p_end))	{
	t <- s;
	while (t<=length(p_end))	{
		p_event <- 1/length(s:t)
		end <- 1+t-s;
		pp <- seq(1/(2*end),1,by=1/end);
#		p_event <- dpois(1,1)
		pla <- accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous(pp,prN=prior,sbN=later);
#		pla <- accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_continuous(pp,p_event,prior,later);
#			pla <- accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete(pp,p_event,prior,later)
		prob_la[s:t] <- prob_la[s:t] + (p_onset[s]*p_end[t]*pla)
		t <- t+1;
		}
	}
prob_la <- prob_la/sum(prob_la)
return(prob_la)
}

# 2020-04-22: routine to eliminate redundant sections introduced
rock_fuzzy_fa_and_la_given_superpositions_full <- function(relv_rock_summary,rock_superposition,rock_to_zones,radiometric_dates,rock_to_isotope_excursion,formations_from_members=T,temporal_precision=0.1)	{
ttl_rocks <- nrow(relv_rock_summary);
rock_need_help <- c();
age_probs_total <- age_probs_ave <- start_probs_ave <- end_probs_ave <- list();
#lnp_end_given_start <- p_end_given_start <- age_probs_total <- list();
tr <- 0;
while (tr < ttl_rocks)	{
	tr <- tr+1;
	print(paste("doing",relv_rock_summary$rock_unit[tr]));
	this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr]));
	if (nrow(this_rock_sequences)==0)	{
		rock_need_help <- c(rock_need_help,relv_rock_summary$rock_unit[tr]);
		if (relv_rock_summary$status[tr]==0)	{
			relv_rock_summary$rock_no[tr] <- relv_rock_summary$formation_no[tr];
			this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr]));
			}
		if (nrow(this_rock_sequences)==0)	{
			this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr-1]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr-2]));
			this_rock_sequences$column <- "Hobbiton";
			this_rock_sequences <- data.frame(column=as.character(c("Hobbiton","Hobbiton")),
											  formation_no_up=as.numeric(c(0,relv_rock_summary$formation_no[tr])),
											  rock_no_sr_up=as.numeric(c(0,relv_rock_summary$rock_no[tr])),
											  rock_no_up=as.numeric(c(0,relv_rock_summary$rock_no[tr])),
											  full_name_up=as.character(c("",relv_rock_summary$rock_unit[tr])),
											  formation_no_lo=as.numeric(c(relv_rock_summary$formation_no[tr],0)),
											  rock_no_sr_lo=as.numeric(c(relv_rock_summary$rock_no[tr],0)),
											  rock_no_lo=as.numeric(c(relv_rock_summary$rock_no[tr],0)),
											  full_name_lo=as.character(c(relv_rock_summary$rock_unit[tr],"")),
											  stringsAsFactors = F);
			this_rock_sequences$formation_no_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$formation_no[tr];
			this_rock_sequences$full_name_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$rock_unit[tr];
			this_rock_sequences$rock_no_sr_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$rock_no[tr];
			fake_older <- match(rock_database$ma_lb[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)],rock_database$ma_ub);
			this_rock_sequences$formation_no_lo[2] <- rock_database$formation_no[fake_older];
			this_rock_sequences$rock_no_sr_lo[2] <- rock_database$rock_no_sr[fake_older];
			this_rock_sequences$rock_no_lo[2] <- rock_database$rock_no[fake_older];
			this_rock_sequences$full_name_lo[2] <- rock_database$full_name[fake_older];

			fake_younger <- match(rock_database$ma_ub[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)],rock_database$ma_lb);
			this_rock_sequences$formation_no_up[1] <- rock_database$formation_no[fake_younger];
			this_rock_sequences$rock_no_sr_up[1] <- rock_database$rock_no_sr[fake_younger];
			this_rock_sequences$rock_no_up[1] <- rock_database$rock_no[fake_younger];
			this_rock_sequences$full_name_up[1] <- rock_database$full_name[fake_younger];
			}
		}
	this_rock_columns <- unique(this_rock_sequences$column);

	# get span of possible ages of collections from rocks in these sections
	poss_age <- seq(rock_database$ma_lb[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)]-temporal_precision/2,
					rock_database$ma_ub[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)]+temporal_precision/2,by=-temporal_precision);

	ttl_sections <- length(this_rock_columns);
	age_probs_by_section <- array(0,dim=c(length(poss_age),1+ttl_sections));
	colnames(age_probs_by_section) <- c("poss_age",this_rock_columns);
	start_probs_by_section <- end_probs_by_section <- age_probs_by_section <- data.frame(age_probs_by_section,stringsAsFactors = F);
	start_probs_by_section$poss_age <- end_probs_by_section$poss_age <- age_probs_by_section$poss_age <- poss_age;
#	p_end_given_start_section <- lnp_end_given_start_section <- list();
#	for (trc in 1:ttl_sections)	{
	trc <- 0;
	while (trc < ttl_sections)	{
		trc <- trc+1;
		section <- subset(rock_superposition,rock_superposition$column==this_rock_columns[trc]);
		if (nrow(section)==0)
			section <- this_rock_sequences;
		section$up_ma_lb <- rock_database$ma_lb[match(section$rock_no_sr_up,rock_database$rock_no_sr)];
		section$up_ma_ub <- rock_database$ma_ub[match(section$rock_no_sr_up,rock_database$rock_no_sr)];
		section$lo_ma_lb <- rock_database$ma_lb[match(section$rock_no_sr_lo,rock_database$rock_no_sr)];
		section$lo_ma_ub <- rock_database$ma_ub[match(section$rock_no_sr_lo,rock_database$rock_no_sr)];
		#		if (sum(is.na(section$lo_ma_lb))>0)
		#			section$lo_ma_lb <- rock_database$ma_lb[match(section$full_name_lo[is.na(section$lo_ma_lb)],rock_database$full_name)]
		#		if (sum(is.na(section$lo_ma_ub))>0)
		#			section$lo_ma_ub <- rock_database$ma_ub[match(section$full_name_lo[is.na(section$lo_ma_ub)],rock_database$full_name)]
		if (relv_rock_summary$status[tr]==0)	{
			rock_no <- relv_rock_summary$rock_no[tr];
#			superposition_info <- accersi_relative_position_within_time_slice(rock_no=relv_rock_summary$rock_no[tr],section,formations_from_members=F,temporal_precision=temporal_precision);
#			                      accersi_relative_position_within_time_slice_full(rock_no,section,radiometric_dates,rock_to_isotope_excursion,formations_from_members=T,temporal_precision=0.1)
			age_probs <- accersi_relative_position_within_time_slice_full(rock_no,section,radiometric_dates,rock_to_isotope_excursion,formations_from_members=F,temporal_precision=temporal_precision);
			} else	{
			rock_no <- relv_rock_summary$rock_no[tr];
			age_probs <- accersi_relative_position_within_time_slice_full(rock_no,section,radiometric_dates,rock_to_isotope_excursion,formations_from_members=T,temporal_precision=temporal_precision);
#			superposition_info <- accersi_relative_position_within_time_slice(rock_no=relv_rock_summary$rock_no[tr],section,formations_from_members=formations_from_members,temporal_precision=temporal_precision);
			}
		# subset superposition_info to include just the relevant time frame
#		age_probs <- superposition_info$age_probs;
		age_probs <- subset(age_probs,round(age_probs$poss_age,3) %in% round(age_probs_by_section$poss_age,3));
		nnn <- match(round(age_probs$poss_age,3),round(age_probs_by_section$poss_age,3));
		age_probs_by_section[nnn,trc+1] <- age_probs$prob_rock_present;
		nnn <- match(round(age_probs$poss_age,3),round(start_probs_by_section$poss_age,3));
		start_probs_by_section[nnn,trc+1] <- age_probs$prob_start;
		nnn <- match(round(age_probs$poss_age,3),round(end_probs_by_section$poss_age,3));
		end_probs_by_section[nnn,trc+1] <- age_probs$prob_end;

#		p_end_given_start_section <- rlist::list.append(p_end_given_start_section,superposition_info$p_end_given_start);
#		logged_info <- expello_na_from_matrix(log(superposition_info$p_end_given_start),0);
#		lnp_end_given_start_section <- rlist::list.append(lnp_end_given_start_section,logged_info);
		#		lnp_end_given_start_section <- rlist::list.append(lnp_end_given_start_section,superposition_info$lnp_end_given_start);
		}
#	names(p_end_given_start_section) <- names(lnp_end_given_start_section) <- this_rock_columns;
	age_probs_total <- rlist::list.append(age_probs_total,age_probs_by_section);
#	p_end_given_start <- rlist::list.append(p_end_given_start,p_end_given_start_section);
#	lnp_end_given_start <- rlist::list.append(lnp_end_given_start,lnp_end_given_start_section);
	if (ttl_sections==1)	{
		colnames(age_probs_by_section) <- colnames(start_probs_by_section) <- colnames(end_probs_by_section) <- c("poss_age","prob_present");
		age_probs_ave <- rlist::list.append(age_probs_ave,age_probs_by_section);
		start_probs_ave <- rlist::list.append(start_probs_ave,start_probs_by_section);
		end_probs_ave <- rlist::list.append(end_probs_ave,end_probs_by_section);
		#		p_end_given_start <- rlist::list.append(p_end_given_start,as.array(p_end_given_start_section));
		#		lnp_end_given_start <- rlist::list.append(lnp_end_given_start,as.array(lnp_end_given_start_section));
		} else	{
		# START HERE 2020-03-20 make sure that all of the sections span the same age range
		redundant_sections <- c();
		for (cn in 2:ttl_sections)	{
			for (nc in (cn+1):(ttl_sections+1))	{
				if (sum(age_probs_by_section[,cn]==age_probs_by_section[,nc])==nrow(age_probs_by_section))	{
					redundant_sections <- unique(c(redundant_sections,nc));
					}
				}
			}
		keepers <- (1:ncol(age_probs_by_section))[!(1:ncol(age_probs_by_section)) %in% redundant_sections];
		age_probs_by_section <- age_probs_by_section[,keepers];
		section_lengths <- c();
		for (trc in 1:ttl_sections)
			section_lengths <- rbind(section_lengths,c(max(age_probs$poss_age),min(age_probs$poss_age)));
#			section_lengths <- rbind(section_lengths,c(max(as.numeric(rownames(p_end_given_start_section[[trc]]))),min(as.numeric(rownames(p_end_given_start_section[[trc]])))));

		start_probs <- data.frame(poss_age=as.numeric(start_probs_by_section$poss_age),
								  prob_present=as.numeric(1-exp(rowSums(log(1-start_probs_by_section[,(1:ttl_sections)+1])))));
		start_probs_ave <- rlist::list.append(start_probs_ave,start_probs);
		end_probs <- data.frame(poss_age=as.numeric(end_probs_by_section$poss_age),prob_present=as.numeric(1-exp(rowSums(log(1-end_probs_by_section[,(1:ttl_sections)+1])))));
		end_probs_ave <- rlist::list.append(end_probs_ave,end_probs);
		cuml_start_probs <- 1-exp(cumsum(log(1-start_probs$prob_present)));
		end_probs_cum <- sort((1-exp(cumsum(log(1-end_probs$prob_present[nrow(end_probs):1])))),decreasing = T);
		ave_probs <- data.frame(poss_age=as.numeric(end_probs_by_section$poss_age),prob_present=as.numeric(exp(log(cuml_start_probs)+log(end_probs_cum))));
		age_probs_ave <- rlist::list.append(age_probs_ave,ave_probs);
		} # end case where the sections do not quite line up.
	names(start_probs_ave) <- names(age_probs_total) <- relv_rock_summary$rock_unit[1:tr];
#	names(lnp_end_given_start) <- names(p_end_given_start) <- relv_rock_summary$rock_unit[1:tr];
	}

names(age_probs_total) <- names(age_probs_ave) <- names(start_probs_ave) <- names(end_probs_ave) <- relv_rock_summary$rock_unit;
#names(p_end_given_start) <- names(lnp_end_given_start) <- relv_rock_summary$rock_unit;

output <- list(age_probs_total,age_probs_ave,start_probs_ave,end_probs_ave);
names(output) <- c("age_probs_total","age_probs_ave","start_probs_ave","end_probs_ave");
return(output);
}

rock_fuzzy_fa_and_la_given_superpositions <- function(relv_rock_summary,rock_superposition,rock_to_zones,formations_from_members=T,temporal_precision=0.1)	{
ttl_rocks <- nrow(relv_rock_summary);
rock_need_help <- c();
lnp_end_given_start <- p_end_given_start <- age_probs_total <- age_probs_ave <- start_probs_ave <- end_probs_ave <- list();
#lnp_end_given_start <- p_end_given_start <- age_probs_total <- list();
tr <- 0;
while (tr < ttl_rocks)	{
	tr <- tr+1;
	print(paste("doing",relv_rock_summary$rock_unit[tr]));
	this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr]));
	if (nrow(this_rock_sequences)==0)	{
		rock_need_help <- c(rock_need_help,relv_rock_summary$rock_unit[tr]);
		if (relv_rock_summary$status[tr]==0)	{
			relv_rock_summary$rock_no[tr] <- relv_rock_summary$formation_no[tr];
			this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr]));
			}
		if (nrow(this_rock_sequences)==0)	{
			this_rock_sequences <- rbind(subset(rock_superposition,rock_superposition$rock_no_sr_up==relv_rock_summary$rock_no[tr-1]),subset(rock_superposition,rock_superposition$rock_no_sr_lo==relv_rock_summary$rock_no[tr-2]));
			this_rock_sequences$column <- "Hobbiton";
			this_rock_sequences <- data.frame(column=as.character(c("Hobbiton","Hobbiton")),
											  formation_no_up=as.numeric(c(0,relv_rock_summary$formation_no[tr])),
											  rock_no_sr_up=as.numeric(c(0,relv_rock_summary$rock_no[tr])),
											  rock_no_up=as.numeric(c(0,relv_rock_summary$rock_no[tr])),
											  full_name_up=as.character(c("",relv_rock_summary$rock_unit[tr])),
											  formation_no_lo=as.numeric(c(relv_rock_summary$formation_no[tr],0)),
											  rock_no_sr_lo=as.numeric(c(relv_rock_summary$rock_no[tr],0)),
											  rock_no_lo=as.numeric(c(relv_rock_summary$rock_no[tr],0)),
											  full_name_lo=as.character(c(relv_rock_summary$rock_unit[tr],"")),
											  stringsAsFactors = F);
			this_rock_sequences$formation_no_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$formation_no[tr];
			this_rock_sequences$full_name_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$rock_unit[tr];
			this_rock_sequences$rock_no_sr_up[this_rock_sequences$rock_no_up==relv_rock_summary$rock_no[tr-1]] <- relv_rock_summary$rock_no[tr];
			fake_older <- match(rock_database$ma_lb[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)],rock_database$ma_ub);
			this_rock_sequences$formation_no_lo[2] <- rock_database$formation_no[fake_older];
			this_rock_sequences$rock_no_sr_lo[2] <- rock_database$rock_no_sr[fake_older];
			this_rock_sequences$rock_no_lo[2] <- rock_database$rock_no[fake_older];
			this_rock_sequences$full_name_lo[2] <- rock_database$full_name[fake_older];

			fake_younger <- match(rock_database$ma_ub[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)],rock_database$ma_lb);
			this_rock_sequences$formation_no_up[1] <- rock_database$formation_no[fake_younger];
			this_rock_sequences$rock_no_sr_up[1] <- rock_database$rock_no_sr[fake_younger];
			this_rock_sequences$rock_no_up[1] <- rock_database$rock_no[fake_younger];
			this_rock_sequences$full_name_up[1] <- rock_database$full_name[fake_younger];
			}
		}
	this_rock_columns <- unique(this_rock_sequences$column);

	poss_age <- seq(rock_database$ma_lb[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)]-temporal_precision/2,
					rock_database$ma_ub[match(relv_rock_summary$rock_no[tr],rock_database$rock_no)]+temporal_precision/2,by=-temporal_precision);

	ttl_sections <- length(this_rock_columns);
	age_probs_by_section <- array(0,dim=c(length(poss_age),1+ttl_sections));
	colnames(age_probs_by_section) <- c("poss_age",this_rock_columns);
	start_probs_by_section <- end_probs_by_section <- age_probs_by_section <- data.frame(age_probs_by_section,stringsAsFactors = F);
	start_probs_by_section$poss_age <- end_probs_by_section$poss_age <- age_probs_by_section$poss_age <- poss_age;
	p_end_given_start_section <- lnp_end_given_start_section <- list();
	for (trc in 1:ttl_sections)	{
		section <- subset(rock_superposition,rock_superposition$column==this_rock_columns[trc]);
		if (nrow(section)==0)
			section <- this_rock_sequences;
		section$up_ma_lb <- rock_database$ma_lb[match(section$rock_no_sr_up,rock_database$rock_no_sr)];
		section$up_ma_ub <- rock_database$ma_ub[match(section$rock_no_sr_up,rock_database$rock_no_sr)];
		section$lo_ma_lb <- rock_database$ma_lb[match(section$rock_no_sr_lo,rock_database$rock_no_sr)];
		section$lo_ma_ub <- rock_database$ma_ub[match(section$rock_no_sr_lo,rock_database$rock_no_sr)];
#		if (sum(is.na(section$lo_ma_lb))>0)
#			section$lo_ma_lb <- rock_database$ma_lb[match(section$full_name_lo[is.na(section$lo_ma_lb)],rock_database$full_name)]
#		if (sum(is.na(section$lo_ma_ub))>0)
#			section$lo_ma_ub <- rock_database$ma_ub[match(section$full_name_lo[is.na(section$lo_ma_ub)],rock_database$full_name)]
		if (relv_rock_summary$status[tr]==0)	{
#			superposition_info <- accersi_relative_position_within_time_slice(rock_no=relv_rock_summary$rock_no[tr],section,formations_from_members=F,temporal_precision=temporal_precision);
			superposition_info <- accersi_relative_position_within_time_slice_full(rock_no=relv_rock_summary$rock_no[tr],section,formations_from_members=F,temporal_precision=temporal_precision);
			} else	{
			superposition_info <- accersi_relative_position_within_time_slice(rock_no=relv_rock_summary$rock_no[tr],section,formations_from_members=formations_from_members,temporal_precision=temporal_precision);
			}
		# subset superposition_info to include just the relevant time frame
		age_probs <- superposition_info$age_probs;
		age_probs <- subset(age_probs,round(age_probs$poss_age,3) %in% round(age_probs_by_section$poss_age,3));
		nnn <- match(round(age_probs$poss_age,3),round(age_probs_by_section$poss_age,3));
		age_probs_by_section[nnn,trc+1] <- age_probs$prob_rock_present;
		start_probs_by_section[nnn,trc+1] <- age_probs$prob_start;
		end_probs_by_section[nnn,trc+1] <- age_probs$prob_end;

		p_end_given_start_section <- rlist::list.append(p_end_given_start_section,superposition_info$p_end_given_start);
		logged_info <- expello_na_from_matrix(log(superposition_info$p_end_given_start),0);
		lnp_end_given_start_section <- rlist::list.append(lnp_end_given_start_section,logged_info);
#		lnp_end_given_start_section <- rlist::list.append(lnp_end_given_start_section,superposition_info$lnp_end_given_start);
		}
	names(p_end_given_start_section) <- names(lnp_end_given_start_section) <- this_rock_columns;
	age_probs_total <- rlist::list.append(age_probs_total,age_probs_by_section);
	p_end_given_start <- rlist::list.append(p_end_given_start,p_end_given_start_section);
	lnp_end_given_start <- rlist::list.append(lnp_end_given_start,lnp_end_given_start_section);
	if (ttl_sections==1)	{
		colnames(age_probs_by_section) <- colnames(start_probs_by_section) <- colnames(end_probs_by_section) <- c("poss_age","prob_present");
		age_probs_ave <- rlist::list.append(age_probs_ave,age_probs_by_section);
		start_probs_ave <- rlist::list.append(start_probs_ave,start_probs_by_section);
		end_probs_ave <- rlist::list.append(end_probs_ave,end_probs_by_section);
#		p_end_given_start <- rlist::list.append(p_end_given_start,as.array(p_end_given_start_section));
#		lnp_end_given_start <- rlist::list.append(lnp_end_given_start,as.array(lnp_end_given_start_section));
		} else	{
		# START HERE 2020-03-20 make sure that all of the sections span the same age range
		section_lengths <- c();
		for (trc in 1:ttl_sections)
			section_lengths <- rbind(section_lengths,c(max(as.numeric(rownames(p_end_given_start_section[[trc]]))),min(as.numeric(rownames(p_end_given_start_section[[trc]])))));

		start_probs <- data.frame(poss_age=as.numeric(start_probs_by_section$poss_age),prob_present=as.numeric(1-exp(rowSums(log(1-start_probs_by_section[,(1:ttl_sections)+1])))));
		start_probs_ave <- rlist::list.append(start_probs_ave,start_probs);
		end_probs <- data.frame(poss_age=as.numeric(end_probs_by_section$poss_age),prob_present=as.numeric(1-exp(rowSums(log(1-end_probs_by_section[,(1:ttl_sections)+1])))));
		end_probs_ave <- rlist::list.append(end_probs_ave,end_probs);
		start_probs_cum <- 1-exp(cumsum(log(1-start_probs$prob_present)));
		end_probs_cum <- sort((1-exp(cumsum(log(1-end_probs$prob_present[nrow(end_probs):1])))),decreasing = T);
		ave_probs <- data.frame(poss_age=as.numeric(end_probs_by_section$poss_age),prob_present=as.numeric(exp(log(start_probs_cum)+log(end_probs_cum))));
		age_probs_ave <- rlist::list.append(age_probs_ave,ave_probs);
		} # end case where the sections do not quite line up.
	names(start_probs_ave) <- names(lnp_end_given_start) <- names(p_end_given_start) <- names(age_probs_total) <- relv_rock_summary$rock_unit[1:tr];
	}

names(age_probs_total) <- names(age_probs_ave) <- names(start_probs_ave) <- relv_rock_summary$rock_unit;
names(end_probs_ave) <- names(p_end_given_start) <- names(lnp_end_given_start) <- relv_rock_summary$rock_unit;

output <- list(age_probs_total,age_probs_ave,start_probs_ave,
			   end_probs_ave,p_end_given_start,lnp_end_given_start);
names(output) <- c("age_probs_total","age_probs_ave","start_probs_ave",
			   "end_probs_ave","p_end_given_start","lnp_end_given_start");
return(output);
}

detritus <- function () {
# work down #
#workdown <- subset(section,section$rock_no_sr_up==rock_no_lowest);
#if (nrow(workdown)>0)	{
#	this_rock_lb <- workdown$up_ma_lb;
#	this_rock_ub <- workdown$up_ma_ub;
#	new_workdown <- subset(section,section$rock_no_sr_up==workdown$rock_no_sr_lo[nrow(workdown)]);
#	while (workdown$lo_ma_ub[nrow(workdown)]==this_rock_ub && nrow(new_workdown)>0)	{
#		next_rock <- workdown$rock_no_sr_lo[nrow(workdown)];
#		new_workdown <- subset(section,section$rock_no_sr_up==next_rock);
#		if (nrow(new_workdown)>1)
#			new_workdown <- subset(new_workdown,new_workdown$formation_no_lo==new_workdown$rock_no_sr_lo);
#		if (nrow(new_workdown)==1)	{
#			workdown <- rbind(workdown,new_workdown);
#			} else if (nrow(new_workdown)>1)	{
#			print(paste(relv_rock_units[tr],"has too many underlying rocks"));
#			}
#		}
#	rock_order_down <- rbind(data.frame(rock_no=as.numeric(workdown$rock_no_sr_up[1]),
#										rock_unit=as.character(workdown$full_name_up[1]),
#										ma_lb=as.numeric(workdown$up_ma_lb[1]),
#										ma_ub=as.numeric(workdown$up_ma_ub[1]),stringsAsFactors = F),
#							 data.frame(rock_no=as.numeric(workdown$rock_no_sr_lo),
#							 		   rock_unit=as.character(workdown$full_name_lo),
#							 		   ma_lb=as.numeric(workdown$lo_ma_lb),
#							 		   ma_ub=as.numeric(workdown$lo_ma_ub),stringsAsFactors = F));
#			workdown_all <- rlist::list.append(workdown_all,rock_order_down);
##			columns_down <- c(columns_down,this_rock_columns[trc]);
#	} else	{
#	rock_order_down <- "";
#	}

## work up #
#workup <- subset(section,section$rock_no_sr_lo==rock_no_highest);
#remainder_overlying <- c();
#if (nrow(workup)>0)	{
#	this_rock_lb <- workup$lo_ma_lb;
#	this_rock_ub <- workup$lo_ma_ub;
#	new_workup <- subset(section,section$rock_no_sr_lo==workup$rock_no_sr_up[nrow(workup)]);
#	while (workup$up_ma_lb[1]==this_rock_lb && nrow(new_workup)>0)	{
#		next_rock <- workup$rock_no_sr_up[1];
#		new_workup <- subset(section,section$rock_no_sr_lo==next_rock);
#		if (nrow(new_workup)==0 && !is.null(remainder_overlying))	{
#			workup[1,] <- remainder_overlying[1,];
#			next_rock <- workup$rock_no_sr_up[1];
#			new_workup <- subset(section,section$rock_no_sr_lo==next_rock);
#			if (nrow(remainder_overlying)>1)	{
#				remainder_overlying <- remainder_overlying[2:nrow(remainder_overlying),];
#				} else	{
#				remainder_overlying <- c();
##				}
#			}
#		if (nrow(new_workup)>1)
#			new_workup <- subset(new_workup,new_workup$formation_no_up!=new_workup$rock_no_sr_up);
#		if (nrow(new_workup)==1)	{
#			workup <- rbind(new_workup,workup);
#			} else if (nrow(new_workup)>1)	{
#			print(paste(new_workup$full_name_up,"has too many overlying rocks"));
#			remainder_overlying <- new_workup[2:nrow(new_workup),];
#			workup <- rbind(new_workup[1,],workup);
#			}
#		}
#	rock_order_up <- rbind(data.frame(rock_no=as.numeric(workup$rock_no_sr_up),
#									  rock_unit=as.character(workup$full_name_up),
#									  ma_lb=as.numeric(workup$up_ma_lb),
#									  ma_ub=as.numeric(workup$up_ma_ub),stringsAsFactors = F),
#						   data.frame(rock_no=as.numeric(workup$rock_no_sr_lo[nrow(workup)]),
#						   		   rock_unit=as.character(workup$full_name_lo[nrow(workup)]),
#						   		   ma_lb=as.numeric(workup$lo_ma_lb[nrow(workup)]),
#						   		   ma_ub=as.numeric(workup$lo_ma_ub[nrow(workup)]),stringsAsFactors = F));
#			#workup_all <- rlist::list.append(workup_all,rock_order_up);
#			#columns_up <- c(columns_up,this_rock_columns[trc]);
#	} else	{	# end tallying overlying rocks
#	rock_order_up <- "";
#	}

#if (is.data.frame(rock_order_up) && is.data.frame(rock_order_down))	{
#	if (nrow(rock_order_down)>1)	{
#		rock_order <- unique(rbind(rock_order_up,rock_order_down));
#		rownames(rock_order) <- 1:nrow(rock_order);
#		} else	{
#		rock_order <- rock_order_up;
#		}
#	} else if (is.data.frame(rock_order_up))	{
#	rock_order <- rock_order_up;
#	} else if (is.data.frame(rock_order_down)) {
#	rock_order <- rock_order_down;
#	} else	{
#	rock_order <- "";
#	}
#this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
#rock_order <- subset(rock_order,rock_order$ma_lb>rock_order$ma_ub[this_rock_hi]);
#this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
#rock_order <- subset(rock_order,rock_order$ma_ub<rock_order$ma_lb[this_rock_lo]);
#this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
#this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);
#seq_rocks <- nrow(rock_order);
#rownames(rock_order) <- 1:seq_rocks;

# look for later appearing rocks that start *after* this rock; in this section, this rock must be older than that.
# if there are multiple rocks that might or might not be in the interval, then make calculations under all combos
# older rocks must end by the latest date of younger rocks #
#ro <- this_rock_hi;
#while (ro > 1)	{
#	ro <- ro-1;
#	older_rocks <- (ro+1):seq_rocks;
#	constrained <- older_rocks[rock_order$ma_ub[older_rocks]<rock_order$ma_ub[ro]];
#	rock_order$ma_ub[constrained] <- rock_order$ma_ub[ro];
#	}
# younger rocks must start by the onset date of older rocks #
#ro <- this_rock_lo;
#while (ro < seq_rocks)	{
#	ro <- ro+1;
#	younger_rocks <- (ro-1):1;
#	restrained <- younger_rocks[rock_order$ma_lb[younger_rocks]>rock_order$ma_lb[ro]];
#	rock_order$ma_lb[restrained] <- rock_order$ma_lb[ro];
#	}

#this_rock_lb_init <- this_rock_lb;
#this_rock_ub_init <- this_rock_ub;
#this_rock_lb <- rock_order$ma_lb[this_rock_lo];
#this_rock_ub <- rock_order$ma_ub[this_rock_hi];
#rock_order <- subset(rock_order,rock_order$ma_lb>this_rock_ub);
#rock_order <- subset(rock_order,rock_order$ma_ub<this_rock_lb);
#this_rock_hi <- match(rock_no_highest,rock_order$rock_no);
#this_rock_lo <- match(rock_no_lowest,rock_order$rock_no);

#ttl_events <- seq_rocks+1;
#if (rock_order$ma_ub[1]<rock_order$ma_ub[this_rock_hi])
#	ttl_events <- ttl_events-1;
#if (rock_order$ma_lb[seq_rocks]>rock_order$ma_lb[this_rock_lo])
#	ttl_events <- ttl_events-1;

#pp <- time_slices <- (seq(this_rock_ub+(temporal_precision/2),this_rock_lb-(temporal_precision/2),by=temporal_precision)-this_rock_ub)/(this_rock_lb-this_rock_ub);
#p_event <- 1/length(time_slices);
#xxx <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event,prN=events_prior_start,sbN=events_after_start);
#yyy <- sapply(pp,accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events,p_event,prN=events_prior_end,sbN=events_after_end);
#plot(time_slices,xxx/sum(xxx),pch=21,bg="red")
#points(time_slices,yyy/sum(yyy),pch=21,bg="blue")
#plot(time_slices,xxx,pch=21,bg="red")
#points(time_slices,yyy,pch=21,bg="blue")

#ttl_time_slices <- length(time_slices);
#lnp_start <- lnp_end <- p_end <- vector(length=ttl_time_slices);
#lnp_end_given_start_matrix <- array(0,dim=c(ttl_time_slices,ttl_time_slices));
#colnames(lnp_end_given_start_matrix) <- rownames(lnp_end_given_start_matrix) <- this_rock_lb-time_slices*(this_rock_lb-this_rock_ub);
#p_now <-ttl_events/ttl_time_slices;
#lnp_now <- log(p_now);
#p_end_given_start_matrix <- lnp_end_given_start_matrix;
#for (ts in 1:ttl_time_slices)	{
#	expected_prior_events <- ttl_events*time_slices[ts];
#	expected_subsequent_events <- ttl_events*(1-time_slices[ts]);
#	lnp_start[ts] <- lnp_now+log(dpois(events_prior_start,expected_prior_events));
#	lnp_end[te] <- lnp_now+log(dpois(events_prior_end,expected_prior_events));
#	lnp_end[ts] <- lnp_now+log(dpois(events_after_end,expected_subsequent_events));
#	subsequent_rate <- events_after_start/(1-time_slices[ts]);
#	for (te in ts:ttl_time_slices)	{
		## rescale time slices by start & get probs of last appearances using accersi_exact_prob_event_at_x_given_n_prior_and_m_later_events_discrete
#		expected_intervening_events <- ttl_events*(time_slices[te]-time_slices[ts]);
#		expected_subsequent_events_dep <- subsequent_rate*(1-time_slices[te]);
#		lnp_end[te] <- log(p_end[te] + exp(lnp_start[ts]+lnp_now+log(dpois(events_after_end,expected_subsequent_events))));
#		lnp_end_given_start_matrix[ts,te] <- lnp_now+log(dpois(0,expected_intervening_events));
#		p_end[te] <- exp(lnp_end[te]);
#		}
#	p_end_given_start_matrix[ts,ts:ttl_time_slices] <- exp(lnp_end_given_start_matrix[ts,ts:ttl_time_slices]);
#	}
#p_start <- exp(lnp_start)/sum(exp(lnp_start));
#points(time_slices,p_start,pch=22,bg="green");
#p_end <- exp(lnp_end)/sum(exp(lnp_end));
#points(time_slices,p_end,pch=22,bg="yellow");
# look for earlier appearing rocks that end *before* this rock; in this section, this rock must be younger than that
#plot(time_slices,p_start);
#points(time_slices,p_end);
#p_started_cum <- cumsum(p_start);
#p_not_ended_cum <- 1-cumsum(p_end);
#p_present <- p_started_cum*p_not_ended_cum;

#par(pin=c(3,3));
#plot(time_slices,p_started_cum,pch=21,bg="red",ylim=c(0,1),cex=0.75,lwd=0.5,col="red4");
#points(time_slices,p_not_ended_cum,pch=21,bg="blue1",cex=0.75,lwd=0.5,col="blue4");
#points(time_slices,p_present,pch=21,bg="purple1",cex=0.75,lwd=0.5,col="purple4");

#output <- list(data.frame(poss_age=as.numeric(time_slices)),prob_rock_present=as.numeric(p_present),prob_start=as.numeric(p_start),
#					 prob_end=as.numeric(p_end),stringsAsFactors = F),
#			   lnp_end_given_start_matrix,
#			   p_end_given_start_matrix);
#names(output) <- c("age_probs","lnp_end_given_start","p_end_given_start");

#return(output)
}

#### Routines for Taxon Ranges ####
# lbs <- c(433.4,431.9,432.0,358.9,355.3); ubs <- c(430.5,430.5,345.0,346.7,348.9);
find_stratigraphic_range_gaps <- function(lbs,ubs)	{
# lbs: vector of oldest possible ages;
# ubs: vector of youngest possible ages;
nranges <- length(lbs);
lbsx <- sort(abs(lbs),decreasing=T)[2:nranges];
ubsx <- sort(abs(ubs),decreasing=T)[1:(nranges-1)];
#cbind(ubsx,lbsx)
gaps <- ubsx-lbsx;
if (!is.null(names(lbsx)))
	names(gaps) <- paste(names(lbsx),"-",names(ubsx),sep="");
return(gaps);
}

accersi_fuzzy_taxon_ranges_from_paleodb_and_rock_database <- function(taxon_list,finds_per_rock,age_probs_ave,all_time)	{
notu <- length(taxon_list);
prob_present <- array(0,dim=c(notu,length(all_time)));
colnames(prob_present) <- round(all_time,3);
rownames(prob_present) <- taxon_list;
nslices <- length(all_time);
ttl_rocks <- ncol(finds_per_rock);
tl <- 0;
p_las <- p_fas <- list();
p_pre_fa <- p_pre_la <- list();
p_obs_range <- good_taxa <- c();
while (tl < notu)	{
	tl <- tl+1;
	relv_finds <- finds_per_rock[tl,finds_per_rock[tl,]>0];
	if (length(relv_finds)>0)	{
		print(paste("Doing",taxon_list[tl]));
		good_taxa <- c(good_taxa,taxon_list[tl]);
		relv_rocks <- (1:ttl_rocks)[finds_per_rock[tl,]>0];
		occupied_rocks <- colnames(finds_per_rock)[finds_per_rock[tl,]>0];
		occupied <- length(relv_finds);
#		strat_age <- p_fa <- p_la <- p_range <- list();
		p_fa <- p_la <- p_range <- array(0,dim=c(length(all_time),occupied));
		names(p_fa) <- names(p_la) <- names(p_range) <- all_time;
		for (oc in 1:occupied)	{
			p_onset_all <- p_end_all <- vector(length=length(all_time));
			this_rock <- relv_rocks[oc];
			### shouldn't this be p_rock_present???
			rock_ages <- age_probs_ave[[this_rock]]$poss_age;
			relv_cells <- match(round(rock_ages,3),round(all_time,3));
			p_onset_all[relv_cells] <- age_probs_ave[[this_rock]]$prob_present;
			p_end_all[relv_cells] <- end_probs_ave[[this_rock]]$prob_present;
			p_fa[relv_cells,oc] <- accersi_p_fa_given_finds_p_start_and_p_end(finds=relv_finds[oc],p_onset=p_onset_all[relv_cells],p_end=p_end_all[relv_cells]);
			p_la[relv_cells,oc] <- accersi_p_la_given_finds_p_start_and_p_end(finds=relv_finds[oc],p_onset=p_onset_all[relv_cells],p_end=p_end_all[relv_cells]);
			}
		p_fas <- rlist::list.append(p_fas,p_fa);
		p_las <- rlist::list.append(p_las,p_la);
		p_appeared <- p_disappeared <- array(0,dim=c(length(all_time),occupied));
		for (oc in 1:occupied)	{
			p_appeared[,oc] <- cumsum(p_fa[,oc]);
			p_disappeared[,oc] <- cumsum(p_la[,oc]);
			}

		p_not_appeared <- 1-p_appeared;
		p_not_disappeared <- 1-p_disappeared;
		p_already_sampled <- 1-rowProds(p_not_appeared);
		p_will_yet_be_sampled <- 1-rowProds(p_disappeared);
		p_strat_range <- round(p_already_sampled*p_will_yet_be_sampled,15);
		rownames(p_not_appeared) <- rownames(p_not_disappeared) <- all_time;
#		p_obs_range <- rlist::list.append(p_obs_range,p_strat_range);
		p_pre_fa <- rlist::list.append(p_pre_fa,p_not_appeared);
		p_pre_la <- rlist::list.append(p_pre_la,p_not_disappeared);
		p_obs_range <- cbind(p_obs_range,p_strat_range);
#		p_pre_fa <- cbind(p_pre_fa,p_not_appeared);
#		p_pre_la <- cbind(p_pre_la,p_not_disappeared);
		}
	}
names(p_fas) <- names(p_las) <- names(p_pre_fa) <- names(p_pre_la) <- colnames(p_obs_range) <- good_taxa;
rownames(p_obs_range) <- all_time;
keepers <- min((1:nrow(p_obs_range))[rowSums(p_obs_range)>0]):nrow(p_obs_range);
p_obs_range <- p_obs_range[keepers,];
new_all_time <- as.numeric(rownames(p_obs_range));

for (gt in 1:length(good_taxa))	{
	p_fas[[gt]] <- p_fas[[gt]][keepers];
	p_las[[gt]] <- p_las[[gt]][keepers];
	p_pre_fa[[gt]] <- p_pre_fa[[gt]][keepers];
	p_pre_la[[gt]] <- p_pre_la[[gt]][keepers];
	}

output <- list(p_fas,p_las,p_pre_fa,p_pre_la,p_obs_range);
names(output) <- c("p_fas","p_las","p_pre_fa","p_pre_la","p_obs_range");
return(output);
}

#poss_ages <- data.frame(lbs=as.numeric(lbs),ubs=as.numeric(ubs),stringsAsFactors = F);
# first written 2020-07-03 for strophomenide divergence project
# returns the probability that the taxon is already present at some point in time.
#poss_ages <- data.frame(lbs=as.numeric(lbs),ubs=as.numeric(ubs),stringsAsFactors = F);
# first written 2020-07-03 for strophomenide divergence project
# returns the probability that the taxon is already present at some point in time.
p_already_first_appeared_given_multiple_finds <- function(poss_ages,prec=0.1)	{
age_sign <- poss_ages$lbs[1]/abs(poss_ages$lbs[1]);
poss_ages <- abs(poss_ages);
poss_ages <- poss_ages[order(-poss_ages$lbs,-poss_ages$ubs),];
oldest_poss <- max(poss_ages$lbs);
youngest_poss <- max(poss_ages$ubs);
rounder <- 1+ceiling(abs(log10(prec/2)));
poss_fa <- round(seq(oldest_poss-(prec/2),youngest_poss,by=-prec),rounder);
poss_fas_ttl <- round(seq(oldest_poss-(prec/2),min(poss_ages$ubs),by=-prec),rounder);
fa_probs <- array(0,dim=c(nrow(poss_ages),length(poss_fas_ttl)));
colnames(fa_probs) <- round(poss_fas_ttl,rounder);
for (n in 1:nrow(poss_ages))	{
	this_rock_poss <- round(seq(poss_ages$lbs[n]-(prec/2),poss_ages$ubs[n],by=-prec),rounder);
	fa_probs[n,match(this_rock_poss,poss_fas_ttl)] <- (poss_ages$lbs[n]-this_rock_poss)/(poss_ages$lbs[n]-poss_ages$ubs[n]);
	done <- match(min(this_rock_poss),poss_fas_ttl);
	if (done < length(poss_fas_ttl))	fa_probs[n,(done+1):length(poss_fas_ttl)] <- 1;
	}
#relv_fa_probs <- fa_probs[,match(poss_fa,poss_fas_ttl)];
#not_fa_probs <- (1-relv_fa_probs)
#ln_not_fa_probs <- log(not_fa_probs);
#ln_prob_not_fa <- colSums(ln_not_fa_probs);
#prob_not_fa <- exp(ln_prob_not_fa);
#prob_fa_overall <- 1-prob_not_fa;
if (nrow(fa_probs)>1)	{
	prob_fa_overall <- 1-exp(colSums(log((1-fa_probs[,match(poss_fa,poss_fas_ttl)]))));
	} else	{
	prob_fa_overall <- 1-exp(log((1-fa_probs[,match(poss_fa,poss_fas_ttl)])));
	}
if (age_sign==-1)
	colnames(fa_probs) <- as.numeric(colnames(fa_probs))/-1;
return(prob_fa_overall);
}

# first written 2020-07-03 for strophomenide divergence project
# returns the probability that the taxon is still present at some point in time.
p_still_to_be_found_again <- function(poss_ages,prec=0.1)	{
age_sign <- poss_ages$lbs[1]/abs(poss_ages$lbs[1]);
poss_ages <- abs(poss_ages);
poss_ages <- poss_ages[order(-poss_ages$lbs,-poss_ages$ubs),];
oldest_poss <- min(poss_ages$lbs);
youngest_poss <- min(poss_ages$ubs);
rounder <- 1+ceiling(abs(log10(prec/2)));
poss_la <- round(seq(oldest_poss-(prec/2),youngest_poss,by=-prec),rounder);
poss_las_ttl <- round(seq(max(poss_ages$lbs)-(prec/2),min(poss_ages$ubs),by=-prec),rounder);
#poss_las_ttl <- round(seq(oldest_poss-(prec/2),min(poss_ages$ubs),by=-prec),rounder);
p_penultimate <- array(0,dim=c(nrow(poss_ages),1+length(poss_las_ttl)));
colnames(p_penultimate) <- round(c(max(poss_las_ttl)+prec,poss_las_ttl),rounder);
for (n in 1:nrow(poss_ages))	{
	# 0: event could not yet have happened;
	# 1: event must already have happened;
	this_rock_poss <- round(seq(poss_ages$lbs[n]-(prec/2),poss_ages$ubs[n],by=-prec),rounder);
	this_rock_span <- round(poss_ages$lbs[n]-poss_ages$ubs[n],rounder)/round(prec,rounder);
	p_penultimate[n,match(this_rock_poss,colnames(p_penultimate))] <- (this_rock_poss-poss_ages$ubs[n])/(poss_ages$lbs[n]-poss_ages$ubs[n]);
	p_penultimate[n,1] <- 1;
	done <- match(max(this_rock_poss),poss_las_ttl);
	if (done > 1)	p_penultimate[n,1:(done-1)] <- 1;
	}

poss_la <- round(c(max(poss_la)+prec,poss_la),rounder);
poss_las_ttl <- round(c(max(poss_las_ttl)+prec,poss_las_ttl),rounder);
if (nrow(p_penultimate)>1)	{
#	1-p_penultimate[,match(poss_la,poss_las_ttl)]
#	log(1-p_penultimate[,match(poss_la,poss_las_ttl)])
#	colSums(log(1-p_penultimate[,match(poss_la,poss_las_ttl)]))
	prob_not_la_overall <- 1-exp(colSums(log(1-p_penultimate[,match(poss_la,poss_las_ttl)])));
	} else if (length(poss_la)==1)	{
	prob_not_la_overall <- 1-exp(sum(log((1-p_penultimate[,match(poss_la,poss_las_ttl)]))));
	names(prob_not_la_overall) <- prob_not_la_overall;
	} else	{
	prob_not_la_overall <- 1-exp(log((1-p_penultimate[,match(poss_la,poss_las_ttl)])));
	}
if (age_sign==-1)
	colnames(prob_not_la_overall) <- as.numeric(colnames(prob_not_la_overall))/-1;
return(prob_not_la_overall);
}

# first written 2020-07-03 for strophomenide divergence project
# returns exact probabilities of first appearances.
p_fa_given_multiple_finds <- function(poss_ages,prec=0.1)	{
# need a loop: p that this is the find x p not already found in another sample for each find;
age_sign <- poss_ages$lbs[1]/abs(poss_ages$lbs[1]);
poss_ages <- abs(poss_ages);
poss_ages <- poss_ages[order(-poss_ages$lbs,-poss_ages$ubs),];
oldest_poss <- max(poss_ages$lbs);
youngest_poss <- max(poss_ages$ubs);
rounder <- 1+ceiling(abs(log10(prec/2)));
poss_fa <- round(seq(oldest_poss-(prec/2),youngest_poss,by=-prec),rounder);
poss_fas_ttl <- round(seq(oldest_poss-(prec/2),min(poss_ages$ubs),by=-prec),rounder);
fa_prob_this_find <- fa_prob_already_found <- array(0,dim=c(nrow(poss_ages),length(poss_fas_ttl)));
colnames(fa_prob_this_find) <- colnames(fa_prob_already_found) <- round(poss_fas_ttl,rounder);
for (n in 1:nrow(poss_ages))	{
	this_rock_poss <- round(seq(poss_ages$lbs[n]-(prec/2),poss_ages$ubs[n],by=-prec),rounder);
	relv_bins <- match(this_rock_poss,poss_fas_ttl);
	fa_prob_this_find[n,relv_bins] <- 1/length(this_rock_poss);
#	fa_prob_this_find[n,relv_bins[1]] <- 0.5/length(this_rock_poss);
	fa_prob_already_found[n,2:(ncol(fa_prob_this_find))] <- cumsum(fa_prob_this_find[n,1:(ncol(fa_prob_this_find)-1)]);
	}
fa_prob_not_already_found <- 1-fa_prob_already_found;
p_fa <- array(0,dim=c(nrow(poss_ages),length(poss_fas_ttl)));
colnames(p_fa) <- round(poss_fas_ttl,rounder);
if (nrow(poss_ages)>1)	{
	for (n in 1:nrow(poss_ages))	{
		m <- (1:nrow(poss_ages))[!(1:nrow(poss_ages)) %in% n];
		if (length(m)>1)	{
			p_not_found_elsewhere <- exp(colSums(log(fa_prob_not_already_found[m,])));
			} else	{
			p_not_found_elsewhere <- fa_prob_not_already_found[m,];
			}
		p_fa[n,] <- exp(log(fa_prob_this_find[n,])+log(p_not_found_elsewhere));
		}
	prob_first <- colSums(p_fa);
	prob_first[!poss_fas_ttl %in% poss_fa] <- 0;
	prob_first <- prob_first/sum(prob_first);
	return(prob_first);
	} else	{
	return(fa_prob_this_find);
	}
}

# returns exact probabilities of first appearances.
#poss_ages <- data.frame(lbs=as.numeric(c(512,512,512,511,511)),ubs=as.numeric(c(510,510,510,509,509)));
p_fa_given_multiple_finds_beta_dist <- function(poss_ages)	{
# need a loop: p that this is the find x p not already found in another sample for each find;
age_sign <- poss_ages$lbs[1]/abs(poss_ages$lbs[1]);
poss_ages <- abs(poss_ages);
poss_ages <- poss_ages[order(-poss_ages$lbs,-poss_ages$ubs),];
p_a <- poss_ages[,colnames(poss_ages) %in% c("lbs","ubs")]
unq_poss_ages <- unique(p_a);
oldest_possible <- max(abs(poss_ages$lbs));
youngest_possible <- max(abs(poss_ages$ubs));
uncertainty <- oldest_possible - youngest_possible;
pp <- seq(0.01,0.99,length=99);
ppt <- round(max(abs(poss_ages$lbs))-pp*uncertainty,4);

unq_fa_probs <- unq_fa_cum_probs <- array(0,dim=c(nrow(unq_poss_ages),99));
colnames(unq_fa_probs) <- colnames(unq_fa_cum_probs) <- round(ppt,4);
for (pf in 1:nrow(unq_poss_ages))	{
	N <- sum(p_a$lbs==unq_poss_ages$lbs[pf] & p_a$ubs==unq_poss_ages$ubs[pf]);
	ppf <- ppt[ppt<unq_poss_ages$lbs[pf]];
	pfspan <- abs(unq_poss_ages$lbs[pf]-unq_poss_ages$ubs[pf]);
	ppff <- (abs(unq_poss_ages$lbs[pf])-ppf)/pfspan;
	unq_fa_probs[pf,match(round(ppf,4),ppt)] <- dbeta(ppff,1,N)*ppff[1];
	unq_fa_cum_probs[pf,] <- cumsum(unq_fa_probs[pf,]);
	}
## now, get P[first find in unique poss span] * P[not found in another span already]
p_not_found_yet <- 1-unq_fa_cum_probs;

p_fa_from_this_possibility <- unq_fa_probs;
for (pf in 1:nrow(unq_poss_ages))	{
	elsewhere <- (1:nrow(unq_poss_ages))[!(1:nrow(unq_poss_ages)) %in% pf]
	if (length(elsewhere)>1)	{
		ln_p_not_found_elsewhere <- colSums(log(p_not_found_yet[elsewhere,]))
		ln_p_not_found_elsewhere <- c(0,ln_p_not_found_elsewhere[1:(length(ln_p_not_found_elsewhere)-1)])
		} else if (length(elsewhere)==1)	{
		ln_p_not_found_elsewhere <- log(p_not_found_yet[elsewhere,]);
		ln_p_not_found_elsewhere <- c(0,ln_p_not_found_elsewhere[1:(length(ln_p_not_found_elsewhere)-1)])
		} else	{
		ln_p_not_found_elsewhere <- rep(0,ncol(p_not_found_yet));
		}
	p_fa_from_this_possibility[pf,] <- exp((log(unq_fa_probs[pf,])+ln_p_not_found_elsewhere));
	}

return(1-exp(colSums(log(1-p_fa_from_this_possibility))));
}

#### Routines to Ordinate Collections ####
tally_range_gaps <- function(comat)	{
#comat: ordered coocurrence matrix;
nspec <- nrow(comat);
gaps <- array(0,dim=nspec);
names(gaps) <- colnames(comat);
for (nn in 1:nspec)	{
	range <- (1:nspec)[comat[nn,]==1];
	gaps[nn] <- (1+max(range)-min(range))-length(range)
#	gaps[nn] <- sum(!(min(range):max(range)) %in% range);
	}
return(gaps)
}

cooccurrence_matrix_add_genera <- function(finds) {
gen_finds <- finds;
gen_finds$accepted_name <- gen_finds$genus;
gen_finds$accepted_no <- gen_finds$genus_no;
gen_finds$accepted_rank <- "genus";
gen_finds$accepted_rank[gen_finds$subgenus_no!=0] <- "subgenus";
finds <- rbind(finds,gen_finds);
finds <- finds[order(finds$collection_no),];
pres_abs_matrix <- generate_taxon_by_locality_matrix(finds);
coccr_mat <- crossprod(t(pres_abs_matrix));
coccr_mat[coccr_mat>1] <- 1;
return(coccr_mat);
}

get_cooccurrence_matrix <- function(pres_abs_matrix) {
# pres_abs_matrix: taxon x locality matrix giving finds
if (is.data.frame(pres_abs_matrix))
	pres_abs_matrix <- as.matrix(pres_abs_matrix);
coccr_mat <- crossprod(t(pres_abs_matrix));
coccr_mat[coccr_mat>1] <- 1;
return(coccr_mat);
}

find_unbroken_ranges <- function(comat)	{
#comat: ordered coocurrence matrix;
nspec <- nrow(comat);
unbroken <- array(0,dim=nspec)
for (nn in 1:nspec)	{
	range <- (1:nspec)[comat[nn,]==1];
	if (max(range)-min(range)==(length(range)-1))	unbroken[nn] <- 1;
	}
return(unbroken)
}

accersi_conjunct_pairs <- function(find_data)	{
conjunct_pairs <- c();
ncoll <- length(unique(find_data$collection_no));
collection_nos <- sort(unique(find_data$collection_no));
taxa <- sort(unique(find_data$accepted_name));

for (c in 1:ncoll)	{
	coll_finds <- subset(find_data,find_data$collection_no==collection_nos[c]);
	find_taxa <- sort(unique(coll_finds$accepted_name));
	find_taxa_no <- match(find_taxa,taxa);
	S <- length(find_taxa);
	if (S>1)	{
		for (n in 1:(S-1))	{
			for (m in (n+1):S)	{
				conjunct_pairs <- rbind(conjunct_pairs,c(find_taxa_no[n],find_taxa_no[m]));
				}
			}
		}
	}
return(conjunct_pairs)
}

disjunct_distance_ordination <- function(find_data)	{
ncoll <- length(unique(find_data$collection_no));
taxon_names <- sort(unique(find_data$accepted_name));
ntaxa <- length(taxon_names);

conjunct_pairs <- accersi_conjunct_pairs(find_data);

avg <- temp <- nc <- vector(length=ntaxa);
#comat <- comat2 <- vector(length=ncoll);
notinclust <- live <- grpsize <- grp <- vector(length=ntaxa);

nco <- nrow(conjunct_pairs);
nc <- vector(length=ntaxa);
for (n in 1:ntaxa)
	nc[n] <- sum(conjunct_pairs[,1]==n) + sum(conjunct_pairs[,2]==n);
ghost <- sum(nc==0);	# taxa without pairs

grp <- array(0,dim=ntaxa);
grpsize <- c();
lnlive <- ngrps <- ngoodgrps <- 0;
for (cp in 1:nco)	{
	x <- conjunct_pairs[cp,1];
	y <- conjunct_pairs[cp,2];
	if (grp[x]==0 && grp[y]==0)	{
		ngrps <- ngrps+1;
		ngoodgrps <- ngoodgrps+1;
		grp[x] <- grp[y] <- ngrps;
		grpsize <- c(grpsize,2);
		} else if (grp[x]==0)	{
		grp[x] <- grp[y];
		grpsize[grp[y]] <- grpsize[grp[y]]+1;
		} else if (grp[y]==0)	{
		grp[y] <- grp[x];
		grpsize[grp[x]] <- grpsize[grp[x]]+1;
		} else if (grp[x]!=grp[y])	{
		ngoodgrps <- ngoodgrps-1;
		q <- min(c(grp[x],grp[y]));
		r <- max(c(grp[x],grp[y]));
		grp[grp==r] <- q;
		grpsize[q] <- grpsize[q]+grpsize[r];
		grpsize[r] <- 0;
		}
	}
group_skree <- sort(grpsize[grpsize>0],decreasing=T);
good_groups <- (1:max(grp))[grpsize>0];
good_groups <- good_groups[order(grpsize[grpsize>0],decreasing = T)]
grp2 <- match(grp,good_groups);
grp2[is.na(grp2)] <- 0;

# if there is more than one good group, then let's focus on the biggest good group.
#if (ngoodgrps > 1)	{
#	x <- y <- 0;
#	}

## initialize variables
stotal <- 0;
hi <- ntaxa;
for (i in 1:ntaxa)	{
	lo <- i;
	avg[i] <- lo/hi;
	temp[i] <- 0;
	stotal <- stotal+avg[i];
	}
#if (length(live)>0)
#	for (i in 1:lnlive)
#		avg[live[i]]

hidif <- MAXNO;
iter<- lastdif <- 0;

while ((hidif>0.00001) && (lastdif!=hidif))	{
	iter <- iter+1;
	if (iter %% 10==0)
		print(paste("iteration: ",iter,"; hidif = ",hidif,"; ",timestamp(prefix="",suffix = "",quiet=T),sep=""));

	for (i in 1:nco)		{
		temp[conjunct_pairs[i,1]] <- temp[conjunct_pairs[i,1]]+avg[conjunct_pairs[i,2]];
		temp[conjunct_pairs[i,2]] <- temp[conjunct_pairs[i,2]]+avg[conjunct_pairs[i,1]];
#		temp[comat[i]] <- temp[comat[i]]+avg[comat2[i]];
#		temp[comat2[i]] <- temp[comat2[i]]+avg[comat[i]];
		}
	hi <- 0;
	lo < - MAXNO;
	for (i in 1:ntaxa)
		if (nc[i]>0)	{
			thang=nc[i]+1;
			temp[i]=(temp[i]+avg[i])/thang;
			if (temp[i]>hi)
				hi=temp[i];
			if (temp[i]<lo)
				lo=temp[i];
			}
	lastdif <- hidif;
	hidif <- stotal <- 0;
	i <- 1;
#	for (i in 1:lnlive;i++)
	while (i <lnlive)	{
		temp[live[i]]=hi;
		i <- i+1;
		}
	for (i in 1:ntaxa)	{
		if (nc[i]>0)	{
			temp[i] <- (temp[i]-lo)/(hi-lo);
			stotal <- stotal+temp[i];
			if (abs(avg[i]-temp[i])>hidif)
				hidif <- abs(avg[i]-temp[i]);
			avg[i] <- temp[i];
			temp[i] <- 0;
			}
		}
	}
names(avg) <- names(grp2) <- taxon_names;
output <- list(avg,group_skree,grp2);
names(output) <- c("ddo","group_skree","linked_group");
return(output);
}

accersi_edge_sequence <- function(ddo,reduced_list,gotliv=F,print_update=F)	{
if (print_update)
	print("Finding edge sequence...");
nedg <- ncimp <- 0;
ntaxa <- length(ddo);
nfinds <- nrow(reduced_list);
collection_nos <- sort(unique(reduced_list$collection_no));
rncoll <- length(collection_nos);
bot <- top <- id<- posit <- mscore <- vector(length=rncoll);
temp <- live <- vector(length=ntaxa);
mat <- vector(length=nfinds);
#mat=(int*)calloc(mknumrex+1,sizeof(int));

for (n in 1:rncoll)	{
	coll_finds <- subset(reduced_list,reduced_list$collection_no==collection_nos[n]);
	S <- sum(ddo[coll_finds$taxon_no]>0);
	mscore[n] <- sum(ddo[coll_finds$taxon_no])/S;
	}
ghosts <- sum(ddo==0);
id <- (1:rncoll)[order(mscore)];
mscore <- sort(mscore);
#plot(1:ncoll,mscore);
posit[id] <- 1:rncoll;
red_collection_info <- data.frame(collection_no=as.numeric(collection_nos),position=as.numeric(posit),stringsAsFactors = F);

incl_taxa <- sort(unique(reduced_list$taxon_no));
rntaxa <- length(incl_taxa);
fst <- lst <- vector(length=rntaxa);
for (tx in 1:rntaxa)	{
	txn <- incl_taxa[tx];
	taxon_finds <- subset(reduced_list,reduced_list$taxon_no==txn);
	if (nrow(taxon_finds)>0)	{
		taxon_colls <- match(taxon_finds$collection_no,collection_nos);
		fst[tx] <- min(posit[taxon_colls]);
		lst[tx] <- max(posit[taxon_colls]);
		}
	}
#for (nc in 1:rncoll)	{
#	coll_finds <- subset(reduced_list,reduced_list$collection_no==collection_nos[n]);
#	coll_taxa <- coll_finds$taxon_no;
#	S <- length(coll_taxa);
#	for (tx in 1:S)	{
#		if (fst[coll_taxa[tx]] > posit[nc])
#			fst[coll_taxa[tx]] <- posit[nc];
#		if (lst[coll_taxa[tx]] < posit[nc])
#			lst[coll_taxa[tx]] <- posit[nc];
#		}
#	}
taxon_position <- data.frame(fst=as.numeric(fst),lst=as.numeric(lst));
output <- list(red_collection_info,taxon_position);
names(output) <- c("ordinated_collections","taxon_ordination_ranges");
return(output);
}

constrained_ddo_optimization <- function(relv_collections)	{
bin_boundaries <- sort(unique(c(relv_collections$ma_lb,relv_collections$ma_ub)),decreasing=T);
relv_collections$faux_bin_lb <- match(relv_collections$ma_lb,bin_boundaries);
relv_collections$faux_bin_ub <- match(relv_collections$ma_ub,bin_boundaries)-1;
mn_bin <- min(relv_collections$faux_bin_lb);
mx_bin <- max(relv_collections$faux_bin_ub);
faux_bins <- unique(data.frame(bin_lb=as.numeric(relv_collections$faux_bin_lb),bin_ub=as.numeric(relv_collections$faux_bin_ub),stringsAsFactors = hell_no));
faux_bins <- faux_bins[order(faux_bins$bin_lb,faux_bins$bin_ub),]
bb <- mn_bin-1;
while (bb < mx_bin)	{
	bb <- bb+1;
	interval_name <- relv_collections$interval_lb[match(bb,relv_collections$faux_bin_lb)];
	if (is.na(interval_name))
		interval_name <- relv_collections$interval_ub[match(bb,relv_collections$faux_bin_ub)];
	poss_bin_colls <- subset(relv_collections,relv_collections$faux_bin_lb<=bb & relv_collections$faux_bin_ub>=bb);
	poss_colls <- nrow(poss_bin_colls);
	if (poss_colls>1)	{
		poss_bin_colls <- poss_bin_colls[order(poss_bin_colls$ddo_pos),];
		ages <- data.frame(lb=as.numeric(poss_bin_colls$faux_bin_lb),ub=as.numeric(poss_bin_colls$faux_bin_ub));
		ages <- constrain_dates_given_superposition(ages);
		poss_bin_colls$faux_bin_lb <- ages$lb;
		poss_bin_colls$faux_bin_ub <- ages$ub;
		ages <- data.frame(lb=as.numeric(poss_bin_colls$ma_lb),ub=as.numeric(poss_bin_colls$ma_ub));
		ages <- constrain_dates_given_superposition(ages,young_low = T);
		poss_bin_colls$ma_lb <- ages$lb;
		poss_bin_colls$ma_ub <- ages$ub;
#		print(poss_bin_colls[c(1,poss_colls),c(94:95,106:110)]);
		age <- poss_bin_colls$ma_lb;
		poss_bin_colls$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
		age <- poss_bin_colls$ma_ub;
		poss_bin_colls$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
		}
#	def_older_colls <- subset(relv_collections,relv_collections$faux_bin_ub<bb);
	older_colls <- sum(relv_collections$faux_bin_ub<bb);
	if (sum(poss_bin_colls$faux_bin_ub==bb)>0)
		poss_bin_colls$ddo_pos_adj <- older_colls+(1:nrow(poss_bin_colls));
	redone <- match(poss_bin_colls$collection_no,relv_collections$collection_no);
	relv_collections[redone,] <- poss_bin_colls;
#	cbind(poss_bin_colls$ddo_pos_adj,poss_bin_colls$faux_bin_lb,poss_bin_colls$faux_bin_ub)
	}
relv_collections$faux_bin_lb <- relv_collections$faux_bin_ub <- NULL;
return(relv_collections);
}

unitary_association_optimization <- function(cooccur_matrix,above_below_matrix,lump=T)	{
# above_below_matrix[i,j]=1 if i first appears above j and =-1 if j first appears above i.
#D <- pairwise_differences_discrete(cooccur_matrix);
#gD <- gower_transform(D);
#init_ord <- eigen(gD);
#first_eigenvector <- init_ord$vectors[,1];
#names(first_eigenvector) <- rownames(cooccur_matrix);
#first_eigenvector <- sort(first_eigenvector);
#f_eigen_order <- match(rownames(cooccur_matrix),names(first_eigenvector));
#write.csv(cooccur_matrix[f_eigen_order,f_eigen_order],"PCO_Ordination.csv");
#mds <- initial_ordination <- vegan::metaMDS(cooccur_matrix);
nn <- match(rownames(above_below_matrix),rownames(cooccur_matrix));
keepers <- (1:nrow(above_below_matrix))[!is.na(nn)];
above_below_matrix <- above_below_matrix[keepers,keepers];
nn <- match(rownames(cooccur_matrix),rownames(above_below_matrix));
above_below_matrix <- above_below_matrix[nn,nn];
nn <- cbind(rownames(cooccur_matrix),rownames(above_below_matrix),colnames(above_below_matrix));

nspec <- nrow(cooccur_matrix);
species_names <- rownames(cooccur_matrix);
print("Finding species pairs with identical co-occurrences");
identical_cooccurrences <- accersi_identical_cooccurrences(cooccur_matrix);
identical_cooccurrences <- unique(identical_cooccurrences);
unique_species_no <- 1:nrow(identical_cooccurrences);
for (sp in 1:nspec)	unique_species_no[sp] <- which(identical_cooccurrences==sp,arr.ind = T)[1,1];

fa_below <- above_below_matrix;
fa_below[fa_below==-1] <- 0;
fa_above <- above_below_matrix;
fa_above[fa_above==1] <- 0;
fa_above <- abs(fa_above);
fa_rank_a <- rowSums(fa_below);
fa_rank_z <- nspec-rowSums(fa_above);

initial_heights_A <- data.frame(taxon_no=as.numeric(1:nspec),clique_no=as.numeric(unique_species_no),cooccrs=rowSums(cooccur_matrix),
							  max_height=as.numeric(fa_rank_z),min_height=as.numeric(fa_rank_a));
#write.csv(initial_heights,paste(area,"_Initial_Rankings.csv",sep=""));
lumpings <- hist(unique_species_no,0:nspec,plot=F)$counts;
for (sp in 1:nspec)	{
	if (lumpings[sp]>1)	{
#		lumping_heights <- subset(initial_heights,initial_heights$clique_no==sp);
		initial_heights_A$max_height[initial_heights_A$clique_no==sp] <- max(initial_heights_A$max_height[initial_heights_A$clique_no==sp]);
		initial_heights_A$min_height[initial_heights_A$clique_no==sp] <- max(initial_heights_A$min_height[initial_heights_A$clique_no==sp]);
		}
	}
initial_heights_A <- initial_heights_A[order(initial_heights_A$min_height,initial_heights_A$max_height,initial_heights_A$clique_no),];
initial_heights_A$clique_no <- match(initial_heights_A$clique_no,unique(initial_heights_A$clique_no));
initial_heights_A <- initial_heights_A[order(initial_heights_A$clique_no,initial_heights_A$min_height,initial_heights_A$max_height),];
#write.csv(initial_heights,paste(area,"_Second_Rankings.csv",sep=""));

clique_info <- initial_heights_A;
clique_info$taxon_no <- clique_info$clique_no;
colnames(clique_info)[1:2] <- c("clique_no","clique_size");
clique_size <- hist(clique_info$clique_no,breaks=c(0:max(clique_info$clique_no)),plot=F)$counts;
clique_info$clique_size <- clique_size[match(clique_info$clique_no,unique(clique_info$clique_no))];
clique_info <- unique(clique_info);
good_cliques <- clique_info[clique_info$clique_size<clique_info$cooccrs,];

initial_heights <- initial_heights_A[initial_heights_A$clique_no %in% good_cliques$clique_no,];
cspec <- nrow(initial_heights);
initial_heights$clique_no <- match(initial_heights$clique_no,unique(initial_heights$clique_no));
good_cliques$clique_no <- match(good_cliques$clique_no,unique(good_cliques$clique_no));
taxon_order <- match(rownames(initial_heights),rownames(cooccur_matrix));
species_names_ordered <- rownames(initial_heights);
lcmat <- cooccur_matrix[taxon_order,taxon_order];
cooccur_matrix_best <- cooccur_matrix_better <- lcmat;

uspec <- nrow(good_cliques);
improvement <- old_best <- best_gaps <- last_best_gaps <- init_gaps <- sum(tally_range_gaps(lcmat));
progress <- round(uspec*(1:100)/100,0);		problems <- c();	attempt <- 1;
while (improvement > 100 && attempt < 5)	{
	order_kept <- cl <- ttl_comps <- 0;
	whatsupdoc <- c();
	while (cl < uspec & sum(is.na(match(rownames(initial_heights),species_names)))==0 & ttl_comps<cspec)	{
		cl <- cl+1;
		this_clique_no <- good_cliques$clique_no[cl];
		ttl_comps <- ttl_comps+1;
		if (!is.na(match(cl,progress)))	{
			if (match(cl,progress)<10)	{
				prop_done <- paste("0",match(cl,progress),"%",sep="");
				} else {
				prop_done <- paste(match(cl,progress),"%",sep="");
				}

			cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',paste("Iteration ",attempt,": ",prop_done,sep=""));
			flush.console()
			}
		clique <- (1:cspec)[initial_heights$clique_no %in% good_cliques$clique_no[cl]];
		clique_species <- rownames(initial_heights)[clique];
		cur_height_cl <- spc <- min(clique);
		clique_size <- length(clique);
		max_height_cl <- cspec;	min_height_cl <- 1;
		for (cls in 1:length(clique))	{
			# note: we are using all of the species names because fa_above & fa_below were not culled
			species_above <- species_names[fa_above[match(clique_species[cls],rownames(fa_above)),]>0];
			species_above_no <- match(species_above,rownames(initial_heights));
			species_above <- species_above[!is.na(species_above_no)];
			if (length(species_above)>0)
				max_height_cl <- min(max_height_cl,min(match(species_above,rownames(initial_heights))));
			species_below <- species_names[fa_below[match(clique_species[cls],rownames(fa_below)),]>0];
			species_below_no <- match(species_below,rownames(initial_heights));
			species_below <- species_below[!is.na(species_below_no)];
			if (length(species_below)>0)
				min_height_cl <- max(min_height_cl,max(match(species_below,rownames(initial_heights))));
			}
		ht <- max(clique)+1;
		cmat <- lcmat <- cooccur_matrix_best;
		clique_comp <- 0;
		if ((attempt %% 1)==0 || (attempt %% 1)==1)	{
			cooccr_species <- colnames(lcmat)[lcmat[spc,]==1];
			cooccr_species <- cooccr_species[!cooccr_species %in% clique_species];
			cooccr_species_no <- match(cooccr_species,colnames(lcmat));
			cooccr_cliques <- unique(initial_heights$clique_no[cooccr_species_no]);
			cooccr_cliques <- cooccr_cliques[!cooccr_cliques<=this_clique_no];
			# eliminate co-occuring species beneath this clique
			cooccr_species <- cooccr_species[initial_heights$clique_no[cooccr_species_no] %in% cooccr_cliques];
			cooccr_species_no <- cooccr_species_no[initial_heights$clique_no[cooccr_species_no] %in% cooccr_cliques]
			#test1 <- match(cooccr_species,rownames(initial_heights));
#			while (ht < max_height && order_kept==0 && length(cooccr_cliques)>0)	{
			while (clique_comp < length(cooccr_cliques))	{
				clique_comp <- clique_comp+1;
				next_clique_no <- cooccr_cliques[clique_comp];
				next_clique <- (1:cspec)[initial_heights$clique_no %in% next_clique_no];
				next_clique_spc <- rownames(initial_heights)[next_clique];
				next_clique_size <- length(next_clique_spc);
				max_height_ncl <- cspec; min_height_ncl <- 0; cur_height_ncl <- min(next_clique);
				for (cls in 1:length(next_clique))	{
					# note: we are using all of the species names because fa_above & fa_below were not culled
					species_above <- species_names[fa_above[match(next_clique_spc[cls],rownames(fa_above)),]>0];
					species_above_no <- match(species_above,rownames(initial_heights));
					species_above <- species_above[!is.na(species_above_no)];
					if (length(species_above)>0)
						max_height_ncl <- min(max_height_ncl,min(match(species_above,rownames(initial_heights))));
					species_below <- species_names[fa_below[match(next_clique_spc[cls],rownames(fa_below)),]>0];
					species_below_no <- match(species_below,rownames(initial_heights));
					species_below <- species_below[!is.na(species_below_no)];
					if (length(species_below)>0)
						min_height_ncl <- max(min_height_ncl,max(match(species_below,rownames(initial_heights))));
					}

#				ht <- max(next_clique+1);
				# remove species in the clique from list of species with which this clique co-occurs;
				if (max_height_ncl < max_height_cl)	{
					# shouldn't be needed, but don't bother if this clique already is above the co-occurring clique
					flipflop1 <- flipflop2 <- 1:cspec;
					}	else if (max_height_cl %in%  min_height_ncl:max_height_ncl && (cur_height_ncl %in%  min_height_ncl:max_height_ncl && cur_height_ncl %in%  min_height_cl:max_height_cl) && max(clique+1)==min(next_clique))	{
					# case where possible ranges of pairs overlap completely for our purposes
					pre <- 1:(cur_height_cl-1);
					pre <- pre[!pre %in% c(0,clique,next_clique)];
					post <- (1:cspec)[!(1:cspec) %in% c(pre,next_clique,clique)];
					flipflop1 <- flipflop2 <- c(pre,next_clique,clique,post);
					lost_spc <- (1:cspec)[!(1:cspec) %in% flipflop1];
#					new_gaps <- c(sum(tally_range_gaps(lcmat[flipflop1,flipflop1])),sum(tally_range_gaps(lcmat[flipflop2,flipflop2])));
					} else if (max_height_cl %in%  min_height_ncl:max_height_ncl && (cur_height_cl %in%  min_height_ncl:max_height_ncl && cur_height_ncl %in%  min_height_cl:max_height_cl))	{
					# here, we can move clique up or next_clique down
					pre1 <- 1:(cur_height_cl-1);
					pre1 <- pre1[!pre1 %in% c(0,clique,next_clique)];
					pre2 <- 1:(cur_height_ncl-1);
					pre2 <- pre2[!pre2 %in% c(0,clique,next_clique)];
					post1 <- (1:cspec)[!(1:cspec) %in% c(pre1,next_clique,clique)];
					post2 <- (1:cspec)[!(1:cspec) %in% c(pre2,next_clique,clique)];
					flipflop1a <- c(pre1,clique,next_clique,post1);
					flipflop1b <- c(pre1,next_clique,clique,post1);
					flipflop2a <- c(pre2,clique,next_clique,post2);
					flipflop2b <- c(pre2,next_clique,clique,post2);
					new_gaps <- c(sum(tally_range_gaps(lcmat[flipflop1a,flipflop1a])),sum(tally_range_gaps(lcmat[flipflop1b,flipflop1b])));
					if (match(min(new_gaps),new_gaps)==1)	{
						flipflop1 <- flipflop1a;
						} else	{
						flipflop1 <- flipflop1b;
						}
					new_gaps <- c(sum(tally_range_gaps(lcmat[flipflop2a,flipflop2a])),sum(tally_range_gaps(lcmat[flipflop2b,flipflop2b])));
					if (match(min(new_gaps),new_gaps)==1)	{
						flipflop2 <- flipflop2a;
						} else	{
						flipflop2 <- flipflop2b;
						}
					lost_spc <- (1:cspec)[!(1:cspec) %in% flipflop1];
					} else if (min_height_ncl < max_height_cl && min_height_ncl>cur_height_cl)	{
					# case where lower species cannot rise as high as higher species, but their possible ranges overlap
					pre <- 1:(min_height_ncl-1);
					pre <- pre[!pre %in% c(next_clique,clique)];
					post <- (1:cspec)[!(1:cspec) %in% c(pre,next_clique,clique)];
					flipflop1 <- c(pre,clique,next_clique,post);
					flipflop2 <- c(pre,next_clique,clique,post);
					lost_spc <- (1:cspec)[!(1:cspec) %in% flipflop1];
					#lost_species <- species_names[!species_names %in% rownames(lcmat)[flipflop1]];
					} else if (min_height_ncl < max_height_cl)	{
					pre <- 1:min(c(clique,next_clique)-1);
					pre <- pre[!pre %in% c(next_clique,clique)];
					post <- (1:cspec)[!(1:cspec) %in% c(pre,next_clique,clique)];
					flipflop1 <- c(pre,next_clique,clique,post);
					flipflop2 <- c(pre,clique,next_clique,post);
					lost_spc <- (1:cspec)[!(1:cspec) %in% flipflop1];
					lost_species <- species_names[!species_names %in% rownames(lcmat)[flipflop1]];
					} else if (min_height_ncl > max_height_cl)	{
					flipflop1 <- flipflop2 <- 1:cspec;
					} else	{
					whatsupdoc <- rbind(whatsupdoc,c(clique_species[1],min_height_cl,max_height_cl,next_clique_spc[1],min_height_ncl,max_height_ncl))
					}
				if (length(lost_spc)==0)	{
					#lost_species2 <- species_names[!species_names %in% rownames(lcmat)[flipflop1]];
					new_gaps <- c(sum(tally_range_gaps(lcmat[flipflop1,flipflop1])),sum(tally_range_gaps(lcmat[flipflop2,flipflop2])));
					if (sum(new_gaps<=best_gaps)>0 && sum(flipflop1!=(1:cspec))>0)	{
						if (match(min(new_gaps),new_gaps)==1){
							lcmat <- cooccur_matrix_better <- lcmat[flipflop1,flipflop1];
							} else	{
							lcmat <- cooccur_matrix_better <- lcmat[flipflop2,flipflop2];
							}
						lost_species3 <- species_names[!species_names %in% rownames(lcmat)]
						best_gaps <- min(new_gaps);
						clique <- match(clique_species,colnames(cooccur_matrix_better));
						initial_heights <- initial_heights[match(rownames(cooccur_matrix_better),rownames(initial_heights)),];
						good_cliques <- good_cliques[match(good_cliques$clique_no,unique(initial_heights$clique_no)),];
						cur_height_cl <- match(clique_species[1],rownames(cooccur_matrix_better));
						} else	{
#						cmat <- lcmat;
						}
					} else	{
					problems <- rbind(problems,c(cl,next_clique_no));
					}
	#			print(clique)
	#			i <- i+1;
				}
			}	else	{
#			(1:cspec)[cooccr_mat_adj[spc,]==1]
			while (ht < max_height_ncl && order_kept==min_height_ncl && (attempt %% 1)==0)	{
	#		init_hts <- match(rownames(initial_heights)[other_coocc_sp],rownames(initial_heights));
				clique_comp <- clique_comp+1;
				this_clique_no <- this_clique_no+1;
				next_clique <- (1:cspec)[initial_heights$clique_no %in% this_clique_no];
				next_clique_spc <- rownames(initial_heights)[next_clique];
				# remove species in the clique from list of species with which this clique co-occurs;
				pre <- 1:min(c(clique,next_clique)-1);
				pre <- pre[!pre %in% c(next_clique,clique)];
				post <- (1:cspec)[!(1:cspec) %in% c(pre,next_clique,clique)];
	#			post <- post[post %in% 1:cspec];
	#		post <- post[!post %in% c(next_clique,clique)];
			flipflop <- c(pre,next_clique,clique,post);
			lost_spc <- (1:cspec)[!(1:cspec) %in% flipflop];
			lost_species <- species_names[!species_names %in% rownames(lcmat)[flipflop]];
			if (length(lost_spc)==0)	{
	#		(1:cspec)[!(1:cspec) %in% flipflop]
	#		length(flipflop)
				next_clique_size <- length(next_clique);
				lost_species2 <- species_names[!species_names %in% rownames(lcmat)[flipflop]];
				cmat <- lcmat[flipflop,flipflop];
				new_gaps <- sum(tally_range_gaps(cmat));
				if (new_gaps<=best_gaps)	{
					lcmat <- cooccur_matrix_better <- cmat;
					lost_species3 <- species_names[!species_names %in% rownames(cmat)]
					best_gaps <- new_gaps;
	#				clique <- adjust;
					} else	{
					cmat <- lcmat;
					}
				} else	{
				problems <- rbind(problems,c(cl,this_clique_no));
				}
			ht <- max(next_clique+1);
	#		print(clique)
	#		i <- i+1;
			}
			}
		if (best_gaps<old_best)	{
			cooccur_matrix_best <- cooccur_matrix_better;
			new_order <- match(rownames(cooccur_matrix_best),rownames(initial_heights));
#			order_kept <- sum(sort(initial_heights$clique_no)!=initial_heights$clique_no)
			lost <- (1:cspec)[!(1:cspec) %in% new_order];
#			initial_heights <- initial_heights[new_order,];
			initial_heights$clique_no <- match(initial_heights$clique_no,unique(initial_heights$clique_no));
			good_cliques$clique_no <- unique(initial_heights$clique_no);
			cl <- cl-1;
			old_best <- best_gaps;
	#		best_gaps <- best_gaps+1;
			}
		}
	write.csv(initial_heights,paste("Try_",attempt,"_Sort_Heights.csv",sep=""));
	write.csv(lcmat,paste("Try_",attempt,"_Cooccurrence_Matrix.csv",sep=""));
	improvement <- last_best_gaps-best_gaps;
	print(paste("This run improved the fit by",improvement,"gaps"));
	attempt <- attempt+1;
	}
return(cooccur_matrix_best);
}

# find age assignments that minimize gaps and range extensions among species
# revised 2020-02-19
# revised 2020-06-15
optimo_paleodb_collection_and_occurrence_stratigraphy <- function(paleodb_finds,paleodb_collections,hierarchical_chronostrat,zone_database,update_search=T)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- remove_duplicate_occurrences_paleodb(occurrences=paleodb_finds);

# make sure that key collections fields are characters, not factors
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fixed_names <- c();
for (fn in 1:length(to_fix))	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]))
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp."
paleodb_finds <- expello_indeterminate_species(paleodb_finds);
taxon_names <- sort(unique(paleodb_finds$accepted_name));

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;

# number taxa
paleodb_finds$taxon_no <- match(taxa,taxon_names);
#paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);

finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
finest_chronostrat <- unique(finest_chronostrat);

### problem is before here: for some reason, there are NAs in the ages.
taxon <- taxon_names;
print(paste("Getting basic stratigraphic range data for",length(taxon),"species..."));
m_r_d <- data.frame(base::t(pbapply::pbsapply(taxon,tally_fuzzy_stratigraphic_ranges_sapply,all_finds=paleodb_finds,hierarchical_chronostrat)));
minimum_range_data <- data.frame(ma_fa_mx=as.numeric(unlist(m_r_d$ma_fa_lb)),
								 ma_fa_mn=as.numeric(unlist(m_r_d$ma_fa_ub)),
								 ma_la_mx=as.numeric(unlist(m_r_d$ma_la_lb)),
								 ma_la_mn=as.numeric(unlist(m_r_d$ma_la_ub)),
								 interval_lb=as.character(unlist(m_r_d$interval_lb)),
								 interval_ub=as.character(unlist(m_r_d$interval_ub)),
								 stringsAsFactors = F);
rownames(minimum_range_data) <- taxon_names;

zone_taxa <- taxon_names[taxon_names %in% zone_database$zone];
z_n <- match(zone_taxa,rownames(minimum_range_data));
these_zones <- zone_database[zone_database$zone %in% zone_taxa,];
these_zones <- these_zones[order(these_zones$zone),];
these_zones_keep <- these_zones[1:length(zone_taxa),];
for (zt in 1:length(zone_taxa))	{
	these_zones_keep[zt,] <- these_zones[match(zone_taxa[zt],these_zones$zone),];
	if (sum(these_zones$zone %in% zone_taxa[zt])>1)	{
		temp_zones <- these_zones[these_zones$zone %in% zone_taxa[zt],];
		these_zones_keep$ma_lb[zt] <- max(temp_zones$ma_lb);
		these_zones_keep$ma_ub[zt] <- min(temp_zones$ma_ub);
		}
	}

zone_rows <- match(these_zones_keep$zone,rownames(minimum_range_data))
print(paste("Double checking minimum stratigraphic ranges of",length(zone_rows),"zone taxa"));
# first, put in ranges for as-yet unranged zone taxa;
blank_zone <- minimum_range_data$ma_fa_mn[zone_rows] %in% 0;
zone_rows_blank <- zone_rows[blank_zone];
age <- these_zones_keep$ma_lb[blank_zone];
minimum_range_data$interval_lb[zone_rows_blank] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
minimum_range_data$ma_fa_mx[zone_rows_blank] <- minimum_range_data$ma_fa_mn[zone_rows_blank] <- age;
age <- these_zones_keep$ma_ub[blank_zone];
minimum_range_data$interval_ub[zone_rows_blank] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
minimum_range_data$ma_la_mx[zone_rows_blank] <- minimum_range_data$ma_la_mn[zone_rows_blank] <- age;
make_older <- minimum_range_data$ma_fa_mx[zone_rows] < these_zones_keep$ma_lb
minimum_range_data$ma_fa_mx[zone_rows[make_older]] <- these_zones_keep$ma_lb[make_older];
make_older <- minimum_range_data$ma_fa_mn[zone_rows] < these_zones_keep$ma_lb;
minimum_range_data$ma_fa_mn[zone_rows[make_older]] <- these_zones_keep$ma_lb[make_older];
make_younger <- minimum_range_data$ma_la_mx[zone_rows] > these_zones_keep$ma_ub;
minimum_range_data$ma_la_mx[zone_rows[make_younger]] <- these_zones_keep$ma_ub[make_younger];
#cbind(minimum_range_data$ma_la_mx[zone_rows[make_younger]],these_zones_keep$ma_ub[make_younger])
make_younger <- minimum_range_data$ma_la_mn[zone_rows] > these_zones_keep$ma_ub;
minimum_range_data$ma_la_mn[zone_rows[make_younger]] <- these_zones_keep$ma_ub[make_younger];
age <- minimum_range_data$ma_fa_mx[zone_rows];
minimum_range_data$interval_lb[zone_rows] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- minimum_range_data$ma_la_mn[zone_rows];
minimum_range_data$interval_ub[zone_rows] <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
#cbind(minimum_range_data$ma_la_mn[zone_rows[make_younger]],these_zones_keep$ma_ub[make_younger])

# now, adjust ranges as needed
#zone_rows
#minimum_range_data$interval_lb[match(zone_taxa,rownames(minimum_range_data))]
#minimum_range_data$interval_ub[match(zone_taxa,rownames(minimum_range_data))]
#zt <- 0;
#while (zt < length(zone_taxa))	{
#for (zt in 1:length(zone_taxa))	{
#	zt <- zt+1;
#	zn <- match(zone_taxa[zt],taxon_names);
#	if (minimum_range_data$ma_fa_mn[zn]==0)
#		print(zt)
#	this_zone <- subset(zone_database,zone_database$zone==zone_taxa[zt]);

#	if (minimum_range_data$ma_fa_mn[zn]==0)	{
#		minimum_range_data$ma_fa_mx[zn] <- minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
#		minimum_range_data$ma_la_mx[zn] <- minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
#		minimum_range_data$interval_lb[zn] <- rebin_collection_with_time_scale(age=minimum_range_data$ma_fa_mx[zn],onset_or_end = "onset",fine_time_scale = finest_chronostrat)
#		minimum_range_data$interval_ub[zn] <- rebin_collection_with_time_scale(age=minimum_range_data$ma_la_mn[zn],onset_or_end = "end",fine_time_scale = finest_chronostrat)
##		minimum_range_data$interval_lb[zn] <- this_zone$interval_lb[match(max(this_zone$ma_lb),this_zone$ma_lb)];
#		minimum_range_data$interval_ub[zn] <- this_zone$interval_ub[match(min(this_zone$ma_ub),this_zone$ma_ub)];
#		} else	{
#		if (minimum_range_data$ma_fa_mx[zn]<max(this_zone$ma_lb))	{
#			minimum_range_data$ma_fa_mx[zn] <- max(this_zone$ma_lb);
#			minimum_range_data$interval_lb[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_lb)<=finest_chronostrat$ma_lb))];
#			}
#		if (minimum_range_data$ma_fa_mn[zn]<max(this_zone$ma_lb)){
#			minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
#			}
#		if (minimum_range_data$ma_la_mx[zn]>min(this_zone$ma_ub)){
#			minimum_range_data$ma_la_mx[zn] <- min(this_zone$ma_ub);
#			minimum_range_data$interval_ub[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_ub)<=finest_chronostrat$ma_ub))];
#			}
#		if (minimum_range_data$ma_la_mn[zn]>min(this_zone$ma_ub))	{
#			minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
#			}
#		}
#	}
#ttt <- (1:ntaxa)[minimum_range_data$ma_fa_mx<minimum_range_data$ma_la_mx];
#which(is.na(minimum_range_data),arr.ind = T)
bin_lb <- hierarchical_chronostrat$bin_first[match(unique(paleodb_collections$interval_lb),hierarchical_chronostrat$interval)];
bin_ub <- hierarchical_chronostrat$bin_last[match(unique(paleodb_collections$interval_ub),hierarchical_chronostrat$interval)];
bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)];
# kluge: eliminate the need for this;
xx <- (1:ncolls)[bin_lb>bin_ub];
if (length(xx)>0)	{
	dummy <- bin_lb[xx];
	bin_lb[xx] <- bin_ub[xx];
	bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[bin_ub[xx]];
	}

#unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
bin_fuzz <- bin_ub-bin_lb;
#hist(bin_fuzz,breaks=-1:max(bin_fuzz))
unfixed_collections <- paleodb_collections$collection_no[bin_fuzz>0];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
min_gap_to_set <- 0;
min_taxa_to_set <- 3;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(bin_ub-bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search)
	print(progress);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- 0;
	uc <- 0;	# uc <- match(problems[12],unfixed_collections)
 	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges <- subset(index_ranges_all,index_ranges_all$ma_fa_mn!=0);
		if (nrow(index_ranges)>0)	{
			ranges <- cbind(hierarchical_chronostrat$bin_first[match(index_ranges$interval_lb,hierarchical_chronostrat$interval)],hierarchical_chronostrat$bin_last[match(index_ranges$interval_ub,hierarchical_chronostrat$interval)]);
			# something is resetting bin_lb & bin_ub to NA; might be ma_lb &/or ma_ub doing it
			poss_range <- bin_lb[coll_no]:bin_ub[coll_no];

			bin_gaps <- array(0,dim=c(length(poss_range)));
			for (rr in 1:nrow(ranges))
				bin_gaps <- bin_gaps+procrustean_binning_one_taxon(range=ranges[rr,],poss_range=c(min(poss_range),max(poss_range)),debug=F);

			# case where some of the previously assigned range induces no gaps in known species
			if ((min(bin_gaps)<=min_gap_to_set && length(poss_range)>sum(bin_gaps<=min_gap_to_set)) && nrow(index_ranges)>=min_taxa_to_set)	{
				bin_lb[coll_no] <- min(poss_range[bin_gaps==min(bin_gaps)]);
				bin_ub[coll_no] <- max(poss_range[bin_gaps==min(bin_gaps)]);
				new_and_improved <- new_and_improved+1;
				paleodb_collections$interval_lb[coll_no] <- finest_chronostrat$interval[match(bin_lb[coll_no],finest_chronostrat$bin_first)];
				paleodb_collections$interval_ub[coll_no] <- finest_chronostrat$interval[match(bin_ub[coll_no],finest_chronostrat$bin_last)];
				paleodb_collections$ma_lb[coll_no] <- min(paleodb_collections$ma_lb[coll_no],finest_chronostrat$ma_lb[match(bin_lb[coll_no],finest_chronostrat$bin_first)]);
				paleodb_collections$ma_ub[coll_no] <- max(paleodb_collections$ma_ub[coll_no],finest_chronostrat$ma_ub[match(bin_ub[coll_no],finest_chronostrat$bin_last)]);
				# case where we've narrowed it down to one bin
				if (sum(bin_gaps==min_gap_to_set)==1)	{
					# the earliest possible first occurrences must be at least as old as the oldest possible first occurrences
					#	for those taxa first known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_fa_mx[sub_index_ranges$ma_fa_mx<max(sub_index_ranges$ma_fa_mn)] <- max(sub_index_ranges$ma_fa_mn);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}

					# the latest possible last occurrences must be at least as young as the youngest possible last occurrences
					#	for those taxa last known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_la_mn[sub_index_ranges$ma_la_mn>min(sub_index_ranges$ma_la_mx)] <- min(sub_index_ranges$ma_la_mx);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}
					indexed_species <- rownames(index_ranges);
					i_s <- match(indexed_species,rownames(minimum_range_data));
					minimum_range_data[i_s,] <- index_ranges;

					if (nrow(index_ranges_all) > nrow(index_ranges))	{
					# update species that do not have constrained finds yet
						unindexed_species <- rownames(index_ranges_all)[index_ranges_all$ma_fa_mn==0];
						u_s <- match(unindexed_species,rownames(minimum_range_data));
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_fa_mx[u_s] <- min(sub_index_ranges$ma_fa_mx);
							minimum_range_data$ma_fa_mn[u_s] <- max(sub_index_ranges$ma_fa_mn);
							} else	{
							minimum_range_data$ma_fa_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_fa_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_la_mx[u_s] <- min(sub_index_ranges$ma_la_mx);
							minimum_range_data$ma_la_mn[u_s] <- max(sub_index_ranges$ma_la_mn);
							} else	{
							minimum_range_data$ma_la_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_la_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						minimum_range_data$interval_lb[u_s] <- paleodb_collections$interval_lb[coll_no];
						minimum_range_data$interval_ub[u_s] <- paleodb_collections$interval_ub[coll_no];
						}
					newly_fixed <- newly_fixed+1;
					}
				}
 			}
 		if (sum(bin_lb==Inf)>0 || sum(is.na(bin_lb))>0)
 			print(uc);
#		print(cbind(uc,minimum_range_data[2309,]));
#		print(c(paleodb_collections$ma_lb[195],paleodb_collections$ma_ub[195]));
		}
#	(1:ncolls)[is.na(as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub))];
	unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
	unfixed <- length(unfixed_collections);
	improved <- new_and_improved;
	fixed <- newly_fixed;
	unfixed <- length(unfixed_collections);
	if (improved==0 && min_taxa_to_set>0)	{
		reboots <- reboots+1;
		if ((reboots %% 2)==0)	{
			min_taxa_to_set <- min_taxa_to_set-1;
			} else	{
			min_gap_to_set <- min_gap_to_set+1;
			}
		if (min_taxa_to_set>0)
			improved <-max(1,improved);
		}
	attempt <- attempt+1;
	fcolls <- hist(bin_ub-bin_lb,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
#	fcolls <- hist(bin_ub-bin_lb,breaks=c(fuzz,max(fuzz)+1),plot=F)$counts;
	progress <- cbind(progress,fcolls);
	if (update_search)
		print(progress);
	}

paleodb_collections$bin_lb <- bin_lb;
paleodb_collections$bin_ub <- bin_ub;
return(paleodb_collections);
}

# written 2020-06-05: makes decisions based on probabilities of filling gaps.
# paleodb_finds <- pbdb_finds; paleodb_collections <- refined_collections_prov;
# paleodb_finds$ma_lb <- paleodb_collections$ma_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]; paleodb_finds$ma_ub <- paleodb_collections$ma_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
# paleodb_finds$interval_lb <- paleodb_collections$interval_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]; paleodb_finds$interval_ub <- paleodb_collections$interval_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
ml_paleodb_collection_and_occurrence_stratigraphy <- function(paleodb_finds,paleodb_collections,hierarchical_chronostrat,zone_database,update_search=T,update_type="Graph",max_rho=0.1,precision=0.1)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- remove_duplicate_occurrences_paleodb(occurrences=paleodb_finds);

paleodb_finds <- paleodb_finds[paleodb_finds$collection_no %in% paleodb_collections$collection_no,];
# make sure that key collections fields are characters, not factors
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fixed_names <- c();
for (fn in 1:length(to_fix))
	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]))
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp."
paleodb_finds <- expello_indeterminate_species(paleodb_finds);
taxon_names <- sort(unique(paleodb_finds$accepted_name));

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;

# number taxa
paleodb_finds$taxon_no <- match(taxa,taxon_names);
#paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
if (is.na(match("ma_lb",colnames(paleodb_finds))))	{
	paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
	paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
	}
if (is.na(match("interval_lb",colnames(paleodb_finds))))	{
	paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
	paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);
	}

finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
finest_chronostrat <- unique(finest_chronostrat);

poss_ages <- sort(unique(c(paleodb_collections$ma_lb,paleodb_collections$ma_ub)),decreasing=T);

if (is.na(match("bin_lb",colnames(paleodb_collections))))	{
#	paleodb_collections$bin_lb <- hierarchical_chronostrat$bin_first[match(unique(paleodb_collections$interval_lb),hierarchical_chronostrat$interval)];
#	paleodb_collections$bin_ub <- hierarchical_chronostrat$bin_last[match(unique(paleodb_collections$interval_ub),hierarchical_chronostrat$interval)];
	paleodb_collections$bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
	paleodb_collections$bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)];
	}
# kluge: eliminate the need for this;
xx <- (1:ncolls)[paleodb_collections$bin_lb>paleodb_collections$bin_ub];
if (length(xx)>0)	{
	dummy <- paleodb_collections$bin_lb[xx];
	paleodb_collections$bin_lb[xx] <- paleodb_collections$bin_ub[xx];
	paleodb_collections$bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[paleodb_collections$bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[paleodb_collections$bin_ub[xx]];
	paleodb_collections$ma_lb[xx] <- finest_chronostrat$ma_lb[paleodb_collections$bin_lb[xx]];
	paleodb_collections$ma_ub[xx] <- finest_chronostrat$ma_ub[paleodb_collections$bin_ub[xx]];
	}

fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
bin_sites <- vector(length=max(paleodb_collections$bin_ub));
bin_sites[1:max(fixed_collections$bin_lb)] <- hist(fixed_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
names(bin_sites)[1:max(paleodb_collections$bin_ub)] <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
bin_sites_cum <- cumsum(bin_sites);

#minimum_range_data <- data.frame(ma_fa_mx=as.numeric(rep(0,ntaxa)),ma_fa_mn=as.numeric(rep(0,ntaxa)),ma_la_mx=as.numeric(rep(0,ntaxa)),ma_la_mn=as.numeric(rep(0,ntaxa)),interval_lb=as.character(rep("",ntaxa)),interval_ub=as.character(rep("",ntaxa)),N=as.numeric(rep(0,ntaxa)),Na=as.numeric(rep(0,ntaxa)),Nz=as.numeric(rep(0,ntaxa)),fa_imp=as.numeric(rep(0,ntaxa)),la_imp=as.numeric(rep(0,ntaxa)),site_lb=as.numeric(rep(0,ntaxa)),site_ub=as.numeric(rep(0,ntaxa)),range=as.numeric(rep(0,ntaxa)),rho=as.numeric(rep(0,ntaxa)),stringsAsFactors = F);
#rownames(minimum_range_data) <- taxon_names;
# this is putting everything in one bin too late: fix it
minimum_range_data <- accersi_minimum_range_information(taxon_names,paleodb_finds,paleodb_collections,bin_sites,finest_chronostrat,max_rho=0.1,precision=precision);
#minimum_range_data_test <- subset(minimum_range_data,minimum_range_data$interval_lb!="");

bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search & tolower(update_type)=="graph")	{
	all_colors <- rainbow(7);
	mxx <- max(progress[,1])+1;
	mxy <- log10(ceiling(nrow(paleodb_collections)/10^floor(log10(nrow(paleodb_collections))))*10^floor(log10(nrow(paleodb_collections))));
	poss_bins <- 0:mxx;
	par(mfrow=c(1,1))
	par(pin=c(4.5,3));
	main_title <- "Refining Chronostratigraphy of Uncertain Sites";
	plot(NA,type='n',axes=FALSE,main=main_title,xlab="Bin Precision",ylab="Collections",xlim=c(0,mxx),ylim=c(0,mxy));
	if (mxx<25)	{
		axis_labels <- 1:mxx
		} else	{
		axis_labels <- seq(2,mxx+1,by=2)-1;
		}
	specified_axis_w_labels(axe=1,max_val=mxx,min_val=0,maj_break=1,med_break=1,min_break=1,label_pos="mid",axis_labels=axis_labels,font_size=0.75);
	numbers <- base_numbers <- c(1,5);
	nn <- 0;
	while (nn < floor(mxy))	{
		nn <- nn+1;
		numbers <- c(numbers,base_numbers*10^nn)
		}
	numbers <- numbers[log10(numbers)<=mxy];
	log10_axes(axe=2,min_ax=0,max_ax=mxy,numbers,font_size=5/6,orient=2);
	for (pp in 1:mxx)	{
		if (progress[pp,2]>0)
			rect(pp-1,log10(0.7),pp,log10(progress[pp,2]),col=all_colors[1],border="gray50",lwd=0.5);
		}
	} else if (update_search)	{
	print(progress);
	}

#collection_nos <- sort(paleodb_collections$collection_no);
list_lengths <- vector(length=length(unfixed_collections));
for (nc in 1:length(unfixed_collections))
	list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];

finest_chronostrat$span <- abs(finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- uc <- 0;	# uc <- match(problems[12],unfixed_collections)
	extensions <- modified_species <- c();
	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges <- subset(index_ranges_all,index_ranges_all$rho!=0);

		if (paleodb_collections$ma_lb[coll_no]==paleodb_collections$ma_ub[coll_no])	{
			paleodb_collections$ma_lb[coll_no] <- finest_chronostrat$ma_lb[paleodb_collections$bin_lb[coll_no]];
			paleodb_collections$ma_ub[coll_no] <- finest_chronostrat$ma_ub[paleodb_collections$bin_ub[coll_no]];
			}
		if (nrow(index_ranges)>1)	{
			site_info <- paleodb_collections[coll_no,colnames(paleodb_collections) %in% c("ma_lb","ma_ub","bin_lb","bin_ub","interval_lb","interval_ub")];
			bl <- site_info$bin_lb;
			prop_lower_bin <- (finest_chronostrat$ma_lb[bl]-site_info$ma_lb)/finest_chronostrat$span[bl];
			if (bl>0)	{
				site_info$site_lb <- bin_sites_cum[bl-1]+(prop_lower_bin*bin_sites[bl]);
				} else	{
				site_info$site_lb <- prop_lower_bin*bin_sites[bl];
				}
			bu <- site_info$bin_ub;
			prop_upper_bin <- abs(finest_chronostrat$ma_ub[bu]-site_info$ma_ub)/finest_chronostrat$span[bu];
			site_info$site_ub <- bin_sites_cum[bu] - (prop_upper_bin*bin_sites[bu]);
			# divide time span into smallest chunks #
			site_info$site_span <- 1+site_info$site_ub-site_info$site_lb;

			# take into account that ma_la_max might be greater than la_imp
#			if (sum(index_ranges$ma_fa_mx<site_info$ma_ub)>0)	extensions <- rbind(extensions,c(uc,"too young"));
#			if (sum(index_ranges$ma_la_mn>site_info$ma_lb)>0)	extensions <- rbind(extensions,c(uc,"too old"));

			# get the basics of this collection #
			poss_ages <- seq(site_info$ma_lb,site_info$ma_ub,by=-precision);
			# the possible finds will span 2+ bins; so, make sure that we get it correct.
			poss_finds <- relv_time <- c();
			for (si in site_info$bin_lb:site_info$bin_ub)	{
				if (si==site_info$bin_lb)	{
					relv_time <- c(relv_time,seq(round(site_info$ma_lb,2),round(finest_chronostrat$ma_ub[si],2),by=-precision));
					slc <- length(relv_time)-1;
					bn_p_slc <- (site_info$site_lb-bin_sites_cum[si-1])/slc;
	#				bn_p_slc <- (bin_sites[si]*(site_info$ma_lb-finest_chronostrat$ma_ub[si])/finest_chronostrat$span[si])/slc;
					poss_finds <- c(poss_finds,rep(bn_p_slc,slc));
					} else if (si==site_info$bin_ub)	{
					if (round(finest_chronostrat$ma_lb[si]-precision,2)==round(site_info$ma_ub,2))	{
						this_time <- round(finest_chronostrat$ma_lb[si]-precision,2);
						} else	{
						this_time <- seq((round(finest_chronostrat$ma_lb[si],2)-precision),round(site_info$ma_ub,2),by=-round(precision,2));
						}
					relv_time <- c(relv_time,this_time);
					slc <- length(this_time);
					bn_p_slc <- (site_info$site_ub-bin_sites_cum[si-1])/slc;
					#bn_p_slc <- (bin_sites[si]*abs(site_info$ma_ub-finest_chronostrat$ma_lb[si])/finest_chronostrat$span[si])/slc;
					poss_finds <- c(poss_finds,rep(bn_p_slc,slc))
					} else	{
					this_time <- seq(finest_chronostrat$ma_lb[si]-precision,finest_chronostrat$ma_ub[si],by=-precision);
					relv_time <- c(relv_time,this_time);
					slc <- length(this_time);
					bn_p_slc <- bin_sites[si]/slc;
					poss_finds <- c(poss_finds,rep(bn_p_slc,slc))
					}
				}

			# possible cases:
			lnl_ages <- array(0,dim=c(nrow(index_ranges),length(poss_ages)-1));
			rownames(lnl_ages) <- rownames(index_ranges);
			names(poss_finds) <- colnames(lnl_ages) <- (poss_ages[1:(length(poss_ages)-1)]+poss_ages[2:length(poss_ages)])/2;

			ir <- 0;
			for (ir in 1:nrow(index_ranges))	{
#			while (ir < nrow(index_ranges))	{
#				ir <- ir+1;
				if (index_ranges$fa_imp[ir]>=site_info$ma_lb & index_ranges$la_imp[ir]<=site_info$ma_ub)	{
					# entirely within range of taxon
	#				lnl_ages[ir,] <- rbind(lnl_ages,rep(0,ncol(lnl_ages)));
					} else if (index_ranges$fa_imp[ir]<=site_info$ma_lb & index_ranges$la_imp[ir] <= site_info$ma_ub) {
					# possibly older than the taxon is known to be
					# use numbers of finds per slice (taking into account different bins) and work up to FA
					relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
					# introduce something here to estimate the gap at lowering rhos with extended earlier bounds
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					} else if (index_ranges$fa_imp[ir]>=site_info$ma_lb & index_ranges$la_imp[ir]<site_info$ma_lb)	{
					# possibly younger than the taxon is known to be
					# use numbers of finds per slice (taking into account different bins) and work up to FA
					relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					} else if (index_ranges$fa_imp[ir]>=site_info$ma_ub & site_info$ma_ub>index_ranges$la_imp[ir])	{
					# possibly older or younger than the taxon is known to be
					relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					# do younger finds
					relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					} else if (index_ranges$ma_fa_lb[ir]<site_info$ma_ub)	{
					# known finds all are older than this collection
					relv_finds_cum <- cumsum(poss_finds)
					lnl_ages[ir,] <- lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
#					print(ir);
					} else if (index_ranges$ma_la_ub[ir]>site_info$ma_lb)	{
					relv_finds_cum <- cumsum(poss_finds[length(poss_finds):1]);
					lnl_ages[ir,] <- relv_finds_cum[length(poss_finds):1]*log(1-index_ranges$rho[ir]);
					# known finds all are younger than this collection
					} else if (index_ranges$fa_imp[ir]<site_info$ma_lb & index_ranges$la_imp[ir]>site_info$ma_ub)	{
					# possibly older finds
					relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					# possibly younger finds
					relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					} else if (index_ranges$la_imp[ir]>site_info$ma_ub)	{
					# possibly younger finds
					relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					} else if (index_ranges$fa_imp[ir]<site_info$ma_lb)	{
					# possibly older finds
					relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
					relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
					lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
					lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
					}
#				print(lnl_ages[ir,])
				}
			age_ln_likelihoods <- colSums(lnl_ages)-max(colSums(lnl_ages));
			new_age_range <- as.numeric(names(age_ln_likelihoods)[age_ln_likelihoods>=-2]);
			if (length(new_age_range)<length(age_ln_likelihoods))	{
				new_and_improved <- new_and_improved+1;
				paleodb_collections$ma_lb[coll_no] <- ceiling(max(new_age_range)/precision)*precision;
				paleodb_collections$ma_ub[coll_no] <- floor(min(new_age_range)/precision)*precision;
				paleodb_collections$interval_lb[coll_no] <- rebin_collection_with_time_scale(max(new_age_range),"onset",finest_chronostrat);
				paleodb_collections$interval_ub[coll_no] <- rebin_collection_with_time_scale(min(new_age_range),"end",finest_chronostrat);
				paleodb_collections$bin_lb[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb[coll_no],finest_chronostrat$interval)];
				paleodb_collections$bin_ub[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_ub[coll_no],finest_chronostrat$interval)];
				bin_fuzz[coll_no] <- paleodb_collections$bin_ub[coll_no]-paleodb_collections$bin_lb[coll_no];
				if (bin_fuzz[coll_no]==0)
					fixed <- fixed+1;
				paleodb_finds$ma_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_lb[coll_no];
				paleodb_finds$ma_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_ub[coll_no];
				paleodb_finds$interval_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_lb[coll_no];
				paleodb_finds$interval_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_ub[coll_no];
				modified_species <- unique(c(modified_species,index_species));
				} # finish updating info
			} # end case of having set species;
 		} # end search of unconstrained_collections

	fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
	bin_sites_old <- bin_sites;
	bin_sites_cum_old <- bin_sites_cum;
	bin_sites <- vector(length=max(paleodb_collections$bin_ub));
	bin_sites[1:max(fixed_collections$bin_lb)] <- hist(fixed_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
	names(bin_sites)[1:max(paleodb_collections$bin_ub)] <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
	bin_sites_cum <- cumsum(bin_sites);

	minimum_range_data <- accersi_minimum_range_information(taxon_names,paleodb_finds,paleodb_collections,bin_sites,finest_chronostrat,max_rho,precision);

	fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
 	bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
	fcolls <- hist(bin_fuzz,breaks=min(progress[,1]):(max(progress[,1])+1),plot=F)$counts;
#
	progress <- cbind(progress,fcolls);
	if (update_search && tolower(update_type)=="graph")	{
		iteration <- ncol(progress);
		if (iteration>(length(all_colors)+1))
			all_colors <- rainbow(iteration);
		interation_cols <- all_colors[1:(iteration-1)];
#		points(3*(1:length(interation_cols)),rep(4,length(interation_cols)),pch=21,bg=interation_cols)
#		interation_cols <- paste("gray",,sep="");
		interation_lwds <- c(0.5,rep(0,iteration-2));
		for (pp in 1:nrow(progress))	{
			print_order <- 1+order(progress[pp,2:iteration],2:iteration,decreasing = T);
			for (ppp in 1:length(print_order))	{
				if (progress[pp,print_order[ppp]]>0)	{
					rect(pp-1,log10(0.7),pp,log10(progress[pp,print_order[ppp]]),col=interation_cols[print_order[ppp]-1],border="gray50",lwd=interation_lwds[ppp]);
#					} else if (max(progress[pp,2:iteration]>0))	{
#					rect(pp-1,log10(0.7),pp,log10(0.7)/2,col="white",border="white",lwd=1.5);
					}
				}
			}
		} else if (update_search) {
		progress <- cbind(progress[1:length(fcolls),],fcolls);
		print(progress);
		}

	unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];
	improved <- unfixed-length(unfixed_collections);
	unfixed <- length(unfixed_collections);
	improved <- sum(abs(progress[,ncol(progress)]-progress[,ncol(progress)-1]));
	list_lengths <- vector(length=unfixed);
	for (nc in 1:length(unfixed_collections))
		list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
	unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];
 	}
paleodb_collections <- put_pbdb_dataframes_into_proper_type(pbdb_data = paleodb_collections);
paleodb_finds <- put_pbdb_dataframes_into_proper_type(pbdb_data = paleodb_finds);
output <- list(paleodb_collections,paleodb_finds);
names(output) <- c("ml_optimized_collections","ml_optimized_paleodb_finds");
return(output);
}

#paleodb_finds=pbdb_finds;paleodb_collections=refined_collections_prov
ml_paleodb_collection_and_occurrence_stratigraphy_taxon_partitioned <- function(paleodb_finds,paleodb_collections,sampling_groups,hierarchical_chronostrat,update_search=T,update_type="Graph",max_rho=0.1,precision=0.1)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- remove_duplicate_occurrences_paleodb(occurrences=paleodb_finds);

# make sure that key collections fields are characters, not factors
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fn <- 0;
fixed_names <- c();
while (fn < length(to_fix))	{
	fn <- fn+1;
	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]));
	}
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp."
paleodb_finds <- expello_indeterminate_species(paleodb_finds);

# bin taxa into sampling groups;
paleodb_finds$sampling_group <- sampling_groups$group[match(paleodb_finds$class,sampling_groups$classes)];
paleodb_finds$sampling_group[is.na(paleodb_finds$sampling_group)][paleodb_finds$genus[is.na(paleodb_finds$sampling_group)] %in% classic_inarticulate_brachiopod_taxa] <- "Inarticulata";
paleodb_finds <- subset(paleodb_finds,!is.na(paleodb_finds$sampling_group));
taxon_names <- sort(unique(paleodb_finds$accepted_name));
taxon_sampling_groups <- paleodb_finds$sampling_group[match(taxon_names,paleodb_finds$accepted_name)];
retained_colls <- sort(unique(paleodb_finds$collection_no));
#retained_colls[retained_colls %in% paleodb_collections$collection_no]
paleodb_collections <- paleodb_collections[paleodb_collections$collection_no %in% paleodb_finds$collection_no,];

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;
ncolls <- nrow(paleodb_collections);

# number taxa
paleodb_finds$taxon_no <- match(taxa,taxon_names);
#paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
if (is.na(match("ma_lb",colnames(paleodb_finds))))	{
	paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
	paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
	}
if (is.na(match("interval_lb",colnames(paleodb_finds))))	{
	paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
	paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);
	}

finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
finest_chronostrat <- unique(finest_chronostrat);

poss_ages <- sort(unique(c(paleodb_collections$ma_lb,paleodb_collections$ma_ub)),decreasing=T);

if (is.na(match("bin_lb",colnames(paleodb_collections))))	{
#	paleodb_collections$bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
#	paleodb_collections$bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval];
	paleodb_collections$bin_lb <- as.numeric(hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)]);
	paleodb_collections$bin_ub <- as.numeric(hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)]);
	}
#ugh <- (1:ncolls)[is.na(paleodb_collections$bin_ub)]
# kluge: eliminate the need for this;
xx <- (1:ncolls)[paleodb_collections$bin_lb>paleodb_collections$bin_ub];
xx <- xx[!is.na(xx)];
if (length(xx)>0)	{
	dummy <- paleodb_collections$bin_lb[xx];
	paleodb_collections$bin_lb[xx] <- paleodb_collections$bin_ub[xx];
	paleodb_collections$bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[paleodb_collections$bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[paleodb_collections$bin_ub[xx]];
	paleodb_collections$ma_lb[xx] <- finest_chronostrat$ma_lb[paleodb_collections$bin_lb[xx]];
	paleodb_collections$ma_ub[xx] <- finest_chronostrat$ma_ub[paleodb_collections$bin_ub[xx]];
	}

fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
control_groups <- sort(unique(sampling_groups$groups));
bin_sites_all <- array(0,dim=c(length(control_groups),length=max(paleodb_collections$bin_ub)));
rownames(bin_sites_all) <- control_groups;
colnames(bin_sites_all) <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
bin_sites_all_cum <- bin_sites_all;
for (cg in 1:length(control_groups))	{
	group_finds <- subset(paleodb_finds,paleodb_finds$sampling_group==control_groups[cg]);
	group_collections <- subset(fixed_collections,fixed_collections$collection_no %in% group_finds$collection_no);
	bin_sites_all[cg,1:max(fixed_collections$bin_lb)] <- hist(group_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
	bin_sites_all_cum[cg,] <- cumsum(bin_sites_all[cg,]);
	group_taxon_names <- taxon_names[taxon_sampling_groups==control_groups[cg]];
	if (cg==1)	{
		minimum_range_data <- accersi_minimum_range_information(taxon_names=group_taxon_names,
																paleodb_finds=group_finds,
																paleodb_collections=paleodb_collections,
																bin_sites=bin_sites_all[cg,],
																finest_chronostrat,max_rho=0.1,precision=precision);
		} else	{
		m_r_d <- accersi_minimum_range_information(taxon_names=group_taxon_names,
												 paleodb_finds=group_finds,
												 paleodb_collections=paleodb_collections,
												 bin_sites=bin_sites_all[cg,],
												 finest_chronostrat,max_rho=0.1,precision=precision);
		minimum_range_data <- rbind(minimum_range_data,m_r_d);
		}
	}
minimum_range_data <- minimum_range_data[order(rownames(minimum_range_data)),];

#minimum_range_data <- data.frame(ma_fa_mx=as.numeric(rep(0,ntaxa)),ma_fa_mn=as.numeric(rep(0,ntaxa)),ma_la_mx=as.numeric(rep(0,ntaxa)),ma_la_mn=as.numeric(rep(0,ntaxa)),interval_lb=as.character(rep("",ntaxa)),interval_ub=as.character(rep("",ntaxa)),N=as.numeric(rep(0,ntaxa)),Na=as.numeric(rep(0,ntaxa)),Nz=as.numeric(rep(0,ntaxa)),fa_imp=as.numeric(rep(0,ntaxa)),la_imp=as.numeric(rep(0,ntaxa)),site_lb=as.numeric(rep(0,ntaxa)),site_ub=as.numeric(rep(0,ntaxa)),range=as.numeric(rep(0,ntaxa)),rho=as.numeric(rep(0,ntaxa)),stringsAsFactors = F);
#rownames(minimum_range_data) <- taxon_names;
# this is putting everything in one bin too late: fix it
#minimum_range_data <- accersi_minimum_range_information(taxon_names,paleodb_finds,paleodb_collections,bin_sites,finest_chronostrat,max_rho=0.1,precision=precision);
#minimum_range_data_test <- subset(minimum_range_data,minimum_range_data$interval_lb!="");

bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search & tolower(update_type)=="graph")	{
	all_colors <- rainbow(7);
	mxx <- max(progress[,1])+1;
	mxy <- log10(ceiling(nrow(paleodb_collections)/10^floor(log10(nrow(paleodb_collections))))*10^floor(log10(nrow(paleodb_collections))));
	poss_bins <- 0:mxx;
	par(mfrow=c(1,1))
	par(pin=c(4.5,3));
	plot(NA,type='n',axes=FALSE,main="",xlab="Bin Precision",ylab="Collections",xlim=c(0,mxx),ylim=c(0,mxy));
	if (mxx<25)	{
		axis_labels <- 1:mxx
		} else	{
		axis_labels <- seq(2,mxx+1,by=2)-1;
		}
	specified_axis_w_labels(axe=1,max_val=mxx,min_val=0,maj_break=1,med_break=1,min_break=1,label_pos="mid",axis_labels=axis_labels,font_size=0.75);
	numbers <- base_numbers <- c(1,5);
	nn <- 0;
	while (nn < floor(mxy))	{
		nn <- nn+1;
		numbers <- c(numbers,base_numbers*10^nn)
		}
	numbers <- numbers[log10(numbers)<=mxy];
	log10_axes(axe=2,min_ax=0,max_ax=mxy,numbers,font_size=5/6,orient=2);
	for (pp in 1:mxx)	{
		if (progress[pp,2]>0)
			rect(pp-1,log10(0.7),pp,log10(progress[pp,2]),col=all_colors[1],border="gray50",lwd=0.5);
		}
	} else if (update_search)	{
	print(progress);
	}

#collection_nos <- sort(paleodb_collections$collection_no);
list_lengths <- vector(length=length(unfixed_collections));
for (nc in 1:length(unfixed_collections))
	list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];

finest_chronostrat$span <- abs(finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- uc <- 0;	# uc <- match(problems[12],unfixed_collections)
	modified_species <- c();
	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges_relv <- subset(index_ranges_all,index_ranges_all$rho!=0);

		if (nrow(index_ranges_relv)>1)	{
			relv_control_groups <- unique(taxon_sampling_groups[match(rownames(index_ranges_relv),taxon_names)]);
			r_c_g <- match(relv_control_groups,control_groups);
			rcg <- length(r_c_g);
			site_info <- s_inf <- paleodb_collections[coll_no,colnames(paleodb_collections) %in% c("ma_lb","ma_ub","bin_lb","bin_ub","interval_lb","interval_ub")];
			bl <- site_info$bin_lb;
			bu <- site_info$bin_ub;
			lnl_ages <- c();
			for (tcg in 1:rcg)	{
				cg <- r_c_g[tcg];	# cg is the control group number
				prop_lower_bin <- (finest_chronostrat$ma_lb[bl]-site_info$ma_lb)/finest_chronostrat$span[bl];		# prop of span from bottom to top
				prop_upper_bin <- abs(finest_chronostrat$ma_ub[bu]-site_info$ma_ub)/finest_chronostrat$span[bu];	# prop of span from top to bottom
				if (bl>0)	{
					site_info$site_lb <- floor(bin_sites_all_cum[cg,bl-1]+(prop_lower_bin*bin_sites_all[cg,bl]));
					} else	{
					site_info$site_lb <- floor(prop_lower_bin*bin_sites_all[cg,bl]);
					}
				site_info$site_ub <- ceiling(bin_sites_all_cum[cg,bu] - (prop_upper_bin*bin_sites_all[cg,bu]));
				bin_sites <- bin_sites_all[cg,];
				index_ranges <- index_ranges_relv[match(taxon_sampling_groups[match(rownames(index_ranges_relv),taxon_names)],control_groups)==cg,];

				lnl_ages <- rbind(lnl_ages,accersi_strat_likelihoods_of_possible_site_ages(site_info=site_info,index_ranges=index_ranges,bin_sites=bin_sites_all[cg,],precision=precision));
				# Now, do likelihoods for each group #
				}

			age_ln_likelihoods <- colSums(lnl_ages)-max(colSums(lnl_ages));
			new_age_range <- as.numeric(names(age_ln_likelihoods)[age_ln_likelihoods>=-2]);
			if (length(new_age_range)<length(age_ln_likelihoods))	{
				new_and_improved <- new_and_improved+1;
				paleodb_collections$ma_lb[coll_no] <- ceiling(max(new_age_range)/precision)*precision;
				paleodb_collections$ma_ub[coll_no] <- floor(min(new_age_range)/precision)*precision;
				paleodb_collections$interval_lb[coll_no] <- rebin_collection_with_time_scale(max(new_age_range),"onset",finest_chronostrat);
				paleodb_collections$interval_ub[coll_no] <- rebin_collection_with_time_scale(min(new_age_range),"end",finest_chronostrat);
				paleodb_collections$bin_lb[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb[coll_no],finest_chronostrat$interval)];
				paleodb_collections$bin_ub[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_ub[coll_no],finest_chronostrat$interval)];
				bin_fuzz[coll_no] <- paleodb_collections$bin_ub[coll_no]-paleodb_collections$bin_lb[coll_no];
				if (bin_fuzz[coll_no]==0)
					fixed <- fixed+1;
				paleodb_finds$ma_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_lb[coll_no];
				paleodb_finds$ma_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_ub[coll_no];
				paleodb_finds$interval_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_lb[coll_no];
				paleodb_finds$interval_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_ub[coll_no];
				modified_species <- unique(c(modified_species,index_species));
				} # finish updating info
			} # end case of having set species;
 		} # end search of unconstrained_collections

	fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
	bin_sites_old <- bin_sites_all;
	bin_sites_cum_old <- bin_sites_all_cum;
#	bin_sites_all <- array(0,dim=c(length(control_groups),length=max(paleodb_collections$bin_ub)));
#	rownames(bin_sites_all) <- control_groups;
#	colnames(bin_sites_all) <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
#	bin_sites_all_cum <- bin_sites_all;
	for (cg in 1:length(control_groups))	{
		group_finds <- subset(paleodb_finds,paleodb_finds$sampling_group==control_groups[cg]);
		group_collections <- subset(fixed_collections,fixed_collections$collection_no %in% group_finds$collection_no);
		bin_sites_all[cg,1:max(fixed_collections$bin_lb)] <- hist(group_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
		bin_sites_all_cum[cg,] <- cumsum(bin_sites_all[cg,]);
		group_taxon_names <- taxon_names[taxon_sampling_groups==control_groups[cg]];
		if (cg==1)	{
			minimum_range_data <- accersi_minimum_range_information(taxon_names=group_taxon_names,
																	paleodb_finds=group_finds,
																	paleodb_collections=paleodb_collections,
																	bin_sites=bin_sites_all[cg,],
																	finest_chronostrat,max_rho=0.1,precision=precision);
			} else	{
			m_r_d <- accersi_minimum_range_information(taxon_names=group_taxon_names,
													   paleodb_finds=group_finds,
													   paleodb_collections=paleodb_collections,
													   bin_sites=bin_sites_all[cg,],
													   finest_chronostrat,max_rho=0.1,precision=precision);
			minimum_range_data <- rbind(minimum_range_data,m_r_d);
			}
		}
	minimum_range_data <- minimum_range_data[order(rownames(minimum_range_data)),];

	fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
 	bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
	fcolls <- hist(bin_fuzz,breaks=min(progress[,1]):(max(progress[,1])+1),plot=F)$counts;
#
	progress <- cbind(progress,fcolls);
	if (update_search && tolower(update_type)=="graph")	{
		iteration <- ncol(progress);
		if (iteration>(length(all_colors)+1))
			all_colors <- rainbow(iteration);
		interation_cols <- all_colors[1:(iteration-1)];
#		points(3*(1:length(interation_cols)),rep(4,length(interation_cols)),pch=21,bg=interation_cols)
#		interation_cols <- paste("gray",,sep="");
		interation_lwds <- c(0.5,rep(0,iteration-2));
		for (pp in 1:nrow(progress))	{
			print_order <- 1+order(progress[pp,2:iteration],2:iteration,decreasing = T);
			for (ppp in 1:length(print_order))	{
				if (progress[pp,print_order[ppp]]>0)	{
					rect(pp-1,log10(0.7),pp,log10(progress[pp,print_order[ppp]]),col=interation_cols[print_order[ppp]-1],border="gray50",lwd=interation_lwds[ppp]);
#					} else if (max(progress[pp,2:iteration]>0))	{
#					rect(pp-1,log10(0.7),pp,log10(0.7)/2,col="white",border="white",lwd=1.5);
					}
				}
			}
		} else if (update_search) {
		progress <- cbind(progress[1:length(fcolls),],fcolls);
		print(progress);
		}

	unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];
	improved <- unfixed-length(unfixed_collections);
	unfixed <- length(unfixed_collections);
	improved <- sum(abs(progress[,ncol(progress)]-progress[,ncol(progress)-1]));
	list_lengths <- vector(length=unfixed);
	for (nc in 1:length(unfixed_collections))
		list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
	unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];
 	}
paleodb_collections <- put_pbdb_dataframes_into_proper_type(pbdb_data = paleodb_collections);
paleodb_finds <- put_pbdb_dataframes_into_proper_type(pbdb_data = paleodb_finds);

output <- list(paleodb_collections,paleodb_finds,progress);
names(output) <- c("ml_optimized_collections","ml_optimized_paleodb_finds","iterations");
return(output);
}

# do this for each sampling group!
accersi_strat_likelihoods_of_possible_site_ages <- function(site_info,index_ranges,bin_sites,precision=0.1)	{
#### add something else to date collections within single bins based on # species with 1…N finds
bin_sites_cum <- cumsum(bin_sites);

poss_ages <- seq(site_info$ma_lb,site_info$ma_ub,by=-precision);
# the possible finds will span 2+ bins; so, make sure that we get it correct.
poss_finds <- relv_time <- c();
for (si in site_info$bin_lb:site_info$bin_ub)	{
	if (si==site_info$bin_lb)	{
		relv_time <- c(relv_time,seq(round(site_info$ma_lb,2),round(finest_chronostrat$ma_ub[si],2),by=-precision));
		slc <- length(relv_time)-1;
		bn_p_slc <- (site_info$site_lb-bin_sites_cum[si-1])/slc;
	#	bn_p_slc <- (bin_sites[si]*(site_info$ma_lb-finest_chronostrat$ma_ub[si])/finest_chronostrat$span[si])/slc;
		poss_finds <- c(poss_finds,rep(bn_p_slc,slc));
		} else if (si==site_info$bin_ub)	{
		if (round(finest_chronostrat$ma_lb[si]-precision,2)==round(site_info$ma_ub,2))	{
			this_time <- round(finest_chronostrat$ma_lb[si]-precision,2);
			} else	{
			this_time <- seq((round(finest_chronostrat$ma_lb[si],2)-precision),round(site_info$ma_ub,2),by=-round(precision,2));
			}
		relv_time <- c(relv_time,this_time);
		slc <- length(this_time);
		bn_p_slc <- (site_info$site_ub-bin_sites_cum[si-1])/slc;
		#bn_p_slc <- (bin_sites[si]*abs(site_info$ma_ub-finest_chronostrat$ma_lb[si])/finest_chronostrat$span[si])/slc;
		poss_finds <- c(poss_finds,rep(bn_p_slc,slc))
		} else	{
		this_time <- seq(finest_chronostrat$ma_lb[si]-precision,finest_chronostrat$ma_ub[si],by=-precision);
		relv_time <- c(relv_time,this_time);
		slc <- length(this_time);
		bn_p_slc <- bin_sites[si]/slc;
		poss_finds <- c(poss_finds,rep(bn_p_slc,slc))
		}
	}

# possible cases:
lnl_ages <- array(0,dim=c(nrow(index_ranges),length(poss_ages)-1));
rownames(lnl_ages) <- rownames(index_ranges);
names(poss_finds) <- colnames(lnl_ages) <- (poss_ages[1:(length(poss_ages)-1)]+poss_ages[2:length(poss_ages)])/2;

ir <- 0;
for (ir in 1:nrow(index_ranges))	{
#while (ir < nrow(index_ranges))	{
#	ir <- ir+1;
	if (index_ranges$fa_imp[ir]>=site_info$ma_lb & index_ranges$la_imp[ir]<=site_info$ma_ub)	{
		# entirely within range of taxon
	#	lnl_ages[ir,] <- rbind(lnl_ages,rep(0,ncol(lnl_ages)));
		} else if (index_ranges$fa_imp[ir]<=site_info$ma_lb & index_ranges$la_imp[ir] <= site_info$ma_ub) {
		# possibly older than the taxon is known to be
		# use numbers of finds per slice (taking into account different bins) and work up to FA
		relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
		# introduce something here to estimate the gap at lowering rhos with extended earlier bounds
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		} else if (index_ranges$fa_imp[ir]>=site_info$ma_lb & index_ranges$la_imp[ir]<site_info$ma_lb)	{
		# possibly younger than the taxon is known to be
		# use numbers of finds per slice (taking into account different bins) and work up to FA
		relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		} else if (index_ranges$fa_imp[ir]>=site_info$ma_ub & site_info$ma_ub>index_ranges$la_imp[ir])	{
		# possibly older or younger than the taxon is known to be
		relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		# do younger finds
		relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		} else if (index_ranges$ma_fa_lb[ir]<site_info$ma_ub)	{
		# known finds all are older than this collection
		relv_finds_cum <- cumsum(poss_finds)
		lnl_ages[ir,] <- lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
#		print(ir);
		} else if (index_ranges$ma_la_ub[ir]>site_info$ma_lb)	{
		relv_finds_cum <- cumsum(poss_finds[length(poss_finds):1]);
		lnl_ages[ir,] <- relv_finds_cum[length(poss_finds):1]*log(1-index_ranges$rho[ir]);
		# known finds all are younger than this collection
		} else if (index_ranges$fa_imp[ir]<site_info$ma_lb & index_ranges$la_imp[ir]>site_info$ma_ub)	{
		# possibly older finds
		relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		# possibly younger finds
		relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		} else if (index_ranges$la_imp[ir]>site_info$ma_ub)	{
		# possibly younger finds
		relv_finds <- poss_finds[as.numeric(names(poss_finds))<index_ranges$la_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir])
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		} else if (index_ranges$fa_imp[ir]<site_info$ma_lb)	{
		# possibly older finds
		relv_finds <- poss_finds[as.numeric(names(poss_finds))>index_ranges$fa_imp[ir]];
		relv_finds_cum <- cumsum(relv_finds[order(names(relv_finds))]);
		lnl_spc <- relv_finds_cum*log(1-index_ranges$rho[ir]);
		lnl_ages[ir,match(names(lnl_spc),colnames(lnl_ages))] <- lnl_spc;
		}
#	print(lnl_ages[ir,])
	}
return(lnl_ages);
}

# score sum of range extensions required to make this the range
procrustean_binning_one_taxon <- function(range,poss_range,debug=FALSE)	{
# range: taxon range in bin numbers (with 1 older than 2)
# poss_range: lower and upper bounds on collection age in bin numbers
if (debug)	print(range);
obs_bins <- range[1]:range[2];
poss_bins <- poss_range[1]:poss_range[2];
#ok_bins <- poss_bins[poss_bins %in% obs_bins];
#if (length(ok_bins)==0)
#	ok_bins <- poss_bins;
bin_score <- c();
for (pb in 1:length(poss_bins))	{
	bin_score <- c(bin_score,min(abs(poss_bins[pb]-obs_bins)));
	}
return(bin_score);
}

#paleodb_finds=pbdb_finds;paleodb_collections=refined_collections_prov
ml_three_rate_optimization_paleodb_collection_and_occurrence_taxon_partitioned <- function(paleodb_finds,paleodb_collections,hierarchical_chronostrat,update_search=T,update_type="Graph",max_rho=0.1,precision=0.1)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- remove_duplicate_occurrences_paleodb(occurrences=paleodb_finds);

# make sure that key collections fields are characters, not factors
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fn <- 0;
fixed_names <- c();
while (fn < length(to_fix))	{
	fn <- fn+1;
	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]));
	}
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp."
paleodb_finds <- expello_indeterminate_species(paleodb_finds);

# bin taxa into sampling groups;
paleodb_finds$sampling_group <- sampling_groups$group[match(paleodb_finds$class,sampling_groups$classes)];
paleodb_finds$sampling_group[is.na(paleodb_finds$sampling_group)][paleodb_finds$genus[is.na(paleodb_finds$sampling_group)] %in% classic_inarticulate_brachiopod_taxa] <- "Inarticulata";
paleodb_finds <- subset(paleodb_finds,!is.na(paleodb_finds$sampling_group));
taxon_names <- sort(unique(paleodb_finds$accepted_name));
taxon_sampling_groups <- paleodb_finds$sampling_group[match(taxon_names,paleodb_finds$accepted_name)];
paleodb_collections <- paleodb_collections[paleodb_collections$collection_no %in% paleodb_finds$collection_no,];

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;
ncolls <- nrow(paleodb_collections);

# number taxa
paleodb_finds$taxon_no <- match(taxa,taxon_names);
#paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
if (is.na(match("ma_lb",colnames(paleodb_finds))))	{
	paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
	paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
	}
if (is.na(match("interval_lb",colnames(paleodb_finds))))	{
	paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
	paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);
	}

finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
finest_chronostrat <- unique(finest_chronostrat);

poss_ages <- sort(unique(c(paleodb_collections$ma_lb,paleodb_collections$ma_ub)),decreasing=T);

if (is.na(match("bin_lb",colnames(paleodb_collections))))	{
	paleodb_collections$bin_lb <- hierarchical_chronostrat$bin_first[match(unique(paleodb_collections$interval_lb),hierarchical_chronostrat$interval)];
	paleodb_collections$bin_ub <- hierarchical_chronostrat$bin_last[match(unique(paleodb_collections$interval_ub),hierarchical_chronostrat$interval)];
	paleodb_collections$bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
	paleodb_collections$bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)];
	}
#ugh <- (1:ncolls)[is.na(paleodb_collections$bin_ub)]
# kluge: eliminate the need for this;
xx <- (1:ncolls)[paleodb_collections$bin_lb>paleodb_collections$bin_ub];
xx <- xx[!is.na(xx)];
if (length(xx)>0)	{
	dummy <- paleodb_collections$bin_lb[xx];
	paleodb_collections$bin_lb[xx] <- paleodb_collections$bin_ub[xx];
	paleodb_collections$bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[paleodb_collections$bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[paleodb_collections$bin_ub[xx]];
	}

fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
control_groups <- sort(unique(sampling_groups$groups));
bin_sites_all <- array(0,dim=c(length(control_groups),length=max(paleodb_collections$bin_ub)));
rownames(bin_sites_all) <- control_groups;
colnames(bin_sites_all) <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
bin_sites_all_cum <- bin_sites_all;
for (cg in 1:length(control_groups))	{
	group_finds <- subset(paleodb_finds,paleodb_finds$sampling_group==control_groups[cg]);
	group_collections <- subset(fixed_collections,fixed_collections$collection_no %in% group_finds$collection_no);
	bin_sites_all[cg,1:max(fixed_collections$bin_lb)] <- hist(group_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
	bin_sites_all_cum[cg,] <- cumsum(bin_sites_all[cg,]);
	group_taxon_names <- taxon_names[taxon_sampling_groups==control_groups[cg]];
	if (cg==1)	{
		minimum_range_data <- accersi_minimum_range_information(taxon_names=group_taxon_names,
																paleodb_finds=group_finds,
																paleodb_collections=paleodb_collections,
																bin_sites=bin_sites_all[cg,],
																finest_chronostrat,max_rho=0.1,precision=precision);
		} else	{
		m_r_d <- accersi_minimum_range_information(taxon_names=group_taxon_names,
												 paleodb_finds=group_finds,
												 paleodb_collections=paleodb_collections,
												 bin_sites=bin_sites_all[cg,],
												 finest_chronostrat,max_rho=0.1,precision=precision);
		minimum_range_data <- rbind(minimum_range_data,m_r_d);
		}
	}
minimum_range_data <- minimum_range_data[order(rownames(minimum_range_data)),];

#minimum_range_data <- data.frame(ma_fa_mx=as.numeric(rep(0,ntaxa)),ma_fa_mn=as.numeric(rep(0,ntaxa)),ma_la_mx=as.numeric(rep(0,ntaxa)),ma_la_mn=as.numeric(rep(0,ntaxa)),interval_lb=as.character(rep("",ntaxa)),interval_ub=as.character(rep("",ntaxa)),N=as.numeric(rep(0,ntaxa)),Na=as.numeric(rep(0,ntaxa)),Nz=as.numeric(rep(0,ntaxa)),fa_imp=as.numeric(rep(0,ntaxa)),la_imp=as.numeric(rep(0,ntaxa)),site_lb=as.numeric(rep(0,ntaxa)),site_ub=as.numeric(rep(0,ntaxa)),range=as.numeric(rep(0,ntaxa)),rho=as.numeric(rep(0,ntaxa)),stringsAsFactors = F);
#rownames(minimum_range_data) <- taxon_names;
# this is putting everything in one bin too late: fix it
#minimum_range_data <- accersi_minimum_range_information(taxon_names,paleodb_finds,paleodb_collections,bin_sites,finest_chronostrat,max_rho=0.1,precision=precision);
#minimum_range_data_test <- subset(minimum_range_data,minimum_range_data$interval_lb!="");

bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search & tolower(update_type)=="graph")	{
	all_colors <- rainbow(7);
	mxx <- max(progress[,1])+1;
	mxy <- log10(ceiling(nrow(paleodb_collections)/10^floor(log10(nrow(paleodb_collections))))*10^floor(log10(nrow(paleodb_collections))));
	poss_bins <- 0:mxx;
	par(mfrow=c(1,1))
	par(pin=c(4.5,3));
	plot(NA,type='n',axes=FALSE,main="",xlab="Bin Precision",ylab="Collections",xlim=c(0,mxx),ylim=c(0,mxy));
	if (mxx<25)	{
		axis_labels <- 1:mxx
		} else	{
		axis_labels <- seq(2,mxx+1,by=2)-1;
		}
	specified_axis_w_labels(axe=1,max_val=mxx,min_val=0,maj_break=1,med_break=1,min_break=1,label_pos="mid",axis_labels=axis_labels,font_size=0.75);
	numbers <- base_numbers <- c(1,5);
	nn <- 0;
	while (nn < floor(mxy))	{
		nn <- nn+1;
		numbers <- c(numbers,base_numbers*10^nn)
		}
	numbers <- numbers[log10(numbers)<=mxy];
	log10_axes(axe=2,min_ax=0,max_ax=mxy,numbers,font_size=5/6,orient=2);
	for (pp in 1:mxx)	{
		if (progress[pp,2]>0)
			rect(pp-1,log10(0.7),pp,log10(progress[pp,2]),col=all_colors[1],border="gray50",lwd=0.5);
		}
	} else if (update_search)	{
	print(progress);
	}

#collection_nos <- sort(paleodb_collections$collection_no);
list_lengths <- vector(length=length(unfixed_collections));
for (nc in 1:length(unfixed_collections))
	list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];

finest_chronostrat$span <- abs(finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- uc <- 0;	# uc <- match(problems[12],unfixed_collections)
	modified_species <- c();
	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges_relv <- subset(index_ranges_all,index_ranges_all$rho!=0);

		if (nrow(index_ranges_relv)>1)	{
			relv_control_groups <- unique(taxon_sampling_groups[match(rownames(index_ranges_relv),taxon_names)]);
			r_c_g <- match(relv_control_groups,control_groups);
			rcg <- length(r_c_g);
			site_info <- s_inf <- paleodb_collections[coll_no,colnames(paleodb_collections) %in% c("ma_lb","ma_ub","bin_lb","bin_ub","interval_lb","interval_ub")];
			bl <- site_info$bin_lb;
			bu <- site_info$bin_ub;
			lnl_ages <- c();
			for (tcg in 1:rcg)	{
				cg <- r_c_g[tcg];	# cg is the control group number
				prop_lower_bin <- (finest_chronostrat$ma_lb[bl]-site_info$ma_lb)/finest_chronostrat$span[bl];		# prop of span from bottom to top
				prop_upper_bin <- abs(finest_chronostrat$ma_ub[bu]-site_info$ma_ub)/finest_chronostrat$span[bu];	# prop of span from top to bottom
				if (bl>0)	{
					site_info$site_lb <- floor(bin_sites_all_cum[cg,bl-1]+(prop_lower_bin*bin_sites_all[cg,bl]));
					} else	{
					site_info$site_lb <- floor(prop_lower_bin*bin_sites_all[cg,bl]);
					}
				site_info$site_ub <- ceiling(bin_sites_all_cum[cg,bu] - (prop_upper_bin*bin_sites_all[cg,bu]));
				bin_sites <- bin_sites_all[cg,];
				index_ranges <- index_ranges_relv[match(taxon_sampling_groups[match(rownames(index_ranges_relv),taxon_names)],control_groups)==cg,];

				lnl_ages <- rbind(lnl_ages,accersi_strat_likelihoods_of_possible_site_ages(site_info=site_info,index_ranges=index_ranges,bin_sites=bin_sites_all[cg,],precision=precision));
				# Now, do likelihoods for each group #
				}

			age_ln_likelihoods <- colSums(lnl_ages)-max(colSums(lnl_ages));
			new_age_range <- as.numeric(names(age_ln_likelihoods)[age_ln_likelihoods>=-2]);
			if (length(new_age_range)<length(age_ln_likelihoods))	{
				new_and_improved <- new_and_improved+1;
				paleodb_collections$ma_lb[coll_no] <- ceiling(max(new_age_range)/precision)*precision;
				paleodb_collections$ma_ub[coll_no] <- floor(min(new_age_range)/precision)*precision;
				paleodb_collections$interval_lb[coll_no] <- rebin_collection_with_time_scale(max(new_age_range),"onset",finest_chronostrat);
				paleodb_collections$interval_ub[coll_no] <- rebin_collection_with_time_scale(min(new_age_range),"end",finest_chronostrat);
				paleodb_collections$bin_lb[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb[coll_no],finest_chronostrat$interval)];
				paleodb_collections$bin_ub[coll_no] <- finest_chronostrat$bin_first[match(paleodb_collections$interval_ub[coll_no],finest_chronostrat$interval)];
				bin_fuzz[coll_no] <- paleodb_collections$bin_ub[coll_no]-paleodb_collections$bin_lb[coll_no];
				if (bin_fuzz[coll_no]==0)
					fixed <- fixed+1;
				paleodb_finds$ma_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_lb[coll_no];
				paleodb_finds$ma_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$ma_ub[coll_no];
				paleodb_finds$interval_lb[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_lb[coll_no];
				paleodb_finds$interval_ub[paleodb_finds$collection_no==unfixed_collections[uc]] <- paleodb_collections$interval_ub[coll_no];
				modified_species <- unique(c(modified_species,index_species));
				} # finish updating info
			} # end case of having set species;
 		} # end search of unconstrained_collections

	fixed_collections <- subset(paleodb_collections,paleodb_collections$bin_lb==paleodb_collections$bin_ub);
	bin_sites_old <- bin_sites_all;
	bin_sites_cum_old <- bin_sites_all_cum;
#	bin_sites_all <- array(0,dim=c(length(control_groups),length=max(paleodb_collections$bin_ub)));
#	rownames(bin_sites_all) <- control_groups;
#	colnames(bin_sites_all) <- finest_chronostrat$interval[1:max(paleodb_collections$bin_ub)];
#	bin_sites_all_cum <- bin_sites_all;
	for (cg in 1:length(control_groups))	{
		group_finds <- subset(paleodb_finds,paleodb_finds$sampling_group==control_groups[cg]);
		group_collections <- subset(fixed_collections,fixed_collections$collection_no %in% group_finds$collection_no);
		bin_sites_all[cg,1:max(fixed_collections$bin_lb)] <- hist(group_collections$bin_lb,breaks=0:max(fixed_collections$bin_lb),plot=F)$counts;
		bin_sites_all_cum[cg,] <- cumsum(bin_sites_all[cg,]);
		group_taxon_names <- taxon_names[taxon_sampling_groups==control_groups[cg]];
		if (cg==1)	{
			minimum_range_data <- accersi_minimum_range_information(taxon_names=group_taxon_names,
																	paleodb_finds=group_finds,
																	paleodb_collections=paleodb_collections,
																	bin_sites=bin_sites_all[cg,],
																	finest_chronostrat,max_rho=0.1,precision=precision);
			} else	{
			m_r_d <- accersi_minimum_range_information(taxon_names=group_taxon_names,
													   paleodb_finds=group_finds,
													   paleodb_collections=paleodb_collections,
													   bin_sites=bin_sites_all[cg,],
													   finest_chronostrat,max_rho=0.1,precision=precision);
			minimum_range_data <- rbind(minimum_range_data,m_r_d);
			}
		}
	minimum_range_data <- minimum_range_data[order(rownames(minimum_range_data)),];

	fuzz <- 0:max(paleodb_collections$bin_ub-paleodb_collections$bin_lb);
 	bin_fuzz <- paleodb_collections$bin_ub-paleodb_collections$bin_lb;
	fcolls <- hist(bin_fuzz,breaks=min(progress[,1]):(max(progress[,1])+1),plot=F)$counts;
#
	progress <- cbind(progress,fcolls);
	if (update_search && tolower(update_type)=="graph")	{
		iteration <- ncol(progress);
		if (iteration>(length(all_colors)+1))
			all_colors <- rainbow(iteration);
		interation_cols <- all_colors[1:(iteration-1)];
#		points(3*(1:length(interation_cols)),rep(4,length(interation_cols)),pch=21,bg=interation_cols)
#		interation_cols <- paste("gray",,sep="");
		interation_lwds <- c(0.5,rep(0,iteration-2));
		for (pp in 1:nrow(progress))	{
			print_order <- 1+order(progress[pp,2:iteration],2:iteration,decreasing = T);
			for (ppp in 1:length(print_order))	{
				if (progress[pp,print_order[ppp]]>0)	{
					rect(pp-1,log10(0.7),pp,log10(progress[pp,print_order[ppp]]),col=interation_cols[print_order[ppp]-1],border="gray50",lwd=interation_lwds[ppp]);
#					} else if (max(progress[pp,2:iteration]>0))	{
#					rect(pp-1,log10(0.7),pp,log10(0.7)/2,col="white",border="white",lwd=1.5);
					}
				}
			}
		} else if (update_search) {
		progress <- cbind(progress[1:length(fcolls),],fcolls);
		print(progress);
		}

	unfixed_collections <- paleodb_collections$collection_no[paleodb_collections$bin_ub!=paleodb_collections$bin_lb];
	improved <- unfixed-length(unfixed_collections);
	unfixed <- length(unfixed_collections);
	improved <- sum(abs(progress[,ncol(progress)]-progress[,ncol(progress)-1]));
	list_lengths <- vector(length=unfixed);
	for (nc in 1:length(unfixed_collections))
		list_lengths[nc] <- sum(paleodb_finds$collection_no==unfixed_collections[nc]);
	unfixed_collections <- unfixed_collections[order(list_lengths,decreasing=T)];
 	}
output <- list(paleodb_collections,paleodb_finds,progress);
names(output) <- c("ml_optimized_collections","ml_optimized_paleodb_finds","iterations");
return(output);
}
