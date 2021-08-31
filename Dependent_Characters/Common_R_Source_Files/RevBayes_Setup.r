#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal
# scribere: write

#library(ggtree);		# install_github("GuangchuangYu/ggtree");
#library(RevGadgets);	# install_github("revbayes/RevGadgets");
ZERO <- 1e-323;
MINEXPN <- 10^-10;
MINNO <- 5e-324;
MAXNO <- 1.797693e+308;
UNKNOWN <- -11;
INAP <- -22;
hell_no <- F;
uncertains <- c("cf.","aff.","ex_gr.");
newick_verbotten <- c(".","?","\"","\'");
letter_states <- LETTERS[!LETTERS %in% c("I","O")];
zzzz <- 0.25;

##### ROUTINES TO WRITE REVBAYES SCRIPTS #######
# divide matrix into multiple matrices with all characters having N=2…? states
#	invariants are placed with 2 state characters
divido_character_matrix_by_state_numbers <- function(nexus_file_name,data_directory,script_directory,output_file_lead="output/",script_file_lead="script/",analysis_name="",no_runs=4,write_rev_bayes_source=TRUE,polymorphs=TRUE, UNKNOWN=-11, INAP=-22)	{
orig_nexus_file_name <- nexus_file_name;
if (analysis_name=="")	{
	j <- strsplit(orig_nexus_file_name,split="",fixed=TRUE)[[1]];
	j <- j[1:(max((1:length(j))[j %in% "."])-1)];
	analysis_name <- paste(j,collapse="");
	}

if (data_directory!="")
	nexus_file_name <- paste(data_directory,nexus_file_name,sep="");
initial_data <- accersi_data_from_nexus_file(nexus_file_name,polymorphs,UNKNOWN,INAP);

n_states <- initial_data$States;
n_states[n_states<2] <- 2;
taxon_names <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
rownames(chmatrix) <- initial_data$OTUs;
n_chars <- length(n_states);
state_orders <- initial_data$State_Order;
if (write_rev_bayes_source)	{
	nstates <- sort(unique(n_states));
	state_numbers <- state_ordering <- matrix_file_names <- c();
	}
#for (st in 1:length(state_numbers))	{
st <- 0;
while (st < length(nstates))	{
	st <- st+1;
	relv_chars <- (1:n_chars)[n_states==nstates[st]];
	orderings <- unique(state_orders[relv_chars]);
	for (sto in 1:length(orderings))	{
		relv_chars <- (1:n_chars)[n_states==nstates[st]];
		relv_chars <- relv_chars[state_orders[relv_chars]==orderings[sto]];
		ch_matrix <- initial_data$Matrix[,relv_chars];
		new_file_name <- paste(analysis_name,"_Matrix_",nstates[st],"_States_",orderings[sto],".nex",sep="");
		if (data_directory!="")
			new_file_name <- paste(data_directory,new_file_name,sep="");
		if (!is.matrix(ch_matrix))	{
			ch_matrix <- matrix(ch_matrix);
			rownames(ch_matrix) <- taxon_names;
			}
		scribio_nexus_file_from_chmatrix(ch_matrix,new_file_name,UNKNOWN,INAP);
		if (write_rev_bayes_source)	{
			matrix_file_names <- c(matrix_file_names,new_file_name);
			state_numbers <- c(state_numbers,nstates[st]);
			state_ordering <- c(state_ordering,orderings[sto]);
			}
		}
	}

if (write_rev_bayes_source)	{
	filename <- paste(paste(script_directory,analysis_name,sep=""),"_Partitioned_Analysis.Rev",sep="");
	revbayes_source <- "clear();"
	revbayes_source <- c(revbayes_source,"");
	file_name_string <- paste(matrix_file_names, collapse = "\", \"");
	file_name_string <- paste("filenames <- v(\"",file_name_string,"\");",sep="");
	file_name_string <-gsub("~/","",file_name_string);
	revbayes_source <- c(revbayes_source,"# Make sure that filenames & directories are correct!!!");
	revbayes_source <- c(revbayes_source,file_name_string);
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,paste("partition_states <- v(",paste(state_numbers, collapse = ","),");",sep=""));
	revbayes_source <- c(revbayes_source,paste("partition_ordering <- v(\"",paste(state_ordering, collapse = "\",\""),"\");",sep=""));
#	state_numbers;
#	state_ordering;
	if(as.numeric(initial_data$Outgroup)!=-1)	{
		if (length(initial_data$Outgroup)>1)	{
			outies <- paste(taxon_names[as.numeric(initial_data$Outgroup)],",");
			outies <- paste("v(",outies,"\")",sep="");
			} else	{
			outies <- taxon_names[as.numeric(initial_data$Outgroup)];
			}
		} else	{
		outies <- "ENTER AN OUTGROUP HERE!"
		}
	revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,"# Set a variable for the number of  partitions");
	revbayes_source <- c(revbayes_source,"n_data_subsets <- filenames.size();");
	revbayes_source <- c(revbayes_source,"dummy <- readDiscreteCharacterData(filenames[1]);");
	revbayes_source <- c(revbayes_source,"taxa <- dummy.taxa();");
	revbayes_source <- c(revbayes_source,"n_species <- dummy.ntaxa();");
	revbayes_source <- c(revbayes_source,"n_branches <- 2 * n_species - 3;");
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,"# Set up tree-search moves");
	revbayes_source <- c(revbayes_source,"topology ~ dnUniformTopology(taxa,outgroup);");
	revbayes_source <- c(revbayes_source,"mvi = 0;  # count moves");
	revbayes_source <- c(revbayes_source,"moves[++mvi] = mvNNI(topology, weight=1.0);   # nearest neighbor interchange");
	revbayes_source <- c(revbayes_source,"moves[++mvi] = mvSPR(topology, weight=1.0);   # subtree pruning");
	revbayes_source <- c(revbayes_source,"for (b in 1:n_branches) {");
	revbayes_source <- c(revbayes_source,"    bl[b] ~ dnExponential(10.0);");
	revbayes_source <- c(revbayes_source,"    moves[++mvi] = mvScale(bl[b]);");
 	revbayes_source <- c(revbayes_source,"   }");
	revbayes_source <- c(revbayes_source,"tau := treeAssembly(topology, bl);  # assign branch lengths to trees");
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,"# Set up appropriate Q-matrices for the partitions");
	revbayes_source <- c(revbayes_source,"#\t(Again, make sure that the directory is OK)");
	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Parameters_for_Analysis_Partitioned_by_States_and_Ordering.Rev\");",sep=""));
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,"# add monitors & commence MCMC'ing");
	revbayes_source <- c(revbayes_source,"#\t(Again, make sure that the source directory is OK)");
	revbayes_source <- c(revbayes_source,"mni = 0;  # count monitors");
#	revbayes_source <- c(revbayes_source,paste("monitors[++mni] = mnModel(filename=\"Documents/RevStudio_Projects/output/",tolower(analysis_name),".log\",printgen=1000,separator=TAB);\t\t#Play with these numbers",sep=""));
#	revbayes_source <- c(revbayes_source,paste("monitors[++mni] = mnFile(tau, filename=\"Documents/RevStudio_Projects/output/",tolower(analysis_name),".trees\",printgen=1000,separator=TAB,tau);\t\t#Play with these numbers",sep=""));
	mcmc_output <- paste(output_file_lead,tolower(analysis_name),sep="");
	revbayes_source <- c(revbayes_source,paste("monitors[++mni] = mnModel(filename=\"",mcmc_output,".log\",printgen=1000,separator=TAB);\t\t#Play with these numbers",sep=""));
	revbayes_source <- c(revbayes_source,paste("monitors[++mni] = mnFile(tau, filename=\"",mcmc_output,".trees\",printgen=1000,separator=TAB,tau);\t\t#Play with these numbers",sep=""));
	revbayes_source <- c(revbayes_source,"monitors[++mni] = mnScreen(printgen=1000);");
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,paste("no_runs=",no_runs,";",sep=""));
	revbayes_source <- c(revbayes_source,"burnin_gens=10000;");
	revbayes_source <- c(revbayes_source,"tuning_int=200;");
	revbayes_source <- c(revbayes_source,"running_gens=1000000;");
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Expecto_MCMC_with_Partitioned_Characters.Rev\");",sep=""));
	revbayes_source <- c(revbayes_source,"");
	
	analyses_outputs <- paste(output_file_lead,tolower(analysis_name),sep="");
	if (no_runs>1)	{
		trees <- maj_rule <- max_prob_tree <- c();
		for (nr in 1:no_runs)	{
			trees <- c(trees,paste("\"",analyses_outputs,"_run_",nr,".tree\"",sep=""));
#			trees <- paste(trees,",",paste("\"",output_file_lead,tolower(analysis_name),"_run_",nr,".tree\"",sep=""));
			maj_rule <- c(maj_rule,paste("\"",analyses_outputs,"_run_",nr,"_maj_rule.tre\"",sep=""));
			max_prob_tree <- c(max_prob_tree,paste("\"",analyses_outputs,"_run_",nr,"_simple_map.tre\"",sep=""));
			}
		trees <- paste(trees,collapse=",");
		maj_rule <- paste(maj_rule,collapse=",");
		max_prob_tree <- paste(max_prob_tree,collapse=",");
		revbayes_source <- c(revbayes_source,paste("tree_files <- v(",trees,");",sep=""));
		revbayes_source <- c(revbayes_source,paste("maj_rule_files <- v(",maj_rule,");",sep=""));
		revbayes_source <- c(revbayes_source,paste("most_probable_files <- v(",max_prob_tree,");",sep=""));
		} else	{
		trees <- paste("\"",analyses_outputs,".trees\"",sep="");
		maj_rule <- paste("\"",analyses_outputs,"_maj_rule.tre\"",sep="");
		max_prob_tree <- paste("\"",analyses_outputs,"_simple_map.tre\"",sep="");
		revbayes_source <- c(revbayes_source,paste("tree_files <- ",trees,";",sep=""));
		revbayes_source <- c(revbayes_source,paste("maj_rule_files <- ",maj_rule,";",sep=""));
		revbayes_source <- c(revbayes_source,paste("most_probable_files <- ",max_prob_tree,";",sep=""));
		}
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"accersi_Consensus.Rev\");",sep=""));
	revbayes_source <- c(revbayes_source,"");
	revbayes_source <- c(revbayes_source,"# Now, go binge-watch Game of Thrones while you wait!");
	write(revbayes_source,file=filename);
#	tree_files <- v("output/acanthoceratoidea_run_1.tree","output/acanthoceratoidea_run_2.tree","output/acanthoceratoidea_run_3.tree","output/acanthoceratoidea_run_4.tree");
#	maj_rule_files <- v("output/acanthoceratoidea_run_1_maj_rule.tre","output/acanthoceratoidea_run_2_maj_rule.tre","output/acanthoceratoidea_run_3_maj_rule.tre","output/acanthoceratoidea_run_4_maj_rule.tre");
#	most_probable_files <- v("output/acanthoceratoidea_run_1_simple_map.tre","output/acanthoceratoidea_run_2_simple_map.tre","output/acanthoceratoidea_run_3_simple_map.tre","output/acanthoceratoidea_run_4_simple_map.tre");
	}
}

divido_character_matrix_by_state_numbers_and_other_partitions <- function(analysis_name="",first_appearances=NULL,write_data_directory="",rate_partition="",trend_partition="",taxa_subset="",data_file_lead="data/",polymorphs=T, UNKNOWN=-11, INAP=-22,coding_threshold=1)	{
# rate_partition: name of acharacter group that should have separate rates (e.g., CHARPARTITION Functional_Partitions  =  Nonfeeding :  1- 45 58- 60, Feeding :  46- 57; in a Nexus file)
# rattrend_partition: name of character group that for which one group ("Driven") has biased change (e.g., Driven_Trend  =  Unbiased :  1- 32 34- 40 43- 50 52- 60, Driven :  33 41- 42 51; in a Nexus file)
# taxa_subset: vector listing the subset of species to include in the final matrices;
# write_data_directory: the directory to which the partitioned matrices should go
# data_file_lead: the directory addendum to paste in front of data files so that RevBayes can find them
if (analysis_name=="")	{
	analysis_name <- "Wombat_Rock";
	}
#initial_data <- accersi_data_from_nexus_file(nexus_file_name=paste(read_data_directory,nexus_file_name,sep=""),polymorphs,UNKNOWN,INAP,character_partitions=rate_partition);
#if (rate_partition!="")	{
print("Your Nexus file probably needs to be converted to RevBayes' preferred format:");
flush.console();
initial_data <- accersi_data_from_chosen_nexus_file(polymorphs,UNKNOWN,INAP,rate_partition,trend_partition);

taxon_names <- initial_data$OTUs;
taxon_names <- gsub("\\(\\?\\)","",taxon_names);
taxon_names <- gsub("\\?","",taxon_names);
taxon_names <- gsub("n\\. sp\\. ","",taxon_names);
taxon_names <- gsub("n\\. gen\\. ","",taxon_names);
taxon_names <- gsub("  "," ",taxon_names);
initial_data$OTUs <- taxon_names;

# added 2020-12-17 to prevent unscored characters from messing up things
# cut underscored characters
scored_taxa <- count_scored_otus_per_character(chmatrix=initial_data$Matrix,UNKNOWN,INAP);
chmatrix <- initial_data$Matrix[,scored_taxa>=coding_threshold];
n_states <- initial_data$States[scored_taxa>=coding_threshold];
n_chars <- length(n_states);
state_types <- initial_data$State_Types[scored_taxa>=coding_threshold];
n_states[n_states<2] <- 2;
rownames(chmatrix) <- initial_data$OTUs;
notu <- nrow(chmatrix);
if (!is.null(first_appearances))	{
	if (is.data.frame(first_appearances) && !is.null(first_appearances$taxon))	{
		taxon_order <- match(gsub("_"," ",first_appearances$taxon),rownames(chmatrix));
#		taxon_order <- match(first_appearances$taxon,rownames(chmatrix));
#		if (sum(is.na(taxon_order))==notu)
			
		if (sum(!is.na(taxon_order))>3 && sum(!is.na(taxon_order)<notu))	{
			tx_ord <- rank(taxon_order[!is.na(taxon_order)]);
			dec <- 0;
			for (i in 2:length(tx_ord))
				if (tx_ord[i]<tx_ord[i-1])
					dec <- dec+1;
			if (dec==0)	taxon_order <- 1:notu;
			}
		if (length(taxon_order)==notu)	{
			otu_fas <- -abs(first_appearances$fa[taxon_order]);
			chmatrix <- rescore_character_matrix_by_first_appearances(chmatrix,otu_fas);
			}
		}
	}

#if (write_rev_bayes_source)	{
nstates <- sort(unique(n_states));
rate_partitions <- initial_data$Rate_Partitions;
trend_partitions <- initial_data$Trend_Partitions;

# separate out all unique combinations that we might use to partition the characters
if (length(unique(rate_partitions))>1 && length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),rate_partitions=as.character(rate_partitions),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else if (length(unique(rate_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),rate_partitions=as.character(rate_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else if (length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	}

if (taxa_subset!="" && length(taxa_subset)<notu)	{
	chmatrix <- chmatrix[match(taxa_subset,rownames(chmatrix)),];
	notu <- nrow(chmatrix);
	}

pc <- 0;
partition_size <- character_rate_partitions <- character_trend_partitions <- matrix_file_names <- coding_bias <- c();
while (pc < nrow(partition_combos))	{
	# check coding_bias here!!!!
	pc <- pc+1;
	# find characters with the appropriate numbers of states
	relv_characters <- (1:n_chars)[n_states %in% partition_combos$n_states[pc]];
	partition_size <- c(partition_size,length(relv_characters));
	if (length(unique(state_types))>1)	{
		# find characters with this ordering
		right_ordering <- (1:n_chars)[state_types %in% partition_combos$type[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_ordering];
		}
	if (length(unique(rate_partitions))>1)	{
		# find characters in this rate partitioning
		right_partition <- (1:n_chars)[rate_partitions %in% partition_combos$rate_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	if (length(unique(trend_partitions))>1)	{
		# find characters in this trend partitioning
		right_partition <- (1:n_chars)[trend_partitions %in% partition_combos$trend_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	# look for invariant & autapomorphic characters
	chmatrix_red <- chmatrix[,relv_characters];
	if (!is.matrix(chmatrix_red))
		chmatrix_red <- as.matrix(chmatrix_red);
	## RevBayes ignores polymorphs; so, if the only cases of derived states are in polymorphic taxa, recode them to derived states;
	observed_states_sans_poly <- count_states(chmatrix_red,UNKNOWN,INAP,include_polymorphs = F);
	poss_polymorph_fixs <- (1:ncol(chmatrix_red))[observed_states_sans_poly<partition_combos$n_states[pc]];
	ppf <- 0;
	while (ppf < length(poss_polymorph_fixs))	{
		ppf <- ppf+1;
		nch <- poss_polymorph_fixs[ppf];
		poly_otus <- (1:notu)[!chmatrix_red[,nch] %in% c(INAP,UNKNOWN,0:partition_combos$n_states[pc])];
		orig_state <- obs_states <- sort(unique(chmatrix_red[chmatrix_red[,nch]>=0,nch]));
		potu <- 0;
		while (potu < length(poly_otus))	{
			potu <- potu+1;
			polystates <- unravel_polymorph(chmatrix_red[poly_otus[potu],nch]);
			if (length(polystates[!polystates %in% obs_states])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- max(polystates[!polystates %in% obs_states]);
				} else if (length(polystates[!polystates %in% orig_state])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- min(polystates[!polystates %in% orig_state]);
				}
			}
		}
	observed_states <- count_states(chmatrix_red,UNKNOWN,INAP);
	
	if (sum(observed_states==1)>0)	{
		coding_bias <- c(coding_bias,"all");
		} else	{
		chstates <- initial_data$States;
		autaps <- list_autapomorphic_characters(chmatrix,chstates,UNKNOWN,INAP);
		if (length(autaps)>0 || partition_combos$n_states[pc]>2)	{
			coding_bias <- c(coding_bias,"variable");
			} else	{
			coding_bias <- c(coding_bias,"informative");
			}
		}
	new_file_name <- paste(analysis_name,"_Matrix_",partition_combos$n_states[pc],"_States",sep="");
	if (length(unique(state_types))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$type[pc],sep="");
	if (length(unique(rate_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$rate_partitions[pc],sep="");
	if (length(unique(trend_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$trend_partitions[pc],sep="");
	orig_file_name <- new_file_name <- paste(new_file_name,".nex",sep="");
	if (write_data_directory!="")
		new_file_name <- paste(write_data_directory,new_file_name,sep="");
	state_symbols <- accersi_state_symbols(n_states=partition_combos$n_states[pc]);
#	state_symbols <- (1:partition_combos$n_states[pc])-1;
#	state_symbols[state_symbols>=10] <- letter_states[state_symbols[state_symbols>=10]-9];
	scribere_rev_bayes_nexus_file_from_character_matrix(ch_matrix=chmatrix_red,state_symbols = state_symbols,new_file_name=new_file_name,UNKNOWN,INAP);
	matrix_file_names <- c(matrix_file_names,paste(data_file_lead,orig_file_name,sep=""));
	}

state_numbers <- partition_combos$n_states;
state_ordering <- partition_combos$type;
if (!is.null(partition_combos$rate_partitions))	{
	character_rate_partitions <- partition_combos$rate_partitions;
	} else	{
	character_rate_partitions <- rep("imagine",nrow(partition_combos));
	}
if (!is.null(partition_combos$trend_partitions))	{
	character_trend_partitions <- partition_combos$trend_partitions;
	} else	{
	character_trend_partitions <- rep("square",nrow(partition_combos));
	}

output <- list(initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,character_rate_partitions,character_trend_partitions,coding_bias);
names(output) <- c("initial_data","matrix_file_names","partition_size","state_numbers","state_ordering","rate_partitions","trend_partitions","coding_bias");
return(output);
}

divido_nexus_file_by_state_numbers_and_other_partitions <- function(nexus_file_name,analysis_name="",first_appearances=NULL,write_data_directory="",rate_partition="",trend_partition="",taxa_subset="",data_file_lead="data/",polymorphs=T, UNKNOWN=-11, INAP=-22,coding_threshold=1)	{
# rate_partition: name of acharacter group that should have separate rates (e.g., CHARPARTITION Functional_Partitions  =  Nonfeeding :  1- 45 58- 60, Feeding :  46- 57; in a Nexus file)
# rattrend_partition: name of character group that for which one group ("Driven") has biased change (e.g., Driven_Trend  =  Unbiased :  1- 32 34- 40 43- 50 52- 60, Driven :  33 41- 42 51; in a Nexus file)
# taxa_subset: vector listing the subset of species to include in the final matrices;
# write_data_directory: the directory to which the partitioned matrices should go
# data_file_lead: the directory addendum to paste in front of data files so that RevBayes can find them
if (analysis_name=="")	{
	analysis_name <- "Wombat_Rock";
	}
#initial_data <- accersi_data_from_nexus_file(nexus_file_name=paste(read_data_directory,nexus_file_name,sep=""),polymorphs,UNKNOWN,INAP,character_partitions=rate_partition);
#if (rate_partition!="")	{
initial_data <- accersi_data_from_nexus_file(nexus_file_name,polymorphs,UNKNOWN,INAP,rate_partition,trend_partition);
#accersi_data_from_chosen_nexus_file(polymorphs,UNKNOWN,INAP,rate_partition,trend_partition);

taxon_names <- initial_data$OTUs;
taxon_names <- gsub("\\(\\?\\)","",taxon_names);
taxon_names <- gsub("\\?","",taxon_names);
taxon_names <- gsub("n\\. sp\\. ","",taxon_names);
taxon_names <- gsub("n\\. gen\\. ","",taxon_names);
taxon_names <- gsub("  "," ",taxon_names);
initial_data$OTUs <- taxon_names;

# added 2020-12-17 to prevent unscored characters from messing up things
# cut underscored characters
scored_taxa <- count_scored_otus_per_character(chmatrix=initial_data$Matrix,UNKNOWN,INAP);
chmatrix <- initial_data$Matrix[,scored_taxa>=coding_threshold];
n_states <- initial_data$States[scored_taxa>=coding_threshold];
n_chars <- length(n_states);
state_types <- initial_data$State_Types[scored_taxa>=coding_threshold];
n_states[n_states<2] <- 2;
rownames(chmatrix) <- initial_data$OTUs;
notu <- nrow(chmatrix);
if (!is.null(first_appearances))	{
	if (is.data.frame(first_appearances) && !is.null(first_appearances$taxon))	{
		taxon_order <- match(gsub("_"," ",first_appearances$taxon),rownames(chmatrix));
#		taxon_order <- match(first_appearances$taxon,rownames(chmatrix));
#		if (sum(is.na(taxon_order))==notu)
			
		if (sum(!is.na(taxon_order))>3 && sum(!is.na(taxon_order)<notu))	{
			tx_ord <- rank(taxon_order[!is.na(taxon_order)]);
			dec <- 0;
			for (i in 2:length(tx_ord))
				if (tx_ord[i]<tx_ord[i-1])
					dec <- dec+1;
			if (dec==0)	taxon_order <- 1:notu;
			}
		if (length(taxon_order)==notu)	{
			otu_fas <- -abs(first_appearances$fa[taxon_order]);
			chmatrix <- rescore_character_matrix_by_first_appearances(chmatrix,otu_fas);
			}
		}
	}

#if (write_rev_bayes_source)	{
nstates <- sort(unique(n_states));
rate_partitions <- initial_data$Rate_Partitions;
trend_partitions <- initial_data$Trend_Partitions;

# separate out all unique combinations that we might use to partition the characters
if (length(unique(rate_partitions))>1 && length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),rate_partitions=as.character(rate_partitions),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else if (length(unique(rate_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),rate_partitions=as.character(rate_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else if (length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	} else	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),type=as.character(state_types),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$type),];
	}

if (taxa_subset!="" && length(taxa_subset)<notu)	{
	chmatrix <- chmatrix[match(taxa_subset,rownames(chmatrix)),];
	notu <- nrow(chmatrix);
	}

pc <- 0;
partition_size <- character_rate_partitions <- character_trend_partitions <- matrix_file_names <- coding_bias <- c();
while (pc < nrow(partition_combos))	{
	# check coding_bias here!!!!
	pc <- pc+1;
	# find characters with the appropriate numbers of states
	relv_characters <- (1:n_chars)[n_states %in% partition_combos$n_states[pc]];
	partition_size <- c(partition_size,length(relv_characters));
	if (length(unique(state_types))>1)	{
		# find characters with this ordering
		right_ordering <- (1:n_chars)[state_types %in% partition_combos$type[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_ordering];
		}
	if (length(unique(rate_partitions))>1)	{
		# find characters in this rate partitioning
		right_partition <- (1:n_chars)[rate_partitions %in% partition_combos$rate_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	if (length(unique(trend_partitions))>1)	{
		# find characters in this trend partitioning
		right_partition <- (1:n_chars)[trend_partitions %in% partition_combos$trend_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	# look for invariant & autapomorphic characters
	chmatrix_red <- chmatrix[,relv_characters];
	if (!is.matrix(chmatrix_red))
		chmatrix_red <- as.matrix(chmatrix_red);
	## RevBayes ignores polymorphs; so, if the only cases of derived states are in polymorphic taxa, recode them to derived states;
	observed_states_sans_poly <- count_states(chmatrix_red,UNKNOWN,INAP,include_polymorphs = F);
	poss_polymorph_fixs <- (1:ncol(chmatrix_red))[observed_states_sans_poly<partition_combos$n_states[pc]];
	ppf <- 0;
	while (ppf < length(poss_polymorph_fixs))	{
		ppf <- ppf+1;
		nch <- poss_polymorph_fixs[ppf];
		poly_otus <- (1:notu)[!chmatrix_red[,nch] %in% c(INAP,UNKNOWN,0:partition_combos$n_states[pc])];
		orig_state <- obs_states <- sort(unique(chmatrix_red[chmatrix_red[,nch]>=0,nch]));
		potu <- 0;
		while (potu < length(poly_otus))	{
			potu <- potu+1;
			polystates <- unravel_polymorph(chmatrix_red[poly_otus[potu],nch]);
			if (length(polystates[!polystates %in% obs_states])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- max(polystates[!polystates %in% obs_states]);
				} else if (length(polystates[!polystates %in% orig_state])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- min(polystates[!polystates %in% orig_state]);
				}
			}
		}
	observed_states <- count_states(chmatrix_red,UNKNOWN,INAP);
	
	if (sum(observed_states==1)>0)	{
		coding_bias <- c(coding_bias,"all");
		} else	{
		chstates <- initial_data$States;
		autaps <- list_autapomorphic_characters(chmatrix,chstates,UNKNOWN,INAP);
		if (length(autaps)>0 || partition_combos$n_states[pc]>2)	{
			coding_bias <- c(coding_bias,"variable");
			} else	{
			coding_bias <- c(coding_bias,"informative");
			}
		}
	new_file_name <- paste(analysis_name,"_Matrix_",partition_combos$n_states[pc],"_States",sep="");
	if (length(unique(state_types))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$type[pc],sep="");
	if (length(unique(rate_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$rate_partitions[pc],sep="");
	if (length(unique(trend_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$trend_partitions[pc],sep="");
	orig_file_name <- new_file_name <- paste(new_file_name,".nex",sep="");
	if (write_data_directory!="")
		new_file_name <- paste(write_data_directory,new_file_name,sep="");
	state_symbols <- accersi_state_symbols(n_states=partition_combos$n_states[pc]);
#	state_symbols <- (1:partition_combos$n_states[pc])-1;
#	state_symbols[state_symbols>=10] <- letter_states[state_symbols[state_symbols>=10]-9];
	scribere_rev_bayes_nexus_file_from_character_matrix(ch_matrix=chmatrix_red,state_symbols = state_symbols,new_file_name=new_file_name,UNKNOWN,INAP);
	matrix_file_names <- c(matrix_file_names,paste(data_file_lead,orig_file_name,sep=""));
	}

state_numbers <- partition_combos$n_states;
state_ordering <- partition_combos$type;
if (!is.null(partition_combos$rate_partitions))	{
	character_rate_partitions <- partition_combos$rate_partitions;
	} else	{
	character_rate_partitions <- rep("imagine",nrow(partition_combos));
	}
if (!is.null(partition_combos$trend_partitions))	{
	character_trend_partitions <- partition_combos$trend_partitions;
	} else	{
	character_trend_partitions <- rep("square",nrow(partition_combos));
	}

output <- list(initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,character_rate_partitions,character_trend_partitions,coding_bias);
names(output) <- c("initial_data","matrix_file_names","partition_size","state_numbers","state_ordering","rate_partitions","trend_partitions","coding_bias");
return(output);
}

#scribere_rev_bayes_script(analysis_name,taxon_names=otu_names,matrix_file_names,state_numbers,state_ordering,outgroup_taxa,unscored_taxa,fbd_parameterization_script,no_runs=4,write_scripts_directory=write_scripts_directory,set_wdir)
#scribere_Rev_Bayes_script_for_partitioned_character_data(            analysis_name,   initial_data,matrix_file_names,state_numbers,state_ordering,coding,           write_scripts_directory,fbd_parameterization_script,character_rate_partitions,fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data <- function(analysis_name="",initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",write_file=T)	{
# analysis_name: name that specificies this particular analysis.  (The clade name is a good choice)
# initial_data: data from nexus file.  
# matrix_file_names: a vector giving the list of all character matrices to be used
# partition_size: vector giving number of characters per partition
# state_numbers: vector giving the number of states for each matrix in matrix_file_names
# state_numbers: vector giving the number of states for each matrix in matrix_file_names
# state_ordering: vector designating unordered or ordered state evolution
# write_scripts_directory: tell R where to send script file
# FBD: if "true", then add script to initiate FBD analyses from a separate file
# set_wdir: tell RevBayes where to set working directory
# output_file_lead: tell RevBayes script where to send/find output (default: "output/")
# script_file_lead: tell RevBayes script where to send/find scripts (default: "script/")
# no_runs: tell RevBayes how many runs to execut (default=4)
#filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Partitioned_Analysis.Rev",sep="");
filename <- paste(write_scripts_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (length(unique(character_trend_partitions))>1)
	filename <- paste(filename,"_Driven_Trend_Test",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Stepping_Stone_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Stepping_Stone_Analysis.Rev",sep="");
	}
revbayes_source <- "clear();"
if (set_wdir=="")	{
	set_wdir <- getwd();
	set_wdir <- strsplit(getwd(),"/")[[1]];
	last_cell <- match("R_Projects",set_wdir);
	set_wdir[last_cell] <- "RevBayes_Projects";
	set_wdir <- paste(set_wdir[1:last_cell],collapse="/");
	}
if (set_wdir!="")
	revbayes_source <- c(revbayes_source,paste("setwd(\"",set_wdir,"\");\t#CHANGE THIS TO THE FOLDER IN WHICH YOU HAVE REVBAYES SCRIPTS & DATA!!!",sep=""));

revbayes_source <- c(revbayes_source,"### This director needs three subdirectories (folders):");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/scripts (additional RevBayes routines that will be used)");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/data (holds data matrices & taxonomic information)");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/output (where trees & logs will be sent)");

revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,"Imperio_Default_Settings.Rev",sep=""),"\");",sep=""));
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"###############################################################################");
revbayes_source <- c(revbayes_source,"# This is (these are) the nexus file(s) that you are using for this analysis  #");
revbayes_source <- c(revbayes_source,"#     Make sure that filenames & directories are correct!!!");                #
revbayes_source <- c(revbayes_source,"###############################################################################");
file_name_string <- paste(matrix_file_names, collapse = "\", \"");
file_name_string <- paste("filenames <- v(\"",file_name_string,"\");",sep="");
file_name_string <-gsub("~/","",file_name_string);
revbayes_source <- c(revbayes_source,file_name_string);
revbayes_source <- c(revbayes_source,paste("partition_chars <- v(",paste(partition_size, collapse = ","),");",sep=""));
revbayes_source <- c(revbayes_source,paste("partition_states <- v(",paste(state_numbers, collapse = ","),");",sep=""));
revbayes_source <- c(revbayes_source,paste("partition_ordering <- v(\"",paste(state_ordering, collapse = "\",\""),"\");",sep=""));
revbayes_source <- c(revbayes_source,paste("coding_bias <- v(\"",paste(coding_bias, collapse = "\",\""),"\");\t## prepare for ascertainment bias in binary characters; 'all': invariant & autapomorphies present; 'variable': all vary & autapomorphies present; 'informative': all vary & no autapomorphies.",sep=""));
if (length(unique(character_rate_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("rate_partitions <- v(\"",paste(character_rate_partitions, collapse = "\",\""),"\");\t# rate partition for each character paritition",sep=""));
	partition_labels <- unique(character_rate_partitions)
	revbayes_source <- c(revbayes_source,paste("rate_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_rate_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	}

if (length(unique(character_trend_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("driven_trend_partitions <- v(\"",paste(tolower(character_trend_partitions), collapse = "\",\""),"\");\t# use 'driven' for characters with biased change",sep=""));
	partition_labels <- unique(tolower(character_trend_partitions));
	revbayes_source <- c(revbayes_source,paste("trend_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_trend_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	}
revbayes_source <- c(revbayes_source,paste("max_age <- ",max_age,";\t\t\t\t\t\t# used if big_bang==TRUE;",sep=""));
revbayes_source <- c(revbayes_source,"");
taxon_names <- initial_data$OTUs;
notu <- length(initial_data$OTUs);
if(as.numeric(initial_data$Outgroup[1])!=-1)	{
	if (length(outgroup_taxa)>1)	{
		outies <- paste(outgroup_taxa,collapse="\",\"");
#			outies <- paste("v(\"",outies,sep="");
		outies <- paste("v(\"",outies,"\")",sep="");
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(",outies,");",sep=""));
		} else	{
		outies <- taxon_names[as.numeric(initial_data$Outgroup)];
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
		}
	hip_crowd <- "\"";
	hip_crowd <- paste(hip_crowd,paste(ingroup_taxa,collapse="\",\""),sep="");
	hip_crowd <- paste(hip_crowd,"\"",sep="");
	hip_crowd <- gsub(" ","_",hip_crowd);
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(",hip_crowd,");",sep=""));
	} else	{
	outies <- "ENTER_AN_OUTGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
	hipsters <- "ENTER_THE_INGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(\"",hipsters,"\");",sep=""));
	}
if (length(initial_data$Unscored_Taxa)>0)
	revbayes_source <- c(revbayes_source,paste("unscored_taxa <- v(",paste(initial_data$Unscored_Taxa,collapse=","),");",sep=""));

revbayes_source <- c(revbayes_source,"among_char_var <- \"lognormal\"\t\t# enter \"gamma\" or \"lognormal\"");
revbayes_source <- c(revbayes_source,"clock_model <- \"strict\";\t# enter \"strict\" for strict clock, or \"lognormal\" for relaxed clock with lognormal variation; we'll add \"dirichlet\ eventually;");

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"#                  Get basic information about the clade                   #");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"n_data_subsets <- filenames.size();");
if (fbd_parameterization_script!="")	{
#	revbayes_source <- c(revbayes_source,paste("intervals = readDataDelimitedFile(file=\"data/",fossil_interval_file,"\",header=true);",sep=""));
	xxx <- strsplit(fossil_interval_file,"/")[[1]];
	fbd_file_name <- xxx[length(xxx)];
	last_cell <- c(match("R_Projects",xxx),match("RevBayes_Projects",xxx))[!is.na(c(match("R_Projects",xxx),match("RevBayes_Projects",xxx)))];
	
	if (length(last_cell)>0 && !is.na(last_cell))	{
		xxx[last_cell+1] <- data_file_lead;
		directory_lead <- paste(xxx[1:(last_cell+1)],collapse="/")
		revbayes_source <- c(revbayes_source,paste("taxa <- readTaxonData(file=\"",data_file_lead,fbd_file_name,"\");",sep=""));
		revbayes_source <- c(revbayes_source,paste("n_taxa <- taxa.size();"));
		} else	{
		revbayes_source <- c(revbayes_source,paste("taxa <- readTaxonData(file=\"",data_file_lead,fbd_file_name,"\");",sep=""));
		revbayes_source <- c(revbayes_source,paste("n_taxa <- taxa.size();"));
		}
	} else {
	revbayes_source <- c(revbayes_source,"dummy <- readDiscreteCharacterData(filenames[1]);");
	revbayes_source <- c(revbayes_source,"taxa <- dummy.taxa();");
	revbayes_source <- c(revbayes_source,"n_taxa <- dummy.ntaxa();");
	}
revbayes_source <- c(revbayes_source,"n_branches <- (2 * n_taxa) - 2;");
revbayes_source <- c(revbayes_source,"");
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up appropriate parameters for speciation, extinction & sampling.     #");
	revbayes_source <- c(revbayes_source,"#      We also set up the tree search here.                                #");
	revbayes_source <- c(revbayes_source,"############################################################################");
#	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,fbd_parameterization_script),"\");",sep="");
	revbayes_source <- c(revbayes_source,"moves = VectorMoves();");
	revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,fbd_parameterization_script,sep=""),"\");",sep=""));
	} else	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up tree-search moves");
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"moves = VectorMoves();");
	revbayes_source <- c(revbayes_source,"topology ~ dnUniformTopology(taxa,outgroup);");
	revbayes_source <- c(revbayes_source,"moves.append = mvNNI(topology, weight=1.0);   # nearest neighbor interchange");
	revbayes_source <- c(revbayes_source,"moves.append = mvSPR(topology, weight=1.0);   # subtree pruning");
	revbayes_source <- c(revbayes_source,"for (b in 1:n_branches) {");
	revbayes_source <- c(revbayes_source,"    bl[b] ~ dnExponential(10.0);");
	revbayes_source <- c(revbayes_source,"    moves.append = mvScale(bl[b]);");
 	revbayes_source <- c(revbayes_source,"   }");
	revbayes_source <- c(revbayes_source,"tau := treeAssembly(topology, bl);  # assign branch lengths to trees");
	revbayes_source <- c(revbayes_source,"");
	}

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Set up appropriate Q-matrices for the partitions");
revbayes_source <- c(revbayes_source,"#   as well as the among-character and among-branch");
revbayes_source <- c(revbayes_source,"#   rate variation models");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the directory is OK)");
revbayes_source <- c(revbayes_source,"############################################################################");
#if (length(unique(character_rate_partitions))>1)	{
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Parameters_for_Analysis_Partitioned_by_States_and_Ordering_and_Class.Rev\");",sep=""));
#	} else	{
#	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Parameters_for_Analysis_Partitioned_by_States_and_Ordering.Rev\");",sep=""));
#	}
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Wrap it all into your model");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"mymodel = model(tau);\t\t# tau should have FBD & character evolution models attached to it");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Add monitors & stone your data");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the source directory is OK)");
revbayes_source <- c(revbayes_source,"# NOTE: the program saves trees once every printgen generations; so, the");
revbayes_source <- c(revbayes_source,"#   lower the number, the more trees you save.");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"monitors = VectorMonitors();");
revbayes_source <- c(revbayes_source,"");

revbayes_source <- c(revbayes_source,"cats=20;              # Number of stepping stones; each will reduce the \"weight\" of the likelihood by ~1/cats)");
revbayes_source <- c(revbayes_source,"burnin_gens=10000;    # Number of generations for the burnin pre-analysis (to tune parameters).");
revbayes_source <- c(revbayes_source,"tuning_int=200;       # Frequency at which burnin analysis will tune parameters (in generations).");
revbayes_source <- c(revbayes_source,"running_gens=100000;	# Number of generations for the real analysis; the bigger the analysis, the more you usually need.");
revbayes_source <- c(revbayes_source,"sampleFreq=1000;      # Frequency of stepping stone samples");

revbayes_source <- c(revbayes_source,"");
stone_output <- paste("output/",tolower(analysis_name),"_",sep="");
revbayes_source <- c(revbayes_source,paste("output_file = \"",stone_output,"\" + among_char_var + \"_char_variation_\" + clock_model;",sep=""));
revbayes_source <- c(revbayes_source,"if (clock_model!=\"strict\")	output_file = output_file + \"_relaxed\";");
revbayes_source <- c(revbayes_source,"output_file = output_file + \"_clock\";");
revbayes_source <- c(revbayes_source,"if (rate_partitions.size()>1)   output_file = output_file + \"_\" + \"_char_rate_partitions_\";");
revbayes_source <- c(revbayes_source,"if (bins>1)	output_file = output_file + \"_\" + bins + \"_bin_skyline\";");
revbayes_source <- c(revbayes_source,"filename1 = output_file + \"_stepping_stoned.log\";");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"monitors.append(mnModel(filename=filename1, printgen=10));");
revbayes_source <- c(revbayes_source,"monitors.append(mnFile(tau,filename=filename1, printgen=10,separator=TAB,tau));");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"filename2 = output_file + \"_stepping_stone_test.out\";");
revbayes_source <- c(revbayes_source,"update = \"Writing to: \" + filename2;");
revbayes_source <- c(revbayes_source,"print(update);");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"pow_p = powerPosterior(mymodel, moves, monitors, filename=filename2, cats=cats, sampleFreq=sampleFreq);  ##Set up your power posterior from everything in completed analysis. Create output for power posterior in quotes");
revbayes_source <- c(revbayes_source,"pow_p.burnin(generations=burnin_gens,tuningInterval=tuning_int);			##Set up power posterior burn in. Should likely be logner than 10000");
revbayes_source <- c(revbayes_source,"pow_p.run(generations=running_gens);");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"#######let run#################");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"ss = steppingStoneSampler(file=filename2, powerColumnName=\"power\", likelihoodColumnName=\"likelihood\");");
revbayes_source <- c(revbayes_source,"ps = pathSampler(file=filename2, powerColumnName=\"power\", likelihoodColumnName=\"likelihood\");");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"ss.marginal();   ##Calculate and display marginal likelihood of stepping stone simulations");
revbayes_source <- c(revbayes_source,"ps.marginal();   ##Calculate and display marginal likelihood of stepping stone simulations");
if (write_file)
	write(revbayes_source,file=filename);
return(revbayes_source);
}

scribere_MCMC_RevBayes_script_for_partitioned_character_data <- function(analysis_name="",initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",write_file=T,no_runs=4)	{
# analysis_name: name that specificies this particular analysis.  (The clade name is a good choice)
# initial_data: data from nexus file.  
# matrix_file_names: a vector giving the list of all character matrices to be used
# partition_size: vector giving number of characters per partition
# state_numbers: vector giving the number of states for each matrix in matrix_file_names
# state_numbers: vector giving the number of states for each matrix in matrix_file_names
# state_ordering: vector designating unordered or ordered state evolution
# write_scripts_directory: tell R where to send script file
# FBD: if "true", then add script to initiate FBD analyses from a separate file
# set_wdir: tell RevBayes where to set working directory
# output_file_lead: tell RevBayes script where to send/find output (default: "output/")
# script_file_lead: tell RevBayes script where to send/find scripts (default: "script/")
# no_runs: tell RevBayes how many runs to execut (default=4)
#filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Partitioned_Analysis.Rev",sep="");
filename <- paste(write_scripts_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (length(unique(character_trend_partitions))>1)
	filename <- paste(filename,"_Driven_Trend_Test",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_MCMC_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_MCMC_Analysis.Rev",sep="");
	}
revbayes_source <- "clear();"
if (set_wdir=="")	{
	set_wdir <- getwd();
	set_wdir <- strsplit(getwd(),"/")[[1]];
	last_cell <- match("R_Projects",set_wdir);
	set_wdir[last_cell] <- "RevBayes_Projects";
	set_wdir <- paste(set_wdir[1:last_cell],collapse="/");
	}
if (set_wdir!="")
	revbayes_source <- c(revbayes_source,paste("setwd(\"",set_wdir,"\");\t#CHANGE THIS TO THE FOLDER IN WHICH YOU HAVE REVBAYES SCRIPTS & DATA!!!",sep=""));

revbayes_source <- c(revbayes_source,"### This director needs three subdirectories (folders):");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/scripts (additional RevBayes routines that will be used)");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/data (holds data matrices & taxonomic information)");
revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/output (where trees & logs will be sent)");

revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,"Imperio_Default_Settings.Rev",sep=""),"\");",sep=""));
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"###############################################################################");
revbayes_source <- c(revbayes_source,"# This is (these are) the nexus file(s) that you are using for this analysis  #");
revbayes_source <- c(revbayes_source,"#     Make sure that filenames & directories are correct!!!");                #
revbayes_source <- c(revbayes_source,"###############################################################################");
file_name_string <- paste(matrix_file_names, collapse = "\", \"");
file_name_string <- paste("filenames <- v(\"",file_name_string,"\");",sep="");
file_name_string <-gsub("~/","",file_name_string);
revbayes_source <- c(revbayes_source,file_name_string);
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,paste("partition_chars <- v(",paste(partition_size, collapse = ","),");",sep=""));
revbayes_source <- c(revbayes_source,paste("partition_states <- v(",paste(state_numbers, collapse = ","),");",sep=""));
revbayes_source <- c(revbayes_source,paste("partition_ordering <- v(\"",paste(state_ordering, collapse = "\",\""),"\");",sep=""));
revbayes_source <- c(revbayes_source,paste("coding_bias <- v(\"",paste(coding_bias, collapse = "\",\""),"\");\t## prepare for ascertainment bias in binary characters; 'all': invariant & autapomorphies present; 'variable': all vary & autapomorphies present; 'informative': all vary & no autapomorphies.",sep=""));
if (length(unique(character_rate_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("rate_partitions <- v(\"",paste(character_rate_partitions, collapse = "\",\""),"\");\t# rate partition for each character paritition",sep=""));
	partition_labels <- unique(character_rate_partitions)
	revbayes_source <- c(revbayes_source,paste("rate_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_rate_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	}

if (length(unique(character_trend_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("driven_trend_partitions <- v(\"",paste(tolower(character_trend_partitions), collapse = "\",\""),"\");\t# use 'driven' for characters with biased change",sep=""));
	partition_labels <- unique(tolower(character_trend_partitions));
	revbayes_source <- c(revbayes_source,paste("trend_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_trend_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	}
revbayes_source <- c(revbayes_source,paste("max_age <- ",max_age,";\t\t\t\t\t\t# used if big_bang==TRUE;",sep=""));
revbayes_source <- c(revbayes_source,"");
taxon_names <- initial_data$OTUs;
notu <- length(initial_data$OTUs);
if(as.numeric(initial_data$Outgroup[1])!=-1)	{
	if (length(outgroup_taxa)>1)	{
		outies <- paste(outgroup_taxa,collapse="\",\"");
#			outies <- paste("v(\"",outies,sep="");
		outies <- paste("v(\"",outies,"\")",sep="");
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(",outies,");",sep=""));
		} else	{
		outies <- taxon_names[as.numeric(initial_data$Outgroup)];
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
		}
	hip_crowd <- "\"";
	hip_crowd <- paste(hip_crowd,paste(ingroup_taxa,collapse="\",\""),sep="");
	hip_crowd <- paste(hip_crowd,"\"",sep="");
	hip_crowd <- gsub(" ","_",hip_crowd);
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(",hip_crowd,");",sep=""));
	} else	{
	outies <- "ENTER_AN_OUTGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
	hipsters <- "ENTER_THE_INGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(\"",hipsters,"\");",sep=""));
	}
if (length(initial_data$Unscored_Taxa)>0)
	revbayes_source <- c(revbayes_source,paste("unscored_taxa <- v(",paste(initial_data$Unscored_Taxa,collapse=","),");",sep=""));

revbayes_source <- c(revbayes_source,"among_char_var <- \"lognormal\"\t\t## ENTER_THE_AMONG-CHARACTER_RATE_DISTRIBUTION_YOU_WISH_TO_USE_HERE\";\t# enter \"gamma\" or \"lognormal\"");
#revbayes_source <- c(revbayes_source,"clock_model <- \"ENTER_THE_CLOCK_MODEL_YOU_WISH_TO_USE_HERE\";\t# enter \"strict\" for strict clock, or \"lognormal\" for relaxed clock with lognormal variation; we'll add \"dirichlet\ eventually;");

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"#                  Get basic information about the clade                   #");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"n_data_subsets <- filenames.size();");
if (fbd_parameterization_script!="")	{
#	revbayes_source <- c(revbayes_source,paste("intervals = readDataDelimitedFile(file=\"data/",fossil_interval_file,"\",header=true);",sep=""));
	xxx <- strsplit(fossil_interval_file,"/")[[1]];
	fbd_file_name <- xxx[length(xxx)];
	last_cell <- match("R_Projects",xxx);
	xxx[last_cell+1] <- data_file_lead;
	directory_lead <- paste(xxx[1:(last_cell+1)],collapse="/")
	revbayes_source <- c(revbayes_source,paste("taxa <- readTaxonData(file=\"",data_file_lead,fbd_file_name,"\");",sep=""));
	revbayes_source <- c(revbayes_source,paste("n_taxa <- taxa.size();"));
	} else {
	revbayes_source <- c(revbayes_source,"dummy <- readDiscreteCharacterData(filenames[1]);");
	revbayes_source <- c(revbayes_source,"taxa <- dummy.taxa();");
	revbayes_source <- c(revbayes_source,"n_taxa <- dummy.ntaxa();");
	}
revbayes_source <- c(revbayes_source,"n_branches <- (2 * n_taxa) - 2;");
revbayes_source <- c(revbayes_source,"");
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up appropriate parameters for speciation, extinction & sampling.     #");
	revbayes_source <- c(revbayes_source,"#      We also set up the tree search here.                                #");
	revbayes_source <- c(revbayes_source,"#                                                                          #");
	revbayes_source <- c(revbayes_source,"# NOTE: This will sometimes freeze; if it does, then edit the script so    #");
	revbayes_source <- c(revbayes_source,"#      origination & extinction are set to 1.0. This usually works!        #");
	revbayes_source <- c(revbayes_source,"############################################################################");
#	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,fbd_parameterization_script),"\");",sep="");
	revbayes_source <- c(revbayes_source,"moves = VectorMoves();");
	revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,fbd_parameterization_script,sep=""),"\");",sep=""));
	} else	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up tree-search moves");
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"moves = VectorMoves();");
	revbayes_source <- c(revbayes_source,"topology ~ dnUniformTopology(taxa,outgroup);");
	revbayes_source <- c(revbayes_source,"moves.append = mvNNI(topology, weight=1.0);   # nearest neighbor interchange");
	revbayes_source <- c(revbayes_source,"moves.append = mvSPR(topology, weight=1.0);   # subtree pruning");
	revbayes_source <- c(revbayes_source,"for (b in 1:n_branches) {");
	revbayes_source <- c(revbayes_source,"    bl[b] ~ dnExponential(10.0);");
	revbayes_source <- c(revbayes_source,"    moves.append = mvScale(bl[b]);");
 	revbayes_source <- c(revbayes_source,"   }");
	revbayes_source <- c(revbayes_source,"tau := treeAssembly(topology, bl);  # assign branch lengths to trees");
	revbayes_source <- c(revbayes_source,"");
	}

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Set up appropriate Q-matrices for the partitions");
revbayes_source <- c(revbayes_source,"#   as well as the among-character and among-branch");
revbayes_source <- c(revbayes_source,"#   rate variation models");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the directory is OK)");
revbayes_source <- c(revbayes_source,"############################################################################");
#if (length(unique(character_rate_partitions))>1)	{
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Parameters_for_Analysis_Partitioned_by_States_and_Ordering_and_Class.Rev\");",sep=""));
#	} else	{
#	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Parameters_for_Analysis_Partitioned_by_States_and_Ordering.Rev\");",sep=""));
#	}
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Wrap it all into your model");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"mymodel = model(tau);\t\t# tau should have FBD & character evolution models attached to it");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Add monitors & commence MCMC'ing");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the source directory is OK)");
revbayes_source <- c(revbayes_source,"# NOTE: the program saves trees once every printgen generations; so, the");
revbayes_source <- c(revbayes_source,"#   lower the number, the more trees you save.");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"monitors = VectorMonitors();");
tree_file_name <- analysis_name;
if (length(unique(character_rate_partitions))>1)
	tree_file_name <- paste(tree_file_name,paste(partition_labels,collapse="_vs_"),sep="_");
mcmc_output <- paste("output/",tree_file_name,sep="");
revbayes_source <- c(revbayes_source,"if (clock_model==\"strict\")\t{");
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnModel(filename=\"",mcmc_output,"_Strict_Clock.log\", printgen=10));",sep=""));
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnFile(tau, filename=\"",mcmc_output,"_Strict_Clock.trees\",printgen=10,separator=TAB,tau));",sep=""));
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,alpha,fbd_p,fbd_q,fbd_r,num_samp_anc,origin_time));");
	} else	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,alpha,origin_time));");
	}
#analyses_outputs <- paste(output_file_lead,tree_file_name,"_Strict_Clock",sep="");
revbayes_source <- c(revbayes_source,"\t} else if (clock_model==\"big_bang\" || clock_model==\"early_burst\") \t{");
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnModel(filename=\"",mcmc_output,"_Early_Burst.log\", printgen=10));",sep=""));
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnFile(tau, filename=\"",mcmc_output,"_Early_Burst.trees\",printgen=10,separator=TAB,tau));",sep=""));
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,rel_bang,alpha,fbd_p,fbd_q,fbd_r,num_samp_anc,origin_time));");
	} else	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,rel_bang,alpha,origin_time));");
	}
#analyses_outputs <- paste(output_file_lead,tree_file_name,"_Early_Burst",sep="");
revbayes_source <- c(revbayes_source,"\t} else if (clock_model==\"uncorrelated\")\t{");
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnModel(filename=\"",mcmc_output,"_Uncorrelated_Relaxed_Clock.log\", printgen=10));",sep=""));
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnFile(tau, filename=\"",mcmc_output,"_Uncorrelated_Relaxed_Clock.trees\",printgen=10,separator=TAB,tau));",sep=""));
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,ucln_var,alpha,fbd_p,fbd_q,fbd_r,num_samp_anc,origin_time));");
	} else	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,ucln_var,alpha,origin_time));");
	}
#analyses_outputs <- paste(output_file_lead,tree_file_name,"_Uncorrelated_Relaxed_Clock",sep="");
revbayes_source <- c(revbayes_source,"\t} else if (clock_model==\"autocorrelated\")\t{");
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnModel(filename=\"",mcmc_output,"_Autocorrelated_Relaxed_Clock.log\", printgen=10));",sep=""));
revbayes_source <- c(revbayes_source,paste("\tmonitors.append(mnFile(tau, filename=\"",mcmc_output,"_Autocorrelated_Relaxed_Clock.trees\",printgen=10,separator=TAB,tau));",sep=""));
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,acln_var,alpha,fbd_p,fbd_q,fbd_r,num_samp_anc,origin_time));");
	} else	{
	revbayes_source <- c(revbayes_source,"\tmonitors.append(mnScreen(printgen=500,mean_rt,acln_var,alpha,origin_time));");
	}
#analyses_outputs <- paste(output_file_lead,tree_file_name,"_Autcorrelated_Relaxed_Clock",sep="");
revbayes_source <- c(revbayes_source,"\t}");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"    ################################################################################");
revbayes_source <- c(revbayes_source,"    # Here are some starting parameters for your MCMC analysis: but use your own!");
revbayes_source <- c(revbayes_source,"    # NOTE: as the number of moves increases, the greater the number of generations");
revbayes_source <- c(revbayes_source,"    #     we need to make a thorough search of parameter space.  So, as taxa and ");
revbayes_source <- c(revbayes_source,"    #     and complexity of character evolution models increases, the greater the ");
revbayes_source <- c(revbayes_source,"    #     number of generations you should use. ");
revbayes_source <- c(revbayes_source,"    ################################################################################");
revbayes_source <- c(revbayes_source,paste("no_runs=",no_runs,";\t\t# Number of independent MCMC analyses. (Even MCMC can get stuck in local optima!)",sep=""));
revbayes_source <- c(revbayes_source,"burnin_gens=10000;\t# Number of generations for the burnin pre-analysis (to tune parameters).");
revbayes_source <- c(revbayes_source,"tuning_int=200;\t\t# Frequency at which burnin analysis will tune parameters (in generations).");
revbayes_source <- c(revbayes_source,"running_gens=1000000;\t# Number of generations for the real analysis; the bigger the analysis, the more you usually need.");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"# Now, go read Anna Karenina.....");
revbayes_source <- c(revbayes_source,paste("source(\"scripts/Expecto_MCMC_with_Partitioned_Characters.Rev\");",sep=""));
revbayes_source <- c(revbayes_source,"# .......");
revbayes_source <- c(revbayes_source,"# Sigh: nobody remembers elementary train safety anymore.  Oh, your trees are done.");
revbayes_source <- c(revbayes_source,"");

revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Prepare MCMC output to get consensus tree(s) and the most probable tree(s)");
revbayes_source <- c(revbayes_source,"#    As always, double check the directories....");
revbayes_source <- c(revbayes_source,"############################################################################");

analyses_outputs <- paste(output_file_lead,tolower(tree_file_name),sep="");
if (no_runs>1)	{
	trees <- maj_rule <- max_prob_tree <- c();
	for (nr in 1:no_runs)	{
		trees <- c(trees,paste("\"",analyses_outputs,"_run_",nr,".trees\"",sep=""));
		maj_rule <- c(maj_rule,paste("\"",analyses_outputs,"_run_",nr,"_maj_rule.tre\"",sep=""));
		max_prob_tree <- c(max_prob_tree,paste("\"",analyses_outputs,"_run_",nr,"_simple_map.tre\"",sep=""));
		}
	trees <- paste(trees,collapse=",");
	maj_rule <- paste(maj_rule,collapse=",");
	max_prob_tree <- paste(max_prob_tree,collapse=",");
	revbayes_source <- c(revbayes_source,paste("tree_files <- v(",trees,");",sep=""));
	revbayes_source <- c(revbayes_source,paste("maj_rule_files <- v(",maj_rule,");",sep=""));
	revbayes_source <- c(revbayes_source,paste("most_probable_files <- v(",max_prob_tree,");",sep=""));
	} else	{
	trees <- paste("\"",analyses_outputs,".trees\"",sep="");
	maj_rule <- paste("\"",analyses_outputs,"_maj_rule.tre\"",sep="");
	max_prob_tree <- paste("\"",analyses_outputs,"_simple_map.tre\"",sep="");
	revbayes_source <- c(revbayes_source,paste("tree_files <- ",trees,";",sep=""));
	revbayes_source <- c(revbayes_source,paste("maj_rule_files <- ",maj_rule,";",sep=""));
	revbayes_source <- c(revbayes_source,paste("most_probable_files <- ",max_prob_tree,";",sep=""));
	}
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Consensus_Tree.Rev\");",sep=""));
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accersi_Most_Probable_Tree.Rev\");",sep=""));
if (write_file)
	write(revbayes_source,file=filename);
return(revbayes_source);
}

# write out script to setup FBD analysis using origination, extinction, preservation and divergence bound dates calculated by some other routine.
scribere_fbd_portion_of_Rev_Bayes_script <- function(analysis_name,write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa="",script_file_lead="scripts/")	{
fbd_script <- c();
fbd_script <- c(fbd_script,"########################################################################");
fbd_script <- c(fbd_script,"# Set up appropriate parameters for speciation, extinction & sampling  #");
fbd_script <- c(fbd_script,"#   \"Seed\" numbers based on analyses of Paleobiology Database data.    #");
fbd_script <- c(fbd_script,"########################################################################");
fbd_script <- c(fbd_script,paste("# Diversification Rates based on ",control_taxon,sep=""));
fbd_script <- c(fbd_script,paste("origination_rate ~ dnExponential(",round(origination,3),");",sep=""));
fbd_script <- c(fbd_script,"moves.append(mvScale(origination_rate, lambda=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvScale(origination_rate, lambda=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvScale(origination_rate, lambda=1.00, weight=1));");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #");
fbd_script <- c(fbd_script,"# NOTE: FBD scripts often allow extinction to vary independently of speciation;     #");
fbd_script <- c(fbd_script,"# However, empirical studies show that these two rates usually are close to equal   #");
fbd_script <- c(fbd_script,"#               and they definitely are not independent.                            #");
fbd_script <- c(fbd_script,"# So, here we'll make turnover (ext/orig) an independent variable and use it        #");
fbd_script <- c(fbd_script,"#               to scale extinction relative to origination                         #");
fbd_script <- c(fbd_script,"# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #");
fbd_script <- c(fbd_script,"turnover ~ dnUnif(0.9, 1.05);");
fbd_script <- c(fbd_script,"moves.append(mvSlide(turnover, delta=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(turnover, delta=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(turnover, delta=1.00, weight=1));");
fbd_script <- c(fbd_script,"extinction_rate := turnover*origination_rate;");

fbd_script <- c(fbd_script,"");
if (sampling_unit=="rock")
	sampling_unit <- "rock units"
fbd_script <- c(fbd_script,paste("# Fossil Sampling Rates based on ",sampling_unit," occupied by ",control_taxon,sep=""));
fbd_script <- c(fbd_script,paste("psi ~ dnExponential(",round(1/psi,3),");",sep=""));
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=1.00, weight=1));");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"# Proportional Taxon Sampling of Youngest Time Slice");
fbd_script <- c(fbd_script,paste("rho <- ",round(rho,3),";	# 'extant' sampling.",sep=""));
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"# Establish Basal Divergence Time");
fbd_script <- c(fbd_script,paste("origin_time ~ dnUnif(",max(divergence_bounds),", ",min(divergence_bounds),");",sep=""));
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=1.00, weight=1));");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"fbd_dist = dnFBDRP(origin=origin_time, lambda=origination_rate, mu=extinction_rate, psi=psi, rho=rho, taxa=taxa);");

fbd_script <- c(fbd_script,"############################################################################");
fbd_script <- c(fbd_script,"#                               Set up tree                                #");
fbd_script <- c(fbd_script,"############################################################################");
fbd_script <- c(fbd_script,"# create the vector of clade constraints");
fbd_script <- c(fbd_script,"constraints = v(ingroup);");
fbd_script <- c(fbd_script,"tau ~ dnConstrainedTopology(fbd_dist,constraints=constraints);");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"moves.append(mvFNPR(tau, weight=n_branches/2));                              # time-tree pruning & grafting");
fbd_script <- c(fbd_script,"moves.append(mvNNI(tau, weight=n_branches/2));                               # nearest-neighbor interchanges");
fbd_script <- c(fbd_script,"moves.append(mvCollapseExpandFossilBranch(tau,origin_time,weight=n_taxa/4)); # consider ancestor-descendant rather than sister species");
fbd_script <- c(fbd_script,"moves.append(mvNodeTimeSlideUniform(tau, weight=n_branches/2));              # adjust divergence times");
fbd_script <- c(fbd_script,"moves.append(mvRootTimeSlideUniform(tau, origin_time, weight=5));            # adjust basal divergence time.");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"num_samp_anc := tau.numSampledAncestors();");
fbd_script <- c(fbd_script,"for (bn in 1:n_branches)\t{");
fbd_script <- c(fbd_script,"\tdivergence_dates[bn]:=tau.nodeAge(bn)                   # this is when a hypothesized ancestor diverges or an OTU is first seen;");
fbd_script <- c(fbd_script,"\tbranch_lengths[bn]:=tau.branchLength(bn);               # this is branch *duration* not expected change!");
fbd_script <- c(fbd_script,"\torigin_dates[bn]:=tau.branchLength(bn)+tau.nodeAge(bn); # this is when a lineage diverged from its ancestor");
fbd_script <- c(fbd_script,"\t}");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"#### Set up deterministic variables for output purposes only  ####");
fbd_script <- c(fbd_script,"fbd_p:=origination_rate;\t\t# origination rate for output");
fbd_script <- c(fbd_script,"fbd_q:=origination_rate;\t\t# extinction rate for output");
fbd_script <- c(fbd_script,"fbd_r:=psi;\t\t# sampling rate for output");
fbd_script <- c(fbd_script,"completeness := psi/(extinction_rate+psi);");
fbd_script <- c(fbd_script,"diversification := origination_rate - extinction_rate;");
fbd_script <- c(fbd_script,"summed_gaps := sum(branch_lengths);");
fbd_script <- c(fbd_script,"");

zombies <- "\"";
zombies <- paste(zombies,paste(extant_taxa,collapse="\",\""),sep="");
zombies <- paste(zombies,"\"",sep="");
zombies <- gsub(" ","_",zombies);
fbd_script <- c(fbd_script,paste("clade_extant = clade(",zombies,");",sep=""));
#fbd_script <- c(fbd_script,"age_extant := tmrca(tau, clade_extant);\t# There is no particularly good reason to keep this!");
fbd_script <- c(fbd_script,"");
if (length(uncoded_taxa)>0 && uncoded_taxa!="")	{
	fbd_script <- c(fbd_script,paste("pruned_tau := fnPruneTree(tau, prune=v(\"",paste(gsub(" ","_",otu_names[uncoded_taxa]),collapse="\",\""),"\"));",sep=""));
#	} else	{
#	fbd_script <- c(fbd_script,paste(fbd_script,"pruned_tau := fnPruneTree(tau, prune=v(\",\"))\t\t#All taxa coded!",sep=""));
	}
output_file <- paste(write_scripts_directory,"Accersi_",analysis_name,"_Range_Based_FBD_Parameterization.Rev",sep="");
write(fbd_script,output_file);
output <- list(fbd_script,paste("Accersi_",analysis_name,"_Range_Based_FBD_Parameterization.Rev",sep=""));
names(output) <- c("script","filename");
#return(paste(study,"_Range_Based_FBD_Parameterization.Rev",sep=""));
return(output);
}

# script that lists "recent" (= latest contemporaneous) taxa
list_faux_extant_taxa <- function(analysis_name,write_scripts_directory,fossil_intervals)	{
extant_intervals <- subset(fossil_intervals,fossil_intervals$min==min(abs(fossil_intervals$min)));
taxon_names_for_file <- gsub(" ","_",extant_intervals$taxon);
clade_extant <- c();
clade_extant <- c(clade_extant,"###################################################################");
clade_extant <- c(clade_extant,"#    Read in the \"Recent\" taxa (i.e., latest co-extant taxa)    #");
clade_extant <- c(clade_extant,"###################################################################");
clade_extant <- c(clade_extant,"");
clade_extant_info <- "clade_extant = clade(";
for (tt in 1:nrow(extant_intervals))	{
	clade_extant_info <- paste(clade_extant_info,"\"",taxon_names_for_file[tt],sep="");
	if (tt < nrow(extant_intervals))
		clade_extant_info <- paste(clade_extant_info,"\",",sep="");
	}
clade_extant_info <- paste(clade_extant_info,"\");",sep="");
clade_extant <- c(clade_extant,clade_extant_info);
output_file <- paste(write_scripts_directory,analysis_name,"_Read_Faux_Extant.Rev",sep="");
write(clade_extant,file=output_file);
return(paste(analysis_name,"_Read_Faux_Extant.Rev",sep=""));
}

#fossil_intervals <- read.table(file=fossil_interval_file,header=T);
#fossil_intervals <- read.table(file=file.choose(),header=T);
faux_extant_taxa <- function(fossil_intervals)	{
return(gsub(" ","_",as.character(fossil_intervals$taxon[fossil_intervals$min==min(abs(fossil_intervals$min))])));
#extant_taxa <- gsub(" ","_",extant_intervals$taxon);
#return(extant_taxa);
}

#scribere_RevBayes_scripts_from_nexus_file_and_PaleoDB_download <- function(analysis_name,nexus_file_name,taxon_subset_file="",rate_partition="",trend_partition="",taxon_level,lump_subgenera=F,species_only=T,bogarted="",rock_unit_databases="",chronostratigraphic_databases="",paleodb_fixes="",control_taxon="",zone_taxa="",onset="Cambrian",end="Holocene",end_FBD="",exclude_uncertain_taxa=T,basic_environments=c("terr","marine","unknown"),sampling_unit="collections",time_scale_stratigraphic_scale="International",temporal_precision=0.1,read_data_directory="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22)	{
scribere_RevBayes_scripts_from_nexus_file_and_PaleoDB_download <- function(analysis_name,taxon_subset_file=F,rate_partition="",trend_partition="",taxon_level="species",lump_subgenera=F,species_only=T,bogarted="",rock_unit_databases="",chronostratigraphic_databases="",paleodb_fixes="",control_taxon="",zone_taxa="",onset="Cambrian",end="Holocene",end_FBD="",exclude_uncertain_taxa=T,basic_environments=c("terr","marine","unknown"),sampling_unit="collections",time_scale_stratigraphic_scale="International",temporal_precision=0.1,write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",polymorphs=T,UNKNOWN=-11,INAP=-22)	{
#### PART 0: Commence ####
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses,");
print("    using data downloaded from the Paleobiology Database (https://www.paleobiodb.org/) for stratigraphic data and then");
print("    for refining/cleaning/updating those data with updated time scales and biozonation information.");
print("");
print("NOTE: The Paleobiology Database should always be considered a STARTING point for these data. This program will also");
print("\toutput summaries of the data that you should check. We encourage you to contribute improvements to these data");
print("\t(to collections, occurrences and/or taxonomy) to the Paleobiology Database.  Improvements to the stratigraphic");
print("\tdatabase used to refine PaleoDB data should be sent to pjwagner@gmail.com");
print("");
print("The program will partition the matrix based on numbers of states and whether states are ordered or unordered. The output");
print("   RevBayes scripts create appropriate Q-matrices for each partition. The program also outputs stratigraphic information");
print("   that RevBayes uses in Fossilized Birth-Death analyses (a type of birth-death-sampling analyses). It provides initial");
print("   estimates of sampling and diversification rates that are used to seed the FBD analyses, although these are varied by");
print("   RevBayes over the MCMC searches. It also provides a \"recent\" sampling rate that reflects the most-likely per-taxon");
print("   sampling rate for the latest interval relevant to the study.  (This pretends that the latest taxa are 'recent'). ");
print("");
if (taxon_level=="genus" && !lump_subgenera)
	taxon_level <- "subgenus";
#### PART 1: GET CHARACTER DATA & BASIC TAXON INFORMATION ####
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	for (i in 1:100)	j <- 1;
	} else	{   
	print("Choose the nexus file that you wish to analyze: ");
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = F,stringsAsFactors=hell_no)[,1];
	verboten <- c("taxon","taxa","species","genus","otu");
	taxa_subset <- taxa_subset[!taxa_subset %in% verboten];
	print("Choose the nexus file that you wish to analyze: ");
	}

basic_data <- divido_character_matrix_by_state_numbers_and_other_partitions(analysis_name=analysis_name,write_data_directory=write_data_directory,rate_partition=rate_partition,trend_partition=trend_partition,taxa_subset=taxa_subset,data_file_lead="data/",polymorphs=polymorphs,UNKNOWN=UNKNOWN,INAP=INAP);
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!
if (species_only)	{
	taxon_names <- ingroup_taxa;
	clade_members <- unique(sapply(taxon_names,divido_genus_names_from_species_names));
	} else	{
	clade_members <- ingroup_taxa;
	}

#### PART 2: LOAD EXTERNAL DATA FOR CLEANING & REFINING PALEODB DATA  ####
fossilworks_collections <- paleodb_fixes$fossilworks_collections;
paleodb_rock_reidentifications <- paleodb_fixes$paleodb_rock_reidentifications;
paleodb_collection_edits <- paleodb_fixes$paleodb_collection_edits;
if (!is.null(paleodb_collection_edits$X))
	paleodb_collection_edits$X <- NULL;
time_scale <- chronostratigraphic_databases$time_scale;
zone_database <- chronostratigraphic_databases$zones;
if (is.list(rock_unit_databases))	{
	rock_database <- rock_unit_databases$rock_unit_database;
	rock_to_zone_database <- rock_unit_databases$rock_to_zone_database;
	rock_to_zone_database$ma_lb <- temporal_precision*round(rock_to_zone_database$ma_lb/temporal_precision,0);
	rock_to_zone_database$ma_ub <- temporal_precision*round(rock_to_zone_database$ma_ub/temporal_precision,0);
	}
time_scale$ma_lb <- temporal_precision*round(time_scale$ma_lb/temporal_precision,0);
time_scale$ma_ub <- temporal_precision*round(time_scale$ma_ub/temporal_precision,0);
zone_database$ma_lb <- temporal_precision*round(zone_database$ma_lb/temporal_precision,0);
zone_database$ma_ub <- temporal_precision*round(zone_database$ma_ub/temporal_precision,0);

zone_database <- subset(zone_database,zone_database$ma_lb<=time_scale$ma_lb[match(onset,time_scale$interval)]+5);
zone_database <- subset(zone_database,zone_database$ma_ub>=time_scale$ma_ub[match(end,time_scale$interval)]-5);

#### PART 3: GET INFORMATION NEEDED TO DOWNLOAD, 'CLEAN' AND ANALYZE STRATIGRAPHIC DATA  ####
compendium <- accersi_updated_taxonomy_for_analyzed_taxa(otu_names=otu_names_used,local_directory=local_directory,study=analysis_name);
if (bogarted)	{
	print("Choose the file with your private stash information: ");
	for (i in 1:100)
		j <- i;
	}
if (bogarted)	{
	bogarted_info <- file.choose();
	print("Reading your private stash now....");
	bogarted_finds <- read.csv(bogarted_info,header = T,stringsAsFactors=hell_no);
	bogarted_finds <- evanesco_na_from_matrix(bogarted_finds,replacement="");
	bogarted_finds <- subset(bogarted_finds,bogarted_finds$identified_name!="");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_name <- bogarted_finds$identified_name;
		bogarted_finds$genus <- as.character(sapply(taxon_name,divido_genus_names_from_species_names));
		if (taxon_level=="subgenus")	{
			genus_name <- bogarted_finds$genus;
			subgenus_results <- sapply(genus_name,divido_subgenus_names_from_genus_names);
#			bogarted_finds$genus <- subgenus_results[1,];
			bogarted_finds$subgenus <- subgenus_results[2,];
			bogarted_finds$subgenus[bogarted_finds$subgenus==""] <- bogarted_finds$genus[bogarted_finds$subgenus==""];
			}
	#	add occurrences
		unique_genera <- unique(bogarted_finds$genus);
		if (!is.null(bogarted_finds$subgenus))
			unique_genera <- unique(c(bogarted_finds$genus,bogarted_finds$subgenus));
		for (u_g in 1:length(unique_genera))	{
			if (!is.na(match(unique_genera[u_g],compendium$taxon_name)))	{
				compendium$n_occs[match(unique_genera[u_g],compendium$taxon_name)] <- length(unique(bogarted_finds$collection_no[unique(which(bogarted_finds==unique_genera[u_g],arr.ind = T)[,1])]));
				}
			}
		
		}
	if (!is.null(bogarted_finds$direct_ma))	{
		bogarted_finds$direct_ma <- as.numeric(bogarted_finds$direct_ma);
		bogarted_finds$direct_ma[is.na(bogarted_finds$direct_ma)] <- 0;
		}
	if (!is.null(bogarted_finds$direct_ma_error))	{
		bogarted_finds$direct_ma_error <- as.numeric(bogarted_finds$direct_ma_error);
		bogarted_finds$direct_ma_error[is.na(bogarted_finds$direct_ma_error)] <- 0;
		}
	bogarted_finds$max_ma <- time_scale$ma_lb[match(bogarted_finds$early_interval,time_scale$interval)];
	bogarted_finds$late_interval[bogarted_finds$late_interval==""] <- bogarted_finds$early_interval[bogarted_finds$late_interval==""];
	bogarted_finds$min_ma <- time_scale$ma_ub[match(bogarted_finds$late_interval,time_scale$interval)];
	bogarted_taxa <- unique(bogarted_finds$identified_name);
	notu <- nrow(compendium);
	compendia_that_are_new <- match(bogarted_finds$identified_name,compendium$taxon_name);
	compendia_that_are_new <- compendia_that_are_new[!is.na(compendia_that_are_new)];
	if (length(compendia_that_are_new)>0)
		compendium$n_occs <- compendium$n_occs+hist(compendia_that_are_new,breaks=(0:notu),plot=F)$counts;
	}

if (sum(compendium$n_occs==0)>0)	{
	missing_taxa <- subset(compendium,compendium$n_occs==0);
	missing_taxa <- subset(missing_taxa,missing_taxa$accepted_name=="?");
	taxon_name <- missing_taxa$taxon_name;
	missing_taxa_rows <- match(taxon_name,compendium$taxon_name);
	taxon_list <- genera <- unique(sapply(taxon_name,divido_genus_names_from_species_names));
	if (length(taxon_list)>0)	{
		taxonomyless_finds <- accersi_occurrences_for_list_of_taxa(taxon_list,paleogeography=paleogeography);
		if (is.data.frame(taxonomyless_finds$collection_compendium))
			for (mt in 1:length(missing_taxa))
				compendium$n_occs[missing_taxa_rows[mt]] <- sum(taxonomyless_finds$occurrences_compendium$identified_name==taxon_name[mt]);
		}
	### insert command for getting occurrences & collections for lists of taxa here.
	}

if (sum(compendium$n_occs==0)>0)	{
#	print(paste("The following taxa have no occurrences:",paste(compendium$taxon_name[compendium$n_occs==0],collapse=", ")));
	print("The following taxa currently have no occurrences entered into the PaleoDB:");
	print(compendium$taxon_name[compendium$n_occs==0]);
#	print(paste("The following taxa are not entered into the PaleoDB:",paste(compendium$taxon_name[compendium$taxon_no==""],collapse=", ")));
	print("");
	if (sum(compendium$taxon_no=="")>0)	{
		print("The following taxa are not entered into the PaleoDB taxonomy tables:");
		print(compendium$taxon_name[compendium$taxon_no==0]);
		}
	print("Enter Data for these into the PaleoDB and try again tomorrow or setup a separate 'bogarted' file with occurrences for these taxa!");
	print("   Make sure the file as formation, member, stage, zonation, etc., information, too. (And consider entering it into the PaleoDB later.)");
	print("Also, make sure that there are no misspellings in your nexus matrix. (Computers do not autocorrect!)");
	return();
	}
otu_names[compendium$accepted_name!="?"] <- compendium$accepted_name[compendium$accepted_name!="?"];
## get paleodb data!!!!
paleodb_data <- accersi_paleodb_data_for_Rev_Bayes(otu_names,analysis_name=analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa,taxon_level,onset,end,basic_environments,time_scale,zone_database,fossilworks_collections,paleodb_rock_reidentifications,paleodb_collection_edits,lump_subgenera,species_only);
control_collections <- paleodb_data$control_collections;
control_occurrences <- paleodb_data$control_occurrences;
control_collections$collection_no <- as.numeric(control_collections$collection_no);
control_occurrences$collection_no <- as.numeric(control_occurrences$collection_no);
control_occurrences$occurrence_no <- as.numeric(control_occurrences$occurrence_no);
if (bogarted)	{
	print("Adding your private stash to the PaleoDB data....");
	bogarted_finds$collection_no <- as.numeric(bogarted_finds$collection_no);
	bogarted_finds$paleodb_collection_no[bogarted_finds$paleodb_collection_no==""] <- 0;
	bogarted_finds$paleodb_collection_no <- as.numeric(bogarted_finds$paleodb_collection_no);
	bogarted_finds$paleodb_collection_no[bogarted_finds$paleodb_collection_no==0] <- bogarted_finds$collection_no[bogarted_finds$paleodb_collection_no==0]+ceiling(max(control_collections$collection_no)/10^(floor(log10(max(control_collections$collection_no)))-1))*10^(floor(log10(max(control_collections$collection_no)))-1);

	colnames(bogarted_finds)[match("collection_no",colnames(bogarted_finds))] <- "my_collection_no"
	colnames(bogarted_finds)[match("paleodb_collection_no",colnames(bogarted_finds))] <- "collection_no"
	column_matches <- match(colnames(bogarted_finds),colnames(control_collections))
	bogarted_coll_info_in_paleodb <- (1:ncol(bogarted_finds))[!is.na(column_matches)];
	column_matches <- column_matches[!is.na(column_matches)];
	new_paleodb_coll <- control_collections[1:length(unique(bogarted_finds$collection_no)),];
	for (nc in 1:ncol(new_paleodb_coll))	{
		if (is.numeric(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- 0;
			} else if (is.character(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- "";
			}
		}
#	new_paleodb_coll <- control_collections[length(unique(bogarted_finds$collection_no)),];
	new_paleodb_coll[,column_matches] <- bogarted_finds[match(unique(bogarted_finds$collection_no),bogarted_finds$collection_no),bogarted_coll_info_in_paleodb];
	control_collections <- rbind(control_collections,new_paleodb_coll);
	
	# set up occcurrences in two steps;
	# first get collection part of occurrences
	column_matches <- match(colnames(bogarted_finds),colnames(control_occurrences))
	bogarted_occr_info_in_paleodb <- (1:ncol(bogarted_finds))[!is.na(column_matches)];
	column_matches <- column_matches[!is.na(column_matches)];
	new_paleodb_occr <- control_occurrences[1:nrow(bogarted_finds),];
	for (nc in 1:ncol(new_paleodb_occr))	{
		if (is.numeric(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.numeric(0);
			} else if (is.character(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.character("");
			new_paleodb_occr[,nc] <- as.character(new_paleodb_occr[,nc]);
			}
		}
	new_paleodb_occr[,column_matches] <- bogarted_finds[,bogarted_occr_info_in_paleodb];
	
	# now get the taxonomy part....
	taxon <- bogarted_taxa;
	for (tt in 1:length(bogarted_taxa))	{
		if (tt==1)	{
			bogarted_taxonomy <- revelare_taxonomy_for_one_taxon(taxon=bogarted_taxa[tt],settle=T);
			} else	{
			bogarted_taxonomy <- rbind(bogarted_taxonomy,revelare_taxonomy_for_one_taxon(taxon=bogarted_taxa[tt],settle=T));
			}
		}
	bogarted_taxonomy$accepted_name[bogarted_taxonomy$taxon_name!=bogarted_taxa] <- bogarted_taxonomy$taxon_name[bogarted_taxonomy$taxon_name!=bogarted_taxa] <- bogarted_taxa[bogarted_taxonomy$taxon_name!=bogarted_taxa];
	bogarted_taxonomy$accepted_rank[match(bogarted_taxonomy$accepted_rank,taxonomic_rank)>match(taxon_level,taxonomic_rank)] <- taxon_level;
	bogarted_taxonomy <- evanesco_na_from_matrix(bogarted_taxonomy,replacement = "");
	bogarted_taxonomy$record_type <- bogarted_taxonomy$flags <- bogarted_taxonomy$early_interval <- bogarted_taxonomy$late_interval <- NULL;
	bogarted_row_to_paledob <- match(bogarted_finds$identified_name,bogarted_taxonomy$taxon_name);
	paleodb_col_to_edit <- match(colnames(bogarted_taxonomy),colnames(control_occurrences));
	paleodb_col_to_edit <- paleodb_col_to_edit[!is.na(paleodb_col_to_edit)];
	bogarted_col_w_fix <- match(colnames(control_occurrences)[paleodb_col_to_edit],colnames(bogarted_taxonomy));
	new_paleodb_occr[,paleodb_col_to_edit] <- bogarted_taxonomy[bogarted_row_to_paledob,bogarted_col_w_fix];
	new_paleodb_occr$record_type <- control_occurrences$record_type[1];
	
	#(1:nrow(control_occurrences))[is.na(as.numeric(control_occurrences$occurrence_no))]
	new_paleodb_occr$occurrence_no <- (1:nrow(new_paleodb_occr))+10^ceiling(log10(max(control_occurrences$occurrence_no)))
	control_occurrences <- rbind(control_occurrences,new_paleodb_occr);
	}

if (taxon_level=="genus" || taxon_level=="subgenus")	
	control_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = control_occurrences);

if (is.data.frame(paleodb_data$zone_occurrences))	{
	zone_occurrences <- paleodb_data$zone_occurrences;
	if (taxon_level=="genus" || taxon_level=="subgenus")
		zone_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = zone_occurrences);
	}
hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(control_collections$early_interval)),unique(as.character(control_collections$late_interval)))),time_scale,regional_scale=time_scale_stratigraphic_scale);
hierarchical_chronostrat$ma_lb <- temporal_precision*round(hierarchical_chronostrat$ma_lb/temporal_precision,0);
hierarchical_chronostrat$ma_ub <- temporal_precision*round(hierarchical_chronostrat$ma_ub/temporal_precision,0);
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
if (is.na(match(end_FBD,finest_chronostrat$interval)))
	end_FBD <- finest_chronostrat$interval[max(which(finest_chronostrat==end_FBD,arr.ind = T)[,1])];
downloaded_collections <- reset_paleodb_intervals_to_desired_time_scale(collections=control_collections,finest_chronostrat = finest_chronostrat,time_scale);
downloaded_collections$min_ma <- round(temporal_precision*round(downloaded_collections$min_ma/temporal_precision,0),floor(-log10(temporal_precision)));
downloaded_collections$max_ma <- round(temporal_precision*round(downloaded_collections$max_ma/temporal_precision,0),floor(-log10(temporal_precision)));
ncolls <- nrow(downloaded_collections);

#### PART 4: REFINE CHRONOSTRATIGRAPHY OF PALEODB DATA  ####
if (!is.null(zone_occurrences))	{
	paleodb_finds <- rbind(control_occurrences,zone_occurrences);
	paleodb_finds <- paleodb_finds[order(paleodb_finds$collection_no,paleodb_finds$occurrence_no),];
	paleodb_finds <- paleodb_finds[match(unique(paleodb_finds$occurrence_no),paleodb_finds$occurrence_no),];
	} else	{
	paleodb_finds <- control_occurrences;
	}
if (is.list(rock_unit_databases))	{
	print("Refining PaleoDB data with rock-unit and biozonation databases...");
	paleodb_data_refined <- refine_collection_dates_with_external_database(study=analysis_name,collections=downloaded_collections,rock_database,zone_database,rock_to_zone_database,time_scale,directory=local_directory);
	refined_collections <- paleodb_data_refined$Recalibrated_Collections;
	chronostrat_units <- unique(c(hierarchical_chronostrat$interval,refined_collections$interval_lb,refined_collections$interval_ub));
	hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units,time_scale,regional_scale=time_scale_stratigraphic_scale);
	finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
	} else {
	if (is.data.frame(zone_database))	{
		print("Refining PaleoDB data with biozonation databases...");
		refined_collections <- refine_paleodb_collection_dates_with_zone_data_only(paleodb_collections=downloaded_collections,paleodb_finds,zone_database,time_scale,hierarchical_chronostrat,finest_chronostrat,examine_finds=T,temporal_precision=0.05);
		} else	{
		refined_collections <- downloaded_collections;
		refined_collections$ma_lb <- refined_collections$max_ma;
		refined_collections$ma_ub <- refined_collections$min_ma;
		refined_collections$interval_lb <- refined_collections$early_interval;
		refined_collections$interval_ub <- refined_collections$late_interval;
		}
	}

age <- temporal_precision*round(refined_collections$ma_lb/temporal_precision,0);
refined_collections$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
age <- temporal_precision*round(refined_collections$ma_ub/temporal_precision,0);
refined_collections$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);

# use quantitative biostratigraphy 101 to refine dates here.
print("Using basic biostratigraphy to minimize gaps for uncertainly aged collections...");
optimized_collections <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds,paleodb_collections=refined_collections,hierarchical_chronostrat,zone_database,update_search=T);
#ddd <- (1:ncolls)[optimized_collections$ma_lb<=optimized_collections$ma_ub]
# use radiometric dates here
if (!is.null(optimized_collections$direct_ma) && sum(optimized_collections$direct_ma>0)>0)	{
	print("Using radiometric data for final ages...");
	optimized_collections <- redate_collections_with_direct_dates(collections=optimized_collections,finest_chronostrat);
	}

# get rock unit numbers if we do not have a stratigraphic database
if (is.null(optimized_collections$rock_no))	{
	print("Putting numbers on unique rock units...");
	optimized_collections <- number_unique_rock_units(paleodb_collections = optimized_collections,zone_database,time_scale);
	}

# for unentered rock units that are unique to their time and location, create dummy numbers.
optimized_collections <- name_unnamed_rock_units(paleodb_collections=optimized_collections,finest_chronostrat);

paleodb_finds <- control_occurrences;
ncolls <- nrow(optimized_collections);
finest_chronostrat$ma_lb <- temporal_precision*round(finest_chronostrat$ma_lb/temporal_precision,0)
finest_chronostrat$ma_ub <- temporal_precision*round(finest_chronostrat$ma_ub/temporal_precision,0)
paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections=optimized_collections,uniform_time_scale = finest_chronostrat);
print(paste("Saving",paste(analysis_name,"_Refined_Collections.csv",sep=""),"..."));
write.csv(paleodb_collections,paste(local_directory,analysis_name,"_Refined_Collections.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

#### PART 5: GET STRATIGRAPHIC DATA THAT REVBAYES CAN USE  ####
# at this point, it becomes a little easier if we have chronostratigraphic data directly tied to finds
#paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections = paleodb_collections,uniform_time_scale = finest_chronostrat);
if (is.null(paleodb_finds$ma_lb))	{
	paleodb_finds$ma_lb <- paleodb_collections$ma_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	paleodb_finds$ma_ub <- paleodb_collections$ma_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	}
if (is.null(paleodb_finds$bin_lb))	{
	if (is.null(paleodb_collections$bin_lb))	{
		paleodb_collections$bin_lb <- match(paleodb_collections$interval_lb,finest_chronostrat$interval);
		paleodb_collections$bin_ub <- match(paleodb_collections$interval_ub,finest_chronostrat$interval);
		}
	paleodb_finds$bin_lb <- paleodb_collections$bin_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	paleodb_finds$bin_ub <- paleodb_collections$bin_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	}
print(paste("Saving",paste(analysis_name,"_Plus_Control_Finds_recalibrated.csv",sep=""),"..."));
write.csv(paleodb_finds,paste(local_directory,analysis_name,"_Plus_Control_Finds_recalibrated.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

## output for only the ingroup
if (taxon_level=="species")	{
	otu_finds <- paleodb_finds[paleodb_finds$accepted_name %in% otu_names_used,];
	} else if (lump_subgenera==T)	{
	otu_finds <- paleodb_finds[paleodb_finds$genus %in% otu_names_used,];
	} else	{
	taxon_names <- otu_names_used;
	genus_names <- sapply(taxon_names,divido_genus_names_from_species_names);
	xxx <- sapply(genus_names,divido_subgenus_names_from_genus_names)
	subgenera_used <- otu_names_used;
	subgenera_used[xxx[2,]!=""] <- xxx[2,xxx[2,]!=""];
	otu_finds <- paleodb_finds[paleodb_finds$subgenus %in% subgenera_used,];
	find_order <- match(otu_finds$subgenus,subgenera_used)
	otu_finds <- otu_finds[order(find_order,-otu_finds$ma_lb),];
	}
print(paste("Saving",paste(analysis_name,"_Finds_recalibrated.csv",sep=""),"..."));
write.csv(otu_finds,paste(local_directory,analysis_name,"_Finds_recalibrated.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

strat_for_Rev_Bayes <- accersi_stratigraphic_information_for_Rev_Bayes(taxa=as.character(otu_names),paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD = end_FBD,temporal_precision=temporal_precision);
end_FBD_b <- rebin_collection_with_time_scale(age=min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_fa),onset_or_end = "end",fine_time_scale = finest_chronostrat);
if (end_FBD!="" && hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)] < match(end_FBD_b,finest_chronostrat$interval))	{
	end_FBD <- end_FBD_b;
	strat_for_Rev_Bayes <- accersi_stratigraphic_information_for_Rev_Bayes(taxa=as.character(otu_names),paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD,temporal_precision);
	}
per_bin_info_ingroup <- accersi_per_stratigraphic_interval_sampling_information_for_Rev_Bayes(taxa=otu_names,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD = end_FBD,temporal_precision=temporal_precision);

# output data for RevBayes & for you to examine 
fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals.tsv",sep="");
fossil_interval_file_FA <- paste(tolower(analysis_name),"_fossil_intervals_FA.tsv",sep="");
print(paste("Saving",fossil_interval_file,"for RevBayes to use and other files for you to examine...."));
fossil_intervals <- strat_for_Rev_Bayes$fossil_intervals;
fossil_intervals$taxon <- gsub(" ","_",otu_names_used);	# make sure that this file uses the same names as the original analysis!
if (taxon_subset_file)	{
	keeper_taxa <- gsub(" ","_",taxa_subset);	# make sure that this file uses the same names as the original analysis!
	keeper_rows <- match(taxa_subset,fossil_intervals$taxon);
	} else	{
	keeper_rows <- 1:nrow(fossil_intervals);
	}
write.table(fossil_intervals[keeper_rows,],file=paste(write_data_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals[keeper_rows,],file=paste(local_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);

if (min(fossil_intervals$min)!=0)	{
	new_min <- min(fossil_intervals$min);
	fossil_intervals$min <- new_min-fossil_intervals$min;
	fossil_intervals$max <- new_min-fossil_intervals$max;
	}
fossil_intervals_FA <- fossil_intervals;
fossil_intervals_FA$min <- fossil_intervals_FA$max-min(fossil_intervals_FA$max);
fossil_intervals_FA$max <- fossil_intervals_FA$max-min(fossil_intervals_FA$max);
print(paste("Saving",fossil_interval_file_FA,"for RevBayes to use and other files for you to examine...."));
write.table(fossil_intervals_FA[keeper_rows,],file=paste(write_data_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_FA[keeper_rows,],file=paste(local_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);

fossil_intervals_fuzzy <- strat_for_Rev_Bayes$fossil_intervals_fuzzy;
fossil_intervals_fuzzy$taxon <- gsub(" ","_",otu_names_used);	# make sure that this file uses the same names as the original analysis!
fuzzy_fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals_fuzzy.tsv",sep="");
print(paste("Saving",fuzzy_fossil_interval_file,"for RevBayes to use and other files for you to examine...."));
write.table(fossil_intervals_fuzzy[keeper_rows,],file=paste(write_data_directory,fuzzy_fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_fuzzy[keeper_rows,],file=paste(local_directory,fuzzy_fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);

detailed_information_file <- paste(tolower(analysis_name),"_detailed_fossil_information.csv",sep="");
detailed_information <- strat_for_Rev_Bayes$fossil_information_detailed;
write.csv(detailed_information[keeper_rows,],file=paste(local_directory,detailed_information_file,sep=""),row.names = F);

write.csv(per_bin_info_ingroup$finds_per_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Finds_per_Bin.csv",sep=""),row.names = T);
write.csv(per_bin_info_ingroup$definite_finds_per_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Definite_Finds_per_Bin.csv",sep=""),row.names = T);
write.csv(per_bin_info_ingroup$sampled_in_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Sampled_in_Bin.csv",sep=""),row.names = T);

# we use the "accepted_name" field for analysis; so, if this is not a species-level analysis, rewrite the field with the appropriate genus or subgenus
if (taxon_level!="species")
	paleodb_finds$accepted_name <- as.character(paleodb_finds[,match(taxon_level,colnames(paleodb_finds))]);

#### PART 6: GET INITIAL SAMPLING ESTIMATES ####
print("Now getting initial estimates of sampling rates per million years.....");
if (taxon_level=="species") {
	taxon_list2 <- sort(unique(paleodb_finds$accepted_name));
	} else	{
	taxon_list2 <- sort(unique(paleodb_finds[,match(taxon_level,colnames(paleodb_finds))]));
	}
end_FBD_z <- finest_chronostrat$interval[nrow(finest_chronostrat)];
per_bin_info_total <- accersi_per_stratigraphic_interval_sampling_information_for_Rev_Bayes(taxa=taxon_list2,paleodb_finds=paleodb_finds,paleodb_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit=sampling_unit,lump_cooccr=T,constrain=T,end_FBD=end_FBD_z,temporal_precision=temporal_precision);
sampled_in_bin_richness <- colSums(per_bin_info_total$finds_per_bin>=0.5);

taxon_ranges <- sepkoskify_paleodb_data(pbdb_finds=paleodb_finds,taxon_names=taxon_list2,interval_names=finest_chronostrat$interval);
taxon_lives <- 1+taxon_ranges$bin_ub-taxon_ranges$bin_lb;
names(taxon_lives) <- rownames(taxon_ranges);
interval_richness <- accersi_synoptic_richness(taxon_ranges=cbind(taxon_ranges$bin_lb,taxon_ranges$bin_ub),interval_names=finest_chronostrat$interval);
if (sampling_unit=="collection" || sampling_unit=="collections")	{
	sample_units_per_bin <- ceiling(tally_collections_occupied_by_subinterval(taxon_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision));
	} else	{
	sample_units_per_bin <- ceiling(tally_rock_units_occupied_by_subinterval(taxon_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision));
	}
if (length(sample_units_per_bin) != length(interval_richness))	{
	sample_units_per_bin <- sample_units_per_bin[names(sample_units_per_bin) %in% names(interval_richness)];
	interval_richness <- interval_richness[names(interval_richness) %in% names(sample_units_per_bin)];
	}
# get the most likely uniform, exponential, beta & lognormal distributions for per-collection or per-rock sampling
print("Estimating best uniform, exponential, beta & lognormal sampling distributions for each interval.....");
fpb <- per_bin_info_total$finds_per_bin;
fpb <- fpb[,(1:ncol(fpb))[colnames(fpb) %in% names(interval_richness)]];
if (ncol(fpb)!=length(sample_units_per_bin))	{
	# put in fix here!
	}

sampling_distributions <- accersi_sampling_distributions_for_RevBayes(finds_per_bin=fpb,sample_units_per_bin,end_FBD="");
# get the expected finds per bin given the distributions found above.
bin_spans <- (finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);
names(bin_spans) <- finest_chronostrat$interval;
bin_spans <- bin_spans[names(bin_spans) %in% names(interval_richness)];
psi_bin <- per_interval_psis(sampling_distributions,sample_units_per_bin);
psi_bin_pma <- per_interval_per_ma_psis(sampling_distributions,sample_units_per_bin,bin_spans);
#psi_bin <- psi_bin*(finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);

# USE THE MEDIAN OVERALL SAMPLING RATE AS THE STARTING SAMPLING RATE PER MILLION YEARS
#psi <- sum(psi_bin*bin_spans)/sum(bin_spans);

# THIS IS OUR SAMPLING PROBABILITY FOR THE LATEST SET ("RECENT") TAXA
if (end_FBD=="")	{
	youngest <- min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_la);
	end_FBD <- rebin_collection_with_time_scale(age=youngest,onset_or_end = onset,fine_time_scale = finest_chronostrat);
	}
if (is.na(match(end_FBD,names(psi_bin))))	{
	end_FBD <- finest_chronostrat$interval[max(which(finest_chronostrat==end_FBD,arr.ind = T)[,1])];
	if (is.na(match(end_FBD,names(psi_bin))))	{
		end_FBD <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb >= min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_fa)),];
		}
	}
last_bin <- hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)];
#psi <- median(psi_bin[1:last_bin]*bin_spans[1:last_bin]);
psi <- sum(psi_bin[1:match(end_FBD,names(psi_bin))])/sum(bin_spans[1:match(end_FBD,names(psi_bin))]);	# total median expected finds divided by total time
print(paste("The median per-ma sampling rate (psi) is: ",round(psi,4),".",sep=""));
faux_recent_bin <- hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)];
if (is.na(psi_bin[faux_recent_bin]))
	faux_recent_bin <- match(end_FBD,names(psi_bin));
rho <- Poisson_rate_to_probability((psi_bin[faux_recent_bin]*bin_spans[faux_recent_bin]));	# per-taxon sampling probability given parameters for the last interval
#fpb_frb <- sort(ceiling(fpb[fpb[,match(end_FBD,colnames(fpb))]>0,match(end_FBD,colnames(fpb))]),decreasing = T);
#rho2 <- chao2(abundance=fpb_frb);
print(paste("The ML per-taxon sampling rate for the final interval (rho) is: ",round(rho,4),".",sep=""));

#### PART 7: GET INITIAL ORIGINATION & EXTINCTION ESTIMATES ####
print("Estimating origination & extinction (given sampling) rates for each interval.....")
sampled_in_bin <- 1*per_bin_info_total$finds_per_bin>0.5;
for (sb in 1:ncol(sampled_in_bin))	{
	sampled_in_bin[sampled_in_bin[,sb]==T,sb] <- 1;
	sampled_in_bin[sampled_in_bin[,sb]==F,sb] <- 0;
	}
if (ncol(sampled_in_bin) != length(interval_richness))
	sampled_in_bin <- sampled_in_bin[,(1:ncol(sampled_in_bin))[colnames(sampled_in_bin) %in% names(interval_richness)]];
synoptic_richness <- interval_richness;
diversification <- accersi_initial_diversification_rates(sampled_in_bin,synoptic_richness=interval_richness,psi_bin=psi_bin,chronostrat=finest_chronostrat);
origination <- diversification[1];
extinction <- diversification[2];
#if (extinction>origination)	extinction <- 0.99*origination;
print(paste("Median ML estimates for origination is: ",round(origination,4)," per myr and extinction is: ",round(extinction,4)," per myr.",sep=""));

# GET INITIAL BOUNDS FOR DIVERGENCE TIMES;
print("Get basic estimate of initial divergence time using Bapst's cal-3 method.....")
phi <- prob_sampling_clade_bapst(p=origination,q=extinction,r=psi);
bound_1 <- max(fossil_intervals_FA$max[keeper_rows]);
bound_2 <- sort(fossil_intervals_FA$max[keeper_rows],decreasing=T)[2];
#initial_divergence <- (simple_probability_divergence(bound_1,bound_2,phi,psi) + simple_likelihood_divergence(bound_1,bound_2,psi))/2;
initial_divergence <- simple_probability_divergence(bound_1,bound_2,phi,psi);
divergence_bounds <- c(bound_1,initial_divergence);
print(paste("Initial divergence bounds are ",round(divergence_bounds[1],2)," to ",round(divergence_bounds[2],2)," Ma before end of study.",sep=""));

#### PART 8: START WRITING SCRIPTS ####
extant_file <- list_faux_extant_taxa(analysis_name,write_scripts_directory,fossil_intervals=strat_for_Rev_Bayes$fossil_intervals[keeper_rows,]);
#extant_taxa <- faux_extant_taxa(fossil_intervals=strat_for_Rev_Bayes$fossil_intervals[keeper_rows,]);
extant_taxa <- faux_extant_taxa(fossil_intervals=fossil_intervals_FA[keeper_rows,]);
#																	 study,              write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_file,otu_names,uncoded_taxa="",script_file_lead="script/"
fbd_parameterization <- scribere_fbd_portion_of_Rev_Bayes_script(analysis_name,write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa=initial_data$Unscored_Taxa,script_file_lead="scripts/");
fbd_parameterization_script <- paste("scripts/",fbd_parameterization$filename,sep="");

max_age <- max(fossil_intervals_FA$max[keeper_rows]);
#scribere_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="script/",no_runs=4);
#                  scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name="",initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",write_file=T)
revbayes_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,character_numbers,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file=fossil_interval_file_FA,set_wdir,output_file_lead="output/",script_file_lead="scripts/");#,no_runs=3);

write(fbd_parameterization$script,file=paste(paste(local_directory,"Accersi_",analysis_name,sep=""),"_Range_Based_FBD_Parameterization.Rev",sep=""));
filename <- paste(local_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_babble,file=filename);
#revbayes_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

scribere_RevBayes_scripts_from_chosen_nexus_file_and_existing_FBD_script_and_data <- function(analysis_name,fa_info=NULL,taxon_subset_file=F,rate_partition="",trend_partition="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",data_file_lead="data/",UNKNOWN=-11,INAP=-22,coding_threshold=1)	{
if (ncol(fa_info)>2)	{
	poss_cols <- unique(which(fa_info==max(fa_info[,2:ncol(fa_info)]),arr.ind = T)[,2]);
	if (length(poss_cols)>1)	{
		col_ages <- colSums(fa_info[,poss_cols]);
		poss_cols <- poss_cols[match(max(col_ages),col_ages)];
		}
	fa_info <- fa_info[,c(1,poss_cols)];
	colnames(fa_info) <- c("taxon","fa");
	}
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses.");
print("   The program will prompt you for (in order!):");
print("      1. The original nexus file;");
print("      2. A .tsv file giving first and last appearance dates of each taxon (in Ma before the youngest taxa);");
print("      3. A RevBayes script setting up the parameters for diversification & sampling;");
print("");
flush.console();
Sys.sleep(zzzz);
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	flush.console();
	Sys.sleep(zzzz);
#	for (i in 1:100)	j <- 1;
	} else	{   
	taxa_subset <- "";
	}
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = T,stringsAsFactors=hell_no)[,1];
	}

#basic_data <- divido_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
basic_data <- divido_character_matrix_by_state_numbers_and_other_partitions(analysis_name,first_appearances=fa_info,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead=data_file_lead,polymorphs=T,UNKNOWN,INAP,coding_threshold = coding_threshold);
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
partition_size <- basic_data$partition_size;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!
print("Select .tsv file with first and last appearance dates for the FBD analysis:");
flush.console();
Sys.sleep(zzzz);
clade_members <- ingroup_taxa;

fossil_interval_file <- file.choose();
print(paste("Using '",fossil_interval_file,"'",sep=""));
flush.console();
fossil_intervals <- read.table(fossil_interval_file,header=T,sep="\t",stringsAsFactors = F);
max_age <- max(fossil_intervals$max);

print("Choose the FBD parameterization script: ");
flush.console();
Sys.sleep(zzzz);
fbd_parameterization_script <- file.choose();
print(paste("Using '",fbd_parameterization_script,"'",sep=""));
flush.console();
Sys.sleep(zzzz);
break_it_down <- strsplit(fbd_parameterization_script,"")[[1]];
if (sum(break_it_down=="/")>0)	{
	script_file_lead <- "scripts/";
	output_file_lead <- "output/";
	pathway <- strsplit(fbd_parameterization_script,"/")[[1]];
	fbd_parameterization_script <- pathway[length(pathway)];
	} else	{
	script_file_lead <- "scripts\\";
	output_file_lead <- "output\\";
	pathway <- strsplit(fbd_parameterization_script,"\\\\")[[1]];
	if (length(pathway)==1)
		pathway <- strsplit(fbd_parameterization_script,"\\")[[1]];
	fbd_parameterization_script <- pathway[length(pathway)];
	}

#                        scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name="",initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",write_file=T)
revbayes_stone_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file,set_wdir,output_file_lead=output_file_lead,script_file_lead=script_file_lead,data_file_lead=data_file_lead,write_file=F);
revbayes_mcmc_babble <- scribere_MCMC_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file,set_wdir,output_file_lead=output_file_lead,script_file_lead=script_file_lead,data_file_lead=data_file_lead,write_file=F,no_runs=3);

filename_ss <- paste(local_directory,analysis_name,"_Stepping_Stone",sep="");
filename_mcmc <- paste(local_directory,analysis_name,"_MCMC",sep="");
if (length(unique(character_rate_partitions))>1)	{
	filename_ss <- paste(filename_ss,"_Rate_Partitioned",sep="");
	filename_mcmc <- paste(filename_mcmc,"_Rate_Partitioned",sep="");
	}
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename_ss <- paste(filename_ss,"_Analysis.Rev",sep="");
	filename_mcmc <- paste(filename_mcmc,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename_ss <- paste(filename_ss,"_FBD_Analysis.Rev",sep="");
	filename_mcmc <- paste(filename_mcmc,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_stone_babble,file=filename_ss);
write(revbayes_mcmc_babble,file=filename_mcmc);
#write(revbayes_babble,file=paste(paste(local_directory,analysis_name,sep=""),"_Partitioned_FBD_Analysis.Rev",sep=""));
#revbayes_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

scribere_RevBayes_scripts_from_nexus_file_and_existing_FBD_script_and_data <- function(analysis_name,fbd_parameterization_script,fossil_interval_file,taxon_subset_file=F,rate_partition="",trend_partition="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",polymorphs=T,UNKNOWN=-11,INAP=-22)	{
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses.");
print("   If conducting FBD analyses, then it relies on the user to provide the name of an FBD parameterization script as.");
print("   well as a file giving fossil intervals. IF you do not have these yet, then you should use another routine:" );
print("\t\tscribere_RevBayes_scripts_from_nexus_file_and_PaleoDB_download()");
print("");
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	for (i in 1:100)	j <- 1;
	} else	{   
#	print("Choose the nexus file that you wish to analyze: ");
#	print("Choose a nexus file to convert to RevBayes preferred format:")
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = T,stringsAsFactors=hell_no)[,1];
#	print("Choose the nexus file that you wish to analyze: ");
#	print("Choose a nexus file to convert to RevBayes preferred format:")
	}

#basic_data <- divido_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
#divido_character_matrix_by_state_numbers_and_other_partitions <- function(analysis_name="",first_appearances=NULL,write_data_directory="",rate_partition="",trend_partition="",taxa_subset="",data_file_lead="data/",polymorphs=T, UNKNOWN=-11, INAP=-22,coding_threshold=1)
basic_data <- divido_character_matrix_by_state_numbers_and_other_partitions(analysis_name,first_appearances=fossil_interval_file,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
print("Select the .tsv file with first and last appearance dates for the FBD analysis:");
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
partition_size <- basic_data$partition_size;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!

#if (species_only)	{
#	taxon_names <- ingroup_taxa;
#	clade_members <- unique(sapply(taxon_names,divido_genus_names_from_species_names));
#	} else	{
clade_members <- ingroup_taxa;
#	}

fossil_intervals <- read.table(file.choose(),header=T);
max_age <- max(fossil_intervals$max);

#                  scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name="",initial_data,matrix_file_names,partition_size,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",write_file=T)
revbayes_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data=initial_data,matrix_file_names=matrix_file_names,partition_size = partition_size,state_numbers=state_numbers,state_ordering=state_ordering,coding_bias=coding_bias,outgroup_taxa=outgroup_taxa,ingroup_taxa=ingroup_taxa,max_age=max_age,write_scripts_directory=write_scripts_directory,fbd_parameterization_script=fbd_parameterization_script,character_rate_partitions=character_rate_partitions,character_trend_partitions=character_trend_partitions,fossil_interval_file=fossil_interval_file,set_wdir=set_wdir,output_file_lead="output/",script_file_lead="scripts/");#,no_runs=4);

filename <- paste(local_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_babble,file=filename);
#write(revbayes_babble,file=paste(paste(local_directory,analysis_name,sep=""),"_Partitioned_FBD_Analysis.Rev",sep=""));
#revbayes_babble <- scribere_Stepping_Stone_RevBayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

# taxon_name <- "'Gyrocystis' platessa";
nexusify_taxon_name <- function(taxon_name)	{
molecularized_name <- strsplit(taxon_name,split="")[[1]];
if (sum(molecularized_name %in% newick_verbotten)==0)	{
	return(gsub(" ","_",taxon_name));
	} else if (molecularized_name[1]!="\'")	{
	return(paste(c("'",molecularized_name,"'"),collapse=""));
	} else	{
	return(taxon_name);
	}
}

accersi_state_symbols <- function(n_states)	{
state_symbols <- (1:n_states)-1;
state_symbols[state_symbols>=10] <- letter_states[state_symbols[state_symbols>=10]-9];
return(state_symbols);
}

#### 
##### ROUTINES TO DOWNLOAD, CLEAN & ORGANIZE PALEODB DATA #######

# get basic stratigraphic information that RevBayes demands
# updated 2020-05-02
accersi_paleodb_data_for_Rev_Bayes <- function(otu_names,analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa=T,taxon_level,onset,end,basic_environments=c("marine","unknown","terrestrial"),paleogeography="scotese",time_scale,zone_database,fossilworks_collections="",paleodb_rock_reidentifications="",paleodb_collection_edits="",lump_subgenera=F,species_only=T)	{
# otu_names: vector giving taxon names (matching those of the original nexus file)	
# study: name of the study (e.g., "Phacophida" or "Ordovician_Bucanids")	
# control_taxon: name of a taxonomic group that can be used as a control for sampling and diversification estimates
# zone_taxa: taxa such as conodonts, forams, trilobites, graptolites, ammonites, etc., that might co-occur with the study clade and that are used for biozonation;
# exclude_uncertain_taxa: if true (default), then questionable assignments are excluded
# onset: oldest collections & occurrences to download (should be older than the oldest members of you clade)
# end: youngest collections & occurrences to download
# basic_environment: environment type to download (defaults to all)
# fossilworks_collections: data.frame downloaded from Fossilworks.org
# paleodb_rock_reidentifications: data.frame providing corrections to the rock assignments of unedittable PaleoDB collections
# paleodb_collection_edits: data.frame providing other corrections to unedittable PaleoDB collections.

print("Getting occurrence & collection data for study taxa....")
data_compendium <- accersi_occurrences_for_list_of_taxa(taxon_list=otu_names,lump_subgenera,species_only,paleogeography=paleogeography);
data_compendium$occurrences_compendium$flags <- simplify2array(data_compendium$occurrences_compendium$flags);
ingroup_finds <- evanesco_na_from_matrix(data=data_compendium$occurrences_compendium,replacement = "");
#ingroup_finds[match(991528,ingroup_finds$occurrence_no),]
if (species_only)	{
	no_species <- c();
	for (ot in 1:length(otu_names))	{
		if((sum(ingroup_finds$identified_rank[unique(which(ingroup_finds==otu_names[ot],arr.ind = T)[,1])]=="species")+sum(ingroup_finds$identified_rank[unique(which(ingroup_finds==otu_names[ot],arr.ind = T)[,1])]=="subspecies"))==0)	{
			no_species <- c(no_species,otu_names[ot]);
			}
		}
	if (length(no_species)>0)	{
		print("Attention: you requested species-level occurrences only, but the following taxa have only genus-level occurrences:");
		print(paste("    ",paste(no_species,collapse=","),sep=""));
		}
	}

if (exclude_uncertain_taxa)	{
	ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain genus, uncertain species");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain genus");
		} else if (taxon_level=="species" || taxon_level=="subspecies")	{
		ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain species");
		}
	}

ingroup_finds$taxon <- as.character(ingroup_finds$taxon);
nfinds <- nrow(ingroup_finds);

ingroup_collections <- data_compendium$collection_compendium[order(data_compendium$collection_compendium$collection_no),];
ingroup_collections <- ingroup_collections[match(unique(ingroup_finds$collection_no),ingroup_collections$collection_no),];
ingroup_collections <- ingroup_collections[order(ingroup_collections$collection_no),];
ncolls <- nrow(ingroup_collections);

# sometimes it returns collections that are too young or too old because of the requested taxa; delete those collections
ingroup_collections <- subset(ingroup_collections,ingroup_collections$max_ma>=time_scale$ma_ub[match(end,time_scale$interval)]);
ingroup_collections <- subset(ingroup_collections,ingroup_collections$min_ma<=time_scale$ma_lb[match(onset,time_scale$interval)]);
ncolls <- nrow(ingroup_collections);
ingroup_finds <- ingroup_finds[(1:nrow(ingroup_finds))[ingroup_finds$collection_no %in% ingroup_collections$collection_no],]

# redate any collections that might be older or younger than the study interval: if we want only Jurassic collections, then Jurassic will be the oldest poss. early interval & youngest possible late interval
ingroup_collections$late_interval[ingroup_collections$min_ma<time_scale$ma_ub[match(end,time_scale$interval)]] <- end;
ingroup_collections$min_ma[ingroup_collections$min_ma<time_scale$ma_ub[match(end,time_scale$interval)]] <- time_scale$ma_ub[match(end,time_scale$interval)];
ingroup_collections$early_interval[ingroup_collections$max_ma>time_scale$ma_lb[match(onset,time_scale$interval)]] <- onset;
ingroup_collections$max_ma[ingroup_collections$min_ma<time_scale$ma_lb[match(onset,time_scale$interval)]] <- time_scale$ma_lb[match(onset,time_scale$interval)];

write.csv(ingroup_collections[order(ingroup_collections$collection_no),],file=paste(local_directory,analysis_name,"_Collections.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(ingroup_finds,file=paste(local_directory,analysis_name,"_Finds.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

if (control_taxon[1]!="")	{
	print("Getting occurrence & collection data for control taxa....");
	control_data <- accersi_data_for_control_groups_to_seed_FBD_analyses(control_taxon,onset,end,basic_environments,species_only = T,paleogeography=paleogeography);
	control_collections <- evanesco_na_from_matrix(control_data$control_collections,"");
	control_occurrences <- evanesco_na_from_matrix(control_data$control_occurrences,"");
	control_occurrences <- expello_indeterminate_species(paleodb_finds = control_occurrences);
	
	if (exclude_uncertain_taxa)	{
		control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain genus, uncertain species");
		if (taxon_level=="genus" || taxon_level=="subgenus")	{
			control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain genus");
			} else if (taxon_level=="species" || taxon_level=="subspecies")	{
			control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain species");
			}
		}
	lost_finds <- ingroup_finds[(1:nrow(ingroup_finds))[!ingroup_finds$occurrence_no %in% control_occurrences$occurrence_no],];
	if (nrow(lost_finds)>0)	{
		lost_finds$taxon <- NULL;
		control_occurrences <- rbind(control_occurrences,lost_finds);
		control_occurrences <- control_occurrences[order(control_occurrences$collection_no,control_occurrences$occurrence_no),];
		}
	control_collections <- control_collections[control_collections$collection_no %in% control_occurrences$collection_no,];

 	# WTF is any of this???
	mislaid_collections <- (1:ncolls)[!ingroup_collections$collection_no %in% control_collections$collection_no];
	if (length(mislaid_collections) > 0 )	{
		control_collections <- rbind(control_collections,ingroup_collections[mislaid_collections,]);
		control_collections <- control_collections[order(control_collections$collection_no),];
		ingroup_finds_temp <- ingroup_finds;
		ingroup_finds_temp$taxon <- NULL;
		control_occurrences <- rbind(control_occurrences,ingroup_finds_temp[ingroup_finds_temp$collection_no %in% ingroup_collections$collection_no[mislaid_collections],]);
		control_occurrences <- control_occurrences[order(control_occurrences$collection_no,control_occurrences$occurrence_no),];
		control_occurrences <- unique(control_occurrences);
		}
	if (length(control_taxon)>0)	{
		control_file_name <- paste(control_taxon,collapse="_&_");
		} else	{
		control_file_name <- control_taxon;
		}
	write.csv(evanesco_na_from_matrix(control_data$control_collections,""),file=paste(local_directory,control_file_name,"_Collections.csv",sep=""),row.names = F);
	write.csv(evanesco_na_from_matrix(control_data$control_occurrences,""),file=paste(local_directory,control_file_name,"_Finds.csv",sep=""),row.names = F);
	} else	{
	control_collections <- ingroup_collections;
	control_occurrences <- ingroup_finds;
	}

#if (!is.null(zone_taxa))	{
if (zone_taxa[1]!="")	{
	print("Getting occurrence data for zone taxa occupying the same collections....")
	zone_taxa_data <- accersi_data_for_control_groups_to_seed_FBD_analyses(control_taxon=zone_taxa,onset,end,basic_environments,species_only=T,paleogeography=paleogeography)
	#zone_collections <- zone_taxa_data$control_collections[zone_taxa_data$control_collections$collection_no %in% control_collections$collection_no];
	zone_occurrences <- zone_taxa_data$control_occurrences[zone_taxa_data$control_occurrences$collection_no %in% control_collections$collection_no,];
	write.csv(zone_taxa_data$control_occurrences[zone_taxa_data$control_occurrences$collection_no %in% control_collections$collection_no,],
			  file=paste(local_directory,analysis_name,"_Zone_Taxa_Finds.csv",sep=""),
			  row.names = F);
	}

# this provides edits to biogeography due to old glitches.
print("Now some basic cleaning of the collections data....")

if (is.data.frame(fossilworks_collections))	{
#	control_collections <- reparo_paleodb_paleogeography_with_fossilworks_data(paleodb_collections=control_collections,fossil_works_geography=fossilworks_collections);
	direct_dates <- data.frame(direct_ma=as.numeric(fossilworks_collections$direct_ma[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   direct_ma_error=as.numeric(fossilworks_collections$direct_ma_error[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   direct_ma_method=as.character(fossilworks_collections$direct_ma_method[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   stringsAsFactors=hell_no);
	direct_dates <- evanesco_na_from_matrix(direct_dates,"");
	control_collections <- cbind(control_collections,direct_dates);
	}

# this provides en masse edits for rock units used in paleodb.

if (is.data.frame(paleodb_rock_reidentifications))
	control_collections <- reparo_unedittable_paleodb_rock_identification(paleodb_collections=control_collections,paleodb_rock_reidentifications=paleodb_rock_reidentifications);

# this provides edits for paleodb collections that cannot currently be edited.
if (is.data.frame(paleodb_collection_edits))
	control_collections <- reparo_unedittable_paleodb_collections(paleodb_collections=control_collections,paleodb_collection_edits=paleodb_collection_edits);

# Correct the age ranges of collections using Gradstein 2012 + addenda (or another time scale)
control_collections <- redate_paleodb_collections_with_time_scale(paleodb_collections=control_collections,time_scale,zone_database);
if (zone_taxa[1]!="")	{
	output <- list(control_collections,control_occurrences,zone_occurrences);
	names(output) <- c("control_collections","control_occurrences","zone_occurrences");
	} else	{
	output <- list(control_collections,control_occurrences);
	names(output) <- c("control_collections","control_occurrences");
	}
return(output);
}

# updated 2020-02-20
# updated 2020-03-05
# updated 2020-04-12
# updated 2020-04-30
# updated 2020-05-04
# updated 2020-06-14
# updated 2020-12-15
accersi_PaleoDB_data_from_chosen_nexus_file <- function(onset,end,rock_unit_databases,chronostratigraphic_databases,paleodb_fixes,control_taxon="",zone_taxa="",taxon_level="species",basic_environments=c("marine","unknown","terrestrial"),paleogeography="scotese",time_scale_stratigraphic_scale="International",temporal_precision=0.05,lump_subgenera=F,analysis_name="",local_directory="",exclude_uncertain_taxa=T,species_only=T,bogarted=F,taxon_subset_file=F)	{
#### PART 0: Commence ####
if (time_scale_stratigraphic_scale=="Standard")
	time_scale_stratigraphic_scale <- "International";
print("This program will read a Nexus file and download collections and occurrences from the Paleobiology Database");
print("   (https://www.paleobiodb.org/) for stratigraphic data and then start refining/cleaning/updating those data");
print("   with updated time scales and biozonation information.");
print("");
print("NOTE: The Paleobiology Database should always be considered a STARTING point for these data. Part of what I have");
print("   designed the output to do is to let you vett occurrence and collection data for your study group.  We encourage");
print("   you to contribute updates to these data (to collections, occurrences and/or taxonomy) to the Paleobiology");
print("   Database (see https://www.youtube.com/channel/UCHxfFXYjYFotJmo_fNTSKJg for tutorials.)  Improvements to the");
print("   stratigraphic database used to refine PaleoDB data should be sent to pjwagner@gmail.com");
print("");
if (taxon_level=="genus" && !lump_subgenera)
	taxon_level <- "subgenus";
#### PART 1: GET TAXON INFORMATION FROM NEXUS FILE ####
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	flush.console();
	for (i in 1:100)	j <- 1;
	} else	{   
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = F,stringsAsFactors=hell_no)[,1];
	verboten <- c("taxon","taxa","species","genus","otu");
	taxa_subset <- taxa_subset[!taxa_subset %in% verboten];
	}

basic_data <- accersi_data_from_chosen_nexus_file();
otu_names_used <- basic_data$OTUs;
taxon_names <- otu_names_used[!tolower(otu_names_used) %in% c("outgroup","out")];
otu_names <- sapply(taxon_names,mundify_taxon_names);
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% basic_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!
if (species_only)	{
	taxon_names <- ingroup_taxa;
	clade_members <- unique(sapply(taxon_names,divido_genus_names_from_species_names));
	} else	{
	clade_members <- ingroup_taxa;
	}

#### PART 2: LOAD EXTERNAL DATA FOR CLEANING & REFINING PALEODB DATA  ####
fossilworks_collections <- paleodb_fixes$fossilworks_collections;
paleodb_rock_reidentifications <- paleodb_fixes$paleodb_rock_reidentifications;
paleodb_collection_edits <- paleodb_fixes$paleodb_collection_edits;
if (!is.null(paleodb_collection_edits$X))
	paleodb_collection_edits$X <- NULL;
time_scale <- chronostratigraphic_databases$time_scale;
zone_database <- chronostratigraphic_databases$zones;
if (is.list(rock_unit_databases))	{
	rock_database <- rock_unit_databases$rock_unit_database;
	rock_to_zone_database <- rock_unit_databases$rock_to_zone_database;
	rock_to_zone_database$ma_lb <- temporal_precision*round(rock_to_zone_database$ma_lb/temporal_precision,0);
	rock_to_zone_database$ma_ub <- temporal_precision*round(rock_to_zone_database$ma_ub/temporal_precision,0);
	}
time_scale$ma_lb <- temporal_precision*round(time_scale$ma_lb/temporal_precision,0);
time_scale$ma_ub <- temporal_precision*round(time_scale$ma_ub/temporal_precision,0);
zone_database$ma_lb <- temporal_precision*round(as.numeric(zone_database$ma_lb)/temporal_precision,0);
zone_database$ma_ub <- temporal_precision*round(zone_database$ma_ub/temporal_precision,0);

zone_database <- subset(zone_database,zone_database$ma_lb<=time_scale$ma_lb[match(onset,time_scale$interval)]+5);
zone_database <- subset(zone_database,zone_database$ma_ub>=time_scale$ma_ub[match(end,time_scale$interval)]-5);

#### PART 3: GET INFORMATION NEEDED TO DOWNLOAD, 'CLEAN' AND ANALYZE STRATIGRAPHIC DATA  ####
compendium <- accersi_updated_taxonomy_for_analyzed_taxa(otu_names=otu_names,local_directory=local_directory,study=analysis_name);

if (bogarted)	{
	print("Choose the file with your private stash information: ");
	flush.console();
	Sys.sleep(zzzz);
	bogarted_info <- file.choose();
	print("Reading your private stash now....");
	flush.console();
	Sys.sleep(zzzz);
#	bogarted_finds <- utils::read.csv(file = read.(bogarted_info), header = TRUE,stringsAsFactors=FALSE,encoding = "UTF-8");
	bogarted_finds <- read.csv(file=bogarted_info,header = T,stringsAsFactors=hell_no,encoding = "UTF-8");
	bogarted_finds <- evanesco_na_from_matrix(bogarted_finds,replacement="");
	bogarted_finds <- subset(bogarted_finds,bogarted_finds$identified_name!="");
	if (!is.na(match("paleodb_collection_no",colnames(bogarted_finds))))	{
		ccc <- colnames(bogarted_finds);
		ccc[match("collection_no",ccc)] <- "my_collection_no";
		ccc[match("paleodb_collection_no",ccc)] <- "collection_no";
		colnames(bogarted_finds) <- ccc;
		}
		
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_name <- bogarted_finds$identified_name;
		bogarted_finds$genus <- as.character(sapply(taxon_name,divido_genus_names_from_species_names));
		if (taxon_level=="subgenus")	{
			genus_name <- bogarted_finds$genus;
			subgenus_results <- sapply(genus_name,divido_subgenus_names_from_genus_names);
#			bogarted_finds$genus <- subgenus_results[1,];
			bogarted_finds$subgenus <- subgenus_results[2,];
			bogarted_finds$subgenus[bogarted_finds$subgenus==""] <- bogarted_finds$genus[bogarted_finds$subgenus==""];
			}
	#	add occurrences
		unique_genera <- unique(bogarted_finds$genus);
		if (!is.null(bogarted_finds$subgenus))
			unique_genera <- unique(c(bogarted_finds$genus,bogarted_finds$subgenus));
		for (u_g in 1:length(unique_genera))	{
			if (!is.na(match(unique_genera[u_g],compendium$taxon_name)))	{
				compendium$n_occs[match(unique_genera[u_g],compendium$taxon_name)] <- length(unique(bogarted_finds$collection_no[unique(which(bogarted_finds==unique_genera[u_g],arr.ind = T)[,1])]));
				}
			}
		}
	if (!is.null(bogarted_finds$direct_ma))	{
		bogarted_finds$direct_ma <- as.numeric(bogarted_finds$direct_ma);
		bogarted_finds$direct_ma[is.na(bogarted_finds$direct_ma)] <- 0;
		}
	if (!is.null(bogarted_finds$direct_ma_error))	{
		bogarted_finds$direct_ma_error <- as.numeric(bogarted_finds$direct_ma_error);
		bogarted_finds$direct_ma_error[is.na(bogarted_finds$direct_ma_error)] <- 0;
		}
	bogarted_finds$max_ma <- time_scale$ma_lb[match(bogarted_finds$early_interval,time_scale$interval)];
	bogarted_finds$late_interval[bogarted_finds$late_interval==""] <- bogarted_finds$early_interval[bogarted_finds$late_interval==""];
	bogarted_finds$min_ma <- time_scale$ma_ub[match(bogarted_finds$late_interval,time_scale$interval)];
	bogarted_finds$accepted_name_orig <- bogarted_finds$accepted_name;
	bogarted_taxa <- unique(bogarted_finds$identified_name);
	for (bt in 1:length(bogarted_taxa))	{
		btn <- match(bogarted_taxa[bt],compendium$taxon_name);
		if (is.na(btn))
			btn <- match(bogarted_taxa[bt],otu_names);
		
		if (!is.na(btn))	{
			this_taxon <- subset(bogarted_finds,bogarted_finds$identified_name==bogarted_taxa[bt]);
#			compendium$n_occs[btn] <- compendium$n_occs[btn] + sum(bogarted_finds$identified_name==bogarted_taxa[bt]);
			compendium$n_occs[btn] <- compendium$n_occs[btn] + nrow(this_taxon);
			compendium$early_interval[btn] <- this_taxon$early_interval[match(max(this_taxon$max_ma),this_taxon$max_ma)];
			if (is.na(this_taxon$late_interval[match(min(this_taxon$min_ma),this_taxon$min_ma)]))	{
				compendium$late_interval[btn] <- this_taxon$early_interval[match(min(this_taxon$min_ma),this_taxon$min_ma)];
				} else	{
				compendium$late_interval[btn] <- this_taxon$late_interval[match(min(this_taxon$min_ma),this_taxon$min_ma)];
				}
			compendium$firstapp_max_ma[btn] <- time_scale$ma_lb[match(compendium$early_interval[btn],time_scale$interval)];
			compendium$firstapp_min_ma[btn] <- time_scale$ma_ub[match(compendium$early_interval[btn],time_scale$interval)];
			compendium$lastapp_max_ma[btn] <- time_scale$ma_lb[match(compendium$late_interval[btn],time_scale$interval)];
			compendium$lastapp_min_ma[btn] <- time_scale$ma_ub[match(compendium$late_interval[btn],time_scale$interval)];
			}
		}
	}

if (sum(compendium$n_occs==0)>0)	{
	missing_taxa <- subset(compendium,compendium$n_occs==0);
	missing_taxa <- subset(missing_taxa,missing_taxa$accepted_name=="?");
	taxon_name <- missing_taxa$taxon_name;
	missing_taxa_rows <- match(taxon_name,compendium$taxon_name);
	taxon_list <- genera <- unique(sapply(taxon_name,divido_genus_names_from_species_names));
	if (length(taxon_list)>0)	{
		taxonomyless_finds <- accersi_occurrences_for_list_of_taxa(taxon_list,paleogeography=paleogeography);
		if (is.data.frame(taxonomyless_finds$collection_compendium))
			for (mt in 1:length(missing_taxa))
				compendium$n_occs[missing_taxa_rows[mt]] <- sum(taxonomyless_finds$occurrences_compendium$identified_name==taxon_name[mt]);
		}
	### insert command for getting occurrences & collections for lists of taxa here.
	}	else	{
	missing_taxa <- "";
	}

if (sum(compendium$n_occs==0)>0)	{
#	print(paste("The following taxa have no occurrences:",paste(compendium$taxon_name[compendium$n_occs==0],collapse=", ")));
	print("The following taxa currently have no occurrences entered into the PaleoDB:");
	print(compendium$taxon_name[compendium$n_occs==0]);
#	print(paste("The following taxa are not entered into the PaleoDB:",paste(compendium$taxon_name[compendium$taxon_no==""],collapse=", ")));
	print("");
	if (sum(compendium$taxon_no=="")>0)	{
		print("The following taxa are not entered into the PaleoDB taxonomy tables:");
		print(compendium$taxon_name[compendium$taxon_no==0]);
		}
	print("Enter Data for these into the PaleoDB and try again tomorrow or setup a separate 'bogarted' file with occurrences for these taxa!");
	print("   Make sure the file as formation, member, stage, zonation, etc., information, too. (And consider entering it into the PaleoDB later.)");
	print("Also, make sure that there are no misspellings in your nexus matrix. (Computers do not autocorrect!)");
	return();
	}
otu_names[compendium$accepted_name!="?"] <- compendium$accepted_name[compendium$accepted_name!="?"];

if (abs(time_scale$ma_lb[match(onset,time_scale$interval)])<max(abs(compendium$firstapp_max_ma)))	{
	stage_info <- accersi_stage_info();
	stage_info <- subset(stage_info,abs(stage_info$onset)>max(abs(compendium$firstapp_max_ma)));
	stage_info <- stage_info[order(abs(stage_info$onset)),];
	if (nrow(stage_info)>1)	{
		onset <- stage_info$interval[2];
		} else	{
		onset <- stage_info$interval[1];
		}
	if (!is.na(match(onset,c("Stage 2","Stage 3","Stage 4"))))	{
		onset <- c("Meishucunian","Atdabanian","Duyunian")[match(end,c("Stage 2","Stage 3","Stage 4"))]
		}
	}

if (abs(time_scale$ma_ub[match(end,time_scale$interval)])>min(abs(compendium$firstapp_max_ma)))	{
	stage_info <- accersi_stage_info();
	stage_info <- subset(stage_info,abs(stage_info$end)<min(abs(compendium$firstapp_max_ma)));
	bb <- 1+sum(stage_info$end>min(abs(compendium$firstapp_max_ma)));
	stage_info <- stage_info[order(-abs(stage_info$onset)),];
	end <- stage_info$interval[bb];
	if (!is.na(match(end,c("Stage 2","Stage 3","Stage 4"))))
		end <- c("Tommotian","Nangaoian","Duyunian")[match(end,c("Stage 2","Stage 3","Stage 4"))]
	}

## get paleodb data!!!!
#                                         function(otu_names,analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa=T,taxon_level,onset,end,basic_environments=c("marine","unknown","terrestrial"),paleogeography="scotese",time_scale,zone_database,fossilworks_collections="",paleodb_rock_reidentifications="",paleodb_collection_edits="",lump_subgenera=F,species_only=T)
paleodb_data <- accersi_paleodb_data_for_Rev_Bayes(otu_names,analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa,taxon_level,onset,end,basic_environments,paleogeography,time_scale,zone_database,fossilworks_collections,paleodb_rock_reidentifications,paleodb_collection_edits,lump_subgenera,species_only);
control_collections <- unique(paleodb_data$control_collections);
control_occurrences <- unique(paleodb_data$control_occurrences);

if (bogarted)	{
	print("Adding your private stash to the PaleoDB data....");
	if (!is.na(match("my_collection_no",colnames(bogarted_finds))))
		bogarted_finds$my_collection_no <- as.numeric(bogarted_finds$my_collection_no);
	bogarted_finds$collection_no[bogarted_finds$collection_no==""] <- 0;
	bogarted_finds$collection_no <- as.numeric(bogarted_finds$collection_no);
	bogarted_finds$collection_no[bogarted_finds$collection_no==0] <- 
		bogarted_finds$my_collection_no[bogarted_finds$collection_no==0]+ceiling(max(control_collections$collection_no)/10^(floor(log10(max(control_collections$collection_no)))-1))*10^(floor(log10(max(control_collections$collection_no)))-1);

	column_matches <- match(colnames(bogarted_finds),colnames(control_collections));
	bogarted_coll_info_in_paleodb <- (1:ncol(bogarted_finds))[!is.na(column_matches)];
	column_matches <- column_matches[!is.na(column_matches)];
	new_paleodb_coll <- control_collections[1:length(unique(bogarted_finds$collection_no)),];
	for (nc in 1:ncol(new_paleodb_coll))	{
		if (is.numeric(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- 0;
			} else if (is.character(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- "";
			}
		}
#	new_paleodb_coll <- control_collections[length(unique(bogarted_finds$collection_no)),];
	new_paleodb_coll[,column_matches] <- bogarted_finds[match(unique(bogarted_finds$collection_no),bogarted_finds$collection_no),bogarted_coll_info_in_paleodb];
	control_collections <- rbind(control_collections,new_paleodb_coll);
	
	# set up occcurrences in two steps;
	# edit already downloaded occurrences
	emended_paleodb_finds <- subset(bogarted_finds,bogarted_finds$occurrence_no %in% control_occurrences$occurrence_no);
	edit_paleodb_rows <- match(emended_paleodb_finds$occurrence_no,control_occurrences$occurrence_no);
	column_matches <- match(colnames(bogarted_finds),colnames(control_occurrences))
	column_matches <- column_matches[!is.na(column_matches)];
	matched_columns <- (1:ncol(emended_paleodb_finds))[!is.na(match(colnames(bogarted_finds),colnames(control_occurrences)))];
	control_occurrences[edit_paleodb_rows,column_matches] <- emended_paleodb_finds[,matched_columns];
	control_occurrences$accepted_name[edit_paleodb_rows] <- control_occurrences$accepted_name_orig[edit_paleodb_rows] <- emended_paleodb_finds$identified_name;
	# add completely new finds
	totally_boggy <- subset(bogarted_finds,!bogarted_finds$occurrence_no %in% control_occurrences$occurrence_no);
	new_paleodb_occr <- control_occurrences[1:nrow(totally_boggy),];
	for (nc in 1:ncol(new_paleodb_occr))	{
		if (is.numeric(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.numeric(0);
			} else if (is.character(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.character("");
			new_paleodb_occr[,nc] <- as.character(new_paleodb_occr[,nc]);
			}
		}
	
	new_paleodb_occr[,column_matches] <- totally_boggy[,matched_columns];

	# now get the taxonomy part....
	totally_bogarted_taxa <- unique(totally_boggy$identified_name);
	for (tt in 1:length(totally_bogarted_taxa))	{
		if (tt==1)	{
			bogarted_taxonomy <- revelare_taxonomy_for_one_taxon(taxon=totally_bogarted_taxa[tt],settle=T);
			} else	{
			bogarted_taxonomy <- rbind(bogarted_taxonomy,revelare_taxonomy_for_one_taxon(taxon=totally_bogarted_taxa[tt],settle=T));
			}
		informal_taxon <- revelare_informal_taxa(taxon_name=totally_bogarted_taxa[tt]);
		if (informal_taxon)	{
			bogarted_taxonomy$taxon_name[tt] <- bogarted_taxonomy$accepted_name[tt] <- totally_bogarted_taxa[tt];
			bogarted_taxonomy$accepted_rank[tt] <- "species";
			}
		}
	bogarted_taxonomy$accepted_name[bogarted_taxonomy$accepted_name!=totally_bogarted_taxa] <- totally_bogarted_taxa[bogarted_taxonomy$taxon_name!=totally_bogarted_taxa];
	bogarted_taxonomy$taxon_name[bogarted_taxonomy$taxon_name!=totally_bogarted_taxa] <- totally_bogarted_taxa[bogarted_taxonomy$taxon_name!=totally_bogarted_taxa];
	bogarted_taxonomy$accepted_rank[match(bogarted_taxonomy$accepted_rank,taxonomic_rank)>match(taxon_level,taxonomic_rank)] <- taxon_level;
	bogarted_taxonomy <- evanesco_na_from_matrix(bogarted_taxonomy,replacement = "");
	bogarted_taxonomy$record_type <- bogarted_taxonomy$flags <- bogarted_taxonomy$early_interval <- bogarted_taxonomy$late_interval <- NULL;
	bogarted_row_to_paledob <- match(totally_boggy$identified_name,bogarted_taxonomy$taxon_name);
	paleodb_col_to_edit <- match(colnames(bogarted_taxonomy),colnames(control_occurrences));
	paleodb_col_to_edit <- paleodb_col_to_edit[!is.na(paleodb_col_to_edit)];
	bogarted_col_w_fix <- match(colnames(control_occurrences)[paleodb_col_to_edit],colnames(bogarted_taxonomy));
	new_paleodb_occr[,paleodb_col_to_edit] <- bogarted_taxonomy[bogarted_row_to_paledob,bogarted_col_w_fix];
	new_paleodb_occr$record_type <- control_occurrences$record_type[1];
	new_paleodb_occr$accepted_name_orig <- new_paleodb_occr$accepted_name;
	
	new_paleodb_occr$occurrence_no <- (1:nrow(new_paleodb_occr))+(2*max(as.numeric(control_occurrences$occurrence_no)));
	control_occurrences <- rbind(control_occurrences,new_paleodb_occr);
	}

this_taxon_rank <- c();
for (tx in 1:length(otu_names))	{
	# if species or subspecies
	if (revelare_informal_taxa(taxon_name = otu_names[tx]))	{
		this_taxon_rank <- c(this_taxon_rank,"species");
		}	else if (length((strsplit(otu_names[tx]," ")[[1]]))==2)	{
		second_name <- strsplit(otu_names[tx]," ")[[1]][2];
		first_character <- strsplit(second_name,"")[[1]][1];
		if (first_character=="(")	{
			this_taxon_rank <- c(this_taxon_rank,"subgenus");
			} else	{
			this_taxon_rank <- c(this_taxon_rank,"species");
			}
		} else if (length((strsplit(otu_names[tx]," ")[[1]]))==1)	{
		this_taxon_rank <- c(this_taxon_rank,"genus");
		} else if (length((strsplit(otu_names[tx]," ")[[1]]))==3)	{
		this_taxon_rank <- c(this_taxon_rank,"subspecies");
		}
	if (this_taxon_rank[tx]=="species" || this_taxon_rank[tx]=="subspecies")	{
		taxon_finds <- unique(rbind(subset(control_occurrences,control_occurrences$accepted_name==otu_names[tx]),subset(control_occurrences,control_occurrences$accepted_name_orig==otu_names[tx])));
		} else	{
		taxon_finds <- subset(control_occurrences,control_occurrences$genus==otu_names[tx]);
		}
	if (nrow(taxon_finds)==0)	{
		taxon_finds <- accersi_occurrence_data(taxa=otu_names[tx],species_only = species_only,save_files=F);
		if (this_taxon_rank=="species" || this_taxon_rank=="subspecies")	{
			taxon_finds$accepted_name <- otu_names[tx];
			} else	{
			taxon_finds$genus <- otu_names[tx];
			}
		emend_these <- (1:nrow(control_occurrences))[control_occurrences$occurrence_no %in% taxon_finds$occurrence_no];
		updates <- (1:nrow(taxon_finds))[taxon_finds$occurrence_no %in% control_occurrences$occurrence_no];
		if (length(emend_these)>0)
			control_occurrences[emend_these,] <- taxon_finds[updates,];
		newbies <- (1:nrow(taxon_finds))[!taxon_finds$occurrence_no %in% control_occurrences$occurrence_no];
		if (length(updates)>0 || sum(taxon_finds$flags[newbies]!="")<length(newbies))	{ 
			newbies <- newbies[taxon_finds$flags[newbies]==""];
			} else {
			# add questionable assignments only as a last resort
			taxon_finds$flags[newbies] <- ""
			}
		if (length(newbies)>0)	{
			control_occurrences <- rbind(control_occurrences,taxon_finds[newbies,]);
			control_occurrences <- control_occurrences[order(control_occurrences$collection_no,control_occurrences$occurrence_no),];
			}
		}
	relv_finds <- match(taxon_finds$occurrence_no,control_occurrences$occurrence_no);
	if (this_taxon_rank[tx]=="species" || this_taxon_rank[tx]=="subspecies")	{
		control_occurrences$accepted_name[relv_finds] <- otu_names[tx];
		} else	{
		control_occurrences$genus[relv_finds] <- otu_names[tx];
		}
	}

### look for occurrence collections not in control collections
occ_colls <- sort(unique(control_occurrences$collection_no));
missing_colls <- occ_colls[!occ_colls %in% control_collections$collection_no];
if (length(missing_colls)>0)	{
	for (mc in 1:length(missing_colls))	{
		xxx <- accersi_single_locality_info(missing_colls[mc],paleogeography=paleogeography);
		if (!"direct_ma" %in% colnames(xxx))	{
			xxx$direct_ma <- 0;
			xxx$direct_ma_error <- xxx$direct_ma_method <- "";
			}
		if (missing_colls[mc] %in% paleodb_collection_edits$collection_no)	{
			fixes <- match(colnames(paleodb_collection_edits),colnames(xxx))
			xxx[,fixes] <- paleodb_collection_edits[match(xxx$collection_no),paleodb_collection_edits$collection_no,];
			}
		control_collections <- rbind(control_collections,xxx);
		}
	control_collections <- control_collections[order(control_collections$collection_no),];
	}

if (taxon_level=="genus" || taxon_level=="subgenus")	
	control_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = control_occurrences);

if (is.data.frame(paleodb_data$zone_occurrences))	{
	zone_occurrences <- paleodb_data$zone_occurrences;
	if (taxon_level=="genus" || taxon_level=="subgenus")
		zone_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = zone_occurrences);
	} else	{
	zone_occurrences <- NULL;
	}

hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(control_collections$early_interval)),unique(as.character(control_collections$late_interval)))),time_scale,regional_scale=time_scale_stratigraphic_scale);
hierarchical_chronostrat$ma_lb <- temporal_precision*round(hierarchical_chronostrat$ma_lb/temporal_precision,0);
hierarchical_chronostrat$ma_ub <- temporal_precision*round(hierarchical_chronostrat$ma_ub/temporal_precision,0);
### If this does not encompass study interval, then modify! ###
if (time_scale$ma_lb[match(onset,time_scale$interval)]>max(hierarchical_chronostrat$ma_lb))	{
	added_scale <- accersi_hierarchical_timescale(chronostrat_units=c(onset,hierarchical_chronostrat$interval[1]),time_scale,regional_scale="International");
	added_end <- match(round(hierarchical_chronostrat$ma_lb[1],3),round(added_scale$ma_ub,3));
	if (!is.na(added_end))	{
		added_scale <- added_scale[1:added_end,];
		} else	{
		ttl_bins <- nrow(added_scale); 
		closest_end <- (1:ttl_bins)[hierarchical_chronostrat$ma_lb[1]-added_scale$ma_ub>0][1];
		added_scale <- added_scale[1:closest_end,];
		added_scale$ma_ub[closest_end] <- hierarchical_chronostrat$ma_lb[1];
		}
	dummy_time_scale <-  rbind(added_scale,hierarchical_chronostrat);
	dummy_time_scale$scale <- time_scale_stratigraphic_scale;
	dummy_time_scale$parent_interval <- dummy_time_scale$bin_first <- dummy_time_scale$bin_last <- NULL;
	hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=dummy_time_scale$interval,time_scale=dummy_time_scale,regional_scale=time_scale_stratigraphic_scale);
	}
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
ranges <- paste(finest_chronostrat$ma_lb,finest_chronostrat$ma_ub,sep="-");
finest_chronostrat <- finest_chronostrat[match(unique(ranges),ranges),];

downloaded_collections <- reset_paleodb_intervals_to_desired_time_scale(collections=control_collections,finest_chronostrat = finest_chronostrat,time_scale);
downloaded_collections$min_ma <- round(temporal_precision*round(downloaded_collections$min_ma/temporal_precision,0),floor(-log10(temporal_precision)));
downloaded_collections$max_ma <- round(temporal_precision*round(downloaded_collections$max_ma/temporal_precision,0),floor(-log10(temporal_precision)));
ncolls <- nrow(downloaded_collections);

#### PART 4: REFINE CHRONOSTRATIGRAPHY OF PALEODB DATA  ####
if (!is.null(zone_occurrences))	{
	paleodb_finds <- rbind(control_occurrences,zone_occurrences);
	paleodb_finds <- paleodb_finds[order(paleodb_finds$collection_no,paleodb_finds$occurrence_no),];
	paleodb_finds <- paleodb_finds[match(unique(paleodb_finds$occurrence_no),paleodb_finds$occurrence_no),];
	} else	{
	paleodb_finds <- control_occurrences;
	}
if (is.list(rock_unit_databases))	{
	print("Refining PaleoDB data with rock-unit and biozonation databases...");
#	time_scale$ma_ub[match(onset,time_scale$interval)];
#	time_scale$ma_ub[match(end,time_scale$interval)];
	lb <- time_scale$ma_lb[match(onset,time_scale$interval)]+round(0.05*time_scale$ma_lb[match(onset,time_scale$interval)],1);
	ub <- time_scale$ma_ub[match(end,time_scale$interval)]-round(0.05*time_scale$ma_ub[match(end,time_scale$interval)],1);
	relv_rock_database <- subset(rock_database,rock_database$ma_lb<=lb);
	relv_rock_database <- subset(relv_rock_database,relv_rock_database$ma_ub>=ub);
	relv_rock_to_zone_database <- subset(rock_to_zone_database,rock_to_zone_database$rock_no %in% relv_rock_database$rock_no);
	relv_zone_database <- subset(zone_database,zone_database$ma_lb>=time_scale$ma_ub[match(end,time_scale$interval)]);
	relv_zone_database <- subset(relv_zone_database,relv_zone_database$ma_ub<=time_scale$ma_lb[match(onset,time_scale$interval)]);
	
	paleodb_data_refined <- refine_collection_dates_with_external_database(study=analysis_name,collections=downloaded_collections,rock_database=relv_rock_database,zone_database=relv_zone_database,rock_to_zone_database=relv_rock_to_zone_database,time_scale,directory=local_directory);
	refined_collections <- paleodb_data_refined$Recalibrated_Collections;
	chronostrat_units <- unique(c(hierarchical_chronostrat$interval,refined_collections$interval_lb,refined_collections$interval_ub));
	hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units,time_scale,regional_scale=time_scale_stratigraphic_scale);
	finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
	} else if (is.data.frame(zone_database))	{
	print("Refining PaleoDB data with biozonation databases...");
	relv_zone_database <- subset(zone_database,zone_database$ma_lb>=time_scale$ma_ub[match(end,time_scale$interval)]);
	relv_zone_database <- subset(relv_zone_database,relv_zone_database$ma_ub<=time_scale$ma_lb[match(onset,time_scale$interval)]);
	refined_collections <- refine_paleodb_collection_dates_with_zone_data_only(paleodb_collections=downloaded_collections,paleodb_finds,zone_database=relv_zone_database,time_scale,hierarchical_chronostrat,finest_chronostrat,examine_finds=T,temporal_precision=0.05);
	} else	{
	refined_collections <- downloaded_collections;
	refined_collections$ma_lb <- refined_collections$max_ma;
	refined_collections$ma_ub <- refined_collections$min_ma;
	refined_collections$interval_lb <- refined_collections$early_interval;
	refined_collections$interval_ub <- refined_collections$late_interval;
	}

if (!is.null(refined_collections$direct_ma) && sum(refined_collections$direct_ma>0)>0)	{
	print("Using radiometric data for final ages...");
	if(max(as.numeric(refined_collections$direct_ma[refined_collections$direct_ma!=""]))>max(finest_chronostrat$ma_lb))	{
		old_dates <- as.numeric(refined_collections$direct_ma[refined_collections$direct_ma!=""]);
		new_old_dates <- old_dates[old_dates>max(finest_chronostrat$ma_lb)];
		new_old_dates <- sort(new_old_dates);
		for (oldies in 1:length(new_old_dates))	{
			relv_intervals <- subset(time_scale,time_scale$ma_ub<new_old_dates[oldies]);
			relv_intervals <- subset(relv_intervals,relv_intervals$ma_lb>=new_old_dates[oldies]);
			interval_spans <- relv_intervals$ma_lb-relv_intervals$ma_ub;
			new_interval <- relv_intervals[match(min(interval_spans),interval_spans),];
			parent_interval <- new_interval$interval;
			new_interval <- tibble::add_column(new_interval, parent_interval, .after = 1);
			bin_first <- bin_last <- min(finest_chronostrat$bin_first)-1;
			new_interval <- tibble::add_column(new_interval, bin_last, .after = 6);
			new_interval <- tibble::add_column(new_interval, bin_first, .after = 6);
			finest_chronostrat <- rbind(new_interval,finest_chronostrat);
			hierarchical_chronostrat <- rbind(new_interval,hierarchical_chronostrat);
			}
		finest_chronostrat <- finest_chronostrat[match(unique(finest_chronostrat$interval),finest_chronostrat$interval),];
		finest_chronostrat$bin_first <- order(finest_chronostrat$bin_first);
		finest_chronostrat$bin_last <- order(finest_chronostrat$bin_last);
		hierarchical_chronostrat <- hierarchical_chronostrat[match(unique(hierarchical_chronostrat$interval),hierarchical_chronostrat$interval),];
		hierarchical_chronostrat$bin_first <- order(hierarchical_chronostrat$bin_first);
		hierarchical_chronostrat$bin_last <- order(hierarchical_chronostrat$bin_last);
		#		new_oldest <- max(as.numeric(refined_collections$direct_ma[refined_collections$direct_ma!=""]));
		}
	refined_collections <- redate_collections_with_direct_dates(collections=refined_collections,finest_chronostrat,temporal_precision = 0.1);
	}

age <- temporal_precision*round(refined_collections$ma_lb/temporal_precision,0);
#age <- sort(unique(temporal_precision*round(refined_collections$ma_lb/temporal_precision,0)));
#xxx <- as.character(sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat));
refined_collections$interval_lb <- as.character(sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat));
bad_colls <- (1:nrow(refined_collections))[refined_collections$interval_lb=="character(0)"];
#refined_collections <- subset(refined_collections,refined_collections$interval_lb!="character(0)");
age <- temporal_precision*round(refined_collections$ma_ub/temporal_precision,0);
#age <- sort(unique(age));
#xxx <- as.character(sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat));
refined_collections$interval_ub <- as.character(sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat));
bad_colls <- sort(unique(c(bad_colls,(1:nrow(refined_collections))[refined_collections$interval_ub=="character(0)"])));
paleodb_finds <- subset(paleodb_finds,!paleodb_finds$collection_no %in% refined_collections$collection_no[bad_colls]);
refined_collections <- refined_collections[refined_collections$collection_no %in% paleodb_finds$collection_no,];
#sum(refined_collections$ma_lb<=refined_collections$ma_ub)

# use quantitative biostratigraphy 101 to refine dates here.
print("Using basic biostratigraphy to minimize gaps for uncertainly aged collections...");
if (is.null(relv_zone_database) && !is.data.frame(zone_database))	{
	relv_zone_database <- "";
	}
# refined_collections lacking some of the collections in finds!!!
optimized_collections <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds=paleodb_finds,paleodb_collections=refined_collections,hierarchical_chronostrat,zone_database=relv_zone_database,update_search=T);
#ddd <- (1:ncolls)[optimized_collections$ma_lb<=optimized_collections$ma_ub]

if (is.null(optimized_collections$bin_lb))
	optimized_collections$bin_lb <- as.numeric(finest_chronostrat$bin_first[match(optimized_collections$interval_lb,finest_chronostrat$interval)]);
if (is.null(optimized_collections$bin_ub))
	optimized_collections$bin_ub <- as.numeric(finest_chronostrat$bin_last[match(optimized_collections$interval_ub,finest_chronostrat$interval)]);
# get rock unit numbers if we do not have a stratigraphic database
if (is.null(optimized_collections$rock_no))	{
	print("Putting numbers on unique rock units...");
	optimized_collections <- number_unique_rock_units(paleodb_collections = optimized_collections,zone_database=relv_zone_database,time_scale=finest_chronostrat);
	}

# for unentered rock units that are unique to their time and location, create dummy numbers.
optimized_collections <- name_unnamed_rock_units(paleodb_collections=optimized_collections,finest_chronostrat);
ncolls <- nrow(optimized_collections);

#### PART 4A: last try to find problem children! ####
optimized_collections_orig <- optimized_collections;
problem_collections <- (1:ncolls)[optimized_collections$ma_lb<=optimized_collections$ma_ub];
pz <- 0;
while (pz < length(problem_collections))	{
	pz <- pz+1;
	if (is.list(rock_database))	{
		prob_formation <- optimized_collections$formation[problem_collections[pz]];
		prob_member <- optimized_collections$member[problem_collections[pz]];
		} else	{
		prob_formation <- prob_member <- "";
		}
	if (is.data.frame(relv_zone_database))	{
		prob_zone <- optimized_collections$zone[problem_collections[pz]];
		} else	{
		prob_zone <- "";
		}
	if (prob_formation!="")	{
		ddd <- data.frame(which(rock_database==prob_formation,arr.ind = T));
		ffff <- sort(unique(ddd$row));
		fff <- ffff[rock_database$member[ffff]==prob_member];
		if (length(fff)==0 && length(ffff)>0)	{
			prob_member <- "";
			fff <- ffff[rock_database$member[ffff]==prob_member];
			if (length(fff)==0)	fff <- ffff[1];
			}
		rock_range <- c(max(rock_database$ma_lb[fff]),min(rock_database$ma_ub[unique(ddd$row)]));
		if (prob_zone!="")	{
			zzz <- data.frame(which(zone_database==prob_zone,arr.ind = T));
			if (nrow(zzz)>0)	{
				zone_range <- c(max(zone_database$ma_lb[zzz$row]),min(zone_database$ma_ub[zzz$row]));
				overlap <- accersi_temporal_overlap (lb1=rock_range[1],lb2=zone_range[1],ub1=rock_range[2],ub2=zone_range[2]);
				if (overlap[1]>0)	{
					optimized_collections$ma_lb[problem_collections[pz]] <- as.numeric(overlap[1]);
					optimized_collections$ma_ub[problem_collections[pz]] <- as.numeric(overlap[2]);
					} else	{
					if (as.numeric(optimized_collections$ref_pubyr[problem_collections[pz]])>=1995)	{
						optimized_collections$ma_lb[problem_collections[pz]] <- as.numeric(zone_range[1]);
						optimized_collections$ma_ub[problem_collections[pz]] <- as.numeric(zone_range[2]);
						} else	{
						optimized_collections$ma_lb[problem_collections[pz]] <- as.numeric(rock_range[1]);
						optimized_collections$ma_ub[problem_collections[pz]] <- as.numeric(rock_range[2]);
						}
					}
				} else	{
				optimized_collections$ma_lb[problem_collections[pz]] <- as.numeric(rock_range[1]);
				optimized_collections$ma_ub[problem_collections[pz]] <- as.numeric(rock_range[2]);
				}
			} else	{
			optimized_collections$ma_lb[problem_collections[pz]] <- rock_range[1];
			optimized_collections$ma_ub[problem_collections[pz]] <- rock_range[2];
			}
		optimized_collections$interval_lb[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=optimized_collections$ma_lb[problem_collections[pz]])];
		optimized_collections$interval_ub[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>optimized_collections$ma_ub[problem_collections[pz]])];
		} else if (prob_zone!="")	{
		zzz <- data.frame(which(zone_database==prob_zone,arr.ind = T));
		if (nrow(zzz)>0)	{
#			print(pz);
			zone_range <- c(max(zone_database$ma_lb[zzz$row]),min(zone_database$ma_ub[zzz$row]));
			if (as.numeric(optimized_collections$ref_pubyr[problem_collections[pz]])>=1995)	{
				optimized_collections$ma_lb[problem_collections[pz]] <- as.numeric(zone_range[1]);
				optimized_collections$ma_ub[problem_collections[pz]] <- as.numeric(rock_range[2]);
				} else	{
				optimized_collections$ma_lb[problem_collections[pz]] <- optimized_collections$max_ma[problem_collections[pz]];
				optimized_collections$ma_ub[problem_collections[pz]] <- optimized_collections$min_ma[problem_collections[pz]];				
				}
			optimized_collections$interval_lb[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=optimized_collections$ma_lb[problem_collections[pz]])];
			optimized_collections$interval_ub[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>optimized_collections$ma_ub[problem_collections[pz]])];
			} else	{
			optimized_collections$ma_lb[problem_collections[pz]] <- optimized_collections$max_ma[problem_collections[pz]];
			optimized_collections$ma_ub[problem_collections[pz]] <- optimized_collections$min_ma[problem_collections[pz]];				
			optimized_collections$interval_lb[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=optimized_collections$ma_lb[problem_collections[pz]])];
			optimized_collections$interval_ub[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>optimized_collections$ma_ub[problem_collections[pz]])];
			}
		} else if (optimized_collections$bin_lb[problem_collections[pz]]<=optimized_collections$bin_ub[problem_collections[pz]])	{
		bin_onset <- finest_chronostrat$ma_lb[optimized_collections$bin_lb[problem_collections[pz]]];
		bin_end <- finest_chronostrat$ma_ub[optimized_collections$bin_lb[problem_collections[pz]]];
		if (optimized_collections$ma_lb[problem_collections[pz]]<bin_onset && optimized_collections$ma_lb[problem_collections[pz]]>bin_end)	{
			optimized_collections$ma_lb[problem_collections[pz]] <- optimized_collections$ma_lb[problem_collections[pz]]+temporal_precision;
			optimized_collections$ma_ub[problem_collections[pz]] <- optimized_collections$ma_ub[problem_collections[pz]]-temporal_precision;
			} else if (optimized_collections$ma_lb[problem_collections[pz]]==bin_onset)	{
			optimized_collections$ma_ub[problem_collections[pz]] <- optimized_collections$ma_ub[problem_collections[pz]]-temporal_precision;
			} else if (optimized_collections$ma_ub[problem_collections[pz]]==bin_end)	{
			optimized_collections$ma_lb[problem_collections[pz]] <- optimized_collections$ma_lb[problem_collections[pz]]+temporal_precision;
			} else	{
			optimized_collections$ma_lb[problem_collections[pz]] <- finest_chronostrat$ma_lb[optimized_collections$bin_lb[problem_collections[pz]]];
			optimized_collections$ma_ub[problem_collections[pz]] <- finest_chronostrat$ma_ub[optimized_collections$bin_lb[problem_collections[pz]]];
			}
		} else	{
		optimized_collections$ma_lb[problem_collections[pz]] <- optimized_collections$max_ma[problem_collections[pz]];
		optimized_collections$ma_ub[problem_collections[pz]] <- optimized_collections$min_ma[problem_collections[pz]];				
		optimized_collections$interval_lb[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=optimized_collections$ma_lb[problem_collections[pz]])];
		optimized_collections$interval_ub[problem_collections[pz]] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>optimized_collections$ma_ub[problem_collections[pz]])];
		}
	}

paleodb_finds <- control_occurrences;
finest_chronostrat$ma_lb <- temporal_precision*round(finest_chronostrat$ma_lb/temporal_precision,0)
finest_chronostrat$ma_ub <- temporal_precision*round(finest_chronostrat$ma_ub/temporal_precision,0)
paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections=optimized_collections,uniform_time_scale = finest_chronostrat);
print(paste("Saving",paste(analysis_name,"_Refined_Collections.csv",sep=""),"..."));
#print(taxon_level);
write.csv(paleodb_collections,paste(local_directory,analysis_name,"_Refined_Collections.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
#print(paste("Saving",paste(analysis_name,"_Plus_Control_Finds.csv",sep=""),"..."));
#write.csv(paleodb_finds,paste(local_directory,analysis_name,"_Plus_Control_Finds.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

#### PART 5: SUMMARIZE STRATIGRAPHIC RANGES #####
if (sum(!unique(this_taxon_rank) %in% c("species","subspecies"))>0)	{
	pbdb_finds <- paleodb_finds;
	noccr <- nrow(pbdb_finds)
	for (tx in 1:length(otu_names))	{
		if (this_taxon_rank[tx]=="genus" || this_taxon_rank[tx]=="subgenus")	{
#			print(c(taxon_level,this_taxon_rank[tx]));
			taxon_occ_nos <- (1:noccr)[pbdb_finds$genus==otu_names[tx]]
#			pbdb_finds$accepted_name[taxon_occ_nos] <- otu_names[tx];
			}
		}
#	old_paleodb_finds <- paleodb_finds;
	xxx <- accersi_stratigraphic_information_for_Rev_Bayes(taxa=otu_names,paleodb_finds = pbdb_finds,paleodb_collections = paleodb_collections,hierarchical_chronostrat = finest_chronostrat,taxon_rank = taxon_level,faux_recent = F);
#	paleodb_finds <- old_paleodb_finds;
	} else	{
	xxx <- accersi_stratigraphic_information_for_Rev_Bayes(taxa=otu_names,paleodb_finds = paleodb_finds,paleodb_collections = paleodb_collections,hierarchical_chronostrat = finest_chronostrat,taxon_rank = taxon_level,faux_recent = F);
	}
taxon_field <- this_taxon_rank;
taxon_field[taxon_field %in% c("species","subspecies")] <- "accepted_name";
this_taxon_field <- match(taxon_field,colnames(paleodb_finds));
for (tx in 1:length(otu_names))	{
	relv_collections <- match(paleodb_finds$collection_no[paleodb_finds[,this_taxon_field[tx]]==otu_names[tx]],paleodb_collections$collection_no);
	occurrence_ages_sites <- data.frame(min=as.numeric(paleodb_collections$ma_ub[relv_collections]),max=as.numeric(paleodb_collections$ma_lb[relv_collections]),stringsAsFactors = hell_no);
	occurrence_ages_sites <- occurrence_ages_sites[order(-occurrence_ages_sites$min,-occurrence_ages_sites$max),];
	new_name <- c();
	for (nn in 1:nrow(occurrence_ages_sites))	{
		if (nrow(occurrence_ages_sites)>1)	{
			new_name <- c(new_name,paste(otu_names[tx],nn,sep=""));
			} else	{
			new_name <- otu_names[tx];
			}
		}
	new_site_occurrences <- data.frame(taxon=as.character(new_name),min=as.numeric(occurrence_ages_sites$min),max=as.numeric(occurrence_ages_sites$max),stringsAsFactors = hell_no);
	if (tx==1)	{
		otu_site_occurrences <- new_site_occurrences;
		} else	{
		otu_site_occurrences <- rbind(otu_site_occurrences,new_site_occurrences);
		}

	occurrence_ages_rocks <- unique(data.frame(min=as.numeric(paleodb_collections$ma_ub[relv_collections]),max=as.numeric(paleodb_collections$ma_lb[relv_collections]),rock_no=as.numeric(paleodb_collections$rock_no_sr[relv_collections]),stringsAsFactors = hell_no));
	occurrence_ages_rocks <- occurrence_ages_rocks[order(-occurrence_ages_rocks$min,-occurrence_ages_rocks$max),];
	new_name <- c();
	for (nn in 1:nrow(occurrence_ages_rocks))	{
		if (nrow(occurrence_ages_rocks)>1)	{
			new_name <- c(new_name,paste(otu_names[tx],nn,sep=""));
			} else	{
			new_name <- otu_names[tx];
			}
		}
	new_rock_occurrences <- data.frame(taxon=as.character(new_name),min=as.numeric(occurrence_ages_rocks$min),max=as.numeric(occurrence_ages_rocks$max),stringsAsFactors = hell_no);
	if (tx==1)	{
		otu_rock_occurrences <- new_rock_occurrences;
		} else	{
		otu_rock_occurrences <- rbind(otu_rock_occurrences,new_rock_occurrences);
		}
	}
fossil_summaries <- xxx$fossil_information_detailed;
fossil_summaries$taxon <- names(otu_names);
#### Wrap it up ####
output <- list(paleodb_finds,paleodb_collections,fossil_summaries,otu_site_occurrences,otu_rock_occurrences,finest_chronostrat,missing_taxa);
names(output) <- c("occurrences","collections","fossil_summaries","otu_site_occurrences","otu_rock_occurrences","relv_time_scale","milk_carton_taxa");

return(output);
}

# organize paleodb data for basic diversification analyses.
accersi_data_for_control_groups_to_seed_FBD_analyses <- function(control_taxon,onset="Cambrian",end="Holocene",basic_environments="terr,marine,unknown",paleogeography="scotese",species_only=T)	{
control_finds <- accersi_occurrence_data(taxa=control_taxon,onset=onset,end=end,basic_environments=basic_environments,species_only=species_only,clean_entered_taxa=T,
					directory="",save_files=F,output_type=".csv");
control_collections <- accersi_collection_data(taxa=control_taxon,onset,end,basic_environments,paleogeography=paleogeography,standardize_members=F,
					directory="",save_files=F,species_only=F,output_type=".csv");

control_collections <- control_collections[control_collections$collection_no %in% unique(control_finds$collection_no),];

output <- list(control_collections,control_finds);
names(output) <- c("control_collections","control_occurrences");
return(output);
}

## Convert PaleoDB data to FBD information ##
# create basic summaries of first and last appearances that RevBayes needs for FBD
# modified 2020-02-24 for more sensible output order
# modified 2020-04-08 for more sensible output order
# modified 2020-04-12 for more sensible output order
accersi_stratigraphic_information_for_Rev_Bayes <- function(taxa,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank,sampling_unit="collection",faux_recent=T,lump_cooccr=T,constrain=F,end_FBD="",temporal_precision=0.1)	{
notu <- length(taxa);
n_intervals <- max(hierarchical_chronostrat$bin_last);
if (taxon_rank=="species")	{
	taxon_rank <- "accepted_name";
	} else if (taxon_rank=="subgenus")	{
	genus_name <- taxa;
	broken_up_names <- sapply(genus_name,divido_subgenus_names_from_genus_names);
	taxa[broken_up_names[2,]!=""] <- broken_up_names[2,broken_up_names[2,]!=""];
	}
if (sampling_unit=="collections" || sampling_unit=="localities" || sampling_unit=="locality" || sampling_unit=="site" || sampling_unit=="sites")
	sampling_unit <- "collection";
taxon_col <- match(taxon_rank,colnames(paleodb_finds));	# different taxonomic ranks are in different columns
if (constrain)	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)]);
	latest_bin <- max(hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)]);
	} else	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first);
	latest_bin <- max(hierarchical_chronostrat$bin_last);
	}
if (end_FBD!="")	{
	latest_bin <- min(latest_bin,hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)]);
	}

if (is.null(paleodb_finds$interval_lb))	{
	paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	}
if (is.null(paleodb_finds$bin_lb))	{
	paleodb_finds$bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_finds$interval_lb,hierarchical_chronostrat$interval)];
	paleodb_finds$bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_finds$interval_ub,hierarchical_chronostrat$interval)];
	}
if (is.null(paleodb_finds$ma_lb))	{
	paleodb_finds$ma_lb <- as.numeric(paleodb_collections$ma_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	paleodb_finds$ma_ub <- as.numeric(paleodb_collections$ma_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	}

fossil_intervals <- data.frame(taxon=as.character(rep("",notu)),
							   max=as.numeric(rep(0,notu)),
							   min=as.numeric(rep(0,notu)),
							   stringsAsFactors=hell_no);
fossil_intervals_fuzzy <- data.frame(taxon=as.character(rep("",notu)),
	earliest_poss_fa=as.numeric(rep(0,notu)),latest_poss_fa=as.numeric(rep(0,notu)),
	earliest_poss_la=as.numeric(rep(0,notu)),latest_poss_la=as.numeric(rep(0,notu)),
	stringsAsFactors=hell_no);
fossil_information_detailed <- data.frame(taxon=as.character(rep("",notu)),
	earliest_poss_fa=as.numeric(rep(0,notu)),latest_poss_fa=as.numeric(rep(0,notu)),
	earliest_poss_la=as.numeric(rep(0,notu)),latest_poss_la=as.numeric(rep(0,notu)),
	total_finds=as.numeric(rep(0,notu)),finds_per_bin=as.character(rep("",notu)),stringsAsFactors=hell_no);
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
tx <- 0;
#for (tx in 1:notu)	{
while (tx < notu)	{
	tx <- tx+1;
	fossil_information_detailed$taxon[tx] <- fossil_intervals_fuzzy$taxon[tx] <- fossil_intervals$taxon[tx] <- taxa[tx];
	# get PaleoDB collection numbers
	#taxon_colls <- as.numeric(paleodb_finds$collection_no[paleodb_finds[,taxon_col]==taxa[tx]]);
	taxon_finds <- subset(paleodb_finds,paleodb_finds[,taxon_col]==taxa[tx]);
	taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus, uncertain species");
	
	if (nrow(taxon_finds)==0)	{
		poss_finds <- which(paleodb_finds==taxa[tx],arr.ind = T);
		if (length(poss_finds)>0)	{
			colls_w_finds <- poss_finds[,1];
			taxon_finds <- paleodb_finds[colls_w_finds,];
			}
		}
	if (taxon_rank=="genus" || taxon_rank=="subgenus")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus");
		} else if (taxon_rank=="species" || taxon_rank=="subspecies")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain species");
		}
#	taxon_finds[order(taxon_finds$ma_lb,decreasing = T),];
	taxon_finds_set <- subset(taxon_finds,taxon_finds$interval_lb==taxon_finds$interval_ub);
	taxon_finds_fzy <- subset(taxon_finds,taxon_finds$interval_lb!=taxon_finds$interval_ub);
	
	# unfixed but definitely older collections
	if (nrow(taxon_finds_set)>0)	{
		if (nrow(taxon_finds_fzy)>0 && (min(taxon_finds_fzy$bin_ub) < min(taxon_finds_set$bin_lb)))	{
	#		print(c(tx,"l"));
			taxon_finds_older <- subset(taxon_finds_fzy,taxon_finds_fzy$bin_ub < min(taxon_finds_set$bin_lb));
			taxon_finds_older$bin_lb <- taxon_finds_older$bin_ub;
			taxon_finds_older$interval_lb <- taxon_finds_older$interval_ub;
			taxon_finds_older$ma_lb <- hierarchical_chronostrat$ma_lb[match(taxon_finds_older$interval_ub,hierarchical_chronostrat$interval)];
			taxon_finds_set <- rbind(taxon_finds_set,taxon_finds_older);
			taxon_finds_set <- taxon_finds_set[order(taxon_finds_set$collection_no),];
			}
		# unfixed but definitely younger collections
		if (nrow(taxon_finds_fzy)>0 && (max(taxon_finds_fzy$bin_lb) > max(taxon_finds_set$bin_ub)))	{
	#		print(c(tx,"u"));
			taxon_finds_younger <- subset(taxon_finds_fzy,taxon_finds_fzy$bin_lb > max(taxon_finds_set$bin_ub));
			taxon_finds_younger$bin_ub <- taxon_finds_younger$bin_lb;
			taxon_finds_younger$interval_ub <- taxon_finds_younger$interval_lb;
			taxon_finds_younger$ma_ub <- hierarchical_chronostrat$ma_ub[match(taxon_finds_younger$interval_lb,hierarchical_chronostrat$interval)];
			taxon_finds_set <- rbind(taxon_finds_set,taxon_finds_younger);
			taxon_finds_set <- taxon_finds_set[order(taxon_finds_set$collection_no),];
			}
		if (lump_cooccr || sampling_unit=="rock" || taxon_rank=="species")	{
			taxon_coll_nos <- sort(unique(taxon_finds_set$collection_no));
			} else	{
			taxon_coll_nos <- sort(taxon_finds_set$collection_no);
			}
		fossil_information_detailed$earliest_poss_fa[tx] <- fossil_intervals_fuzzy$earliest_poss_fa[tx] <- fossil_intervals$max[tx] <- max(taxon_finds_set$ma_lb);
		fossil_information_detailed$latest_poss_fa[tx] <- fossil_intervals_fuzzy$latest_poss_fa[tx] <- max(taxon_finds_set$ma_ub);
		fossil_information_detailed$earliest_poss_la[tx] <- fossil_intervals_fuzzy$earliest_poss_la[tx] <- min(taxon_finds_set$ma_lb);
		fossil_information_detailed$latest_poss_la[tx] <- fossil_intervals_fuzzy$latest_poss_la[tx] <- fossil_intervals$min[tx] <- min(taxon_finds_set$ma_ub);
		} else if (nrow(taxon_finds_fzy)>0)	{
		fossil_information_detailed$earliest_poss_fa[tx] <- fossil_intervals_fuzzy$earliest_poss_fa[tx] <- fossil_intervals$max[tx] <- max(taxon_finds$ma_lb);
		fossil_information_detailed$latest_poss_fa[tx] <- fossil_intervals_fuzzy$latest_poss_fa[tx] <- max(taxon_finds$ma_ub);
		fossil_information_detailed$earliest_poss_la[tx] <- fossil_intervals_fuzzy$earliest_poss_la[tx] <- min(taxon_finds$ma_lb);
		fossil_information_detailed$latest_poss_la[tx] <- fossil_intervals_fuzzy$latest_poss_la[tx] <- fossil_intervals$min[tx] <- min(taxon_finds$ma_ub);
		}
	
	if (nrow(taxon_finds)>0)	{
		coll_no <- match(taxon_finds$collection_no,paleodb_collections$collection_no);
		if (lump_cooccr || sampling_unit=="rock")
			coll_no <- unique(coll_no);
		fossil_information_detailed$total_finds[tx] <- length(coll_no);
		taxon_collections <- paleodb_collections[coll_no,];
		if (sampling_unit=="collection")	{
#			early_interval <- as.character(paleodb_collections$interval_lb[coll_no]);
#			late_interval <- as.character(paleodb_collections$interval_ub[coll_no]);
			f_p_b <- tally_collections_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=constrain,temporal_precision=temporal_precision);
			} else	{
#			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat,constrain=constrain,temporal_precision=temporal_precision);
			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat);
			f_p_b <- tally_rock_units_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
			}
		names(f_p_b) <- finest_chronostrat$interval[earliest_bin:latest_bin];
		fossil_information_detailed$finds_per_bin[tx] <- paste(sprintf("%5.2f", f_p_b[earliest_bin:latest_bin]),collapse=";");
		} else	{
		f_p_b <- rep(0,1+(latest_bin-earliest_bin));
		fossil_information_detailed$finds_per_bin[tx] <- paste(sprintf("%5.2f", f_p_b),collapse=";");
		}
#	print(paleodb_collections$ma_lb[123])
	}
#fossil_intervals[31,]
#fossil_intervals_fuzzy[31,]
if (faux_recent)	{
	ingroup_ranges <- sepkoskify_paleodb_data(paleodb_finds,taxon_names=taxa,interval_names=finest_chronostrat$interval);
	if (end_FBD=="")	{
		youngest <- min(ingroup_ranges$ma_min);
		end_FBD <- rebin_collection_with_time_scale(age=youngest,onset_or_end = onset,fine_time_scale = finest_chronostrat);
		faux_recent <- min(fossil_intervals$min);
		fossil_intervals$min[fossil_intervals$max<youngest] <- youngest-0.01;
		fossil_intervals$min <- round(fossil_intervals$min-faux_recent,3);
		fossil_intervals$max <- round(fossil_intervals$max-faux_recent,3);

		fossil_intervals_fuzzy$latest_poss_la[fossil_intervals_fuzzy$latest_poss_la<youngest] <- youngest-0.01;
		fossil_intervals_fuzzy$earliest_poss_la[fossil_intervals_fuzzy$earliest_poss_la<youngest] <- youngest-0.01;
		fossil_intervals_fuzzy$latest_poss_la <- fossil_intervals_fuzzy$latest_poss_la-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_la <- fossil_intervals_fuzzy$earliest_poss_la-faux_recent;
		fossil_intervals_fuzzy$latest_poss_fa <- fossil_intervals_fuzzy$latest_poss_fa-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_fa <- fossil_intervals_fuzzy$earliest_poss_fa-faux_recent;
		} else	{
		if (min(fossil_intervals_fuzzy$latest_poss_fa)<finest_chronostrat$ma_lb[latest_bin])	{
			latest_bin <- 1+finest_chronostrat$bin_first[match(rebin_collection_with_time_scale(age=min(fossil_intervals_fuzzy$latest_poss_fa),onset_or_end = "onset",fine_time_scale = finest_chronostrat),finest_chronostrat$interval)];
			if (latest_bin>max(finest_chronostrat$bin_first))	{
				latest_bin <- max(finest_chronostrat$bin_first);
				end_FBD <- finest_chronostrat$interval[latest_bin]; 
				} else {
				end_FBD <- finest_chronostrat$interval[latest_bin]; 
				}
			}
		faux_recent <- youngest <- finest_chronostrat$ma_lb[match(end_FBD,finest_chronostrat$interval)];
		fossil_intervals$min <- fossil_intervals$min-youngest;
		fossil_intervals$min[fossil_intervals$min<0] <- 0;
		fossil_intervals$max <- fossil_intervals$max-youngest;
		fossil_intervals_fuzzy$latest_poss_la[fossil_intervals_fuzzy$latest_poss_la<youngest] <- youngest;
		fossil_intervals_fuzzy$earliest_poss_la[fossil_intervals_fuzzy$earliest_poss_la<youngest] <- youngest;
		fossil_intervals_fuzzy$latest_poss_la <- fossil_intervals_fuzzy$latest_poss_la-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_la <- fossil_intervals_fuzzy$earliest_poss_la-faux_recent;
		fossil_intervals_fuzzy$latest_poss_fa <- fossil_intervals_fuzzy$latest_poss_fa-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_fa <- fossil_intervals_fuzzy$earliest_poss_fa-faux_recent;
		}
	}

fossil_intervals$min <- temporal_precision*round(fossil_intervals$min/temporal_precision,0);
fossil_intervals$max <- temporal_precision*round(fossil_intervals$max/temporal_precision,0);
fossil_information_detailed$latest_poss_la <- fossil_intervals_fuzzy$latest_poss_la <- temporal_precision*round(fossil_intervals_fuzzy$latest_poss_la/temporal_precision,0);
fossil_information_detailed$latest_poss_fa <- fossil_intervals_fuzzy$latest_poss_fa <- temporal_precision*round(fossil_intervals_fuzzy$latest_poss_fa/temporal_precision,0);
fossil_information_detailed$earliest_poss_la <- fossil_intervals_fuzzy$earliest_poss_la <- temporal_precision*round(fossil_intervals_fuzzy$earliest_poss_la/temporal_precision,0);
fossil_information_detailed$earliest_poss_fa <- fossil_intervals_fuzzy$earliest_poss_fa <- temporal_precision*round(fossil_intervals_fuzzy$earliest_poss_fa/temporal_precision,0);

finds_per_bin_mat <- as.numeric(strsplit(fossil_information_detailed$finds_per_bin[1],";")[[1]]);
names(finds_per_bin_mat) <- finest_chronostrat$interval[1:length(finds_per_bin_mat)];
for (tx in 2:notu)	{
	finds_per_bin_mat <- rbind(finds_per_bin_mat,as.numeric(strsplit(fossil_information_detailed$finds_per_bin[tx],";")[[1]]));
	}
keepers <- (1:ncol(finds_per_bin_mat))[colSums(finds_per_bin_mat)>0];
finds_per_bin_mat <- finds_per_bin_mat[,keepers];

finds_per_bin_new <- c();
for (tx in 1:notu)	{
	this_row <- paste(sprintf("%5.2f", as.numeric(finds_per_bin_mat[tx,])),collapse=";");
	this_row <- gsub(" ","",this_row);
	finds_per_bin_new <- rbind(finds_per_bin_new,this_row);
	}
fossil_information_detailed$finds_per_bin <- finds_per_bin_new;
new_header <- paste(colnames(finds_per_bin_mat),collapse=";")
colnames(fossil_information_detailed)[ncol(fossil_information_detailed)] <- new_header;

output <- list(fossil_intervals,fossil_intervals_fuzzy,fossil_information_detailed);
names(output) <- c("fossil_intervals","fossil_intervals_fuzzy","fossil_information_detailed");
return(output);
}

# matrices of finds per bin in different ways
accersi_per_stratigraphic_interval_sampling_information_for_Rev_Bayes <- function(taxa,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank,sampling_unit="collection",lump_cooccr=T,constrain=F,end_FBD="",temporal_precision=0.1)	{
# taxa: vector of taxon names
# paleodb_finds: data.frame of Paleobiology Database finds
# paleodb_collections: data.frame of collections from the PaleoDB corresponding to finds & linked by collection_no
# hierarchical_chronostrat: chronostratigraphic scheme with information aboust subintervals
# taxon_rank: rank (genus, species, etc.) at which analysis is conducted
# sampling_unit: "collection" to count collections occupied, "rock" to count rock-units occupied.
# lump_cooccr: if TRUE, the only one occurrence per colleciton
# constrain: if TRUE, then do not extend possible range of find counts beyond the latest possible FA & earliest possible LA
if (constrain)	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)]);
	latest_bin <- max(hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)]);
	} else	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first);
	latest_bin <- max(hierarchical_chronostrat$bin_last);
	}
if (end_FBD!="")
	latest_bin <- min(latest_bin,hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)]);
notu <- length(taxa);
n_intervals <- max(hierarchical_chronostrat$bin_last);
if (taxon_rank=="species")	{
	taxon_rank <- "accepted_name";
	} else if (taxon_rank=="subgenus")	{
	genus_name <- taxa;
	broken_up_names <- sapply(genus_name,divido_subgenus_names_from_genus_names);
	taxa[broken_up_names[2,]!=""] <- broken_up_names[2,broken_up_names[2,]!=""];
	}
taxon_col <- match(taxon_rank,colnames(paleodb_finds));
if (sampling_unit=="collections" || sampling_unit=="localities" || sampling_unit=="locality" || sampling_unit=="site" || sampling_unit=="sites")
	sampling_unit <- "collection";
finds_per_bin <- array(0,dim=c(notu,n_intervals));
colnames(finds_per_bin) <- hierarchical_chronostrat$interval[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last];
rownames(finds_per_bin) <- taxa;
sampled_in_bin <- definite_finds_per_bin <- finds_per_bin <- data.frame(finds_per_bin);
if (sampling_unit=="collection")	{
	sample_col_coll <- match("collection_no",colnames(paleodb_collections));
	sample_col_occr <- match("collection_no",colnames(paleodb_finds));
	} else	{
	sample_col_coll <- match("rock_no_sr",colnames(paleodb_collections));
	if (is.null(paleodb_finds$rock_no_sr))
		paleodb_finds$rock_no_sr <- paleodb_collections$rock_no_sr[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	sample_col_occr <- match("rock_no_sr",colnames(paleodb_finds));
	}

tx <- 0;
while (tx < notu)	{			# for debugging.....
	tx <- tx+1;					# for debugging.....
#for (tx in 1:notu)	{
	# get PaleoDB collection numbers
#	print(taxa[tx]);
	# note: we want to use collections initially just because some collections within the same rock
	#	unit can have different ages;
	taxon_finds <- subset(paleodb_finds,paleodb_finds[,taxon_col]==taxa[tx]);
	taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus, uncertain species");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus");
		} else if (taxon_level=="species" || taxon_level=="subspecies")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain species");
		}
	if (nrow(taxon_finds)>0)	{
		occupation <- taxon_finds$collection_no;
		if (lump_cooccr || sampling_unit!="collection")
			occupation <- unique(occupation);
#	taxon_colls <- as.numeric(paleodb_finds$collection_no[paleodb_finds[,taxon_col]==taxa[tx]]);
	# get which collection in this data set the PaleoDB colleciton is
		coll_no <- match(occupation,paleodb_collections$collection_no);
		early_interval <- as.character(paleodb_collections$interval_lb[coll_no]);
		late_interval <- as.character(paleodb_collections$interval_ub[coll_no]);
		taxon_collections <- paleodb_collections[coll_no,];
		if (sampling_unit=="collection")	{
			finds_per_bin[tx,] <- tally_collections_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
#			finds_per_bin[tx,] <- tally_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr,constrain);
			definite_finds_per_bin[tx,] <- tally_definite_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
			} else	{
			## START HERE!!!!
			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat=subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last));
			finds_per_bin[tx,] <- tally_rock_units_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
			definite_finds_per_bin[tx,] <- tally_rock_units_definitely_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F);
			}
		sampled_in_bin[tx,] <- tally_presence_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
		}
	# latest possible first appearance would be the earliest late interval
	}
colnames(sampled_in_bin) <- colnames(definite_finds_per_bin) <- colnames(finds_per_bin) <- gsub("\\."," ",colnames(finds_per_bin));
output <- list(finds_per_bin[earliest_bin:latest_bin],definite_finds_per_bin[earliest_bin:latest_bin],sampled_in_bin[earliest_bin:latest_bin]);
names(output) <- c("finds_per_bin","definite_finds_per_bin","sampled_in_bin");
return(output);
}

#### ROUTINES TO GET BASIC SAMPLING & DIVERSIFICATION NUMBERS INFORMATION ####
# routine to test out best uniform, exponential, beta & lognormal distributions
# updated 2020-04-13
accersi_sampling_distributions_for_RevBayes <- function(finds_per_bin,sample_units_per_bin,end_FBD="",update_search=T,minS=4)	{
ok <- 0;
good_bins <- c();
sampled_in_bin <- colSums(round(finds_per_bin+0.001,0)>0);
# get fractional presences for all rocks
all_bin_sampling <- round(sample_units_per_bin,0);
if (end_FBD!="")	{
	last_bin <- match(end_FBD,colnames(finds_per_bin));
	} else	{
	last_bin <- ncol(finds_per_bin);
	}
sampled_in_bin <- colSums(finds_per_bin>=0.5);
interval_richness <- accersi_synoptic_richness_from_sampled_in_bin(finds_per_bin);

classic_richness <- setup_three_timer_analysis(samples_per_interval = finds_per_bin);
synoptic_richness <- classic_richness$sepkoski_richness;

#colMax(finds_per_bin)
bn <- 0;
while (bn < last_bin)	{
	bn <- bn+1;
#for (bn in 1:last_bin)	{
	interval_name <- names(interval_richness)[bn];
	Sb <- synoptic_richness[bn];			# bin richness
	bn_coll <- round(sample_units_per_bin[bn],0);
	bin_occcurences <- sort(round(finds_per_bin[,bn]+0.01,0),decreasing = T);
	if (sum(bin_occcurences>0)>=minS)	{
		if (bn_coll<=max(bin_occcurences))
			bn_coll <- 1+max(bin_occcurences);
		bin_no <- match(interval_name,colnames(finds_per_bin));
		tod <- sort(bin_occcurences,decreasing = T)[1:Sb];
		if (tod[1]==1)
			tod[1] <- 2;
		if (max(tod)>1)	{
			good_bins <- c(good_bins,interval_name);
			if (update_search)
				print(paste("finding sampling distributions for",interval_name));
			if (ok==1)	{
				sampling_uni_bin <- rbind(sampling_uni_bin,data.frame(base::t(optimize_uniform_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_exp_bin <- rbind(sampling_exp_bin,data.frame(base::t(optimize_exponential_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_bta_bin <- rbind(sampling_bta_bin,data.frame(base::t(optimize_beta_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_lgn_bin <- rbind(sampling_lgn_bin,data.frame(base::t(optimize_lognormal_occupancy(finds=tod,ncoll=bn_coll))));
				} else	{
				sampling_uni_bin <- data.frame(base::t(optimize_uniform_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_exp_bin <- data.frame(base::t(optimize_exponential_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_bta_bin <- data.frame(base::t(optimize_beta_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_lgn_bin <- data.frame(base::t(optimize_lognormal_occupancy(finds=tod,ncoll=bn_coll)));
				ok <- 1;
				}
			}
		}
	}
#names(all_bin_sampling) <- finest_chronostrat$interval;
#names(bin_sampling) <- good_bins;
rownames(sampling_uni_bin) <- rownames(sampling_exp_bin) <- rownames(sampling_bta_bin) <- rownames(sampling_lgn_bin) <- good_bins;
#rownames(sampling_uni_bin) <- rownames(sampling_exp_bin) <- rownames(sampling_lgn_bin) <- good_bins;
output <- list(sampling_uni_bin,sampling_exp_bin,sampling_bta_bin,sampling_lgn_bin);
#output <- list(sampling_uni_bin,sampling_exp_bin,sampling_lgn_bin);
names(output) <- c("uniform","exponential","beta","lognormal");
#names(output) <- c("sampling_uni_bin","sampling_exp_bin","sampling_lgn_bin");
return(output);
}

# routine to get the typical per-interval sampling rate given the best fit distributions for individual intervals
accersi_median_per_unit_sampling_probability_from_best_distributions <- function(sampling_distributions)	{

diff_dists <-names(sampling_distributions);
diff_dists <- gsub("best_","",diff_dists);
diff_dists <- gsub("_sampling","",diff_dists);
aiccs <- c();
for (dd in 1:length(diff_dists))
	aiccs <- cbind(aiccs,sampling_distributions[[dd]]$AICc);
rownames(aiccs) <- rownames(sampling_distributions[[1]]);
colnames(aiccs) <- diff_dists;
	
bin_spans <- psi_binA <- psi_binB <- psi_binC <- rho_binA <- rho_binB <- c();
per_unit_sampling <- per_unit_sampling_A <- c();
for (bb in 1:nrow(aiccs))	{
#	bn <- match(rownames(aiccs)[bb],chronostrat$interval);
#	bin_span <- chronostrat$ma_lb[bn]-chronostrat$ma_ub[bn];
#	bin_spans <- c(bin_spans,bin_span);
	this_dist <- match(min(aiccs[bb,]),aiccs[bb,]);
	dist_info <- sampling_distributions[[this_dist]][bb,];
	priors <- exp(-aiccs[bb,]/2)/sum(exp(-aiccs[bb,]/2));
	if (diff_dists[this_dist]=="uniform")	{
		per_unit_sampling <- c(per_unit_sampling,dist_info$scale);
		} else if (diff_dists[this_dist]=="exponential")	{
		best_distribution <- scaled_exponential_distribution(sc=dist_info$scale,decay=dist_info$decay);
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution[1:dist_info$richness]));
		} else if (diff_dists[this_dist]=="beta")	{
		best_distribution <- beta_distribution(shape1 = dist_info$alpha,shape2=dist_info$beta,S=dist_info$richness)
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution));
		} else if (diff_dists[this_dist]=="lognormal")	{
		best_distribution <- scaled_lognormal_distribution(sc=dist_info$scale,mag=dist_info$mag_var,S=dist_info$richness);
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution));
		}
	}
names(per_unit_sampling) <- rownames(sampling_distributions[[1]]);
return(per_unit_sampling);
}

# routine to get the typical per-collection or per-formation sampling rate given the best fit distributions for individual intervals
accersi_median_prob_find_per_unit_finds_given_multiple_distributions <- function(sampling_distributions)	{
diff_dists <-names(sampling_distributions);
diff_dists <- gsub("best_","",diff_dists);
diff_dists <- gsub("_sampling","",diff_dists);
aiccs <- c();
for (dd in 1:length(diff_dists))
	aiccs <- cbind(aiccs,sampling_distributions[[dd]]$AICc);
aiccs <- data.frame(aiccs,stringsAsFactors=hell_no);
per_unit_sampling <- data.frame(array(0,dim=dim(aiccs)),stringsAsFactors=hell_no);
rownames(per_unit_sampling) <- rownames(aiccs) <- rownames(sampling_distributions[[1]]);
colnames(per_unit_sampling) <- colnames(aiccs) <- diff_dists;

for (bb in 1:nrow(aiccs))	{
#	bn <- match(rownames(aiccs)[bb],chronostrat$interval);
#	bin_span <- chronostrat$ma_lb[bn]-chronostrat$ma_ub[bn];
#	bin_spans <- c(bin_spans,bin_span);
	priors <- exp(-aiccs[bb,]/2)/sum(exp(-aiccs[bb,]/2));
	if (!is.na(match("uniform",diff_dists)))	{
		md <- match("uniform",diff_dists);
		per_unit_sampling$uniform[bb] <- sampling_distributions[[md]]$scale[bb];
		}
	if (!is.na(match("exponential",diff_dists)))	{
		md <- match("exponential",diff_dists);
		best_distribution <- scaled_exponential_distribution(sc=sampling_distributions[[md]]$scale[bb],decay=sampling_distributions[[md]]$decay[bb]);
		per_unit_sampling$exponential[bb] <- median(best_distribution[1:sampling_distributions[[md]]$richness[bb]]);
		}
	if (!is.na(match("beta",diff_dists)))	{
		md <- match("beta",diff_dists);
		best_distribution <- beta_distribution(shape1=sampling_distributions[[md]]$alpha[bb],shape2=sampling_distributions[[md]]$beta[bb],S=sampling_distributions[[md]]$richness[bb])
		per_unit_sampling$beta[bb] <- median(best_distribution);
		}
	if (!is.na(match("lognormal",diff_dists)))	{
		md <- match("lognormal",diff_dists);
		best_distribution <- scaled_lognormal_distribution(sc=sampling_distributions[[md]]$scale[bb],mag=sampling_distributions[[md]]$mag_var[bb],S=sampling_distributions[[md]]$richness[bb]);
		per_unit_sampling$lognormal[bb] <- median(best_distribution);
		}
	}
output <- list(per_unit_sampling,aiccs);
names(output) <- c("per_unit_sampling","AICc")
return(output);
}

weighted_median_sampling_probability <- function(sampling_distributions,sampling_opportunities)	{
per_bin_sampling_basics <- accersi_median_prob_find_per_unit_finds_given_multiple_distributions(sampling_distributions);
sampling_priors <- exp(per_bin_sampling_basics$AICc/-2)/rowSums(exp(per_bin_sampling_basics$AICc/-2));
per_unit_sampling <- per_bin_sampling_basics$per_unit_sampling;
per_unit_weighted_prob <- c();
for (bb in 1:nrow(per_unit_sampling))	{
	per_unit_weighted_prob <- c(per_unit_weighted_prob,sum(sampling_priors[bb,]*(1-(1-per_unit_sampling[bb,]))));
	}
names(per_unit_weighted_prob) <- rownames(per_unit_sampling);
return(per_unit_weighted_prob);
}

# expected finds per million years per interval
per_interval_per_ma_psis <- function(sampling_distributions,sample_units_per_bin,bin_spans)	{
per_unit_sampling_probs <- weighted_median_sampling_probability(sampling_distributions,sampling_opportunities=sample_units_per_bin[names(sample_units_per_bin) %in% rownames(sampling_distributions[[1]])]);
bin_psi_pma <- c();
for (bn in 1:length(sample_units_per_bin))	{
	if (!is.na(match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs))))	{
		bb <- match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs));
		bin_psi_pma <- c(bin_psi_pma,probability_to_Poisson_rate(1-(1-per_unit_sampling_probs[bb])^sample_units_per_bin[bn])/bin_spans[bn]);
		} else	{
		bin_psi_pma <- c(bin_psi_pma,probability_to_Poisson_rate(1-median((1-per_unit_sampling_probs)^sample_units_per_bin[bn]))/bin_spans[bn]);
		}
	}
names(bin_psi_pma) <- names(sample_units_per_bin);
return(bin_psi_pma)
}

# expected finds per interval
per_interval_psis <- function(sampling_distributions,sample_units_per_bin)	{
# per_unit_sampling_probs gives the average probability of sampling a taxon per collection or rock unit
per_unit_sampling_probs <- weighted_median_sampling_probability(sampling_distributions,sampling_opportunities=sample_units_per_bin[names(sample_units_per_bin) %in% rownames(sampling_distributions[[1]])]);
bin_psi <- c();
for (bn in 1:length(sample_units_per_bin))	{
	if (!is.na(match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs))))	{
		bb <- match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs));
		bin_psi <- c(bin_psi,probability_to_Poisson_rate(1-(1-per_unit_sampling_probs[bb])^sample_units_per_bin[bn]));
		} else	{
		bin_psi <- c(bin_psi,probability_to_Poisson_rate(1-median((1-per_unit_sampling_probs)^sample_units_per_bin[bn])));
		}
	}
names(bin_psi) <- names(sample_units_per_bin);
return(bin_psi)
}

# first pass estimate of diversification rates given sampling
accersi_initial_diversification_rates <- function(sampled_in_bin,synoptic_richness,psi_bin,chronostrat)	{

gap_taxa <- synoptic_richness-colSums(sampled_in_bin);
if (length(psi_bin)>length(gap_taxa))
	psi_bin <- psi_bin[names(psi_bin) %in% names(gap_taxa)];
chronostrat <- subset(chronostrat,chronostrat$interval %in% names(psi_bin));
bin_spans <- chronostrat$ma_lb-chronostrat$ma_ub;
names(bin_spans) <- chronostrat$interval;
ttl_taxa <- nrow(sampled_in_bin);
per_bin_origination <- per_bin_extinction <- data.frame(best_rate=rep(0,length(psi_bin)),log_likelihood=rep(0,length(psi_bin)));
for (bn in 1:length(psi_bin))	{
	# do origination
	S1 <- sum(sampled_in_bin[,bn]);
	if (bn>2)	{
		pmiss <- dpois(0,psi_bin[bn-1]*bin_spans[bn-1]);
		if (!is.na(pmiss))	{
			two_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn-1]);
			two_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn-1]);
			three_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn-2]);
			three_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn-2]);
			gap_filler <- length(three_timers[!three_timers %in% two_timers]);
			per_bin_origination[bn,] <- accersi_best_diversification_given_sampling(pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=T);
			} else	{
			per_bin_origination[bn,] <- c(0,0);
			}
		}  else	{
		per_bin_origination[bn,] <- c(0,0);
		}
	# do extinction
	if (bn<(ncol(sampled_in_bin)-2))	{
		pmiss <- dpois(0,psi_bin[bn+1]*bin_spans[bn+1]);
		if (!is.na(pmiss))	{
			two_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn+1]);
			two_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn+1]);
			three_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn+2]);
			three_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn+2]);
			gap_filler <- length(three_timers[!three_timers %in% two_timers]);
			per_bin_extinction[bn,] <- accersi_best_diversification_given_sampling(pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=T);
			} else	{
			per_bin_extinction[bn,] <- c(0,0);
			}
		}	else	{
		per_bin_extinction[bn,] <- c(0,0);
		}
	}
rownames(per_bin_origination) <- rownames(per_bin_extinction) <- names(gap_taxa);
per_bin_origination$best_rate <- per_bin_origination$best_rate/bin_spans;
per_bin_extinction$best_rate <- per_bin_extinction$best_rate/bin_spans;
per_bin_origination <- per_bin_origination[per_bin_origination$best_rate>0,];
per_bin_extinction <- per_bin_extinction[per_bin_extinction$best_rate>0,];

bin_spans <- chronostrat$ma_lb-chronostrat$ma_ub;
names(bin_spans) <- chronostrat$interval;
origination <- sum(bin_spans[match(rownames(per_bin_origination),names(bin_spans))]*per_bin_origination$best_rate)/sum(bin_spans[match(rownames(per_bin_origination),names(bin_spans))]);
extinction <- sum(bin_spans[match(rownames(per_bin_extinction),names(bin_spans))]*per_bin_extinction$best_rate)/sum(bin_spans[match(rownames(per_bin_extinction),names(bin_spans))]);
output <- c(origination,extinction);
names(output) <- c("origination","extinction");
return(output);
}

# altered version of Liow & Starrfelt's TRIPPs estimate
modified_TRIPPs <- function(psi,taxon_intervals_finds,taxon_intervals_spans)	{
dp <- c();
for (i in 1:length(taxon_intervals_finds))	{
	if (is.integer(taxon_intervals_finds[i]))	{
		dp <- c(dp,dpois(x=taxon_intervals_finds[i],lambda=psi*taxon_intervals_spans[i]))
		} else	{
		dp1 <- dpois(x=floor(taxon_intervals_finds[i]),lambda=psi*taxon_intervals_spans[i])
		dp2 <- dpois(x=ceiling(taxon_intervals_finds[i]),lambda=psi*taxon_intervals_spans[i])
		wt2 <- (ceiling(taxon_intervals_finds[i])-taxon_intervals_finds[i]);
		wt1 <- 1-wt2;
		dp <- c(dp,(wt1*dp1)+(wt2*dp2));
		}
	}
dp[dp<MINNO] <- MINNO;
lndp <- log(dp);
return(sum(lndp));
}

#### ROUTINES TO GET BASIC TAXONOMIC INFORMATION ####
accersi_updated_taxonomy_for_analyzed_taxa <- function(otu_names,local_directory="",study="")	{
notu <- length(otu_names);
taxon <- otu_names;
print("Getting taxonomic data....")
#dummy <- revelare_taxonomy_for_one_taxon(taxon="Lophospira");

unentered_taxa <- c();
for (n in 1:notu)	{
	this_taxon_info <- revelare_taxonomy_for_one_taxon(taxon=otu_names[n],settle=F);
	if (ncol(this_taxon_info)<ncol(dummy))	{
		unentered_taxa <- c(unentered_taxa,otu_names[n]);
		this_taxon_info <- revelare_taxonomy_for_one_taxon(taxon=otu_names[n],settle=T);
		}	
	if (n==1)	{
#		initial_compendium <- revelare_taxonomy_for_one_taxon(taxon=otu_names[n]);
		initial_compendium <- this_taxon_info;
		} else	{
		initial_compendium <- rbind(initial_compendium,this_taxon_info);
		}
	}
#initial_compendium <- data.frame(base::t(sapply(taxon,revelare_taxonomy_for_one_taxon)),stringsAsFactors=hell_no);
#revelare_taxonomy_for_one_taxon(taxon="Proetus latifrons")
initial_compendium <- evanesco_na_from_matrix(initial_compendium,"");
compendium_headers <- colnames(initial_compendium);
taxon_compendium <- data.frame(array("",dim=dim(initial_compendium)),stringsAsFactors=hell_no);
colnames(taxon_compendium) <- compendium_headers;

for (cn in 1:ncol(initial_compendium))	{
	old_info <- unlist(initial_compendium[,cn]);
	header_words <- strsplit(compendium_headers[cn],"_")[[1]];
	if (header_words[length(header_words)] %in% paleodb_numeric_fields)	{
		old_info[old_info %in% missing_data_assignment] <- 0;
		taxon_compendium[,cn] <- as.numeric(old_info);
		} else	{
		taxon_compendium[,cn] <- as.character(old_info);
		taxon_compendium[is.na(taxon_compendium[,cn]),cn] <- "";
		}
	}
#for (rn in 1:nrow(initial_compendium))
#	for (cn in 1:ncol(initial_compendium))
#		taxon_compendium[rn,cn] <- as.character(initial_compendium[rn,cn]);
#taxon_compendium <- evanesco_na_from_matrix(taxon_compendium,replacement="");
#colnames(taxon_compendium)
if (length(unentered_taxa)>0)	{
	missing_taxa_info <- taxon_compendium[match(unentered_taxa,otu_names),];
	missing_taxa_info$taxon_name <- unentered_taxa;
	} else	{
	missing_taxa_info <- taxon_compendium[!1:length(otu_names),];
	}
#subset(taxon_compendium,taxon_compendium$accepted_name=="?");
found_taxa_info <- subset(taxon_compendium,taxon_compendium$n_occs>0);

if (nrow(missing_taxa_info)>0)	{
	print("Not all of your taxa are entered into the PaleoDB: you should fix that!");
	if (study=="")
		study <- paste(otu_names[1],"Clade");
	file_name <- paste(local_directory,study,"_Unentered_Taxa.csv",sep="");
	write.csv(missing_taxa_info,file=file_name,row.names = FALSE);
	}
taxon_compendium$n_occs <- as.numeric(taxon_compendium$n_occs);
taxon_compendium$orig_no <- as.numeric(taxon_compendium$orig_no);
taxon_compendium$taxon_no <- as.numeric(taxon_compendium$taxon_no);
taxon_compendium$accepted_no <- as.numeric(taxon_compendium$accepted_no);
taxon_compendium$parent_no <- as.numeric(taxon_compendium$parent_no);
taxon_compendium$immpar_no <- as.numeric(taxon_compendium$immpar_no);
taxon_compendium$phylum_no <- as.numeric(taxon_compendium$phylum_no);
taxon_compendium$order_no <- as.numeric(taxon_compendium$order_no);
taxon_compendium$family_no <- as.numeric(taxon_compendium$family_no);
taxon_compendium$genus_no <- as.numeric(taxon_compendium$genus_no);
taxon_compendium$subgenus_no <- as.numeric(taxon_compendium$subgenus_no);
taxon_compendium <- evanesco_na_from_matrix(taxon_compendium,replacement=0);

#taxon_list <- simplify2array(found_taxa_info$accepted_name);
#output <- list(found_taxa_info,missing_taxa_info);
return(taxon_compendium);
}

##### CAL3 #####
### Routines for Phylogeny Prior Probabilities
prob_sampling_clade_bapst <- function(p,q,r)	{
# modified from Bapst (2013) equation 2
# p: origination (lambda in Bapst 2013)
# q: extinction  (mu in Bapst 2013)
# r: sampling (psi in Bapst 2013)
s <- 1		# richness of clade (K in Bapst's paper)
pmiss <- mxp <- 0
num <- den <- 1
while ((mxp/100000)<(num/den))	{
	num <- (p^(s-1))*(q^s)*choose((2*s)-2,s-1)
	den <- s*((p+q+r)^((2*s)-1))
	pmiss <- pmiss+(num/den)
	if (mxp<(num/den))	mxp <- (num/den)
#	print(c(s,num,den,num/den,pmiss))
	s <- s+1
	}
return(1-pmiss)
}

##### ROUTINES TO WRITE CHARACTER MATRIX INFORMATION IN NEXUS FILE #######
scribere_rev_bayes_nexus_file_from_character_matrix <- function(ch_matrix,state_symbols,new_file_name,UNKNOWN=-11,INAP=-22)	{
# ch_matrix: character matrix being printed to file.
# new_file_name: name of nexus file, including the directory to which it should be printed.
# no_states: number of states; note that if only a portion of the taxa are used, then sometimes no otus have some states
if (!is.matrix(ch_matrix))	{
	ch_matrix <- data.frame(ch=ch_matrix);
	}
notu <- nrow(ch_matrix);
taxon_names <- rownames(ch_matrix);
nchars <- ncol(ch_matrix);

nexus_file_content <- c();
nexus_file_content <- rbind("#NEXUS","","BEGIN DATA;")
nexus_file_content <- rbind(nexus_file_content,paste("	DIMENSIONS  NTAX=",notu," NCHAR=",nchars,";",sep=""));

nexus_file_content <- rbind(nexus_file_content,paste("	FORMAT DATATYPE = STANDARD RESPECTCASE GAP = - MISSING = ? SYMBOLS = \"",paste(state_symbols,collapse=" "),"\";"));
nexus_file_content <- rbind(nexus_file_content,"	MATRIX");

string_to_count <- taxon_names;
name_lengths <- sapply(string_to_count,count_characters_in_string);
max_name_length <- max(name_lengths);
need_quotes <- c(".","(",")","[","]");
for (nn in 1:notu)	{
	test_name <- strsplit(taxon_names[nn],split="",fixed=TRUE)[[1]]
	if (sum(test_name %in% need_quotes)==0)	{
		taxon <- gsub(" ","_",taxon_names[nn]);
		} else	{
		taxon <- paste("\"",taxon_names[nn],"\"",sep="");
		name_lengths[nn] <- name_lengths[nn]+2;
		}
	this_line <- paste("\t",taxon,paste(rep(" ",(5+(max_name_length-name_lengths[nn]))),collapse=""),sep="");
	otu_code <- c();
	for (ch in 1:nchars)	{
		if (ch_matrix[nn,ch]>=0 && ch_matrix[nn,ch]<=9)	{
			otu_code <- paste(otu_code,ch_matrix[nn,ch],sep="");
			} else if (ch_matrix[nn,ch]>9)	{
			otu_code <- paste(otu_code,all_states[1+ch_matrix[nn,ch]],sep="");	# note: we need +1 because of state 0
			} else if (ch_matrix[nn,ch]==UNKNOWN)	{
			otu_code <- paste(otu_code,"?",sep="");
			} else if (ch_matrix[nn,ch]==INAP)	{
			otu_code <- paste(otu_code,"-",sep="");
			} else if (ch_matrix[nn,ch]<0)	{
			polystates <- strsplit(as.character(ch_matrix[nn,ch]),split="",fixed=TRUE)[[1]];
			polystates <- as.numeric(polystates[polystates!="-"]);
			otu_code <- paste(otu_code,ravel_polymorph_for_file(polystates),sep="");
			}
		}
	nexus_file_content <- rbind(nexus_file_content,paste(this_line,otu_code,sep=""));
	}

nexus_file_content <- rbind(nexus_file_content,";");
nexus_file_content <- rbind(nexus_file_content,"END;");
#nexus_file_content <- rbind(nexus_file_content,"begin mrbayes;");
#nexus_file_content <- rbind(nexus_file_content,"	set autoclose=yes nowarn=yes;");
#nexus_file_content <- rbind(nexus_file_content,"	lset nst=6 rates=invgamma;");
#nexus_file_content <- rbind(nexus_file_content,"	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); ");
#nexus_file_content <- rbind(nexus_file_content,"	prset applyto=(all) ratepr=variable;");
#nexus_file_content <- rbind(nexus_file_content,"	mcmcp ngen= 100000000 relburnin=yes burninfrac=0.25 printfreq=10000  samplefreq=10000 nchains=4 savebrlens=yes;");
#nexus_file_content <- rbind(nexus_file_content,"	mcmc;");
#nexus_file_content <- rbind(nexus_file_content,"	sumt;");
#nexus_file_content <- rbind(nexus_file_content,"end;");
write(nexus_file_content,file=new_file_name);
}

##### ELEMENTARY DIVERGENCE TIME ESTIMATION	#####
simple_likelihood_divergence <- function(bound_1,bound_2,psi,cutoff=exp(-4),precision=0.01)	{
# bound_1: oldest taxon;
# bound_2: 2nd oldest taxon;
# psi: probability of sampling a taxon per time unit 1;
# cutoff: area under likelihood curve at which we fix the divergence (null: 1 unit support)
# precision: temporal precision for estimating divergence, with expected finds = psi*precision
bound_1 <- abs(bound_1);
bound_2 <- abs(bound_2);
min_divergence_date <- divergence_date <- max(bound_1,bound_2);
gap_ldf <- c();				# likelihood density function for gap
pgap <- 1;
while (pgap > 10^-10)	{
	pgap <- dpois(0,psi*((divergence_date-bound_1)+(divergence_date-bound_2)));
	gap_ldf <- c(gap_ldf,pgap);
	divergence_date <- divergence_date+precision;
	}
return(min_divergence_date+(precision*(1+sum(gap_ldf>cutoff))));
#return(min_divergence_date+(precision*(1+sum(cumsum(gap_ldf/sum(gap_ldf))<0.5))));
#return(min_divergence_date+(precision*(1+sum(cumsum(gap_pdf)<0.5))));
}

simple_probability_divergence <- function(bound_1,bound_2,phi,psi,cutoff=exp(-8),precision=0.01)	{
# bound_1: oldest taxon;
# bound_2: 2nd oldest taxon;
# phi: probability of sampled sister taxon (lineage or clade) arising;
# psi: probability of sampling a taxon per time unit 1;
# cutoff: area under likelihood curve at which we fix the divergence.  (Default: logP = -1)
# precision: temporal precision for estimating divergence, with expected finds = psi*precision
bound_1 <- abs(bound_1);
bound_2 <- abs(bound_2);
min_divergence_date <- divergence_date <- max(bound_1,bound_2);
gap_pdf <- c();
pgap <- lgap <- 1;
while (exp(log(pgap)+log(lgap)) > 10^-10)	{
#	gap_pdf <- c(gap_pdf,dpois(0,phi*((divergence_date-bound_1)+(divergence_date-bound_2))));
	pgap <- (1-phi)^((divergence_date-bound_1)+(divergence_date-bound_2));
	lgap <- dpois(0,psi*((divergence_date-bound_1)+(divergence_date-bound_2)));
	gap_pdf <- c(gap_pdf,exp(log(pgap)+log(lgap)));
	divergence_date <- divergence_date+precision;
	}
return(min_divergence_date+(precision*(1+(sum(gap_pdf>cutoff)))));
}

##### ROUTINES FOR A PRIOR EVALUATION OF POSSIBLE CHARACTER MODELS #####
revelare_possible_driven_trends <- function(analysis_name,nexus_file_name=NULL,nexus_file_directory=NULL,outgroup="",strat_data_file=NULL,UNKNOWN=-11,INAP=-22)	{
if (is.null(nexus_file_name))	{
	nexus_file_name <- file.choose();
	basic_data <- accersi_data_from_nexus_file(nexus_file_name);
	} else	{
	new_file_name <- paste(nexus_file_directory,nexus_file_name,sep="");
	basic_data <- accersi_data_from_nexus_file(nexus_file_name=new_file_name);
	}
ch_matrix <- basic_data$Matrix;
otus <- basic_data$OTUs;
if (outgroup=="")
	outgroup <- otus[as.numeric(basic_data$Outgroup)];

ch_matrix_red <- ch_matrix[!rownames(ch_matrix) %in% outgroup,];

if (is.null(strat_data_file))
	strat_data_file <- file.choose();
finds_per_bin_data <- read.csv(strat_data_file,header=T,stringsAsFactors=hell_no);
if (is.na(match("Species",colnames(finds_per_bin_data))))	{
	ccns <- colnames(finds_per_bin_data);
	ccns[1] <- "Species";
	colnames(finds_per_bin_data) <- ccns;
	}
finds_per_bin_data <- subset(finds_per_bin_data,!finds_per_bin_data$Species %in% outgroup);
rounded_finds_per_bin <- finds_per_bin_data;
rownames(rounded_finds_per_bin) <- finds_per_bin_data$Species;
rounded_finds_per_bin$Species <- NULL;
notu <- nrow(rounded_finds_per_bin);
for (i in 1:notu)
	rounded_finds_per_bin[i,] <- round(cumsum(as.numeric(rounded_finds_per_bin[i,]))+0.001,0);
bin_names <- colnames(rounded_finds_per_bin);
bin_names <- bin_names[colSums(rounded_finds_per_bin)>0]
rounded_finds_per_bin <- rounded_finds_per_bin[,colnames(rounded_finds_per_bin) %in% bin_names];
nbins <- length(bin_names);
fa_bin_no <- fa_bin <- c();
for (i in 1:notu)
	fa_bin_no <- c(fa_bin_no,min((1:nbins)[rounded_finds_per_bin[i,]>0]));
fa_bin <- bin_names[fa_bin_no];
#fa_bin_no[fa_bin_no==5] <- 4;
nchars <- ncol(ch_matrix_red);
nstates <- basic_data$States;
all_char_time_state_matrix <- assessed_characters <- significant_cases <- kendall_pval <- kendall_tau <- c();

for (nch in 1:nchars)	{
	char_states <- ch_matrix_red[,nch];
	char_binning <- fa_bin_no*((char_states != UNKNOWN)*(char_states != INAP));	# first appearances of states in different taxa
	char_binning <- char_binning[char_binning!=0];								# cull out species with unknown or inapplicable
	unique_bins <- sort(unique(char_binning));
	char_states <- char_states[char_states!=UNKNOWN];							# cull out species with unknown
	char_states <- char_states[char_states!=INAP];								# cull out species with inapplicable

	if (sum(char_states<0)>0)	{
		poly_cases <- (1:length(char_states))[char_states<0];
		poly_fas <- char_binning[char_states<0];
		for (pf in 1:length(poly_fas))	{
			polystates <- unravel_polymorph(char_states[poly_cases[pf]]);
			char_states[poly_cases[pf]] <- polystates[1];
			char_states <- c(char_states,polystates[2]);
			char_binning <- c(char_binning,poly_fas[pf]);
			}
		}
	names(char_states) <- names(char_binning);
	state_nos <- sort(unique(char_states[char_states>=0]));
	this_chars_fas <- fa_bin_no[match(names(char_states),rownames(ch_matrix_red))];
	state_fas <- c();
	for (nst in 1:nstates[nch])
		state_fas <- c(state_fas,min(char_binning[char_states==state_nos[nst]]));
	state_nos <- state_nos[order(state_fas)];
	
	if (length(unique(char_states))==1)	{
		kendall_pval <- c(kendall_pval,1.0);
		kendall_tau <- c(kendall_tau,0.00);
		} else	{
		kendall_pval <- c(kendall_pval,cor.test(x=char_binning,y=char_states,method="kendall")$p.value);
		kendall_tau <- c(kendall_tau,cor.test(x=char_binning,y=char_states,method="kendall")$estimate);
		}
	
	time_state_matrix <- array(0,dim=c(nbins,max(nstates)));
	for (nn in 1:length(char_states))	{
		if (char_states[nn]<0)	{
			polystates <- unravel_polymorph(poly=char_states[nn]);
			time_state_matrix[this_chars_fas[nn],match(polystates,state_nos)] <- time_state_matrix[this_chars_fas[nn],match(polystates,state_nos)]+(1/length(polystates));
			} else	{
			time_state_matrix[this_chars_fas[nn],match(char_states[nn],state_nos)] <- time_state_matrix[this_chars_fas[nn],match(char_states[nn],state_nos)]+1;
			}
		}
	rownames(time_state_matrix) <- bin_names;
	time_state_matrix <- subset(time_state_matrix,rowSums(time_state_matrix)>0);
	
	all_char_time_state_matrix <- rbind(all_char_time_state_matrix,cbind(rep(nch,nrow(time_state_matrix)),rownames(time_state_matrix),time_state_matrix));

	bin_cases <- rowSums(time_state_matrix);
	mann_whitney_results <- c();
	comparisons <- 0;
	for (bn in 1:(nrow(time_state_matrix)-1))	{
		if (bin_cases[bn] > 0 && bin_cases[bn+1] > 0)	{
			data_v <- category <- c();
			for (bnb in 0:1)	{
				for (bnst in 1:nstates[nch])	{
					data_v <- c(data_v,rep(state_nos[bnst],round(time_state_matrix[bn+bnb,bnst])));
					category <- c(category,rep(bnb,round(time_state_matrix[bn+bnb,bnst])));
					}
				}
			comparisons <- comparisons + 1;
			mann_whitney_results <- rbind(mann_whitney_results,mann_whitney(data_v,category));
			}
		}
	if (comparisons>0)	{
		assessed_characters <- c(assessed_characters,nch);
		significant_cases <- c(significant_cases,sum(mann_whitney_results$pval<=0.10));
		} else	{
		significant_cases <- c(significant_cases,0);
		}
	}
kendall_summary <- data.frame(nch=as.numeric(1:nchars),tau=as.numeric(kendall_tau),pval=as.numeric(kendall_pval),stringsAsFactors=hell_no);
kendall_summary <- kendall_summary[order(kendall_summary$pval),];
trend_summary <- data.frame(nch=as.numeric(1:nchars),states=as.numeric(nstates),tau=as.numeric(kendall_tau),pval=as.numeric(kendall_pval),MW_sign=as.numeric(significant_cases),stringsAsFactors=hell_no);
#trend_summary[order(trend_summary$states,trend_summary$pval),];
state_distributions <- data.frame(all_char_time_state_matrix,stringsAsFactors=hell_no);
colnames(state_distributions) <- c("char","interval",paste("st_",((1:max(nstates))-1),sep=""));

output <- list(kendall_summary,trend_summary,state_distributions);
names(output) <- c("Kendall_Summary","Trend_Summary","State_Histograms");

return(output);
#hist(kendall_summary$pval,breaks=((0:100)/100))

#cbind(assessed_characters[significant_cases==max(significant_cases)],significant_cases[significant_cases==max(significant_cases)])
#(1:nchars)[significant_cases==max(significant_cases)]
}

#median(c(0.136, 0.008, 0.078, 0.054, 0.156, 0.235, 0.855, 0.133, 0.065, 0.016, 0.169, 0.903, 0.013, 0.160, 0.591, 0.021, 0.286, 0.124, 0.006, 0.727, 0.142, 0.297,0.115, 0.034, 0.025, 0.074, 0.810, 0.100, 0.373, 0.019, 0.095, 0.590, 0.168, 0.137, 0.088, 0.021, 0.339, 0.350, 0.158, 0.294, 0.132, 0.163, 0.084, 0.034,0.112, 0.058, 0.731, 0.162, 0.359, 0.182, 1.093, 0.066, 0.074, 0.359, 0.488, 0.166, 0.385, 0.450, 0.504, 0.300, 0.293, 0.261, 0.654, 0.025, 0.219, 0.012,0.282, 0.099, 0.184, 0.126, 0.072, 0.082, 0.287, 0.684, 0.020, 0.398, 0.298, 0.444, 0.153, 0.126, 0.048, 0.032, 0.305, 0.015, 0.232, 0.328, 0.121, 0.031,0.748, 0.001, 0.023, 0.113, 0.312, 0.018, 0.636, 0.032, 0.260, 0.116, 0.033, 0.642))
# added 2020-02-25
rescore_character_matrix_by_first_appearances <- function(chmatrix,otu_fas,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
for (nch in 1:nchars)	{
	otu_states <- chmatrix[,nch];
	if (sum(!otu_states %in% c(UNKNOWN,INAP))>0)
		chmatrix[,nch] <- rescore_states_by_first_appearances(otu_states,otu_fas,UNKNOWN,INAP);
	}
return(chmatrix);
}

rescore_states_by_first_appearances <- function(otu_states,otu_fas,UNKNOWN=-11,INAP=-22)	{
# remove taxa with inapplicable & unknown conditions
otus <- (1:length(otu_states))[!otu_states %in% c(UNKNOWN,INAP)];
otu_fas_relv <- otu_fas[!otu_states %in% c(UNKNOWN,INAP)];
otu_states_relv <- otu_states[!otu_states %in% c(UNKNOWN,INAP)];

# deal with polymorphics:
poly_otus <- (1:length(otu_states_relv))[otu_states_relv<0];
poly_otus_orig <- otus[otu_states_relv<0];
poly_scores <- otu_states_relv[otu_states_relv<0];
if (length(poly_otus)>0)	{
	for (po in 1:length(poly_otus))	{
		psts <- unravel_polymorph(poly_scores[po]);	# observed_polymorphic states
		otu_states_relv[poly_otus[po]] <- psts[1];
		otu_states_relv <- c(otu_states_relv,psts[2:length(psts)]);
		otu_fas_relv <- c(otu_fas_relv,rep(otu_fas_relv[poly_otus[po]],length(psts)-1));
		}
	}
relv_states <- sort(unique(otu_states_relv));
state_appearance <- c();
for (rs in 1:length(relv_states))
	state_appearance <- c(state_appearance,min(otu_fas_relv[otu_states_relv==relv_states[rs]]));

new_coding <- otu_states;
if (sum(state_appearance==state_appearance[order(state_appearance)]) < length(state_appearance))	{
	states_by_appearances <- relv_states[order(state_appearance)];
	for (rs in 1:length(relv_states))
		new_coding[otu_states==relv_states[rs]] <- match(relv_states[rs],states_by_appearances)-1;
#	plot(new_coding[new_coding>=0],otu_states[otu_states>=0])
	po <- 0;
	while (po < length(poly_otus))	{
		po <- po+1;
		new_coding[poly_otus_orig[po]] <- ravel_polymorph(relv_states[match(unravel_polymorph(otu_states[poly_otus_orig[po]]),states_by_appearances)])
		}
	}
return(new_coding);
}

##### PROBABILITY 101 ####
probability_to_Poisson_rate <- function(pn)	{
return(-1*log(1-pn))
}

Poisson_rate_to_probability <- function(expectation)	{
return(1-exp(-expectation))
}

##### KLUGES! #####
# clean NA from matrix
evanesco_na_from_matrix <- function(data, replacement="")  {
for (i in 1:ncol(data))	{
	if(sum(is.na(data[,i]))>0)	{
		duds <- (1:nrow(data))[is.na(data[,i])]
		data[duds,i] <- replacement
		}
	}
return(data)
}

# clean NA from vector
evanesco_na_from_vector <- function(data, replacement="")	{
if(sum(is.na(data))>0)	{
	duds <- (1:length(data))[is.na(data)]
	data[duds] <- replacement
	}
return(data)
}

# Change String To Title Case 
transmogrify_to_title_case <- function(name) {
name_part <- strsplit(name, " ")[[1]]
return (paste(toupper(substring(name_part, 1,1)), substring(name_part, 2),sep="", collapse=" "));
}

#find_me <- paleodb_clean_member_basic_no_rock[need_info_member[3]]
matrix_row_match <- function(find_me,examined_matrix)	{
xxx <- which(examined_matrix==find_me,arr.ind=TRUE);
#examined_matrix[unique(xxx[,1]),]
if (nrow(xxx)==0)	{
	return(-1);
	} else if (length(unique(xxx[,1]))==1)	{
	return(xxx[1,1])
	} else {
	vvv <- unique(xxx[,1]);
	mx <- 0;
	mj <- 0;
	for (j in 1:length(vvv))	{
		if (sum(vvv[j]==xxx[,1])>mx)	{
			mj <- j;
			mx <- sum(vvv[j]==xxx[,1])
			} else if (sum(vvv[j]==xxx[,1])==mx)	{
			mj <- c(mj,j)
			}
		}
	if (length(mj)==1)	{
		return(vvv[mj]);
		} else {
		return(-2);
		}
	}
}

#character_array <- strsplit(rock_unit_name,split = "");
simplify2vector <- function(character_array)	{
return(simplify2array(character_array)[,1])
}

revelare_file_type <- function(filename)	{
file_info <- strsplit(filename,split="")[[1]];
f_i <- length(file_info);
f_i_s <- 1+(1:f_i)[file_info %in% "."];
file_type <- c();
for (i in f_i_s:f_i)	file_type <- paste(file_type,file_info[i],sep="");
return(file_type);
}

insert_cell_into_vector_x <- function(x,new_value,cell_no)	{
# routine to put take vector x and add a new value at some point (cell_no)
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (cell_no==1)	{
	return(c(new_value,x));
	} else if (cell_no==(length(x)+1))	{
	return(c(x,new_value));
	} else	{
	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

insert_row_into_matrix_x <- function(x,new_row,row_no)	{
# routine to put take matrix x and add a new row at that bumps other rows down one
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (row_no==1)	{
	xx <- rbind(new_row,x);
	colnames(xx) <- colnames(x);
	return(xx);
	} else if (row_no==(nrow(x)+1))	{
	return(rbind(x,new_row));
	} else	{
	return(rbind(x[1:(row_no-1),],new_row,x[(row_no:nrow(x)),]));
#	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

accersi_stage_info <- function()	{
interval <- c("Ediacaran","Fortunian","Stage 2","Stage 3","Stage 4","Wuliuan","Drumian","Guzhangian","Paibian","Jiangshanian","Stage 10","Tremadoc","Floian","Dapingian","Darriwilian","Sandbian","Katian","Hirnantian","Rhuddanian","Aeronian","Telychian","Sheinwoodian","Homerian","Gorstian","Ludfordian","Pridoli","Lochkovian","Pragian","Emsian","Eifelian","Givetian","Frasnian","Famennian","Tournaisian","Visean","Serpukhovian","Bashkirian","Asselian","Sakmarian","Artinskian","Kungurian","Roadian","Wordian","Capitanian","Wuchiapingian","Changhsingian");
interval <- c(interval,"Induan","Olenekian","Anisian","Ladinian","Carnian","Norian","Rhaetian","Hettangian","Sinemurian","Pliensbachian","Toarcian","Aalenian","Bajocian","Bathonian","Callovian","Oxfordian","Kimmeridgian","Tithonian","Berriasian","Valanginian","Hauterivian","Barremian","Aptian","Albian","Cenomanian","Turonian","Coniacian","Santonian","Campanian","Maastrichtian","Paleocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Holocene");
onset <- c(635.0,538.7,529.0,521.0,514.0,509.0,504.5,500.5,497.0,494.0,489.5,485.4,477.7,470.0,467.3,458.4,453.0,445.2,443.8,440.8,438.5,433.4,430.5,427.4,425.6,423.0,419.2,410.8,407.6,393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,298.9,295.5,290.1,279.3,272.3,268.8,265.1,259.9,254.1);
onset <- c(onset,252.2,251.2,247.2,242,237,228,208.5,201.3,199.3,190.8,182.7,174.1,170.3,168.3,166.1,163.5,157.3,152.1,145,139.8,132.9,129.4,125,113,100.5,93.9,89.8,86.3,83.6,72.1,66,56,33.9,23.03,5.333,2.588,0.0117);
end <- c(538.7,529.0,521.0,514.0,509.0,504.5,500.5,497.0,494.0,489.5,485.4,477.7,470.0,467.3,458.4,453.0,445.2,443.8,440.8,438.5,433.4,430.5,427.4,425.6,423.0,419.2,410.8,407.6,393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,315.5,295.5,290.1,279.3,272.3,268.8,265.1,259.9,254.1,252.2);
end <- c(end,251.2,247.2,242,237,228,208.5,201.3,199.3,190.8,182.7,174.1,170.3,168.3,166.1,163.5,157.3,152.1,145,139.8,132.9,129.4,125,113,100.5,93.9,89.8,86.3,83.6,72.1,66,56,33.9,23.03,5.333,2.588,0.0117,0);
color <- c("#FED96A","#99B575","#A6BA80","#A6C583","#B3CA8E","#B3D492","#BFD99D","#CCDFAA","#CCEBAE","#D9F0BB","#E6F5C9","#33A97E","#41B087","#66C092","#74C69C","#8CD094","#99D69F","#A6DBAB","#A6DCB5","#B3E1C2","#BFE6CF","#BFE6C3","#CCEBD1","#CCECDD","#D9F0DF","#E6F5E1","#E5B75A","#E5C468","#E5D075","#F1D576","#F1E185","#F2EDAD","#F2EDC5","#8CB06C","#A6B96C","#BFC26B","#99C2B5","#E36350","#E36F5C","#E37B68","#E38776","#FB8069","#FB8D76","#FB9A85","#FCB4A2","#FCC0B2");
color <- c(color,"#A4469F","#B051A5","#BC75B7","#C983BF","#C99BCB","#D6AAD3","#E3B9DB","#4EB3D3","#67BCD8","#80C5DD","#99CEE3","#9AD9DD","#A6DDE0","#B3E2E3","#BFE7E5","#BFE7F1","#CCECF4","#D9F1F7","#8CCD60","#99D36A","#A6D975","#B3DF7F","#BFE48A","#CCEA97","#B3DE53","#BFE35D","#CCE968","#D9EF74","#E6F47F","#F2FA8C","#FDA75F","#FDB46C","#FDC07A","#FFFF00","#FFFF99","#FFF2AE","#FEF2E0");
return(data.frame(interval=as.character(interval),onset=as.numeric(onset),end=as.numeric(end),color=as.character(color),stringsAsFactors = hell_no));
}
