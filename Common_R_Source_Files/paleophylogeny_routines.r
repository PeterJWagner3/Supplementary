#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal

### udpated with routines written for first phylogenetics course at UNL
UNKNOWN <- -11; 
INAP <- -22;

#### ROUTINES RELATED TO BASIC INFORMATION ABOUT TREE TOPOLOGY ####
accersi_clade_reunion_given_paleodb_data <- function(clade_members)	{
paleodb_taxonomic_data <- accersi_taxonomic_data_for_list_of_taxa(taxon_list=clade_members);
basic_taxonomic_data <- evanesco_na_from_matrix(data=paleodb_taxonomic_data$entered_taxa,replacement="");
paleodb_stored_ranks <- c("phylum","class","order","family","genus");
paleodb_not_entered_default <- c("NO_PHYLUM_SPECIFIED","NO_CLASS_SPECIFIED","NO_ORDER_SPECIFIED ","NO_FAMILY_SPECIFIED");
relv_columns <- (1:ncol(basic_taxonomic_data))[colnames(basic_taxonomic_data) %in% paleodb_stored_ranks];
basic_taxonomic_data[,relv_columns];
}

accersi_branch_order_with_stratigraphy <- function(mtree,FAs)	{
branch_order <- accersi_patristic_distance_from_base(mtree);	# order branches from base of the tree
used_order <- branch_order + (1-(1:length(branch_order))/100);	# for ties, put nodes before otus
strat_order <- date_taxa_on_tree_simple(mtree,FAs);

tosort <- cbind(strat_order,used_order);
return(order(tosort[,1],tosort[,2],decreasing=FALSE))
}

accersi_patristic_distance_from_base <- function(atree)	{
ttus <- max(atree)
pat_dist_from_base <- vector(length=ttus)
if (length(dim(as.array(atree)))==1)	{
	base <- min(atree[atree>0])
	Nnode <- 1+max(atree)-base
	pat_dist_from_base[base] <- 0
	for (n in (base+1):ttus)	pat_dist_from_base[n] <- 1+pat_dist_from_base[atree[n]]
	for (s in 1:(base-1))		pat_dist_from_base[s] <- 1+pat_dist_from_base[atree[s]]
	}	else 	{
	Nnode <- dim(atree)[1]
	otus <- ttus-Nnode
	base <- otus+1
	pat_dist_from_base[base] <- 0
	for (tx in 1:Nnode)	{
		f1 <- sum(atree[tx,]>0)
		ht <- tx+otus
		for (f in 1:f1)	{
			f2 <- atree[tx,f]
			pat_dist_from_base[f2] <- pat_dist_from_base[ht]+1
			}
		}
	}
return(pat_dist_from_base)
}

transfigure_hard_polytomies <- function(rmtree,apos,rfas)	{
Nnode <- dim(rmtree)[1]
node_div <- vector(length=Nnode)
for (n in 1:Nnode)	{
	node_div[n] <- sum(rmtree[n,]>0)
	}
mxtmy <- max(node_div)
polys <- (1:Nnode)[node_div>2]
for (p in 1:length(polys))	{
	pn <- polys[p]
	f1 <- rmtree[pn,1:node_div[pn]]
	if (!is.na(match(0,apos[f1])))	{
		anc <- f1[match(0,apos[f1])]
		ds <- f1[!f1 %in% anc]
		ds <- ds[order(rfas[ds])]
		for (dd in 2:length(ds))	{
			rmtree <- rbind(rmtree,c(anc,ds[dd],rep(0,(mxtmy-2))))
			Nnode <- Nnode + 1
			}
		rmtree[pn,] <- c(anc,ds[1],rep(0,(mxtmy-2)))
		}
	}
return(rmtree[,1:2])
}

accersi_node_richness_from_vector_tree <- function(vector_tree)	{
vntree <- transform_vector_tree_to_venn_tree(vector_tree);
node_richness <- c();
for (nd in 1:nrow(vntree))
	node_richness <- c(node_richness,sum(vntree[nd,]>0));
return(node_richness);
}

tally_node_richness_from_vector_tree <- function(vector_tree)	{
vntree <- transform_vector_tree_to_venn_tree(vector_tree);
node_richness <- c();
for (nd in 1:nrow(vntree))
	node_richness <- c(node_richness,sum(vntree[nd,]>0));
return(node_richness);
}

tally_node_richness_from_venn_tree <- function(venn_tree)	{
node_richness <- c();
for (nd in 1:nrow(venn_tree))	node_richness <- c(node_richness,sum(venn_tree[nd,]>0));
return(node_richness);
}

# matrix_tree <- mat_tree;  divergence_times <- branch_dates_re;
stratify_matrix_tree_and_divergence_times <- function(matrix_tree,divergence_times,ancestral_taxa) {
# written 2022-09-22
# matrix_tree: as in all the other routines.
# divergence_times: $origin gives when branch starts, $divergence gives when branch ends; origin[i] should be <= divergence[]
colnames(divergence_times) <- c("origin","divergence");
if (!is.data.frame(divergence_times)) divergence_times <- data.frame(divergence_times);
if (divergence_times$origin[1]>divergence_times$divergence[1]) divergence_times <- -1*divergence_times;
n_Nodes <- nrow(matrix_tree);
notu <- max(matrix_tree)-n_Nodes;
# this is giving nodes the wrong number!!! 
rownames(matrix_tree) <- notu+(1:n_Nodes);
descendant_htu <- accersi_descendant_nodes_in_matrix_tree(matrix_tree,reduce_nodes=F); # list descendant clades as htu numbers
node_heights <- vector(length=n_Nodes);
for (nd in 2:n_Nodes)  node_heights[nd] <- sum(descendant_htu==nd+notu);
node_info <- data.frame(ntu=as.numeric(1:n_Nodes),htu=as.numeric(notu+1:n_Nodes),
                        height=as.numeric(node_heights),origin=divergence_times$origin[notu+1:n_Nodes])
node_info <- node_info[order(node_info$height,-node_info$origin),];
venn_tree <- transform_matrix_tree_to_venn_tree(matrix_tree);
rownames(venn_tree) <- rownames(matrix_tree);
matrix_tree2 <- transform_venn_tree_to_matrix_tree(venn_tree[node_info$ntu,]);
rownames(matrix_tree2) <- node_info$htu;
divergence_times[notu+1:n_Nodes,] <- divergence_times[node_info$htu,];

nd <- 0;
while (nd < n_Nodes)  {
  nd <- nd+1;
  f1 <- matrix_tree2[nd,matrix_tree2[nd,]>0];
  s1 <- f1[f1<=notu];
  d1 <- f1[f1>notu];
  if (length(s1)>0) {
    new_order <- f1[order(divergence_times$origin[f1],decreasing=T)];
    matrix_tree2[nd,matrix_tree2[nd,]>0] <- c(new_order[new_order %in% ancestral_taxa],new_order[!new_order %in% ancestral_taxa]);
    }
  }
venn_tree2 <- transform_matrix_tree_to_venn_tree(matrix_tree2);
output <- list(matrix_tree2,divergence_times,venn_tree2[1,]);
names(output) <- c("matrix_tree","divergence_times","good_tree_draw_order");
return(output);
}

# matrix_tree <- mat_tree
accersi_descendant_nodes_in_matrix_tree <- function(matrix_tree,reduce_nodes=T)  {
nNodes <- nrow(matrix_tree);
notu <- max(matrix_tree)-nNodes;
rownames(matrix_tree) <- notu+(1:nNodes);
tree_nodes <- matrix_tree;
tree_nodes <- cbind(matrix_tree,array(0,dim=c(nNodes,nNodes-ncol(matrix_tree))));
for (nd in 1:nNodes)  {
  tree_nodes[nd,tree_nodes[nd,]<=notu] <- 0;
  tree_nodes[nd,] <- sort(tree_nodes[nd,],decreasing=T)
  tree_nodes[nd,tree_nodes[nd,]>0] <- sort(tree_nodes[nd,tree_nodes[nd,]>0]);
  }
if (!is.null(rownames(matrix_tree)))  {
  rownames(tree_nodes) <- rownames(matrix_tree);
  } else  {
  rownames(tree_nodes) <- 1:nNodes;
  }
if (reduce_nodes) {
  tree_nodes <- tree_nodes[rowSums(tree_nodes)>0,];
  nNodes2 <- nrow(tree_nodes);
  } else  {
  nNodes2 <- nNodes;
  }
nd <- nNodes2;
while (nd > 1)  {
  htu <- nd+notu;
  gramps <- which(tree_nodes==htu,arr.ind = T);
  grampsr <- gramps[,1];
  grampsc <- gramps[,2];
  dcl <- tree_nodes[nd,tree_nodes[nd,]>0];
  d1 <- length(dcl);
  if (d1>0) {
    tree_nodes[grampsr,((grampsc+d1+1):nNodes)] <- tree_nodes[grampsr,(grampsc+1):(nNodes-d1)];
    tree_nodes[grampsr,(grampsc+1):(grampsc+d1)] <- dcl;
    }
  nd <- nd-1;
  }
tree_nodes <- tree_nodes[,colSums(tree_nodes)>0];

return(tree_nodes);
}

whatitomy <- function(vector_tree)	{
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
tomy <- vector(length=nrow(mtree));
for (i in 1:nrow(mtree))	tomy[i] <- sum(mtree[i,]>0)
return(tomy);
}

#rptree <- vector_tree
#rfas <- c(-10,-9,-6,-5)
#apos <- c(1,0,1,1)
#rmtree <- transform_vector_tree_to_matrix_tree(vector_tree=rptree)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Character reconstruction routines (likelihood, parsimony, etc.) ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
set_marginals <- function(taxon_names,chmatrix,nstates,UNKNOWN=-11,INAP=-22)	{
notu <- length(taxon_names);
nttu <- (2*notu)-1;
ncharss <- length(nstates)
mxstates <- max(nstates)
marginals <- array(0,dim=c((2*notu)-1,ncharss,mxstates))
node_names <- vector(length=notu-1)
for (n in 1:(notu-1))	{
	if (n<10)	{
		node_names[n] <- paste("Node_0",n,sep="")
		}	else {
		node_names[n] <- paste("Node_",n,sep="")
		}
	}
names(marginals[1:nttu,,]) <- c(taxon_names,node_names)
char_names <- vector(length=ncharss)
for (ch in 1:ncharss)	{
	if (ch<10)	{
		if (ncharss<100)	{
			char_names[ch] <- paste("Char_0",ch,sep="")
			}	else	{
			char_names[ch] <- paste("Char_00",ch,sep="")	
			}
		}	else if (ch<100)	{
		if (ncharss<100)	{
			char_names[ch] <- paste("Char_0",ch,sep="")
			}	else	{
			char_names[ch] <- paste("Char_",ch,sep="")
			}
		}	else	char_names[ch] <- paste("Char_",ch,sep="")
	}
names(marginals[,1:ncharss,]) <- char_names;
state_names <- vector(length=mxstates);
for (st in 1:mxstates)	{
	state_names[st] <- paste("State_",st,sep="")
	}
names(marginals[,,1:mxstates]) <- state_names;

state_list <- list_states(chmatrix)

for (ch in 1:ncharss)	{
  for (sp in 1:notu)	{
		if (chmatrix[sp,ch]>=0)	{
			marginals[sp,ch,match(chmatrix[sp,ch],state_list[[ch]])] <- 1;
			}	else if (chmatrix[sp,ch]==UNKNOWN || chmatrix[sp,ch]==INAP)	{
			for (st in 1:nstates[ch])	marginals[sp,ch,st] <- 1/nstates[ch];
			}	else {
			polym <- unravel_polymorph(chmatrix[sp,ch]);
			marginals[sp,ch,match(polym,state_list[[ch]])] <- 1/length(polym);
			}
		}
	marginals[(notu+1):((2*notu)-1),ch,1:nstates[ch]] <- 1;
	}

#rownames(marginals[,3,]) <- c(taxon_names,node_names);
return(marginals);
}

# get minimum possible divergence times with time separating each taxon & htu from its ancestral node
accersi_minimum_divergences <- function(vector_tree,strat_ranges,apos,min_prec=0.1)	{
# vector_tree: vector giving node from which lineage evolved
# strat_ranges: matrix giving first & last appearances
# apos: apomorphies along each branch
# min_prec: minimum
colnames(strat_ranges)[colnames(strat_ranges) %in% c("fa","FAD","fad")] <- "FA";
colnames(strat_ranges)[colnames(strat_ranges) %in% c("la","LAD","lad")] <- "LA";
notu <- min(vector_tree[vector_tree>0])-1;
gaps <- accersi_gaps_from_vector_tree(vector_tree,strat_ranges[1:notu,],apos);
node_dates <- date_clades_from_vector_tree(vector_tree,FA=strat_ranges$FA[1:notu]);

# divergence$onset will give appearance time
# divergence$end will give FA for otus and divergences of daughters for htu
divergences <- data.frame(onset=as.numeric(node_dates),end=as.numeric(node_dates),stringsAsFactors = F);
otu_ranges <- abs(strat_ranges[,2]-strat_ranges[,1])
mtree <- transform_vector_tree_to_matrix_tree(vector_tree)
Nnode <- nrow(mtree);
notu <- nrow(strat_ranges);
ohtu <- length(vector_tree);
node_poss_anc <- accersi_poss_ancestors_for_nodes(vector_tree,FA=strat_ranges$FA,apos)[(notu+1):ohtu];
for (nn in Nnode:1)	{
	# if two sister species, then use nu & sum(apos) to get expected time
	#	separating them.  Remember, the time is divided!
	htu <- nn+notu;
	f1 <- mtree[nn,mtree[nn,]>0];
	if (node_poss_anc[nn]==0)	{
#		divergences$end[htu] <- divergences$onset[f1] <- min(divergences$end[f1])-min_prec;
		divergences$end[htu] <- divergences$onset[f1] <- min(divergences$end[f1])-min_prec;
		if (divergences$onset[htu] >= divergences$end[htu])
			divergences$onset[htu] <- divergences$end[htu]-min_prec;
		} else	{
		anc <- node_poss_anc[nn];
		f1a <- f1[f1!=anc];
    divergences$onset[f1a[divergences$onset[f1a]==divergences$end[f1a]]] <- divergences$onset[f1a[divergences$onset[f1a]==divergences$end[f1a]]]-min_prec;
		latest_divergence <- max(divergences$onset[f1a]);
		baby_sister <- f1a[match(latest_divergence,divergences$onset[f1a])];
		big_sisters <- f1a[f1a!=baby_sister];
		if (divergences$onset[baby_sister]==divergences$end[baby_sister])
		  latest_divergence <- divergences$onset[baby_sister] <- divergences$onset[baby_sister]-min_prec;
		if (strat_ranges$LA[anc] < latest_divergence)
		  # make sure upper bound of node ranges up to latest descendant
		  divergences$end[htu] <- latest_divergence;
		earliest_divergence <- min(c(divergences$onset[anc],divergences$onset[big_sisters]));
		if (earliest_divergence < divergences$onset[anc])
		  divergences$onset[anc] <- earliest_divergence;
    if (earliest_divergence < divergences$onset[htu])
#		divergences$onset[f1] <- divergences$end[anc] <- min(divergences$onset[f1a])-min_prec; # inlude ancestor here!
#		divergences$end[htu] <- min(divergences$end[htu],min(divergences$onset[f1a]))-min_prec;
		if (divergences$end[htu] <= divergences$onset[htu])
			divergences$onset[htu] <- divergences$end[htu]-min_prec;
		### this will need safeguards for cases where 1+ descendant is as old as or older than the ancestor
		}
#	nn <- nn-1
	}
#divergences$onset[notu+1] <- divergences$end[notu+1];
return(divergences)
}

accersi_simple_clock_divergences <- function(vector_tree,strat_ranges,apos,dstates,min_prec=0.1)	{
# vector_tree: vector giving node from which lineage evolved
# strat_ranges: matrix giving first & last appearances
# apos: apomorphies along each branch
# dstate: numbers of derived n_states
# min_prec: minimum
notu <- min(vector_tree[vector_tree>0])-1;
gaps <- accersi_gaps_from_vector_tree(vector_tree,strat_ranges[1:notu,],apos);
dates <- date_clades_from_vector_tree(vector_tree,FA=strat_ranges$FA[1:notu]);
# divergence$onset will give appearance time
# divergence$end will give FA for otus and divergences of daughters for htu
divergences <- data.frame(onset=as.numeric(dates),end=as.numeric(dates),stringsAsFactors = F);
otu_ranges <- abs(strat_ranges[,2]-strat_ranges[,1])
inu <- (sum(apos)/dstates)/(sum(gaps)+sum(otu_ranges))
mtree <- transform_vector_tree_to_matrix_tree(vector_tree)
Nnode <- nrow(mtree);
notu <- nrow(strat_ranges);
ohtu <- length(vector_tree);
node_poss_anc <- accersi_poss_ancestors_for_nodes(vector_tree,FA=strat_ranges$FA,apos)[(notu+1):ohtu]
for (nn in Nnode:1)	{
	# if two sister species, then use nu & sum(apos) to get expected time
	#	separating them.  Remember, the time is divided!
	htu <- nn+notu;
	if (node_poss_anc[nn]==0)	{
		tictoc <- max(min_prec,(sum(apos[mtree[nn,]])/dstates)/inu);
		f1 <- mtree[nn,mtree[nn,]>0]
		ghost <- sum(divergences$end[f1]-min(divergences$end[f1]));
		if (tictoc > ghost && sum(f1<=notu)==length(f1))	{
			divergences$end[htu] <- divergences$onset[f1] <- round((min(divergences$end[f1])-(tictoc-ghost)/2),1)
			if (divergences$onset[htu] > divergences$end[htu]) divergences$onset[htu] <- divergences$end[htu];
			} else if (var(divergences$onset[f1])!=0)	{
			bds <- tictocs <- c()
			for (f in 1:length(f1))	{
				tictocs <- c(tictocs,(apos[f1[f]]/dstates)/inu)
				bds <- c(bds,divergences$onset[f1[f]]-min(divergences$onset[f1]))
				}
			ext <- min_prec
			ble <- 0
			lbd <- 1
			while (ble <= lbd)	{
				lbd <- 1
				for (f in 1:length(f1))
					lbd <- lbd*dpois(apos[f1[f]],inu*dstates*(bds[f]+ext));
				if (ble < lbd)	{
					ble <- lbd
					ext <- ext + min_prec;
					}
				#print(c(ble,lbd))
				}
			ext <- ext-min_prec;
			divergences$end[htu] <- divergences$onset[f1] <- min(divergences$onset[f1])-ext;
			if (divergences$onset[htu] > divergences$end[htu])
				divergences$onset[htu] <- divergences$end[htu]-min_prec;
			}
		# end case of sister-taxa
		} else	{
			# we have an ancestor: see if taxa can diverge within alotted time
			anc <- node_poss_anc[nn];
			# list descendants (1 for bifurcation)
			f1 <- mtree[nn,mtree[nn,]!=node_poss_anc[nn]];
			f1 <- f1[f1>0];
			for (f in 1:length(f1))	{
				tictoc <- (apos[f1[f]]/dstates)/inu
				if (divergences[f1[f],2] > strat_ranges[anc,2])	{
				# possible case of anagenesis: put f1 divergence as ancestral LA
					divergences[f1[f],1] <- strat_ranges[anc,2]
					}	else	{
					divergences[f1[f],1] <- divergences[f1[f],1]-tictoc
					# if this pushes descendant divergence below ancestral FA
					#	then modify divergence date of ancestor & node
					if (divergences[anc,2] > divergences[f1[f],1])	{
						# remember, ancestral morphology must be present by this
						#	time: so, we need to reset the upper bound on how long
						#	the ancestor had to accumulate change
						divergences[anc,2] <- divergences[f1[f],1]
						divergences[anc,1] <- min(divergences[anc,])
						}
					}
#				divergences[f1[f],1] <- min(strat_ranges[anc,2],divergences[f1[f],1]-tictoc)
#				if (divergences[anc,2]>divergences[f1[f],1])
#					divergences[htu,2] <- divergences[anc,1] <- divergences[f1[f],1]-min_prec
				} # end calibration of descendants
			# Make sure that all divergence time for node is that of 
			#	ancestral species && that ancestral species has no time for
			#	accumulating any change.
			divergences[htu,2] <- divergences[anc,2]
			divergences[htu,1] <- min(divergences[htu,1],divergences[anc,1])
			divergences[anc,1] <- divergences[anc,2]
		# end case of ancestor-descendant nodes
		}
#	nn <- nn-1
	}
return(divergences);
}

# ROUTINES RELATED TO MINIMUM STEPS PARSIMONY GIVEN A TREE
# updated 2020-03-13
# updated 2020-08-31
accersi_exhaustion_curve <- function(mtree,cdata,types,FAs=1,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# based on routine used in Wagner (2000).  Uses simple parsimony optimization
#	and (if available) stratigraphic data to order branches and show how
#	many novel n_states appear with how many changes along each branch working
#	up the tree.
# mtree: matrix tree, where each row gives the branches stemming from a node
# cdata: character data
# types: 1 for unordered, 0 for ordered
# UNKNOWN: numeric representation for "?"
# INAP: numeric representation for inapplicable or gap
# outgroup: the number of the outgroup taxon
if (length(outgroup)>1)
	outgroup <- outgroup[match(min(FAs[outgroup]),FAs[outgroup])];
	
n_states <- count_states(chmatrix=cdata);
char_state_combs <- c();
for (ch in 1:length(n_states))	{
	for (st in 1:n_states[ch])
		char_state_combs <- c(char_state_combs,paste(ch,"_",st,sep=""));
	}

char_history <- accersi_minimum_steps_history(mtree=mtree,cdata=cdata,types=types,UNKNOWN=UNKNOWN,INAP=INAP,outgroup=outgroup);
#char_history <- accersi_minimum_steps_history(mtree,cdata,types,UNKNOWN,INAP,outgroup);
initial_states <- char_history$Ancestral_Reconstructions[1,];
character_history <- char_history$Changes;
if (!is.data.frame(character_history))
	character_history <- data.frame(character_history,stringsAsFactors = F);
if (length(FAs)==1)	{
	branch_order <- accersi_patristic_distance_from_base(mtree)	# order branches from base of the tree
	used_order <- branch_order + (1-(1:length(branch_order))/100)	# for ties, put nodes before otus
	exhaust_order <- order(used_order,decreasing=FALSE)		# list branches in order
	}	else	{
	### insert ordering with stratigraphy here.	
	exhaust_order <- used_order <- accersi_branch_order_with_stratigraphy(mtree,FAs)
	}
exhaust_matrix <- matrix(0,max(character_history$dch),max(character_history$dst));
found <- total_steps <- novel_states <- 0;
exhaustion <- data.frame(branch_no=as.numeric(exhaust_order),steps=as.numeric(rep(0,length(exhaust_order))),ttl_steps=as.numeric(rep(0,length(exhaust_order))),novel_states=as.numeric(rep(0,length(exhaust_order))),stringsAsFactors = F);
branch_names <- rep("branch_",length(exhaust_order));
if (max(mtree)<100)	{
	branch_names[exhaust_order<10] <- paste(branch_names[exhaust_order<10],"0",sep="");
	} else if (max(mtree)<100)	{
	branch_names[exhaust_order<10] <- paste(branch_names[exhaust_order<10],"0",sep="");
	branch_names[exhaust_order<100] <- paste(branch_names[exhaust_order<100],"0",sep="");
	}
branch_names <- paste(branch_names,exhaust_order,sep="");
state_histories <- array(0,dim=c(length(branch_names),length(char_state_combs)));
rownames(state_histories) <- branch_names;
colnames(state_histories) <- char_state_combs;
for (b in 1:length(exhaust_order))	{
	br_ch <- sum(character_history$dbr==exhaust_order[b]);
	if (b>1)
		state_histories[b,] <- state_histories[b-1,];
	if (br_ch>0)	{
		local_change <- subset(character_history,character_history$dbr==exhaust_order[b])
		for (d in 1:br_ch)	{
#			if (br_ch>1)	{
			ch <- as.numeric(local_change$dch[d]);
			st <- as.numeric(local_change$dst[d]);
#				}	else	{
#				ch <- as.numeric(local_change$dch);
#				st <- as.numeric(local_change$dst);
#				}
			combo <- paste(ch,st,sep="_");
			cs <- match(combo,colnames(state_histories));
			if (st==initial_states[ch] && exhaust_matrix[ch,st]==0)	{
				state_histories[b,cs] <- exhaust_matrix[ch,st] <- 2;
				} else	{
				state_histories[b,cs] <- exhaust_matrix[ch,st] <- exhaust_matrix[ch,st]+1;
				}
			if (exhaust_matrix[ch,st]==1)	novel_states <- novel_states+1;
			}
		total_steps <- total_steps+br_ch;
		exhaustion$branch_no[b] <- exhaust_order[b];
		exhaustion$steps[b] <- br_ch;
		exhaustion$ttl_steps[b] <- total_steps;
		exhaustion$novel_states[b] <- novel_states;
		} else	{
		exhaustion$branch_no[b] <- exhaust_order[b];
		exhaustion$ttl_steps[b] <- total_steps;
		exhaustion$novel_states[b] <- novel_states;
		}
	}
	 #b <- b+1}
rownames(exhaustion) <- branch_names;

exhaust_output <- list(exhaustion,exhaust_matrix,character_history[order(match(character_history$dbr,exhaust_order)),],state_histories);
names(exhaust_output) <- c("exhaustion","state_derivations","character_history","state_histories");
return(exhaust_output);
}

accersi_minimum_steps_history <- function(mtree,cdata,types,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# mtree: matrix tree, where each row gives the branches stemming from a node
# cdata: character data
# types: 1 for unordered, 0 for ordered
# UNKNOWN: numeric representation for "?"
# INAP: numeric representation for inapplicable or gap
# outgroup: the number of the outgroup taxon
otus <- nrow(cdata);
nNodes <- nrow(mtree);
nchars <- ncol(cdata);
steps <- vector(length=nchars);
full_matrix <- changes_matrix <- array(0,dim=c(otus+nNodes,nchars));
for (c in 1:nchars)	{
	cvector <- cdata[,c];
	type <- types[c];
#	sort(unique(cvector));
	if (sum(cvector>=0)>0)	{
		char_evolution <- Sankoff_character(mtree,cvector,type,UNKNOWN,INAP,outgroup);
		steps[c] <- char_evolution$Steps;
		full_matrix[,c] <- char_evolution$n_states;
		changes_matrix[,c] <- char_evolution$Derivation;
#		if (c==1)	{
#			full_matrix <- char_evolution$n_states;
#			changes_matrix <- char_evolution$Derivation;
#			}	else	{
#			full_matrix <- cbind(full_matrix,char_evolution$n_states);
#			changes_matrix <- cbind(changes_matrix,char_evolution$Derivation);
#			}
		} else	{
		full_matrix[,c] <- rep(-11,max(mtree));
#		dummy_states[1:length(cvector)] <- cvector;
		changes_matrix[,c] <- rep(0,max(mtree));
#		if (c==1)	{
#			full_matrix <- dummy_states;
#			changes_matrix <- dummy_changes;
#			} else	{
#			full_matrix <- cbind(full_matrix,dummy_states);
#			changes_matrix <- cbind(changes_matrix,dummy_changes);
#			}
		}
	}
#which(is.na(changes_matrix),arr.ind = T)
char_labels <- vector(length=nchars);
for (c in 1:nchars)	{
	if (c<10)	{
		if (nchars<10)	{
			char_labels[c] <- paste("ch_",c,sep="")
			}	else if (nchars>9 && nchars<100)	{
			char_labels[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_labels[c] <- paste("ch_00",c,sep="")
			}
		# end case with < 10 characters
		} else if (c<100)	{
		if (nchars<100)	{
		char_labels[c] <- paste("ch_",c,sep="");
			}	else	{
			char_labels[c] <- paste("ch_0",c,sep="");
			}
		} else	{
		char_labels[c] <- paste("ch_",c,sep="");
		}
	}
#which(is.na(changes_matrix),arr.ind = T)
ttl_br <- otus+nNodes;
nonzero <- rbr <- 1;
for (br in 1:ttl_br)	{
	nonzero <- sum(changes_matrix[br,]>0);
#	(1:ncol(cdata))[is.na(changes_matrix[br,])]
	if (nonzero > 0)	{
		dch <- (1:nchars)[changes_matrix[br,]>0];
		dst <- changes_matrix[br,changes_matrix[br,]>0];
		dbr <- rep(br,length(dch));
		if (rbr==1)	{
			branch_changes <- cbind(dbr,dch,dst);
			} else	{
			branch_changes <- rbind(branch_changes,cbind(dbr,dch,dst));
			}
		rbr <- rbr+1;
		}
	}
rownames(branch_changes) <- rep("",nrow(branch_changes));
branch_changes <- as.data.frame(branch_changes,stringsAsFactors = F);
colnames(full_matrix) <- colnames(changes_matrix) <- char_labels
output <- list(steps,full_matrix[((otus+1):(otus+nNodes)),],branch_changes)
names(output) <- c("Steps","Ancestral_Reconstructions","Changes")
return(output);
}

Sankoff_character <- function(mtree,cvector,type=1,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# method for reconstructing ancestral conditions.
#obs_states <- sort(unique(cvector[cvector>=0]));
#cvector <- cvector[1:otus]
obs_states <- minimum_state(matrix(cvector,nrow=length(cvector),ncol=1)):maximum_state(matrix(cvector,nrow=length(cvector),ncol=1));
adj <- min(obs_states)-1;
ttl_states <- length(obs_states);
nNodes <- nrow(mtree);
otus <- length(cvector);
scored <- (1:otus);
sankoff_matrix <- matrix(1,otus+nNodes,(length=ttl_states));
for (s in 1:otus)	{
	if (cvector[s] >= 0)	{
		st <- match(cvector[s],obs_states)
		sankoff_matrix[s,st] <- 0;
		} else if (cvector[s]!=UNKNOWN && cvector[s]!=INAP) {
		combo <- as.character(cvector[s]);
		sankoff_matrix[s,] <- 1;
		for (os in max(obs_states):min(obs_states))	{
			oss <- os-adj;
		  ost <- obs_states[oss];
			ostt <- match(ost,obs_states);
			combo2 <- gsub(ost,"•",combo);
			if (combo2!=combo)	{
				sankoff_matrix[s,ostt] <- 0;
				combo <- combo2;
				}
			}
		} else	sankoff_matrix[s,] <- 0
	}
cvector <- c(cvector,rep(UNKNOWN,nNodes));
node_rich <- vector(length=nNodes);
for (n in nNodes:1)	{
	missing <- gap <- 0;
	f1s <- mtree[n,mtree[n,]>0];
	node_rich[n] <- f1 <- length(f1s);
	# if all taxa are scored
	# mtree[n,]
#	if (sum(mtree[n,(1:f1)] %in% scored)==f1)	{
	ht <- n+otus;
	sankoff_matrix[ht,] <- 0;
#	if (ncol(sankoff_matrix)>1)	{
#		sankoff_matrix[ht,] <- colSums(sankoff_matrix[f1s,]);
#		} else	{
#		sankoff_matrix[ht,] <- sum(sankoff_matrix[f1s,]);
#		}
	for (s in 1:f1)	{
		#### add something to deal with inapplicables here.
		# list n_states that demand more than minimal change
		sp <- mtree[n,s]
		if (cvector[sp]!=INAP && cvector[sp]!=UNKNOWN)	{
			#	we will add a step to each of these because if sankoff is:
			#		0 1 1 for n_states 0, 1 & 2
			#	then we need 1 step from either 1 or 2
			n_s <- obs_states[sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
			a_s <- obs_states[!sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
			no_step <- match(n_s,obs_states)
			add_step <- match(a_s,obs_states)
			sankoff_matrix[ht,no_step] <- sankoff_matrix[ht,no_step]+min(sankoff_matrix[sp,])
			if (length(no_step) < length(obs_states))	{
				sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+min(sankoff_matrix[sp,])+1
	#			sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+sankoff_matrix[sp,add_step]
	#			sankoff_matrix[ht,]
				}
			}	else if (cvector[sp]==INAP)	{
			gap <- gap+1;
			}	else if (cvector[sp]==UNKNOWN)	{
			missing <- missing+1;
			}
		}
	if (missing==f1 || (missing>0 && (missing+gap)==f1))	{
		cvector[ht] <- UNKNOWN
		scored <- c(scored,ht)
		}	else if (gap==f1)	{
		cvector[ht] <- INAP
		scored <- c(scored,ht)
		}	else if (sum(sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,]))==1)	{
		# see if we can fix a state
		ast <- match(min(sankoff_matrix[ht,]),sankoff_matrix[ht,])
		cvector[ht] <- obs_states[ast]
		scored <- c(scored,ht)
		}	else	{
		cvector[ht] <- ravel_polymorph(obs_states[sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,])])
		}
#	n <- n-1;
	}

# this somehow missed the one change.
base <- otus + 1;
ch_steps <- min(sankoff_matrix[base,]);

if (cvector[base]<0 && (cvector[base]!=UNKNOWN && cvector[base]!=INAP))	{
	#poss_starts <- (1:ttl_states)[sankoff_matrix[base,] %in% min(sankoff_matrix[base,])]
	poss_starts <- unravel_polymorph(cvector[base])
	# IF the designated outgroup is attached to the first node AND if it 
	#	has one of the possible nodal n_states, then assign that to basal node
	if (!is.na(match(outgroup,mtree[1,])) && !is.na(match(cvector[outgroup],poss_starts)))	{
		cvector[base] <- cvector[outgroup]
		}	else	{
		# otherwise, just grab one of them at random: it really doesn't matter
		grab <- ceiling(runif(1)/(1/length(poss_starts)))
		cvector[base] <- poss_starts[grab]
		}
	}

### start here: work up the tree, using cvector[htu] to set the state
changes_above <- vector(length=nNodes)
for (n in 1:nNodes)	{
	ht <- n+otus;
	if (cvector[ht]!=INAP && cvector[ht]!=UNKNOWN)	{
		f1 <- node_rich[n]
		# if all taxa are scored in mtree[n,]
		if (sum(mtree[n,(1:f1)] %in% scored)<f1)	{
			anc_st <- cvector[ht]
			unscored <- mtree[n,!mtree[n,(1:f1)] %in% scored]
			unscored <- unscored[unscored>0]
			for (u in 1:length(unscored))	{
				if (length(unscored)>1)	{
					f2 <- unscored[u]
					}	else	{
					f2 <- unscored
					}
				# if parent value matches one of the possible values for the daughter node
				#	then go with that value
				if(anc_st>=0 && (!is.na(match(anc_st,obs_states[sankoff_matrix[f2,] %in% min(sankoff_matrix[f2,])]))))	{
					cvector[f2] <- anc_st
					}	else if (anc_st<0 && (anc_st!=UNKNOWN && anc_st!=INAP))	{
					anc_poss <- unravel_polymorph(anc_st)
					f2_poss <- unravel_polymorph(cvector[f2])
					reduced_poss <- f2_poss[f2_poss %in% anc_poss]
					if (length(reduced_poss)==1)	{
						cvector[ht] <- cvector[f2] <- reduced_poss
						}	else	{
	#					cvector[ht] <- cvector[f2] <- ravel_polymorph(reduced_poss)
						# if it is still up in the air at this point, then you might as
						#	well roll the dice!
						cvector[f2] <- reduced_poss[ceiling(runif(1)/(1/length(reduced_poss)))]
						}
					}	else	{
				# it really does not matter at this point what state we pick: it will all be the same.
					f2_poss <- unravel_polymorph(cvector[f2])
					cvector[f2] <- f2_poss[ceiling(runif(1)/(1/length(f2_poss)))]
					}
					# end case where we can assign ancestral condition to daughter node
				}	# end search of unscored nodes
			}	# I could add a routine here to compare sets of equally parsimonious n_states
		changes_above[n] <- min(sankoff_matrix[n+otus,])
		}		# if node has (0,1) and daughter node has (1,2), then go with 1
	}
cchanges <- rep(0,length(cvector))
for (n in 1:nNodes)	{
	ht <- otus+n
	if (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP)	{
		f1 <- node_rich[n]
		for (f in 1:f1)	{
			f2 <- mtree[n,f]
			if (cvector[f2] != cvector[ht])	{
				if ((cvector[f2]!=UNKNOWN && cvector[f2]!=INAP) && (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP))	{
					# if both taxa are resolved, then this is easy
					if (cvector[f2]<0)	{
						ply <- unravel_polymorph(cvector[f2])
						if (is.na(match(cvector[ht],ply)))	cchanges[f2] <- match(min(ply),obs_states)
						}	else if (length(cvector[f2])==1 && length(cvector[ht])==1)	{
						cchanges[f2] <- match(cvector[f2],obs_states)
						}
					}
				}
			}
		}
	}
output <- list(ch_steps,cvector,cchanges)
names(output) <- c("Steps","n_states","Derivation")
return(output)
}

Sankoff_character_effed <- function(mtree,cvector,type=1,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# 2020-09-01: updated to accommodate n_states appearing only in polymorphics
# method for reconstructing ancestral conditions.
# cvector <- cvector[1:otus]
obs_states <- 1:maximum_state(chmatrix=array(cvector,dim=c(length(cvector),1)))
max_st <- max(obs_states);
min_st <- min(obs_states);
ttl_states <- length(obs_states);
nNodes <- nrow(mtree);
otus <- length(cvector);
scored <- (1:otus);
sankoff_matrix <- matrix(1,otus+nNodes,(length=ttl_states));
for (s in 1:otus)	{
	if (cvector[s] >= 0)	{
		st <- match(cvector[s],obs_states)
		sankoff_matrix[s,st] <- 0
		}	else if (cvector[s]==UNKNOWN || cvector[s]==INAP)	{
		sankoff_matrix[s,] <- 0
		} else	{
		if (max_st<10)	{
			sankoff_matrix[s,unravel_polymorph(cvector[s])] <- 0;
			} else	{
			sankoff_matrix[s,unravel_polymorph_badass(cvector[s],minst = min_st)] <- 0;
			}
		}
	}
cvector <- c(cvector,rep(UNKNOWN,nNodes));
node_rich <- vector(length=nNodes);
#for (nn in nNodes:(87-otus))	{
for (nn in nNodes:1)	{
	ht <- nn+otus;
	missing <- gap <- 0;
	node_rich[nn] <- f1 <- length(mtree[nn,mtree[nn,]>0]);
	# if all taxa are scored
	# mtree[n,]
#	if (sum(mtree[n,(1:f1)] %in% scored)==f1)	{
#	sankoff_matrix[ht,] <- 0;
	if (ncol(sankoff_matrix)>1)	{
		sankoff_matrix[ht,] <- colSums(sankoff_matrix[mtree[nn,mtree[nn,]>0],]);
		} else	{
		sankoff_matrix[ht,] <- sum(sankoff_matrix[mtree[nn,mtree[nn,]>0],]);
		}
	gap <- sum(cvector[mtree[nn,mtree[nn,]>0]]==INAP);
	missing <- sum(cvector[mtree[nn,mtree[nn,]>0]]==UNKNOWN);
#	for (s in 1:f1)	{
		#### add something to deal with inapplicables here.
		# list n_states that demand more than minimal change
#		sp <- mtree[nn,s];
#		if (cvector[sp]!=INAP && cvector[sp]!=UNKNOWN)	{
			#	we will add a step to each of these because if sankoff is:
			#		0 1 1 for n_states 0, 1 & 2
			#	then we need 1 step from either 1 or 2
#			n_s <- obs_states[sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
#			a_s <- obs_states[!sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
#			no_step <- match(n_s,obs_states);
#			add_step <- match(a_s,obs_states);
#			sankoff_matrix[ht,no_step] <- sankoff_matrix[ht,no_step]+min(sankoff_matrix[sp,])
#			if (length(no_step) < length(obs_states))
#				sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+min(sankoff_matrix[sp,])+1
	#			sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+sankoff_matrix[sp,add_step]
	#			sankoff_matrix[ht,]
	#			}
#			}	else if (cvector[sp]==INAP)	{
#			gap <- gap+1
#			}	else if (cvector[sp]==UNKNOWN)	{
#			missing <- missing+1
#			}
#		}
	if (missing==f1 || (missing >0 && (missing+gap)==f1))	{
		cvector[ht] <- UNKNOWN;
		scored <- c(scored,ht);
		}	else if (gap==f1)	{
		cvector[ht] <- INAP;
		scored <- c(scored,ht);
		}	else if (sum(sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,]))==1)	{
		# see if we can fix a state
		ast <- match(min(sankoff_matrix[ht,]),sankoff_matrix[ht,]);
		cvector[ht] <- obs_states[ast];
		scored <- c(scored,ht);
		}	else	{
		cvector[ht] <- ravel_polymorph(obs_states[sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,])]);
		}
	}

base <- otus+1;
ch_steps <- min(sankoff_matrix[base,]);

if (cvector[base]<0 && (cvector[base]!=UNKNOWN && cvector[base]!=INAP))	{
	#poss_starts <- (1:ttl_states)[sankoff_matrix[base,] %in% min(sankoff_matrix[base,])]
	if (max_st < 10)	{
		poss_starts <- unravel_polymorph(cvector[base]);
		} else	{
		poss_starts <- unravel_polymorph_badass(poly=cvector[base],minst=min_st);
		}
	# IF the designated outgroup is attached to the first node AND if it 
	#	has one of the possible nodal n_states, then assign that to basal node
	if (!is.na(match(outgroup,mtree[1,])) && !is.na(match(cvector[outgroup],poss_starts)))	{
		cvector[base] <- cvector[outgroup]
		}	else	{
		# otherwise, just grab one of them at random: it really doesn't matter
		grab <- ceiling(runif(1)/(1/length(poss_starts)))
		cvector[base] <- poss_starts[grab]
		}
	}

### start here: work up the tree, using cvector[htu] to set the state
changes_above <- vector(length=nNodes);
for (nn in 1:nNodes)	{
	ht <- nn+otus;
	if (cvector[ht]!=INAP && cvector[ht]!=UNKNOWN)	{
		f1 <- node_rich[nn]
		# if all taxa are scored in mtree[n,]
		if (sum(mtree[nn,(1:f1)] %in% scored)<f1)	{
			anc_st <- cvector[ht];
			unscored <- mtree[nn,!mtree[nn,(1:f1)] %in% scored]
			unscored <- unscored[unscored>0]
			for (u in 1:length(unscored))	{
				if (length(unscored)>1)	{
					f2 <- unscored[u]
					}	else	{
					f2 <- unscored
					}
				# if parent value matches one of the possible values for the daughter node
				#	then go with that value
				if(anc_st>=0 && (!is.na(match(anc_st,obs_states[sankoff_matrix[f2,] %in% min(sankoff_matrix[f2,])]))))	{
					cvector[f2] <- anc_st
					}	else if (anc_st<0 && (anc_st!=UNKNOWN && anc_st!=INAP))	{
					anc_poss <- unravel_polymorph(anc_st)
					f2_poss <- unravel_polymorph(cvector[f2])
					reduced_poss <- f2_poss[f2_poss %in% anc_poss]
					if (length(reduced_poss)==1)	{
						cvector[ht] <- cvector[f2] <- reduced_poss
						}	else	{
	#					cvector[ht] <- cvector[f2] <- ravel_polymorph(reduced_poss)
						# if it is still up in the air at this point, then you might as
						#	well roll the dice!
						cvector[f2] <- reduced_poss[ceiling(runif(1)/(1/length(reduced_poss)))]
						}
					}	else	{
				# it really does not matter at this point what state we pick: it will all be the same.
					f2_poss <- unravel_polymorph(cvector[f2])
					cvector[f2] <- f2_poss[ceiling(runif(1)/(1/length(f2_poss)))]
					}
					# end case where we can assign ancestral condition to daughter node
				}	# end search of unscored nodes
			}	# I could add a routine here to compare sets of equally parsimonious n_states
		changes_above[nn] <- min(sankoff_matrix[nn+otus,])
		}		# if node has (0,1) and daughter node has (1,2), then go with 1
	}
cchanges <- rep(0,length(cvector));
#for (nn in 1:nNodes)	{
for (nn in 1:nNodes)	{
	ht <- otus+nn;
	if (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP)	{
		f1 <- node_rich[nn]
		for (f in 1:f1)	{
			f2 <- mtree[nn,f];
			if (cvector[f2] != cvector[ht])	{
				if ((cvector[f2]!=UNKNOWN && cvector[f2]!=INAP) && (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP))	{
					# if both taxa are resolved, then this is easy
					if (cvector[f2]<0)	{
						ply <- unravel_polymorph(cvector[f2])
						if (is.na(match(cvector[ht],ply)))	cchanges[f2] <- match(min(ply),obs_states)
						}	else if (length(cvector[f2])==1 && length(cvector[ht])==1)	{
						cchanges[f2] <- match(cvector[f2],obs_states);
						if (is.na(cchanges[f2]))	{
							#print(paste(nn,f2,cvector[f2]));
							# this shouldn't be necessary;
							# somehow, cvector[f2] sometimes is zero.
							cchanges[f2] <- 0;
							cvector[f2] <- cvector[ht];
							}
						}
					} # end case with coded characters
				} # end case where descendant & node mismatch
			} # end search of descendants
		} # end case of coded character
	}
output <- list(ch_steps,cvector,cchanges);
names(output) <- c("Steps","n_states","Derivation")
return(output)
}

# get the maximum number of times a character can change on a tree
#	given minimum steps.  
accersi_maximum_parsimony_steps_per_character <- function(chmatrix,chstates=NULL)	{
# routine to get the maximum number of steps for a character on a minimum steps tree.
#	This assumes a star-phylogeny (zero resolution) where one state is ancestral
#		and all other n_states are derived n times given n otus with that state
if (is.null(chstates))	chstates <- count_states(chmatrix);
nchars <- ncol(chmatrix);
max_steps <- vector(length=nchars);
for (c in 1:nchars)	{
	scored <- sum(chmatrix[,c]>=0)
	max_shared <- 0;
	for (s in 0:chstates[c])	{
		x <- sum(chmatrix[,c]==s);
		if (max_shared < x)	max_shared <- x;
		}
	max_steps[c] <- scored - max_shared;
	}
return(max_steps);
}

# get proportion of coded characters that are homoplastic
# homoplasy2 <- c();
proportional_homoplasy <- function(chmatrix,types,vector_tree,outgroup=1,UNKNOWN=-11,INAP=-22,reps=50)	{
nchars <- ncol(chmatrix);
ttl_coded <- homoplastics <- vector(length=nchars);
for (nch in 1:nchars)	ttl_coded[nch] <- notu-(sum(chmatrix[,nch]==UNKNOWN)+sum(chmatrix[,nch]==INAP));
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
notu <- nrow(chmatrix);
prob_homoplastic <- vector(length=reps);
for (rr in 1:reps)	{
	pars_estimates <- accersi_minimum_steps_history(mtree=mtree,cdata=chmatrix,types=types,UNKNOWN=UNKNOWN,INAP=INAP,outgroup=outgroup);
	ancestral_states <- pars_estimates$Ancestral_Reconstructions;
#	print(ancestral_states[1,1]);
	all_changes <- pars_estimates$Changes;
	#max(all_changes$dbr);
	branch_order <- c(1:notu,notu+(nrow(venn_tree):1));
	pr_homoplastic <- vector(length=nchars);
	for (nch in 1:nchars)	{
		char_changes <- all_changes[all_changes$dch==nch,];
		while (nrow(char_changes)==0 && nch<nchars)	{
			nch <- nch+1;
			char_changes <- all_changes[all_changes$dch==nch,];
			}
		if (nch>nchars)
			break;
		# order character changes so that we look at nodes first
		char_changes <- char_changes[order(match(char_changes$dbr,branch_order)),];
		orig_state <- ancestral_states[1,nch];
		state_richness <- sum(chmatrix[,nch]==orig_state);
		other_states <- (1:n_states[nch])[!(1:n_states[nch]) %in% orig_state];
		revert_to_orig <- subset(char_changes,char_changes$dst==orig_state);
		originals <- (1:notu)[chmatrix[,nch]==orig_state];
		ro <- 0;
		copycats <- copycat_sizes <- c();
		while (ro < nrow(revert_to_orig))	{
			ro <- ro+1;
			tbr <- revert_to_orig$dbr[ro];
			if (tbr>notu)	{
				nbr <- tbr-notu;
				tsp <- venn_tree[nbr,venn_tree[nbr,]>0];
				copycat_sizes <- c(copycat_sizes,sum(chmatrix[tsp,nch]==orig_state));
				copycats <- c(copycats,tsp[chmatrix[tsp,nch]==orig_state]);
				} else	{
				copycats <- c(copycats,tbr);
				copycat_sizes <- c(copycat_sizes,1);
				}
			}
		originals <- originals[!originals %in% copycats];	# species retaining original derivation
		copycat_sizes <- c(length(originals),copycat_sizes)
#		copy_cat_props <- copycat_sizes/sum(copycat_sizes);	# prob of grabbing state from each derivation
#		copy_cat_props_rem <- (copycat_sizes-1)/(ttl_coded[nch]-1);	# prob of regrabbing state from each derivation
		p_homoplasy <- (state_richness/ttl_coded[nch])*(1-(sum(copycat_sizes*(copycat_sizes-1))/max(1,(state_richness*(state_richness-1)))));
	#	p_homoplasy <- (state_richness/ttl_coded[nch])*(1-sum(copy_cat_props*copy_cat_props_rem));
		
		derived_states <- subset(char_changes,char_changes$dst!=orig_state);
		for (os in 1:length(other_states))	{
			this_state <- subset(derived_states,derived_states$dst==other_states[os]);
			state_richness <- c(state_richness,sum(chmatrix[,nch]==other_states[os]));
			if (nrow(this_state)>1)	{
				ro <- 0;
				spc_w_state <- c();
				derivation_richness <- c();
				while (ro < nrow(this_state))	{
					ro <- ro+1;
					tbr <- this_state$dbr[ro];
					if (tbr>notu)	{
						nbr <- tbr - notu;
						# get descendants of node with this state
						this_case <- venn_tree[nbr,venn_tree[nbr,]>0][chmatrix[venn_tree[nbr,venn_tree[nbr,]>0],nch]==other_states[os]];
						# remove descendants of node that got this state later
						this_case <- this_case[!this_case %in% spc_w_state];
						derivation_richness <- c(derivation_richness,length(this_case));
						} else	{
						derivation_richness <- c(derivation_richness,1);
						spc_w_state <- c(spc_w_state,tbr);
						}
					}
				p_homoplasy <- c(p_homoplasy,
								 (state_richness[os+1]/ttl_coded[nch])*(1-(sum(derivation_richness*(derivation_richness-1))/max(1,(state_richness[os+1]*(state_richness[os+1]-1))))));
				
#				derivation_richness_prop <- derivation_richness/sum(derivation_richness);
#				derivation_richness_prop_rem <- (derivation_richness-1)/max(1,sum((derivation_richness-1)));
#				p_homoplasy <- c(p_homoplasy,(state_richness[os+1]/ttl_coded[nch])*(1-sum(derivation_richness_prop*derivation_richness_prop_rem)));
				}
			}
		
#		pr_homoplastic <- c(pr_homoplastic,sum(p_homoplasy));
		pr_homoplastic[nch] <- sum(p_homoplasy);
		}
	prob_homoplastic[rr] <- mean(pr_homoplastic);
#	prop_homoplastic <- c(prop_homoplastic,sum(homoplastics)/sum(ttl_coded));
	} # end rr loop

#homoplasy2 <- cbind(homoplasy2,homoplastics);
#print(c(sum(homoplastics),sum(ttl_coded)));
return(median(prob_homoplastic));
}

proportional_homoplasy_old <- function(chmatrix,types,vector_tree,outgroup=1,UNKNOWN=-11,INAP=-22,reps=50)	{
nchars <- ncol(chmatrix);
ttl_coded <- homoplastics <- vector(length=nchars);
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
notu <- nrow(chmatrix);
prop_homoplastic <- c();
for (rr in 1:reps)	{
	pars_estimates <- accersi_minimum_steps_history(mtree=mtree,cdata=chmatrix,types=types,UNKNOWN=UNKNOWN,INAP=INAP,outgroup=outgroup);
	ancestral_states <- pars_estimates$Ancestral_Reconstructions;
#	print(ancestral_states[1,1]);
	all_changes <- pars_estimates$Changes;
	#max(all_changes$dbr);
	branch_order <- c(1:notu,notu+(nrow(venn_tree):1));
	for (nch in 1:nchars)	{
		ttl_coded[nch] <- notu-(sum(chmatrix[,nch]==UNKNOWN)+sum(chmatrix[,nch]==INAP));
		orig_state <- ancestral_states[1,nch];
		other_states <- (1:n_states[nch])[!(1:n_states[nch]) %in% orig_state];
		char_changes <- all_changes[all_changes$dch==nch,];
		char_changes <- char_changes[order(match(char_changes$dbr,branch_order)),];
		revert_to_orig <- subset(char_changes,char_changes$dst==orig_state);
		ro <- 0;
		copycats <- c();
		while (ro < nrow(revert_to_orig))	{
			ro <- ro+1;
			tbr <- revert_to_orig$dbr[ro];
			if (tbr>notu)	{
				nbr <- tbr-notu;
				tsp <- venn_tree[nbr,venn_tree[nbr,]>0];
				copycats <- c(copycats,tsp[chmatrix[tsp,nch]==orig_state]);
				} else	{
				copycats <- c(copycats,tbr);
				}
			}
	
		derived_states <- subset(char_changes,char_changes$dst!=orig_state);
		for (os in 1:length(other_states))	{
			this_state <- subset(derived_states,derived_states$dst==other_states[os]);
			if (nrow(this_state)>1)	{
				ro <- 0;
				spc_w_state <- c();
				derivation_richness <- c();
				while (ro < nrow(this_state))	{
					ro <- ro+1;
					tbr <- this_state$dbr[ro];
					if (tbr>notu)	{
						nbr <- tbr - notu;
						# get descendants of node with this state
						this_case <- venn_tree[nbr,venn_tree[nbr,]>0][chmatrix[venn_tree[nbr,venn_tree[nbr,]>0],nch]==other_states[os]];
						# remove descendants of node that got this state later
						this_case <- this_case[!this_case %in% spc_w_state];
						derivation_richness <- c(derivation_richness,length(this_case));
						} else	{
						derivation_richness <- c(derivation_richness,1);
						spc_w_state <- c(spc_w_state,tbr);
						}
					}
				best_case <- match(max(derivation_richness),derivation_richness);
				for (ro in 1:nrow(this_state))	{
					if (ro!=best_case)	{
						tbr <- this_state$dbr[ro];
						if (tbr>notu)	{
							nbr <- tbr - notu;
							# get descendants of node with this state
							this_case <- venn_tree[nbr,venn_tree[nbr,]>0][chmatrix[venn_tree[nbr,venn_tree[nbr,]>0],nch]==other_states[os]];
							# remove descendants of node that got this state later
							copycats <- unique(c(copycats,this_case));
							} else	{
							copycats <- c(copycats,tbr);
							}
						}
					} # finish looking at derivations of this state n_states
				} # finish examining other n_states
	#		if (nrow(this_state)>1)	{
	#			print(nch);
	#			}
			}
		homoplastics[nch] <- length(copycats);
		}
	prop_homoplastic <- c(prop_homoplastic,sum(homoplastics)/sum(ttl_coded));
	} # end rr loop

#homoplasy2 <- cbind(homoplasy2,homoplastics);
#print(c(sum(homoplastics),sum(ttl_coded)));
return(median(prop_homoplastic));
}

# get proportion of coded characters that are homoplastic
proportional_homoplasy_per_character <- function(chmatrix,types,vector_tree,outgroup=1,UNKNOWN=-11,INAP=-22)	{
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
notu <- nrow(chmatrix);
pars_estimates <- accersi_minimum_steps_history(mtree=mtree,cdata=chmatrix,types=types,UNKNOWN=UNKNOWN,INAP=INAP,outgroup=outgroup);
ancestral_states <- pars_estimates$Ancestral_Reconstructions;
all_changes <- pars_estimates$Changes;
#max(all_changes$dbr);
ttl_coded <- homoplastics <- vector(length=nchars);
branch_order <- c(1:notu,notu+(nrow(venn_tree):1));
for (nch in 1:nchars)	{
	ttl_coded[nch] <- notu-(sum(chmatrix[,nch]==UNKNOWN)+sum(chmatrix[,nch]==INAP));
	orig_state <- ancestral_states[1,nch];
	other_states <- (1:n_states[nch])[!(1:n_states[nch]) %in% orig_state];
	char_changes <- all_changes[all_changes$dch==nch,];
	char_changes <- char_changes[order(match(char_changes$dbr,branch_order)),];
	revert_to_orig <- subset(char_changes,char_changes$dst==orig_state);
	ro <- 0;
	copycats <- c();
	while (ro < nrow(revert_to_orig))	{
		ro <- ro+1;
		tbr <- revert_to_orig$dbr[ro];
		if (tbr>notu)	{
			nbr <- tbr-notu;
			tsp <- venn_tree[nbr,venn_tree[nbr,]>0];
			copycats <- c(copycats,tsp[chmatrix[tsp,nch]==orig_state]);
			} else	{
			copycats <- c(copycats,tbr);
			}
		}

	derived_states <- subset(char_changes,char_changes$dst!=orig_state);
	for (os in 1:length(other_states))	{
		this_state <- subset(derived_states,derived_states$dst==other_states[os]);
		if (nrow(this_state)>1)	{
			ro <- 0;
			spc_w_state <- c();
			derivation_richness <- c();
			while (ro < nrow(this_state))	{
				ro <- ro+1;
				tbr <- this_state$dbr[ro];
				if (tbr>notu)	{
					nbr <- tbr - notu;
					# get descendants of node with this state
					this_case <- venn_tree[nbr,venn_tree[nbr,]>0][chmatrix[venn_tree[nbr,venn_tree[nbr,]>0],nch]==other_states[os]];
					# remove descendants of node that got this state later
					this_case <- this_case[!this_case %in% spc_w_state];
					derivation_richness <- c(derivation_richness,length(this_case));
					} else	{
					derivation_richness <- c(derivation_richness,1);
					spc_w_state <- c(spc_w_state,tbr);
					}
				}
			best_case <- match(max(derivation_richness),derivation_richness);
			for (ro in 1:nrow(this_state))	{
				if (ro!=best_case)	{
					tbr <- this_state$dbr[ro];
					if (tbr>notu)	{
						nbr <- tbr - notu;
						# get descendants of node with this state
						this_case <- venn_tree[nbr,venn_tree[nbr,]>0][chmatrix[venn_tree[nbr,venn_tree[nbr,]>0],nch]==other_states[os]];
						# remove descendants of node that got this state later
						copycats <- unique(c(copycats,this_case));
						} else	{
						copycats <- c(copycats,tbr);
						}
					}
				}
			}
#		if (nrow(this_state)>1)	{
#			print(nch);
#			}
		}
	homoplastics[nch] <- length(copycats);
	}
output <- list(ttl_coded,homoplastics);
names(output) <- c("coded_otus","homoplastic_otus");
return(output);
}

# likelihood routines
seed_Mk_rate <- function(vector_tree,chmatrix,steps,durations,UNKNOWN=-11,INAP=-22)	{
# steps: vector with changes on branches leading to otus;
# add_myr deleted  
notu <- min(vector_tree[vector_tree>0]) - 1;
ghosts <- durations$FA[1:notu]-durations$LA[vector_tree[1:notu]];
st_p_ch_p_myr <- poss_ch <- c();
for (gh in 1:notu)	{
	poss_ch <- c(poss_ch,sum(!chmatrix[gh,] %in% c(UNKNOWN,INAP)));
		if (ghosts[gh]>0)	{
		st_p_ch_p_myr <- c(st_p_ch_p_myr,((steps[gh]/poss_ch[gh])/ghosts[gh]));
		}
	}
init_rate <- data.frame(mean=as.numeric((sum(steps[ghosts>0])/sum(poss_ch))/sum(ghosts)),median=as.numeric(median(st_p_ch_p_myr)),stringsAsFactors = F);
return(init_rate);
}

Mk_continuous <- function(bd,k,alpha)	{
# bd: time (branch duration)
# k: n_states
# alpha: instantaneous rate
# from Lewis 2001
pstasis <- (1/k)+((k-1)*exp(-k*alpha*bd)/k)
pchange <- (1/k)-(exp(-k*alpha*bd)/k)
return(c(pstasis,pchange))
}

Mk_punctuated_background <- function(bd,k,epsilon,lambda)	{
# bd: time
# k: n_states
# epsilon: pulsed rate
# lambda: speciation rate
# from Wagner & Marcot 2001
pstasis <- (1/k)+((k-1)*exp(-k*epsilon*(lambda*bd))/k)
pchange <- (1/k)-(exp(-k*epsilon*(lambda*bd))/k)
return(c(pstasis,pchange))
}

Mk_continuous_and_punctuated <- function(k,alpha,epsilon,lambda,m,bd)	{
pstasis <- (1/k)+((k-1)*exp(-k*(epsilon*(m+lambda*bd)+(alpha*bd)))/k)
pchange <- (1/k)-(exp(-k*(epsilon*(m+lambda*bd)+(alpha*bd)))/k)
return(c(pstasis,pchange))
}

Mk_punctuated <- function(bd,k,epsilon,lambda,m)	{
# bd: time over which speciation & change can accrue
# k: number of n_states
# epsilon: probability of deriving a state at branching
# lambda: origination rate
# m: minimum number of branches (0, 1 or 0.5 for there had to be one here or for the sister-taxon; 0.67 for trichotomy)
# modified from Wagner & Marcot 2001
# initial punctuation follows a binomial
pstasis1 <- 1-((k-1)*epsilon)
pchange1 <- epsilon
# subsequent change follows a Poisson with lambda*bd subbing for time
pstasis2 <- (1/k)+((k-1)*exp(-k*epsilon*lambda*bd)/k)
pchange2 <- (1/k)-(exp(-k*epsilon*lambda*bd)/k)
pstasisA <- pchangeA <- 0
if (m>0)	{
	# if there might have been branching, then get binomial change/stasis followed by Poisson change/stasis
	pstasisA <- (pstasis1*pstasis2)+((k-1)*(pchange1*pchange2))
	pchangeA <- (pchange1*pstasis2)+(pstasis1*pchange2)+((k-2)*pchange1*pchange2)
	}
pstasis <- (m*pstasisA)+((1-m)*pstasis2)
pchange <- (m*pchangeA)+((1-m)*pchange2)
return(c(pstasis,pchange))
}

Mk_punctuated_plus_continuous <- function(bd,k,alpha,epsilon,lambda,m)	{
# bd: time over which speciation & change can accrue
# k: number of n_states
# alpha: instantaneous rate of deriving a state
# epsilon: probability of deriving a state at branching
# lambda: origination rate
# m: minimum number of branches (0, 1 or 0.5 for there had to be one here or for the sister-taxon; 0.67 for trichotomy)
# modified from Wagner & Marcot 2001
# initial punctuation follows a binomial
pstasis1 <- 1-((k-1)*epsilon)
pchange1 <- epsilon
# subsequent change follows a Poisson with lambda*bd subbing for time
pstasis2 <- (1/k)+((k-1)*exp(-k*(alpha+(epsilon*lambda))*bd)/k)
pchange2 <- (1/k)-(exp(-k*(alpha+(epsilon*lambda))*bd)/k)
#-k*(alpha+(epsilon*lambda))*bd
#-k*((alpha*bd)+(epsilon*lambda*bd))
pstasisA <- pchangeA <- 0
if (m>0)	{
	# if there might have been branching, then get binomial change/stasis followed by Poisson change/stasis
	pstasisA <- (pstasis1*pstasis2)+((k-1)*(pchange1*pchange2))
	pchangeA <- ((pchange1*pstasis2)+(pstasis1*pchange2)+((k-2)*pchange1*pchange2))
	}
pstasis <- (m*pstasisA)+((1-m)*pstasis2)
pchange <- (m*pchangeA)+((1-m)*pchange2)
return(c(pstasis,pchange))
}

accersi_tree_log_likelihood_simple <- function(m_alpha,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#print(m_alpha);  this is the good one!
alpha_sts <- vector(length=max(n_states))
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
ttu <- max(vector_tree)
nNodes <- ttu-notu
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
#nn <- nNodes+1;
base <- notu+1;
for (nn in nNodes:1)	{
#	nn <- nn-1;
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k];
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- exp(log(mk[1,d])+log(marginals[dd,ch,st]));
					# prob net change
					mg_ch <- exp(log(mk[2,d])+log(marginals[dd,ch,(1:k)[!(1:k) %in% st]]));
					marginals[htu,ch,st] <- exp(log(marginals[htu,ch,st])+log(mg_st+sum(mg_ch)));
#					if (ch==1) print(marginals[htu,ch,]);
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,]);
return(sum(log(basal_char_likelihoods)))
}

accersi_tree_log_likelihood_simple_dud <- function(m_alpha,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#print(m_alpha);
alpha_sts <- vector(length=max(n_states));
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2);
ttu <- max(vector_tree);
nNodes <- ttu-notu;
mxstates <- max(n_states);
tnchars <- length(n_states);
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
#nn <- nNodes+1;
base <- notu+1;
for (nn in nNodes:1)	{
	nn <- nn-1;
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k];
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- exp(log(mk[1,d])+log(marginals[dd,ch,st]));
					# prob net change
					mg_ch <- exp(log(mk[2,d])+log(marginals[dd,ch,(1:k)[!(1:k) %in% st]]));
					marginals[htu,ch,st] <- exp(log(marginals[htu,ch,st])+log(mg_st+sum(mg_ch)));
#					if (ch==1) print(marginals[htu,ch,]);
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k]);
			}
		}
#	marginals[,1,];
	}
basal_char_likelihoods <- rowSums(marginals[base,,]);
return(sum(log(basal_char_likelihoods)))
}

accersi_tree_marginals_simple <- function(m_alpha,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#
alpha_sts <- vector(length=max(n_states))
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
ttu <- max(vector_tree)
nNodes <- ttu-notu
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
#nn <- nNodes+1;
base <- notu+1;
for (nn in nNodes:1)	{
#	nn <- nn-1;
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k];
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- exp(log(mk[1,d])+log(marginals[dd,ch,st]));
					# prob net change
					mg_ch <- exp(log(mk[2,d])+log(marginals[dd,ch,(1:k)[!(1:k) %in% st]]));
					marginals[htu,ch,st] <- exp(log(marginals[htu,ch,st])+log(mg_st+sum(mg_ch)));
#					if (ch==1) print(marginals[htu,ch,]);
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
return(log(basal_char_likelihoods));
}

accersi_tree_bds_probability_simple <- function(m_alpha,divergence_times,vector_tree,notu,marginals,n_states,BDS)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#
alpha_sts <- vector(length=max(n_states))
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
ttu <- max(vector_tree)
nNodes <- ttu-notu
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
nn <- nNodes
base <- notu+1
gap_lgprobs <- 0.0
for (nn in nNodes:1)	{
	htu <- nn+notu
	daughters <- (1:ttu)[vector_tree==htu]
	for (d in 1:length(daughters))	{
		stgap <- (1:nrow(BDS))[BDS$ma > divergence_times[htu]]
		stgap <- stgap[BDS$ma[stgap] < divergence_times[daughters[d]]]
		gap_lgprobs <- gap_lgprobs+sum(log(BDS$pgap[stgap]))
		}
	tomy <- length(daughters)
	bd <- abs(divergence_times[daughters]-divergence_times[htu])
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k]
		rchars <- (1:tnchars)[n_states==k]
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
log_prob <- sum(log(basal_char_likelihoods))+gap_lgprobs
return(log_prob)
}

accersi_tree_bds_probability_given_basal_divergence_and_rate <- function(divergence_times,onset,m_alpha,vector_tree,notu,marginals,n_states,BDS)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#
alpha_sts <- vector(length=max(n_states))
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
ttu <- max(vector_tree)
nNodes <- ttu-notu
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
nn <- nNodes
base <- notu+1
gap_lgprobs <- 0.0
for (nn in nNodes:1)	{
	htu <- nn+notu
	daughters <- (1:ttu)[vector_tree==htu]
	for (d in 1:length(daughters))	{
		stgap <- (1:nrow(BDS))[BDS$ma > divergence_times[htu]]
		stgap <- stgap[BDS$ma[stgap] < divergence_times[daughters[d]]]
		gap_lgprobs <- gap_lgprobs+sum(log(BDS$pgap[stgap]))
		}
	tomy <- length(daughters)
	bd <- abs(divergence_times[daughters]-divergence_times[htu])
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k]
		rchars <- (1:tnchars)[n_states==k]
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
log_prob <- sum(log(basal_char_likelihoods))+gap_lgprobs
return(log_prob)
}

accersi_tree_log_likelihood_simple_check_nodes <- function(node_ages,m_alpha,vector_tree,notu,fas,marginals,n_states)	{
# starting_node_ages: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
ttu <- max(vector_tree)
nNodes <- ttu-notu
changes <- 1
while (changes>0)	{
	changes <- 0
	for (n in nNodes:2)	{
		htu <- n+notu
		anc <- vector_tree[htu]-notu
		if (node_ages[anc] > node_ages[n])	{
			da <- node_ages[n]
			node_ages[n] <- node_ages[anc]
			node_ages[anc] <- da
			changes <- changes+1
			}
		}
	}
starting_node_ages <- c(fas,node_ages)
for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
nn <- nNodes
base <- notu+1
for (nn in nNodes:1)	{
	htu <- nn+notu
	daughters <- (1:ttu)[vector_tree==htu]
	tomy <- length(daughters)
	bd <- abs(starting_node_ages[daughters]-starting_node_ages[htu])
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k]
		rchars <- (1:tnchars)[n_states==k]
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
return(sum(log(basal_char_likelihoods)))
}

# do another version of this in which we optimize the rate!
accersi_tree_likelihood_divergences_scaled_to_base <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving onset of branches
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
first_sample <- min(divergence_times[1:notu])
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree)
nNodes <- ttu-notu
predates <- (1:ttu)[divergence_times<first_sample]
orig_base <- divergence_times[ur_b]
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))

for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha/(1+st-2)
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
nn <- nNodes
base <- notu+1
for (nn in nNodes:1)	{
	htu <- nn+notu
	daughters <- (1:ttu)[vector_tree==htu]
	tomy <- length(daughters)
	bd <- abs(divergence_times[daughters]-divergence_times[htu])
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k]
		rchars <- (1:tnchars)[n_states==k]
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
return(sum(log(basal_char_likelihoods)))
}

# do another version of this in which we optimize the rate!
accersi_tree_likelihood_divergences_scaled_to_base_w_two_rates <- function(base,m_alpha_1,m_alpha_2,branch_rates,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
first_sample <- min(divergence_times[1:notu])
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree)
nNodes <- ttu-notu
predates <- (1:ttu)[divergence_times<first_sample]
orig_base <- divergence_times[ur_b]
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))

for (st in 2:max(n_states))	alpha_sts[st] <- m_alpha_1/(1+st-2)
mxstates <- max(n_states)
tnchars <- length(n_states)
#marginals <- array(0,dim=c(ttu,tnchars,mxstates))
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates))
nn <- nNodes
base <- notu+1
for (nn in nNodes:1)	{
	htu <- nn+notu
	daughters <- (1:ttu)[vector_tree==htu]
	tomy <- length(daughters)
	bd <- abs(divergence_times[daughters]-divergence_times[htu])
	for (k in 2:mxstates)	{
		alpha <- alpha_sts[k]
		rchars <- (1:tnchars)[n_states==k]
		mk <- sapply(bd,Mk_continuous,k,alpha)	# mk model for each node, with first value being P[net stasis]
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[base,,])
return(sum(log(basal_char_likelihoods)))
}

# m_alpha <- 0.01307749
# m_alpha <- 0.01082389
# marginals <- allocate_node_marginals(chmatrix_2,nttu,n_states)
optimo_rate_divergences_scaled_to_base_old <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree)
nNodes <- ttu-notu
predates <- (1:ttu)[divergence_times<first_sample]
orig_base <- divergence_times[ur_b]
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
cl <- list(fnscale=-1);
ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
ml_alpha <- ddd$par;
mlgln <- ddd$value
results <- list(ml_alpha,mlgln,divergence_times)
names(results) <- c("ml_Rate","lnL","divergence_times")
return(results)
}

#optimo_rate_and_tree_divergences_posterior_scaled_to_base
optimo_rate_and_tree_divergences_posterior_scaled_to_base <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times)
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times")
return(results)
}

#base=base_rep[bb];m_alpha_1=m_alpha_bang_rep;m_alpha_2=m_alpha_post_rep;divergence_times=init_divergence_times_rep;branch_rates=branch_rates_rep;vector_tree=vector_tree_ord;notu=notu_ord;marginals=marginals;n_states=n_states;BDS=BDS;rate_quants=rate_quants
optimo_rate_and_tree_divergences_posterior_scaled_to_base_two_rates_new <- function(base,m_alpha_1,m_alpha_2,divergence_times,branch_rates,vector_tree,notu,marginals,n_states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
#ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
m_alphas <- c(m_alpha_1,m_alpha_2);
ddd <- optim(m_alphas,fn=likelihood_of_early_and_late_alphas_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,branch_rates=branch_rates,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times);
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times");
return(results)
}

#rate_shift <- bang_end
# written 2020-07-25
optimo_rate_and_tree_divergences_posterior_scaled_to_base_rate_shift <- function(base,m_alpha_1,m_alpha_2,divergence_times,rate_shift,vector_tree,notu,marginals,n_states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
#ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
m_alphas <- c(m_alpha_1,m_alpha_2);
ddd <- optim(m_alphas,fn=likelihood_of_rate_shift_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,rate_shift=rate_shift,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
#ddd <- optim(m_alphas,fn=likelihood_of_early_and_late_alphas_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,branch_rates=branch_rates,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times);
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times");
return(results)
}

optimo_rate_and_tree_divergences_posterior_scaled_to_base_two_rates <- function(base,m_alpha_1,m_alpha_2,divergence_times,branch_rates,vector_tree,notu,marginals,n_states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
#ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
m_alphas <- c(m_alpha_1,m_alpha_2);
ddd <- optim(m_alphas,fn=likelihood_of_early_and_late_alphas_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,branch_rates=branch_rates,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times);
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times");
return(results)
}

optimo_rate_and_tree_divergences_likelihood_scaled_to_base <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree)
nNodes <- ttu-notu
predates <- (1:ttu)[divergence_times<first_sample]
orig_base <- divergence_times[ur_b]
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,marginals=marginals,n_states=n_states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value
results <- list(ml_alpha,mlgln,divergence_times)
names(results) <- c("ml_Rate","lnL","divergence_times")
return(results)
}

optimo_rate_divergences_scaled_to_base_bds_redone <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#original_divergence_times <- divergence_times;
#original_basal_divergence <- min(original_divergence_times);
latest <- max(-abs(divergence_times));
node_divergences <- divergences[(notu+2):length(vector_tree)];
best_nodal_divergences <- optim(par=node_divergences,
						  fn=accersi_tree_bds_probability_given_divergences_and_rate,
						  m_alpha=m_alpha,
						  m_alpha=m_alpha,
						  onset=base,
						  vector_tree=vector_tree,
						  fas=fas,
						  marginals=marginals,
						  n_states=n_states,
						  BDS=BDS,
						  min=base,max=latest);

first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))

### change routine here: redate all of the branches to get the best overall rate....
(abs(base-first_sample)/abs(orig_base-first_sample))


#sapply(m_alpha,accersi_tree_log_likelihood_simple,)								divergence_times,                 vector_tree,      notu,     marginals,          n_states,       BDS)
ddd <- optim(m_alpha,fn=accersi_tree_bds_probability_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,BDS=BDS,lower=m_alpha/10,upper=1.0,control=cl)
ml_alpha <- ddd$par
mlgln <- ddd$value
results <- list(ml_alpha,mlgln)
names(results) <- c("mp_Rate","lnP")
return(results)
}

optimo_rate_divergences_scaled_to_base_bds <- function(base,m_alpha,divergence_times,vector_tree,notu,marginals,n_states,BDS)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#original_divergence_times <- divergence_times;
#original_basal_divergence <- min(original_divergence_times);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))

### change routine here: redate all of the branches to get the best overall rate....
#(abs(base-first_sample)/abs(orig_base-first_sample))
#sapply(m_alpha,accersi_tree_log_likelihood_simple,)								divergence_times,                 vector_tree,      notu,     marginals,          n_states,       BDS)
ddd <- optim(m_alpha,fn=accersi_tree_bds_probability_simple,method="L-BFGS-B",divergence_times=divergence_times,vector_tree=vector_tree,notu=notu,marginals=marginals,n_states=n_states,BDS=BDS,lower=m_alpha/10,upper=1.0,control=cl)
ml_alpha <- ddd$par
mlgln <- ddd$value
results <- list(ml_alpha,mlgln)
names(results) <- c("mp_Rate","lnP")
return(results)
}

# rescale vector giving divergence time of each taxon & node
rescale_early_divergences <- function(base,divergence_times,vector_tree,notu)	{
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
first_sample <- min(divergence_times[1:notu])
ur_b <- notu+1	# htu number of basal node
ttu <- max(vector_tree)
nNodes <- ttu-notu
predates <- (1:ttu)[divergence_times<first_sample]
orig_base <- divergence_times[ur_b]
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample))
return(divergence_times)
}

# rescale array giving onset & end time of each "ghost" branch on a phylogeny
rescale_basal_divergences <- function(new_base,divergence_dates,highest_base)	{
if (!is.data.frame(divergence_dates))
	divergence_dates <- data.frame(onset=as.numeric(divergence_dates[,1]),end=as.numeric(divergence_dates[,2]),stringsAsFactors = F);
colnames(divergence_dates) <- c("onset","end");
divergence_dates <- -abs(divergence_dates);
old_base <- min(divergence_dates$onset);
old_span <- highest_base-old_base;
new_span <- highest_base-new_base;
for (dd in 1:nrow(divergence_dates))	{
	if (divergence_dates$onset[dd]<highest_base)
		divergence_dates$onset[dd] <- highest_base-new_span*((highest_base-divergence_dates$onset[dd])/old_span);
	if (divergence_dates$end[dd]<highest_base)
		divergence_dates$end[dd] <- highest_base-new_span*((highest_base-divergence_dates$end[dd])/old_span);
	}
#sum(divergence_dates$end-divergence_dates$onset)
return(divergence_dates);
}

# likelihood routines
#t1 <- 1
#t2 <- 4
#k <- 2
#alpha <- 0.05
Mk_continuous_species_stasis <- function(k,alpha,t1,t2)	{
#stasis_marg <- Mk_continuous(k,alpha,bd=t1)*Mk_continuous(k,alpha,bd=t2)
return(Mk_continuous(k,alpha,bd=abs(t1-t2))[1])
#stasis_prob <- stasis_marg/sum(stasis_marg)
#return(stasis_prob)
}

Mk_continuous_varying_rates <- function(k,alphas_c,bd)	{
bins <- length(alphas_c)
rate <- 0
for (b in 1:bins)
	rate <- rate+(alphas_c[b]*bd[b])

pstasis <- (1/k)+((k-1)*exp(-k*rate)/k)
pchange <- (1/k)-(exp(-k*rate)/k)
return(pstasis,pchange)
}

Mk_punctuated_varying_rates <- function(k,epsilons,origination_rates,m,bin_durations)	{
nbins <- length(epsilons);
rate <- 0;
for (b in 1:nbins)	{;
	if (b==1)	{;
		mb <- m	# minimum number of branchings happened at outset	
		}	else mb <- 0
	rate <- rate+(epsilons[b]*(mb+origination_rates[b]*bin_durations[b]))
	}

pstasis <- (1/k)+((k-1)*exp(-k*rate)/k)
pchange <- (1/k)-(exp(-k*rate)/k)
return(pstasis,pchange)
}

Mk_mixed_varying_rates <- function(k,alphas_c,alphas_p,lambdas,m,bd)	{
bins <- length(alphas_p)
rate <- 0
for (b in 1:bins)	{
	if (b==1)	{
		mb <- m	# minimum number of branchings happened at outset	
		}	else mb <- 0
	rate <- rate+(alphas_c[b]*bd[b])+(alphas_p[b]*(mb+labdas[b]*bd[b]))
	}

pstasis <- (1/k)+((k-1)*exp(-k*rate)/k)
pchange <- (1/k)-(exp(-k*rate)/k)
return(pstasis,pchange)
}

get_stages_and_rates_given_onset_and_end <- function(chronology,rates,onset,end)	{
# if just one rate, then use rates <- matrix(0,length(x),1) & rates[,1] <- x
b1 <- match(FALSE,chronology>=onset)-1
b2 <- match(FALSE,chronology>end)
return(rates[b1:b2,])
}

accersi_branch_likelihoods_continuous <- function(desc_states,anc_states,marginals,k,alphas_c,bd)	{
loutc <- Mk_continuous_varying_rates(k,alphas_c,bd)
for (as in 1:k)	{
	for (ds in 1:k)	{
		if (as==ds)	{
			marginals[as] <- marginals[as]*loutc[1]
			} else {
			marginals[as] <- marginals[as]*loutc[2]
			}
		}
	}
return(marginals)
}

## send character
accersi_branch_likelihoods_punctuated <- function(desc_states,anc_states,char_marginals,k,alphas_c,alphas_p,lambdas,m,bd)	{
if (max(alphas_c)==0)	{
	loutc <- Mk_punctuated_varying_rates(k,alphas_p,lambdas,m,bd)
	}	else {
	loutc <- Mk_mixed_varying_rates(k,alphas_c,alphas_p,lambdas,m,bd)
	}
for (as in 1:k)	{
	for (ds in 1:k)	{
		if (as==ds)	{
			char_marginals[as] <- char_marginals[as]*loutc[1]
			}	else {
			char_marginals[as] <- char_marginals[as]*loutc[2]
			}
		}
	}
return(marginals)
}

get_node_likelihood <- function(descendants,ancestor,marginals,otu,strat_info,mode)	{
return(descendants)
}

allocate_node_marginals <- function(chmatrix,nttu,n_states)	{
# chmatrix: original character matrix;
# nttu: total taxa + nodes;
# n_states: n_states per character
mxstates <- max(n_states);
tnchars <- ncol(chmatrix);
notu <- nrow(chmatrix);
marginals <- array(0,dim=c(nttu,tnchars,mxstates));
n_states[n_states<2] <- 2;
for (n in 1:notu)	{
	for (c in 1:tnchars)	{
		if (chmatrix[n,c]>=0)	{
			marginals[n,c,chmatrix[n,c]] <- 1
			} else if (chmatrix[n,c]==UNKNOWN || chmatrix[n,c]==INAP)	{
			marginals[n,c,1:n_states[c]]	<- 1/n_states[c];
			} else {
			pc <- -1*chmatrix[n,c];
			polys <- vector(length=ceiling(log10(pc+0.5)));
			ex <- length(polys);
			for (i in 1:length(polys))	{
				polys[i] <- floor(pc/(10^(ex-1)));
				pc <- pc-(polys[i]*10^(ex-1));
				ex <- ex-1;
				}
#			polys <- polys+1;
			marginals[n,c,polys] <- 1/length(polys)
			}
		}
	}
return(marginals)
}

accersi_matrix_tree_from_ape_tree <- function(ape)	{
### ape is the output of read.nexus
#	ape$edge gives apes' version of the tree, which alters the numbers
#	This changes them back to the numbers in the original (1,(2,3)) format
#	The information for doing this is in ape$tip.label
base <- ape$edge[1,1]
notu <- base-1	# the lowest htu number is 1 more than maximum otu number
yyy <- ape$edge[order(ape$edge[,2]),]
yyy[1:notu,2] <- as.numeric(ape$tip.label)
zzz <- yyy[order(yyy[,2]),]
vector_tree <- zzz[order(zzz[,2]),1]
vector_tree <- c(vector_tree[1:(base-1)],-1,vector_tree[base:length(vector_tree)])
mtree <- transform_vector_tree_to_matrix_tree(vector_tree)
}

accersi_starter_divergences_min <- function(vector_tree,simple_clock,strat_ranges,apos,prec=0.1)	{
# vector_tree: vector giving htu from which each otu & htu evovled
# simple_clock: vector giving expected branch durations given overall rates
# strat_ranges: first & last appearances of species
if (is.matrix(strat_ranges))
	strat_ranges <- data.frame(strat_ranges,stringsAsFactors = F);
if (sum(colnames(strat_ranges) %in% c("fa","FA","fad","FAD","onset"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("fa","FA","fad","FAD","onset")];
	col_titles[fix_col] <- "ma_lb";
	colnames(strat_ranges) <- col_titles;
	}
if (sum(colnames(strat_ranges) %in% c("la","LA","Lad","LAD","end"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("la","LA","Lad","LAD","end")];
	col_titles[fix_col] <- "ma_ub";
	colnames(strat_ranges) <- col_titles;
	}

mtree_to_date <- transform_vector_tree_to_matrix_tree(vector_tree=vector_tree);
notu <- min(vector_tree[vector_tree>0])-1;
poss_anc <- (1:notu)[apos[1:notu]==0];
nttu <- max(vector_tree)
nNodes <- nttu-notu
node_ages_starters <- c(-abs(strat_ranges$ma_lb[1:notu]),rep(0,nNodes));
nn <- nNodes;
#while (nn>0)	{
earliest_onset <- -max(abs(strat_ranges$ma_lb))-prec;
for (nn in nNodes:1)	{
	htu <- notu+nn;							# get vector number for the node's age
	daughters <- (1:nttu)[vector_tree==htu];		# get descendants of this node
#	if (max(simple_clock[daughters])==0)
#		simple_clock[daughters[length(daughters)]] <-prec;
	local_clock <- node_ages_starters[daughters]-simple_clock[daughters];	# get divergence times implied by simple clock
	# make sure that local clock is at least as old minimum divergence times
	local_clock[local_clock>min(node_ages_starters[daughters])] <- min(node_ages_starters[daughters])
	node_ages_starters[htu] <- prec*round(mean(local_clock)/prec,0);
#	if (min(strat_ranges$ma_lb[daughters])<node_ages_starters[htu])
#		node_ages_starters[htu] <- min(strat_ranges$ma_lb[daughters])-0.1;
#	nn <- nn-1
#	daughters
#	node_ages_starters[daughters]
#	node_ages_starters
	}
return(node_ages_starters)
}

accersi_starter_divergences <- function(vector_tree,simple_clock,strat_ranges,prec=0.1)	{
# vector_tree: vector giving htu from which each otu & htu evovled
# simple_clock: vector giving expected branch durations given overall rates
# strat_ranges: first & last appearances of species
if (is.matrix(strat_ranges))
	strat_ranges <- data.frame(strat_ranges,stringsAsFactors = F);
if (sum(colnames(strat_ranges) %in% c("fa","FA","fad","FAD","Ma_LB","onset"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("fa","FA","fad","FAD","Ma_LB","onset")];
	col_titles[fix_col] <- "ma_lb";
	colnames(strat_ranges) <- col_titles;
	}
if (sum(colnames(strat_ranges) %in% c("la","LA","Lad","LAD","Ma_UB","end"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("la","LA","Lad","LAD","Ma_UB","end")];
	col_titles[fix_col] <- "ma_ub";
	colnames(strat_ranges) <- col_titles;
	}

notu <- min(vector_tree[vector_tree>0])-1;
nttu <- max(vector_tree);
nNodes <- nttu-notu;
node_ages_starters <- c(-abs(strat_ranges$ma_lb[1:notu]),rep(0,nNodes));
nn <- nNodes;
#while (nn>0)	{
for (nn in nNodes:1)	{
	htu <- notu+nn;							# get vector number for the node's age
	daughters <- (1:nttu)[vector_tree==htu];		# get descendants of this node
	if (max(simple_clock[daughters])==0)
		simple_clock[daughters[length(daughters)]] <- prec;
	local_clock <- node_ages_starters[daughters]-simple_clock[daughters];	# get divergence times implied by simple clock
	# make sure that local clock is at least as old minimum divergence times
	local_clock[local_clock>min(node_ages_starters[daughters])] <- min(node_ages_starters[daughters]);
	node_ages_starters[htu] <- prec*round(mean(local_clock)/prec,0);
#	nn <- nn-1
#	daughters
#	node_ages_starters[daughters]
#	node_ages_starters
	}
return(node_ages_starters);
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### basic "stratocladistic" routines ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# get minimum implied gap along each branch
accersi_gaps_from_vector_tree <- function(vector_tree,strat_ranges,apos=NULL)	{
if (is.matrix(strat_ranges))	strat_ranges <- data.frame(strat_ranges,stringsAsFactors = F);
if (ncol(strat_ranges)==2)	{
	colnames(strat_ranges) <- c("FA","LA");
	} else if (ncol(strat_ranges)==3)	{
	colnames(strat_ranges) <- c("taxon","FA","LA");
	}
ohtu <- length(vector_tree);
notu <- nrow(strat_ranges);
if (is.null(apos))
	apos <- rep(1,ohtu);
base <- notu+1;
dates <- date_clades_from_vector_tree(vector_tree,FA=strat_ranges$FA);
gaps <- vector(length=ohtu);
node_poss_anc <- accersi_poss_ancestors_for_nodes(vector_tree,FA=strat_ranges$FA,apos);

for (tx in 1:ohtu)	{
	if (tx==base)	tx <- tx+1
	if (node_poss_anc[vector_tree[tx]]==0)	{
		gaps[tx] <- abs(dates[tx]-dates[vector_tree[tx]])
		} else	{
		gaps[tx] <- max(0,dates[tx]-strat_ranges$LA[node_poss_anc[vector_tree[tx]]])
		}
	}
return(gaps);
}

# get possible ancestors given apomorphies
accersi_poss_ancestors_for_nodes <- function(vector_tree,FA,apos)	{
ohtu <- length(vector_tree);
#notu <- length(FA);
notu <- match(-1,vector_tree)-1;
node_poss_anc <- rep(0,ohtu);
if (!is.na(apos[1]))
	for (n in 1:notu)
		if (apos[n]==0)	
			if (node_poss_anc[vector_tree[n]]==0)	{
				node_poss_anc[vector_tree[n]] <- n;
				} else if (FA[node_poss_anc[vector_tree[n]]] > FA[n])
				FA[node_poss_anc[vector_tree[n]]] <- FA[n];
return(node_poss_anc)
}

# get node ages using minimum necessary ages
date_taxa_on_tree_simple <- function(vector_tree,FAs)	{
ttus <- max(vector_tree)	# total htus+otus
dates <- vector(length=ttus)
notu <- length(FAs)	# number of otus
if (length(dim(as.array(vector_tree)))==1)	{
	# routine if tree given as a single vector with ancestral
	base <- notu+1
	end <- max(FAs)
	for (nd in base:ttus)	dates[nd] <- end
	for (sp in 1:notu)	{
		if (vector_tree[sp]>0)	{
			dates[sp] <- FAs[sp]
			anc <- vector_tree[sp]
			if (dates[anc]>dates[sp])	dates[anc] <- dates[sp]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	for (nd in ttus:(base+1))	{
		if (vector_tree[nd]>0)	{
			anc <- vector_tree[nd]
			if (dates[anc]>dates[nd])	dates[anc] <- dates[nd]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	} else	{
	Nnode <- dim(vector_tree)[1]
	dates[1:notu] <- FAs
	for (n in Nnode:1)	{
		f1 <- sum(vector_tree[n,]>0)
		ht <- n+notu
		dates[ht] <- min(dates[vector_tree[n,(1:f1)]])
		}
	}
return(dates)
}

# collect node ages with plotting in mind
date_taxa_on_tree_for_plotting <- function(vector_tree,FAs,tie_break=1/3)	{
ttus <- max(vector_tree)	# total htus+otus
dates <- vector(length=ttus)
notu <- length(FAs)	# number of otus
if (length(dim(as.array(vector_tree)))==1)	{
	# routine if tree given as a single vector with ancestral
	base <- notu+1
	end <- max(FAs)
	for (nd in base:ttus)	dates[nd] <- end
	
	for (sp in 1:notu)	{
		if (vector_tree[sp]>0)	{
			dates[sp] <- FAs[sp]
			anc <- vector_tree[sp]
			if (dates[anc]>dates[sp])	dates[anc] <- dates[sp]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	for (nd in ttus:(base+1))	{
		if (vector_tree[nd]>0)	{
			anc <- vector_tree[nd]
			if (dates[anc]>dates[nd])	dates[anc] <- dates[nd]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	grandparents <- sort(unique(vector_tree[(base+1):ttus]),decreasing=TRUE)
	for (g in 1:length(grandparents))	{
		htu <- grandparents[g]
		dn <- ((base+1):ttus)[vector_tree[(base+1):ttus]==htu]
		if (length(dn)>0)	{
			if (min(dates[dn]) <= dates[htu])	dates[htu] <- min(dates[dn])-tie_break
			}
		}
	} else	{
	Nnode <- dim(vector_tree)[1]
	dates[1:notu] <- FAs
	for (n in Nnode:1)	{
		f1 <- sum(vector_tree[n,]>0)
		ht <- n+notu
		dates[ht] <- min(dates[vector_tree[n,(1:f1)]])
		}
	}
return(dates)
}

accersi_node_durations <- function(vector_tree,FAs)	{
divergences <- date_clades_from_vector_tree(vector_tree,FAs);
notu <- min(vector_tree[vector_tree>0])-1;
htu <- -sort(unique(-vector_tree[vector_tree>0]));

node_durations <- c();
for (ht in 1:length(htu))	{
	if (vector_tree[htu[ht]]>0)	{
		node_durations <- rbind(c(divergences[vector_tree[htu[ht]]],divergences[htu[ht]]),node_durations);
		} else	{
		node_durations <- rbind(c(divergences[htu[ht]],divergences[htu[ht]]),node_durations);
		}
	}
node_durations <- data.frame(fa=as.numeric(node_durations[,1]),la=as.numeric(node_durations[,2]),stringsAsFactors = F);
rownames(node_durations) <- paste(rep("htu_",nrow(node_durations)),sort(htu),sep="");
return(node_durations)
}

# collect node ages from vector tree
date_clades_from_vector_tree <- function(vector_tree,FAs)	{
ohtu <- length(vector_tree)
notu <- length(FAs)
base <- notu+1
dates <- vector(length=ohtu)
end <- max(FAs)
for (nd in base:ohtu)	dates[nd] <- end
for (sp in 1:notu)	{
	if (vector_tree[sp]>0)	{
		dates[sp] <- FAs[sp]
		anc <- vector_tree[sp]
		if (dates[anc]>dates[sp])	dates[anc] <- dates[sp]
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
for (nd in ohtu:(base+1))	{
	if (vector_tree[nd]>0)	{
		anc <- vector_tree[nd]
		if (dates[anc]>dates[nd])	dates[anc] <- dates[nd]
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
return(dates)
}

date_clades_from_vector_tree_with_constraints <- function(vector_tree,FAs,constraints)	{
ohtu <- length(vector_tree);
notu <- length(FAs);
base <- notu+1;
htu <- base:max(vector_tree);
dates <- vector(length=ohtu);
end <- max(FAs);
for (nd in base:ohtu)	dates[nd] <- end;
dates[htu[constraints!=0]] <- constraints[constraints!=0];
for (sp in 1:notu)	{
	if (vector_tree[sp]>0)	{
		dates[sp] <- FAs[sp];
		anc <- vector_tree[sp];
		if (dates[anc]>dates[sp])	dates[anc] <- dates[sp];
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
for (nd in ohtu:(base+1))	{
	if (vector_tree[nd]>0)	{
		anc <- vector_tree[nd];
		if (dates[anc]>dates[nd])	dates[anc] <- dates[nd];
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
return(dates)
}

date_notu_divergences_from_vector_tree_with_constraints <- function(vector_tree,FAs,constraints,apomorphies)	{
ohtu <- length(vector_tree);
notu <- length(FAs);
nNodes <- ohtu-notu;
base <- notu+1;
#htu <- base:max(vector_tree);
dates <- date_clades_from_vector_tree_with_constraints(vector_tree,FAs,constraints);
sampled_node <- vector(length=nNodes);
obs_ancestors <- (1:notu)[apomorphies==0];
for (nd in 1:nNodes)	{
	htu <- nd+notu;
	f1 <- (1:ohtu)[vector_tree==htu];
	if (sum(f1 %in% obs_ancestors)>0)	{
		sampled_node[nd] <- 1;
		} else	{
		dates[f1] <- max(dates[f1]);
		}
	}
return(dates)
}

date_onsets_and_ends_from_vector_tree_with_constraints <- function(vector_tree,FAs,constraints,apomorphies,min_dur=1)	{
ohtu <- length(vector_tree);
notu <- length(FAs);
nNodes <- ohtu-notu;
base <- notu+1;
#htu <- base:max(vector_tree);
dates <- date_clades_from_vector_tree_with_constraints(vector_tree,FAs,constraints);
sampled_node <- vector(length=nNodes);
obs_ancestors <- (1:notu)[apomorphies==0];
dates2 <- c(FAs,rep(0,nNodes));
for (nd in nNodes:1)	{
	htu <- nd+notu;
	f1 <- (1:ohtu)[vector_tree==htu];
	if (sum(f1 %in% obs_ancestors)>0)	{
		sampled_node[nd] <- 1;
		} else	{
		dates[f1] <- dates[htu];
		}
	dates2[htu] <- max(dates[f1]);
	if (dates2[htu]==dates[htu])	dates[htu] <- dates[htu]-abs(min_dur);
	}
durations <- data.frame(onset=as.numeric(dates),end=as.numeric(dates2));
return(durations)
}

# routine to get oldest member within each node 
accersi_youngest_daughter_in_node <- function(tree,FAs)	{
notu <- length(FAs)
oldest_daughter <- date_nodes(tree,FAs)
if (is.matrix(tree))	{
	Nnode <- dim(tree)[1]
	smiths <- vector(length=Nnode)
	for (n in Nnode:1)	{
		if (sum(m_tree[n,]<=notu)==length(m_tree[n,]))	{
			smiths[n] <- max(FAs[m_tree[n,]])
			}	else	{
			dc <- m_tree[n,m_tree[n,]>notu]-notu
			smiths[n] <- min(oldest_daughter[dc])
			}
		}
	}	else	{
	ttu <- length(tree)
	Nnode <- ttu-notu
	smiths <- rep(min(FAs),Nnode)
	for (n in Nnode:1)	{
		htu <- n+notu
		f1 <- (1:ttu)[tree==htu]
		if (sum(f1<=notu)==length(f1))	{
			smiths[n] <- max(FAs[f1])
			}	else	{
			dc <- f1[f1>notu]-notu
			smiths[n] <- min(oldest_daughter[dc])
			}
		}
	}
return(smiths)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Mk Model Fun (including Q-Matrix Mania() ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

identity_matrix <- function(n)	{
idm <- array(0,c(n,n))
for (i in 1:n)	idm[i,i] <- 1
return(idm)
}

multiply_square_matrices <- function(B1,B2)	{
B1 %*% B2
n <- nrow(B2)
BB <- array(0,c(n,n))
for (i in 1:n)	{
	for (j in 1:n)	{
		BB[i,j] <- 0
		for (k in 1:n)	BB[i,j] <- BB[i,j]+B1[i,k]*B2[k,j]
		}
	}
return(BB)
}

construct_Q_matrix_unordered <- function(k)	{
Q <- matrix(1,k,k)
for (i in 1:k)	Q[i,i] <- 1-k
return(Q);
}

construct_Q_matrix_ordered <- function(k)	{
Q <- matrix(0,k,k)
for (i in 2:(k-1))	{
	Q[i,i] <- -2;
	Q[i,i-1] <- Q[i,i+1] <- 1;
	}
Q[1,1] <- Q[k,k] <- -1;
Q[1,2] <- Q[k,k-1] <- 1;
return(Q);
}

construct_Q_matrix_for_divergent_hierarchical_multistate <- function(all_poss_combos,char_dependencies,UNKNOWN=-11,INAP=-22) {
all_chars <- as.numeric(colnames(all_poss_combos));
if (nrow(all_poss_combos)>100)	{
	print("Getting basic distances among state combinations")
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=T);
	} else	{
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=F);
	}
Q[Q>1] <- 0;
colnames(Q) <- rownames(Q) <- rownames(all_poss_combos);
	# now, go through and adjust transitions to "new" morphospace to be 1/poss. landing spots.
keychar_nos <- unique(char_dependencies);
for (kn in length(keychar_nos):1)	{
	kc <- match(keychar_nos[kn],all_chars);
	istates <- sort(unique(all_poss_combos[,kc]),decreasing=T);
	for (ks in 1:nstates[kc])	{
		istate <- istates[ks];
		state_combos <- subset(all_poss_combos,all_poss_combos[,kc]==istate);
		inap_pairs <- unique(which(state_combos==INAP,arr.ind = T)[,2]);
		if (length(inap_pairs)==1)	{
			other_ind_states <- unique(all_poss_combos[which(all_poss_combos[,inap_pairs]==INAP,arr.ind = T),kc]);
			} else	{
			other_ind_states <- unique(all_poss_combos[unique(which(all_poss_combos[,inap_pairs]==INAP,arr.ind = T)[,1]),kc]);
			}
		this_state <- (1:a_p_c)[all_poss_combos[,kc] == istate];
		not_other_ind_states <- istates[!istates %in% other_ind_states];
		for (noi in 1:length(not_other_ind_states))	{
			one_step_away <- (1:a_p_c)[all_poss_combos[,kc] %in% not_other_ind_states[noi]];
			Q[this_state,one_step_away] <- 1/length(one_step_away);
			}
		}
	}
nn <- rowSums(Q);
Q <- Q/nn;
#for (qq in 1:nrow(Q))	Q[qq,qq] <- stasis[qq];
for (qq in 1:nrow(Q))	Q[qq,qq] <- -1;
return(Q);
}

# char_dependencies <- secondary_dependencies
construct_Q_matrix_for_simple_hierarchical_multistate <- function(all_poss_combos,char_dependencies=c(),UNKNOWN=-11,INAP=-22)	{
orig_chars <- as.numeric(colnames(all_poss_combos));
ind_char_orig <- min(orig_chars[orig_chars %in% char_dependencies]);
ind_char <- match(ind_char_orig,orig_chars);
dep_dependencies <- char_dependencies[!all_chars==char_dependencies];
semi_indies <- unique(char_dependencies[char_dependencies!=ind_char_orig]);
nchars <- ncol(all_poss_combos);
dchars <- length(char_dependencies)-1;
nstates <- count_states(all_poss_combos);
secondaries <- unique(which(all_poss_combos==INAP,arr.ind = T)[,2]);
nstates[secondaries] <- nstates[secondaries]+1;
	#combos <- combos[order(combos[,1],combos[,2],combos[,3]),];
ind_states <- sort(unique(all_poss_combos[,ind_char][all_poss_combos[,ind_char]>=0]))
keystates <- accersi_key_states_for_independent_character_in_hierarchical(all_poss_combos,char_dependencies,INAP,UNKNOWN);
wrong_states <- ind_states[!ind_states %in% keystates[dep_dependencies %in% ind_char_orig,]];
if (nrow(all_poss_combos)>100)	{
	print("Getting basic distances among state combinations")
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=T);
	} else	{
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=F);
	}
	# if there are secondaries, the figure out how to weight them here!!!
colnames(Q) <- rownames(Q) <- rownames(all_poss_combos);
ttl_states <- nrow(Q);
rcombos <- ttl_states-length(wrstates);
Q[match(wrstates,all_poss_combos[,ind_char]),(length(wrstates)+1):ttl_states] <- 1/rcombos;
Q[Q>1] <- 0;
si <- 0;
while (si < length(semi_indies))	{
	si <- si+1;
	sc <- match(semi_indies[si],c(ind_char,dep_chars));
	cs <- (1:length(c(ind_char,dep_chars)))[secondary_dependencies==semi_indies[si]];
	relv_combos <- unique(all_poss_combos[,c(sc,cs)]);
	relv_combos <- relv_combos[(1:nrow(relv_combos))[!(1:nrow(relv_combos)) %in% unique(which(relv_combos==INAP,arr.ind=T)[,1])],];
	relv_combos_all <- all_poss_combos[,c(sc,cs)];
	if (length(cs)>1)	{
		key_combos <- which(all_poss_combos[,cs]==INAP,arr.ind=T)[,1]
		} else	{
		key_combos <- as.numeric(which(all_poss_combos[,cs]==INAP,arr.ind=T))
		}
	# fix revl_combos to get all of the right matches & not just those for one possibility
	key_combos <- key_combos[!key_combos %in% key_combos[all_poss_combos[,sc]==INAP]];
	relv_combos_2 <- all_poss_combos[,c(sc,cs)];
	relv_combos_2 <- relv_combos_2[!(1:ttl_states) %in% which(relv_combos_2==INAP,arr.ind=T)[,1],];
	combos_key <- match(rownames(relv_combos_2),rownames(Q));
	#	Q[key_combos,combos_key]==1;
	for (i in 1:length(key_combos))
		Q[key_combos[i],combos_key][Q[key_combos[i],combos_key]==1] <- 1/nrow(relv_combos);
	#	1/nrow(relv_combos)
	}
stasis <- rowSums(Q);
for (qq in 1:nrow(Q))	Q[qq,qq] <- -stasis[qq];


}

construct_Q_matrix_for_hierarchical_multistate <- function(all_poss_combos,char_dependencies,UNKNOWN=-11,INAP=-22) {
orig_char_nos <- as.numeric(colnames(all_poss_combos));
nchars <- ncol(all_poss_combos);
ncombos <- nrow(all_poss_combos);
nstates <- count_states(all_poss_combos);
char_dependencies_2 <- match(char_dependencies,orig_char_nos);
names(char_dependencies) <- names(char_dependencies_2) <- orig_char_nos;
if (nrow(all_poss_combos)>100)	{
	print("Getting basic distances among state combinations")
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=T);
	} else	{
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=F);
	}
Q[Q>1] <- 0;
colnames(Q) <- rownames(Q) <- rownames(all_poss_combos);

keystates <- accersi_key_states_for_independent_character_in_hierarchical(all_poss_combos,char_dependencies,INAP,UNKNOWN);
# find transitions that can jump to inapplicables
implosions <- vector(length=(nrow(Q)));
for (nn in 1:nrow(Q))	implosions[nn] <- sum(all_poss_combos[nn,]==INAP);

obs_combos <- all_poss_combos[!(1:nrow(all_poss_combos)) %in% unique(which(all_poss_combos==INAP,arr.ind = T)[,1]),];
nstates2 <- count_states(obs_combos);
if (sum(nstates2!=nstates)>0)	{
	for (nc in 1:nchars){
		if (nstates[nc]!=nstates2[nc])	{
			dummy <- obs_combos[1,];
			add_states <- unique(all_poss_combos[,nc][!all_poss_combos[,nc] %in% obs_combos[,nc]]);
			add_states <- add_states[add_states!=INAP];
			for (as in 1:length(add_states))	{
				dummy[nc] <- add_states[as];
				obs_combos <- rbind(obs_combos,dummy);
				}
			}
		}
	}
all_poss_no_inaps <- accersi_all_theoretical_character_state_combinations(obs_combinations=obs_combos);
all_poss_dists <- pairwise_differences_discrete(all_poss_no_inaps,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=F);
rownames(all_poss_dists) <- colnames(all_poss_dists) <- rownames(all_poss_no_inaps);

implosions <- sum_identicals <- vector(length=ncombos);
for (nc in 1:ncombos)	{
	dd <- accersi_distance_for_one_taxon_from_every_other(single_taxon=all_poss_combos[nc,],all_others=all_poss_no_inaps);
	sum_identicals[nc] <- sum(dd==0);
	implosions[nc] <- sum(all_poss_combos[nc,]==INAP);
	}
for (nc in 1:ncombos)	{
	if (implosions[nc]>0)	{
		voids <- (1:nchars)[all_poss_combos[nc,]==INAP];
		ind_chars <- unique(match(char_dependencies[voids],orig_char_nos));
		ind_chars_orig <- as.numeric(colnames(all_poss_combos))[ind_chars];
		if (length(ind_chars)==1)	{
			jumps <- (1:ncombos)[Q[nc,]==1 & all_poss_combos[,voids[1]]!=INAP];
			Q[nc,jumps] <- 1/length(jumps);
			} else	{
#			relv_combos <- all_poss_combos[Q[nc,]==1,voids];
#			cbind(Q[nc,],sum_identicals,all_poss_combos[,voids]);
			Q[nc,Q[nc,]==1] <- sum_identicals[Q[nc,]==1]/sum(sum_identicals[Q[nc,]==1]);
#			jumps <- ic_combos <- vector(length=length(ind_chars));
#			for (ic in 1:length(ind_chars))	{
#				ic_combos[ic] <- prod(nstates[voids][char_dependencies_2[voids]==ind_chars[ic]]);
#				if (ic==1)	{
#					jumps[ic] <- 1/ic_combos[ic];
#					} else	{
#					jumps[ic] <- jumps[ic-1]/ic_combos[ic];	# make sure that this is always ic-1
#					}
#				}
			
			}
		}
	}

#nn <- rowSums(Q);
#Q <- Q/nn;
stasis <- -rowSums(Q);
#rowsums(stasis)
for (nc in 1:ncombos)	Q[nc,nc] <- stasis[nc];
#for (nc in 1:ncombos)	Q[nc,nc] <- -1;
return(Q);
}

construct_Q_matrix_random <- function(k)	{
Q <- matrix(0,k,k);
for (i in 1:k)	{
	for (j in 1:k)	{
		if (i!=j)	{
			Q[i,j] <- runif(1);
			}
		}
	Q[i,] <- Q[i,]/sum(Q[i,]);
	}
for (i in 1:k)	Q[i,i] <- -1;
return(Q)
}

construct_Q_matrix_beta <- function(k)	{
Q <- matrix(0,k,k);
for (i in 1:k)	{
	for (j in 1:k)	{
		if (i!=j)	{
			Q[i,j] <- runif(1);
			}
		}
	Q[i,] <- Q[i,]/sum(Q[i,]);
	}
for (i in 1:k)	Q[i,i] <- -1;
return(Q)
}

construct_Q_matrix_F81 <- function(state_rates) {
k <- length(state_rates);
Q <- matrix(0,k,k);
for (i in 1:k)	{
	jj <- (1:k)[(1:k)!=i];
	Q[jj,i] <- state_rates[i];
	}
for (i in 1:k)	Q[i,i] <--sum(Q[i,]);
return(Q);
}

random_Q_matrix_F81 <- function(k, shape1=1, shape2=1) {
pi_st <- rbeta(n=k,shape1=shape1,shape2=shape2);
Q <- matrix(0,k,k);
for (i in 1:k)	{
	jj <- (1:k)[(1:k)!=i];
	Q[jj,i] <- pi_st[i];
	}
for (i in 1:k)	Q[i,i] <--sum(Q[i,]);
return(Q);
}


construct_Q_matrix_F81_driven <- function(k, shape1=1, shape2=1) {
pi_st <- sort(rbeta(n=k,shape1=shape1,shape2=shape2));
Q <- matrix(0,k,k);
for (i in 1:k)	{
	jj <- (1:k)[(1:k)!=i];
	Q[jj,i] <- pi_st[i];
	}
for (i in 1:k)	Q[i,i] <--sum(Q[i,]);
return(Q);
}


matrix_exponentiation <- function(Q)	{
n <- nrow(Q)
transition <- identity_matrix(n)
i <- 1
last_tran <- round(transition[1,1],8)
u <- 0
while (u==0)	{
	qq <- Q
	j <- 2
	while (j<=i)	{
#		qq <- qq %*% Q
		qq <- multiply_square_matrices(qq,Q)
		j <- j+1
		}
	transition <- transition+(qq/factorial(i))
	qq
	transition
	if (last_tran==round(transition[1,1],8))	{
		u <- 1;
		} else	{
		last_tran <- round(transition[1,1],8);
		}
	i <- i+1
	}
return(transition)
}

Mk_stasis_under_punctuation_ts <- function(bd,epsilon,k,lambda,mu)	{
# epsilon: per-branch probability of deriving a state
# k: number of n_states
# lambda: branching rate
# mu: extinction rate
# bd: duration of stasis
p_same_species <- exp(-mu*bd)/prob_paraclade_survival(a=1,p=lambda,q=mu,bd)
p_diff_species <- 1-p_same_species
p_outcomes <- Mk_punctuated(bd,k,epsilon,lambda,m=1)
return(p_same_species +(p_diff_species*p_outcomes[1]))
}

Mk_stasis_under_continuous_bd <- function(bd,alpha,k)	{
# alpha: instantaneous rate of state derivation
# k: number of n_states
# ts_ interval of time
return(Mk_continuous(bd,k,alpha)[1])
}

prob_paraclade_ancestor_extant <- function(lambda,mu,bd)	{
pexx <- prob_paraclade_survival(a=1,p=lambda,q=mu,tm=bd)
pex <- exp(-mu*bd)
return(pex/pexx)
}

pairwise_disimilarity_discrete <- function(chmatrix,n_states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
#	uchmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#uchmatrix <- unique(chmatrix)
uchmatrix <- chmatrix
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
					## different with single scored n_states
					if (uchmatrix[sp1,ch]>=0 && uchmatrix[sp2,ch]>=0)	{
						if (types[ch]==0 || n_states[ch]==2)	{
							diss <- diss+1
							} else {
							if (weight_ordered==FALSE)	{
								x <- abs(uchmatrix[sp1,ch]-uchmatrix[sp2,ch])
								} else if (weight_ordered==TRUE)	{
								x <- abs(uchmatrix[sp1,ch]-uchmatrix[sp2,ch])/(n_states[ch]-1)
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
							x <- sum(cvec2 %in% cvec1)/length(cvec2);
							diss <- diss+x;
#							x <- match(cvec1,cvec2)
#							if (length(cvec1)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec1)>1)	{
#								x <- mundus_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1);
							diss <- diss+x;
#							x <- match(cvec2,cvec1)		# pick up here!
#							if (length(cvec2)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec2)>1)	{
#								x <- mundus_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}			
						} # end case of 1 or 2 polymorphsics
					} # end case of disagreement
				} # end case of two coded characters
			}	# end going through characters
		if (num==0)	{
			diss <- 1
			num <- 2
			}
		udis_matrix[sp2,sp1] <- udis_matrix[sp1,sp2] <- diss/num
		}	# end comparison between sp2 & sp1
	}	#end going throuch species

#notu <- nrow(chmatrix)
#dis_matrix <- matrix(0,notu,notu)
#xxx <- prodlim::row.match(chmatrix,list(uchmatrix)[[1]])
#for (u in 1:(unotu-1))	{
#	yyy <- (1:notu)[xxx %in% u]
#	for (u2 in 2:unotu)	{
#		zzz <- (1:notu)[xxx %in% u2]
#		dis_matrix[zzz,yyy] <- dis_matrix[yyy,zzz] <- udis_matrix[u,u2]
#		}
#	}
return (udis_matrix)
}

pairwise_similarity_discrete <- function(chmatrix,n_states,types,weight_ordered=TRUE,polymorphs=TRUE,UNKNOWN=-11,INAP=-22)	{
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
					## different with single scored n_states
					if (chmatrix[sp1,ch]>=0 && chmatrix[sp2,ch]>=0)	{
						if (types[ch]==0 || n_states[ch]==2)	{
							diss <- diss+1
							} else {
							if (weight_ordered==FALSE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])
								} else if (weight_ordered==TRUE)	{
								x <- abs(chmatrix[sp1,ch]-chmatrix[sp2,ch])/(n_states[ch]-1)
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
#								x <- mundus_na_from_vector(x,-1)
#								for (q in 1:length(x))	if (x[q]==-1)	diss <- diss+1/(length(x))
#								}
							}	else {
							# numerator is # matches; denominator is # poss. matches
							x <- sum(cvec1 %in% cvec2)/length(cvec1)
#							x <- match(cvec2,cvec1)		# pick up here!
#							if (length(cvec2)==1 && is.na(x))	{
#								diss <- diss+1
#								} else if (length(cvec2)>1)	{
#								x <- mundus_na_from_vector(x,-1)
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

old_outgroup_search <- function()	{
	while (outgroup==-1 && nexusfile[tx_pt,1]=="\t" && tolower(nexusfile[tx_pt,2])=="t" && tolower(nexusfile[tx_pt,3])=="a" && tolower(nexusfile[tx_pt,4])=="x" && tolower(nexusfile[tx_pt,5])=="s" && tolower(nexusfile[tx_pt,6])=="e")	{
		j <- strsplit(nexus[tx_pt],split=" ",fixed=TRUE)[[1]]
		if (!is.na(match("outgroup",tolower(j))))	{
			out <- j[length(j)];
			out <- strsplit(out,split="",fixed=TRUE)[[1]]
			while (out[length(out)]==";" || out[length(out)]==" " || out[length(out)]==",")
				out <- out[1:((length(out)-1))]
#			otg <- match("outgroup",tolower(j))+1
#			while (j[otg]==":" || j[otg]==" " || j[otg]=="")	otg <- otg+1
			outg <- 0;
			outgroup <- c();
			through <- FALSE;
#			while (out[a]!="," && out[a]!=";")	{
			for (a in 1:length(out))	{
				if (out[a]=="-")	{
					# add taxon to list of outgroups
					outgroup <- c(outgroup,outg);
					through <- TRUE;
					outg <- 0
					} else if (out[a]==" ")	{
					if (through)	{
						outgroup <- c(outgroup,seq(1+outgroup[length(outgroup)],outg,by=1));
						} else	{
						outgroup <- c(outgroup,outg);
						}
					outg <- 0
					through <- FALSE;
					} else	{
					outg <- (10*outg)+as.numeric(out[a]);
					}
#				if (out[a]!="=")
#					outgroup <- (10*outgroup)+as.numeric(out[a])
#				a <- a+1
#				outgroup
				}
			if (through)	{
				outgroup <- c(outgroup,seq(1+outgroup[length(outgroup)],outg,by=1));
				} else	{
				outgroup <- c(outgroup,outg);
				}
			}
		tx_pt <- tx_pt+1
		}
	}

#### Routines for Mark & Curtis Project ####
prob_divergence_that_is_sampled <- function(lambdas,phis,min_exp=10^-10)	{
# lambdas: expected per-interval branchings (with expectation per cell in vector)
# phis: probability that a species diverging at this point in time every has a sampled successor
# min_exp: lowest probability of richnes that might be observed to use as a cutoff.
p_sampled_paraclade <- c()	# probability of 1…K paraclades evolving x 1-P[missing all K | phis]
for (b in 1:length(lambdas))	{
	latest_p <- 1;
	K <- 1;
	p_sampled_paraclade <- c(p_sampled_paraclade,0)
	while (latest_p>=min_exp)	{
		latest_p <- dpois(K,lambdas[b])
		p_sampled_paraclade[b] <- p_sampled_paraclade[b]+(1-(1-phis[b])^K)*latest_p
		K <- K+1
		}
	}
return(p_sampled_paraclade)
}

prob_predecessor_is_missed_over_intervals_continuous <- function(psis)	{
# psis: expected finds per interval, with each cell representing an interval
return(exp(-sum(psis)))
}


prob_no_sampled_sisters_over_intervals <- function(lambdas,phis)	{
# lambdas: expected per-interval branchings (with expectation per cell in vector)
# phis: probability that a species diverging at this point in time every has a sampled successor
return(exp(sum(log(1-prob_divergence_that_is_sampled(lambdas,phis)))))
}

prob_gap_in_phylogeny <- function(lambdas,phis,psis)	{
# lambdas: expected per-interval branchings (with expectation per cell in vector)
# phis: probability that a species diverging at this point in time every has a sampled successor
# psis: expected finds per interval, with each cell representing an interval
return(prob_no_sampled_sisters_over_intervals(lambdas,phis)*prob_predecessor_is_missed_over_intervals_continuous(psis))
}


# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
likelihood_of_alpha_given_divergences <- function(m_alpha,divergence_times,vector_tree,marginals,n_states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
if (m_alpha < 0)
  m_alpha <- abs(m_alpha);
#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha <- m_alpha*rate_quants;
log_like <- sapply(alpha,calculate_rate_likelihood_given_tree_divergences_and_characters,vector_tree,divergence_times,marginals,n_states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
likelihood_of_divergences_given_alpha <- function(divergence_times,m_alpha,vector_tree,marginals,n_states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
if (m_alpha < 0)	m_alpha <- abs(m_alpha);
#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha <- m_alpha*rate_quants;
log_like <- sapply(alpha,calculate_rate_likelihood_given_tree_divergences_and_characters,vector_tree,divergence_times,marginals,n_states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
probability_of_divergences_given_alpha_diversification_and_sampling <- function(divergence_times,m_alpha,vector_tree,marginals,n_states,BDS,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
if (m_alpha < 0)	m_alpha <- abs(m_alpha);
#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha <- m_alpha*rate_quants;
log_prob <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
log_like <- sapply(alpha,calculate_rate_likelihood_given_tree_divergences_and_characters,vector_tree,divergence_times,marginals,n_states);
#log_prob <- ;
return(BayesFactor::logMeanExpLogs(log_prob)+BayesFactor::logMeanExpLogs(log_like))
}

probability_of_divergences_given_epsilon_diversification_and_sampling <- function(divergence_times,m_epsilon,vector_tree,marginals,strat_ranges,n_states,BDS,rate_quants=1)	{
# m_epsilon: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
if (m_epsilon < 0)	m_epsilon <- abs(m_epsilon);
#alpha_sts <- alpha*1/((1:mxstates)-1);
epsilon <- m_epsilon*rate_quants;
log_prob <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vector_tree,BDS);
log_like <- sapply(epsilon,calculate_rate_likelihood_given_tree_divergences_and_characters,vector_tree,divergence_times,marginals,n_states);
#log_prob <- ;
return(BayesFactor::logMeanExpLogs(log_prob)+BayesFactor::logMeanExpLogs(log_like))
}

# routine to get likelihood of one rate given tree & divergence times.  Good for optimizing one rate or calculating rates over lognormal/gamma distributions
calculate_rate_likelihood_given_tree_divergences_and_characters <- function(alpha,vector_tree,divergence_times,marginals,n_states)	{
# divergence times: vector giving divergence dates for notus & htus
# marginals: likelihoods of each state at each branch/tip; = 1 for tips if state present (1/n if polymorphic with n different n_states)
# n_states: n_states for each character;
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
n_Nodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(n_Nodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- alpha*1/((1:mxstates)-1);
for (nn in n_Nodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha=alpha_sts[k]);	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
#		for (rc in 1:length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

calculate_rate_likelihood_given_tree_divergences_characters_punctuated_change_and_origination_rates <- function(epsilon,vector_tree,divergence_times,strat_ranges,marginals,n_states,origination_rates)	{
# divergence times: vector giving divergence dates for notus & htus
# marginals: likelihoods of each state at each branch/tip; = 1 for tips if state present (1/n if polymorphic with n different n_states)
# n_states: n_states for each character;
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
n_Nodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(n_Nodes,tnchars,mxstates));
gap_lgprobs <- 0;
epsilon_sts <- epsilon*1/((1:mxstates)-1);
nbins <- nrow(origination_rates);
bin_durations <- c(origination_rates$age[2:nbins]-origination_rates$age[1:(nbins-1)],0);
for (nn in n_Nodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
#	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	ave_lambda <- bd <- vector(length=tomy);
	min_branchings <- rep(1,tomy);
	d_dates <- unique(divergences[daughters]);
	poss_ancestor <- daughters[daughters<=notu];
	assumed_ancestor <- 0;
	if (length(d_dates)>1 && length(poss_ancestor)>0) {
	  poss_ancestor <- poss_ancestor[match(min(divergences[poss_ancestor]),divergences[poss_ancestor])];
	  if (divergences[poss_ancestor]==min(divergences[daughters]))  {
	    assumed_ancestor <- poss_ancestor;
	    min_branchings[match(assumed_ancestor,daughters)] <- 0;
	    }
	  }
	
	for (dt in 1:tomy)  {
	  if (daughters[dt]>notu) {
	    bna <- max(sum(divergence_times[htu]>origination_rates$age),1);
	    bnz <- min(sum(divergence_times[daughters[dt]]>origination_rates$age),nrow(origination_rates)-1);
	    strt <-divergence_times[htu];
	    end <-divergence_times[daughters[dt]];
	    } else  {
	    bna <- max(sum(divergence_times[daughters[dt]]>origination_rates$age),1);
	    bnz <- min(sum(observed_ranges$FAD[daughters[dt]]>origination_rates$age),nrow(origination_rates));
	    strt <- divergence_times[daughters[dt]];
	    end <-strat_ranges$FAD[daughters[dt]];
	    }
	  if (end==strt) {
	    ave_lambda[dt] <- 0;
	    bd[dt] <- 0;
	    } else  {
  	  top_offset <- abs(strt-origination_rates$age[bna]);
  	  bottom_offset <- abs(end-origination_rates$age[bnz]);
  	  relv_bin_lengths <- bin_durations[bna:bnz];
	    bd[dt] <- end-strt;
  	  if (bna!=bnz) {
  	    relv_bin_lengths[1] <- origination_rates$age[bna+1]-strt;
  	    relv_bin_lengths[length(relv_bin_lengths)] <- end-origination_rates$age[bnz];
  	    } else  {
  	    relv_bin_lengths <- end-strt;
  	    }
  	  ave_lambda[dt] <- sum(relv_bin_lengths*origination_rates$lambda[bna:bnz])/bd[dt];
	    }
	  }
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
#		mk <- sapply(bd,Mk_continuous,k,alpha=epsilon_sts[k]);	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
#		for (rc in 1:length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
    		if (daughters[d]!=assumed_ancestor){
			    mk <- Mk_punctuated(bd[d],k,epsilon=epsilon_sts[k],ave_lambda[d],m=min_branchings[d]);	# mk model for each node, with first value being P[net stasis]
				  dd <- daughters[d]	# dth taxon attached to node
				  for (st in 1:k)	{
					# prob net stasis
				  	mg_st <- mk[1]*marginals[dd,ch,st]
					# prob net change
				  	mg_ch <- mk[2]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
				  	marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
				  	}	
    		  }
			  }
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

  calculate_rate_likelihood_given_tree_divergences_characters_and_different_branch_rates <- function(alpha,vector_tree,divergence_times,marginals,n_states,branch_rates)	{
# divergence times: vector giving divergence dates for notus & htus
# marginals: likelihoods of each state at each branch/tip; = 1 for tips if state present (1/n if polymorphic with n different n_states)
# n_states: n_states for each character;
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- alpha*1/((1:mxstates)-1);
for (nn in nNodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha=alpha_sts[k]);	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
#		for (rc in 1:length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

prior_probability_of_evolutionary_history_given_basal_divergence <- function(divergence_times,vector_tree,BDS)	{
# modified 2019-08-26
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
#
notu <- match(-1,vector_tree)-1;
nNodes <- length(vector_tree)-notu;
gap_lgprobs <- 0;
ttu <- length(vector_tree);
for (nn in nNodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vector_tree==htu];
	for (d in 1:length(daughters))	{
		stgap <- (1:nrow(BDS))[BDS$ma > divergence_times[htu]];		# all intervals after divergence
		stgap <- stgap[BDS$ma[stgap] < divergence_times[daughters[d]]];	# intervals between divergence & appearance
		gap_lgprobs <- gap_lgprobs+sum(log(BDS$pgap[stgap]));
		}
	}
return(gap_lgprobs);
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
# created 2020-07-11
likelihood_of_rate_shift_given_divergences <- function(m_alphas,divergence_times,rate_shift,vector_tree,marginals,n_states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
m_alphas <- abs(m_alphas);
m_alpha_1 <- m_alphas[1];
m_alpha_2 <- m_alphas[2];

#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha_1 <- m_alpha_1*rate_quants;
alpha_2 <- m_alpha_2*rate_quants;
alphas <- cbind(m_alpha_1*rate_quants,m_alpha_2*rate_quants);
#log_like_e <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_1,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
#log_like_u <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_2,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);

log_like <- apply(alphas,MARGIN=1,FUN=calculate_rate_shift_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,rate_shift=rate_shift,marginals=marginals,n_states=n_states);
#lll <- vector(length=length(rate_quants));
#for (ii in 1:length(rate_quants))	{
#	lll[ii] <- calculate_rate_shift_likelihood_given_tree_divergences_and_characters(alphas=alphas[ii,],vector_tree,divergence_times,rate_shift,marginals,n_states);
#	}

#log_like <- c();
#for (rq in 1:nrow(alphas))
#	log_like <- c(log_like,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alphas[rq,],vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states));
#	log_like <- sapply(alphas,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
# modified 2020-07-11
likelihood_of_early_and_late_alphas_given_divergences <- function(m_alphas,divergence_times,branch_rates,vector_tree,marginals,n_states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
m_alphas <- abs(m_alphas);
m_alpha_1 <- m_alphas[1];
m_alpha_2 <- m_alphas[2];

#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha_1 <- m_alpha_1*rate_quants;
alpha_2 <- m_alpha_2*rate_quants;
alphas <- cbind(m_alpha_1*rate_quants,m_alpha_2*rate_quants);
#log_like_e <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_1,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
#log_like_u <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_2,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);

log_like <- apply(alphas,MARGIN=1,FUN=calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
#log_like <- c();
#for (rq in 1:nrow(alphas))
#	log_like <- c(log_like,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alphas[rq,],vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states));
#	log_like <- sapply(alphas,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# created 2020-07-11
likelihood_of_early_and_late_alphas_given_divergences_shift <- function(m_alphas,divergence_times,boaty_mcboatface="",vector_tree,marginals,n_states,shift,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vector_tree: vector giving htu number of node from which each htu & otu evolved
# n_states: number fo state per character
#
m_alphas <- abs(m_alphas);
m_alpha_1 <- m_alphas[1];
m_alpha_2 <- m_alphas[2];

#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha_1 <- m_alpha_1*rate_quants;
alpha_2 <- m_alpha_2*rate_quants;
alphas <- cbind(m_alpha_1*rate_quants,m_alpha_2*rate_quants);
#log_like <- apply(alphas,MARGIN=1,FUN=calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=boaty_mcboatface,marginals=marginals,n_states=n_states);
log_like <- apply(alphas,MARGIN=1,FUN=calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters_shift,vector_tree=vector_tree,divergence_times=divergence_times,shift=shift,marginals=marginals,n_states=n_states);
#log_like <- c();
#for (rq in 1:nrow(alphas))
#	log_like <- c(log_like,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alphas[rq,],vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states));
#	log_like <- sapply(alphas,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vector_tree=vector_tree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,n_states=n_states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# created 2020-07-25
calculate_rate_shift_likelihood_given_tree_divergences_and_characters <- function(alphas,vector_tree,divergence_times,rate_shift,marginals,n_states,rate_quants=1)	{
# alphas: the first and second rate
# NOTE: Do this separately for all 4 quartiles if we have rate variation!!!
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- c();
for (aa in 1:length(alphas))
	alpha_sts <- rbind(alpha_sts,alphas[aa]*1/((1:mxstates)-1));

for (nn in nNodes:1)	{
	htu <- nn+notu;
#	b_rate <- branch_rates[htu];
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	node_onset <- divergence_times[htu];
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
		mk <- array(0,dim=c(2,length(daughters)))
		for (dd in 1:length(daughters))	{
			f1 <- daughters[dd];
			f1_onset <- divergence_times[f1];
			if (node_onset < rate_shift && f1_onset > rate_shift)	{
				r1 <- abs(rate_shift-node_onset)/(f1_onset-node_onset);
				r2 <- 1-r1;
				} else if (f1_onset <= rate_shift)	{
				r1 <- 1;
				r2 <- 0;
				} else	{
				r1 <- 0;
				r2 <- 1;
				}
			br_alphas <- (r1*alpha_sts[1,k]) + (r2*alpha_sts[2,k]);
#			mk[,dd] <- 0
			mk[,dd] <- Mk_continuous(bd=bd[dd],k,br_alphas);
#				mk[,dd] <- (mk[,dd] + Mk_continuous(bd=bd[dd],k,br_alphas[bra])/length(br_alphas));
			}
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc + 1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k]);
#			print(logmarginals[nn,ch,1:k])
			}
		}
#	print(sum(logmarginals[nn,,]));
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like);
}

calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters <- function(alphas,vector_tree,divergence_times,branch_rates,marginals,n_states,rate_quants=1)	{
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- c();
for (aa in 1:length(alphas))
	alpha_sts <- rbind(alpha_sts,alphas[aa]*1/((1:mxstates)-1));

for (nn in nNodes:1)	{
	htu <- nn+notu;
	b_rate <- branch_rates[htu];
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha=alpha_sts[b_rate,k]);	# mk model for each node, with first value being P[net stasis]
#		print(dim(mk));
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc + 1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k]);
#			print(logmarginals[nn,ch,1:k])
			}
		}
	print(sum(logmarginals[nn,,]))
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

# redone to use the shift date instead.
# created 2020-07-11
calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters_shift <- function(alphas,vector_tree,divergence_times,shift,marginals,n_states,rate_quants=1)	{
notu <- match(-1,vector_tree)-1;
tree_base <- notu+1;
ttu <- length(vector_tree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- c();
state_options <- 1/1:(mxstates-1);
for (aa in 1:length(alphas))
	alpha_sts <- rbind(alpha_sts,as.numeric(alphas[aa])/(state_options));
#	alpha_sts <- rbind(alpha_sts,alphas[aa]*(1/((2:mxstates)-1)));
colnames(alpha_sts) <- paste("k",2:mxstates,sep="")
rownames(alpha_sts) <- c(paste("R1_",1:4,sep=""),paste("R2_",1:4,sep=""))

for (nn in nNodes:1)	{
	htu <- nn+notu;
#	b_rate <- branch_rates[htu];
	daughters <- (1:ttu)[vector_tree==htu];
	tomy <- length(daughters);
	br_dr <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (ch in 1:tnchars)
		marginals[htu,ch,1:n_states[ch]] <- 1;	

	for (d in 1:tomy)	{
		# set up rate quantiles, with branches spanning shifts getting weighted rates
		dd <- daughters[d];
		if (divergence_times[htu]<shift)	{
			if (divergence_times[dd]<=shift)	{
				rate_quants <- alpha_sts[1:4,];
				} else	{
				p_bng <- abs(shift-divergence_times[htu])/abs(divergence_times[dd]-divergence_times[htu]);
				p_nbng <- 1-p_bng;
				rate_quants <- (p_bng*alpha_sts[1:4,])+((p_nbng)*alpha_sts[5:8,]);
				}
			} else	{
			rate_quants <- alpha_sts[5:8,];
			}
		for (k in 2:mxstates)	{
			rchars <- (1:tnchars)[n_states==k];
			mk <- c();
			for (qr in 1:4)
				mk <- rbind(mk,Mk_continuous(bd=br_dr[d],k=k,alpha=rate_quants[qr,k-1]))
			for (rc in 1:length(rchars))	{
				ch <- rchars[rc];
				for (st in 1:k)	{
					# prob net stasis
					#mg_st <- mk[1,d]*marginals[dd,ch,st]
					mg_st <- mean(mk[1:4,1]*marginals[dd,ch,st]);
					# prob net change
					mg_ch <- mean(mk[1:4,2]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]);
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch));
					}
				}
			}
		} # end examining daughter taxa;
	
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[n_states==k];
		mk <- sapply(br_dr,Mk_continuous,k,alpha=alpha_sts[b_rate,k-1]);	# mk model for each node, with first value being P[net stasis]
#		print(dim(mk));
		for (rc in 1:length(rchars))	{
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch));
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k]);
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

center_budding_phylogeny_sc <- function(vector_tree,durations,sampled_ancestors)	{
# vector_tree: tree as a vector with ancestral node HTU given for each taxon (-1 for basal node)
# durations: origins and ends
# function to get relative positions of lineages on a phylogeny;
venntree <- transform_vector_tree_to_venn_tree(vector_tree);
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
node_richness <- accersi_node_richness_from_vector_tree(vector_tree = vector_tree);
nNode <- nrow(mtree);
notu <- length(vector_tree)-nNode;
node_ages <- c();
if (nrow(durations)==(notu+nNode))	{
	branch_ages <- durations[,1];
	} else	{
	for (nd in 1:nNode)	node_ages <- c(node_ages,min(durations[venntree[nd,venntree[nd,]>0],1]));
	branch_ages <- c(durations[1:notu,1],node_ages);
	}
if (length(sampled_ancestors) < (notu+nNode))
	sampled_ancestors <- c(rep(0,notu),sampled_ancestors);

ttl_richness <- c(rep(1,notu),node_richness);
##patristic_distances <- accersi_patristic_distance_from_base(atree=mtree);#max_nodes <- max(patristic_distances);

last_left <- "left";	# move up the axis
last_right <- "right";	# move down the axis
accounted <- c();
nd <- 0;
phy_pos <- rep(0,nNode+notu);
#for (nd in 1:nNode)	{
while (nd < nNode)	{
	nd <- nd+1;
	htu <- nd+notu;			# htu number of node;
	tf1 <- sum(mtree[nd,]>0);
	if (sampled_ancestors[htu]!=0)	{
		tf1 <- tf1-1;
		phy_pos[sampled_ancestors[htu]] <- phy_pos[notu+nd];
		}  
	f1 <- mtree[nd,!mtree[nd,] %in% c(sampled_ancestors[htu],0)];
	f1 <- f1[order(-ttl_richness[f1])];
	if (length(f1)>2)	{
		right <- left <- 0;
		prop_richness <- ttl_richness[f1]/sum(ttl_richness[f1]);
		f1cc <- length(f1);
		left <- f1[1]; 
		sum_prop <- prop_richness[1];
		while (sum_prop <= 0.45)	{
			sum_prop <- sum_prop+prop_richness[f1cc];
			left <- c(left,f1[f1cc]);
			f1cc <- f1cc-1;
			}
		right <- f1[!f1 %in% left];
		right <- right[order(-abs(branch_ages[right]))];
		left <- left[order(-abs(branch_ages[left]))];
		# shift rest of the tree away from ancestral node
#		phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]]-sum(ttl_richness[right]);
#		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+sum(ttl_richness[left]);
		rr <- 1;
		phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]]-ttl_richness[right[rr]]
		phy_pos[right[rr]] <- phy_pos[htu]-ttl_richness[right[rr]];
		while (rr < length(right))	{
			rr <- rr+1;
			if (sum(phy_pos<phy_pos[htu])>0)
				phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]]-(2*ttl_richness[right[rr]]);
#			phy_pos[right[rr]] <- (phy_pos[right[rr-1]]-ttl_richness[right[rr-1]])-ttl_richness[right[rr]];
			phy_pos[right[rr]] <- phy_pos[htu]-ttl_richness[right[rr]];
#			phy_pos[c(htu,right)];
			}
		ll <- 1;
		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+ttl_richness[left[ll]]
#		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+sum(ttl_richness[left[ll]]);
		phy_pos[left[ll]] <- phy_pos[htu] + ttl_richness[left[ll]];
		while (ll < length(left))	{
			ll <- ll+1;
			if (sum(phy_pos>phy_pos[htu])>0)
				phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+(2*ttl_richness[left[ll]]);
#			phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]]+2*ttl_richness[left[ll]]);
#			phy_pos[left[ll]] <- (phy_pos[left[ll-1]]+ttl_richness[left[ll-1]])+ttl_richness[left[ll]];
			phy_pos[left[ll]] <- phy_pos[htu]+ttl_richness[left[ll]];
			}
		}	else if (length(f1)==2)	{
		f1 <- f1[order(-ttl_richness[f1])];
		if (phy_pos[htu]<phy_pos[1+notu])	{
			# going left is positive, so shift everything above this up by this amount
			if (last_right=="right")	{
				right <- f1[2];
				left <- f1[1];
				last_right <- "left"
				} else	{
				right <- f1[1];
				left <- f1[2];
				last_right <- "right"
				}
			} else	{
			# going right is negative, so shift everything below this down by this amount
			if (last_left=="left")	{
				right <- f1[1];
				left <- f1[2];
				last_left <- "right";
				} else	{
				right <- f1[2];
				left <- f1[1];
				last_left <- "left";
				}
			}
		# shift rest of the tree away from ancestral node
		phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[left];
		phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[right];
		phy_pos[left] <- phy_pos[htu] + ttl_richness[left];
		phy_pos[right] <- phy_pos[htu] - ttl_richness[right];
		}	else if (length(f1)==1)	{
		if (phy_pos[htu]<phy_pos[1+notu])	{
			if (last_right=="right")	{
				# going left is positive, so shift everything above this up by this amount
#				phy_y[phy_y>(phy_y[htu]+ttl_richness[f1])] <- phy_y[phy_y>phy_y[htu]] + ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] + ttl_richness[f1];
				last_right <- "left";
				} else {
				# going right is negative, so shift everything below this down by this amount
#				phy_y[phy_y<(phy_y[htu]-ttl_richness[f1])] <- phy_y[phy_y>ttl_richness[f1]] - ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] - ttl_richness[f1];
				last_right <- "right";
				}
			} else	{
			if (last_left=="right")	{
#				phy_y[phy_y>(phy_y[htu]+ttl_richness[f1])] <- phy_y[phy_y>phy_y[htu]] + ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos>phy_pos[htu]] <- phy_pos[phy_pos>phy_pos[htu]] + ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] + ttl_richness[f1];
				last_left <- "left";
				} else	{
				# going right is negative, so shift everything below this down by this amount
#				phy_y[phy_y<(phy_y[htu]-ttl_richness[f1])] <- phy_y[phy_y>ttl_richness[f1]] - ttl_richness[f1];
				# shift rest of the tree away from ancestral node
				phy_pos[phy_pos<phy_pos[htu]] <- phy_pos[phy_pos<phy_pos[htu]] - ttl_richness[f1];
				phy_pos[f1] <- phy_pos[htu] - ttl_richness[f1];
				last_left <- "right";
				}
			}
		} 
	
	phy_pos <- phy_pos-phy_pos[notu+1];	# recenter around the base of the tree
#	print(c(nd,phy_pos));
	# now do species
	accounted <- c(accounted,mtree[nd,mtree[nd,]>0]);
	}
final_pos <- match(phy_pos,sort(unique(phy_pos)));
needed_edits <- hist(final_pos,breaks=0:max(final_pos),plot=F)$counts
too_many_here <- (1:length(needed_edits))[needed_edits>2];
while (length(too_many_here)>0) {
	ttl_minions <- 1:(notu+nNode);
	tmh <- 0;
	while (tmh < length(too_many_here))	{
		tmh <- tmh+1;
		problems <- ttl_minions[final_pos==too_many_here[tmh]];	# taxa overlapping each other
		these_nodes <- c();
		for (pp in 1:length(problems))
			these_nodes <- c(these_nodes,which(mtree==problems[pp],arr.ind = T)[1]);	# get the nodes containing problem cases;
		problem_ancestors <- problems[(notu+these_nodes) %in% problems];				# separate out sampled ancestors
		problem_ancestors_htu <- notu+these_nodes[match(problem_ancestors,problems)];	# keep track of the htu to which they belong, however!
		these_nodes <- these_nodes[!problems %in% problem_ancestors];
		problems <- problems[!problems %in% problem_ancestors];	# remove sampled ancestors for now
		starting_points <- final_pos[notu+these_nodes];			# positions of ancestral nodes/taxa
		adjust2 <- adjust <- starting_points-too_many_here[tmh];
		adjust2[adjust<0]<- -(length(adjust[adjust<0]):1);
		adjust2[adjust>0]<- 1:length(adjust[adjust>0]);
		final_pos[final_pos<too_many_here[tmh]] <- final_pos[final_pos<too_many_here[tmh]]+min(adjust2);
		final_pos[final_pos>too_many_here[tmh]] <- final_pos[final_pos>too_many_here[tmh]]+max(adjust2);
		phy_pos[phy_pos<phy_pos[problems[1]]] <- phy_pos[phy_pos<phy_pos[problems[1]]]+min(adjust2);
		phy_pos[phy_pos>phy_pos[problems[1]]] <- phy_pos[phy_pos>phy_pos[problems[1]]]+max(adjust2);
		final_pos[problems] <- final_pos[problems]+adjust;
		phy_pos[problems] <- phy_pos[problems]+adjust2;
		final_pos[problem_ancestors] <- final_pos[problem_ancestors_htu];
		phy_pos[problem_ancestors] <- phy_pos[problem_ancestors_htu];
		too_many_here <- too_many_here+max(adjust2);
		}
	final_pos <- match(phy_pos,sort(unique(phy_pos)));
	needed_edits <- hist(final_pos,breaks=0:max(final_pos),plot=F)$counts
	too_many_here <- (1:length(needed_edits))[needed_edits>2]
	}
return(final_pos);
}

#(1,(2,4,(9,(7,(10,((18,(26,(33,36)),(34,40)),(19,27))),(13,25,(14,16)))),(3,15,(5,(11,17,24,(21,29),(22,31)),(28,(20,(30,32))))),(6,(8,39,(12,23,35,37,38)))))
ladderize_vector_tree <- function(vector_tree)	{
  # modified 2022-09-05
  lvector_tree <- lvector_tree_old <- vector_tree;
  venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
  nNodes <- nrow(venn_tree);
  notu <- max(venn_tree);
  ttu <- length(vector_tree);
  matrix_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
  node_rich <- tally_node_richness_from_vector_tree(vector_tree);
  node_rich_all <- c(rep(0,notu),node_rich);
  node_daughters <- vector(length=nNodes);
  #for (nd in 1:nNodes)  {
  nd <- 0;
  while (nd < nNodes)  {
    #sort(matrix_tree[nd,matrix_tree[nd,]>0])
    nd <- nd+1;
    node_daughters[nd] <-  sum(matrix_tree[nd,]>0);
    daughter_clades <- matrix_tree[nd,matrix_tree[nd,]>notu];
    
    if (length(daughter_clades)>0) {
      daughter_clades_ranked <- daughter_clades[order(node_rich_all[daughter_clades])];
      daughter_nodes <- daughter_clades-notu;
      daughter_nodes_ranked <- daughter_nodes[order(node_rich_all[daughter_clades])];
      lower_nodes <- (1:max(1,(nd-1)))[!(1:max(1,(nd-1))) %in% c(nd,daughter_nodes)];
      upper_nodes <- (nd:nNodes)[!(nd:nNodes) %in% c(nd,daughter_nodes)];
      venn_tree <- venn_tree[c(lower_nodes,nd,daughter_nodes_ranked,upper_nodes),];
      matrix_tree <- transform_venn_tree_to_matrix_tree(venn_tree);
      node_rich <- tally_node_richness_from_venn_tree(venn_tree);
      node_rich_all <- c(rep(0,notu),node_rich);
      lvector_tree <- transform_venn_tree_to_vector_tree(venn_tree);
      #    (notu+(1:nNodes))[!(notu+(1:nNodes)) %in% unique(lvector_tree[lvector_tree>0])]
    }
  }
  return(lvector_tree);
}

#(1,(2,4,(9,(7,(10,((18,(26,(33,36)),(34,40)),(19,27))),(13,25,(14,16)))),(3,15,(5,(11,17,24,(21,29),(22,31)),(28,(20,(30,32))))),(6,(8,39,(12,23,35,37,38)))))
ladderize_vector_tree_new <- function(vector_tree)	{
# modified 2022-09-05
lvector_tree <- lvector_tree_old <- vector_tree;
venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
nNodes <- nrow(venn_tree);
notu <- max(venn_tree);
ttu <- length(vector_tree);
matrix_tree_orig <- matrix_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
node_rich <- tally_node_richness_from_vector_tree(vector_tree);
node_rich_all <- c(rep(0,notu),node_rich);
node_daughters <- vector(length=nNodes);
#for (nd in 1:nNodes)  {
nd <- 0;
#venn_tree_new <- venn_tree[1,];
while (nd < nNodes) {
  nd <- nd+1;
  f1 <- matrix_tree[nd,matrix_tree[nd,]>0];
  dcl <- f1[f1>notu];
  if (length(dcl)>1)  {
    dcl2 <- dcl[order(node_rich_all[dcl])];
    ncl <- dcl-notu;
    ncl2 <- dcl2-notu;
    matrix_tree[ncl,] <- matrix_tree[ncl2,];
    venn_tree[ncl,] <- venn_tree[ncl2,];
    node_rich_all[dcl] <- node_rich_all[dcl2];
    node_rich[ncl] <- node_rich[ncl2];
    }
  }
ladder <- transform_matrix_tree_to_vector_tree(matrix_tree);
return(lvector_tree);
}

# base_alpha <- 0.05;
# burst <- 3;
# burst_start <- 66.0;
# burst_end <- 61.6;
# t1 <- 68
# t2 <- 62;
#t1 <- 40;
#t2 <- 50;
#t1 <- 64 ;
#t2 <- 60 ;
#proportion_overlap(t1,t2,b1,b2)
proportion_overlap <- function(t1,t2,b1,b2)	{
return((max(0,(b1-t2))+max(0,(b2-t1))-max(0,(b2-t2))+max(0,(b1-t1)))/(t1-t2));
}

delayed_burst_at_time_tt <- function(base_alpha,burst,burst_start,burst_end,t1,t2)	{
#if (t1<=burst_start && t2>=burst_end)	{
#	return(base_alpha*burst);
#	} else if (t1>burst_start && t2>=burst_end)	{
#	tt <- t1-t2;
#	return((((burst_start-t2)/tt)*burst*base_alpha)+(((t1-burst_start)/tt)*base_alpha))
#	} else if (t1>=burst_end && t2<=burst_end)	{
#	tt <- t1-t2;
#	return((((t1-burst_end)/tt)*burst*base_alpha)+(((t2-burst_end)/tt)*base_alpha))
#	}
return(burst*base_alpha*((max(0,(b1-t2))+max(0,(b2-t1))-max(0,(b2-t2))+max(0,(b1-t1)))/(t1-t2))+base_alpha*(1-((max(0,(b1-t2))+max(0,(b2-t1))-max(0,(b2-t2))+max(0,(b1-t1)))/(t1-t2))));
}

early_burst_rate_at_time_tt_old <- function(final_alpha,bang_mag,max_age,tt)	{
rel_bang <- 2^bang_mag;
final_alpha*(rel_bang^(tt/max_age))
#(rel_bang*(((rel_bang^(tt1/max_age))/max_age) - ((rel_bang^(tt2/max_age))/max_age))/log(rel_bang))
return(final_alpha*(rel_bang^(tt/max_age)));
}
# final_alpha <- 0.1;
# bang_mag <- 2;
# tt1 <- 6
# tt2 <- 5
# max_age <- 10;
early_burst_rate_at_time_tt <- function(init_alpha,bang_mag,max_age,tt)	{
rel_bang <- 2^bang_mag;
#final_alpha*(rel_bang^(tt/max_age));
init_alpha*rel_bang^((tt-max_age)/max_age);

#(rel_bang*(((rel_bang^(tt1/max_age))/max_age) - ((rel_bang^(tt2/max_age))/max_age))/log(rel_bang))
#return(final_alpha*(rel_bang^(tt/max_age)));
return(init_alpha*rel_bang^((tt-max_age)/max_age));
}

# t1 <- 11
# t2 <- 10
# final_alpha <- 0.05
# init_alpha <- 0.20
early_burst_rate_at_from_t1_to_t2 <- function(final_alpha,bang_mag,max_age,t1,t2)	{
rel_bang <- 2^bang_mag;
tt <- abs(t2-t1);	# branch duration

#((rel_bang^((t1-max_age)/max_age))-(rel_bang^((t2-t1)/t1)))/log(rel_bang)
#final_alpha*(((rel_bang^(t1/max_age)-rel_bang^(t2/max_age))/log(rel_bang))/(tt/max_age))
#answer <- init_alpha*(((rel_bang^((t1-max_age)/max_age)-rel_bang^((t2-max_age)/max_age))/log(rel_bang))/(tt/max_age))
answer <- init_alpha*(((rel_bang^((t1-max_age)/max_age)-rel_bang^((t2-max_age)/max_age))/log(rel_bang))*(max_age/tt));

#(rel_bang*(((rel_bang^(t1/max_age))/max_age) - ((rel_bang^(t2/max_age))/max_age))/log(rel_bang))
#final_alpha*(((rel_bang^(t1/t2)-rel_bang^(t2/max_age))/log(rel_bang))/(tt/max_age));
#init_alpha/(((rel_bang^(t2/t1)-rel_bang^((t2-max_age)/max_age))/log(rel_bang))/(tt/max_age))
#	(((rel_bang^(t1/t2)-rel_bang^(t2/max_age))/log(rel_bang))/(tt/max_age))
return(answer);
}

#(init_alpha*rel_bang^((t1-max_age)/max_age)*init_alpha*rel_bang^((t2-max_age)/max_age))^0.5
#final_alpha*
#init_alpha*(((rel_bang^(t1/max_age)-rel_bang^(t2/max_age))/log(rel_bang))/(tt/max_age))/rel_bang

#for (tt1 in 9:1)	{
#	tt2 <- tt1-2;
#	print(early_burst_rate_at_time_tt(final_alpha,bang_mag,max_age,tt1,tt2));
#	}

#### Return to the 1990s ####
calculate_SCI <- function(vector_tree, fad)	{
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
tree_dates <- date_taxa_on_tree_simple(vector_tree=vector_tree,FAs=fad);
notu <- length(fad)
clades <- nrow(mtree);

tomy <- daughter_clades <- vector(length=clades);
sci_numerator <- sci_denominator <- 0;
for (nd in clades:1)	{
	tomy[nd] <- sum(mtree[nd,]>0);
	daughter_clades[nd] <- sum(mtree[nd,]>notu);
	f1 <- mtree[nd,mtree[nd,]!=0];
	if (daughter_clades[nd]>0)	{
		daughters <- mtree[nd,mtree[nd,]>0];
		daughters_w_daughters <- mtree[nd,mtree[nd,]>notu];
		consistency <- !tree_dates[daughters_w_daughters] %in% min(tree_dates[daughters]);
		sci_numerator <- sci_numerator+sum(consistency);
		sci_denominator <- sci_denominator+length(consistency);
		}
	}
return(sci_numerator/sci_denominator);
}

calculate_RCI <- function(vector_tree,strat_ranges,apomorphies)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
atu <- length(vector_tree);
tree_dates <- date_taxa_on_tree_simple(vector_tree=vector_tree,FAs=strat_ranges$FAD);
ghost_lineages <- naive_ghost_lineages <- tree_dates[1:notu]-tree_dates[vector_tree[1:notu]]
ghost_taxa <- naive_ghost_taxa <- tree_dates[(notu+2):atu]-tree_dates[vector_tree[(notu+2):atu]]
ghosts <- naive_ghosts <- c(naive_ghost_lineages,0,naive_ghost_taxa);
MIG_nv <- sum(naive_ghost_lineages) + sum(naive_ghost_taxa)
SRL <- sum(strat_ranges$LAD-strat_ranges$FAD);

ancestral_nodes <- sort(vector_tree[apomorphies==0],decreasing=T);
an <- 0;
while (an < length(ancestral_nodes))	{
	an <- an+1;
	f0 <- (1:atu)[vector_tree==ancestral_nodes[an]];
	f1 <- f0[apomorphies[f0]>0];
	anc <- f0[!f0 %in% f1];
	if (strat_ranges$FAD[anc]==tree_dates[ancestral_nodes[an]])	{
		ghosts[f1] <- tree_dates[f1]-strat_ranges$LAD[anc];
		ghosts[f1][ghosts[f1]<0] <- 0;
		}
	}
MIG <- sum(ghosts);
RCI <- data.frame(cladistic=as.numeric(1-(MIG_nv/SRL)),phylogenetic=as.numeric(1-(MIG/SRL)),stringsAsFactors = F)
return(RCI)
}

calculate_GER <- function(vector_tree,strat_ranges,apomorphies)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
atu <- length(vector_tree);
tree_dates <- date_taxa_on_tree_simple(vector_tree=vector_tree,FAs=strat_ranges$FAD);
ghost_lineages <- naive_ghost_lineages <- tree_dates[1:notu]-tree_dates[vector_tree[1:notu]]
ghost_taxa <- naive_ghost_taxa <- tree_dates[(notu+2):atu]-tree_dates[vector_tree[(notu+2):atu]]
ghosts <- naive_ghosts <- c(naive_ghost_lineages,0,naive_ghost_taxa);
MIG_cl <- sum(naive_ghost_lineages) + sum(naive_ghost_taxa)

ancestral_nodes <- sort(vector_tree[apomorphies==0],decreasing=T);
an <- 0;
while (an < length(ancestral_nodes))	{
	an <- an+1;
	f0 <- (1:atu)[vector_tree==ancestral_nodes[an]];
	f1 <- f0[apomorphies[f0]>0];
	anc <- f0[!f0 %in% f1];
	if (strat_ranges$FAD[anc]==tree_dates[ancestral_nodes[an]])	{
		ghosts[f1] <- tree_dates[f1]-strat_ranges$LAD[anc];
		ghosts[f1][ghosts[f1]<0] <- 0;
		}
	}
MIG_ph <- sum(ghosts);
max_gaps <- maximum_gaps(strat_ranges);
min_gaps_cladistic <- cladistic_minimum_gaps(strat_ranges);
min_gaps_phylogenetic <- phylogenetic_minimum_gaps(strat_ranges);

GER <- data.frame(cladistic=as.numeric((MIG_cl-min_gaps_cladistic)/(max_gaps-MIG_cl)),phylogenetic=as.numeric((MIG_ph-min_gaps_phylogenetic)/(max_gaps-MIG_ph)),stringsAsFactors = F);
return(GER);
}

calculate_MSM <- function(vector_tree,strat_ranges,apomorphies)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
mtree <- transform_vector_tree_to_matrix_tree(vector_tree);
atu <- length(vector_tree);
tree_dates <- date_taxa_on_tree_simple(vector_tree=vector_tree,FAs=strat_ranges$FAD);
ghost_lineages <- naive_ghost_lineages <- tree_dates[1:notu]-tree_dates[vector_tree[1:notu]]
ghost_taxa <- naive_ghost_taxa <- tree_dates[(notu+2):atu]-tree_dates[vector_tree[(notu+2):atu]]
ghosts <- naive_ghosts <- c(naive_ghost_lineages,0,naive_ghost_taxa);
MIG_cl <- sum(naive_ghost_lineages) + sum(naive_ghost_taxa)

ancestral_nodes <- sort(vector_tree[apomorphies==0],decreasing=T);
an <- 0;
while (an < length(ancestral_nodes))	{
	an <- an+1;
	f0 <- (1:atu)[vector_tree==ancestral_nodes[an]];
	f1 <- f0[apomorphies[f0]>0];
	anc <- f0[!f0 %in% f1];
	if (strat_ranges$FAD[anc]==tree_dates[ancestral_nodes[an]])	{
		ghosts[f1] <- tree_dates[f1]-strat_ranges$LAD[anc];
		ghosts[f1][ghosts[f1]<0] <- 0;
		}
	}
MIG_ph <- sum(ghosts);
min_gaps_cladistic <- cladistic_minimum_gaps(strat_ranges);
min_gaps_phylogenetic <- phylogenetic_minimum_gaps(strat_ranges);

MSM <- data.frame(cladistic=as.numeric(min_gaps_cladistic/MIG_cl),phylogenetic=as.numeric(min_gaps_phylogenetic/MIG_ph),stringsAsFactors = F);
}

cladistic_minimum_gaps <- function(FAs)	{
return(max(FAs)-min(FAs));
}

phylogenetic_minimum_gaps <- function(strat_ranges)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
notu <- nrow(strat_ranges);
rownames(strat_ranges) <- 1:notu
gap_pre <- gap_pst <- vector(length=notu);
for (n in 1:notu)	{
	other_ranges <- strat_ranges[!(1:notu) %in% n,];
	# get taxa known after taxon first appears
	overlapping_set <- subset(other_ranges,other_ranges$LAD>=strat_ranges$FAD[n]);
	overlapping_set <- subset(overlapping_set,overlapping_set$FAD<=strat_ranges$LAD[n]);
	if (nrow(overlapping_set)==0)	{
		older_set <- subset(other_ranges,other_ranges$FAD<=strat_ranges$LAD[n]);
		younger_set <- subset(other_ranges,other_ranges$LAD>=strat_ranges$FAD[n]);
		if (nrow(younger_set)>0 && min(younger_set$FAD) > strat_ranges$LAD[n])	{
			gap_pst[n] <- min(younger_set$FAD) - strat_ranges$LAD[n]
#			min_gap <- min_gap + min(younger_set$FAD) - strat_ranges$LAD[n];
			}
		if (nrow(older_set)>0 && strat_ranges$FAD[n] > max(older_set$LAD))	{
			gap_pre[n] <- strat_ranges$FAD[n] - max(older_set$LAD);
#			min_gap <- min_gap + (strat_ranges$FAD[n] - max(older_set$LAD));
			}
		}
	}
#min_gap <- sum(gap_pre) + sum(gap_pst);
min_gap <- sum(gap_pre);
return(min_gap);
}

maximum_gaps <- function(strat_ranges)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
return(sum(strat_ranges$FAD-min(strat_ranges$FAD)))
}