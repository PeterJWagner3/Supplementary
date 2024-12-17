#library(ape)
#library(paleotree)
#library(phangorn)
#library(phylobase)
#library(gtools)

clear_na_from_matrix <- function(data, replacement)  {
size <- dim(data)
for (i in 1:size[1])	{
	for (j in 1:size[2]) if (is.na(data[i,j]))	data[i,j] <- replacement
	}
return(data)
}

clear_matrix_na_with_another_cell_value <- function(data,j, k)	{
size <- dim(data)
for (i in 1:size[1])	{
	if (is.na(data[i,j]))	data[i,j] <- data[i,k]
	}
return(data)
}

clear_na_from_vector <- function(data, replacement)	{
size <- length(data)
for (i in 1:size[1])	if (is.na(data[i]))	data[i] <- replacement
return(data)
}

remove_elements_from_vector <- function(to_remove,orig_vector)	{
return(orig_vector[! orig_vector %in% to_remove])
}

replace_elements_from_vector <- function(to_remove,to_replace,orig_vector,keep_length=TRUE)	{
reduced <- remove_elements_from_vector(to_remove,orig_vector)
i <- sum(reduced>0)+1
j <- sum(reduced>0)+length(to_replace)
reduced[i:j] <- to_replace
if (keep_length)	{
	fortify <- length(orig_vector)-length(reduced)
	reduced <- c(reduced,rep(0,fortify))
	}
return(reduced)
}

accersi_nodes_from_base_for_paleotree <- function(record) {
# describe how far above the base of the tree each evolved taxon is
ANCESTORID <- match("ancestor.id",colnames(record))
ANCESTOR <- match("ancestral",colnames(record))
if (is.na(ANCESTOR))	{
	ancestral <- vector(length=length(record[,1]))
	for (i in 2:length(record[,1]))	ancestral[record[i,2]] <- ancestral[record[i,2]]+1
	record <- cbind(record,ancestral)
	ANCESTOR <- match("ancestral",colnames(record))
	}
raw_taxa <- dim(record)[1]
#record[,ANCESTORID]
noderanks <- vector(length=raw_taxa)
for (i in 1:raw_taxa)	noderanks[record[i,ANCESTORID]] <- 1
# this gives the node rank for taxa descended from this species
for (i in 2:raw_taxa)	noderanks[i] <- noderanks[i]+noderanks[record[i,ANCESTORID]]
return(noderanks)
}

#nodestart <- rank(record[,2],ties.method=c("min"))	# this gives the ranked numbers for raw_ancestral species
accersi_ape_tree_from_paleo_tree <- function(record,sampled_only=T)	{
#ANCESTORID <- match("ancestor.id",colnames(record))
#ORIG <- match("orig.time",colnames(record))
#EXTN <- match("ext.time",colnames(record))
#FAD <- match("FAD",colnames(record))
#x <- rank(record[,FAD],na.last=FALSE)
x <- rank(record$FAD,na.last=FALSE)
raw_taxa <- length(x)
LAD <- match("LAD",colnames(record))
RANK_FAD <- match("rank_FAD",colnames(record))
sampled_taxa <- sum(record$FAD>0)
record$RANK_FAD <- 1+nrow(record)-rank(record$FAD,na.last=FALSE);

OBSRANGE <- match("obsrange",colnames(record))
if (is.na(OBSRANGE))
	obsrange <- abs(record$FAD-record$LAD);

for (sp in 1:nrow(record))
	ancestral <- c(ancestral,sum(record$ancestor.id==sp));	

noderanks <- vector(length=raw_taxa);
for (i in 1:dim(record)[1])		noderanks[record[i,ANCESTORID]] <- 1;
for (i in 2:length(noderanks))	noderanks[i] <- noderanks[i]+noderanks[i-1];	# this gives the node rank for taxa descended from this species
ond <- raw_taxa+1;
noderanks <- noderanks+raw_taxa;
nodeanc <- vector(length=max(noderanks));
ancnode <- vector(length=raw_taxa);
for (i in 1:nrow(record))	{
	if (record$ancestor.id[i]>0)	{
		nodeanc[noderanks[i]] <- i
		ancnode[i] <- noderanks[i]
		}
	}

# now, used rank order of appearances to reduce the tree
appearances <- (1+raw_taxa)-rank(record[,FAD])
ss <- 1
sampled <- vector(length=sampled_taxa)
for (i in 1:raw_taxa)  {
	if (appearances[i]<=sampled_taxa)	{
		sampled[ss] <- i
		ss <- ss+1
		}
	}

venn_ape <- matrix(0,length(nodeanc),sampled_taxa)
node_rich <- vector(length=length(nodeanc))
mxnd <- 0
for (s in 1:sampled_taxa)  venn_ape[(raw_taxa+1),s] <- sampled[s]
node_rich[raw_taxa+1] <- sampled_taxa
for (s in 1:sampled_taxa)  {
	sp <- sampled[s]
	if (record[sp,ANCESTORID]>0) {
		anc_sp <- record[sp,ANCESTORID]      # get ancestors_sampled
		anc_nd <- ancnode[anc_sp]   # get ancestral node
		if (anc_nd>mxnd)  mxnd <- anc_nd
		# this loop will add a sampled species to each node where it has an ancestor.  If we miss multiple
		# ancestors, then it will be in multiple rows.
		if (anc_nd>ond)	{
			node_rich[anc_nd] <- node_rich[anc_nd]+1
			venn_ape[anc_nd,node_rich[anc_nd]] <- sp
			while (anc_nd>(ond+1)) {
				anc_sp <- record[anc_sp,ANCESTORID]
				anc_nd <- ancnode[anc_sp]
				if (anc_nd>(raw_taxa+1))	{
					node_rich[anc_nd] <- node_rich[anc_nd]+1
					venn_ape[anc_nd,node_rich[anc_nd]] <- sp
					}
				}	# add it to lower nodes
			}		# add taxon to node
		if (ancnode[sp]!=0)	{
			new_nd <- ancnode[sp]
			node_rich[new_nd] <- 1
			venn_ape[new_nd,1] <- sp
			if (new_nd>mxnd)	mxnd <- new_nd
			}	# end establishing new node if we find an ancestor
		}	# if this a case about which we need to worry
	}

ttl_nodes <- 1+(mxnd-ond)
node_desc <- node_rich
raw_ape <- venn_ape
# reduce daughter clades to htu number.
for (i in mxnd:(ond+1))  {
	if (node_desc[i]>1)	{
		NTU <- i	# NTU is node taxon unit
		for (j in (i-1):1)	{
			if (node_desc[j]>=node_desc[i])	{
				if (all(match(subset(raw_ape[i,], raw_ape[i,] > 0), raw_ape[j,],nomatch=FALSE)))	{   
					row.temp <- subset(raw_ape[j, ], !raw_ape[j, ] %in% subset(raw_ape[i,], raw_ape[i,] >0) )
					row.temp.2 <- c(subset(row.temp, row.temp > 0), NTU )
					row.final <- c(row.temp.2, rep(0, sampled_taxa - length(row.temp.2)) )
					raw_ape[j, ] <- row.final
					node_desc[j] <- 1+node_desc[j]-node_desc[i]
					}	# end case where derived clade members are found in an older clade
				}
			}
		}
	}

# get taxon numbers for nodes to rank and renumber
live_nd <- 0	#number of relevant nodes
for (n1 in (raw_taxa+1):mxnd)	if (node_desc[n1]>1)	live_nd <- live_nd+1
node_nums <- vector(length=live_nd)
samp_anc <- vector(length=max(sampled))
n <- 1
for (n1 in (raw_taxa+1):mxnd)	{
	if (node_desc[n1]>1)	{
		node_nums[n] <- n1
		n <- n+1
		}
	anc <- nodeanc[n1]
	if (record[anc,FAD]>0)	samp_anc[anc] <- n1	# this sampled species is the ancestor of node n1
	}   # get raw node numbers that have 2+ descendants

dummy <- c(sampled,node_nums)
#order(dummy)
taxon_ranks <- vector(length=mxnd)
taxon_ranks[dummy] <- order(dummy)
nodes <- length(node_nums)

ape <- vector(length=max(taxon_ranks))
ape[sampled_taxa+1] <- 0
for (n in (raw_taxa+1):mxnd) {
	if (node_desc[n]>1) {
		anc <- taxon_ranks[n]
		for (f in 1:node_desc[n]) {
			sp <- raw_ape[n,f]
			rnk <- taxon_ranks[sp]
			ape[rnk] <- anc
			}
		}
	}

# first time we do raw_ancestral.2
final_ancestral <- vector(length=sampled_taxa)
node_ancestors <- vector(length=length(ape))
for (c in ond:mxnd)	{
	if (node_desc[c]>1 && raw_ape[c,1]<=sampled_taxa)	{
		sp <- raw_ape[c,1]
		s <- taxon_ranks[sp]
		if (ancnode[sp]==c)	{
			final_ancestral[s] <- 1
			node_ancestors[ape[s]] <- s
			}
		}
	}

venn_tree <- matrix(0,live_nd,sampled_taxa)
c <- 0
for (i in ond:mxnd)	{
	if (node_desc[i]>1)	{
		c <- c+1
		for (j in 1:node_rich[i])	venn_tree[c,j] <- taxon_ranks[venn_ape[i,j]]
		}
	}
# now, get branch lengths.  Then renumber things
raw_branchings <- vector(length=(mxnd+1))
raw_divergences <- vector(length=(mxnd+1))
divergences <- matrix(0,(sampled_taxa+live_nd),2)
branchings <- vector(length=(sampled_taxa+live_nd))
# get the divergence times of the nodes by getting the origin of the last common ancestor
for (n1 in 1:live_nd)	{
	branchings[sampled_taxa+n1] <- raw_branchings[node_nums[n1]] <- 1
	divergences[(sampled_taxa+n1),1] <- raw_divergences[node_nums[n1]] <- record[nodeanc[node_nums[n1]],ORIG]
	}

#now, set aside observed ranges, which form the basis for divergences for otus
ranges <- matrix(0,sampled_taxa,2)
for (s in 1:sampled_taxa)	{
	sp <- sampled[s]
	divergences[s,2] <- ranges[s,1] <- record[sp,FAD]
	ranges[s,2] <- record[sp,LAD]
	if (final_ancestral[s]==0)	{
		branchings[s] <- raw_branchings[sp] <- 1
		divergences[s,1] <- raw_divergences[sp] <- record[sp,ORIG]
		} else {
		divergences[ape[s],2] <- ranges[s,1]
		raw_divergences[sp] <- record[sp,ORIG]
		divergences[s,1] <- divergences[s,2]
		}
	}

# add further unsampled ancestors 
max_obs_sp <- max(sampled)
for (n1 in mxnd:(raw_taxa+1)) {
	if (node_desc[n1]==1) {
		sp <- raw_ape[n1,1]
		s <- taxon_ranks[sp]
		if (sp>max_obs_sp || (sp<=max_obs_sp && ancnode[sp]!=n1))	{
			branchings[s] <- raw_branchings[sp] <- raw_branchings[sp]+1		# add branch
			divergences[s,1] <- raw_divergences[sp] <- record[nodeanc[n1],ORIG]		# push divergence times down
			}
		}
	}

#make sure that nodes without sampled ancestors have upper divergence based on latest divergence time of a descendant
fred <- matrix(0,2,length(ape))
for (i in 1:length(ape)) fred[1,i] <- i
fred[2,] <- ape
for (htu in max(ape):(sampled_taxa+1))	{
#	htu <- i+sampled_taxa
	if (divergences[htu,2]==0)	{
		scions <- subset(fred[1,],fred[2,]==htu)
#		scions <- subset(venn_tree[nd,],venn_tree[nd,]>0)
#		divergences[htu,2] <- max(divergences[,1])
		divergences[htu,2] <- max(divergences[scions,1])
#		for (s in 1:sampled_taxa)	{
#			if (ape[s]==htu && divergences[s,1]<divergences[htu,2])
#				divergences[htu,2] <- divergences[s,1]
#			}
		if (htu<max(ape))	{
			for (s in (htu+1):max(ape))	{
				if (ape[s]==htu && divergences[s,1]<divergences[htu,2])
					divergences[htu,2] <- divergences[s,1]
				}
			}
		}
	}

branchings[sampled_taxa+1] <- 0	# don't bother with a branch length for the basal most node
output <- list(ape,branchings,divergences,venn_tree,ranges,final_ancestral,node_ancestors)
return(output)
}

get_ape_tree_from_paleo_tree_set_taxa <- function(record,sampled_taxa)	{
ANCESTORID <- match("ancestor.id",colnames(record))
ORIG <- match("orig.time",colnames(record))
EXTN <- match("ext.time",colnames(record))
FAD <- match("FAD",colnames(record))
x <- rank(record[,FAD],na.last=FALSE)
raw_taxa <- length(x)
LAD <- match("LAD",colnames(record))
RANK_FAD <- match("rank_FAD",colnames(record))
if (is.na(RANK_FAD))	{
	x <- rank(record[,FAD],na.last=FALSE)
	raw_taxa <- length(x)
	rank_FAD <- vector(length=raw_taxa)
	for (i in 1:raw_taxa)  rank_FAD[i] <- 1+raw_taxa-x[i]
	record <- cbind(record,rank_FAD)
	RANK_FAD <- match("rank_FAD",colnames(record))
	}
OBSRANGE <- match("obsrange",colnames(record))
if (is.na(OBSRANGE))  {
	obsrange <- vector(length=length(record[,1]))
	for (i in 1:length(record[,1]))	obsrange[i] <- abs(record[i,FAD]-record[i,LAD])
	}
ANCESTOR <- match("ancestral",colnames(record))
if (is.na(ANCESTOR))	{
	ancestral <- vector(length=length(record[,1]))
	for (i in 2:length(record[,1]))	ancestral[record[i,2]] <- ancestral[record[i,2]]+1
	record <- cbind(record,ancestral)
	ANCESTOR <- match("ancestral",colnames(record))
	}

noderanks <- vector(length=raw_taxa)
for (i in 1:dim(record)[1])		noderanks[record[i,ANCESTORID]] <- 1
for (i in 2:length(noderanks))	noderanks[i] <- noderanks[i]+noderanks[i-1]	# this gives the node rank for taxa descended from this species
ond <- raw_taxa+1
noderanks <- noderanks+raw_taxa
nodeanc <- vector(length=max(noderanks))
ancnode <- vector(length=raw_taxa)
for (i in 1:length(record[,ANCESTORID]))	{
	if (record[i,ANCESTOR]>0)	{
		nodeanc[noderanks[i]] <- i
		ancnode[i] <- noderanks[i]
		}
	}

# now, used rank order of appearances to reduce the tree
appearances <- (1+raw_taxa)-rank(record[,FAD])
ss <- 1
sampled <- vector(length=sampled_taxa)
for (i in 1:raw_taxa)  {
	if (appearances[i]<=sampled_taxa)	{
		sampled[ss] <- i
		ss <- ss+1
		}
	}

venn_ape <- matrix(0,length(nodeanc),sampled_taxa)
node_rich <- vector(length=length(nodeanc))
mxnd <- 0
for (s in 1:sampled_taxa)  venn_ape[(raw_taxa+1),s] <- sampled[s]
node_rich[raw_taxa+1] <- sampled_taxa
for (s in 1:sampled_taxa)  {
	sp <- sampled[s]
	if (record[sp,ANCESTORID]>0) {
		anc_sp <- record[sp,ANCESTORID]      # get ancestors_sampled
		anc_nd <- ancnode[anc_sp]   # get ancestral node
		if (anc_nd>mxnd)  mxnd <- anc_nd
		# this loop will add a sampled species to each node where it has an ancestor.  If we miss multiple
		# ancestors, then it will be in multiple rows.
		if (anc_nd>ond)	{
			node_rich[anc_nd] <- node_rich[anc_nd]+1
			venn_ape[anc_nd,node_rich[anc_nd]] <- sp
			while (anc_nd>(ond+1)) {
				anc_sp <- record[anc_sp,ANCESTORID]
				anc_nd <- ancnode[anc_sp]
				if (anc_nd>(raw_taxa+1))	{
					node_rich[anc_nd] <- node_rich[anc_nd]+1
					venn_ape[anc_nd,node_rich[anc_nd]] <- sp
					}
				}	# add it to lower nodes
			}		# add taxon to node
		if (ancnode[sp]!=0)	{
			new_nd <- ancnode[sp]
			node_rich[new_nd] <- 1
			venn_ape[new_nd,1] <- sp
			if (new_nd>mxnd)	mxnd <- new_nd
			}	# end establishing new node if we find an ancestor
		}	# if this a case about which we need to worry
	}

ttl_nodes <- 1+(mxnd-ond)
node_desc <- node_rich
raw_ape <- venn_ape
# reduce daughter clades to htu number.
for (i in mxnd:(ond+1))  {
	if (node_desc[i]>1)	{
		NTU <- i	# NTU is node taxon unit
		for (j in (i-1):1)	{
			if (node_desc[j]>=node_desc[i])	{
				if (all(match(subset(raw_ape[i,], raw_ape[i,] > 0), raw_ape[j,],nomatch=FALSE)))	{   
					row.temp <- subset(raw_ape[j, ], !raw_ape[j, ] %in% subset(raw_ape[i,], raw_ape[i,] >0) )
					row.temp.2 <- c(subset(row.temp, row.temp > 0), NTU )
					row.final <- c(row.temp.2, rep(0, sampled_taxa - length(row.temp.2)) )
					raw_ape[j, ] <- row.final
					node_desc[j] <- 1+node_desc[j]-node_desc[i]
					}	# end case where derived clade members are found in an older clade
				}
			}
		}
	}

# get taxon numbers for nodes to rank and renumber
live_nd <- 0	#number of relevant nodes
for (n1 in (raw_taxa+1):mxnd)	if (node_desc[n1]>1)	live_nd <- live_nd+1
node_nums <- vector(length=live_nd)
samp_anc <- vector(length=max(sampled))
n <- 1
for (n1 in (raw_taxa+1):mxnd)	{
	if (node_desc[n1]>1)	{
		node_nums[n] <- n1
		n <- n+1
		}
	anc <- nodeanc[n1]
	if (record[anc,FAD]>0)	samp_anc[anc] <- n1	# this sampled species is the ancestor of node n1
	}   # get raw node numbers that have 2+ descendants

dummy <- c(sampled,node_nums)
#order(dummy)
taxon_ranks <- vector(length=mxnd)
taxon_ranks[dummy] <- order(dummy)
nodes <- length(node_nums)

ape <- vector(length=max(taxon_ranks))
ape[sampled_taxa+1] <- 0
for (n in (raw_taxa+1):mxnd) {
	if (node_desc[n]>1) {
		anc <- taxon_ranks[n]
		for (f in 1:node_desc[n]) {
			sp <- raw_ape[n,f]
			rnk <- taxon_ranks[sp]
			ape[rnk] <- anc
			}
		}
	}

# first time we do raw_ancestral.2
final_ancestral <- vector(length=sampled_taxa)
node_ancestors <- vector(length=length(ape))
for (c in ond:mxnd)	{
	if (node_desc[c]>1 && raw_ape[c,1]<=sampled_taxa)	{
		sp <- raw_ape[c,1]
		s <- taxon_ranks[sp]
		if (ancnode[sp]==c)	{
			final_ancestral[s] <- 1
			node_ancestors[ape[s]] <- s
			}
		}
	}

venn_tree <- matrix(0,live_nd,sampled_taxa)
c <- 0
for (i in ond:mxnd)	{
	if (node_desc[i]>1)	{
		c <- c+1
		for (j in 1:node_rich[i])	venn_tree[c,j] <- taxon_ranks[venn_ape[i,j]]
		}
	}
# now, get branch lengths.  Then renumber things
raw_branchings <- vector(length=(mxnd+1))
raw_divergences <- vector(length=(mxnd+1))
divergences <- matrix(0,(sampled_taxa+live_nd),2)
branchings <- vector(length=(sampled_taxa+live_nd))
# get the divergence times of the nodes by getting the origin of the last common ancestor
for (n1 in 1:live_nd)	{
	branchings[sampled_taxa+n1] <- raw_branchings[node_nums[n1]] <- 1
	divergences[(sampled_taxa+n1),1] <- raw_divergences[node_nums[n1]] <- record[nodeanc[node_nums[n1]],ORIG]
	}

#now, set aside observed ranges, which form the basis for divergences for otus
ranges <- matrix(0,sampled_taxa,2)
for (s in 1:sampled_taxa)	{
	sp <- sampled[s]
	divergences[s,2] <- ranges[s,1] <- record[sp,FAD]
	ranges[s,2] <- record[sp,LAD]
	if (final_ancestral[s]==0)	{
		branchings[s] <- raw_branchings[sp] <- 1
		divergences[s,1] <- raw_divergences[sp] <- record[sp,ORIG]
		} else {
		divergences[ape[s],2] <- ranges[s,1]
		raw_divergences[sp] <- record[sp,ORIG]
		divergences[s,1] <- divergences[s,2]
		}
	}

# add further unsampled ancestors 
max_obs_sp <- max(sampled)
for (n1 in mxnd:(raw_taxa+1)) {
	if (node_desc[n1]==1) {
		sp <- raw_ape[n1,1]
		s <- taxon_ranks[sp]
		if (sp>max_obs_sp || (sp<=max_obs_sp && ancnode[sp]!=n1))	{
			branchings[s] <- raw_branchings[sp] <- raw_branchings[sp]+1		# add branch
			divergences[s,1] <- raw_divergences[sp] <- record[nodeanc[n1],ORIG]		# push divergence times down
			}
		}
	}

#make sure that nodes without sampled ancestors have upper divergence based on latest divergence time of a descendant
fred <- matrix(0,2,length(ape))
for (i in 1:length(ape)) fred[1,i] <- i
fred[2,] <- ape
for (htu in max(ape):(sampled_taxa+1))	{
#	htu <- i+sampled_taxa
	if (divergences[htu,2]==0)	{
		scions <- subset(fred[1,],fred[2,]==htu)
#		scions <- subset(venn_tree[nd,],venn_tree[nd,]>0)
#		divergences[htu,2] <- max(divergences[,1])
		divergences[htu,2] <- max(divergences[scions,1])
#		for (s in 1:sampled_taxa)	{
#			if (ape[s]==htu && divergences[s,1]<divergences[htu,2])
#				divergences[htu,2] <- divergences[s,1]
#			}
		if (htu<max(ape))	{
			for (s in (htu+1):max(ape))	{
				if (ape[s]==htu && divergences[s,1]<divergences[htu,2])
					divergences[htu,2] <- divergences[s,1]
				}
			}
		}
	}

branchings[sampled_taxa+1] <- 0	# don't bother with a branch length for the basal most node
output <- list(ape,branchings,divergences,venn_tree,ranges,final_ancestral,node_ancestors)
return(output)
}

scramble_multistates <- function(nstates)	{
dstates <- vector(length=nstates)
for (i in 1:nstates)	dstates[i] <- i-1
for (i in 1:nstates)	{
	p <- i+ceiling((nstates-i)*runif(1))
	q <- dstates[i]
	dstates[i] <- dstates[p]
	dstates[p] <- q
	}
return(dstates)
}

flip_binaries <- function(ch_vector)	{
new_vector <- ch_vector;
new_vector[ch_vector==0] <- 1;
new_vector[ch_vector==1] <- 0;
return(new_vector);
}

evolve_X_discrete_steps <- function(venn_tree,sampled_taxa,nodes,branchings,nchars,states,types,steps,UNKNOWN,INAP)	{
simchmatrix <- matrix(0,sampled_taxa,nchars)
branches <- sum(branchings)
ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=sampled_taxa+nodes)
char_changes <- vector(length=nchars)
br <- 1		#for debugging
for (b in 1:length(branchings))	{
	bb <- br
	if (branchings[b]>0)	for (br in bb:(bb+branchings[b]))	ab[br] <- b
	}
fb <- length(unique(ab))		# free branches (i.e., those that can change 1+ times)

# first, make sure that all states appear
for (ch in 1:nchars)	{
	a <- 0
	while (a<(states[ch]-1))	{
		brdelts <- ceiling(fb*runif(4*states[ch]-1))
		brdelts <- unique(brdelts)
		a <- length(brdelts)
		}		# make sure that we have enough unique branches to evolve every state
	for (c in 1:(states[ch]-1))	{
		b <- brdelts[c]
		br <- ab[b]
		if (br<=sampled_taxa)	{
#			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch!=INAP)	simchmatrix[br,ch] <- dstates[simchmatrix[br,ch]]
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	simchmatrix[br,ch] <- c
			} else {
			ndb <- br-sampled_taxa
			for (d in 1:length(subset(venn_tree[ndb,],venn_tree[ndb,]>0)))	{
				sp <- venn_tree[ndb,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- c
				}	
			}
		char_changes[ch] <- 1+char_changes[ch]
		branch_changes[br] <- branch_changes[br]+1
		}
	}

for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=sampled_taxa)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndc <- br-sampled_taxa
				for (d in 1:length(subset(venn_tree[ndc,],venn_tree[ndc,]>0)))	{
					sp <- venn_tree[ndc,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		}	# end case where branch needed changes
	}
delta <- sum(char_changes)
for (d in delta:steps)	{
	b <- ceiling(fb*runif(1))
	br <- ab[b]
	ch <- ceiling(nchars*runif(1))
	if (states[ch]==2)	{
		dstates <- c(1,0)
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	
		} 
	if (br<=sampled_taxa)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			nde <- br-sampled_taxa
			for (d in 1:length(subset(venn_tree[nde,],venn_tree[nde,]>0)))	{
				sp <- venn_tree[nde,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	}

output <- list(simchmatrix,branch_changes,char_changes)
}

evolve_1_discrete_step <- function(simchmatrix,branch_changes,char_changes,venn_tree,sampled_taxa,nodes,branchings,nchars,states,types,steps,UNKNOWN,INAP)	{
# char_changes: number of changes for each character
# branch_changes: number of changes on each branch
branches <- sum(branchings)
ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
br <- 1		#for debugging
for (b in 1:length(branchings))	{
	bb <- br
	if (branchings[b]>0)	for (br in bb:(bb+branchings[b]))	ab[br] <- b
	}
fb <- length(unique(ab))		# free branches

# first, make sure that all states appear
b <- ceiling(fb*runif(1))
br <- ab[b]
ch <- ceiling(nchars*runif(1))
if (states[ch]==2)	{
	dstates <- c(1,0)
	} else if (types[ch]==0)	{
	dstates <- scramble_multistates(states[ch])	
	} 
if (br<=sampled_taxa)	{
	if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
		simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
		} else {
		ndf <- br-sampled_taxa
		for (d in 1:length(subset(venn_tree[ndf,],venn_tree[ndf,]>0)))	{
			sp <- venn_tree[ndf,d]
			if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
			}	
		}
	#change_character_on_branch(ch,branchings[b])
	}
branch_changes[br] <- branch_changes[br]+1
char_changes[ch] <- char_changes[ch]+1

output <- list(simchmatrix,branch_changes,char_changes)
}

compatibilty_over_range_of_steps <- function(venn_tree,sampled_taxa,nodes,branchings,nchars,states,types,minsteps,maxsteps,UNKNOWN,INAP)	{
simchmatrix <- matrix(0,sampled_taxa,nchars)
branches <- sum(branchings)
ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=sampled_taxa+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=maxsteps)
br <- 0		#for debugging
for (b in 1:length(branchings))	{
	if (branchings[b]>0)	{
	    for (bb in 1:branchings[b])	{
	      br <- br+1
	      ab[br] <- b
	      }
	  }
	}

tb <- length(ab)
unique_ab <- unique(ab)
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:sampled_taxa)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

nodes <- length(branchings)-sampled_taxa
rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(chmatrix,nchars,sampled_taxa,branchings,rich,ab) # possible available branches matrix
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
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=sampled_taxa)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-sampled_taxa
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

# second, make sure that all states appear
for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=sampled_taxa)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-sampled_taxa
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
simcompmat <- compatibility_matrix(simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
character_compats <- vector(length=nchars)
for (c1 in 1:nchars)	for (c2 in 1:nchars) if (c1!=c2)	character_compats[c1] <- character_compats[c1]+simcompmat[c1,c2]
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
#print(c(delta,simcompat[delta]))
for (d in (delta+1):maxsteps)	{
	ch <- ceiling(nchars*runif(1))
	while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
  c <- char_changes[ch]+1
	br <- pabm[ch,c]

	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	if (br<=sampled_taxa)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-sampled_taxa
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
	# update compatibility matrix
	simcompmat[ch,] <- new_compat
	simcompmat[,ch] <- new_compat
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	#print(c(d,simcompat[d]))

	# update character compatibilities that might have been altered
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
return (simcompat)
}

evolve_to_particular_compatibilty_old <- function(ttl_compat,venn_tree,sampled_taxa,nodes,branchings,nchars,states,types,minsteps,maxsteps,UNKNOWN,INAP,rep)	{
simchmatrix <- matrix(0,sampled_taxa,nchars)
branches <- sum(branchings)
ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=sampled_taxa+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=maxsteps)
br <- 0		#for debugging
for (b in 1:length(branchings))	{
	if (branchings[b]>0)	{
		for (bb in 1:branchings[b])	{
			br <- br+1
			ab[br] <- b
			}
		}
	}

tb <- length(ab)
unique_ab <- unique(ab)
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:sampled_taxa)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

nodes <- length(branchings)-sampled_taxa
rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(chmatrix,nchars,sampled_taxa,branchings,rich,ab) # possible available branches matrix
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
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=sampled_taxa)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-sampled_taxa
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

# second, make sure that all states appear
for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=sampled_taxa)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-sampled_taxa
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
simcompmat <- compatibility_matrix(simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
character_compats <- vector(length=nchars)
for (c1 in 1:nchars)	for (c2 in 1:nchars) if (c1!=c2)	character_compats[c1] <- character_compats[c1]+simcompmat[c1,c2]
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
if (rep>0)	print(c(rep,delta,simcompat[delta]))

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
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	if (br<=sampled_taxa)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-sampled_taxa
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,nchars,states,types,sampled_taxa,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (rep>0 && counter%%10==0)	print(c(rep,d,simcompat[d]))

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
output <- list(d,simcompmat)
return (output)
}

get_possible_branches_for_character <- function(chmatrix,ch,sampled_taxa,branchings,rich,ab)	{
poss_brs <- length(ab)
pab <- vector(length=poss_brs)
a <- 1
for (b in 1:poss_brs)	{
	if (ab[b]<=sampled_taxa)	{
		s <- ab[b]
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	{
#			for (z in 1:branchings[s])	{
				pab[a] <- s
				a <- a+1
#				}	# add each speciation for the branch
			}	# case where species is scored
#		b <- b+(branchings[s]-1)
		} else {
		ndx <- ab[b]-sampled_taxa
		accept <- 0
		for (r in 1:rich[ndx])	{
			s <- venn_tree[ndx,r]
			if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	{
				accept <- 1
				r <- rich[ndx]
				}	# case where a clade member is scored
			}
		if (accept==1)	{
#			for (z in 1:branchings[nd+sampled_taxa])	{
				pab[a] <- ndx+sampled_taxa
				a <- a+1
#				}	# add each speciation
			}
#		b <- b+(branchings[nd+sampled_taxa]-1)
		}	# end case of node
	}
return(pab)
}

get_possible_branches_for_all_characters_old <- function(chmatrix,nchars,sampled_taxa,branchings,rich,ab)	{
poss_brs <- length(ab)
pabm <- matrix(0,nchars,poss_brs)
for (ch in 1:nchars)	{
	a <- 1
	for (b in 1:poss_brs)	{
		if (ab[b]<=sampled_taxa)	{
			s <- ab[b]
			if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	{
				pabm[ch,a] <- s
				a <- a+1
			  }	# case where species is scored
			 } else {
			  nda <- ab[b]-sampled_taxa
			  accept <- 0
			  for (r in 1:rich[nda])	{
				  s <- venn_tree[nda,r]
				  if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	{
				  	accept <- 1
				  	r <- rich[nda]
				  	}	# case where a clade member is scored
				  }
			  if (accept==1)	{
			  	pabm[ch,a] <- nda+sampled_taxa
			  	a <- a+1
			  }
		  }	# end case of node
    }
  }
return(pabm)
}

gower_transformation <- function(distmat, n)	{
dmat <- matrix(0,n,n)
for (i in 1:n)  for (j in 1:n)  dmat[i,j]=-0.5*distmat[i,j]
  
gower <- matrix(0,n,n)
average <- vector(length=n)
  
for (i in 1:n)	average[i] <- mean(dmat[i,])

ave=mean(dmat)
  
for (i in 1:(n-1))	for (j in (i+1):n)	gower[i,j] <- gower[j,i] <- distmat[i,j]-(average[i]+average[j])+ave
  
return (gower)
}

evolve_clade_for_limited_time <- function(lambda,mu,limit)	{
ttl <- total_progeny <- 1;
history <- matrix(0,1,5);
history[1,2] <- min(-log(runif(1))/mu,limit);
history[1,3] <- anc <- 0;
clade <- vector(length=1);
clade[1] <- 1;
while (anc<total_progeny)	{
	anc <- anc+1
	lifetime <- history[anc,2]-history[anc,1]
	history[anc,4] <- daughters <- qpois(runif(1),lambda*lifetime)
	history[anc,5] <- dpois(daughters,lambda*lifetime)
	if (daughters>0)	{
		newbies <- vector(length=daughters)
		node_hist <- matrix(0,daughters,5)
		for (n in 1:daughters)	{
			newbies[n] <- ttl+n
			node_hist[n,1] <- history[anc,1]+(runif(1)*(history[anc,2]-history[anc,1]))
			life <- -log(runif(1))/mu
			node_hist[n,2] <- min((node_hist[n,1]+life),limit)
			node_hist[n,3] <- anc
			}
		history <- rbind(history,node_hist)
		total_progeny <- total_progeny+daughters
		}
	}
colnames(history) <- c("birth","death","ancestor","descendants","lambda[no desc]");
return(data.frame(history));
}

Wagner_sampling <- function(durations,fr,end_date=0,missed_fa=-1,missed_la=-2,cutoff=0)	{
taxa <- nrow(durations);
durations[durations[,2]<cutoff,2] <- cutoff
dur <- abs(durations[,1]-durations[,2])
finds <- vector(length=taxa)
wag_samples <- matrix(missed_fa,taxa,2)
if (missed_fa!=missed_la)	wag_samples[,2] <- missed_la
colnames(wag_samples) <- c("FAD","LAD")
mx_finds <- 1;
all_finds <- array(0,dim=c(taxa,mx_finds));
for (s in 1:taxa)	{
	finds[s] <- rpois(1,fr*dur[s])
	if (finds[s]>0)	{
		if (finds[s]>mx_finds)	{
			added <- finds[s] - mx_finds;
			dummy <- array(0,dim=c(taxa,added));
			all_finds <- cbind(all_finds,dummy);
			mx_finds <- finds[s];
			}
		dist_finds <- sort(runif(finds[s]));
		all_finds[s,1:finds[s]] <- dist_finds;
		wag_samples[s,1] <- durations[s,1]-(min(dist_finds)*dur[s]);
		wag_samples[s,2] <- durations[s,1]-(max(dist_finds)*dur[s]);
		}	
	}
wag_ranges <- data.frame(FAD=as.numeric(wag_samples[,1]),LAD=as.numeric(wag_samples[,2]),stringsAsFactors = F);
rownames(wag_ranges) <- 1:taxa;
output <- list(wag_ranges,all_finds);
names(output) <- c("strat_ranges_all","finds_all");
return(output)
}

accersi_otu_numbers_from_paleo_tree_sampled_taxa <- function(record)	{
ORIGID <- match("taxon.id",colnames(record))
ANCESTORID <- match("ancestor.id",colnames(record))
FAD <- match("FAD",colnames(record))
otu_orig_nos <- subset(record,record[,FAD]>0)[,ORIGID]
names(otu_orig_nos) <- rep(NULL,length(otu_orig_nos))
return(otu_orig_nos)
}

# 2017-01-12 rewritten yet again....
accersi_rawest_venn_tree_from_paleo_tree <- function(record)	{
# Again, remember that we are jumping between raw evolved number & OTU numbers
ANCESTORID <- match("ancestor.id",colnames(record))
FAD <- match("FAD",colnames(record))
# number of sampled taxa
notus <- dim(subset(record,record[,FAD]>0))[1]

# get the number of descendants for each taxon & number of taxa with descendants
# 		success of simulation numbered species
ancestral_success <- accersi_evolutionary_success(record)
# simulation numbered ancestral species
true_ancestors <- as.numeric(names(ancestral_success))
# number of nodes in raw tree
raw_nodes <- length(true_ancestors)
raw_venn <- matrix(0,raw_nodes,notus)
raw_venn_richness <- vector(length=raw_nodes)
# now, get the OTU numbers corresponding to simulation numbers
# 	THIS IS OUR ROSETTA STONE!!!
taxon_numbering_rosetta <- accersi_otu_numbers_from_paleo_tree_sampled_taxa(record)
for (n in 1:notus)	{
	# make sure that ancestors are in their own nodes
	spc <- taxon_numbering_rosetta[n]	# sim number corresponding to OTU number
	if (is.na(match(spc,true_ancestors)))	{
		anc <- record[spc,ANCESTORID]
		}	else	{
		anc <- spc	
		}
	node <- match(anc,true_ancestors)
	raw_venn[node,1+sum(raw_venn[node,]>0)] <- taxon_numbering_rosetta[n]
	raw_venn_richness[node] <- raw_venn_richness[node]+1
	while (node>1)	{
		spc <- anc
		anc <- record[spc,ANCESTORID]	# get next node down
		node <- match(anc,true_ancestors)
		raw_venn[node,1+sum(raw_venn[node,]>0)] <- taxon_numbering_rosetta[n]
		raw_venn_richness[node] <- raw_venn_richness[node]+1
		}
	}
# Name Rows by the simulated species number of the ancestor
rownames(raw_venn) <- true_ancestors
# get the relevant nodes: i.e., those with 2+ observed species
relevant_raw_nodes <- (1:raw_nodes)[raw_venn_richness>0]
rawest_venn_tree <- raw_venn[relevant_raw_nodes,]
rawest_venn_richness <- raw_venn_richness[relevant_raw_nodes]
#write.table(rawest_venn_tree,file="Rawest_Venn_Tree.xls",sep="\t",row.names = TRUE)
# now, eliminate redundant branches
return(rawest_venn_tree)
}

# 2017-01-12 rewritten yet again....
accersi_evolutionary_success <- function(record)	{
ANCESTORID <- match("ancestor.id",colnames(record))
max_anc <- max(record[,ANCESTORID])
true_ancestors <- sort(unique(record[,ANCESTORID]))[2:length(unique(record[,ANCESTORID]))]
# table(x) gives number of times a value appears in array x
#success <- table(record[,ANCESTORID])[2:length(true_ancestors)]
success <- table(record[,ANCESTORID])
# cull the zero!
return(success[as.numeric(names(success))>0])
}

convert_venn_tree_to_noded_tree <- function (venn_tree)	{
ttl_ancestors <- dim(venn_tree)[1]
notus <- dim(venn_tree)[2]
max_otus <- max(venn_tree)
noded_tree <- venn_tree
venn_richness <- as.numeric(rowSums(venn_tree>0))
# get the nodes with at least 2 but not all species.  (Skip basal one!)
relevant_nodes <- sort((1:length(venn_richness))[venn_richness>1],decreasing=TRUE)
for (rn in 1:(length(relevant_nodes)-1))	{
	ta <- relevant_nodes[rn]	# node number
	token <- as.numeric(noded_tree[ta,1])	# first species
	sibs <- sum(noded_tree[ta,]>0)	# branches coming from node ta
	sisters <- as.numeric(noded_tree[ta,(1:sibs)])
	nests <- which(noded_tree==token,arr.ind=TRUE)
	if (length(nests)==2)	{
		xx <- rownames(nests)
		nests <- matrix(0,1,2)
		nests[1,] <- which(noded_tree==token,arr.ind=TRUE)
		rownames(nests) <- xx
		}	else	{
		nests <- nests[order(nests[,1],decreasing=TRUE),]
		}
	replacement <- as.numeric(rownames(noded_tree)[ta])+max_otus		# htu number
	# 2017-01-14: redo this 
	start <- match(ta,nests[,1])+1
	for (nd in start:dim(nests)[1])	{
		lnd <- nests[nd,1]	# lower node
		anc_node_rich <- sum(noded_tree[lnd,]>0)
		# take members of clade ta out of clade lnd
		rewrite <- remove_elements_from_vector(to_remove=sisters,orig_vector=noded_tree[lnd,])
		# add back 0's for eliminated elements
		rewrite <- c(rewrite,(rep(0,sibs)))
		marker <- 1+sum(rewrite>0)
		# put in HTU number for node
		rewrite[marker] <- replacement
		noded_tree[lnd,] <- rewrite
		}
#	ta <- ta-1
	}
maxtomy <- sum(colSums(noded_tree)>0)
return(noded_tree[,1:maxtomy])
# write.table(noded_tree,"Wagner_Trad_tree.xls",row.names=TRUE,sep="\t")
}

transform_venn_tree_to_vector_tree <- function (venn_tree)	{
Nnodes <- dim(venn_tree)[1]
notus <- dim(venn_tree)[2]
max_otus <- max(venn_tree)
base <- max_otus+1
vtree <- vector(length=(max_otus+Nnodes))
otus <- sort(venn_tree[1,])
for (s in 1:notus)	{
	spc <- otus[s]
	vtree[spc] <- max_otus+sort(which(venn_tree==spc,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
vtree[base] <- -1
vtree[base+1] <- base
for (n in 3:Nnodes)	{
	htu <- max_otus+n
	lead <- venn_tree[n,1]
	vtree[htu] <- max_otus+sort(which(venn_tree[1:(n-1),]==lead,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
return(vtree)
}

accersi_branch_divergence <- function(record)	{
ANCESTORID <- match("ancestor.id",colnames(record))
ORIG <- match("orig.time",colnames(record))
EXTN <- match("ext.time",colnames(record))
FAD <- match("FAD",colnames(record))
x <- rank(record[,FAD],na.last=FALSE)
raw_taxa <- length(x)
sampled_taxa <- sum(record[,FAD]>0)
ANCESTOR <- match("ancestral",colnames(record))
if (is.na(ANCESTOR))	{
	ancestral <- vector(length=length(record[,1]))
	for (i in 2:length(record[,1]))	ancestral[record[i,2]] <- ancestral[record[i,2]]+1
	record <- cbind(record,ancestral)
	ANCESTOR <- match("ancestral",colnames(record))
	}
otus <- accersi_otu_numbers_from_paleo_tree_sampled_taxa(record)
raw_venn_tree <- accersi_rawest_venn_tree_from_paleo_tree(record)
raw_venn_richness <- as.numeric(rowSums(raw_venn_tree>0))
raw_noded_tree <- convert_venn_tree_to_noded_tree(raw_venn_tree)
raw_noded_richness <- as.numeric(rowSums(raw_noded_tree>0))
notu <- length(otus)
max_otu <- max(otus)
htus <- max_otu+as.numeric(rownames(raw_noded_tree))
atus <- c(otus,htus)
divergences <- matrix(0,length(atus),2)
colnames(divergences) <- c("origin","end")
rownames(divergences) <- atus
raw_nodes <- dim(raw_venn_tree)[1]
for (n in notu:1)	{
	spc <- otus[n]
	nest <- which(raw_venn_tree==spc,arr.ind=TRUE)
	if (length(nest)==2)	{
		xx <- rownames(nest)
		nest <- matrix(0,1,2)
		nest[1,] <- which(raw_venn_tree==spc,arr.ind=TRUE)
		rownames(nest) <- xx
		}	else	{
		nest <- nest[order(nest[,1],decreasing=TRUE),]
		}
	divergences[n,2] <- record[spc,FAD]
	depth <- sum(raw_venn_richness[nest[,1]]==1)
	if (depth==0)	{
		divergences[n,1] <- record[spc,ORIG]
		}	else	{
		anc <- as.numeric(rownames(nest)[depth])
		divergences[n,1] <- record[anc,ORIG]
		}
	}

branch_type <- c(rep("sampled",notu),rep("unsampled",length(htus)))
for (nd in raw_nodes:2)	{
	while (raw_venn_richness[nd]==1 && nd>2)	nd <- nd-1
#	dummy_tree <- raw_venn_tree
	if (nd>1)	{
		token <- as.numeric(raw_venn_tree[nd,1])
		nest <- which(raw_venn_tree==token,arr.ind=TRUE)
		nest <- nest[order(nest[,1],decreasing=TRUE),]
		# if node is redundant, then the lead species will be repeated
		#	and diversity will remain the same
		depth_1 <- sum(nest[,2]==1)
		depth_2 <- sum(raw_venn_richness[as.numeric(nest[,1])]==raw_venn_richness[nd])
		depth <- min(depth_1,depth_2)
		# rownames gives the simulated species number of the node
		node_anc <- as.numeric(rownames(nest)[depth])
		tu <- notu+nd
		# see if common ancestor is sampled.
		if (node_anc==token)	{
			ot <- match(token,otus)
			divergences[tu,] <- divergences[ot,]
			divergences[ot,1] <- divergences[ot,2]	# set ancestral ghost lineage to zero
			branch_type[tu] <- "sampled"
			}	else	{
			divergences[tu,1] <- record[node_anc,ORIG]
			sisters <- raw_noded_tree[nd,(1:raw_noded_richness[nd])]
			sister_nos <- match(sisters,atus)
			divergences[tu,2] <- max(divergences[sister_nos,1])
			}	# 2017-01-13 HERE!!!!
		nd <- nd-1
		divergences[(notu+1):length(atus),]
		}
	}
return(divergences)
}

accersi_branch_durations <- function(record)	{
ANCESTORID <- match("ancestor.id",colnames(record))
ORIG <- match("orig.time",colnames(record))
EXTN <- match("ext.time",colnames(record))
FAD <- match("FAD",colnames(record))
x <- rank(record[,FAD],na.last=FALSE)
raw_taxa <- length(x)
sampled_taxa <- sum(record[,FAD]>0)
ANCESTOR <- match("ancestral",colnames(record))
if (is.na(ANCESTOR))	{
	ancestral <- vector(length=length(record[,1]))
	for (i in 2:length(record[,1]))	ancestral[record[i,2]] <- ancestral[record[i,2]]+1
	record <- cbind(record,ancestral)
	ANCESTOR <- match("ancestral",colnames(record))
	}
otus <- accersi_otu_numbers_from_paleo_tree_sampled_taxa(record)
raw_venn_tree <- accersi_rawest_venn_tree_from_paleo_tree(record)
raw_venn_richness <- as.numeric(rowSums(raw_venn_tree>0))
raw_noded_tree <- convert_venn_tree_to_noded_tree(venn_tree=raw_venn_tree)
raw_noded_richness <- as.numeric(rowSums(raw_noded_tree>0))
notu <- length(otus)
max_otu <- max(otus)
htus <- max_otu+as.numeric(rownames(raw_noded_tree))
atus <- c(otus,htus)
divergences <- matrix(0,length(atus),2)
colnames(divergences) <- c("origin","end")
rownames(divergences) <- atus
raw_nodes <- dim(raw_venn_tree)[1]
for (n in notu:1)	{
	spc <- otus[n]
	nest <- which(raw_venn_tree==spc,arr.ind=TRUE)
	if (length(nest)==2)	{
		xx <- rownames(nest)
		nest <- matrix(0,1,2)
		nest[1,] <- which(raw_venn_tree==spc,arr.ind=TRUE)
		rownames(nest) <- xx
		}	else	{
		nest <- nest[order(nest[,1],decreasing=TRUE),]
		}
	divergences[n,2] <- record[spc,FAD]
	depth <- sum(raw_venn_richness[nest[,1]]==1)
	if (depth==0)	{
		divergences[n,1] <- record[spc,ORIG]
		}	else	{
		anc <- as.numeric(rownames(nest)[depth])
		divergences[n,1] <- record[anc,ORIG]
		}
	}

branch_type <- c(rep("sampled",notu),rep("unsampled",length(htus)))
for (nd in raw_nodes:2)	{
	while (raw_venn_richness[nd]==1 && nd>2)	nd <- nd-1
#	dummy_tree <- raw_venn_tree
	if (nd>1)	{
		token <- as.numeric(raw_venn_tree[nd,1])
		nest <- which(raw_venn_tree==token,arr.ind=TRUE)
		nest <- nest[order(nest[,1],decreasing=TRUE),]
		# if node is redundant, then the lead species will be repeated
		#	and diversity will remain the same
		depth_1 <- sum(nest[,2]==1)
		depth_2 <- sum(raw_venn_richness[as.numeric(nest[,1])]==raw_venn_richness[nd])
		depth <- min(depth_1,depth_2)
		# rownames gives the simulated species number of the node
		node_anc <- as.numeric(rownames(nest)[depth])
		tu <- notu+nd
		# see if common ancestor is sampled.
		if (node_anc==token)	{
			ot <- match(token,otus)
			divergences[tu,] <- divergences[ot,]
			divergences[ot,1] <- divergences[ot,2]	# set ancestral ghost lineage to zero
			branch_type[tu] <- "sampled_ancestor"
			}	else	{
			divergences[tu,1] <- record[node_anc,ORIG]
			sisters <- raw_noded_tree[nd,(1:raw_noded_richness[nd])]
			sister_nos <- match(sisters,atus)
			divergences[tu,2] <- max(divergences[sister_nos,1])
			}	# 2017-01-13 HERE!!!!
		nd <- nd-1
		divergences[(notu+1):length(atus),]
		}
	}
branch_durations <- abs(divergences[,1]-divergences[,2])
names(branch_durations) <- branch_type
return(branch_durations[branch_durations>0])
}

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


#### Routines from 2019 Quantitative Course ####
accersi_fossil_record_summary_from_paleotree_output <- function(fossil_record,exclude_recent=F)	{
true_history <- fossilRecord2fossilTaxa(fossil_record);	# extract origins, ancestors, durations & other information;
if (!is.data.frame(true_history))
	true_history <- data.frame(true_history,stringsAsFactors = F);
nTotalTaxa <- nrow(true_history);
sampling_history <- matrix(0,nrow(true_history),1);
sepkoski <- data.frame(FAD=as.numeric(rep(0,nrow(true_history))),LAD=as.numeric(rep(0,nrow(true_history))));
for (i in 1:length(fossil_record))	{
	if (length(fossil_record[[i]]$sampling.times)>ncol(sampling_history))	{
		dummy <- matrix(0,nTotalTaxa,length(fossil_record[[i]]$sampling.times)-ncol(sampling_history))
		sampling_history <- cbind(sampling_history,dummy);
		}
	if (length(fossil_record[[i]]$sampling.times)>0)	{
		sampling_history[i,1:length(fossil_record[[i]]$sampling.times)] <- fossil_record[[i]]$sampling.times;
		}
	if (max(sampling_history[i,])>0)	{
		sepkoski$FAD[i] <- max(sampling_history[i,]);
		if (true_history$still.alive[i]==0)	{
			sepkoski$LAD[i] <- min(sampling_history[i,sampling_history[i,]>0]);
			} else	{
			sepkoski$LAD[i] <- 0;
			}
		}
	}
output <- list(sepkoski,sampling_history);
names(output) <- c("strat_ranges","sampling_history");
return(output);
}

accersi_vector_tree_from_paleotree_output <- function(simulated_history)	{
node_ancestors <- sort(unique(simulated_history$ancestor.id));
notu <- nrow(simulated_history);
nNodes <- length(node_ancestors);
htus <- (1:nNodes)+notu;
vector_tree <- vector(length=notu+nNodes);
vector_tree[1:notu] <- notu+match(simulated_history$ancestor.id,node_ancestors)
vector_tree[(notu+1):(notu+nNodes)] <- notu+match(simulated_history$ancestor.id[node_ancestors],node_ancestors);
vector_tree[node_ancestors] <- notu+match(node_ancestors,node_ancestors);
vector_tree[notu+1] <- -1;
return(vector_tree);
}

accersi_vector_tree_and_ancestors_from_paleotree_output <- function(fossil_record)	{
# returns 2 vectors:
# vector_tree: ancestor of each simulated species
# node_ancestor: which species each node represents
simulated_history <- data.frame(fossilRecord2fossilTaxa(fossil_record),stringsAsFactors = F);
node_ancestors <- sort(unique(simulated_history$ancestor.id));
notu <- nrow(simulated_history);
nNodes <- length(node_ancestors);
htus <- (1:nNodes)+notu;
vector_tree <- vector(length=notu+nNodes);
vector_tree[1:notu] <- notu+match(simulated_history$ancestor.id,node_ancestors)
vector_tree[(notu+1):(notu+nNodes)] <- notu+match(simulated_history$ancestor.id[node_ancestors],node_ancestors);
vector_tree[node_ancestors] <- notu+match(node_ancestors,node_ancestors);
vector_tree[notu+1] <- -1;
#vector_tree[1:notu]
output <- list(vector_tree,node_ancestors);
names(output) <- c("vector_tree","node_ancestors");
return(output);
}

accersi_sampled_only_vector_tree_from_paleo_tree <- function(fossil_record_plus)	{
# basic procedure: create a Venn-Diagram tree and eliminate the unsampled taxa, then reduce it to just the sampled taxa; extract vector tree from that.

sampled_taxa <- (1:nrow(fossil_record_plus))[fossil_record_plus$FAD!=0];
sampled_richness <- length(sampled_taxa);
vector_tree <- accersi_vector_tree_from_paleotree_output(simulated_history=fossil_record_plus)
venn_tree_orig <- transform_vector_tree_to_venn_tree(vector_tree);

sampled_original_nodes <- original_ancestor <- sampled_node_richness <- node_names <- c();
venn_tree_sampled <- matrix(0,nrow=0,ncol=sampled_richness);
for (nn in 1:nrow(venn_tree_orig))	{
	fossil_reps <- venn_tree_orig[nn,venn_tree_orig[nn,] %in% sampled_taxa];
	if (length(fossil_reps)>0)	{
		sampled_original_nodes <- c(sampled_original_nodes,nn);
		original_ancestor <- c(original_ancestor,venn_tree_orig[nn,1]);
		sampled_node_richness <- c(sampled_node_richness,length(fossil_reps));
		node_names <- c(node_names,(paste(fossil_reps,collapse=",")));
		fossil_reps <- c(fossil_reps,rep(0,sampled_richness-length(fossil_reps)));
		venn_tree_sampled <- rbind(venn_tree_sampled,fossil_reps);
#		print(c(sampled_node_richness[nrow(venn_tree_sampled)],venn_tree_sampled[nrow(venn_tree_sampled),1:10]));
		}
	nn <- nn+1;
	}
ancestral_origination <- fossil_record_plus$orig.time[original_ancestor];	# divergence of the original ancestor of a node;
ancestral_extinction <- fossil_record_plus$ext.time[original_ancestor];		# extinction of the original ancestor of a node;
sampled_taxon_divergences <- orig_divergences <- fossil_record_plus$orig.time[sampled_taxa];	# divergences of sampled taxa
sampled_taxon_extinctions <- fossil_record_plus$ext.time[sampled_taxa];		# extinctions of sampled taxa

# two step reduction. First, use the node_names string to identify the same rows.  Tally the lumped ancestors as you do.
unique_node_names <- unique(node_names);
keepers <- match(unique_node_names,node_names);
srepeek <- 1+length(node_names)-match(unique_node_names,node_names[length(node_names):1]);	# match things in reverse to get latest examples
total_branches_2 <- rep(1,sampled_richness);	# give all species a branch length of 1;
total_branches_2[(1:sampled_richness)[sampled_taxa %in% original_ancestor]] <- 0;	# make ancestors 0: the node will take over.
for (un in 1:length(unique_node_names))
	total_branches_2 <- c(total_branches_2,sum(node_names==unique_node_names[un]));

venn_tree_sampled_2 <- venn_tree_sampled[keepers,];
original_ancestor_2 <- original_ancestor[keepers];
ancestral_origination_2 <- fossil_record_plus$orig.time[original_ancestor_2];
ancestral_extinction_2 <- fossil_record_plus$ext.time[original_ancestor[srepeek]];
sampled_original_nodes_2 <- sampled_original_nodes[keepers];
node_names_2 <- node_names[keepers];
node_divergences <- fossil_record_plus$orig.time[original_ancestor_2];	# set the nodes to their earliest member after compression
sampled_node_richness_2 <- sampled_node_richness[keepers];

# now, compress singletons: they get the origin time of their first "node"
compressed_ancestors_singletons <- original_ancestor_2[sampled_node_richness_2==1]	# original ancestors of nodes with only one sampled species
sampled_singletons <- match(venn_tree_sampled_2[sampled_node_richness_2==1,1],sampled_taxa)	# connect to sampled taxa;
# put divergence times of singletons with their unsampled ancestor's divergence times
sampled_taxon_divergences[sampled_singletons] <- fossil_record_plus$orig.time[compressed_ancestors_singletons];
# remove rows with just one species
venn_tree_sampled_3 <- venn_tree_sampled_2[sampled_node_richness_2>1,];
sampled_nodes <- nrow(venn_tree_sampled_3);
original_ancestor_3 <- original_ancestor_2[sampled_node_richness_2>1];
ancestral_origination_3 <- ancestral_origination_2[sampled_node_richness_2>1];
ancestral_extinction_3 <- ancestral_extinction_2[sampled_node_richness_2>1];
sampled_original_nodes_3 <- sampled_original_nodes_2[sampled_node_richness_2>1];

# use the name strings to find out how many unsampled ancestors single species have.
species_only_names <- node_names_2[sampled_node_richness_2==1];
unique_species_only_names <- unique(species_only_names);
species_only_numbers <- venn_tree_sampled_2[sampled_node_richness_2==1,1];
for (us in 1:length(unique_species_only_names))
	total_branches_2[match(species_only_numbers[us],sampled_taxa)] <- total_branches_2[match(species_only_numbers[us],sampled_taxa)]+sum(species_only_names==unique_species_only_names[us]);
# eliminate nodes containing only one sampled species
total_branches_3 <- total_branches_2[sampled_node_richness_2>1];

sampled_ancestors <- sampled_taxa[sampled_taxa %in% original_ancestor_3];							# these are the OTU numbers for the sampled tree
unobserved_nodes <- sampled_richness+(1:sampled_nodes)[!original_ancestor_3 %in% sampled_taxa];		# these are the HTU numbers for the sampled tree
observed_nodes <- sampled_richness+(1:sampled_nodes)[original_ancestor_3 %in% sampled_taxa];		# these are the HTU numbers for the sampled tree
htu_nos <- sampled_richness+(1:sampled_nodes);														# 
unobserved_htus <- htu_nos[!original_ancestor_3 %in% sampled_taxa];
observed_htus <- htu_nos[original_ancestor_3 %in% sampled_taxa];
unsampled_ancestor_orig_nos <- original_ancestor_3[!original_ancestor_3 %in% sampled_taxa];
node_original_species_nos <- original_ancestor_3;

venn_tree <- matrix(0,nrow=0,ncol=sampled_richness);
sampled_dummy <- c(0,sampled_taxa);
for (nn in 1:sampled_nodes)
	venn_tree <- rbind(venn_tree,match(venn_tree_sampled_3[nn,],sampled_dummy)-1);
vector_tree <- transform_venn_tree_to_vector_tree(venn_tree);
mxtree <- transform_vector_tree_to_matrix_tree(vector_tree);
	
all_divergence <- c(sampled_taxon_divergences,ancestral_origination_3);
all_extinction <- c(sampled_taxon_extinctions,ancestral_extinction_3);

# now, because we have budding models sometimes, let's make sure that sister species with unsampled ancestors have common divergence times.
node_nos <- sort(c(observed_htus,unobserved_htus));
for (nn in nrow(venn_tree):1)	{
	htu_number <- nn+sampled_richness;
	if (node_nos[nn] %in% unobserved_htus)	{
		f1 <- mxtree[nn,mxtree[nn,]!=0];
		all_extinction[htu_number] <- all_divergence[f1] <- max(all_divergence[f1]);
		}
	}

sampled_durations <- data.frame(origination=as.numeric(all_divergence),extinction=as.numeric(all_extinction),stringsAsFactors = F);
#cbind(node_original_species_nos,ancestral_origination_2,ancestral_extinction_2)

output <- list(vector_tree,venn_tree,sampled_durations,sampled_taxa,sampled_ancestors,match(sampled_ancestors,sampled_taxa),observed_nodes,unobserved_nodes,observed_htus,unobserved_htus,unsampled_ancestor_orig_nos,node_original_species_nos,sampled_original_nodes_3,total_branches_3);
#output <-     list(vector_tree,venn_tree,sampled_durations,sampled_taxa,sampled_ancestors,match(sampled_ancestors,sampled_taxa),unsampled_nodes,unsampled_ancestor_orig_nos,node_original_species_nos,sampled_original_nodes_3);
names(output) <- c("vector_tree","venn_tree","durations","sampled_taxa","sampled_ancestors_orig_no","sampled_ancestral_otus","observed_nodes_orig","unobserved_nodes_orig","observed_nodes_red","unobserved_nodes_red","unsampled_ancestors_orig_no","node_orig_species_no","original_node_no","total_branches");
return(output);
}


	##### Divergence Times #######
accersi_unsampled_branch_durations <- function(reduced_history,strat_ranges_red)	{
ctree <- reduced_history$vector_tree;								# vector tree giving cladogram for sampled species only
sotu <- match(-1,ctree)-1;
mtree <- transform_vector_tree_to_matrix_tree(ctree);
ancestral_otus <- reduced_history$sampled_ancestral_otus;			# original taxon numbers of sampled ancestors
observed_nodes_red <- reduced_history$observed_nodes_red;
unobserved_nodes_red <- reduced_history$unobserved_nodes_red;
sampled_ancestors_red <- c(ancestral_otus,rep(0,length(unobserved_nodes_red)));	# sampled otu numbers of nods
sampled_ancestors_red <- sampled_ancestors_red[order(c(observed_nodes_red,unobserved_nodes_red))];
durations_red <- -abs(reduced_history$durations);					# durations including those of unsampled ancestors
#min_node_age <- vector(length=nrow(mtree));
#dummy <- data.frame(FAD=as.numeric(rep(0,length(ctree)-sotu)),LAD=as.numeric(rep(0,length(ctree)-sotu)))
dummy <- matrix(0,nrow=length(ctree)-sotu,ncol=2);
colnames(dummy) <- colnames(strat_ranges_red);
strat_ranges_red <- rbind(strat_ranges_red,dummy);
#strat_ranges_red[observed_nodes_red,] <- strat_ranges_red[ancestral_otus,];
range_extension <- signor_lipps <- jaanusson <- vector(length=length(ctree));
nNodes <- nrow(mtree);
for (nd in nNodes:1)	{
	f1 <- mtree[nd,!mtree[nd,] %in% c(0,ancestral_otus)];
	jaanusson[f1] <- strat_ranges_red[f1,1]-durations_red[f1,1];
	if (sampled_ancestors_red[nd]!=0)	{
		sll <- max(durations_red[f1,1]-strat_ranges_red[sampled_ancestors_red[nd],2]);
		signor_lipps[sampled_ancestors_red[nd]] <- max(c(0,sll));
		strat_ranges_red[nd+sotu,] <- strat_ranges_red[sampled_ancestors_red[nd],];
		}
	}
output <- list(jaanusson,signor_lipps);
names(output) <- c("jaanusson_lineages","signor_lipps_lineages");
return(output);
}
