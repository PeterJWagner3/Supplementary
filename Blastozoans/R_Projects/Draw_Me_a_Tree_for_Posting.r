#### Setup Program ####
this_directory <- "~/Documents/R_Projects/Tree_Drawing/";
setwd(this_directory);
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
franky <- "Franklin Gothic Medium";  #Or whatever other R-compatible font you'd like to use...

#install.packages("png")
library(png);
library(zoo);
source(paste(common_source_folder,"General_Plot_Templates.r",sep=""));
source(paste(common_source_folder,"Nexus_File_Routines.r",sep=""));

data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";
load(paste(data_for_R_folder,"Gradstein_2020_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations

# set parameters & files for printing tree ####
tree_to_print <- "1A";

if (tree_to_print %in% c("1A","1B"))	{
	strat_rank <- "Epoch";		# I recommend either "Stage" or "Stage Slice"
	rank_to_plot <- "order";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Diploporita_Fuzzy_Ranges.xlsx",sep="");
	otu_information_file <- paste(this_directory,"taxon_info/Diploporita_Information.xlsx",sep="");
	if (tree_to_print=="1A")	{
		mctree_file <- paste(this_directory,"trees/dip_mcc.tre",sep="");
		} else if (tree_to_print=="1B")	{
		mctree_file <- paste(this_directory,"trees/diploploranNoOutGroupNoJump.mcc.tre",sep="");
		}
	} else if (tree_to_print %in% c("2A","2B"))	{
	strat_rank <- "Epoch";		# I recommend either "Stage" or "Stage Slice"
	rank_to_plot <- "order";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Blastoidea_Fuzzy_Ranges.xlsx",sep="");
	if (tree_to_print=="2A")	{
		mctree_file <- paste(this_directory,"trees/blastoid_constrained.mcc.tre",sep="");
		} else if (tree_to_print=="2B")	{
		}
	} else if (tree_to_print %in% c("3A","3B"))	{
	if (tree_to_print=="3A")	{
		} else if (tree_to_print=="3B")	{
		}
	strat_rank <- "Stage Slice";		# I recommend either "Stage" or "Stage Slice"
	rank_to_plot <- "family";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Paracrinoidea_Fuzzy_Ranges.xlsx",sep="");
	otu_information_file <- paste(this_directory,"taxon_info/Paracrinoidea_Information.xlsx",sep="");
	}

# load stratigraphic data & set it up for figures ####
if (gsub(".csv","",otu_information_file)!=otu_information_file)	{
	otu_information <- read.csv(otu_information_file,header=T,stringsAsFactors = F);
	} else if (gsub(".xlsx","",otu_information_file)!=otu_information_file)	{
	otu_information <- as.data.frame(readxl::read_xlsx(otu_information_file));
	}
colnames(otu_information)[colnames(otu_information)=="species"] <- "taxon";
colnames(otu_information) <- tolower(colnames(otu_information));

if (gsub(".csv","",range_file)!=range_file)	{
	fuzzy_ranges <- read.csv(range_file,header=T,stringsAsFactors = F);
	} else if (gsub(".xlsx","",range_file)!=range_file)	{
	fuzzy_ranges <- as.data.frame(readxl::read_xlsx(range_file));
	}
colnames(fuzzy_ranges)[colnames(fuzzy_ranges)=="species"] <- "taxon";
colnames(fuzzy_ranges) <- tolower(colnames(fuzzy_ranges));
fuzzy_ranges$taxon <-gsub("_"," ",fuzzy_ranges$taxon);
#fuzzy_ranges$taxon <- gsub("Crassatellites","Crassatella",fuzzy_ranges$taxon);
faux_recent <- min(fuzzy_ranges$fa_ub);
oldest_fa <- -(max(fuzzy_ranges$fa_ub)+max(fuzzy_ranges$fa_lb))/2;
rescaled_fa_bounds <- fuzzy_ranges;
rescaled_fa_bounds$la_lb <- rescaled_fa_bounds$la_ub <- NULL;
rescaled_fa_bounds$fa_lb <- rescaled_fa_bounds$fa_lb - faux_recent;
rescaled_fa_bounds$fa_ub <- rescaled_fa_bounds$fa_ub - faux_recent;

finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$chronostratigraphic_rank==strat_rank)
finest_chronostrat <- subset(finest_chronostrat,finest_chronostrat$scale=="International")
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];

stage_slice_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
stage_slice_chronostrat <- stage_slice_chronostrat[order(-abs(stage_slice_chronostrat$ma_lb)),];

if (strat_rank=="Stage Slice")	finest_chronostrat <- stage_slice_chronostrat;

#### Get Tree Information for Plotting ####
mctree <- scan(file=mctree_file,what=character(),sep="\n");
mctree <- gsub("\t","",mctree)
begin_taxa <- match("begin taxa;",tolower(mctree));
notu <- as.numeric(strsplit(gsub(";","",mctree[begin_taxa+1]),split = "=")[[1]][2]);

taxlabels <- match("taxlabels",tolower(mctree));

otu_names <- fuzzy_ranges$taxon[fuzzy_ranges$taxon %in% gsub("_"," ",mctree[(taxlabels+1):(taxlabels+notu)])];
otu_names_nex <- gsub(" ","_",otu_names);
tree_line <- 1+match("begin trees;",tolower(mctree));
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\]","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&U\\]","",newick_string_full);

tree_output <- read_newick_string_mcmc(newick_string_full,otu_names);
clade_posteriors <- tree_output$clade_posteriors;
sampled_ancestors <- tree_output$ancestral;

basal_taxon <- fuzzy_ranges$taxon[match(max(fuzzy_ranges$fa_lb),fuzzy_ranges$fa_lb)];
basal_date_raw <- max(((tree_output$hpd$lb+tree_output$hpd$ub)/2)[1:notu]);	# base from Bayesian tree before rescaling 0 to actual upper bound
basal_taxon_no <- match(basal_date_raw,((tree_output$hpd$lb+tree_output$hpd$ub)/2)[1:notu]);

date_rescale <- min(fuzzy_ranges$fa_ub);
fuzzy_ranges$hpd_lb <- tree_output$hpd_age$lb[1:notu];
fuzzy_ranges$hpd_ub <- tree_output$hpd_age$ub[1:notu];
#write.csv(fuzzy_ranges,"WTF.csv",row.names = FALSE);
hpd <- -(tree_output$hpd+date_rescale);

#hpd <- data.frame(-min(fuzzy_ranges$fa_lb)-tree_output$hpd);
hpd$md <- (hpd$lb+hpd$ub)/2;
branch_durations <- abs(tree_output$branch_durations);
#branch_durations[(1:notu)[sampled_ancestors==1]] <- 0;
vector_tree <- tree_output$vector_tree;
notu <- match(-1,tree_output$vector_tree)-1;
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree=tree_output$vector_tree);
nNodes <- nrow(mat_tree);
lowest_taxon <- min(which(mat_tree %in% (1:notu),arr.ind = T));
outgp <- mat_tree[lowest_taxon,mat_tree[lowest_taxon,]<notu];
#-max(abs(hpd$md[outgp]))

### Construct proper phylogeny ####
root_ma <- min(hpd$md);
branch_ranges <- data.frame(start=as.numeric(vector_tree),
														finish=as.numeric(vector_tree));
rownames(branch_ranges) <- names(branch_durations)[!is.na(names(branch_durations))];
branch_ranges$finish <- hpd$md;

branch_ranges[notu+1,] <- root_ma;
branch_ranges$start[1:notu] <- hpd$md[vector_tree[1:notu]];
branch_ranges$start[(notu+2):(notu+nNodes)] <-hpd$md[vector_tree[(notu+2):(notu+nNodes)]];

ancestral_spc <- (1:notu)[sampled_ancestors==1];
names(ancestral_spc) <- names(sampled_ancestors)[sampled_ancestors==1];
anagenetic_ancestors <- rep(0,notu);
names(anagenetic_ancestors) <- names(sampled_ancestors);
if (length(ancestral_spc)>0)	{
	pbdb_finds <- pbdb_data_list$pbdb_finds;
	pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
	venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
	mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
#	cinctan_finds <- read.csv(paste(this_directory,"Cincta_Finds.csv",sep=""),header=T,stringsAsFactors = F);
	for (an in 1:length(ancestral_spc))	{
		anc <- ancestral_spc[an];
		daughter_node <- max(which(venn_tree==anc,arr.ind = T)[,1]);
		its_node <- which(mat_tree==anc,arr.ind = T)[,1]
		desc <- mat_tree[its_node,mat_tree[its_node,]!=anc];
	#	hpd[anc,]
	#	hpd[desc,]
	#	fuzzy_ranges[anc,]
	#	fuzzy_ranges[desc,]
		# prior routine using venn_tree
		progeny <- venn_tree[daughter_node,venn_tree[daughter_node,]>0];
		progeny <- progeny[progeny!=anc];
		anc_find_data <- subset(pbdb_finds,pbdb_finds$accepted_name==otu_names[anc])
		if(sum(abs(fuzzy_ranges$la_lb[anc])<=abs(fuzzy_ranges$fa_ub[progeny]))==0)	{
#			find_data <- subset(pbdb_finds,pbdb_finds$accepted_name %in% otu_names[c(anc,progeny)]);
			accepted_names <- unique(pbdb_taxonomy$accepted_name[pbdb_taxonomy$taxon_name %in% otu_names[c(anc,progeny)]]);
			accepted_names <- accepted_names[!is.na(accepted_names)];
			poss_names <- pbdb_taxonomy$taxon_name[pbdb_taxonomy$accepted_name %in% accepted_names];
			poss_taxon_nos <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$accepted_name %in% accepted_names];
			find_data <- pbdb_finds[pbdb_finds$accepted_no %in% poss_taxon_nos,]
			find_data$accepted_name <- pbdb_taxonomy$accepted_name[match(find_data$accepted_no,pbdb_taxonomy$taxon_no)];
#			pbdb_finds[pbdb_finds$identified_name=="Crassatella tumidula",]
			coccur_mat <- accersi_cooccurence_matrix(find_data);
			if (sum(coccur_mat[match(otu_names[anc],rownames(coccur_mat)),])==1)
#			if (coccur_mat[1,2]==0)
				anagenetic_ancestors[anc] <- 1;
			}
		}
	}

#nNodes <- nrow(mat_tree)
phylo_axis <- get_phylo_axis_from_newick_string_w_anagenesis(newick_string = tree_output$newick,sampled_ancestors=sampled_ancestors,anagenetic_ancestors=anagenetic_ancestors);
write.csv(cbind(sampled_ancestors,anagenetic_ancestors),"Ancestors2.csv",row.names = TRUE)
names(phylo_axis) <- c((otu_names),paste("nd_",1:nNodes,sep=""));


### use this to identify which sequences get lumped ####
# get anagenetic series of species
anagenetic_anc_no <- (1:notu)[anagenetic_ancestors==1];
names(anagenetic_anc_no) <- names(anagenetic_ancestors)[anagenetic_anc_no];
axis_ranks <- rank(phylo_axis[1:notu]);
axis_ranks <- axis_ranks-min(axis_ranks)+1;
anagenetic_overlaps <- hist(axis_ranks[1:notu],breaks=(min(axis_ranks[1:notu])-1):ceiling(max(axis_ranks[1:notu])),plot=F)$counts;
#anagenetic_overlaps <- hist(phylo_axis[1:notu],breaks=((min(phylo_axis[1:notu])-1):max(phylo_axis[1:notu])),plot=F)$counts;
anagenetic_segments <- names(anagenetic_ancestors[anagenetic_ancestors==1]);
#anagenetic_segments <- otu_names[anagenetic_ancestors==1];
anagenetic_segments_nex <- gsub(" ","_",anagenetic_segments);
#anagenetic_segments_nex <- otu_names_nex[anagenetic_ancestors==1];
names(phylo_axis) <- all_names <- rownames(branch_ranges);
names(anagenetic_overlaps) <- sort(unique(phylo_axis[1:notu]));
anagenetic_overlaps <- anagenetic_overlaps[anagenetic_overlaps>1];
anagenetic_fuzzies <- data.frame(ancestor=as.character(),descendant=as.character(),
								 fa_lb=as.numeric(),fa_ub=as.numeric(),la_lb=as.numeric(),la_ub=as.numeric(),
								 anc_no=as.numeric(),dsc_no=as.numeric(),
								 phylo_axis=as.numeric(),taxon=as.character(),stringsAsFactors = F);
an <- 0;
while (an<length(anagenetic_overlaps))	{
	an <- an+1;
	dummy <- data.frame(ancestor=as.character("a"),descendant=as.character("z"),
						fa_lb=as.numeric(0),fa_ub=as.numeric(0),la_lb=as.numeric(0),la_ub=as.numeric(0),
						anc_no=as.numeric(0),dsc_no=as.numeric(0),phylo_axis=as.numeric(0),taxon=as.character("tx"),stringsAsFactors = F);
	anagenetic_fuzzies <- rbind(anagenetic_fuzzies,dummy);
	morphospc <- (1:notu)[phylo_axis[1:notu]==as.numeric(names(anagenetic_overlaps)[an])];
#	if (length(unique(fuzzy_ranges$fa_lb[morphospc]))==1 && length(unique(fuzzy_ranges$fa_ub[morphospc])) && length(unique(fuzzy_ranges$la_lb[morphospc])))	{
	anc_no <- morphospc[morphospc %in% anagenetic_anc_no];
	dsc_no <- morphospc[!morphospc %in% anagenetic_anc_no];
	#	all_names[morphospc][!all_names[morphospc] %in% anagenetic_segments_nex]
	anagenetic_fuzzies$ancestor[an] <- all_names[anc_no];
	anagenetic_fuzzies$descendant[an] <- all_names[dsc_no];
	anagenetic_fuzzies$fa_lb[an] <- fuzzy_ranges$fa_lb[anc_no];
	anagenetic_fuzzies$fa_ub[an] <- max(fuzzy_ranges$fa_ub[morphospc]);
	anagenetic_fuzzies$la_lb[an] <- min(fuzzy_ranges$la_lb[morphospc]);
	anagenetic_fuzzies$la_ub[an] <- fuzzy_ranges$la_ub[dsc_no];
	anagenetic_fuzzies$anc_no[an] <- anc_no;
	anagenetic_fuzzies$dsc_no[an] <- dsc_no;
	anagenetic_fuzzies$taxon[an] <- paste(otu_names[anc_no],"->",otu_names[dsc_no],sep="");
	anagenetic_fuzzies$phylo_axis[an] <- phylo_axis[anc_no];
#		}
	}

#### plot time scale ####
st <- sum(max(abs(branch_ranges$start)) <= abs(finest_chronostrat$ma_lb));	# starts!
en <- sum(min(abs(fuzzy_ranges$la_ub)) <= abs(finest_chronostrat$ma_ub));	# endsd!

standard_time_scale <- finest_chronostrat[st:en,];
standard_time_scale$ma_lb <- -abs(standard_time_scale$ma_lb);
standard_time_scale$ma_ub <- -abs(standard_time_scale$ma_ub);
standard_time_scale <- standard_time_scale[standard_time_scale$interval_sr=="",]

# Add names/symbols for plotting
strat_names <- standard_time_scale$st;
# Get the colors for time units
strat_colors <- standard_time_scale$color;
# Set the oldest and youngest intervals by names on whatever chronostratigraphic scale you used in the analysis
oldest_interval <- standard_time_scale$interval[1];
youngest_interval <- standard_time_scale$interval[nrow(standard_time_scale)];
# get the time scale that will be plotted
time_scale_to_plot <- unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub));
# here, I want to adjust the earliest date: I don't want the whole Cambrian
time_scale_to_plot[1] <- -abs(stage_slice_chronostrat$ma_lb[sum(-abs(stage_slice_chronostrat$ma_lb) < (-ceiling(max(abs(branch_ranges$start)))))]);

# Now, set the onset and end of the x-axis
onset <- min(time_scale_to_plot);
end <- max(time_scale_to_plot);
# Finally, set up breaks for x-axis
#yearbreaks <- c(5,25,50);					# set breaks for x-axis (minor, medium & major)
yearbreaks <- as.numeric(set_axis_breaks(max_no=end,min_no=onset));
# now, set up the y-axis: this will reflect your data

use_strat_labels <- T;						# if T, then strat_names will be plotted on X-axis inside boxes
alt_back <- F;								# if T, then the background will alternat shades between major intervals
plot_title <- "";							# Name of the plot; enter "" for nothing
ordinate <- "";								# Label of Y-axis
hues <- "T";								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
ysize <- 6*4.285714285/6;
ysize <- 5;
xsize <- 5; 
mxy <- max(phylo_axis)
mny <- min(phylo_axis);									# set maximum y value
#mny <- 0;
yearbreaks <- sort(as.numeric(set_axis_breaks_new(max(abs(time_scale_to_plot))-min(abs(time_scale_to_plot)))));
#yearbreaks <- c(1,0.5,0.1);
yearbreaks <- sort(yearbreaks);

end_extra <- abs(end-onset)/4;
ysize <-  min(6,4.285714285*notu/20);
#ysize <- 5;
xsize <- 1.25*5; 
mxy <- max(phylo_axis)
mny <- min(phylo_axis);									# set maximum y value
ysize <- xsize/1.25;
print_names <- TRUE;
print_pix <- FALSE;
strat_label_size <- 2/3;					# size of the labels for chronostratigraphic units
myr_size <- 1;
Phanerozoic_Timescale_Plot_Flexible(onset,end+end_extra,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma                   ",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky,xlab_size=1.5,myr_size=myr_size);
#print_names <- T;
#print_pix <- F;
#strat_label_size <- 2;					# size of the labels for chronostratigraphic units
#Phanerozoic_Timescale_Plot_Flexible(onset,end+20,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky);
#### Plot Tree ####
#Phanerozoic_Timescale_Plot_Flexible(onset,end+end_extra,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma                   ",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky,xlab_size=1.5,myr_size=myr_size);
otu_colors <- otu_information$color;
names(otu_colors) <- otu_information$taxon;
# plot stratigraphic ranges
for (sp in 1:notu)	{
#	sp <- sp+1;
	if (!sp %in% anagenetic_fuzzies$anc_no & !sp %in% anagenetic_fuzzies$dsc_no)	{
		# upper limit exceeds lower limit, but no definite range
		if (fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp] && fuzzy_ranges$fa_ub[sp]>fuzzy_ranges$la_ub[sp])	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else if (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp] && fuzzy_ranges$fa_lb[sp]>fuzzy_ranges$la_lb[sp])	{
		# lower limit exceeds exceeds limit, but no definite range
			rect(-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			}
		if (abs(fuzzy_ranges$fa_ub[sp])>abs(fuzzy_ranges$la_lb[sp]))
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp],lwd=0);
#		if ((fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp]) || (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp]))	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			} else	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp]);
#			}
		} else if (sp %in% anagenetic_fuzzies$dsc_no)	{
		fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
		rect(-abs(anagenetic_fuzzies$fa_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
		if (abs(anagenetic_fuzzies$fa_ub[fuzzy_dude])>abs(anagenetic_fuzzies$la_lb[fuzzy_dude]))
			rect(-abs(anagenetic_fuzzies$fa_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=otu_colors[sp],lwd=0);
#		if ()
		}
	}

observed_nodes_all <- rep(0,nNodes);
an <- 0;
# draw out ancestors specially
while (an < length(ancestral_spc))	{
	an <- an+1;
	observed_nodes_all[which(mat_tree==ancestral_spc[an],arr.ind = T)[1]] <- ancestral_spc[an];
	}
# draw branches linking taxa to nodes and nodes to deeper nodes
for (nn in 1:length(phylo_axis))	{
#	nn <- nn+1;
	if (nn<=notu)	{
		segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=1,lty=2)
		} else if (nn>(notu+1))	{
		nd <- nn-notu;
		if (observed_nodes_all[nd]==0)	{
			segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=5*clade_posteriors[nd]);
			} else	{
			segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=5*clade_posteriors[nd],lty=3);
			}
		} # end case for clade
	}
# link nodes & branches
for (nd in 1:nNodes)	{
	htu <- nd+notu;
	f1 <- mat_tree[nd,mat_tree[nd,]>0];
	segments(branch_ranges$finish[htu],min(phylo_axis[f1]),branch_ranges$finish[htu],max(phylo_axis[f1]));
	}
# draw reconstructed FAs
for (sp in 1:notu)	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp],pch=8,cex=1.0);
top <- mxy+1;
taxon_font_size <- 0.6;
if (print_names)	{
	sp <- 0;
#	for (sp in 1:notu)	{
	while (sp < notu)	{
		sp <- sp+1;
		if (!sp %in% anagenetic_fuzzies$anc_no && !sp %in% anagenetic_fuzzies$dsc_no)	{
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],otu_names[sp],pos=4,cex=taxon_font_size,family=franky);
			} else	if (sp %in% anagenetic_fuzzies$dsc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
#		lineage_name <- paste(anagenetic_fuzzies$ancestor[fuzzy_dude],"->",anagenetic_fuzzies$descendant[fuzzy_dude],sep="")
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],anagenetic_fuzzies$taxon[fuzzy_dude],pos=4,cex=taxon_font_size,family=franky);
			}
		}
	}
