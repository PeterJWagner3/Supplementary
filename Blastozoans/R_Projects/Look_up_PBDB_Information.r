#### Setup Program ####
# Setup Your Analysis
pbdb_directory <- "~/Documents/R_Projects/PaleoDB_for_RevBayes/";	# directory to folder where you keep common source
setwd(pbdb_directory);

#devtools::install_github("kassambara/r2excel");
library(readxl);
source('~/Documents/R_Projects/Common_R_Source_Files/Chronos.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Historical_Diversity_Metrics.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Nexus_File_Routines.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Occurrence_Data_Routines.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Sampling_and_Occupancy_Distributions.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Stratigraphy.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_Stats_and_Probability_101.r'); 		#

franky <- "Franklin Gothic Medium";
load("~/Documents/R_Projects/Data_for_R/Gradstein_2020_Augmented.RData"); # refined Gradstein 2020 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/Rock_Unit_Database.RData"); # Rock Unit information
load("~/Documents/R_Projects/Data_for_R/PaleoDB_Edits.RData"); # Edits of PBDB Data
load("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData"); # PBDB Data
load("~/Documents/R_Projects/Data_for_R/PBDB_References.RData"); # PBDB Data
load("~/Documents/R_Projects/Data_for_R/Sepkoski_Genus_Compendium.RData"); # Jack's Compendium

temporal_precision <- 0.05;
time_scale <- gradstein_2020_emended$time_scale;
time_scale$ma_lb <- round(time_scale$ma_lb,2);
time_scale$ma_ub <- round(time_scale$ma_ub,2);
time_scale$interval_sr[is.na(time_scale$interval_sr)] <- "";

all_stage_scale <- time_scale[time_scale$chronostratigraphic_rank=="Stage" & time_scale$scale=="International",];
all_stage_scale <- all_stage_scale[order(-all_stage_scale$ma_lb),];
all_stage_scale <- all_stage_scale[all_stage_scale$interval_sr=="",];

stage_slice_scale <- time_scale[time_scale$scale=="Stage Slice",]
stage_slice_scale <- stage_slice_scale[order(-stage_slice_scale$ma_lb),];
stage_slice_scale_2 <- adjacent_stage_slice_lumping(stage_slice_scale);

rock_database <- rock_unit_data$rock_unit_database;
pbdb_sites <- pbdb_data_list$pbdb_sites_refined;
pbdb_finds <- pbdb_data_list$pbdb_finds;
#pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Trematocystis magniporatus"],]
nfinds <- nrow(pbdb_finds);
pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_opinions <- pbdb_data_list$pbdb_opinions;
paleodb_taxonomy_edits <- paleodb_fixes$paleodb_taxonomy_edits;
pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% paleodb_taxonomy_edits$taxon_no,] <- paleodb_taxonomy_edits[paleodb_taxonomy_edits$taxon_no %in% pbdb_taxonomy$taxon_no,];

# Setup Specifics for This Summary ####
#taxon_data <- accersi_data_from_nexus_file("Diploporita_Sheffield_&_Sumrall_2019.nex");
#analysis_name <- "Diploporita";
#taxon_data <- accersi_data_from_nexus_file("Blastoidea_Bauer_2018.nex");
#analysis_name <- "Blastoidea";
#taxon_data <- accersi_data_from_nexus_file("Paracrinoidea_Limbeck_et_al_2023.nex");
analysis_name <- "Proetoidea";
taxon_data_file <- "Desmoceratidae.xlsx";
taxon_data_file <- "/Users/peterjwagner/Documents/R_Projects/Devonian_Trilobites/Proetida_final.nex";
if (gsub(".xlsx","",taxon_data_file)!=taxon_data_file)	{
	taxon_data <- as.data.frame(read_xlsx(taxon_data_file))
	otu_names <- taxon_data$species;
	} else if (gsub(".csv","",taxon_data_file)!=taxon_data_file)	{
	taxon_data <- read.csv(taxon_data_file,header=T);
	otu_names <- taxon_data$species;
	} else if (gsub(".nex","",taxon_data_file)!=taxon_data_file)	{
	taxon_data <- accersi_data_from_nexus_file(taxon_data_file);
	otu_names <- taxon_data$OTUs;
	}

otu_names <- gsub("_"," ",otu_names);
otu_names <- sapply(otu_names,mundify_taxon_names);
otu_names <- unique(otu_names);
notu <- length(otu_names);

otus_in_pbdb <- match(otu_names,pbdb_taxonomy$taxon_name);
otus_in_pbdb <- otus_in_pbdb[!is.na(otus_in_pbdb)];
otu_orders <- unique(pbdb_taxonomy$order[otus_in_pbdb]);
otu_orders <- otu_orders[otu_orders!=""];
norders <- length(otu_orders); 

otu_taxonomy <- pbdb_taxonomy[pbdb_taxonomy$order %in% otu_orders,];
otu_taxonomy <- add_higher_taxon_to_pbdb_taxonomy("subfamily",otu_taxonomy);
otu_taxonomy <- add_higher_taxon_to_pbdb_taxonomy("superfamily",otu_taxonomy);

output_info <- data.frame(taxon=as.character(otu_names),accepted_name=as.character(otu_names),
						  max_age=as.numeric(rep(0,notu)),min_age=as.numeric(rep(0,notu)),
						  FA_ub=as.numeric(rep(0,notu)),FA_lb=as.numeric(rep(0,notu)),
						  LA_ub=as.numeric(rep(0,notu)),LA_lb=as.numeric(rep(0,notu)),
						  taxon_attr=as.character(rep("",notu)),
						  subfamily=as.character(rep("",notu)),family=as.character(rep("",notu)),superfamily=as.character(rep("",notu)));
species_collection_info <- data.frame(taxon=as.character(),rock_unit_senior=as.character(),
									  zone=as.character(),ma_lb=as.numeric(),ma_ub=as.numeric());
for (i in 1:notu)	{
	nn <- match(otu_names[i],otu_taxonomy$taxon_name);
	taxon_nos <- otu_taxonomy$taxon_no[otu_taxonomy$accepted_no %in% otu_taxonomy$accepted_no[nn]];
	if (is.na(nn))	{
		otu_genus <- divido_genus_names_from_species_names(otu_names[i]);
		if (!is.subgenus(otu_genus))	{
			try_this <- paste(divido_genus_names_from_species_names(otu_names[i])," (",divido_genus_names_from_species_names(otu_names[i]),") ",divido_species_epithets(otu_names[i]),sep="");
			nn <- match(try_this,otu_taxonomy$taxon_name);
			} else	{
			otu_epithet <- divido_species_epithets(otu_names[i]);
			otu_subgenera <- divido_subgenus_names_from_genus_names(otu_genus);
			nn <- match(paste(otu_subgenera[2],otu_epithet),otu_taxonomy$taxon_name);
			if (is.na(nn))
				nn <- match(paste(otu_subgenera[1],otu_epithet),otu_taxonomy$taxon_name);
			}
		}
	colls <- c();
	if (!is.na(nn))	{
		output_info$accepted_name[i] <- otu_taxonomy$accepted_name[nn];
		output_info$taxon_attr[i] <- otu_taxonomy$taxon_attr[nn];
		output_info$subfamily[i] <- otu_taxonomy$subfamily[nn];
		output_info$family[i] <- otu_taxonomy$family[nn];
		output_info$superfamily[i] <- otu_taxonomy$superfamily[nn];
		colls <- pbdb_finds$collection_no[(1:nfinds)[pbdb_finds$accepted_no %in% taxon_nos]];
		if (length(colls)==0)
			colls <- pbdb_finds$collection_no[(1:nfinds)[pbdb_finds$accepted_name_orig %in% otu_taxonomy$accepted_name[nn]]];
		} else	{
		colls <- pbdb_finds$collection_no[(1:nfinds)[pbdb_finds$identified_name %in% otu_names[i]]];
		}
	if (length(colls)>0)	{
		output_info$FA_lb[i] <- max(pbdb_sites$ma_lb[pbdb_sites$collection_no %in% colls]);
		output_info$FA_ub[i] <- max(pbdb_sites$ma_ub[pbdb_sites$collection_no %in% colls]);
		output_info$LA_lb[i] <- min(pbdb_sites$ma_lb[pbdb_sites$collection_no %in% colls]);
		output_info$LA_ub[i] <- min(pbdb_sites$ma_ub[pbdb_sites$collection_no %in% colls]);
		taxon_sites <- pbdb_sites[pbdb_sites$collection_no %in% colls,c("rock_unit_senior","zone","ma_lb","ma_ub")];
#		write.csv(pbdb_sites[pbdb_sites$collection_no %in% colls,],"pomum.csv",row.names = F)
		taxon_sites <- cbind(taxon=as.character(rep(otu_names[i],nrow(taxon_sites))),taxon_sites);
		taxon_sites <- taxon_sites[order(-taxon_sites$ma_lb,-taxon_sites$ma_ub),];
		species_collection_info <- rbind(species_collection_info,unique(taxon_sites));
		}
	}
output_file_name3 <- paste(analysis_name,"_Finds.csv",sep="");
write.csv(species_collection_info,output_file_name3,row.names = F);

output_file_name1 <- paste(analysis_name,"_Info.csv",sep="");
output_info$max_age <- (floor(output_info$FA_lb*10)/10)-(floor(min(output_info$FA_ub)*10)/10);
output_info$min_age <- (floor(output_info$FA_ub*10)/10)-(floor(min(output_info$FA_ub)*10)/10);
nexus_taxa <- output_info$taxon;
output_info$taxon <- gsub("_"," ",output_info$taxon);
write.csv(output_info,output_file_name1,row.names = F,fileEncoding = "UTF-8");
writexl::write_xlsx(output_info,gsub("csv","xlsx",output_file_name1));
output_info$taxon <- nexus_taxa;
#pbdb_finds[pbdb_finds$collection_no==138713,]
# More Stuff ####
#write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name=="Desmoceras (Pseudoughligella) bearskinense",],"bearskinense.csv",row.names = F,fileEncoding = "UTF-8");
#pbdb_sites[pbdb_sites$collection_no==86497,]

relv_periods <- time_scale[time_scale$chronostratigraphic_rank=="Period" & time_scale$scale=="International",];
relv_periods <- relv_periods[relv_periods$ma_lb>=min(output_info$FA_ub) & relv_periods$ma_ub<=max(output_info$FA_lb+5),];
relv_stages <- all_stage_scale[all_stage_scale$ma_lb>=min(output_info$FA_ub) & all_stage_scale$ma_ub<=max(output_info$FA_lb+5),];
relv_stages <- relv_stages[order(-relv_stages$ma_lb),];
relv_slices <- stage_slice_scale_2[stage_slice_scale_2$ma_lb>=min(output_info$FA_ub) & stage_slice_scale_2$ma_ub<=max(output_info$FA_lb+5),];
range1 <- paste(round(sort(relv_stages$ma_lb-(floor(min(output_info$FA_ub)*10)/10)),2),collapse=",");
range2 <- paste(round(sort(relv_slices$ma_lb-(floor(min(output_info$FA_ub)*10)/10)),2),collapse=",");
paste("timeline2 <- v(",range1,");",sep="");
paste("timeline3 <- v(",range2,");",sep="");
range3 <- paste(round(sort(relv_periods$ma_lb-(floor(min(output_info$FA_ub)*10)/10)),2),collapse=",");
paste("timeline4 <- v(",range3,");",sep="");

write.table(output_info,output_file_name2,sep="\t",row.names = F,quote=F)
pbdb_opinions[pbdb_opinions$taxon_name=="Glyptosphaerites leuchtenbergi",]
pbdb_taxonomy[pbdb_taxonomy$orig_no==381814,];
pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Macurdablastus uniplicatus"],]
pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Sphaeronites pomum"],]
#c(max(pbdb_sites$ma_lb[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Sphaeronites pomum"]]),max(pbdb_sites$ma_ub[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Sphaeronites pomum"]]))
pbdb_sites[pbdb_sites$collection_no==87101,]
write.csv(pbdb_taxonomy[pbdb_taxonomy$accepted_name %in% c("Cribroblastus melonoides","Decemoblastus melonoides"),],"Decemoblastus_melonoides.csv",row.names = F);
i <- match("Holocystites cylindricus",analysis_taxa$Taxon)
if (length(ingroup)<3)	{
	ingroup_title <- paste(ingroup,collapse="_&_");
	} else	{
	ingroup_title <- paste(ingroup[1],"_et_al",sep="");
	}
ingroup_rank <- unique(pbdb_taxonomy$taxon_rank[match(ingroup,pbdb_taxonomy$taxon_name)]);
total_group_rank <- pbdb_taxonomy$taxon_rank[match(total_group,pbdb_taxonomy$taxon_name)];
if (!total_group_rank %in% colnames(pbdb_finds))	{
	total_group_daughters <- accersi_daughter_taxa(total_group,pbdb_taxonomy);
	}	else	{
	total_group_daughters <- data.frame(daughter=total_group,daughter_rank=total_group_rank);
	}
if (!ingroup_rank %in% colnames(pbdb_finds))	{
	ingroup_daughters <- accersi_daughter_taxa(total_group,pbdb_taxonomy);
	}	else	{
	ingroup_daughters <- data.frame(daughter=ingroup,daughter_rank=ingroup_rank);
	}

pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name_orig=="Pachythaerus vindinnensis"],];
pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name_orig=="Pachythaerus similoides"],];

study_stage_scale <- time_scale[time_scale$interval %in% study_stages,];
study_stage_scale <- study_stage_scale[order(-study_stage_scale$ma_lb),];
interval_lb <- study_stages[1];
interval_ub <- study_stages[length(study_stages)];
myr_lb <- time_scale$ma_lb[match(interval_lb,time_scale$interval)];
myr_ub <- time_scale$ma_ub[match(interval_ub,time_scale$interval)];

biogeography_key <- as.data.frame(readxl::read_xlsx(biogeography_key_file));
biogeography_key$region[is.na(biogeography_key$region)] <- "";
biogeography_key$province[is.na(biogeography_key$province)] <- "";

excluded_species <- c();
# if you are excluding species, then upload a file for doing so.
# In the example below, species known only from pygidia are excluded: cranidia must be present to be included in the analyses
#excluded_species <- as.data.frame(readxl::read_xlsx ("Pygidium_Only.xlsx"));

# Start Main Program ####
if (length(total_group_daughters)==1 && total_group_daughters[1]==total_group)	{
	rank_col <- match(total_group_rank,colnames(pbdb_finds));
	control_finds <- pbdb_finds[pbdb_finds$class %in% total_group,];
	rank_col2 <- match(daughter_ranks[dr],colnames(pbdb_taxonomy));
	control_taxonomy <- pbdb_taxonomy[pbdb_taxonomy[,rank_col2]==total_group,];
	} else	{
	control_taxonomy <- pbdb_taxonomy[pbdb_taxonomy$taxon_name==total_group,];
	control_taxonomy <- rbind(control_taxonomy,pbdb_taxonomy[pbdb_taxonomy$parent_name %in% total_group,]);
	daughter_ranks <- unique(total_group_daughters$daughter_rank);
	control_finds <- pbdb_finds[pbdb_finds$collection_no==1,];
	control_finds <- control_finds[control_finds$collection_no>1,];
	for (dr in 1:length(daughter_ranks))	{
		rank_col <- match(daughter_ranks[dr],colnames(pbdb_finds));
		rank_col2 <- match(daughter_ranks[dr],colnames(pbdb_taxonomy));
		these_daughters <- total_group_daughters[total_group_daughters$daughter_rank==daughter_ranks[dr],];
		for (td in 1:nrow(these_daughters))	{
			control_finds <- rbind(control_finds,pbdb_finds[pbdb_finds[,rank_col] %in% total_group_daughters$daughter[td],]);
			control_taxonomy <- rbind(control_taxonomy,pbdb_taxonomy[pbdb_taxonomy[,rank_col2] %in% total_group_daughters$daughter[td],])
			}
		}
	}
control_sites <- pbdb_sites[pbdb_sites$collection_no %in% control_finds$collection_no,];
dud_taxa <- control_taxonomy$accepted_name[control_taxonomy$taxon_name %in% excluded_species$species]
control_taxonomy <- control_taxonomy[!control_taxonomy$accepted_name %in% dud_taxa,];
control_finds <- control_finds[!control_finds$accepted_name %in% dud_taxa,];

control_species <- control_taxonomy[control_taxonomy$accepted_rank %in% c("species","subspecies"),];
control_species$epithets <- sapply(control_species$accepted_name,divido_species_epithets);
control_species$accepted_genus <- sapply(control_species$accepted_name,divido_genus_names_from_species_names);
alt_combo1 <- alt_combo2 <- alt_combo3 <- paste(control_species$accepted_genus,control_species$epithets);
subgenus_rows <- (1:nrow(control_species))[sapply(control_species$accepted_genus,is.subgenus)];
notsubgenus_rows <- (1:nrow(control_species))[!(1:nrow(control_species)) %in% subgenus_rows];
#nn <- base::t(sapply(control_species$accepted_genus[sapply(control_species$accepted_genus,is.subgenus)],divido_subgenus_names_from_genus_names));
#mm <- base::t(sapply(control_species$genus[sapply(control_species$genus,is.subgenus)],divido_subgenus_names_from_genus_names));
nn <- base::t(sapply(control_species$accepted_genus,divido_subgenus_names_from_genus_names));
alt_combo2[subgenus_rows] <- paste(nn[subgenus_rows,1],control_species$epithets[subgenus_rows]);
alt_combo2[notsubgenus_rows] <- paste(nn[notsubgenus_rows,1]," (",nn[notsubgenus_rows,1],") ",control_species$epithets[notsubgenus_rows],sep="");
alt_combo3[subgenus_rows] <- paste(nn[subgenus_rows,2],control_species$epithets[subgenus_rows]);
alt_combo3[notsubgenus_rows] <- paste(nn[notsubgenus_rows,1]," (",nn[notsubgenus_rows,1],") ",control_species$epithets[notsubgenus_rows],sep="");
alt_combo4 <- alt_combo5 <- alt_combo6 <- paste(control_species$genus,control_species$epithets);
subgenus_rows_2 <- (1:nrow(control_species))[sapply(control_species$genus,is.subgenus)];
notsubgenus_rows_2 <- (1:nrow(control_species))[!(1:nrow(control_species)) %in% subgenus_rows_2];
mm <- base::t(sapply(control_species$genus,divido_subgenus_names_from_genus_names));
alt_combo5[subgenus_rows_2] <- paste(mm[subgenus_rows_2,1],control_species$epithets[subgenus_rows_2]);
alt_combo5[notsubgenus_rows_2] <- paste(mm[notsubgenus_rows_2,1]," (",mm[notsubgenus_rows_2,1],") ",control_species$epithets[notsubgenus_rows_2],sep="");
alt_combo6[subgenus_rows_2] <- paste(mm[subgenus_rows_2,2],control_species$epithets[subgenus_rows_2]);
alt_combo6[notsubgenus_rows_2] <- paste(mm[notsubgenus_rows_2,1]," (",mm[notsubgenus_rows_2,1],") ",control_species$epithets[notsubgenus_rows_2],sep="");

control_species <- cbind(control_species,alt_combo1,alt_combo2,alt_combo3,alt_combo4,alt_combo5,alt_combo6);
#which(control_species %in% pygidial_duds,arr.ind = T)[,1]

#control_finds$accepted_name[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no]
control_finds$accepted_name[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- control_taxonomy$accepted_name[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)];
control_finds$subgenus_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- as.numeric(control_taxonomy$subgenus_no[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)]);
control_finds$genus[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- control_taxonomy$genus[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)];
control_finds$genus_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- as.numeric(control_taxonomy$genus_no[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)]);
control_finds$family[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- control_taxonomy$family[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)];
control_finds$family_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- as.numeric(control_taxonomy$family_no[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)]);
control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no] <- as.numeric(control_taxonomy$accepted_no[match(control_finds$accepted_no[control_finds$accepted_no %in% paleodb_taxonomy_edits$taxon_no],control_taxonomy$taxon_no)]);

control_taxonomy <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank="superfamily",pbdb_taxonomy=control_taxonomy);
control_taxonomy <- add_higher_taxon_to_pbdb_taxonomy("subfamily",control_taxonomy);

# add faux subfamilies to taxa assigned only to a family
if (!is.null(control_taxonomy$subfamily))	{
	for (pt in 1:nrow(control_taxonomy))	{
		if (control_taxonomy$subfamily[pt]=="" && control_taxonomy$accepted_rank[pt] %in% c("genus","subgenus"))	{
			taxon <- gsub(" ","%20",control_taxonomy$taxon_name[pt]);
			httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
			new_info <- read.csv(httpTt,header=T,stringsAsFactors = F);
			new_info <- new_info[new_info$accepted_name==taxon,];
			if (nrow(new_info)>0)	{
				if (!is.na(match(new_info$parent_name[1],control_taxonomy$taxon_name)) & control_taxonomy$accepted_rank[match(new_info$parent_name[1],control_taxonomy$taxon_name)]=="subfamily")	{
					control_taxonomy$subfamily[pt] <- new_info$parent_name[1];
					control_taxonomy$subfamily_no[pt] <- new_info$parent_no[1];
					control_taxonomy$parent_name[pt] <- new_info$parent_name[1];
					control_taxonomy$parent_no[pt] <- new_info$parent_no[1];
					} else if (gsub("idae","inae",control_taxonomy$family[pt]) %in% control_taxonomy$accepted_name)	{
					faux_subfamily <- gsub("idae","inae",control_taxonomy$family[pt]);
					control_taxonomy$subfamily[pt] <- faux_subfamily;
					control_taxonomy$subfamily_no[pt] <- control_taxonomy$taxon_no[control_taxonomy$accepted_name==faux_subfamily][1];
					}
				}
			}
		}
	}

# make sure that we use either Lophospira or Lophospira (Lophospira) consistently.
control_taxonomy$taxon_name_orig <- control_taxonomy$taxon_name;
control_subgenera <- control_taxonomy[control_taxonomy$taxon_rank=="subgenus" & control_taxonomy$difference=="",];
subgenera <- sort(control_taxonomy$taxon_name[control_taxonomy$taxon_rank=="subgenus" & control_taxonomy$difference==""]);
check <- c();
for (sg in 1:length(subgenera))	{
	taxon_names <- divido_subgenus_names_from_genus_names(subgenera[sg]);
	taxon_name <- taxon_names[2];
	if (!taxon_name %in% control_taxonomy$taxon_name)	{
		dummy <- control_taxonomy[control_taxonomy$taxon_name==subgenera[sg],];
		dummy$taxon_name <- taxon_name;
		dummy$taxon_rank <- "genus";
		dummy$difference <- "recombined as";
		control_taxonomy <- rbind(control_taxonomy,dummy);
		} else if (taxon_names[1]==taxon_names[2] && sum(control_subgenera$parent_name %in% taxon_names[2])>1)	{
#		proetide_subgenera[proetide_subgenera$parent_name %in% taxon_names[2],]
		control_taxonomy$genus[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==taxon_name] <- subgenera[sg];
		control_taxonomy$taxon_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==taxon_name] <- gsub(taxon_name,subgenera[sg],control_taxonomy$taxon_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==taxon_name]);
		control_taxonomy$accepted_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==taxon_name] <- gsub(taxon_name,subgenera[sg],control_taxonomy$accepted_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==taxon_name]);
		control_taxonomy$difference[control_taxonomy$taxon_name==taxon_name] <- "see subgenus";
		control_taxonomy$accepted_name[control_taxonomy$taxon_name==taxon_name] <- subgenera[sg];
		control_finds$accepted_name[control_finds$genus %in% taxon_name] <- gsub(taxon_name,subgenera[sg],control_finds$accepted_name[control_finds$genus %in% taxon_name]);
		control_finds$genus[control_finds$genus %in% taxon_name] <- gsub(taxon_name,subgenera[sg],control_finds$genus[control_finds$genus %in% taxon_name]);
		} else if (taxon_names[1]==taxon_names[2])	{
		subdummy <- paste(taxon_name," \\(",taxon_name,"\\)",sep="")
		control_taxonomy$genus[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==subgenera[sg]] <- taxon_name;
		control_taxonomy$taxon_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==subgenera[sg]] <- gsub(subdummy,taxon_name,control_taxonomy$taxon_name[control_taxonomy$taxon_rank %in% c("species","subspecies") & control_taxonomy$parent_name==subgenera[sg]]);
		control_taxonomy$difference[control_taxonomy$taxon_name==subgenera[sg]] <- "obviated subgenus";
#		proetide_taxonomy$accepted_name[proetide_taxonomy$taxon_name==subgenera[sg]] <- taxon_name;
		control_taxonomy$accepted_name[control_taxonomy$genus %in% subgenera[sg]] <- gsub(subdummy,taxon_name,control_finds$accepted_name[control_finds$genus %in% subgenera[sg]]);
		control_taxonomy$accepted_rank[control_taxonomy$taxon_name==subgenera[sg]] <- "genus";
		control_finds$accepted_name[control_finds$genus %in% subgenera[sg]] <- gsub(subdummy,taxon_name,control_finds$accepted_name[control_finds$genus %in% subgenera[sg]]);
		control_finds$genus[control_finds$genus %in% subgenera[sg]] <- taxon_name;
		} else if (control_taxonomy$difference[match(subgenera[sg],control_taxonomy$taxon_name)]=="" && control_taxonomy$difference[match(taxon_name,control_taxonomy$taxon_name)]=="") {
		check <- c(check,taxon_name,subgenera[sg]);
		}
	}

# make sure that we can find the requested higher taxa in PBDB data!
total_group_order <- control_taxonomy$order[control_taxonomy$taxon_name==total_group];
jackobites <- sepkoski_compendium[sepkoski_compendium$Order %in% total_group_order,];
if (length(ingroup_daughters)==2 && ingroup_daughters$daughter[1]==ingroup[length(ingroup)])	{
	rank_col <- match(total_group_rank,colnames(pbdb_finds));
#	control_finds <- pbdb_finds[pbdb_finds$class %in% total_group,];
	rank_col2 <- match(daughter_ranks[dr],colnames(pbdb_taxonomy));
	ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy[,rank_col2]==total_group,];
	} else	{
	ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% ingroup,];
	ingroup_taxonomy <- rbind(ingroup_taxonomy,pbdb_taxonomy[pbdb_taxonomy$parent_name %in% total_group,]);
	daughter_ranks <- unique(total_group_daughters$daughter_rank);
#	control_finds <- pbdb_finds[pbdb_finds$collection_no==1,];
#	control_finds <- control_finds[control_finds$collection_no>1,];
	for (dr in 1:length(daughter_ranks))	{
		rank_col <- match(daughter_ranks[dr],colnames(pbdb_finds));
		rank_col2 <- match(daughter_ranks[dr],colnames(pbdb_taxonomy));
		these_daughters <- total_group_daughters[total_group_daughters$daughter_rank==daughter_ranks[dr],];
		for (td in 1:nrow(these_daughters))	{
#			control_finds <- rbind(control_finds,pbdb_finds[pbdb_finds[,rank_col] %in% total_group_daughters$daughter[td],]);
			ingroup_taxonomy <- rbind(ingroup_taxonomy,pbdb_taxonomy[pbdb_taxonomy[,rank_col2] %in% total_group_daughters$daughter[td],])
			}
		}
	}

pbdb_control_taxonomy <- control_taxonomy[control_taxonomy$accepted_rank %in% c("genus","subgenus"),];
pbdb_control_valid <- pbdb_control_taxonomy[pbdb_control_taxonomy$difference %in% "",];
ctaxa <- nrow(pbdb_control_valid);

# Get Character Data ####
phylogenetic_info <- accersi_data_from_chosen_nexus_file();
analyzed_species <- phylogenetic_info$OTUs;
analyzed_species <- unique(gsub("redo ","",analyzed_species));
notu <- length(analyzed_species);
analyzed_species <- gsub("_"," ",analyzed_species);
analyzed_genera <- analyzed_genera_orig <- sapply(analyzed_species,divido_genus_names_from_species_names);
analyzed_genera[!is.na(match(analyzed_genera_orig,pbdb_control_taxonomy$taxon_name))] <- pbdb_control_taxonomy$accepted_name[match(analyzed_genera_orig[!is.na(match(analyzed_genera,pbdb_control_taxonomy$taxon_name))],pbdb_control_taxonomy$taxon_name)];
analyzed_families <- unique(control_taxonomy$family[control_taxonomy$taxon_name %in% analyzed_genera]);
for (ks in 1:length(analyzed_species))	{
	nn <- control_species$accepted_name[match(analyzed_species[ks],control_species$taxon_name)];
	if (is.na(nn))	nn <- control_species$accepted_name[match(analyzed_species[ks],control_species$taxon_name_orig)];
	if (is.na(nn))	nn <- control_species$accepted_name[match(analyzed_species[ks],control_species$accepted_name)];
	if (is.na(nn))	{
		nn <- gsub(analyzed_genera_orig[ks],analyzed_genera[ks],analyzed_species[ks]);
		nn <- control_species$accepted_name[match(nn,control_species$taxon_name)];
		}
	if (is.na(nn))	nn <- unique(control_species$accepted_name[which(control_species==analyzed_species[ks],arr.ind = T)[,1]]);
	if (is.na(nn) || length(nn)==0)	{
		genus_ks <- divido_genus_names_from_species_names(analyzed_species[ks]);
		species_ks <- divido_species_epithets(analyzed_species[ks]);
		genus_subgenus <- divido_subgenus_names_from_genus_names(genus_ks);
		if (genus_subgenus[2]!="")	{
			combo_1 <- paste(genus_subgenus[1],species_ks);
			combo_2 <- paste(genus_subgenus[2],species_ks);
			nn <- unique(c(unique(control_species$accepted_name[unique(which(control_species==combo_1,arr.ind=T)[,1])]),
										 unique(control_species$accepted_name[unique(which(control_species==combo_2,arr.ind=T)[,1])])));
			} else	{
			combo_1 <- paste(genus_subgenus[1]," (",genus_subgenus[1],") ",species_ks,sep="");
			nn <- control_species$accepted_name[match(analyzed_species[ks],control_species$taxon_name)];
			if (is.na(nn))
				nn <- unique(control_species$accepted_name[unique(which(control_species==combo_1,arr.ind=T)[,1])]);
			}
		}
	if (!is.na(nn) && length(nn)==1)	analyzed_species[ks] <- nn;
	}

# get biogeographic information 
study_sites <- control_sites[control_sites$ma_ub<myr_lb & control_sites$ma_lb>myr_ub,];
study_sites$realm <- rock_database$realm[match(study_sites$rock_no_sr,rock_database$rock_no)];
study_sites$region <- rock_database$region[match(study_sites$rock_no_sr,rock_database$rock_no)];
study_sites$region[study_sites$region==""] <- study_sites$realm[study_sites$region==""]
study_sites$region[study_sites$region=="All"] <- "";
study_plates <- sort(unique(study_sites$geoplate));
for (sp in 1:length(study_plates))	{
	site_regions <- study_sites$region[study_sites$geoplate==study_plates[sp]];
	plate_regions <- unique(study_sites$region[study_sites$geoplate==study_plates[sp]]);
	if ("" %in% plate_regions)	{
		plate_regions <- plate_regions[plate_regions!=""];
		site_regions <- site_regions[site_regions!=""];
		if (length(plate_regions)>0)	{
			region_counts <- hist(match(site_regions,plate_regions),breaks=0:length(plate_regions),plot=F)$counts;
			names(region_counts) <- plate_regions;
			best_region <- names(region_counts)[match(max(region_counts),region_counts)];
			best_realm <- biogeography_key$realm[match(best_region,biogeography_key$region)];
			study_sites$realm[study_sites$geoplate==study_plates[sp] & study_sites$region==""] <- best_realm;
			study_sites$region[study_sites$geoplate==study_plates[sp] & study_sites$region==""] <- best_region;
			}
		}
	}

# Get summary info for total group
control_group_output <- data.frame(genus=as.character(pbdb_control_valid$taxon_name),
								   subfamily=as.character(pbdb_control_valid$subfamily),family=as.character(pbdb_control_valid$family),
								   superfamily=as.character(pbdb_control_valid$superfamily),
								   nspc=as.numeric(rep(0,ctaxa)),nspc_unentered=as.numeric(rep(0,ctaxa)),
								   coded_spc=as.numeric(rep(0,ctaxa)),nfinds=as.numeric(rep(0,ctaxa)),
								   fa_lb=as.numeric(rep(0,ctaxa)),fa_ub=as.numeric(rep(0,ctaxa)),
								   la_lb=as.numeric(rep(0,ctaxa)),la_ub=as.numeric(rep(0,ctaxa)),
								   fa_sepkoski=as.numeric(rep(0,ctaxa)),la_sepkoski=as.numeric(rep(0,ctaxa)),
								   fa=as.numeric(rep(0,ctaxa)),la=as.numeric(rep(0,ctaxa)));
target_control_taxa <- control_group_output;
missing_species <- missing_target_species <- c();
dummy_info <- data.frame(genus=as.character(),species=as.character(),author=as.character(),
												 nfinds=as.numeric(),nrocks=as.numeric(),
												 fa=as.character(),la=as.character(),
												 geoplates=as.character(),region=as.character());
uncoded_target_species <- dummy_info;
coded_target_species <- data.frame(genus=as.character(analyzed_genera),species=as.character(analyzed_species),author=as.character(rep("",notu)),
								   nfinds=as.numeric(rep(0,notu)),nrocks=as.numeric(rep(0,notu)),
								   fa=as.character(rep("",notu)),la=as.character(rep("",notu)),
								   geoplates=as.character(rep("",notu)),region=as.character(rep("",notu)));
for (np in 1:ctaxa)	{
	genus_finds <- control_finds[control_finds$genus %in% pbdb_control_valid$taxon_name[np],];
	if (nrow(genus_finds)==0 && !is.subgenus(pbdb_control_valid$taxon_name[np]) && pbdb_control_valid$n_occs[np]>0)	{
		eponymous <- make_eponymous_subgenus(pbdb_control_valid$taxon_name[np]);
		if (eponymous %in% control_taxonomy$accepted_name)	{
			pbdb_control_valid[np,] <- control_taxonomy[control_taxonomy$accepted_name==eponymous,];
			} else	{
			pbdb_control_valid$taxon_name[np] <- eponymous;
			}
		genus_finds <- control_finds[control_finds$genus %in% pbdb_control_valid$taxon_name[np],];
		}
	genus_finds <- genus_finds[genus_finds$identified_rank %in% c("species","subspecies"),];
	genus_finds <- expello_indeterminate_species(genus_finds,"identified_name");
	genus_finds$identified_name <- sapply(genus_finds$identified_name,mundify_taxon_names);
	genus_finds$accepted_name <- sapply(genus_finds$accepted_name,mundify_taxon_names);
	unique_id_accepted_combos <- unique(genus_finds[,colnames(genus_finds) %in% c("identified_name","accepted_name")]);
	unique_id_accepted_combos$species_epithet <- sapply(unique_id_accepted_combos$identified_name,divido_species_epithets);
	unique_id_accepted_combos <- unique_id_accepted_combos[order(unique_id_accepted_combos$species_epithet),];
	entered_species <- unique_id_accepted_combos$identified_name[!unique_id_accepted_combos$accepted_name %in% pbdb_control_valid$taxon_name[np]];
#	entered_species <- unique_id_accepted_combos$identified_name[!unique_id_accepted_combos$accepted_name %in% paste(pbdb_proetoids_valid$taxon_name[np]," ",sep="")]
	hoffas <- unique_id_accepted_combos$identified_name[unique_id_accepted_combos$accepted_name %in% pbdb_control_valid$taxon_name[np]];
	missing_species <-rbind(missing_species,cbind(rep(pbdb_control_valid$taxon_name[np],length(hoffas)),hoffas,rep(pbdb_control_valid$family[np],length(hoffas))));
	unique_species <- unique(c(hoffas,entered_species));
	control_group_output$coded_spc[np] <- sum(analyzed_genera %in% pbdb_control_valid$taxon_name[np]);
	control_group_output$nspc[np] <- length(unique_species);
	control_group_output$nspc_unentered[np] <- length(hoffas);
	genus_sites <- study_sites[study_sites$collection_no %in% genus_finds$collection_no,];
	genus_sites$realm[is.na(genus_sites$realm)] <- genus_sites$region[is.na(genus_sites$realm)];
	control_group_output$nfinds[np] <- nrow(genus_sites);
	if (nrow(genus_sites)>0)	{
		control_group_output$fa_lb[np] <- max(genus_sites$ma_lb);
		control_group_output$fa_ub[np] <- max(genus_sites$ma_ub);
		control_group_output$la_lb[np] <- min(genus_sites$ma_lb);
		control_group_output$la_ub[np] <- min(genus_sites$ma_ub);
		}
	if (!is.na(match(pbdb_control_valid$taxon_name[np],jackobites$Genus)))	{
		control_group_output$fa_sepkoski[np] <- time_scale$ma_lb[match(jackobites$FA[match(pbdb_control_valid$taxon_name[np],jackobites$Genus)],time_scale$interval)];
		control_group_output$la_sepkoski[np] <- time_scale$ma_ub[match(jackobites$LA[match(pbdb_control_valid$taxon_name[np],jackobites$Genus)],time_scale$interval)];
		} else if (!is.na(match(pbdb_control_valid$accepted_name[np],jackobites$Genus)))	{
		control_group_output$fa_sepkoski[np] <- time_scale$ma_lb[match(jackobites$FA[match(pbdb_control_valid$taxon_name[np],jackobites$Genus)],time_scale$interval)];
		control_group_output$la_sepkoski[np] <- time_scale$ma_ub[match(jackobites$LA[match(pbdb_control_valid$taxon_name[np],jackobites$Genus)],time_scale$interval)];
		}
	# case where we don't have any PBDB dates
	if (control_group_output$fa_lb[np]==0)	{
		control_group_output$fa[np] <- control_group_output$fa_sepkoski[np];
		control_group_output$la[np] <- control_group_output$la_sepkoski[np];
		} else	{
		control_group_output$fa[np] <- control_group_output$fa_lb[np];
		control_group_output$la[np] <- control_group_output$la_ub[np];
		}
	if (control_group_output$fa[np]>myr_ub & control_group_output$la[np]<myr_lb)	{
		target_sites <- genus_sites[genus_sites$ma_lb>myr_ub & genus_sites$ma_ub<myr_lb,];
		if (nrow(target_sites)>0)	{
			target_finds <- genus_finds[genus_finds$collection_no %in% target_sites$collection_no,];
			unique_id_accepted_combos <- unique(target_finds[,colnames(target_finds) %in% c("identified_name","accepted_name")]);
			unique_id_accepted_combos$species_epithet <- sapply(unique_id_accepted_combos$identified_name,divido_species_epithets);
			unique_id_accepted_combos <- unique_id_accepted_combos[order(unique_id_accepted_combos$species_epithet),];
			entered_species <- unique(unique_id_accepted_combos$accepted_name[!unique_id_accepted_combos$accepted_name %in% pbdb_control_valid$taxon_name[np]]);
			hoffas <- unique_id_accepted_combos$identified_name[unique_id_accepted_combos$accepted_name %in% pbdb_control_valid$taxon_name[np]]
			unique_species <- unique(c(hoffas,entered_species));
			missing_target_species <- rbind(missing_target_species,cbind(rep(pbdb_control_valid$taxon_name[np],length(hoffas)),hoffas,rep(pbdb_control_valid$family[np],length(hoffas))));
			target_control_taxa$nspc[np] <- length(unique_species);
			target_control_taxa$nspc_unentered[np] <- length(hoffas);
			target_control_taxa$coded_spc[np] <- control_group_output$coded_spc[np];
			target_control_taxa$nfinds[np] <- nrow(target_finds);
			target_control_taxa$fa[np] <- target_control_taxa$fa_lb[np] <- max(target_sites$ma_lb);
			target_control_taxa$fa_ub[np] <- max(target_sites$ma_ub);
			target_control_taxa$la_lb[np] <- min(target_sites$ma_lb);
			target_control_taxa$la[np] <- target_control_taxa$la_ub[np] <- min(target_sites$ma_ub);
			target_control_taxa$fa_sepkoski[np] <- control_group_output$fa_sepkoski[np];
			target_control_taxa$la_sepkoski[np] <- control_group_output$la_sepkoski[np];
			
			add_these <- dummy_info;
			us <- 0;
			while (us < length(unique_species))	{
				us <- us+1;
				add_me <- data.frame(genus=as.character(pbdb_control_valid$taxon_name[np]),
														 species=as.character(unique_species[us]),
														 author=as.character(""),
														 nfinds=as.numeric(0),nrocks=as.numeric(0),
														 fa=as.character(""),la=as.character(""),
														 geoplates=as.character(""),region=as.character(""));
				taxon_row <- unique(which(control_taxonomy==unique_species[us],arr.ind=T)[,1]);
				if (length(taxon_row)==0 & control_taxonomy$accepted_rank[match(pbdb_control_valid$taxon_name[np],control_taxonomy$accepted_name)]=="genus")	{
					dummy_subgenus <- paste(pbdb_control_valid$taxon_name[np]," (",pbdb_control_valid$taxon_name[np],")",sep="");
					dummy_species <- gsub(pbdb_control_valid$taxon_name[np],dummy_subgenus,unique_species[us])
					taxon_row <- unique(which(control_taxonomy==dummy_species,arr.ind=T)[,1]);
					}
				if (length(taxon_row)>1)
					taxon_row <- taxon_row[control_taxonomy$taxon_no[taxon_row]==control_taxonomy$accepted_no[taxon_row]];
				if (length(taxon_row)>1)
					taxon_row <- taxon_row[control_taxonomy$accepted_rank[taxon_row] %in% c("species","subspecies")];
				if (length(taxon_row)>1)	{
					if (!is.subspecies(unique_species[us]))
						taxon_row <- taxon_row[control_taxonomy$accepted_rank[taxon_row] %in% "species"];
					}
				if (length(taxon_row)==1)	{
					add_me$genus <- control_taxonomy$genus[taxon_row];
					add_me$species <- control_taxonomy$accepted_name[taxon_row];
					add_me$author <- control_taxonomy$taxon_attr[taxon_row];
					taxon_sites <- target_sites[target_sites$collection_no %in% unique(target_finds$collection_no[target_finds$accepted_no %in% control_taxonomy$taxon_no[taxon_row]]),];
					if (nrow(taxon_sites)==0)
						taxon_sites <- target_sites[target_sites$collection_no %in% unique(target_finds$collection_no[target_finds$identified_name %in% unique_species[us]]),];
					if (nrow(taxon_sites)>0)	{
						add_me$nfinds <- nrow(taxon_sites);
						add_me$fa <- taxon_sites$interval_lb[match(max(taxon_sites$ma_lb),taxon_sites$ma_lb)];
						add_me$la <- taxon_sites$interval_ub[match(min(taxon_sites$ma_ub),taxon_sites$ma_ub)];
						add_me$nrocks <- length(unique(taxon_sites$rock_no_sr));
						if (nrow(taxon_sites)==1)	{
							add_me$geoplates <- paste(unique(taxon_sites$geoplate),collapse=", ");
							regions <- unique(rock_database$region[match(taxon_sites$rock_no_sr[taxon_sites$rock_no_sr>0],rock_database$rock_no)]);
							add_me$region <- paste(regions,collapse=", ");
							} else	{
							oldest_sites <- taxon_sites[taxon_sites$ma_ub>=max(taxon_sites$ma_ub),];
							add_me$geoplates <- paste(unique(oldest_sites$geoplate),collapse=", ");
							regions <- unique(rock_database$region[match(oldest_sites$rock_no_sr[oldest_sites$rock_no_sr>0],rock_database$rock_no)]);
							add_me$region <- paste(regions,collapse=", ");
							}
#						if (sum(taxon_sites$rock_no_sr>0)>0)	{
#							regions <- unique(rock_database$region[match(taxon_sites$rock_no_sr[taxon_sites$rock_no_sr>0],rock_database$rock_no)]);
#							add_me$region <- paste(regions,collapse=", ");
#							}
						}
					} else if (length(taxon_row)==0)	{
					taxon_sites <- target_sites[target_sites$collection_no %in% unique(target_finds$collection_no[target_finds$identified_name %in% unique_species[us]]),];
					if (nrow(taxon_sites)==0) { #& ingroup_taxonomy$accepted_rank[match(pbdb_proetoids_valid$taxon_name[np],ingroup_taxonomy$accepted_name)]=="genus")	{
						dummy_subgenus <- paste(pbdb_control_valid$taxon_name[np]," (",pbdb_control_valid$taxon_name[np],")",sep="");
						dummy_species <- gsub(pbdb_control_valid$taxon_name[np],dummy_subgenus,unique_species[us])
						dummy_species_2 <- paste(pbdb_control_valid$taxon_name[np],divido_species_epithets(unique_species[us]));
						take_these <- unique(c(target_finds$collection_no[target_finds$identified_name %in% c(dummy_species,dummy_species_2)],
																	 target_finds$collection_no[target_finds$accepted_name %in% c(dummy_species,dummy_species_2)]));
						taxon_sites <- target_sites[target_sites$collection_no %in% take_these,];
						}
					add_me$nfinds <- nrow(taxon_sites);
					add_me$fa <- taxon_sites$interval_lb[match(max(taxon_sites$ma_lb),taxon_sites$ma_lb)];
					add_me$la <- taxon_sites$interval_ub[match(min(taxon_sites$ma_ub),taxon_sites$ma_ub)];
					oldest_sites <- taxon_sites[taxon_sites$ma_ub>=max(taxon_sites$ma_ub),];
					add_me$geoplates <- paste(unique(oldest_sites$geoplate),collapse=", ");
					if (sum(taxon_sites$rock_no_sr>0)>0)	{
						add_me$nrocks <- length(unique(taxon_sites$rock_no_sr));
						regions <- unique(rock_database$region[match(oldest_sites$rock_no_sr[oldest_sites$rock_no_sr>0],rock_database$rock_no)]);
						add_me$region <- paste(regions,collapse=", ");
						}
					} else if (length(taxon_row)>1)	{
					print(paste(np,us));
					}
				if (!add_me$species %in% analyzed_species)	{
					add_these <- rbind(add_these,add_me);
					} else	{
					coded_target_species[match(add_me$species,analyzed_species),] <- add_me;
#					coded_target_species <- rbind(coded_target_species,add_these);
					}
				}
			if (nrow(add_these)>0)	{
				add_these$nrocks[add_these$nrocks==0] <- 1;
				add_these <- add_these[order(-(rank(add_these$nrocks)+rank(add_these$nfinds))/2,add_these$species),];
				uncoded_target_species <- rbind(uncoded_target_species,add_these);
				}
			}
		if (control_group_output$fa[np]==0)	{
			target_control_taxa$fa[np] <- min(myr_lb,control_group_output$fa[np]);
			target_control_taxa$la[np] <- max(myr_ub,control_group_output$la[np]);
			}
		}
	}

# make sure that we have information for coded species
for (cs in 1:notu)	{
	taxon_row <- c();
	if (coded_target_species$author[cs]=="" || is.na(coded_target_species$author[cs]))	{
		taxon_row <- unique(which(ingroup_taxonomy==coded_target_species$species[cs],arr.ind=T)[,1]);
		if (length(taxon_row)==0)	{
			if (is.subgenus(divido_genus_names_from_species_names(coded_target_species$species[cs])))	{
				this_genus <- divido_genus_names_from_species_names(coded_target_species$species[cs]);
				genus_names <- divido_subgenus_names_from_genus_names(divido_genus_names_from_species_names(coded_target_species$species[cs]));
				dummy_genus <- paste(genus_names[1]," \\(",genus_names[2],"\\)",sep="")
				if (length(unique(genus_names))==1)	{
					alt_names <- gsub(dummy_genus,genus_names[1],coded_target_species$species[cs]);
					taxon_row <- unique(which(ingroup_taxonomy==alt_names,arr.ind=T)[,1]);
					} else	{
					alt_names <- gsub(dummy_genus,genus_names[1],coded_target_species$species[cs]);
					taxon_row <- unique(which(ingroup_taxonomy==alt_names,arr.ind=T)[,1]);
					alt_names <- c(alt_names,gsub(dummy_genus,genus_names[2],coded_target_species$species[cs]));
					taxon_row <- c(taxon_row,unique(which(ingroup_taxonomy==alt_names[2],arr.ind=T)[,1]));
					}
				}
			if (length(taxon_row)>1)	taxon_row <- taxon_row[ingroup_taxonomy$accepted_rank[taxon_row] %in% c("species","subspecies")];
			if (length(taxon_row)>1)	taxon_row <- taxon_row[ingroup_taxonomy$accepted_name[taxon_row]==ingroup_taxonomy$taxon_name[taxon_row]];
			if (length(taxon_row)==1) coded_target_species$author[cs] <- ingroup_taxonomy$taxon_attr[taxon_row];
			}
		if (length(taxon_row==0))	{
			taxon_row <- unique(which(control_taxonomy==coded_target_species$species[cs],arr.ind=T)[,1]);
			if (length(taxon_row)>0)	{
				if (length(unique(control_taxonomy$accepted_name[taxon_row]))>1)	{
					if (control_taxonomy$accepted_rank[taxon_row] %in% c("species","subspecies"))	{
						taxon_row <- taxon_row[control_taxonomy$accepted_rank[taxon_row] %in% "species"];
						}
					}
				if (length(unique(control_taxonomy$accepted_name[taxon_row]))==1)	{
					coded_target_species$author[cs] <- control_taxonomy$taxon_attr[taxon_row[1]];
					}
				}
			}
		}
	if (coded_target_species$nfinds[cs]==0)	{
		taxon_row <- unique(which(ingroup_taxonomy==coded_target_species$species[cs],arr.ind=T)[,1]);
		taxon_row <- taxon_row[ingroup_taxonomy$accepted_rank[taxon_row] %in% c("species","subspecies")];
		ingroup_taxonomy[taxon_row,]
		if (sum(c("species","subspecies") %in% ingroup_taxonomy$accepted_rank[taxon_row])==2 && length(unique(ingroup_taxonomy$accepted_name[taxon_row]))>1)	{
			keeper <- ingroup_taxonomy$accepted_name[taxon_row[ingroup_taxonomy$taxon_name[taxon_row]==coded_target_species$species[cs]]];
			if (length(keeper)==1)
				taxon_row <- taxon_row[ingroup_taxonomy$accepted_name[taxon_row] %in% keeper]
			}
		if (length(taxon_row)==0)	{
			if (is.subgenus(divido_genus_names_from_species_names(coded_target_species$species[cs])))	{
				this_genus <- divido_genus_names_from_species_names(coded_target_species$species[cs]);
				genus_names <- divido_subgenus_names_from_genus_names(divido_genus_names_from_species_names(coded_target_species$species[cs]));
				dummy_genus <- paste(genus_names[1]," \\(",genus_names[2],"\\)",sep="")
				if (length(unique(genus_names))==1)	{
					alt_names <- gsub(dummy_genus,genus_names[1],coded_target_species$species[cs]);
					taxon_row <- unique(which(ingroup_taxonomy==alt_names,arr.ind=T)[,1]);
					} else	{
					alt_names <- gsub(dummy_genus,genus_names[1],coded_target_species$species[cs]);
					taxon_row <- unique(which(ingroup_taxonomy==alt_names,arr.ind=T)[,1]);
					alt_names <- c(alt_names,gsub(dummy_genus,genus_names[2],coded_target_species$species[cs]));
					taxon_row <- c(taxon_row,unique(which(ingroup_taxonomy==alt_names[2],arr.ind=T)[,1]));
					}
				}
			}
		otu_finds <- control_finds[control_finds$accepted_no %in% ingroup_taxonomy$taxon_no[taxon_row],];
		coded_target_species$author[cs] <- ingroup_taxonomy$taxon_attr[taxon_row[1]];
		if (nrow(otu_finds)>0)	{
			otu_sites <- study_sites[unique(study_sites$collection_no) %in% otu_finds$collection_no,];
			otu_sites$realm <- otu_sites$region <- rep("",nrow(otu_sites));
			otu_sites$realm[otu_sites$rock_no_sr!=0] <- rock_database$realm[match(otu_sites$rock_no_sr[otu_sites$rock_no_sr!=0],rock_database$rock_no)];
			otu_sites$region[otu_sites$rock_no_sr!=0] <- rock_database$region[match(otu_sites$rock_no_sr[otu_sites$rock_no_sr!=0],rock_database$rock_no)];
			otu_sites$region[otu_sites$region==""] <- otu_sites$realm[otu_sites$region==""];
			coded_target_species$nfinds[cs] <- nrow(otu_sites);
			coded_target_species$nrocks[cs] <- length(unique(otu_sites$rock_no_sr));
			coded_target_species$fa[cs] <- otu_sites$interval_lb[match(max(otu_sites$ma_lb),otu_sites$ma_lb)];
			coded_target_species$la[cs] <- otu_sites$interval_ub[match(min(otu_sites$ma_ub),otu_sites$ma_ub)];
			coded_target_species$geoplates[cs] <- paste(unique(otu_sites$geoplate),collapse=",");
			if (coded_target_species$nrocks[cs]>0)	{
				coded_target_species$region[cs] <- paste(unique(otu_sites$region[!otu_sites$region %in% c("","All")]),collapse=", ");
				} else	{
				coded_target_species$nrocks[cs] <- 1;
				}
			}
		}
	}

# organize data for output
control_group_output <- control_group_output[order(control_group_output$fa,control_group_output$la,control_group_output$superfamily,control_group_output$family,control_group_output$subfamily),];
target_control_taxa <- target_control_taxa[order(target_control_taxa$fa,target_control_taxa$la,target_control_taxa$superfamily,target_control_taxa$family,target_control_taxa$subfamily),];
target_control_taxa <- target_control_taxa[target_control_taxa$fa>0,];

strict_target_control_taxa <- target_control_taxa[target_control_taxa[,ingroup_rank] %in% ingroup,];
strict_target_control_taxa <- strict_target_control_taxa[order(-strict_target_control_taxa$nspc),];
strict_uncoded_target_species <- uncoded_target_species[uncoded_target_species$genus %in% strict_target_control_taxa$genus,];
strict_uncoded_target_species$genus_no <- match(strict_uncoded_target_species$genus,strict_target_control_taxa$genus);
strict_uncoded_target_species <- strict_uncoded_target_species[order(strict_uncoded_target_species$genus_no,-strict_uncoded_target_species$nfinds,strict_uncoded_target_species$species),];
strict_uncoded_target_species$status <- rep("uncoded",nrow(strict_uncoded_target_species));
strict_uncoded_target_species$accepted_no <- control_species$accepted_no[match(strict_uncoded_target_species$species,control_species$taxon_name)]
strict_uncoded_target_species$accepted_no[is.na(strict_uncoded_target_species$accepted_no)] <- control_species$accepted_no[match(strict_uncoded_target_species$species[is.na(strict_uncoded_target_species$accepted_no)],control_species$taxon_name_orig)];
strict_uncoded_target_species <- strict_uncoded_target_species[!strict_uncoded_target_species$species %in% excluded_species$species,]
	#control_species[control_species$taxon_name=="Macroblepharum perstes",]

not_informals <- (1:nrow(strict_uncoded_target_species))[!sapply(strict_uncoded_target_species$species,revelare_informal_taxa)];
strict_uncoded_target_species <- strict_uncoded_target_species[not_informals,];
unauthored <- (1:nrow(strict_uncoded_target_species))[strict_uncoded_target_species$author==""];
#strict_uncoded_target_species$species[unauthored]

uncspecies <- nrow(strict_uncoded_target_species)
uc <- 0;
while (uc < uncspecies)	{
	uc <- uc+1;
	if (strict_uncoded_target_species$author[uc]=="" || is.na(strict_uncoded_target_species$author[uc]))	{
		taxon_row <- unique(which(control_species==strict_uncoded_target_species$species[uc],arr.ind=T)[,1]);
		if (length(taxon_row)==0)	{
			if (is.subgenus(divido_genus_names_from_species_names(strict_uncoded_target_species$species[uc])))	{
				this_genus <- divido_genus_names_from_species_names(strict_uncoded_target_species$species[uc]);
				genus_names <- divido_subgenus_names_from_genus_names(divido_genus_names_from_species_names(strict_uncoded_target_species$species[uc]));
				dummy_genus <- paste(genus_names[1]," \\(",genus_names[2],"\\)",sep="");
				dummy_subgenus <- paste(genus_names[1]," \\(",genus_names[1],"\\)",sep="");
				if (length(unique(genus_names))==1)	{
					alt_names <- gsub(dummy_genus,genus_names[1],strict_uncoded_target_species$species[uc]);
					taxon_row <- unique(which(control_species==alt_names,arr.ind=T)[,1]);
					} else if (genus_names[2]=="")	{
					alt_names <- gsub(dummy_subgenus,genus_names[1],strict_uncoded_target_species$species[uc]);
					} else	{
					alt_names <- gsub(dummy_genus,genus_names[1],strict_uncoded_target_species$species[uc]);
					taxon_row <- unique(which(control_species==alt_names,arr.ind=T)[,1]);
					alt_names <- c(alt_names,gsub(dummy_genus,genus_names[2],strict_uncoded_target_species$species[uc]));
					taxon_row <- c(taxon_row,unique(which(control_species==alt_names[2],arr.ind=T)[,1]));
					}
				}
			if (length(taxon_row)>1)	taxon_row <- taxon_row[control_species$accepted_rank[taxon_row] %in% c("species","subspecies")];
			if (length(taxon_row)>1)	taxon_row <- taxon_row[control_species$accepted_name[taxon_row]==control_species$taxon_name[taxon_row]];
			if (length(taxon_row)>1)	taxon_row <- taxon_row[control_species$accepted_name[taxon_row] %in% alt_names];
			}
    if (length(taxon_row)>1)
      taxon_row <- taxon_row[control_species$taxon_name[taxon_row]==control_species$accepted_name[taxon_row]];
    if (length(taxon_row)>1)  {
      authority_test <- gsub("\\(","",control_species$taxon_attr[taxon_row]);
      authority_test <- gsub("\\)","",authority_test);
      authority_test <- unique(authority_test);
      if (length(authority_test)==1)  {
        strict_uncoded_target_species$author[uc] <- paste("(",authority_test,")",sep="");
        } else  {
        print(uc);
        }
      }
		if (length(taxon_row)==1)
      strict_uncoded_target_species$author[uc] <- control_species$taxon_attr[taxon_row];
		}	
	if (is.na(strict_uncoded_target_species$accepted_no[uc]))	{
		poss_matches <- unique(which(control_species==strict_uncoded_target_species$species[uc],arr.ind=T)[,1]);
		if (length(poss_matches)==0)	{
			if (is.subgenus(divido_genus_names_from_species_names(strict_uncoded_target_species$species[uc])))	{
				genus_subgenus <- divido_subgenus_names_from_genus_names(divido_genus_names_from_species_names(strict_uncoded_target_species$species[uc]));
				epithet <- divido_species_epithets(strict_uncoded_target_species$species[uc]);
				combos <- paste(genus_subgenus[1],epithet);
				combos <- c(combos,paste(genus_subgenus[2],epithet));
				combos <- c(combos,paste(genus_subgenus[2]," (",genus_subgenus[2],") ",epithet,sep=""));
				poss_matches <- (1:nrow(control_species))[control_species$taxon_name %in% combos]
				}
			}
		if (length(poss_matches)>0 && length(unique(control_species$accepted_no[poss_matches]))==1)
			strict_uncoded_target_species$accepted_no[uc] <- control_species$accepted_no[poss_matches[1]];
		}
	}

coded_target_species$genus_no <- match(coded_target_species$genus,strict_target_control_taxa$genus);
nag <- sum(is.na(coded_target_species$genus_no));
ng <- 0;
nags <- (1:notu)[is.na(coded_target_species$genus_no)];
while (ng < nag)	{
	ng <- ng+1;
	genus_subgenus <- divido_subgenus_names_from_genus_names(coded_target_species$genus[ng]);
	if (genus_subgenus[2]!="")
		coded_target_species$genus_no[ng] <- match(genus_subgenus[2],strict_target_control_taxa$genus);
	}
coded_target_species$status <- rep("coded",nrow(coded_target_species));
coded_target_species$accepted_no <- control_species$accepted_no[match(coded_target_species$species,control_species$taxon_name)]
coded_target_species$accepted_no[is.na(coded_target_species$accepted_no)] <- control_species$accepted_no[match(coded_target_species$species[is.na(coded_target_species$accepted_no)],control_species$taxon_name_orig)];

# strict_uncoded_target_species lost!
all_target_species <- rbind(strict_uncoded_target_species,coded_target_species);
all_target_species <- all_target_species[!is.na(all_target_species$genus_no),];
target_control_taxa$nspc[target_control_taxa$nspc==0] <- 1;

# eliminate duplicates
all_target_species$accepted_no[is.na(all_target_species$accepted_no)] <- -max(all_target_species$accepted_no[!is.na(all_target_species$accepted_no)])+1000+(1:sum(is.na(all_target_species$accepted_no)));
all_target_species <- all_target_species[order(abs(all_target_species$accepted_no)),];
a_t_s <- nrow(all_target_species);
if (length(unique(all_target_species$accepted_no)) < a_t_s)	{
	unique_nos <- unique(all_target_species$accepted_no);
	max_no <- max(match(all_target_species$accepted_no,unique_nos));
	duplicate_test <- hist(match(all_target_species$accepted_no,unique_nos),breaks=0:max_no,plot=F)$counts;
	names(duplicate_test) <- unique_nos;
	duplicate_test <- duplicate_test[duplicate_test>1];
	dup_nos <- as.numeric(names(duplicate_test))
	dt <- 0;
	while (dt < length(dup_nos))	{
		dt <- dt+1;
		dupc_cases <- all_target_species[all_target_species$accepted_no %in% dup_nos[dt],];
		if (length(unique(c("uncoded","coded") %in% dupc_cases$status))==1)	{
			all_target_species[all_target_species$accepted_no %in% dup_nos[dt],] <- dupc_cases[dupc_cases$status=="coded",];
			} else	{
			print(dt)
			}
		}
	all_target_species <- unique(all_target_species);
	a_t_s <- nrow(all_target_species);
	}

# Output Coding Progress ####
write.csv(control_group_output,paste(total_group,"_Genera_Summary.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(target_control_taxa,paste("Target_",total_group,"_Genera_Summary.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(missing_species,paste("Unentered_",total_group,"_Species.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(missing_target_species,paste("Unentered_Target_",total_group,"_Species.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

write.csv(strict_target_control_taxa,paste("Target_",ingroup_title,"_Species.csv",sep=""),row.names = F,fileEncoding = "UTF-8");
write.csv(strict_uncoded_target_species,paste("Uncoded_",ingroup_title,"_Species_Information.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
#xlsx::write.xlsx(strict_uncoded_target_species,"Uncoded_Species_Information.xlsx");
writexl::write_xlsx(strict_uncoded_target_species,paste("Uncoded_",ingroup_title,"_Species_Information.xlsx",sep=""));
write.csv(coded_target_species,paste(ingroup_title,"_Coded_Species_Information.csv",sep=""),row.names = F,fileEncoding = "UTF-8");

write.csv(all_target_species,paste("All_Relevant_",total_group,"_Species_Info.csv",sep=""),row.names = F);
writexl::write_xlsx(all_target_species,paste("All_Relevant_",total_group,"_Species_Info.xlsx",sep=""));

# Illustrate Coding Progress ####
ingroup_title2 <- gsub("_"," ",ingroup_title);
nctaxa <- nrow(target_control_taxa);
mxspc <- max(target_control_taxa$nspc);
nsctaxa <- nrow(strict_target_control_taxa);
strict_target_control_taxa$subfamily[strict_target_control_taxa$subfamily==""] <- strict_target_control_taxa$family[strict_target_control_taxa$subfamily==""];
strict_target_control_taxa$subfamily <- gsub("idae","inae",strict_target_control_taxa$subfamily);
#strict_target_control_taxa$family[strict_target_control_taxa$subfamily=="Pseudotissotiinae"] <- "Pseudotissotiidae"
family_subfamily_combos <- strict_target_control_taxa[,colnames(strict_target_control_taxa) %in% c("subfamily","family")];
family_subfamily_combos <- unique(family_subfamily_combos);
family_subfamily_combos <- family_subfamily_combos[order(family_subfamily_combos$subfamily),];
family_subfamily_combos <- family_subfamily_combos[order(family_subfamily_combos$family),];
family_subfamily_combos$color <- rainbow(nrow(family_subfamily_combos)+2)[1:nrow(family_subfamily_combos)];
ingroup_families <- unique(family_subfamily_combos$family);
nfamilies <- length(ingroup_families)
nsubfamilies <- hist(match(family_subfamily_combos$family,unique(family_subfamily_combos$family)),breaks=0:nfamilies,plot=F)$count;

xsize <- 6.0;
ysize <- xsize*(3.5/9);
mxspc <- max(strict_target_control_taxa$nspc);
specify_basic_plot(mxx=nsctaxa,mnx=0,mxy=mxspc*(3.5/3),mny=0,abcissa=paste(ingroup_title2,"Genera Ranked"),ordinate="No. Species     ",xsize=xsize,ysize=ysize,font=franky);
specified_axis(axe=1,max_val=nsctaxa,min_val=0,maj_break = 10,med_break=5,min_break = 1,font=franky);
specified_axis(axe=2,max_val=mxspc,min_val=0,maj_break = 5,med_break=1,min_break = 1,font=franky);
for (ntp in 1:nsctaxa)	{
	taxon_color <- family_subfamily_combos$color[match(strict_target_control_taxa$subfamily[ntp],family_subfamily_combos$subfamily)]
	rect(ntp-1,0,ntp,strict_target_control_taxa$nspc[ntp],col=makeTransparent(taxon_color,100),lwd=0.5);
#	rect(ntp-1,0,ntp,strict_target_control_taxa$coded_spc[ntp],col=makeTransparent(taxon_color,100),lwd=0.5);
	rect(ntp-1,0,ntp,strict_target_control_taxa$coded_spc[ntp],col="black",lwd=0.5);
	text(ntp-0.5,strict_target_control_taxa$nspc[ntp]+1,labels=strict_target_control_taxa$genus[ntp],family=franky,cex=0.4,srt=45,pos=4,offset=0);
	if (strict_target_control_taxa$coded_spc[ntp]>=strict_target_control_taxa$nspc[ntp]/3)	{
		segments(ntp-1,strict_target_control_taxa$nspc[ntp]/3,ntp,strict_target_control_taxa$nspc[ntp]/3,col="white",lwd=2);
		} else	{
		segments(ntp-1,strict_target_control_taxa$nspc[ntp]/3,ntp,strict_target_control_taxa$nspc[ntp]/3,col="black",lwd=2);
		}
	}

if (sum(nsubfamilies)>5)	{
	taxon_break <- sum(nsubfamilies)/2;
	brk <- match(min(abs(cumsum(nsubfamilies)-taxon_break)),abs(cumsum(nsubfamilies)-taxon_break));
	init_families <- ingroup_families[1:brk];
	last_families <- ingroup_families[(brk+1):nfamilies];
	i_f <- 0.575*nsctaxa;
	i_sf <- 0.6*nsctaxa;
	l_f <- 0.825*nsctaxa;
	l_sf <- 0.85*nsctaxa;
	txn <- 1;
	yy <- mxspc*33/35;
	dy <- mxspc*1.5/35
	for (f1 in 1:length(init_families))	{
		scions <- sum(init_families[f1]==family_subfamily_combos$family);
		scion_names <- family_subfamily_combos$subfamily[family_subfamily_combos$family==init_families[f1]];
		scion_colors <- family_subfamily_combos$color[family_subfamily_combos$family==init_families[f1]];
		if (scions>1)	{
			text(i_f,yy,init_families[f1],pos=4,family=franky,cex=0.67);
			yy <- yy-dy;
			for (f2 in 1:scions)	{
#				points(l_f,yy,pch=22,bg=scion_colors[f2],cex=0.75)
				points(i_f,yy,pch=22,bg=makeTransparent(scion_colors[f2],100),cex=0.75);
				text(i_f,yy,paste(" ",scion_names[f2]),pos=4,family=franky,cex=0.67);
				yy <- yy-dy;
				txn <- txn+1;
				}
			} else	{
#			points(l_f,yy,pch=22,bg=family_subfamily_combos$color[txn],cex=0.75)
			points(i_f,yy,pch=22,bg=makeTransparent(family_subfamily_combos$color[txn],100),cex=0.75);
			text(i_f,yy,init_families[f1],pos=4,family=franky,cex=0.67);
			yy <- yy-dy;
			}
		}
	yy <- mxspc*33/35;
	dy <- mxspc*1.5/35;
	for (f1 in 1:length(last_families))	{
		scions <- sum(last_families[f1]==family_subfamily_combos$family);
		scion_names <- family_subfamily_combos$subfamily[family_subfamily_combos$family==last_families[f1]];
		scion_colors <- family_subfamily_combos$color[family_subfamily_combos$family==last_families[f1]];
		if (scions>1)	{
			text(l_f,yy,last_families[f1],pos=4,family=franky,cex=0.67);
			yy <- yy-dy;
			for (f2 in 1:scions)	{
#				points(l_f,yy,pch=22,bg=scion_colors[f2],cex=0.75)
				points(l_f,yy,pch=22,bg=makeTransparent(scion_colors[f2],100),cex=0.75);
				text(l_f,yy,paste(" ",scion_names[f2]),pos=4,family=franky,cex=0.67);
				yy <- yy-dy;
				txn <- txn+1;
				}
			} else	{
#			points(l_f,yy,pch=22,bg=family_subfamily_combos$color[txn],cex=0.75)
			points(l_f,yy,pch=22,bg=makeTransparent(family_subfamily_combos$color[txn],100),cex=0.75);
			text(l_f,yy,last_families[f1],pos=4,family=franky,cex=0.67);
			yy <- yy-dy;
			}
		}
	}

all_relv_target_species <- all_target_species[all_target_species$fa!="",];
a_t_s <- nrow(all_relv_target_species);
all_relv_target_species$fa_lb <- time_scale$ma_lb[match(all_relv_target_species$fa,time_scale$interval)];
all_relv_target_species$fa_ub <- time_scale$ma_ub[match(all_relv_target_species$fa,time_scale$interval)];
study_stage_scale <- study_stage_scale[order(-study_stage_scale$ma_lb),];
nstages <- nrow(study_stage_scale);
stage_summaries <- data.frame(stage=as.character(study_stage_scale$interval),
															all_species=as.numeric(rep(0,nstages)),
															coded_species=as.numeric(rep(0,nstages)));
for (ss in 1:nstages)	{
	stage_species <- all_relv_target_species[all_relv_target_species$fa_lb>study_stage_scale$ma_ub[ss] & all_relv_target_species$fa_ub<study_stage_scale$ma_lb[ss],]
	stage_summaries$all_species[ss] <- nrow(stage_species);
	stage_summaries$coded_species[ss] <- sum(stage_species$status=="coded");
	}

mxy <- ceiling(max(stage_summaries$all_species/10))*10;	# set maximum y value
mny <- 0;							# set maximum y value
ybreaks <- as.numeric(set_axis_breaks_new(max_no=mxy));
ysize <- 4*4.285714285/6;
xsize <- 4.5;
use_strat_labels <- T;						# if T, then strat_names will be plotted on X-axis inside boxes
strat_label_size <- 2;					# size of the labels for chronostratigraphic units
alt_back <- F;								# if T, then the background will alternat shades between major intervals
plot_title <- ""			# Name of the plot; enter "" for nothing
ordinate <- "First Sampled Species";							# Label of Y-axis
hues <- T;								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
strat_colors <- study_stage_scale$color;
strat_names <- study_stage_scale$st;
strat_names[strat_names=="Pra"] <- "Pr";
onset <- -max(abs(study_stage_scale$ma_lb));
end <- -min(abs(study_stage_scale$ma_ub));
time_scale_to_plot <- c(-study_stage_scale$ma_lb,end);
yearbreaks <- as.numeric(set_axis_breaks_new(max_no=abs(onset),min_no=abs(end)));
Phanerozoic_Timescale_Plot_Flexible(onset,end,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky);
specified_axis(axe=2,max_val = mxy,min_val=0,maj_break = ybreaks[1],med_break=ybreaks[2],min_break=ybreaks[3],orient=2,font=franky);

for (ss in 1:nrow(study_stage_scale))	{
	rect(-study_stage_scale$ma_lb[ss],0,-study_stage_scale$ma_ub[ss],stage_summaries$all_species[ss],col=makeTransparent(study_stage_scale$color[ss],100));
	rect(-study_stage_scale$ma_lb[ss],0,-study_stage_scale$ma_ub[ss],stage_summaries$coded_species[ss],col="black");
	segments(-study_stage_scale$ma_lb[ss],stage_summaries$all_species[ss]/3,-study_stage_scale$ma_ub[ss],stage_summaries$all_species[ss]/3,lty=2,col="purple",lwd=2);
	}

nslices <- nrow(stage_slice_scale_2);
stage_slice_summaries <- data.frame(stage_slice=as.character(stage_slice_scale_2$interval),
																		all_species=as.numeric(rep(0,nslices)),
																		coded_species=as.numeric(rep(0,nslices)));
for (ss in 1:nslices)	{
	stage_species <- all_relv_target_species[all_relv_target_species$fa_lb>stage_slice_scale_2$ma_ub[ss] & all_relv_target_species$fa_ub<stage_slice_scale_2$ma_lb[ss],]
	stage_slice_summaries$all_species[ss] <- nrow(stage_species);
	stage_slice_summaries$coded_species[ss] <- sum(stage_species$status=="coded");
	}

mxy <- ceiling(max(stage_slice_summaries$all_species/5))*5;	# set maximum y value
mny <- 0;							# set maximum y value
ybreaks <- as.numeric(set_axis_breaks_new(max_no=mxy));
Phanerozoic_Timescale_Plot_Flexible(onset,end,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky);
specified_axis(axe=2,max_val = mxy,min_val=0,maj_break = ybreaks[1],med_break=ybreaks[2],min_break=ybreaks[3],orient=2,font=franky);
for (ss in 1:nslices)	{
	rect(-stage_slice_scale_2$ma_lb[ss],0,-stage_slice_scale_2$ma_ub[ss],stage_slice_summaries$all_species[ss],col=makeTransparent(stage_slice_scale_2$color[ss],100));
	rect(-stage_slice_scale_2$ma_lb[ss],0,-stage_slice_scale_2$ma_ub[ss],stage_slice_summaries$coded_species[ss],col="black");
	segments(-stage_slice_scale_2$ma_lb[ss],stage_slice_summaries$all_species[ss]/3,-stage_slice_scale_2$ma_ub[ss],stage_slice_summaries$all_species[ss]/3,lty=1,col="purple",lwd=2);
	}

all_relv_target_species$region <- gsub("East Laurentia","Eastern Americas",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Midcontinent North America","Eastern Americas",all_relv_target_species$region);
all_relv_target_species$region <- gsub("West Laurentia","Cordilleran",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Western Canada","Cordilleran",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Victoria","Tasman",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Australia","Tasman",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Japan","South China",all_relv_target_species$region);
all_relv_target_species$region <- gsub("North China","South China",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Alborz","Uralian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("East Baltica","Uralian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Mongolia","Uralian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Siberian Platform","Uralian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Mongolia","Uralian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("Midland Valley","Rhenish-Bohemian",all_relv_target_species$region);
all_relv_target_species$region <- gsub("East Avalonia","Rhenish-Bohemian",all_relv_target_species$region);
study_sites$region <- gsub("East Avalonia","Rhenish-Bohemian",study_sites$region);

all_relv_target_species <- all_relv_target_species[!all_relv_target_species$species %in% excluded_species$species,];
a_t_s <- nrow(all_relv_target_species);
relv_regions <- relv_plates <- c();
for (sp in 1:a_t_s)	{
	relv_regions <- unique(c(relv_regions,strsplit(all_relv_target_species$region[sp],", ")[[1]]));
	all_relv_target_species$geoplates[sp] <- gsub(", ",",",all_relv_target_species$geoplates[sp]);
	all_geoplates <- strsplit(all_relv_target_species$geoplates[sp],",")[[1]]
	if (!is.na(as.numeric(all_geoplates[1])))	{
		all_geoplates <- as.numeric(all_geoplates)
		relv_plates <- unique(c(relv_plates,all_geoplates));
		} else	{
		relv_plates <- unique(c(relv_plates,as.numeric(strsplit(all_relv_target_species$geoplates[sp],",")[[1]])));
		}
	}
#relv_plates <- sort(relv_plates);
relv_plates <- sort(unique(study_sites$geoplate));
relv_regions <- unique(study_sites$region[study_sites$region!=""])
#relv_regions[relv_regions!=""];

nregions <- length(relv_regions);
nplates <- length(relv_plates);
species_by_region <- data.frame(region=as.character(relv_regions),all_species=rep(0,nregions),coded_species=rep(0,nregions));
species_by_plate <- data.frame(plate=as.numeric(relv_plates),all_species=rep(0,nplates),coded_species=rep(0,nplates));
name_cols <- match(c("taxon_name","accepted_name","alt_combo1","alt_combo2","alt_combo3","alt_combo4","alt_combo5","alt_combo6"),colnames(control_species));
uglies <- c();
# set aside separate information for the first stage in which species appear
first_stage_info <- all_relv_target_species;
age <- first_stage_info$fa_lb;
first_stage_info$first_stage <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",fine_time_scale=all_stage_scale);
tspecies <- nrow(control_species);
for (sp in 1:a_t_s)	{
#	spc_rows <- unique(which(control_species==all_relv_target_species$species[sp],arr.ind=T)[,1]);
	spc_rows <- (1:tspecies)[control_species$taxon_no==all_relv_target_species$accepted_no[sp]];
	if (length(spc_rows)>0 && length(unique(control_species$accepted_rank[spc_rows]))==2)	{
		if (is.subspecies(all_relv_target_species$species[sp]))	{
			spc_rows <- spc_rows[control_species$accepted_rank[spc_rows]=="subspecies"];
			} else	{
			spc_rows <- spc_rows[control_species$accepted_rank[spc_rows]=="species"];
			}
		}
#	control_species[spc_rows,]
	other_combos <- possible_names <- c();
	if (length(spc_rows)==0)	{
		epithet <- divido_species_epithets(all_relv_target_species$species[sp]);
		poss_spc_rows <- unique(which(control_species==epithet,arr.ind=T)[,1]);
		control_species[poss_spc_rows,]
		epithet_atomized <- strsplit(epithet,"")[[1]];
		l_e_a <- length(epithet_atomized);
		this_genus <- divido_genus_names_from_species_names(all_relv_target_species$species[sp]);
		if (epithet_atomized[l_e_a]=="a")	{
			alt_combos <- c(paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"um"),collapse="")),paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"us"),collapse="")));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[1],arr.ind=T)[,1]));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[2],arr.ind=T)[,1]));
			possible_names <- alt_combos[alt_combos %in% control_species[spc_rows,]];
			} else if (epithet_atomized[l_e_a]=="m" && epithet_atomized[l_e_a-1]=="u")	{
			alt_combos <- c(paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"a"),collapse="")),paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"us"),collapse="")));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[1],arr.ind=T)[,1]));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[2],arr.ind=T)[,1]));
			possible_names <- alt_combos[alt_combos %in% control_species[spc_rows,]];
			} else if (epithet_atomized[l_e_a]=="s" && epithet_atomized[l_e_a-1]=="u")	{
			alt_combos <- c(paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"a"),collapse="")),paste(this_genus,paste(c(epithet_atomized[1:(l_e_a-1)],"um"),collapse="")));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[1],arr.ind=T)[,1]));
			spc_rows <- c(spc_rows,unique(which(control_species==alt_combos[2],arr.ind=T)[,1]));
			possible_names <- alt_combos[alt_combos %in% control_species[spc_rows,]];
			}
		if (length(possible_names)==0)	{
			genus_subgenus <- divido_subgenus_names_from_genus_names(this_genus);
			if (genus_subgenus[2]!="")	{
				spc_rows <- unique(c(spc_rows,unique(which(control_species==paste(genus_subgenus[1],epithet),arr.ind = T)[,1])));
				spc_rows <- unique(c(spc_rows,unique(which(control_species==paste(genus_subgenus[2],epithet),arr.ind = T)[,1])));
				spc_rows <- unique(c(spc_rows,unique(which(control_species==paste(genus_subgenus[1]," (",genus_subgenus[1],") ",epithet,sep=""),arr.ind = T)[,1])));
				spc_rows <- unique(c(spc_rows,unique(which(control_species==paste(genus_subgenus[2]," (",genus_subgenus[2],") ",epithet,sep=""),arr.ind = T)[,1])));
				}
			}
		}
#		control_species[spc_rows,]
	if (length(spc_rows)==0)	uglies <- c(uglies,sp);
	sr <- 0;
#	possible_names <- c();
	while(sr < length(spc_rows))	{
		sr <- sr+1;
		for (i in 1:length(name_cols))	possible_names <- unique(c(possible_names,control_species[spc_rows[sr],name_cols[i]]));
		possible_names <- unique(c(possible_names,all_relv_target_species$species[sp]),other_combos);
		}
	
#	control_species$taxon_name[proetide_taxonomy$accepted_name %in% proetide_taxonomy$accepted_name[proetide_taxonomy$taxon_name %in% all_relv_target_species$species[sp]]];
#	control_species$taxon_name[proetide_taxonomy$accepted_name %in% proetide_taxonomy$accepted_name[proetide_taxonomy$taxon_name %in% all_relv_target_species$species[sp]]];
	if (length(possible_names)>0)	{
		taxon_sites <- study_sites[study_sites$collection_no %in% control_finds$collection_no[control_finds$accepted_no %in% all_relv_target_species$accepted_no[sp]],];
		if (nrow(taxon_sites)==0)
			taxon_sites <- study_sites[study_sites$collection_no %in% control_finds$collection_no[control_finds$accepted_name %in% possible_names],];
		if (nrow(taxon_sites)==0)
			taxon_sites <- study_sites[study_sites$collection_no %in% control_finds$collection_no[control_finds$accepted_name_orig %in% possible_names],];
		if (nrow(taxon_sites)==0)	{
			taxon_sites <- trilo_sites[trilo_sites$collection_no %in% control_finds$collection_no[control_finds$accepted_name %in% possible_names],];
			if (nrow(taxon_sites)>0)	{
				taxon_sites$realm <- rock_database$realm[match(taxon_sites$rock_no_sr,rock_database$rock_no)];
				taxon_sites$region <- rock_database$region[match(taxon_sites$rock_no_sr,rock_database$rock_no)];
				taxon_sites$region[taxon_sites$region==""] <- taxon_sites$realm[taxon_sites$region==""];
				}
			}
		if (nrow(taxon_sites)==0)	uglies <- c(uglies,sp);

		oldest_sites <- taxon_sites[taxon_sites$ma_lb>=max(taxon_sites$ma_ub),];
		spc_regions <- unique(oldest_sites$region);
		spc_regions <- spc_regions[spc_regions!=""];
		species_by_region$all_species[match(spc_regions,species_by_region$region)] <- species_by_region$all_species[match(spc_regions,species_by_region$region)]+1;
#	all_relv_target_species$geoplates[sp] <- gsub(", ",",",all_relv_target_species$geoplates[sp]);
		spc_plates <- unique(oldest_sites$geoplate);
		species_by_plate$all_species[match(spc_plates,species_by_plate$plate)] <- species_by_plate$all_species[match(spc_plates,species_by_plate$plate)]+1;
		if (all_relv_target_species$status[sp]=="coded")	{
			species_by_region$coded_species[match(spc_regions,species_by_region$region)] <- species_by_region$coded_species[match(spc_regions,species_by_region$region)]+1;
			species_by_plate$coded_species[match(spc_plates,species_by_plate$plate)] <- species_by_plate$coded_species[match(spc_plates,species_by_plate$plate)]+1;
			}
		
		first_stage <- study_stage_scale$interval[sum(study_stage_scale$ma_lb>=max(oldest_sites$ma_lb))];
		first_stage_no <- sum(study_stage_scale$ma_lb>=max(oldest_sites$ma_lb));
		# Start Here. Collect info on sampling from first stage only
		first_stage_sites <- taxon_sites[taxon_sites$ma_lb>study_stage_scale$ma_ub[first_stage_no],];
		if (nrow(first_stage_sites)==0)	{
			first_stage_sites <- rbind(oldest_sites,oldest_sites);
			if (nrow(first_stage_sites)>1)	{
				interval_lbs <- first_stage_sites $interval_lb;
				for (i in 0:9)	interval_lbs <- gsub(i,"",interval_lbs);
				if (length(unique(interval_lbs))>1)	{
					uglier <- data.frame(interval=interval_lbs,ma_lb=first_stage_sites$ma_lb,ma_ub=first_stage_sites$ma_ub)
					uglier <- uglier[order(-uglier$ma_lb),];
					unique_intervals <- unique(uglier$interval);
					cutoff <- uglier$ma_lb[match(unique_intervals[2],uglier$interval)];
					first_stage_sites <- first_stage_sites[first_stage_sites$ma_lb>cutoff,];
					}
				}
			}
		first_stage_info$nfinds[sp] <- nrow(first_stage_sites);
		first_stage_info$nrocks[sp] <- max(1,length(unique(first_stage_sites$rock_no_sr[first_stage_sites$rock_no_sr!=0])));
		first_stage_info$geoplates[sp] <- paste(unique(first_stage_sites$geoplate),collapse=", ");
		first_stage_sites$region[first_stage_sites$region==""] <- first_stage_sites$realm[first_stage_sites$region==""];
		first_stage_sites$region[first_stage_sites$region=="All"] <- "";
		first_stage_info$region[sp] <- paste(unique(first_stage_sites$region[first_stage_sites$region!=""]),collapse=",");
		}
	}

xsize <- 6;
ysize <- xsize*(3.5/xsize);
mxspc <- ceiling(max(species_by_plate$all_species)/10)*10;
specify_basic_plot(mxx=nplates,mnx=0,mxy=mxspc*(3.5/3),mny=0,abcissa="GPlate",
									 ordinate="No. Species     ",xsize=xsize,ysize=ysize,font=franky);
xbreaks <- as.numeric(set_axis_breaks_new(nplates));
xbreaks[xbreaks<1] <- 1;
specified_axis_w_labels(axe=1,max_val=nplates,min_val=0,maj_break=1,med_break=1,min_break=1,axis_labels=as.character(species_by_plate$plate),label_pos = "mid",font=franky,font_size=0.5,orient=2);
ybreaks <- as.numeric(set_axis_breaks_new(mxspc));
specified_axis(axe=2,max_val=mxspc,min_val=0,maj_break=ybreaks[1],med_break=ybreaks[2],min_break=ybreaks[3],font=franky);
continent_colors <- rainbow(9);
for (np in 1:nplates)	{
	plate_color <- continent_colors[floor(species_by_plate$plate[np]/100)]
	rect(np-1,0,np,species_by_plate$all_species[np],col=makeTransparent(plate_color,100),lwd=0.5);
	rect(np-1,0,np,species_by_plate$coded_species[np],col="black",lwd=0.5);
	if (species_by_plate$coded_species[np]>species_by_plate$all_species[np]/3){
		segments(np-1,species_by_plate$all_species[np]/3,np,species_by_plate$all_species[np]/3,lty=3,col="white",lwd=2);
		} else	{
		segments(np-1,species_by_plate$all_species[np]/3,np,species_by_plate$all_species[np]/3,lty=3,col="black",lwd=2);
		}
	}

# get region : realm key from somewhere!!!
species_by_region <- cbind(realm=as.character(biogeography_key$realm[match(species_by_region$region,biogeography_key$biogeographic_unit)]),species_by_region);
species_by_region <- species_by_region[order(species_by_region$realm),]
region_colors <- rainbow(nregions+1);

xsize <- 6.0;
ysize <- xsize*(3.5/9);
mxspc <- ceiling(max(species_by_region$all_species)/10)*10;
specify_basic_plot(mxx=nregions,mnx=0,mxy=mxspc*(3.5/3),mny=0,
									 abcissa="",
									 ordinate="No. Species     ",xsize=xsize,ysize=ysize,font=franky);
specified_axis_w_labels(axe=1,max_val=nregions,min_val=0,maj_break=1,med_break=1,min_break=1,axis_labels=species_by_region$region,orient=2,label_pos = "mid",font=franky,font_size =0.75);
ybreaks <- as.numeric(set_axis_breaks_new(mxspc));
specified_axis(axe=2,max_val=mxspc,min_val=0,maj_break=ybreaks[1],med_break=ybreaks[2],min_break=ybreaks[3],font=franky);

for (nr in 1:nregions)	{
	rect(nr-1,0,nr,species_by_region$all_species[nr],col=makeTransparent(region_colors[nr],100));
	rect(nr-1,0,nr,species_by_region$coded_species[nr],col="black");
	if (species_by_region$coded_species[nr]>species_by_region$all_species[nr]/3){
		segments(nr-1,species_by_region$all_species[nr]/3,nr,species_by_region$all_species[nr]/3,lty=3,col="white",lwd=2);
		} else	{
		segments(nr-1,species_by_region$all_species[nr]/3,nr,species_by_region$all_species[nr]/3,lty=3,col="black",lwd=2);
		}
	}

write.csv(proetoid_output,"Proetoid_Genera_Summary.csv",row.names=F,fileEncoding = "UTF-8");

all_target_species$ma_lb <- stage_slice_scale$ma_lb[match(all_target_species$fa,stage_slice_scale$interval)];
all_target_species$ma_ub <- stage_slice_scale$ma_ub[match(all_target_species$la,stage_slice_scale$interval)];
all_target_species$ma_ub[is.na(all_target_species$ma_ub) & !is.na(all_target_species$ma_lb)] <- min(stage_slice_scale$ma_ub);
dataless_species <- all_target_species[is.na(all_target_species$ma_lb),];
all_target_species <- all_target_species[!is.na(all_target_species$ma_lb),];
ngenera_st <- vector(length=nstages);
mxspc <- 0;
for (st in 1:nstages)	{
	stage_lb <- study_stage_scale$ma_lb[st];
	stage_ub <- study_stage_scale$ma_ub[st];
	stage_sites <- study_sites[study_sites$ma_ub<stage_lb & study_sites$ma_lb>stage_ub,]
	stage_finds <- control_finds[control_finds$collection_no %in% stage_sites$collection_no,];
	stage_target_species <- first_stage_info[first_stage_info$first_stage %in% study_stage_scale$interval[st],];
#	stage_target_species$species[stage_target_species$nrocks==0] <- 1;
	stage_target_species$sampling <- (rank(stage_target_species$nfinds)+rank(stage_target_species$nrocks))/2;
	stage_target_species <- stage_target_species[order(stage_target_species$genus_no,-stage_target_species$sampling),];
	nspecies_st <- nrow(stage_target_species);
	ngenera_st <- length(unique(stage_target_species$genus_no));
	genus_names_st <- unique(stage_target_species$genus);
	
	# now, sort these by species richness, then sampling.
	genus_richness <- hist(match(stage_target_species$genus_no,sort(unique(stage_target_species$genus_no))),breaks=0:ngenera_st,plot=F)$counts;
	if (mxspc < max(genus_richness))	mxspc <- max(genus_richness);
	genus_finds <- genus_rocks <- genus_sampling <- vector(length=ngenera_st);
	for (gn in 1:ngenera_st)	{
		genus_aliases <- c(genus_names_st[gn],divido_subgenus_names_from_genus_names(genus_names_st[gn])[2]);
		genus_aliases <- genus_aliases[genus_aliases!=""];
		genus_sites <- study_sites[study_sites$collection_no %in% unique(stage_finds$collection_no[stage_finds$genus %in% genus_aliases]),];
		genus_finds[gn] <- nrow(genus_sites);
		genus_forms <- unique(c(genus_sites$rock_no_sr[genus_sites$rock_no_sr>0],genus_sites$rock2_no_sr[genus_sites$rock2_no_sr>0]));
		if (length(genus_forms)>0)	{
			genus_rocks[gn] <- length(genus_forms);
			} else	{
			genus_rocks[gn] <- length(unique(genus_sites$geoplate));
			}
		}
	genus_sampling <- (rank(genus_finds)+rank(genus_rocks))/2;
	sorty_mcsortface <- data.frame(genus_no=unique(stage_target_species$genus_no),genus_richness=genus_richness,genus_finds=genus_finds,genus_rocks=genus_rocks,genus_sampling=genus_sampling);
	sorty_mcsortface <- sorty_mcsortface[order(-sorty_mcsortface$genus_richness,-sorty_mcsortface$genus_sampling,-sorty_mcsortface$genus_rocks),]
	stage_target_species$genus_rank <- match(stage_target_species$genus_no,sorty_mcsortface$genus_no);
	stage_target_species <- stage_target_species[order(stage_target_species$genus_rank,-stage_target_species$sampling),]
	output_csv <- paste(study_stage_scale$interval[st],"_Target_Species.csv",sep="");
	output_xcl <- paste(study_stage_scale$interval[st],"_Target_Species.xlsx",sep="");
	write.csv(stage_target_species,output_csv,row.names = F,fileEncoding = "UTF-8");
	writexl::write_xlsx(stage_target_species,output_xcl,);
	}

for (st in 1:nstages)  {
	stage_lb <- study_stage_scale$ma_lb[st];
	stage_ub <- study_stage_scale$ma_ub[st];
	stage_sites <- study_sites[study_sites$ma_ub<stage_lb & study_sites$ma_lb>stage_ub,]
	stage_finds <- control_finds[control_finds$collection_no %in% stage_sites$collection_no,];
	stage_target_species <- all_target_species[all_target_species$ma_lb<=stage_lb & all_target_species$ma_lb>stage_ub,];
	nspecies_st <- nrow(stage_target_species);
	for (sp in 1:nspecies_st)	{
		
		}
	multiplates <- (1:nspecies_st)[!stage_target_species$geoplates==gsub(",","",stage_target_species$geoplates)];
	mp <- 0;
	while (mp < length(multiplates))	{
		mp <- mp+1;
		species_sites <- study_sites[study_sites$collection_no %in% unique(control_finds$collection_no[which(control_finds==stage_target_species$species[multiplates[mp]],arr.ind=T)[,1]]),]
		oldest_sites <- species_sites[species_sites$ma_lb>=max(species_sites$ma_ub),];
		stage_target_species$geoplates[multiplates[mp]] <- paste(unique(oldest_sites$geoplate),collapse=",");
		}
	multiregions <- (1:nspecies_st)[!stage_target_species$region==gsub(",","",stage_target_species$region)];
	mr <- 0;
	while (mr < length(multiregions))	{
		mr <- mr+1;
		species_sites <- study_sites[study_sites$collection_no %in% unique(control_finds$collection_no[which(control_finds==stage_target_species$species[multiregions[mr]],arr.ind=T)[,1]]),]
		oldest_sites <- species_sites[species_sites$ma_lb>=max(species_sites$ma_ub),];
		stage_target_species$region[multiregions[mr]] <- paste(unique(oldest_sites$region),collapse=", ");
		}
	noregions <- (1:nspecies_st)[stage_target_species$region==""];
	nr <- 0;
	while (nr < length(noregions))	{
		nr <- nr+1;
		species_sites <- study_sites[study_sites$collection_no %in% unique(control_finds$collection_no[which(control_finds==stage_target_species$species[noregions[nr]],arr.ind=T)[,1]]),]
		oldest_sites <- species_sites[species_sites$ma_lb>=max(species_sites$ma_ub),];
		stage_target_species$region[noregions[nr]] <- paste(unique(oldest_sites$region),collapse=", ");
		}
	stage_genera <- unique(stage_target_species$genus);
	ngenera_st[st] <- length(stage_genera);
	stage_genus_species <- hist(match(stage_target_species$genus,stage_genera),breaks=0:ngenera_st[st],plot=F)$counts;
	stage_target_genus_nos <- unique(stage_target_species$genus_no);
	names(stage_genus_species) <- stage_target_genus_nos;
	genus_finds <- genus_rocks <- vector(length=length(stage_genera));
	species_sampling <- vector(length=nspecies_st);
	for (g in 1:length(stage_genera))	{
#		(rank(stage_target_species$nfinds[stage_target_species$genus==stage_genera[g]])+rank(stage_target_species$nrocks[stage_target_species$genus==stage_genera[g]]))/2;
		genus_finds[g] <- sum(stage_target_species$nfinds[stage_target_species$genus==stage_genera[g]]);
		genus_rocks[g] <- sum(stage_target_species$nrocks[stage_target_species$genus==stage_genera[g]]);
		genus_species <- stage_target_species[stage_target_species$genus==stage_genera[g],];
		}
	genus_sampling <- (rank(genus_finds)+rank(genus_rocks))/2;
	stage_target_genus_nos <- stage_target_genus_nos[order(-stage_genus_species,-genus_sampling)];
	stage_target_species <- stage_target_species[order(match(stage_target_species$genus_no,stage_target_genus_nos)),];
	output_csv <- paste(study_stage_scale$interval[st],"_Target_Species.csv",sep="");
	output_xcl <- paste(study_stage_scale$interval[st],"_Target_Species.xlsx",sep="");
	write.csv(stage_target_species,output_csv,row.names = F,fileEncoding = "UTF-8");
	writexl::write_xlsx(stage_target_species,output_xcl,);
  # order genera by # species & output them
	if (mxspc < max(stage_genus_species))
		mxspc <- max(stage_genus_species);
  }

xsize <- 6;
ysize <- 3.5;
mx_gen <- max(ngenera_st);
for (st in 1:nstages)  {
  stage_lb <- study_stage_scale$ma_lb[st];
  stage_ub <- study_stage_scale$ma_ub[st];
  stage_target_species <- all_target_species[all_target_species$ma_lb<=stage_lb & all_target_species$ma_lb>stage_ub,];
  stage_genera <- unique(stage_target_species$genus);
  ngenera_st[st] <- length(stage_genera);
  stage_genus_species <- hist(match(stage_target_species$genus,stage_genera),breaks=0:ngenera_st[st],plot=F)$counts;
  names(stage_genus_species) <- stage_genera;
  coded_stage_target_species <- stage_target_species[stage_target_species$status=="coded",];
  stage_genus_species_coded <- stage_genus_species;
  stage_genus_species_coded[stage_genus_species_coded>0] <- 0;
  if (nrow(coded_stage_target_species)>0) {
    stage_genus_species_coded <- hist(match(coded_stage_target_species$genus,stage_genera),breaks=0:ngenera_st[st],plot=F)$counts;
    names(stage_genus_species_coded) <- stage_genera;
    }
  stage_genus_species <- sort(stage_genus_species,decreasing = T);
  stage_genus_species_coded <- stage_genus_species_coded[match(names(stage_genus_species),names(stage_genus_species_coded))]
#  mxspc <- max(stage_genus_species);
  stage_families <- unique(proetide_taxonomy$family[match(stage_genera,proetide_taxonomy$taxon_name)]);
  if ("Phillipsiidae" %in% stage_families)  {
    abcissa <- paste(study_stage_scale$interval[st],"Proetid and Phillipsid Genera Ranked");
    } else  {
    abcissa <- paste(study_stage_scale$interval[st],"Proetid Genera Ranked");
    }
  specify_basic_plot(mxx=ngenera_st[st],mnx=0,mxy=mxspc*ysize/3,mny=0,
  									 abcissa=abcissa,ordinate="No. Species",
  									 xsize=xsize*ngenera_st[st]/mx_gen,ysize=ysize,font=franky);
  specified_axis(axe=1,max_val=ngenera_st[st],min_val=0,maj_break = 10,med_break=5,min_break = 1,font=franky);
  specified_axis(axe=2,max_val=mxspc,min_val=0,maj_break = 5,med_break=1,min_break = 1,font=franky);
  for (ng in 1:ngenera_st[st])  {
    rect(ng-1,0,ng,stage_genus_species[ng],col=study_stage_scale$col[st]);
    if (stage_genus_species_coded[ng]>0)
      rect(ng-1,0,ng,stage_genus_species_coded[ng],col="black");
		text(ng-0.5,stage_genus_species[ng]+1,labels=names(stage_genus_species)[ng],family=franky,cex=0.5,srt=45,pos=4,offset=0);
    segments(ng-1,stage_genus_species[ng]/3,ng,stage_genus_species[ng]/3,col="purple",lwd=2);
  	}
  }

# 'While I'm at it' data checking, etc., routines.... ####
# Look for post-Hastarian proetids
new_ub <- time_scale$ma_ub[match("Pennsylvanian",time_scale$interval)];
later_sites <- trilo_sites[trilo_sites$ma_ub<myr_ub & trilo_sites$ma_lb>new_ub,];
later_sites$realm <- rock_database$realm[match(later_sites$rock_no_sr,rock_database$rock_no)]
later_sites$region <- rock_database$region[match(later_sites$rock_no_sr,rock_database$rock_no)]
later_sites$region[later_sites$region==""] <- later_sites$realm[later_sites$region==""]
later_finds <- control_finds[control_finds$collection_no %in% later_sites$collection_no,];
last_proetid_finds <- later_finds[later_finds$family %in% "Proetidae",];
last_proetid_finds <- last_proetid_finds[last_proetid_finds$identified_rank %in% c("species","subspecies"),];
later_proetid_genera <- sort(unique(last_proetid_finds$genus));
later_proetid_species <- sort(unique(last_proetid_finds$identified_name));
nn <- unique(data.frame(species=last_proetid_finds$identified_name[last_proetid_finds$identified_name %in% later_proetid_species], 
												accepted_species=last_proetid_finds$accepted_name[last_proetid_finds$identified_name %in% later_proetid_species], 
												genus=last_proetid_finds$genus[last_proetid_finds$identified_name %in% later_proetid_species]));
nn$subfamily <- proetide_taxonomy$subfamily[match(nn$genus,proetide_taxonomy$taxon_name)];
nn$family <- ingroup_taxonomy$family[match(nn$genus,ingroup_taxonomy$taxon_name)];
nn$family <- proetide_taxonomy$family[match(nn$subfamily,proetide_taxonomy$taxon_name)]


later_sites$paleolng[1];
later_sites$paleolat[1];
reconstruct(c(later_sites$lng[1],later_sites$lat[1]),round(later_sites$ma_ub[1],0));
reconstruct(c(later_sites$lng[1],later_sites$lat[1]),round(later_sites$ma_ub[1],0));
xx <- read.table("https://gws.gplates.org/reconstruct/reconstruct_points/?points=-98.1,31&time=352&model=MERDITH2021")
xx$V6
xx$V7
ingroup_taxonomy[ingroup_taxonomy$taxon_name %in% c("Liobole","Liobole (Panibole)","Panibole"),]
lat_lng_age <- site_lng_lat_age[1,]

#write.csv(ingroup_taxonomy[ingroup_taxonomy$taxon_no %in% 431358,],"B_spinosus.csv",row.names=F,fileEncoding="UTF-8");
#control_finds[control_finds$taxon_name %in% "Xenocybe",]
#ingroup_taxonomy[ingroup_taxonomy$accepted_name %in% "Xenocybe",]
#control_finds[control_finds$accepted_name %in% "Brachymetopus maccoyi spinosus",]

control_finds <- control_finds[control_finds$identified_rank %in% c("species","subspecies"),];
authorizers <- unique(control_finds$authorizer);
authorizer_counts <- hist(match(control_finds$authorizer,authorizers),breaks=(0:length(authorizers)),plot=F)$counts
names(authorizer_counts) <- authorizers;
authorizer_counts <- -sort(-authorizer_counts);
reference_nos <- sort(unique(control_finds$reference_no))
reference_no_counts <- hist(control_finds$reference_no,breaks=(0:max(control_finds$reference_no)),plot=F)$counts
names(reference_no_counts) <- 1:max(reference_nos);
reference_no_counts <- reference_no_counts[reference_no_counts>0];
reference_no_counts <- sort(reference_no_counts,decreasing = T);
pbdb_alexandria$pbdb_references$distinct_citation[pbdb_alexandria$pbdb_references$reference_no %in% names(reference_no_counts)[1:11]];
main_ref_nos <- as.numeric(names(reference_no_counts)[1:10]);

# The End ####
write.csv(pbdb_taxonomy[pbdb_taxonomy$parent_name=="Acanthoceratoidea" & pbdb_taxonomy$taxon_rank=="family"& pbdb_taxonomy$difference=="",],"Acanthoceratoidea.csv",row.names=F);
write.csv(pbdb_taxonomy[pbdb_taxonomy$taxon_name=="Pseudotissotiinae",],"Pseudotissotiinae.csv",row.names = F)

Kamerunoceras_turoniense <- pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_name %in% "Kamerunoceras turoniense"],];

Kamerunoceras_turoniense$collection_no[4:5]
{}

