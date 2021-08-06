#### Setup Program ####
source('~/Documents/R_Projects/Common_R_Source_Files/Chronos.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Historical_Diversity_Metrics.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Occurrence_Data_Routines.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Sampling_and_Occupancy_Distributions.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Stratigraphy.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r'); 		#

taxon_fields <- c("phylum","phylum_no","class","class_no","order","order_no","suborder","suborder_no","superfamily","superfamily_no","family","family_no","subfamily","subfamily_no","genus","genus_no","subgenus","subgenus_no")
prior_clam_directory <- "~/Documents/R_Projects/Prior_Clams/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";	# directory to folder where you keep common source
setwd(prior_clam_directory);
prec <- 0.1;
dummy_parent <- accersi_taxonomic_data_for_one_taxon_no(1);
parent_nos <- c("orig_no","taxon_no","accepted_no","parent_no","immpar_no","ref_pubyr","reference_no","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma",taxonomic_field[!taxonomic_field %in% taxonomic_rank],"type_taxon_no","n_orders","n_families","n_genera","n_species");
dummy_parent[,colnames(dummy_parent) %in% parent_nos] <- 0;
dummy_parent[,!colnames(dummy_parent) %in% parent_nos] <- "";

load("~/Documents/R_Projects/Data_for_R/Gradstein_2012_Augmented.RData"); # refined Gradstein 2012 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/Rock_Unit_Database.RData"); # refined Gradstein 2012 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/PaleoDB_Edits.RData"); # refined Gradstein 2012 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData"); # refined Gradstein 2012 timescale & biozonations

# load data #
time_scale_stratigraphic_scale <- "Stage Slice";
temporal_precision <- 0.05;
rock_unit_databases <- rock_unit_data;						# information about rock ages, biozonations and sequence stratigraphy
rock_database <- rock_unit_databases$rock_unit_database;
rock_to_zone_database <- rock_unit_databases$rock_to_zone_database;
chronostratigraphic_databases <- gradstein_2012_emended;	# information about time scales, zone ages, etc.
time_scale <- chronostratigraphic_databases$time_scale;
zone_database <- chronostratigraphic_databases$zones;
rock_superposition <- rock_unit_data$rock_superposition_database;
fossilworks_collections <- paleodb_fixes$fossilworks_collections;
paleodb_rock_reidentifications <- paleodb_fixes$paleodb_rock_reidentifications;
paleodb_collection_edits <- paleodb_fixes$paleodb_collection_edits;
if (!is.na(match("X",colnames(paleodb_collection_edits))))
	paleodb_collection_edits[,match("X",colnames(paleodb_collection_edits)):ncol(paleodb_collection_edits)] <- NULL;
#finest_chronostrat <- read.csv(paste(data_for_R_folder,"Slices_Final.csv",sep=""),header=T,stringsAsFactors = F,fileEncoding = "UTF-8")
finest_chronostrat <- pbdb_data_list$time_scale;

#### Download Data ####
initial_data <- read.csv(paste(prior_clam_directory,"FO_and_LO_Wagner.csv",sep=""),header=T,stringsAsFactors = F)
initial_data <- initial_data[initial_data$Role!="",];
initial_data$ma_lb <- gradstein_2012_emended$time_scale$ma_lb[match(initial_data$First.Occurrence,gradstein_2012_emended$time_scale$interval)];
initial_data$ma_ub <- gradstein_2012_emended$time_scale$ma_ub[match(initial_data$Last.Occurrence,gradstein_2012_emended$time_scale$interval)];

pbdb_data_list$pbdb_finds$accepted_name <- gsub("Crassatella \\(Crassatella\\)","Crassatella",pbdb_data_list$pbdb_finds$accepted_name)
pbdb_data_list$pbdb_finds$accepted_name[pbdb_data_list$pbdb_finds$identified_name=="Crassatellites carolinensis"] <- "Crassatella carolinensis";
#pbdb_data_list$pbdb_finds$flags[pbdb_data_list$pbdb_finds$occurrence_no %in% c(97978,45522,45542,45562,45667,100392,100607,100471)] <- "uncertain species";

clade_finds <- rbind(pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$accepted_name_orig %in% initial_data$Species,],
					 pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$accepted_name %in% initial_data$Species,]);
clade_finds <- clade_finds[match(unique(clade_finds$occurrence_no),clade_finds$occurrence_no),];
#missing_species <- initial_data$Species[!initial_data$Species %in% c(clade_finds$accepted_name,clade_finds$accepted_name_orig)];
#other_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$identified_name %in% missing_species,];
#pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$collection_no==97077,]

clade_sites <- pbdb_data_list$pbdb_sites_refined[pbdb_data_list$pbdb_sites_refined$collection_no %in% clade_finds$collection_no,];
#write.csv(clade_sites,"Rowans_Clams.csv",row.names = F,fileEncoding = "UTF-8");
rockless_clade_sites <- clade_sites[clade_sites$rock_no_sr==0 & clade_sites$formation!="",]
#sort(unique(rockless_clade_sites$formation[rockless_clade_sites$formation!=""]))

rockless_clade_sites_2 <- subset(clade_sites,clade_sites$formation=="");
rockless_clade_sites_2 <- subset(rockless_clade_sites_2,rockless_clade_sites_2$member=="");
rockless_clade_sites_2 <- subset(rockless_clade_sites_2,rockless_clade_sites_2$stratgroup=="");
paste(rockless_clade_sites_2$collection_no,collapse=",")

#clade_sites <- read.csv("Rowans_Clams.csv",header = T,fileEncoding = "UTF-8",stringsAsFactors = hell_no);
clade_finds$ma_lb <- clade_sites$ma_lb[match(clade_finds$collection_no,clade_sites$collection_no)];
clade_finds$ma_ub <- clade_sites$ma_ub[match(clade_finds$collection_no,clade_sites$collection_no)];
dummy_find <- clade_finds[clade_finds$collection_no==221123,];
dummy_find$identified_name <- dummy_find$accepted_name <- dummy_find$accepted_name_orig <- "Crassatella carolinensis";
clade_finds <- rbind(clade_finds,dummy_find);
clade_finds <- clade_finds[order(clade_finds$collection_no),];

crassatellids <- initial_data$Species;
nspec <- length(crassatellids);
crassatellid_appearances <- data.frame(species=as.character(crassatellids),fa_lb=as.numeric(1:nspec),fa_ub=as.numeric(1:nspec),la_lb=as.numeric(1:nspec),la_ub=as.numeric(1:nspec),
									   all_finds=as.numeric(1:nspec),fa_finds=as.numeric(1:nspec),la_finds=as.numeric(1:nspec),discarded_finds=as.numeric(1:nspec),stringsAsFactors = F);

#crassatella2_finds <- accersi_occurrence_data("Crassatella",save_file=F);
#crassatella2_sites <- accersi_collection_data("Crassatella",save_file=F);
#write.csv(crassatella2_sites,"More_Sites.csv",row.names = F,fileEncoding = "UTF-8");
#more_finds <- read.csv("Crassatella_Occurrences_Hodgei.csv",header=T,fileEncoding = "UTF-8",stringsAsFactors = F);
#clade_finds <- rbind(clade_finds,more_finds)
#write.csv(clade_finds,"Crassatella_PBDB_Finds.csv",row.names = F,fileEncoding = "UTF-8")
#clade_finds <- read.csv("Crassatella_PBDB_Finds.csv",header = T,fileEncoding = "UTF-8",stringsAsFactors = F)
#
for (sp in 1:nspec)	{
#	which(clade_finds==crassatellids[sp],arr.ind==T)
	species_finds <- clade_finds[sort(unique(which(clade_finds==crassatellids[sp],arr.ind = T)[,1])),];
	species_finds <- species_finds[species_finds$flags!="uncertain species",];
	species_finds_orig <- species_finds <- species_finds[!is.na(species_finds$ma_lb),];
	species_finds <- species_finds[species_finds$ma_lb<=initial_data$ma_lb[sp] & species_finds$ma_ub>=initial_data$ma_ub[sp],];
	if (nrow(species_finds)==0)	{
		species_finds <- clade_finds[sort(unique(which(clade_finds==crassatellids[sp],arr.ind = T)[,1])),];
		species_finds <- species_finds[!is.na(species_finds$ma_lb),];
		species_finds <- species_finds[species_finds$flags!="uncertain species",];
		}
	species_finds <- species_finds[order(-species_finds$ma_lb,-species_finds$ma_ub),];
	sort(unique(clade_sites$formation[clade_sites$collection_no %in% species_finds$collection_no][clade_sites$formation[clade_sites$collection_no %in% species_finds$collection_no]!=""]))
	crassatellid_appearances$fa_lb[sp] <- max(species_finds$ma_lb);
	crassatellid_appearances$fa_ub[sp] <- max(species_finds$ma_ub);
	crassatellid_appearances$la_lb[sp] <- min(species_finds$ma_lb);
	crassatellid_appearances$la_ub[sp] <- min(species_finds$ma_ub);
	crassatellid_appearances$all_finds[sp] <- nrow(species_finds);
	crassatellid_appearances$fa_finds[sp] <- sum(species_finds$ma_lb>crassatellid_appearances$fa_ub[sp]);
	crassatellid_appearances$la_finds[sp] <- sum(species_finds$ma_ub<crassatellid_appearances$la_lb[sp]);
	crassatellid_appearances$discarded_finds[sp] <- nrow(species_finds_orig)-nrow(species_finds);
	}

crassatellid_appearances$fa_lb <- ceiling(10*crassatellid_appearances$fa_lb)/10;
crassatellid_appearances$fa_ub <- floor(10*crassatellid_appearances$fa_ub)/10;
crassatellid_appearances$la_lb <- ceiling(10*crassatellid_appearances$la_lb)/10;
crassatellid_appearances$la_ub <- floor(10*crassatellid_appearances$la_ub)/10;
write.csv(crassatellid_appearances,"Crassatellid FAs & LAs.csv",row.names = F,fileEncoding = "UTF-8");

sort(unique(rockless_clade_sites$formation[rockless_clade_sites$formation!=""]))
rocks <- unique(cbind(rockless_clade_sites$formation,rockless_clade_sites$member));
rocks <- rocks[order(rocks[,1],rocks[,2]),];

max_ma <- time_scale$ma_lb[match(onset,time_scale$interval)];
min_ma <- time_scale$ma_ub[match(end,time_scale$interval)];

basic_environments <- c("marine","unknown");	# choices: marine; terr; unknown
pbdb_sites <- pbdb_data_list$pbdb_sites_refined[pbdb_data_list$pbdb_sites_refined$ma_lb>min_ma & pbdb_data_list$pbdb_sites_refined$ma_ub<max_ma,];
pbdb_finds <- pbdb_data_list$pbdb_finds[pbdb_data_list$pbdb_finds$collection_no %in% pbdb_sites$collection_no,];
pbdb_finds <- expello_indeterminate_species(pbdb_finds);
pbdb_finds <- pbdb_finds[!pbdb_finds$phylum %in% c("Problematica",plant_phyla,protist_phyla),];
pbdb_finds$order[pbdb_finds$order_no==0] <- "";
nfinds <- nrow(pbdb_finds);
{}
