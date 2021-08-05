#### Setup Program ####
source('~/Documents/R_Projects/Common_R_Source_Files/Chronos.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Data_Downloading_v4.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Historical_Diversity_Metrics.r');	#
source('~/Documents/R_Projects/Common_R_Source_Files/Occurrence_Data_Routines.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Sampling_and_Occupancy_Distributions.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Stratigraphy.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r'); 		#
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_Stats_and_Probability_101.r'); 		#
therapsid_directory <- "~/Documents/R_Projects/Therapsid_Mania/";	# directory to folder where you keep common source
setwd(therapsid_directory);

load("~/Documents/R_Projects/Data_for_R/Gradstein_2012_Augmented.RData"); # refined Gradstein 2012 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/Rock_Unit_Database.RData"); # refined Gradstein 2012 timescale & biozonations
load("~/Documents/R_Projects/Data_for_R/PaleoDB_Edits.RData"); # refined Gradstein 2012 timescale & biozonations
#load("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData"); # refined Gradstein 2012 timescale & biozonations

temporal_precision <- 0.05;
time_scale <- gradstein_2012_emended$time_scale;
time_scale$ma_lb <- floor(time_scale$ma_lb/temporal_precision)*temporal_precision;
time_scale$ma_ub <- ceiling(time_scale$ma_ub/temporal_precision)*temporal_precision;
zone_database <- gradstein_2012_emended$zones;
zone_database$ma_lb <- floor(zone_database$ma_lb/temporal_precision)*temporal_precision;
zone_database$ma_ub <- ceiling(zone_database$ma_ub/temporal_precision)*temporal_precision;

# Download Data ####
therapsid_finds <- accersi_occurrence_data("Therapsida",onset="Paleozoic",end="Triassic");
therapsid_sites <- accersi_collection_data("Therapsida",onset="Paleozoic",end="Triassic");
therapsid_sites$late_interval[therapsid_sites$late_interval==""] <- therapsid_sites$early_interval[therapsid_sites$late_interval==""];

hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(therapsid_sites$early_interval)),unique(as.character(therapsid_sites$late_interval)))),time_scale,regional_scale="International");
post_paleozoic <- subset(hierarchical_chronostrat,hierarchical_chronostrat$ma_lb<242);
finest_post_paleozoic <- post_paleozoic[post_paleozoic$bin_first==post_paleozoic$bin_last,];
hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(therapsid_sites$early_interval)),unique(as.character(therapsid_sites$late_interval)))),time_scale,regional_scale="Stage Slice");
finest_paleozoic <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- rbind(finest_paleozoic,finest_post_paleozoic);
finest_chronostrat <- finest_chronostrat[order(finest_chronostrat$ma_lb),]
finest_chronostrat <- finest_chronostrat[unique(match(finest_chronostrat$ma_lb,finest_chronostrat$ma_lb)),];
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];

therapsid_sites <- reparo_unedittable_paleodb_collections(therapsid_sites,paleodb_fixes$paleodb_collection_edits);

therapsid_taxa <- sort(unique(therapsid_finds$accepted_name));
th_sites <- nrow(therapsid_sites);

therapsid_otus <- read.csv(paste(therapsid_directory,"Therapsid_Ages.csv",sep=""),header=T,fileEncoding = "UTF-8",stringsAsFactors = hell_no);
therapsid_otus$species <- gsub("_"," ",therapsid_otus$species);
taxonomic_information <- accersi_taxonomic_data_for_list_of_species(species_list=therapsid_otus$species);
current_names <- taxonomic_information$entered_species$accepted_name[match(therapsid_otus$species,taxonomic_information$entered_species$taxon_name)];
therapsid_otus$species[!is.na(current_names)] <- current_names[!is.na(current_names)];
#cbind(taxonomic_information$entered_species$taxon_name,taxonomic_information$entered_species$accepted_name)
nspecies <- nrow(therapsid_otus);
clade_sites <- therapsid_sites[1,];
clade_sites <- clade_sites[clade_sites$max_ma>1000,];
for (nn in 1:nspecies)	{
	species_sites <- therapsid_sites[therapsid_sites$collection_no %in% unique(therapsid_finds$collection_no[unique(which(therapsid_finds==therapsid_otus$species[nn],arr.ind=T)[,1])]),];
	clade_sites <- rbind(clade_sites,species_sites);
	}
clade_sites <- unique(clade_sites);
clade_sites <- clade_sites[order(clade_sites$collection_no),];

clade_rock_info_2 <- organize_pbdb_rock_data(paleodb_collections = clade_sites);
pbdb_rocks <- clade_rock_info_2$pbdb_rocks;
lower_pbdb_rocks <- tolower_dataframe(pbdb_rocks);
rock_database <- rock_unit_data$rock_unit_database;
lower_wagner_rocks <- tolower_dataframe(rock_database);
geoplate_list <- clade_rock_info_2$geoplate_list;
rock_no <- 1:nrow(pbdb_rocks);
rock_no_sr <- pbapply::pbsapply(rock_no,match_organized_pbdb_rock_to_external_database,pbdb_rocks=lower_pbdb_rocks,wagner_rocks=lower_wagner_rocks,geoplate_list=geoplate_list);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no_sr=as.numeric(rock_no_sr), .before = 1);
write.csv(pbdb_rocks,"Analyzed_Therapsids_Rock.csv",row.names = F);

rock_to_zone_database <- rock_unit_data$rock_to_zone_database;
clade_sites_refinement <- refine_pbdb_collections_w_external_databases(paleodb_collections=clade_sites,rock_database,zone_database,rock_to_zone_database,finest_chronostrat);
clade_sites_refined <- clade_sites_refinement$refined_collections;
write.csv(clade_sites_refined,paste(therapsid_directory,"Therapsid_Sites_Refined.csv",sep=""),row.names = F)

# find upper and lower bounds of possible FA & LA for each species ##
for (nn in 1:nspecies)	{
	species_sites <- clade_sites_refined[clade_sites_refined$collection_no %in% unique(therapsid_finds$collection_no[unique(which(therapsid_finds==therapsid_otus$species[nn],arr.ind=T)[,1])]),];
	if (nn==1)	{
		therapsid_info <- data.frame(species=as.character(rep(therapsid_otus$species[nn],nrow(species_sites))),
									 formation=as.character(species_sites$formation),
									 member=as.character(species_sites$member),
									 group=as.character(species_sites$stratgroup),
									 ma_lb=as.numeric(species_sites$ma_lb),ma_ub=as.numeric(species_sites$ma_ub),
									 early_interval=as.character(species_sites$early_interval),
									 late_interval=as.character(species_sites$late_interval),
									 zone=as.character(species_sites$zone));
		} else	{
		dummy <- data.frame(species=as.character(rep(therapsid_otus$species[nn],nrow(species_sites))),
									 formation=as.character(species_sites$formation),
									 member=as.character(species_sites$member),
									 group=as.character(species_sites$stratgroup),
									 ma_lb=as.numeric(species_sites$ma_lb),ma_ub=as.numeric(species_sites$ma_ub),
									 early_interval=as.character(species_sites$early_interval),
									 late_interval=as.character(species_sites$late_interval),
									 zone=as.character(species_sites$zone));
			
		therapsid_info <- rbind(therapsid_info,dummy);
		}
	}

therapsid_info <- unique(therapsid_info);
therapsid_rock_info <- therapsid_info;
therapsid_rock_info$species <- NULL;
therapsid_rock_info <- unique(therapsid_rock_info);
write.csv(therapsid_info,paste(therapsid_directory,"Therapsid_PBDB_Summaries.csv",sep=""),row.names = F);
write.csv(therapsid_rock_info,paste(therapsid_directory,"Therapsid_Rocks.csv",sep=""),row.names = F);

clade_intervals <- data.frame(taxon=as.character(therapsid_otus$species),
							  min_ub=as.numeric(rep(0,nspecies)),min_lb=as.numeric(rep(0,nspecies)),
							  max_ub=as.numeric(rep(0,nspecies)),max_lb=as.numeric(rep(0,nspecies)),
							  stringsAsFactors = hell_no);

for (nn in 1:nspecies)	{
	clade_intervals$min_ub[nn] <- min(therapsid_info$ma_ub[therapsid_info$species==therapsid_otus$species[nn]]);
	clade_intervals$min_lb[nn] <- min(therapsid_info$ma_lb[therapsid_info$species==therapsid_otus$species[nn]]); 
	clade_intervals$max_ub[nn] <- max(therapsid_info$ma_ub[therapsid_info$species==therapsid_otus$species[nn]]);
	clade_intervals$max_lb[nn] <- max(therapsid_info$ma_lb[therapsid_info$species==therapsid_otus$species[nn]]);
	}
missing_taxa <- clade_intervals$taxon[is.infinite(clade_intervals$min_ub)];	# this should be empty!!
write.csv(clade_intervals,"Therapsid_Fossil_Info.csv",row.names = F);

# Get Info for Synapsids & Archosaurs ####
all_finds <- accersi_occurrence_data("Synapsida,Archosauromorpha",onset="Permian",end="Triassic",directory = therapsid_directory);
all_sites <- accersi_collection_data("Synapsida,Archosauromorpha",onset="Permian",end="Triassic",directory = therapsid_directory);
all_sites <- reparo_unedittable_paleodb_collections(all_sites,paleodb_fixes$paleodb_collection_edits);

ttl_sites <- nrow(all_sites);
all_taxa <- sort(unique(all_finds$accepted_name));

#all_sites <- unique(rbind(therapsid_sites,archosaur_sites));
all_sites <- all_sites[order(all_sites$collection_no),];
all_sites <- all_sites[all_sites$max_ma>min(therapsid_otus$min),];

all_sites_refinement <- refine_pbdb_collections_w_external_databases(paleodb_collections=all_sites,rock_database,zone_database,rock_to_zone_database,finest_chronostrat);
all_sites_refined <- all_sites_refinement$refined_collections;

# Separate "Youngest" Finds ####
youngest <- (1:nspecies)[clade_intervals$max_lb %in% min(clade_intervals$max_lb)];
latest_bounds <- c(clade_intervals$max_lb[youngest[1]],max(clade_intervals$max_ub[youngest]))
rho_sites <- all_sites_refined[all_sites_refined$ma_lb>min(latest_bounds),];
rho_sites <- rho_sites[rho_sites$ma_ub<max(latest_bounds),];
rho_finds <- all_finds[all_finds$collection_no %in% rho_sites$collection_no,];
rho_finds <- rho_finds[rho_finds$accepted_rank %in% c("subspecies","species"),];
rho_taxa <- sort(unique(rho_finds$accepted_name));
rho_taxonomic_information <- accersi_taxonomic_data_for_list_of_species(species_list=rho_taxa);

# get rid of taxa based on footprints or other trace fossils
httpS <- "https://paleobiodb.org/data1.2/occs/taxa.csv?base_name=Synapsida&pres=form,ichno&interval=Triassic,Triassic&show=full";
httpA <- "https://paleobiodb.org/data1.2/occs/taxa.csv?base_name=Archosauromorpha&pres=form,ichno&interval=Triassic,Triassic&show=full";
ichno_taxa <-rbind(read.csv(httpS,header=T,stringsAsFactors=hell_no,encoding="UTF-8"),read.csv(httpA,header=T,stringsAsFactors=hell_no,encoding="UTF-8"));
rho_taxa <- rho_taxa[!rho_taxa %in% ichno_taxa$accepted_name];

rho_sites <- rho_sites[rho_sites$collection_no %in% rho_finds$collection_no,];
write.csv(rho_sites,"Last_Call_Sites.csv",row.names = F,fileEncoding = "UTF-8");
write.csv(rho_finds,"Last_Call_Finds.csv",row.names = F,fileEncoding = "UTF-8");

rtaxa <- length(rho_taxa);
taxon_sites <- taxon_rocks <- vector(length=rtaxa);
rho_sites <- name_unnamed_rock_units(paleodb_collections = rho_sites,finest_chronostrat = finest_chronostrat);
for (rr in 1:rtaxa)	{
	rr_sites <- rho_sites[rho_sites$collection_no %in% rho_finds$collection_no[rho_finds$accepted_name==rho_taxa[rr]],];
	rr_sites <- unique(rho_finds$collection_no[rho_finds$accepted_name==rho_taxa[rr]]);
	rr_rocks <- unique(rho_sites$rock_no_sr[match(rr_sites,rho_sites$collection_no)]);
	if (length(rr_rocks)>1)	rr_rocks <- rr_rocks[rr_rocks>0];
	taxon_sites[rr] <- length(rr_sites);
	taxon_rocks[rr] <- length(unique(rr_rocks));
	}
control_taxa_info <- data.frame(species=as.character(rho_taxa),
								sites=as.numeric(taxon_sites),
								rock_units=as.numeric(taxon_rocks));
write.csv(control_taxa_info,"Rho_Taxa_Sampling.csv",row.names = F);

# Now, get sampling! ####
par(pin=c(3,3));
hist(taxon_sites,breaks=0:max(taxon_sites),xlab="Fossiliferous Sites Occupied",ylab="No. Large Vertebrate Species",main="");
hist(taxon_rocks,breaks=0:max(taxon_rocks),xlab="Fossiliferous Rock-Units Occupied",ylab="No. Large Vertebrate Species",main="");

# contrast best uniform, exponential, beta & lognormal distributions
# 	First by sites (= collections or localities)
rsites <- nrow(rho_sites);
taxon_sites <- sort(taxon_sites,decreasing=T);
best_uniform_sites <- optimize_uniform_occupancy(finds=taxon_sites,ncoll=rsites);
best_exponential_sites <- optimize_exponential_occupancy(finds=taxon_sites,ncoll=rsites);
best_beta_sites <- optimize_beta_occupancy(finds=taxon_sites,ncoll=rsites);
best_lognormal_sites <- optimize_lognormal_occupancy(finds=taxon_sites,ncoll=rsites);
best_possible_sites <- saturated_occupancy_model(finds=taxon_sites);

# 	Second by rock-units (formations, members &/or groups)
rrocks <- length(unique(rho_sites$rock_no_sr));
taxon_rocks <- sort(taxon_rocks,decreasing=T);
best_uniform_rocks <- optimize_uniform_occupancy(finds=taxon_rocks,ncoll=rrocks);
best_exponential_rocks <- optimize_exponential_occupancy(finds=taxon_rocks,ncoll=rrocks);
best_beta_rocks <- optimize_beta_occupancy(finds=taxon_rocks,ncoll=rrocks);
best_lognormal_rocks <- optimize_lognormal_occupancy(finds=taxon_rocks,ncoll=rrocks);
best_possible_rocks <- saturated_occupancy_model(finds=taxon_rocks);

best_uniforms <- rbind(best_uniform_sites,best_uniform_rocks);
rownames(best_uniforms) <- c("Sites","Rocks");
best_exponentials <- rbind(best_exponential_sites,best_exponential_rocks);
best_betas <- rbind(best_beta_sites,best_beta_rocks);
best_lognormals <- rbind(best_lognormal_sites,best_lognormal_rocks);
best_possible <- rbind(best_possible_sites,best_possible_rocks);

best_distributions <- cbind(best_uniforms,best_exponentials,best_betas,best_lognormals,best_possible);
write.csv(best_distributions,"Best_Sampling_Occupancy_Distributions.csv",row.names = T);

aics <- (1:ncol(best_distributions))[colnames(best_distributions) %in% "AICc"];
richness <- (1:ncol(best_distributions))[colnames(best_distributions) %in% "richness"];
best_sites_s <- best_distributions[1,richness[match(min(best_distributions[1,aics]),best_distributions[1,aics])]];
best_rocks_s <- best_distributions[2,richness[match(min(best_distributions[2,aics]),best_distributions[2,aics])]];
rho <- ((rtaxa/best_sites_s)+(rtaxa/best_rocks_s))/2;

print(rho)

{}