reupdate_opinions <- FALSE;
# setup functions ####
update_opinions <- function(taxon=c("Life","Ichnofossils"),ops_modified_after="1990-01-01",inc_children=TRUE)	{
#	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod")
options(warn=-1);
taxon <- paste(taxon, collapse = ",");
if (inc_children)	{
	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTO <- paste("https://paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&ops_modified_after=",ops_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}

pbdb_data <- read.csv(httpTO,header=T,stringsAsFactors = F);
new_opinions <- put_pbdb_dataframes_into_proper_type(pbdb_data);
options(warn=1);
return(new_opinions)
}

update_taxonomy <- function(taxon=c("Life","Ichnofossils"),taxa_modified_after="1990-01-01",inc_children=TRUE)	{
#httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,",&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,",&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="")
taxon <- paste(taxon, collapse = ",");
options(warn=-1);
if (inc_children)	{
	if (taxa_modified_after!="1990-01-01")	{
		httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&variant=all&rel=all_children&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
		} else {
		httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&variant=all&rel=all_children&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
		}
#	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#	httpTt <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&ops_modified_after=",taxa_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	if (taxa_modified_after!="1990-01-01")	{
		httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
		} else	{
		httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&private&variant=all&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,refattr,ent,entname,crmod",sep="");
		}
#	httpTt <- "https://paleobiodb.org/data1.2/taxa/list.csv?base_name=Deuterostomia&variant=all&private&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,seq,img,ref,ent,entname,crmod";
#	httpTt <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&taxa_modified_after=",taxa_modified_after,"&private&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr&all_records",sep="");
#	httpTt <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&ops_modified_after=",taxa_modified_after,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}
#fetch <- RCurl::getURL(httpTt,timeout(60*60));
#pbdb_data <- read.csv(httpTt,header=T,stringsAsFactors = F,timeout(60*60));
pbdb_data <- read.csv(httpTt,header=TRUE,stringsAsFactors = FALSE);
new_taxa <- put_pbdb_dataframes_into_proper_type(pbdb_data);
options(warn=1);
return(new_taxa)
}

clean_pbdb_fields <- function(pbdb_data)	{
dirty_fields <- c("collection_name","collection_aka","state","county","geogcomments","formation","member","stratgroup","county","localsection","regionalsection","stratcomments","lithdescript","geology_comments","author","ref_author","primary_reference","collection_comments","taxonomy_comments","primary_reference","authorizer","enterer","modifier");
d_f <- (1:ncol(pbdb_data))[colnames(pbdb_data) %in% dirty_fields];
for (df in 1:length(d_f))	{
	web_text <- pbdb_data[,d_f[df]];
	if (sum(web_text!="")>0)	pbdb_data[,d_f[df]][web_text!=""] <- pbapply::pbsapply(web_text[web_text!=""],mundify_web_text_dull);
	}
return(pbdb_data);
}

#pbdb_data <- updated_taxonomy;
clean_the_bastards <- function(pbdb_data)	{
pbdb_data$flags[pbdb_data$flags=="FALSE"] <- "F";
pbdb_data$flags <- gsub("B","",pbdb_data$flags);
pbdb_data$phylum_no[pbdb_data$phylum=="NO_PHYLUM_SPECIFIED"] <- 0;
pbdb_data$phylum[pbdb_data$phylum=="NO_PHYLUM_SPECIFIED"] <- "";
pbdb_data$phylum_no[is.na(pbdb_data$phylum_no)] <- 0;
pbdb_data$phylum[is.na(pbdb_data$phylum_no)] <- "";
pbdb_data$class_no[pbdb_data$class=="NO_CLASS_SPECIFIED"] <- 0;
pbdb_data$class[pbdb_data$class=="NO_CLASS_SPECIFIED"] <- "";
pbdb_data$class_no[is.na(pbdb_data$class_no)] <- 0;
pbdb_data$class[is.na(pbdb_data$class_no)] <- "";
pbdb_data$order_no[pbdb_data$order=="NO_ORDER_SPECIFIED"] <- 0;
pbdb_data$order[pbdb_data$order=="NO_ORDER_SPECIFIED"] <- "";
pbdb_data$order_no[is.na(pbdb_data$order_no)] <- 0;
pbdb_data$order[is.na(pbdb_data$order_no)] <- "";
pbdb_data$family_no[pbdb_data$family=="NO_FAMILY_SPECIFIED"] <- 0;
pbdb_data$family[pbdb_data$family=="NO_FAMILY_SPECIFIED"] <- "";
pbdb_data$family_no[is.na(pbdb_data$family_no)] <- 0;
pbdb_data$family[is.na(pbdb_data$family_no)] <- "";
pbdb_data$genus_no[pbdb_data$genus=="NO_GENUS_SPECIFIED"] <- 0;
pbdb_data$genus[pbdb_data$genus=="NO_GENUS_SPECIFIED"] <- "";
pbdb_data$phylum_no <- as.numeric(pbdb_data$phylum_no);
pbdb_data$class_no <- as.numeric(pbdb_data$class_no);
pbdb_data$order_no <- as.numeric(pbdb_data$order_no);
pbdb_data$family_no <- as.numeric(pbdb_data$family_no);
pbdb_data$genus_no <- as.numeric(pbdb_data$genus_no);
pbdb_data$genus_no[is.na(pbdb_data$genus_no)] <- 0;
if (!is.null(pbdb_data$subgenus_no))	{
	pbdb_data$subgenus_no[pbdb_data$subgenus_no==""] <- 0;
	pbdb_data$subgenus_no[is.na(pbdb_data$subgenus_no)] <- 0;
	pbdb_data$subgenus_no <- as.numeric(pbdb_data$subgenus_no);
	}
return(pbdb_data)
}

# main function ####
update_pbdb_RData <- function(pbdb_data_list,rock_unit_data,gradstein_2020_emended,paleodb_fixes,do_not_forget_taxa=c(),do_not_forget_start="Archean",do_not_forget_end="Phanerozoic",reupdate_opinions=FALSE)	{
start_time <- date();
# reload current RData & setup relevant databases ####
pbdb_finds <- put_pbdb_dataframes_into_proper_type(paleodb_data=pbdb_data_list$pbdb_finds);
pbdb_sites <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_sites);
#pbdb_sites_refined <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_sites_refined);
pbdb_taxonomy <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_taxonomy);
pbdb_opinions <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_opinions);
pbdb_rocks <- pbdb_data_list$pbdb_rocks;
pbdb_rocks_sites <- pbdb_data_list$pbdb_rocks_sites;

time_scale <- gradstein_2020_emended$time_scale;
finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
zone_database <- gradstein_2020_emended$zones;

rock_database <- rock_unit_data$rock_unit_database;

fossilworks_collections <- paleodb_fixes$fossilworks_collections;
pbdb_taxonomy_fixes <- paleodb_fixes$paleodb_taxonomy_edits;

pbdb_taxonomy_fixes <- clear_na_from_matrix(pbdb_taxonomy_fixes,"");
pbdb_taxonomy_fixes <- clear_na_from_matrix(pbdb_taxonomy_fixes,0);
pbdb_taxonomy_species_fixes <- pbdb_taxonomy_fixes[pbdb_taxonomy_fixes$taxon_rank %in% c("species","subspecies"),];
pbdb_taxonomy_fixes$created <- gsub(" UTC","",pbdb_taxonomy_fixes$created);
pbdb_taxonomy_fixes$modified <- gsub(" UTC","",pbdb_taxonomy_fixes$modified);

# I hate that I have to do this!!! ####
occurrence_no <- pbdb_finds$occurrence_no[is.na(pbdb_finds$modified)];
fixed_finds <- data.frame(base::t(pbapply::pbsapply(occurrence_no,update_occurrence_from_occurrence_no)));
new_modified <- unlist(fixed_finds$modified);
if (sum(new_modified=="0000-00-00 00:00:00")>0)	{
	new_created <- unlist(fixed_finds$created);
	new_modified[new_modified=="0000-00-00 00:00:00"] <- new_created[new_modified=="0000-00-00 00:00:00"];
	}

pbdb_finds$modified[is.na(pbdb_finds$modified)] <- as.character.Date(new_modified);
#pbdb_finds$modified[pbdb_finds$occurrence_no %in% occurrence_no];

# get new/revised data ####
print("Getting occurrence & locality data entered/modified since the last download");
occs_modified_after <- as.Date(max(pbdb_finds$modified[!is.na(pbdb_finds$modified)]))-2;
print(paste("Looking for data modified after",occs_modified_after));
new_finds <- update_occurrence_data(occs_modified_after = occs_modified_after,basic_environments = c("terrestrial","marine","unknown"),species_only=FALSE,save_files = FALSE);
print(paste("Adding data modified as recently as",max(new_finds$modified)));
if (length(do_not_forget_taxa)>0 && sum(do_not_forget_taxa=="Metazoa")!=length(do_not_forget_taxa))	{
	print("Adding additional finds from specified taxonomic group(s)");
	other_finds <- accersi_occurrence_data(taxa = do_not_forget_taxa,onset=do_not_forget_start,end=do_not_forget_end,species_only=F,save_files = F);
#	other_finds <- accersi_occurrence_data(taxa = do_not_forget,onset="Cretaceous",end="Cenozoic",species_only=F,save_files = F);
	new_finds <- rbind(new_finds,other_finds);
	} else if (!do_not_forget_start %in% c("","Archean") && !do_not_forget_end %in% c("","Phanerozoic"))	{
	if (do_not_forget_start==do_not_forget_end)	{
		print(paste("Adding additional finds from the ",do_not_forget_start,".",sep=""))
		} else	{
		print(paste("Adding additional finds from the ",do_not_forget_start,"-",do_not_forget_end,".",sep=""))
		}
#	print("Adding additional finds from specified interval(s)");
	other_finds <- accersi_occurrence_data(taxa = "Life",onset=do_not_forget_start,end=do_not_forget_end,save_files = F,species_only = F);
	new_finds <- rbind(new_finds,other_finds);
	}
new_finds <- put_pbdb_dataframes_into_proper_type(new_finds);
new_finds <- clean_pbdb_fields(new_finds);
new_finds <- new_finds[match(unique(new_finds$occurrence_no),new_finds$occurrence_no),];
#pbdb_data_list$pbdb_finds <- tibble::add_column(pbdb_data_list$pbdb_finds,updater=as.character(rep("",nrow(pbdb_data_list$pbdb_finds))), .before = match("abund_value",colnames(pbdb_data_list$pbdb_finds)));
#pbdb_data_list$pbdb_finds <- tibble::add_column(pbdb_data_list$pbdb_finds,updated=pbdb_data_list$pbdb_finds$modified, .before = match("accepted_name_orig",colnames(pbdb_data_list$pbdb_finds)));
#cbind(c(colnames(pbdb_data_list$pbdb_finds),"",""),colnames(new_finds))
colls_modified_after <- as.Date(min(occs_modified_after,max(pbdb_sites$modified[!is.na(pbdb_sites$modified)])));
options(warn=-1);
new_sites <- update_collection_data(colls_modified_after = occs_modified_after,basic_environments = c("terrestrial","marine","unknown"),save_files = F);
if (length(do_not_forget_taxa)>0 && sum(do_not_forget_taxa=="Metazoa")!=length(do_not_forget_taxa))	{
	print("Adding additional sites from specified taxonomic group(s)");
	other_sites <- accersi_collection_data(taxa = do_not_forget_taxa,onset=do_not_forget_start,end=do_not_forget_end,species_only=F,save_files = F);
#	other_sites <- accersi_collection_data(taxa = do_not_forget,onset="Cretaceous",end="Cenozoic",save_files = F);
	new_sites <- rbind(new_sites,other_sites);
	} else if (!do_not_forget_start %in% c("","Archean") && !do_not_forget_end %in% c("","Phanerozoic"))	{
	print("Adding additional sites from specified interval(s)");
	other_sites <- accersi_collection_data(taxa = "Life",onset=do_not_forget_start,end=do_not_forget_end,save_files = F);
	new_sites <- rbind(new_sites,other_sites);
	}

# find "lost" collections and get their information
coll_id <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites$collection_no]);
coll_id <- sort(c(coll_id,pbdb_sites$collection_no[!pbdb_sites$collection_no %in% pbdb_data_list$pbdb_sites_refined$collection_no]));
coll_id <- coll_id[!coll_id %in% new_sites$collection_no];
if (length(coll_id)>0)	{
	print(paste("Adding",length(coll_id),"lost sites"));
	ns <- 1;
	lost_sites <- accersi_data_for_one_collection(coll_id[ns]);
	while (ns < length(coll_id))	{
		ns <- ns+1;
		lost_sites <- rbind(lost_sites,accersi_data_for_one_collection(coll_id[ns]));
		}
#	colnames(new_sites)
#	colnames(lost_sites)[!colnames(lost_sites) %in% colnames(new_sites)]
	new_sites <- rbind(new_sites,lost_sites);
	new_sites <- new_sites[order(new_sites$modified,decreasing = T),];
	}
#other_sites <- accersi_collection_data(taxa = "Ichthyosauromorpha,Plesiosauria,Mosasauria",onset="Permian",end="Cenozoic",save_files = F);
new_sites <- new_sites[match(unique(new_sites$collection_no),new_sites$collection_no),];
new_sites <- new_sites[order(new_sites$collection_no),];
new_sites <- put_pbdb_dataframes_into_proper_type(new_sites);
new_sites <- clean_pbdb_fields(new_sites);
new_sites <- new_sites[,colnames(new_sites) %in% colnames(pbdb_sites)];
#match(colnames(pbdb_sites),colnames(new_sites))
#print(nrow(new_sites))

# added 2023-03-15
#sum(coll_id %in% new_sites$collection_no)
# added 2023-03-15
#sum(coll_id %in% new_sites$collection_no)
#print("Redoing paleolng & paleolat for new collections.");
#if (sum(new_sites$paleomodel %in% c("gp_mid","scotese"))>0)	{
#	old_sites <- new_sites[new_sites$paleomodel %in% c("gp_mid","scotese"),];
#	old_sites <- old_sites[old_sites$max_ma>1,];
#	old_sites <- mundify_paleolatitude_and_paleolongitude(old_sites=old_sites,old_model=c("gp_mid","scotese"));
#	old_sites <- old_sites[order(old_sites$collection_no),];
#	new_sites$paleolat[new_sites$collection_no %in% old_sites$collection_no] <- old_sites$paleolat;
#	new_sites$paleolng[new_sites$collection_no %in% old_sites$collection_no] <- old_sites$paleolng;
#	new_sites$paleomodel[new_sites$collection_no %in% old_sites$collection_no] <- old_sites$paleomodel;
#	}
new_coll_numbers <- sort(unique(new_sites$collection_no));
#print(nrow(new_sites))
#cd <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites$collection_no]);
options(warn=1);

# get taxonomic updates ####
#### Get new & updated opinions ####
print("Getting taxonomic opinions & data entered/modified since the last download");
ops_modified_after <- as.Date(max(pbdb_opinions$modified[!is.na(pbdb_opinions$modified)]))-1;
new_opinions <- update_opinions(ops_modified_after = ops_modified_after);
#new_opinions <- rbind(new_opinions,update_opinions("Ichnofossils",ops_modified_after = occs_modified_after));
new_opinions <- put_pbdb_dataframes_into_proper_type(new_opinions);
new_opinions <- clean_pbdb_fields(new_opinions);

if (reupdate_opinions)	{
# get older opinions for taxa with new opinions
	opinionated <- gsub(" ","%20",unique(new_opinions$taxon_name));
#parents <- unique(new_opinions$parent_name); # ? use this
#nparents <- length(parents);
	print("Collecting older opinions for taxa with new opinions");
#for (p in 1:nparents)
#	old_opinions <- rbind(old_opinions,update_opinions(taxon=parents[p],inc_children = T));

#old_opinions <- update_opinions();
	old_opinions <- pbdb_opinions;
	old_opinions <- old_opinions[old_opinions$opinion_no==0,];
	i <- 0;
	while (i < length(opinionated))	{
		i <- i+1;
		old_opinions <- rbind(old_opinions,update_opinions(taxon=opinionated[i],inc_children = FALSE));
		}
#old_opinions <- base::t(pbapply::pbsapply(opinionated,update_opinions,inc_children = F));
# eliminate duplicate opinions;
	opinion_nos <- pbdb_opinions$opinion_no;
	opinion_no_counts <- hist(opinion_nos,breaks=0:max(opinion_nos),plot=F)$counts;
	names(opinion_no_counts) <- 1:max(opinion_nos);
	dup_opinions <- opinion_no_counts[opinion_no_counts>1];
	dp <- 0;
	while (dp < length(dup_opinions))	{
		dp <- dp+1;
		fix_me <- as.numeric(names(dup_opinions)[dp]);
		prob_rows <- (1:nrow(pbdb_opinions))[pbdb_opinions$opinion_no %in% fix_me];
		last_modified <- max(pbdb_opinions$modified[pbdb_opinions$opinion_no %in% fix_me]);
		keep_me <- prob_rows[pbdb_opinions$modified[pbdb_opinions$opinion_no %in% fix_me] %in% last_modified]
		kill_me <- prob_rows[!prob_rows %in% keep_me[1]];
		pbdb_opinions <- pbdb_opinions[(1:nrow(pbdb_opinions)) %in% kill_me,];
		}
	#writexl::write_xlsx(new_finds,"new_finds.xlsx");
	#writexl::write_xlsx(new_sites,"new_sites.xlsx");
	#writexl::write_xlsx(new_opinions,"new_opinions.xlsx");
	#old_opinions <- pbdb_opinions;
	fix_me <- match(old_opinions$opinion_no,pbdb_opinions$opinion_no);
	old_opinions <- old_opinions[!is.na(fix_me),];
	fix_me <- fix_me[!is.na(fix_me)];
	pbdb_opinions[fix_me,] <- old_opinions;
	}

#### Get new & updated taxa ####
print("Collecting updated taxonomic information.");
taxa_modified_after <- ops_modified_after;
#taxa_modified_after <- as.Date(max(pbdb_taxonomy$modified[!is.na(pbdb_taxonomy$modified)]))-1;
#taxa_modified_after <- "2010-04-30";
new_taxa <- update_taxonomy(taxa_modified_after = taxa_modified_after,inc_children = TRUE);
#new_taxa <- rbind(new_taxa,update_taxonomy(taxon = "Ichnofossils",taxa_modified_after = taxa_modified_after,inc_children = T));
#new_taxa <- data.frame(readxl::read_xlsx("new_taxa.xlsx"));
new_taxa <- put_pbdb_dataframes_into_proper_type(new_taxa);
new_taxa <- clean_pbdb_fields(new_taxa);
#writexl::write_xlsx(new_taxa,"new_taxa.xlsx");
#new_taxa[new_taxa$taxon_name=="Agnostida",]
new_opinions_spc <- new_opinions[new_opinions$taxon_rank %in% c("species","subspecies"),];
new_opinions_spc <- new_opinions_spc[new_opinions_spc$opinion_type=="class",];

#### Get finds for taxa with updated opinions including new taxa ####
taxon <- new_opinions_spc$taxon_name;
taxon_id <- unique(new_opinions_spc$orig_no[!is.na(new_opinions_spc$orig_no)]);
#pbdb_list <- data.frame(base::t(pbapply::pbsapply(taxon_id,accersi_occurrences_for_one_taxon_no)));
#updated_finds <- list_to_dataframe_for_pbdb_data(pbdb_list);
#updated_finds <- data.frame(readxl::read_xlsx("updated_finds.xlsx"));
#updated_finds <- put_pbdb_dataframes_into_proper_type(updated_finds);
#updated_finds <- clean_pbdb_fields(updated_finds);
#td <- 1+max((1:length(taxon_id))[taxon_id %in% updated_finds$accepted_no]);
if (length(taxon_id)>0)	updated_finds <- accersi_occurrences_for_one_taxon_no(taxon_id[1]);
td <- 1;
options(warn=2);
#nn <- base::t(pbapply::pbsapply(taxon_id,accersi_occurrences_for_one_taxon_no));
updates <- round(length(taxon_id)*(1:19)/20,0);
post_update <- FALSE; # for debugging...
while (td < length(taxon_id))	{
	if (post_update)	writexl::write_xlsx(new_taxa[match(taxon_id[1:td],new_taxa$taxon_no),],"taxon_ids.xlsx");
	if (td %in% updates & post_update)	{
		print(paste(round(100*td/length(taxon_id),0),"% done.",sep=""));
		writexl::write_xlsx(updated_finds,"updated_finds.xlsx");
		}
#while (td <= 832)	{
#	taxon[td]
	td <- td+1;
	if(!taxon_id[td] %in% new_taxa$taxon_no || !new_taxa$difference[new_taxa$taxon_no==taxon_id[td]] %in% nomens)	{
		if (taxon_id[td]%%1e+05==0)	{
#		pbdb_taxonomy$taxon_name[pbdb_taxonomy$taxon_no==50000]
#		dummy <- paste((taxon_id[td]/1e+05),"0000",sep="");
			updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(dummy));
			} else	{
			updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(taxon_id[td]));
			}
#		updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(taxon_id[td]));
		}
	}
options(warn=1);
updated_finds <- unique(updated_finds);
updated_finds$modified[is.na(updated_finds$modified)] <- updated_finds$created[is.na(updated_finds$modified)];
updated_finds <- put_pbdb_dataframes_into_proper_type(updated_finds);
updated_finds <- clean_pbdb_fields(updated_finds);
new_finds <- rbind(new_finds,updated_finds[!updated_finds$occurrence_no %in% new_finds$occurrence_no,]);
new_finds <- new_finds[match(unique(new_finds$occurrence_no),new_finds$occurrence_no),];
#writexl::write_xlsx(new_finds,"new_finds.xlsx");

new_opinions_gen <- new_opinions[new_opinions$taxon_rank %in% c("genus","subgenus"),];
new_opinions_gen <- new_opinions_gen[new_opinions_gen$opinion_type=="class",];
taxon <- new_opinions_gen$taxon_name;
taxon_id <- unique(new_opinions_gen$orig_no[!is.na(new_opinions_gen$orig_no)]);

#### update occurrence data for name changes ####
print("Updating occurrence data for taxa with new taxonomic information.");
td <- 1;
if (length(taxon_id)>0)	updated_finds <- accersi_occurrences_for_one_taxon_no(taxon_id[1]);
options(warn=2);
updates <- round(length(taxon_id)*(1:19)/20,0);
while (td < length(taxon_id))	{
	if (td %in% updates & post_update)	{
		print(paste(round(100*td/length(taxon_id),0),"% done.",sep=""));
		writexl::write_xlsx(updated_finds,"updated_finds2.xlsx");
		}
	td <- td+1;
	if (taxon_id[td]%%1e+05==0)	{
#		pbdb_taxonomy$taxon_name[pbdb_taxonomy$taxon_no==50000]
		dummy <- paste((taxon_id[td]/1e+05),"0000",sep="");
		updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(dummy));
		} else	{
		updated_finds <- rbind(updated_finds,accersi_occurrences_for_one_taxon_no(taxon_id[td]));
		}
	}
options(warn=1);

#### Clean up finds after taxonomic updates ####
updated_finds$modified[is.na(updated_finds$modified)] <- updated_finds$created[is.na(updated_finds$modified)];
updated_finds <- put_pbdb_dataframes_into_proper_type(updated_finds);
updated_finds <- clean_pbdb_fields(updated_finds);
new_finds <- rbind(new_finds,updated_finds[!updated_finds$occurrence_no %in% new_finds$occurrence_no,]);
new_finds <- new_finds[match(unique(new_finds$occurrence_no),new_finds$occurrence_no),];

taxon_name <- new_finds$accepted_name;
informals <- pbapply::pbsapply(taxon_name,revelare_informal_taxa,keep_author_specific=TRUE);
new_finds_spc <- new_finds[!informals,];
dummy_finds <- unique(rbind(pbdb_data_list$pbdb_finds,new_finds_spc));

coll_id2 <- dummy_finds$collection_no[!dummy_finds$collection_no %in% c(pbdb_sites$collection_no,new_sites$collection_no)];
new_sites_spc <- unique(new_sites[new_sites$collection_no %in% unique(dummy_finds$collection_no),]);
cd <- unique(pbdb_finds$collection_no[!new_sites_spc$collection_no %in% pbdb_sites$collection_no]);
dummy_finds <- NULL;

# repair & redate new sites ####
print("Repairing & Redating the new & updated localities.");
new_sites_spc <- unique(new_sites_spc); # new_sites_spc <- unique(new_sites);
new_sites_spc <- new_sites_spc[order(new_sites_spc$collection_no),];
new_sites_spc <- reparo_unedittable_paleodb_collections(new_sites_spc,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);

if (sum(new_sites_spc$direct_ma>0)>0)	{
	new_sites_redated <- redate_paleodb_collections_with_direct_dates(paleodb_collections=new_sites_spc,finest_chronostrat);
	new_sites_redated <- redate_paleodb_collections_with_time_scale(paleodb_collections=new_sites_redated,time_scale=time_scale,zone_database = zone_database);
	} else	{
	new_sites_redated <- redate_paleodb_collections_with_time_scale(paleodb_collections=new_sites_spc,time_scale=time_scale,zone_database = zone_database);
	}
if (sum(new_sites_redated$zone!="")>0)	new_sites_redated <- redate_paleodb_collections_with_zones(paleodb_collections = new_sites_redated,zone_database = zone_database,time_scale = time_scale,emend_paleodb=T);

coll_no <- unique(sort(c(pbdb_sites$collection_no[is.na(pbdb_sites$paleolat)],pbdb_sites$collection_no[pbdb_sites$geoplate==0])));
#if (length(coll_no)>0)	{
#	revised_paleogeography <- data.frame(base::t(pbapply::pbsapply(coll_no,paleogeographic_orphanarium)));
#	for (cc in 1:ncol(revised_paleogeography))	revised_paleogeography[,cc] <- as.vector(unlist(revised_paleogeography[,cc]));
#	revised_paleogeography <- put_pbdb_dataframes_into_proper_type(revised_paleogeography);
#	leelas <- match(revised_paleogeography$collection_no,pbdb_sites$collection_no);
#	mutants <- match(colnames(revised_paleogeography),colnames(pbdb_sites));
#	pbdb_sites[leelas,mutants] <- revised_paleogeography;
#	}

if (is.list(pbdb_sites$zone))	pbdb_sites$zone <- as.character(unlist(pbdb_sites$zone));

# separate new finds & sites from recently edited finds & sites ####
print("Integrating modified finds & new finds.");
edited_finds <- new_finds_spc[(1:nrow(new_finds_spc))[new_finds_spc$occurrence_no %in% pbdb_finds$occurrence_no],];
added_finds <- new_finds_spc[(1:nrow(new_finds_spc))[!new_finds_spc$occurrence_no %in% pbdb_finds$occurrence_no],];
pbdb_finds[match(edited_finds$occurrence_no,pbdb_finds$occurrence_no),] <- edited_finds;
pbdb_finds <- rbind(pbdb_finds,added_finds);
pbdb_finds <- pbdb_finds[order(pbdb_finds$collection_no,pbdb_finds$occurrence_no),];

# update taxonomy that PBDB will not
pbdb_finds$accepted_name[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_name[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_name_orig[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_name[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_rank[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_rank[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$accepted_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$accepted_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$subgenus_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$subgenus_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$genus_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$genus_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$genus[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$genus[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$family_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$family_no[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];
pbdb_finds$family[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no] <- pbdb_taxonomy_species_fixes$family[match(pbdb_finds$identified_no[pbdb_finds$identified_no %in% pbdb_taxonomy_species_fixes$taxon_no],pbdb_taxonomy_species_fixes$taxon_no)];

age <- new_sites_redated$max_ma;
new_sites_redated$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
age <- new_sites_redated$min_ma;
new_sites_redated$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
new_sites_redated$created <- as.Date(new_sites_redated$created);
new_sites_redated$modified <- as.Date(new_sites_redated$modified);

if (!is.null(pbdb_sites$rock_no))		pbdb_sites$rock_no <- NULL;
if (!is.null(pbdb_sites$formation_no))	pbdb_sites$formation_no <- NULL;
edited_sites <- new_sites_redated[new_sites_redated$collection_no %in% pbdb_sites$collection_no,];
added_sites <- new_sites_redated[(1:nrow(new_sites_redated))[!new_sites_redated$collection_no %in% pbdb_sites$collection_no],];
edited_sites <- put_pbdb_dataframes_into_proper_type(edited_sites);
added_sites <- put_pbdb_dataframes_into_proper_type(added_sites);
nn <- match(edited_sites$collection_no,pbdb_sites$collection_no);

edited_sites$created <- as.Date(edited_sites$created);
edited_sites$modified <- as.Date(edited_sites$modified);
added_sites$created <- as.Date(added_sites$created);
added_sites$modified <- as.Date(added_sites$modified);
if (is.null(pbdb_sites$interval_lb))	{
	age <- pbdb_sites$max_ma;
	pbdb_sites$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
	age <- pbdb_sites$min_ma;
	pbdb_sites$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
	}
# 2023-03-20: something changed in the output here!!
edited_sites <- edited_sites[,colnames(edited_sites) %in% colnames(pbdb_sites)];
added_sites <- added_sites[,colnames(added_sites) %in% colnames(pbdb_sites)];
pbdb_sites[nn,] <- edited_sites;
pbdb_sites <- rbind(pbdb_sites,added_sites);
pbdb_sites <- pbdb_sites[order(pbdb_sites$collection_no),];
nsites <- nrow(pbdb_sites);
pbdb_sites <- reparo_unedittable_paleodb_collections(pbdb_sites,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);
time_scale <- gradstein_2020_emended$time_scale;
time_scale <- rbind(time_scale,finest_chronostrat[!finest_chronostrat$interval %in% time_scale$interval,]);
pbdb_sites <- redate_paleodb_collections_with_time_scale(paleodb_collections=pbdb_sites,time_scale,zone_database);
age <- pbdb_sites$max_ma;
pbdb_sites$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
age <- pbdb_sites$min_ma;
pbdb_sites$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
if (sum(is.na(pbdb_sites$modified))>0)	{
	fix_these1 <- (1:nsites)[is.na(pbdb_sites$created)];
	fix_these2 <- (1:nsites)[is.na(pbdb_sites$modified)];
	fix_these <- sort(unique(c(fix_these1,fix_these2)));
	collection_nos <- pbdb_sites$collection_no[fix_these];
	effed_up_sites <- accersi_collection_data_for_list_of_collection_nos(collection_nos);
	pbdb_sites$created[fix_these] <- effed_up_sites$created;
	pbdb_sites$modified[fix_these] <- effed_up_sites$modified;
	}

#pbdb_taxonomy <- tibble::add_column(pbdb_data_list$pbdb_finds,updater=as.character(rep("",nrow(pbdb_data_list$pbdb_finds))), .before = match("abund_value",colnames(pbdb_data_list$pbdb_finds)));
#pbdb_taxonomy <- tibble::add_column(pbdb_taxonomy,updated=pbdb_taxonomy$modified,.after = match("modified",colnames(pbdb_taxonomy)));
#pbdb_taxonomy <- tibble::add_column(pbdb_taxonomy,updater=pbdb_taxonomy$modifier,.after = match("modifier",colnames(pbdb_taxonomy)));
#pbdb_taxonomy <- tibble::add_column(pbdb_taxonomy,updater_no=pbdb_taxonomy$modifier_no,.after = match("modifier_no",colnames(pbdb_taxonomy)));

pbdb_taxonomy <- rbind(pbdb_taxonomy[!pbdb_taxonomy$taxon_no %in% new_taxa$taxon_no,],new_taxa);
colnames(new_taxa)
pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% pbdb_taxonomy_fixes$taxon_no,] <- pbdb_taxonomy_fixes[pbdb_taxonomy_fixes$taxon_no %in% pbdb_taxonomy$taxon_no,];
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];
pbdb_taxonomy_fixes <- pbdb_taxonomy_fixes[order(pbdb_taxonomy_fixes$taxon_no),];
corrected_phyla_nos <- pbdb_taxonomy_fixes$taxon_no[pbdb_taxonomy_fixes$accepted_rank %in% "phylum"];
ct <- 0;
while (ct < length(corrected_phyla_nos))	{
	ct <- ct+1;
	new_no <- pbdb_taxonomy_fixes$accepted_no[match(corrected_phyla_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$phylum[pbdb_taxonomy$phylum_no %in% corrected_phyla_nos[ct]] <- pbdb_taxonomy_fixes$accepted_name[match(corrected_phyla_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$phylum_no[pbdb_taxonomy$phylum_no %in% corrected_phyla_nos[ct]] <- new_no;
	}
corrected_class_nos <- pbdb_taxonomy_fixes$taxon_no[pbdb_taxonomy_fixes$accepted_rank %in% "class"];
ct <- 0;
while (ct < length(corrected_class_nos))	{
	ct <- ct+1;
	new_no <- pbdb_taxonomy_fixes$accepted_no[match(corrected_class_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$class[pbdb_taxonomy$class_no %in% corrected_class_nos[ct]] <- pbdb_taxonomy_fixes$accepted_name[match(corrected_class_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$class_no[pbdb_taxonomy$class_no %in% corrected_class_nos[ct]] <- new_no;
	}
corrected_order_nos <- pbdb_taxonomy_fixes$taxon_no[pbdb_taxonomy_fixes$accepted_rank %in% "order"];
ct <- 0;
while (ct < length(corrected_order_nos))	{
	ct <- ct+1;
	new_no <- pbdb_taxonomy_fixes$accepted_no[match(corrected_order_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$order[pbdb_taxonomy$order_no %in% corrected_order_nos[ct]] <- pbdb_taxonomy_fixes$accepted_name[match(corrected_order_nos[ct],pbdb_taxonomy_fixes$taxon_no)];
	pbdb_taxonomy$order_no[pbdb_taxonomy$order_no %in% corrected_order_nos[ct]] <- new_no;
	}
#pbdb_taxonomy_fixes$accepted_no[match(corrected_order_nos,pbdb_taxonomy_fixes$taxon_no)]
#pbdb_taxonomy$phylum[pbdb_taxonomy$phylum_no %in% pbdb_taxonomy_fixes$taxon_no] <- pbdb_taxonomy_fixes$accepted_name[pbdb_taxonomy_fixes$taxon_no %in% ]

#pbdb_taxonomy$order_no[pbdb_taxonomy$order %in% "Trepostomata"] <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$taxon_name %in% "Trepostomida"]
#pbdb_taxonomy$order[pbdb_taxonomy$order %in% "Trepostomata"] <- "Trepostomida";
#missing_genera <- pbdb_finds$genus[pbdb_finds$genus %in% ]

# added 2022-11-04
pbdb_finds$subgenus_no[is.na(pbdb_finds$subgenus_no)] <- 0;
pbdb_finds$genus_no[is.na(pbdb_finds$genus_no)] <- 0;
pbdb_finds$genus[is.na(pbdb_finds$genus)] <- "";
missing_genera <- unique(pbdb_finds$genus[!pbdb_finds$genus %in% c("","NO_GENUS_SPECIFIED")][is.na(match(pbdb_finds$genus[!pbdb_finds$genus %in% c("","NO_GENUS_SPECIFIED")],pbdb_taxonomy$taxon_name))]);
missing_genera <- gsub(" ","%20",missing_genera);
mg <- 0;
while (mg < length(missing_genera))	{
	mg <- mg+1;
	milk_carton <- update_taxonomy(taxon=missing_genera[mg],inc_children = F);
	if (ncol(milk_carton)==ncol(pbdb_taxonomy))
		pbdb_taxonomy <- rbind(pbdb_taxonomy,milk_carton);
	}
#genus_rows <- (1:nrow(pbdb_taxonomy))[pbdb_taxonomy$accepted_rank %in% c("genus","subgenus") & pbdb_taxonomy$flags==""]
pbdb_finds$family_no[!pbdb_finds$genus %in% ""] <- as.numeric(pbdb_taxonomy$family_no[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)]);
pbdb_finds$family[!pbdb_finds$genus %in% ""] <- pbdb_taxonomy$family[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)];
pbdb_finds$order_no[!pbdb_finds$genus %in% ""] <- as.numeric(pbdb_taxonomy$order_no[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)]);
pbdb_finds$order[!pbdb_finds$genus %in% ""] <- pbdb_taxonomy$order[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)];
pbdb_finds$class_no[!pbdb_finds$genus %in% ""] <- as.numeric(pbdb_taxonomy$class_no[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)]);
pbdb_finds$class[!pbdb_finds$genus %in% ""] <- pbdb_taxonomy$class[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)];
pbdb_finds$phylum_no[!pbdb_finds$genus %in% ""] <- as.numeric(pbdb_taxonomy$phylum_no[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)]);
pbdb_finds$phylum[!pbdb_finds$genus %in% ""] <- pbdb_taxonomy$phylum[match(pbdb_finds$genus_no[!pbdb_finds$genus %in% ""],pbdb_taxonomy$taxon_no)];
pbdb_finds$family_no[is.na(pbdb_finds$family_no)] <- 0;
pbdb_finds$family[is.na(pbdb_finds$family)] <- "";
pbdb_finds$order_no[is.na(pbdb_finds$order_no)] <- 0;
pbdb_finds$order[is.na(pbdb_finds$order)] <- "";
pbdb_finds$class_no[is.na(pbdb_finds$class_no)] <- 0;
pbdb_finds$class[is.na(pbdb_finds$class)] <- "";
pbdb_finds$phylum_no[is.na(pbdb_finds$phylum_no)] <- 0;
pbdb_finds$phylum[is.na(pbdb_finds$phylum)] <- "";
pbdb_finds <- clean_the_bastards(pbdb_finds);
pbdb_finds$subgenus_no[is.na(pbdb_finds$subgenus_no)] <- 0;
#pbdb_finds$family_no[pbdb_finds$family %in% c(NA,"NO_FAMILY_SPECIFIED")] <- 0;
#pbdb_finds$family[pbdb_finds$family %in% c(NA,"NO_FAMILY_SPECIFIED")] <- "";
#pbdb_finds$order_no[pbdb_finds$order %in% c(NA,"NO_ORDER_SPECIFIED")] <- 0;
#pbdb_finds$order[pbdb_finds$order %in% c(NA,"NO_ORDER_SPECIFIED")] <- "";
#pbdb_finds$class_no[pbdb_finds$class %in% c(NA,"NO_CLASS_SPECIFIED")] <- 0;
#pbdb_finds$class[pbdb_finds$class %in% c(NA,"NO_CLASS_SPECIFIED")] <- "";
#pbdb_finds$phylum_no[pbdb_finds$phylum %in% c(NA,"NO_PHYLUM_SPECIFIED")] <- 0;
#pbdb_finds$phylum[pbdb_finds$phylum %in% c(NA,"NO_PHYLUM_SPECIFIED")] <- "";

edited_opinions <- new_opinions[(1:nrow(new_opinions))[new_opinions$opinion_no %in% pbdb_opinions$opinion_no],];
added_opinions <- new_opinions[(1:nrow(new_opinions))[!new_opinions$opinion_no %in% pbdb_opinions$opinion_no],];
#colnames(pbdb_opinions)
#pbdb_opinions <- tibble::add_column(pbdb_opinions,updated=pbdb_opinions$modified,.after = match("modified",colnames(pbdb_opinions)));
#pbdb_opinions <- tibble::add_column(pbdb_opinions,updater=pbdb_opinions$modifier,.after = match("modifier",colnames(pbdb_opinions)));
#pbdb_opinions <- tibble::add_column(pbdb_opinions,updater_no=pbdb_opinions$modifier_no,.after = match("modifier_no",colnames(pbdb_opinions)));

pbdb_opinions[match(edited_opinions$opinion_no,pbdb_opinions$opinion_no),] <- edited_opinions;
pbdb_opinions <- rbind(pbdb_opinions,added_opinions);
pbdb_opinions <- pbdb_opinions[order(pbdb_opinions$opinion_no),];

# generate PBDB rock unit database ####
print("Generating a database of rock units in the PBDB.");
print("Cleaning stupid rock names yet AGAIN...")
named_rock_unit <- pbdb_sites$formation[pbdb_sites$formation!=""];
pbdb_sites$formation[pbdb_sites$formation!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$member[pbdb_sites$member!=""];
pbdb_sites$member[pbdb_sites$member!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""];
pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$formation_alt[pbdb_sites$formation_alt!=""];
pbdb_sites$formation_alt[pbdb_sites$formation_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$member_alt[pbdb_sites$member_alt!=""];
pbdb_sites$member_alt[pbdb_sites$member_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
named_rock_unit <- pbdb_sites$stratgroup_alt[pbdb_sites$stratgroup_alt!=""];
pbdb_sites$stratgroup_alt[pbdb_sites$stratgroup_alt!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);

# some legitimate rock names get eradicated this way; put them back if you can!
pbdb_sites <- reparo_unedittable_paleodb_collections(pbdb_sites,paleodb_collection_edits = paleodb_fixes$paleodb_collection_edits);

pbdb_sites_w_rocks <- unique(rbind(subset(pbdb_sites,pbdb_sites$formation!=""),
								   subset(pbdb_sites,pbdb_sites$stratgroup!=""),
								   subset(pbdb_sites,pbdb_sites$member!="")));
pbdb_sites_w_rocks <- pbdb_sites_w_rocks[order(pbdb_sites_w_rocks$collection_no),];
pbdb_rock_info <- organize_pbdb_rock_data(paleodb_collections = pbdb_sites_w_rocks,geosplit = T);
pbdb_rocks <- pbdb_rock_info$pbdb_rocks;
#write.csv(pbdb_rock_info$pbdb_rocks,"PBDB_Rocks2.csv",row.names = F)
rock_names <- pbdb_rocks$rock_unit_clean_no_rock;
nrocks <- nrow(pbdb_rocks);
group_only <- (1:nrocks)[pbdb_rocks$formation=="" & pbdb_rocks$member==""];
last_group <- (min((1:nrocks)[pbdb_rocks$formation!=""])-1);
formation_names <- c(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],
					 sort(pbdb_rocks$formation_clean_no_rock[pbdb_rocks$formation_clean_no_rock!=""]));
rock_no_orig <- rock_no <- 1:nrocks;
rock_no_sr <- match(pbdb_rocks$rock_unit_clean_no_rock,rock_names);
formation_no <- rock_no_sr[match(pbdb_rocks$formation_clean_no_rock,formation_names)];
formation_no[1:last_group] <- rock_no_sr[match(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],formation_names)];
member_only <- (1:nrocks)[pbdb_rocks$formation=="" & pbdb_rocks$member!=""];
formations_and_members <- member_only[pbdb_rocks$member_clean_no_rock[member_only] %in% formation_names]
formation_no[formations_and_members] <- formation_no[match(pbdb_rocks$member_clean_no_rock[formations_and_members],formation_names)]
pbdb_rocks <- tibble::add_column(pbdb_rocks, formation_no=as.numeric(formation_no), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no_sr=as.numeric(rock_no_sr), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no=as.numeric(rock_no), .before = 1);
pbdb_rocks <- tibble::add_column(pbdb_rocks, rock_no_orig=as.numeric(rock_no_orig), .before = 1);
formation_no_dummy <- match(pbdb_rocks$formation_clean_no_rock,rock_names);
#formation_no_dummy <- match(pbdb_rocks$formation_clean_no_rock,pbdb_rocks$rock_unit_clean_no_rock);
formation_no_dummy[1:last_group] <- match(pbdb_rocks$rock_unit_clean_no_rock[1:last_group],formation_names);
#
add_here <- (1:nrow(pbdb_rocks))[is.na(formation_no_dummy)];
add_formations <- pbdb_rocks$formation_clean_no_rock[is.na(formation_no_dummy)];
#match("Ablakoskovolgy",pbdb_rocks$)
kluge_city <- cbind(add_here,add_formations);#add_formations <- unique(pbdb_rocks$formation_clean_no_rock[is.na(formation_no_dummy)])
kluge_city <- kluge_city[match(unique(add_formations),add_formations),];
add_here <- as.numeric(kluge_city[,1]);
a_f <- length(add_here);
pbdb_rocks_old <- pbdb_rocks;
for (af in a_f:1)	{
	if (af>1)
		while (af>1 & add_here[af]==1+add_here[af-1] & pbdb_rocks$formation_clean_basic[add_here[af]]==pbdb_rocks$formation_clean_basic[add_here[af-1]])
			af <- af-1;
	nn <- add_here[af];
	add_rock <- pbdb_rocks[nn,];
	add_rock$rock_no_orig <- -1;
	add_rock$rock_no <- add_rock$rock_no-0.5;
	add_rock$rock_no_sr <- add_rock$rock_no_sr-0.5;
	add_rock$member <- add_rock$member_clean_basic <- add_rock$member_clean_no_rock <- add_rock$member_clean_no_rock_formal <- "";
	add_rock$full_name <- add_rock$formation;
	add_rock$rock_unit_clean_basic <- add_rock$formation_clean_basic;
	add_rock$rock_unit_clean_no_rock <- add_rock$formation_clean_no_rock;
	add_rock$rock_unit_clean_no_rock_formal <- add_rock$formation_clean_no_rock_formal;
	yin <- 1:(nn-1);
	yang  <- nn:nrow(pbdb_rocks);
	pbdb_rocks <- rbind(pbdb_rocks[yin,],add_rock,pbdb_rocks[yang,]);
#	pbdb_rocks[(nn-5):(nn+2),];
	}
nrocks <- nrow(pbdb_rocks);
pbdb_rocks_redone <- pbdb_rocks;
#nn <- match("Ablakoskovolgy",pbdb_rocks$formation_clean_no_rock)
new_rock_no <- 1:nrocks;
new_rock_no_sr <- new_rock_no[match(pbdb_rocks$rock_no_sr,pbdb_rocks$rock_no_sr)];
new_formation_no <- new_rock_no_sr[match(pbdb_rocks$formation_no,pbdb_rocks$formation_no)];
pbdb_rocks_redone$rock_no <- new_rock_no;
pbdb_rocks_redone$rock_no_sr <- new_rock_no_sr;
pbdb_rocks_redone$formation_no <- new_formation_no;

pbdb_rocks_site_lists <- pbdb_rock_info$site_list;

# eliminate accidental duplicate sites 1####
print("We will now pause to remove replicates from the data...");
dup_colls <- hist(pbdb_sites$collection_no,breaks=0:max(pbdb_sites$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites[pbdb_sites$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
#	if (nrow(dup_sites)>1 && max(dup_sites$pbdb_rock_no)>0)
#		dup_sites <- dup_sites[dup_sites$pbdb_rock_no>0,];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites[pbdb_sites$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites <- unique(pbdb_sites);

pbdb_sites_refined <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_sites_refined);
new_sites_modified <- pbdb_sites[pbdb_sites$collection_no %in% new_sites$collection_no,];
new_fields <- colnames(pbdb_sites_refined)[!colnames(pbdb_sites_refined) %in% colnames(new_sites_modified)];
dummy <- pbdb_sites_refined[1:nrow(new_sites_modified),match(new_fields,colnames(pbdb_sites_refined))];
nf <- 0;
while (nf < length(new_fields))	{
	nf <- nf+1;
	psrc <- match(new_fields[nf],colnames(dummy));
	if (is.numeric(dummy[,psrc]))	{
		dummy[,psrc] <- 0;
		} else	{
		dummy[,psrc] <- "";
		}
	}

new_sites_modified <- cbind(new_sites_modified,dummy);
added_sites_refn <- new_sites_modified[!new_sites_modified$collection_no %in% pbdb_sites_refined$collection_no,];
edited_sites_refn <- new_sites_modified[new_sites_modified$collection_no %in% pbdb_sites_refined$collection_no,];
rr <- match(edited_sites_refn$collection_no,pbdb_sites_refined$collection_no);
#match(colnames(new_sites_modified),colnames(pbdb_sites_refined))
pbdb_sites_refined[rr,] <- edited_sites_refn;
pbdb_sites_refined <- rbind(pbdb_sites_refined,added_sites_refn);

if (!is.null(pbdb_sites_refined))	{
dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[(dup_sites$ma_lb-dup_sites$ma_ub)==min(dup_sites$ma_lb-dup_sites$ma_ub),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);
}

cd <- unique(pbdb_finds$collection_no[!pbdb_finds$collection_no %in% pbdb_sites_refined$collection_no]);
# Refine Rocks by Period ####
print("Refining information about rock units in each geological period.");
effed <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
#sum(pbdb_sites_refined$interval_lb[effed]==pbdb_sites_refined$interval_ub[effed])
effed_int <- pbdb_sites_refined$interval_lb[effed];
pbdb_sites_refined$interval_lb[effed] <- pbdb_sites_refined$interval_ub[effed];
pbdb_sites_refined$interval_ub[effed] <- effed_int;
pbdb_sites_refined$ma_lb[effed] <- finest_chronostrat$ma_lb[match(pbdb_sites_refined$interval_lb[effed],finest_chronostrat$interval)];
pbdb_sites_refined$ma_ub[effed] <- finest_chronostrat$ma_ub[match(pbdb_sites_refined$interval_ub[effed],finest_chronostrat$interval)];

age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
zone_database <- gradstein_2020_emended$zones;
chronostrat_units <- unique(c(pbdb_sites$early_interval,pbdb_sites$late_interval));
hierarchical_chronostrat <- accersi_hierarchical_timescale(chronostrat_units=chronostrat_units,time_scale=gradstein_2020_emended$time_scale,regional_scale="Stage Slice");
myrs <- sort(unique(round(c(hierarchical_chronostrat$ma_lb,hierarchical_chronostrat$ma_ub),5)),decreasing=T)
hierarchical_chronostrat$bin_first <- match(hierarchical_chronostrat$ma_lb,myrs);
hierarchical_chronostrat$bin_last <- match(hierarchical_chronostrat$ma_ub,myrs)-1;
#write.csv(hierarchical_chronostrat,"Hierarchical_Chronostrat.csv",row.names = F);
#hierarchical_chronostrat$bin_last <- match(hierarchical_chronostrat$bin_last,sort(unique(c(hierarchical_chronostrat$bin_first,hierarchical_chronostrat$bin_last))));
old_ma_mids <- round((pbdb_sites_refined$ma_lb+pbdb_sites_refined$ma_ub)/2,0);
opt_periods <- c("Orosirian","Statherian","Calymmian","Ectasian","Stenian","Tonian","Cryogenian","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene");
#hierarchical_chronostrat$interval[hierarchical_chronostrat$interval!=hierarchical_chronostrat$st] <- hierarchical_chronostrat$st[hierarchical_chronostrat$interval!=hierarchical_chronostrat$st];

for (b1 in 1:(length(opt_periods)))	{
	print(paste("doing the",opt_periods[b1]));
	max_ma <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b1],gradstein_2020_emended$time_scale$interval)];
	if (b1<length(opt_periods))	{
		min_ma <- gradstein_2020_emended$time_scale$ma_ub[match(opt_periods[b1],gradstein_2020_emended$time_scale$interval)];
		} else	{
		min_ma <- 0;
		}
	interval_sites <- pbdb_sites[pbdb_sites$max_ma>=min_ma & pbdb_sites$min_ma<=max_ma,];

	rock_database <- rock_unit_data$rock_unit_database[rock_unit_data$rock_unit_database$ma_ub<(max_ma+25) & rock_unit_data$rock_unit_database$ma_lb>(min_ma-25),];
	rock_to_zone_database <- rock_unit_data$rock_to_zone_database[rock_unit_data$rock_to_zone_database$rock_no %in% rock_database$rock_no,];
	zone_database <- gradstein_2020_emended$zone[gradstein_2020_emended$zone$ma_ub<=(max_ma+25) & gradstein_2020_emended$zone$ma_lb>=(min_ma-25),];
	zone_database <- zone_database[!is.na(zone_database$ma_lb),];

	refined_info <- refine_pbdb_collections_w_external_databases(paleodb_collections=interval_sites,rock_database,zone_database,rock_to_zone_database,finest_chronostrat);
	#colnames(refined_info$refined_collections)[!colnames(refined_info$refined_collections) %in% colnames(pbdb_sites)]
	new_sites_ref <- refined_info$refined_collections[!refined_info$refined_collections$collection_no %in% pbdb_sites_refined$collection_no,];
	old_sites_ref <- refined_info$refined_collections[refined_info$refined_collections$collection_no %in% pbdb_sites_refined$collection_no,];
	old_sites_ref <- old_sites_ref[order(old_sites_ref$collection_no),];
#	old_sites$rock_no[old_sites$formation=="Simeh–Kuh"] <- old_sites$rock_no_sr[old_sites$formation=="Simeh–Kuh"] <- old_sites$formation_no[old_sites$formation=="Simeh–Kuh"] <- 8995;
#		now_rocked <- old_sites$collection_no[old_sites$rock_unit_senior!=""];
#		was_rocked <- pbdb_sites_refined$collection_no[pbdb_sites_refined$rock_unit_senior!=""];
#		rock_change <- now_rocked[!now_rocked %in% was_rocked]
#	now_rockless <- old_sites$collection_no[!old_sites$rock_unit_senior!=""];
#	was_rockless <- pbdb_sites_refined$collection_no[!pbdb_sites_refined$rock_unit_senior!=""];
#	change_rock <- now_rockless[!now_rockless %in% was_rockless]
#	old_sites <- old_sites[!old_sites$collection_no %in% change_rock,];
#		pbdb_sites_refined$rock_unit_senior[pbdb_sites_refined$collection_no %in% change_rock]
		# update sites
#	old_sites <- old_sites[,colnames(old_sites) %in% colnames(pbdb_sites_refined)];
	# only update if there is now rock info!!!!
#	newly_rocked_collections <- old_sites$collection_no[old_sites$rock_no_sr>0];
#	old_sites_rocked <- subset(old_sites,old_sites$rock_no>0);
	# update old sites
	pbdb_sites_refined[match(old_sites_ref$collection_no,pbdb_sites_refined$collection_no),match(colnames(old_sites_ref),colnames(pbdb_sites_refined))] <- old_sites_ref[,colnames(old_sites_ref) %in% colnames(pbdb_sites_refined)];
	pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
#	pbdb_sites_refined[match(old_sites_rocked$collection_no,pbdb_sites_refined$collection_no),match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))] <- old_sites_rocked[,colnames(old_sites_rocked) %in% colnames(pbdb_sites_refined)];
#	pbdb_sites_refined[match(old_sites_rocked$collection_no,pbdb_sites_refined$collection_no),match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))] <- old_sites_rocked[,match(colnames(old_sites_rocked),colnames(pbdb_sites_refined))];
#	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% old_sites$collection_no,] <- old_sites;

	# add new sites
	new_sites_ref <- new_sites_ref[,colnames(new_sites_ref) %in% colnames(pbdb_sites_refined)];
	if (nrow(new_sites_ref)>0)	{
		if (ncol(new_sites_ref)<ncol(pbdb_sites_refined))	{
			lost_fields <- colnames(pbdb_sites_refined)[!colnames(pbdb_sites_refined) %in% colnames(new_sites_ref)];
			dummy <- array(0,dim=c(nrow(new_sites_ref),length(lost_fields)));
			colnames(dummy) <- lost_fields;
			new_sites_ref <- cbind(new_sites_ref,dummy);
			}
#		if (ncol(pbdb_sites_refined)<=ncol(new_sites))	{
		pbdb_sites_refined <- rbind(pbdb_sites_refined,new_sites_ref[,match(colnames(new_sites_ref),colnames(pbdb_sites_refined))]);
#			}
		}
	pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
	}

dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1 && max(dup_sites$rock_no)>0 && !is.infinite(abs(dup_sites$rock_no)))
		dup_sites <- dup_sites[dup_sites$rock_no>0,];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);
beepr::beep("wilhelm");

write.csv(pbdb_sites_refined,"PBDB_Sites_Refined.csv",row.names = F,fileEncoding='UTF-8');
writexl::write_xlsx(pbdb_sites_refined,"PBDB_Sites_Refined.xlsx");

# update collection ages based on rocks if there is disagreement ####
#pbdb_sites_refined <- as.data.frame(readxl::read_xlsx("PBDB_Sites_Refined.xlsx"));
#pbdb_sites_refined <- clean_pbdb_fields(pbdb_data=pbdb_sites_refined);
pbdb_sites_refined <- put_pbdb_dataframes_into_proper_type(pbdb_sites_refined);

nsites <- nrow(pbdb_sites_refined);
effed <- (1:nsites)[is.na(pbdb_sites_refined$ma_lb)];
if (length(effed)>0)	{
	pbdb_sites_refined$ma_lb[effed] <- time_scale$ma_lb[match(pbdb_sites_refined$early_interval[effed],time_scale$interval)];
	age <- pbdb_sites_refined$ma_lb[effed];
	pbdb_sites_refined$interval_lb[effed] <- sapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
	pbdb_sites_refined$ma_ub[effed] <- time_scale$ma_ub[match(pbdb_sites_refined$late_interval[effed],time_scale$interval)];
	age <- pbdb_sites_refined$ma_ub[effed];
	pbdb_sites_refined$interval_ub[effed] <- sapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
	}
effed <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];

#sum(pbdb_sites_refined$interval_lb[effed]==pbdb_sites_refined$interval_ub[effed])
effed_int <- pbdb_sites_refined$interval_lb[effed];
pbdb_sites_refined$interval_lb[effed] <- pbdb_sites_refined$interval_ub[effed];
pbdb_sites_refined$interval_ub[effed] <- effed_int;
pbdb_sites_refined$ma_lb[effed] <- finest_chronostrat$ma_lb[match(pbdb_sites_refined$interval_lb[effed],finest_chronostrat$interval)];
pbdb_sites_refined$ma_ub[effed] <- finest_chronostrat$ma_ub[match(pbdb_sites_refined$interval_ub[effed],finest_chronostrat$interval)];

rock_database <- rock_unit_data$rock_unit_database;
zone_database <- gradstein_2020_emended$zones;
#cn <- cbind(pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$rock_no>0],rock_database$rock_no_sr[match(pbdb_sites_refined$rock_no[pbdb_sites_refined$rock_no>0],rock_database$rock_no)]);
#write.csv(cn,"Rock_Numbers.csv",row.names = F);
#pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$rock_no>0] <- rock_database$rock_no_sr[match(pbdb_sites_refined$rock_no[pbdb_sites_refined$rock_no>0],rock_database$rock_no)];
#pbdb_sites_refined$rock_no_sr[pbdb_sites_refined$collection_no %in% c(155912,155918)] <- pbdb_sites_refined$rock_no[pbdb_sites_refined$collection_no %in% c(155912,155918)] <- 3461;
#pbdb_sites_refined$formation[pbdb_sites_refined$collection_no %in% c(155912,155918)] <- "Dundee Flagstone";
#pbdb_sites_refined$zone[pbdb_sites_refined$collection_no %in% c(155912,155918)] <- "Emphanisporites micronatus - Streelispora newportensis";
#pbdb_sites_refined$stratgroup[pbdb_sites_refined$collection_no %in% c(155912,155918)] <- "Arbuthnott-Garvock";
#pbdb_sites$formation[pbdb_sites$collection_no %in% c(155912,155918)] <- "Dundee Flagstone";
#pbdb_sites$zone[pbdb_sites$collection_no %in% c(155912,155918)] <- "Emphanisporites micronatus - Streelispora newportensis";
#pbdb_sites$stratgroup[pbdb_sites$collection_no %in% c(155912,155918)] <- "Arbuthnott-Garvock";

rocked_sites <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$rock_no_sr>0];
# update sites with overly old lower bounds
effed1 <- rocked_sites[pbdb_sites_refined$ma_lb[rocked_sites]>rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[rocked_sites],rock_database$rock_no)]];
ef <- 0;
#rock_database[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed1[ef]],];
#rock_database$ma_lb[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed1[ef]]] <- 61.66;
#rock_database$ma_ub[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed1[ef]]] <- 57.55;
while (ef < length(effed1))	{
	ef <- ef+1;
	overlap <- accersi_temporal_overlap(lb1=pbdb_sites_refined$ma_lb[effed1[ef]],
										ub1=pbdb_sites_refined$ma_ub[effed1[ef]],
										lb2=rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed1[ef]],rock_database$rock_no)],
										ub2=rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[effed1[ef]],rock_database$rock_no)]);
	if (overlap$ma_lb>0 || overlap$ma_ub>0)	{
		pbdb_sites_refined$ma_lb[effed1[ef]] <- overlap$ma_lb;
		pbdb_sites_refined$ma_ub[effed1[ef]] <- overlap$ma_ub;
		}
	}
effed1 <- rocked_sites[pbdb_sites_refined$ma_lb[rocked_sites]>rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[rocked_sites],rock_database$rock_no)]];
test1 <- pbdb_sites_refined$ma_lb[effed1]-rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed1],rock_database$rock_no)];
too_long <- 5;
if (length(test1[test1<too_long])>0)	{
	# correct outdated dates in pbdb_sites_refined
	pbdb_sites_refined$ma_lb[effed1] <- time_scale$ma_lb[match(pbdb_sites_refined$early_interval[effed1],time_scale$interval)];
	effed1_up <- effed1[pbdb_sites_refined$ma_ub[effed1]>=time_scale$ma_lb[match(pbdb_sites_refined$early_interval[effed1],time_scale$interval)]];
	pbdb_sites_refined$ma_ub[effed1_up] <- time_scale$ma_ub[match(pbdb_sites_refined$late_interval[effed1_up],time_scale$interval)];
	}

effed1 <- rocked_sites[pbdb_sites_refined$ma_lb[rocked_sites]>rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[rocked_sites],rock_database$rock_no)]];
test1 <- pbdb_sites_refined$ma_lb[effed1]-rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed1],rock_database$rock_no)];
if (length(test1)>0)	{
	dummy <- pbdb_sites_refined[effed1,];
	dummy[,!colnames(dummy) %in% c("collection_no","rock_unit_senior","rock_no_sr","ma_lb","ma_ub")] <- NULL;
	dummy <- dummy[order(dummy$rock_no_sr),];
	dummy$ma_lb_rock <- rock_database$ma_lb[match(dummy$rock_no_sr,rock_database$rock_no)];
	dummy$ma_ub_rock <- rock_database$ma_ub[match(dummy$rock_no_sr,rock_database$rock_no)];
	write.csv(dummy,"Collections_w_Rocks_Needing_Older_Dates.csv",row.names = F);
	}
#pbdb_sites_refined$collection_no[effed1[8]]
#rock_database[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed1[8]],]
pbdb_sites_refined$ma_lb[effed1] <- rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed1],rock_database$rock_no)];
effed1_up <- effed1[pbdb_sites_refined$ma_ub[effed1]>=time_scale$ma_lb[match(pbdb_sites_refined$early_interval[effed1],time_scale$interval)]];
pbdb_sites_refined$ma_ub[effed1_up] <- time_scale$ma_ub[match(pbdb_sites_refined$late_interval[effed1_up],time_scale$interval)];

# update sites with overly young upper bounds
effed2 <- rocked_sites[pbdb_sites_refined$ma_ub[rocked_sites]<rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[rocked_sites],rock_database$rock_no)]];
#rock_database[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed2[ef]],];
rock_database[rock_database$rock_no %in% pbdb_sites_refined$rock_no_sr[effed2],];
#rock_database$ma_lb[rock_database$rock_no_sr==5402] <- 473.5;
#rock_database$ma_ub[rock_database$rock_no_sr==5402] <- 470.5;
#rock_database$ma_lb[rock_database$rock_no_sr==21195] <- 59.24;
#rock_database$ma_ub[rock_database$rock_no_sr==21195] <- 55.1;
#rock_database$ma_lb[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed2[ef]]] <- 384.1;
#rock_database$ma_ub[rock_database$rock_no==pbdb_sites_refined$rock_no_sr[effed2[ef]]] <- 382.1;
ef <- 0;
while (ef < length(effed2))	{
	ef <- ef+1;
	overlap <- accersi_temporal_overlap(lb1=pbdb_sites_refined$ma_lb[effed2[ef]],
										ub1=pbdb_sites_refined$ma_ub[effed2[ef]],
										lb2=rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed2[ef]],rock_database$rock_no)],
										ub2=rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[effed2[ef]],rock_database$rock_no)]);
	if (overlap$ma_lb>0 || overlap$ma_ub>0)	{
		pbdb_sites_refined$ma_lb[effed2[ef]] <- overlap$ma_lb;
		pbdb_sites_refined$ma_ub[effed2[ef]] <- overlap$ma_ub;
		}
	}
effed2 <- rocked_sites[pbdb_sites_refined$ma_ub[rocked_sites]<rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[rocked_sites],rock_database$rock_no)]];
#cbind(pbdb_sites_refined$ma_ub[effed2],rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[effed2],rock_database$rock_no)]);
test2 <- rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[effed2],rock_database$rock_no)]-pbdb_sites_refined$ma_ub[effed2];
#paste(pbdb_sites_refined$collection_no[effed2],collapse=",")
if (length(test2[test2>too_long])>0)	{
	dummy <- pbdb_sites_refined[effed2,];
	dummy[,!colnames(dummy) %in% c("collection_no","rock_unit_senior","rock_no_sr","ma_lb","ma_ub")] <- NULL;
	dummy <- dummy[order(dummy$rock_no_sr),];
	dummy$ma_lb_rock <- rock_database$ma_lb[match(dummy$rock_no_sr,rock_database$rock_no)];
	dummy$ma_ub_rock <- rock_database$ma_ub[match(dummy$rock_no_sr,rock_database$rock_no)];
	write.csv(dummy,"Collections_w_Rocks_Needing_Younger_Dates.csv",row.names = F);
	}
effed2 <- effed2[test2<=too_long];
pbdb_sites_refined$ma_ub[effed2] <- rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[effed2],rock_database$rock_no)];
effed2_lo <- effed2[pbdb_sites_refined$ma_lb[effed2]<=pbdb_sites_refined$ma_ub[effed2]];
ef <- 0;
#rock_database[rock_database$rock_no_sr %in% pbdb_sites_refined$rock_no_sr[effed2_lo],]
#rock_database$ma_lb[rock_database$rock_no_sr %in% pbdb_sites_refined$rock_no_sr[effed2_lo]] <- 276.9;
while (ef < length(effed2_lo))	{
	ef <- ef+1;
	poss_fas <- c(time_scale$ma_lb[match(pbdb_sites_refined$late_interval[effed2_lo[ef]],time_scale$interval)],
				  rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[effed2_lo[ef]],rock_database$rock_no)]);
	poss_fas <- poss_fas[poss_fas>pbdb_sites_refined$ma_ub[effed2_lo[ef]]];
	pbdb_sites_refined$ma_lb[effed2_lo[ef]] <- min(poss_fas);
	}
#hist(test2)

effed <- (1:nrow(pbdb_sites_refined))[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
#sum(pbdb_sites_refined$interval_lb[effed]==pbdb_sites_refined$interval_ub[effed])
effed_int <- pbdb_sites_refined$interval_lb[effed];
pbdb_sites_refined$interval_lb[effed] <- pbdb_sites_refined$interval_ub[effed];
pbdb_sites_refined$interval_ub[effed] <- effed_int;
pbdb_sites_refined$ma_lb[effed] <- finest_chronostrat$ma_lb[match(pbdb_sites_refined$interval_lb[effed],finest_chronostrat$interval)];
pbdb_sites_refined$ma_ub[effed] <- finest_chronostrat$ma_ub[match(pbdb_sites_refined$interval_ub[effed],finest_chronostrat$interval)];

# add PBDB Rock numbers ####
print("Add numbers to rock units from recently created database.");
p_r_s_l <- length(pbdb_rocks_site_lists);
nsites <- nrow(pbdb_sites_refined);
pbdb_sites_refined$pbdb_formation_no <- pbdb_sites_refined$pbdb_rock_no_sr <- pbdb_sites_refined$pbdb_rock_no <- rep(0,nsites);
#pbdb_rocks_refined$rock_no <- 1:nrow(pbdb_rocks);
for (pr in 1:p_r_s_l)	{
	rd <- match(pr,pbdb_rocks_redone$rock_no_orig);
	pbdb_rows <- match(pbdb_rocks_site_lists[[pr]],pbdb_sites_refined$collection_no);
	pbdb_sites_refined$pbdb_rock_no[pbdb_rows] <- pbdb_rocks_redone$rock_no[rd];
	pbdb_sites_refined$pbdb_rock_no_sr[pbdb_rows] <- pbdb_rocks_redone$rock_no_sr[rd];
	pbdb_sites_refined$pbdb_formation_no[pbdb_rows] <- pbdb_rocks_redone$formation_no[rd];
	}
beepr::beep("wilhelm");
# deal with duds ####
nsites <- nrow(pbdb_sites_refined);
dud_sites <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];

# Optimize by Period ####
print("Optimize localities using limited Unitary Association.");
pbdb_sites_refined_orig <- pbdb_sites_refined;
single_binners <- pbdb_data_list$single_binners;
single_binners_a <- single_binners_z <- c();
age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);
opt_periods <- opt_periods[match("Ediacaran",opt_periods):length(opt_periods)];
single_binners_a <- vector(length=length(opt_periods)-1);
for (b2 in 1:(length(opt_periods)-1))	{
	print(paste("Optimizing the",opt_periods[b2],"+",opt_periods[b2+1]));
	max_ma <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b2],gradstein_2020_emended$time_scale$interval)];
	min_ma <- gradstein_2020_emended$time_scale$ma_ub[match(opt_periods[b2+1],gradstein_2020_emended$time_scale$interval)];
#	max_ma_b <- gradstein_2020_emended$time_scale$ma_lb[match(opt_periods[b2-1],gradstein_2020_emended$time_scale$interval)];
	interval_sites <- pbdb_sites_refined[pbdb_sites_refined$ma_lb<=max_ma & pbdb_sites_refined$ma_ub>=min_ma,];
	single_binners_a[b2] <- sum(interval_sites$interval_lb==interval_sites$interval_ub);
#	interval_sites_b <- pbdb_sites_refined[pbdb_sites_refined$ma_lb<=max_ma_b & pbdb_sites_refined$ma_ub>=min_ma,];
#	interval_finds <- pbdb_finds[pbdb_finds$collection_no %in% interval_sites_b$collection_no,];
	interval_finds <- pbdb_finds[pbdb_finds$collection_no %in% interval_sites$collection_no,];
	interval_finds <- interval_finds[interval_finds$identified_rank %in% c("species","subspecies"),];
	if (is.null(interval_sites$ma_lb))	interval_sites$ma_lb <- interval_sites$max_ma;
	if (is.null(interval_sites$ma_ub))	interval_sites$ma_ub <- interval_sites$min_ma;
	interval_zones <- gradstein_2020_emended$zone[gradstein_2020_emended$zone$ma_ub<=(max_ma+25) & gradstein_2020_emended$zone$ma_lb>=(min_ma-25),];
	options(warn=-1);
	#write.csv(interval_sites,"Ediacaran-Cambrian_Sites.csv",row.names = F)
	interval_finds$early_interval <- interval_sites$early_interval[match(interval_finds$collection_no,interval_sites$collection_no)];
	interval_finds$late_interval <- interval_sites$late_interval[match(interval_finds$collection_no,interval_sites$collection_no)];
	stages <- sort(unique(c(interval_finds$early_interval,interval_finds$late_interval)));
	#match(stages,time_scale$interval)
	# paleodb_finds=interval_finds;paleodb_collections=interval_sites;hierarchical_chronostrat=hierarchical_chronostrat;zone_database=interval_zones;
	zone_database <- gradstein_2020_emended$zones;
	optimized_sites <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds=interval_finds,paleodb_collections=interval_sites,hierarchical_chronostrat=hierarchical_chronostrat,zone_database=interval_zones);
	#write.csv(optimized_sites,"Ediacaran-Cambrian_Sites_Refined.csv",row.names = F)
#	keep_these <- match(colnames(optimized_sites)[colnames(optimized_sites) %in% colnames(pbdb_sites)],colnames(optimized_sites));
#	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no,] <- optimized_sites[,keep_these];
	single_binners_z <- c(single_binners_z,sum(optimized_sites$interval_lb==optimized_sites$interval_ub));
#	pbdb_sites_refined$ma_lb[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no] <- optimized_sites$ma_lb;
	pbdb_sites_refined[pbdb_sites_refined$collection_no %in% optimized_sites$collection_no,] <- optimized_sites;
	options(warn=1);
	}
if (ncol(single_binners)-1 < 10)	{
	new_singles <- paste("single_binners_0",ncol(single_binners)-1,sep="");
	} else	{
	new_singles <- paste("single_binners_",ncol(single_binners)-1,sep="");
	}
#pbdb_sites_refined[pbdb_sites_refined$collection_no==1441,]

age <- pbdb_sites_refined$ma_lb;
pbdb_sites_refined$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"onset",finest_chronostrat);
age <- pbdb_sites_refined$ma_ub;
pbdb_sites_refined$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,"end",finest_chronostrat);

single_binners <- cbind(single_binners,single_binners_z);
colnames(single_binners)[ncol(single_binners)] <- new_singles;
write.csv(single_binners,"Single_Binners.csv",row.names = F);
beepr::beep("wilhelm");

# deal with duds 2 ####
nsites <- nrow(pbdb_sites_refined);
dud_sites <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];

# eliminate accidental duplicate sites 2####
print("Eliminate any duplicate information.");
dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
names(dup_colls) <- 1:max(pbdb_sites_refined$collection_no);
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
#	dup_sites$created <- dup_sites$modified <- as.Date("2021-03-07")
	if (nrow(dup_sites)>1)	dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1 && max(dup_sites$rock_no)>0)	dup_sites <- dup_sites[dup_sites$rock_no>0,];
	if (nrow(dup_sites)>1)	dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

dup_colls <- hist(pbdb_sites_refined$collection_no,breaks=0:max(pbdb_sites_refined$collection_no),plot=F)$counts;
#names(dup_colls) <- (1:max(pbdb_sites_refined$collection_no));
dup_colls <- dup_colls[dup_colls>1];
dc <- 0;
while (dc < length(dup_colls))	{
	dc <- dc+1;
	dup_sites_init <- pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),];
	dup_sites <- unique(dup_sites_init);
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$modified==max(dup_sites$modified),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[(dup_sites$ma_lb-dup_sites$ma_ub)==min(dup_sites$ma_lb-dup_sites$ma_ub),];
	if (nrow(dup_sites)>1)
		dup_sites <- dup_sites[dup_sites$n_occs==max(dup_sites$n_occs),];
	if (nrow(dup_sites)==1)	{
		for (ds in 1:nrow(dup_sites_init))	dup_sites_init[ds,] <- dup_sites;
		pbdb_sites_refined[pbdb_sites_refined$collection_no %in% as.numeric(names(dup_colls)[dc]),] <- dup_sites_init;
		} else	{
		(1:ncol(dup_sites))[dup_sites[1,]!=dup_sites[2,]];
		print(dc);
		}
	}
if (length(dup_colls)>0)	pbdb_sites_refined <- unique(pbdb_sites_refined);

# I hate that I have to do this 2!!! ####
print("Make sure that dates are dates & not numbers.");
occurrence_no <- pbdb_finds$occurrence_no[is.na(pbdb_finds$modified)];
fixed_finds <- data.frame(base::t(pbapply::pbsapply(occurrence_no,update_occurrence_from_occurrence_no)));
new_modified <- unlist(fixed_finds$modified);
if (sum(new_modified=="0000-00-00 00:00:00")>0)	{
	new_created <- unlist(fixed_finds$created);
	new_modified[new_modified=="0000-00-00 00:00:00"] <- new_created[new_modified=="0000-00-00 00:00:00"];
	}

pbdb_finds$modified[is.na(pbdb_finds$modified)] <- as.character.Date(new_modified);
#pbdb_finds$modified[pbdb_finds$occurrence_no %in% occurrence_no];

# put a bow on it and gift it to yourself ####
print("You deserve this....");
pbdb_data_list$pbdb_finds <- pbdb_finds;
pbdb_data_list$pbdb_sites <- pbdb_sites;
pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
pbdb_data_list$pbdb_opinions <- pbdb_opinions;
pbdb_data_list$time_scale <- finest_chronostrat;
pbdb_data_list$pbdb_rocks <- pbdb_rocks_redone;
pbdb_data_list$pbdb_rocks_sites <- pbdb_rocks_sites;
pbdb_data_list$single_binners <- single_binners;
#save(pbdb_data_list,file=paste(getwd(),"/data/Paleobiology_Database.RData",sep=""));
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Marine_Tetrapod_Diversification/data/Paleobiology_Database.RData",sep=""));

pbdb_data_list_smol <- list(pbdb_data_list$pbdb_sites_refined,pbdb_data_list$pbdb_finds,pbdb_data_list$pbdb_taxonomy,pbdb_data_list$time_scale[pbdb_data_list$time_scale$scale=="Stage Slice",]);
names(pbdb_data_list_smol) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
#pbdb_data_list_for_class <- list(pbdb_sites_refined,pbdb_finds,pbdb_taxonomy,finest_chronostrat);
#names(pbdb_data_list_for_class) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
#save(pbdb_data_list_for_class,file=paste("~/Documents/R_Projects/Data_for_R/PBDB_Data_for_Invert_Paleo.RData",sep=""));
# 2023-15-03: Update paleogeography ####
#print("Remapping collections with new dates…")
#pbdb_sites_refined_orig <- pbdb_sites_refined_orig[order(pbdb_sites_refined_orig$collection_no),];
#pbdb_sites_refined <- pbdb_sites_refined[order(pbdb_sites_refined$collection_no),];
#new_ma_mids <- round((pbdb_sites_refined$ma_lb+pbdb_sites_refined$ma_ub)/2,0);
#old_ma_mids <- round((pbdb_sites_refined_orig$ma_lb+pbdb_sites_refined_orig$ma_ub)/2,0);
#length(new_ma_mids); length(old_ma_mids)
#collections_to_remap <- (1:nrow(pbdb_sites_refined))[new_ma_mids!=old_ma_mids];
#sum(is.na(pbdb_sites_refined$lat[collections_to_remap]));
#if (length(collections_to_remap)>0)	{
#	roh_rows <- collections_to_remap[abs(pbdb_sites_refined$lat[collections_to_remap])>90];
#	if (length(roh_rows)>0)	{
#		redone <- accersi_collection_data_for_list_of_collection_nos(pbdb_sites_refined$collection_no[roh_rows]);
#		pbdb_sites_refined$lat[roh_rows] <- redone$lat;
#		pbdb_sites_refined$lng[roh_rows] <- redone$lng;
#		}
#	roh_rows <- collections_to_remap[abs(pbdb_sites_refined$lat[collections_to_remap])>90];
#	collections_to_remap <- collections_to_remap[!collections_to_remap %in% roh_rows];
#	redone_sites <- mundify_paleolatitude_and_paleolongitude(old_sites=pbdb_sites_refined[collections_to_remap,]);
#	pbdb_sites_refined$paleolat[match(redone_sites$collection_no,pbdb_sites_refined$collection_no)] <- redone_sites$paleolat;
#	pbdb_sites_refined$paleolng[match(redone_sites$collection_no,pbdb_sites_refined$collection_no)] <- redone_sites$paleolng;
#	pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
#	pbdb_data_list_smol$pbdb_sites <- pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
#	save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#	save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
#	}
# The Gray Havens ####
epilogue <- paste("Completing run started at",start_time,"at",date());
if (length(do_not_forget_taxa)>1)	{
	epilogue <- paste(epilogue,"with complete download of");
	if (length(do_not_forget_taxa)==1)	{
		epilogue <- paste(epilogue,do_not_forget_taxa);
		} else if (length(do_not_forget_taxa)==2)	{
		do_not_forget_taxa[2] <- paste(" &",do_not_forget_taxa[2]);
		epilogue <- paste(epilogue,paste(do_not_forget_taxa,collapse=""));
		} else if (length(do_not_forget_taxa)>2)	{
		do_not_forget_taxa[length(do_not_forget_taxa)] <- paste(" &",do_not_forget_taxa[length(do_not_forget_taxa)]);
		do_not_forget_taxa[2:(length(do_not_forget_taxa)-1)] <- paste(",",do_not_forget_taxa[2:(length(do_not_forget_taxa)-1)]);
		epilogue <- paste(epilogue,paste(do_not_forget_taxa,collapse=""));
		}
	if (do_not_forget_start!="")	{
		epilogue <- paste(epilogue,paste("from the"));
		if (do_not_forget_start!=do_not_forget_end)	{
			epilogue <- paste(epilogue,paste(do_not_forget_start,"-",do_not_forget_end));
			} else	{
			epilogue <- paste(epilogue,do_not_forget_start);
			}
		}
	} else if (do_not_forget_start!="")	{
	epilogue <- paste(epilogue,paste("with complete download of data from the"));
	if (do_not_forget_start!=do_not_forget_end)	{
		epilogue <- paste(epilogue,paste(do_not_forget_start,"-",do_not_forget_end));
		} else	{
		epilogue <- paste(epilogue,do_not_forget_start);
		}
	}
print(epilogue);
beepr::beep("wilhelm");
return(pbdb_data_list);
}

fix_sites_based_on_rock_data_only <- function(pbdb_data_list,rock_database)	{
pbdb_sites_refined <- pbdb_data_list$pbdb_sites_refined;
nsites <- nrow(pbdb_sites_refined);
sites_w_rocks_in_database <- (1:nsites)[pbdb_sites_refined$rock_no_sr>0];
too_old <- sites_w_rocks_in_database[pbdb_sites_refined$ma_lb[sites_w_rocks_in_database]-rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[sites_w_rocks_in_database],rock_database$rock_no)]>0];
too_young <- sites_w_rocks_in_database[pbdb_sites_refined$ma_ub[sites_w_rocks_in_database]-rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[sites_w_rocks_in_database],rock_database$rock_no)]<0];
pbdb_sites_refined$ma_lb[too_old] <- rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[too_old],rock_database$rock_no)]
pbdb_sites_refined$ma_ub[too_young] <- rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[too_young],rock_database$rock_no)]
effed <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
pbdb_sites_refined$ma_lb[effed] <- pbdb_data_list$pbdb_sites_refined$ma_lb[effed];
pbdb_sites_refined$ma_ub[effed] <- pbdb_data_list$pbdb_sites_refined$ma_ub[effed];
wonky_sites <- pbdb_sites_refined[effed,];
wonky_sites <- wonky_sites[order(wonky_sites$rock_no_sr),];
rock_lb <- rock_database$ma_lb[match(wonky_sites$rock_no_sr,rock_database$rock_no)];
rock_ub <- rock_database$ma_ub[match(wonky_sites$rock_no_sr,rock_database$rock_no)];
wonky_sites <- tibble::add_column(wonky_sites, rock_ub=as.numeric(rock_ub), .after = match("ma_ub",colnames(wonky_sites)));
wonky_sites <- tibble::add_column(wonky_sites, rock_lb=as.numeric(rock_lb), .after = match("ma_ub",colnames(wonky_sites)));
write.csv(wonky_sites,"Wonky_Sites.csv",row.names=F,fileEncoding = 'UTF-8');

pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
pbdb_data_list_smol <- list(pbdb_data_list$pbdb_sites_refined,pbdb_data_list$pbdb_finds,pbdb_data_list$pbdb_taxonomy,pbdb_data_list$time_scale[pbdb_data_list$time_scale$scale=="Stage Slice",]);
names(pbdb_data_list_smol) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
return(pbdb_data_list);
}

#rock_to_zone_database <- rock_unit_data$rock_to_zone_database;
#rock_database <- rock_unit_data$rock_unit_database;
fix_sites_based_on_rock_data_and_rock_to_zone_data <- function(pbdb_data_list,rock_database,rock_to_zone_database)	{
pbdb_sites_refined <- pbdb_data_list$pbdb_sites_refined;
nsites <- nrow(pbdb_sites_refined);
sites_w_rocks_in_database <- (1:nsites)[pbdb_sites_refined$rock_no_sr>0];
too_old <- sites_w_rocks_in_database[pbdb_sites_refined$ma_lb[sites_w_rocks_in_database]-rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[sites_w_rocks_in_database],rock_database$rock_no)]>0];
too_young <- sites_w_rocks_in_database[pbdb_sites_refined$ma_ub[sites_w_rocks_in_database]-rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[sites_w_rocks_in_database],rock_database$rock_no)]<0];
pbdb_sites_refined$ma_lb[too_old] <- rock_database$ma_lb[match(pbdb_sites_refined$rock_no_sr[too_old],rock_database$rock_no)]
pbdb_sites_refined$ma_ub[too_young] <- rock_database$ma_ub[match(pbdb_sites_refined$rock_no_sr[too_young],rock_database$rock_no)]
effed <- (1:nsites)[pbdb_sites_refined$ma_lb<=pbdb_sites_refined$ma_ub];
pbdb_sites_refined$ma_lb[effed] <- pbdb_data_list$pbdb_sites_refined$ma_lb[effed];
pbdb_sites_refined$ma_ub[effed] <- pbdb_data_list$pbdb_sites_refined$ma_ub[effed];
wonky_sites <- pbdb_sites_refined[effed,];
wonky_sites <- wonky_sites[order(wonky_sites$rock_no_sr),];
rock_lb <- rock_database$ma_lb[match(wonky_sites$rock_no_sr,rock_database$rock_no)];
rock_ub <- rock_database$ma_ub[match(wonky_sites$rock_no_sr,rock_database$rock_no)];
wonky_sites <- tibble::add_column(wonky_sites, rock_ub=as.numeric(rock_ub), .after = match("ma_ub",colnames(wonky_sites)));
wonky_sites <- tibble::add_column(wonky_sites, rock_lb=as.numeric(rock_lb), .after = match("ma_ub",colnames(wonky_sites)));
wsites <- nrow(wonky_sites);
write.csv(wonky_sites,"Wonky_Sites.csv",row.names=F,fileEncoding = 'UTF-8');

for (i in 1:wsites)	{
	if (wonky_sites$ma_ub[i]==wonky_sites$rock_lb[i])	{

		} else if (wonky_sites$ma_lb[i]==wonky_sites$rock_ub[i])	{

		}
	}

pbdb_data_list$pbdb_sites_refined <- pbdb_sites_refined;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
pbdb_data_list_smol <- list(pbdb_data_list$pbdb_sites_refined,pbdb_data_list$pbdb_finds,pbdb_data_list$pbdb_taxonomy,pbdb_data_list$time_scale[pbdb_data_list$time_scale$scale=="Stage Slice",]);
names(pbdb_data_list_smol) <- c("pbdb_sites","pbdb_finds","pbdb_taxonomy","time_scale");
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
}

#taxon <- "Cephalaspidomorpha";
update_taxonomic_data_in_pbdb_RData <- function(taxon,pbdb_data_list,pbdb_taxonomy_edits)	{
print(paste("Getting taxonomic information for",taxon));
reloaded_taxonomy <- update_taxonomy(taxon);
reloaded_taxonomy <- reloaded_taxonomy[order(reloaded_taxonomy$taxon_no,reloaded_taxonomy$orig_no),];
pbdb_taxonomy_edits <- pbdb_taxonomy_edits[order(pbdb_taxonomy_edits$taxon_no,pbdb_taxonomy_edits$orig_no),];

reloaded_taxonomy[reloaded_taxonomy$taxon_no %in% pbdb_taxonomy_edits$taxon_no,] <- pbdb_taxonomy_edits[pbdb_taxonomy_edits$taxon_no %in% reloaded_taxonomy$taxon_no,];
reloaded_taxonomy <- reloaded_taxonomy[order(reloaded_taxonomy$taxon_no,reloaded_taxonomy$orig_no),];
reloaded_taxonomy <- clean_the_bastards(pbdb_data=reloaded_taxonomy);

pbdb_data_list$pbdb_taxonomy <- rbind(pbdb_data_list$pbdb_taxonomy,pbdb_taxonomy_edits[!pbdb_taxonomy_edits$taxon_no %in% pbdb_data_list$pbdb_taxonomy$taxon_no,]);

pbdb_taxonomy <- clean_the_bastards(pbdb_data_list$pbdb_taxonomy);
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];
updated_taxonomy <- reloaded_taxonomy[reloaded_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no,];
new_taxonomy <- reloaded_taxonomy[!reloaded_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no,];

pbdb_finds <- pbdb_data_list$pbdb_finds;
print(paste("Getting taxonomic information for taxa formerly classified in",taxon));
this_taxon_info <- reloaded_taxonomy[reloaded_taxonomy$taxon_name %in% taxon,];
if (nrow(this_taxon_info)==0)
	this_taxon_info <- update_taxonomy(taxon,inc_children = F);
if (nrow(this_taxon_info)>0)
	this_taxon_info <- this_taxon_info[this_taxon_info$flags %in% c("","B"),];
if (nrow(this_taxon_info)==1)	{
	taxon_rank <- this_taxon_info$taxon_rank;
	} else	{
	if (nrow(this_taxon_info)==0)	{
		this_taxon_info <- updated_taxonomy[updated_taxonomy$taxon_name %in% taxon,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$flags %in% c("","B","BV"),];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$n_occs>0,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$orig_no==this_taxon_info$taxon_no,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$modified==max(this_taxon_info$modified),];
		if (nrow(this_taxon_info)==1)
			taxon_rank <- this_taxon_info$taxon_rank;
		} else	{
		attributes <- gsub("\\(","",this_taxon_info$taxon_attr);
		attributes <- gsub("\\)","",attributes);
		if (length(unique(attributes))==1)	{

			}
		print("Pete, we got problems!");
		return(pbdb_data_list);
		}
	}
if (taxon_rank %in% colnames(pbdb_finds))	{
	candidates <- data.frame(taxon_no=this_taxon_info$taxon_no,taxon=taxon,taxon_rank=taxon_rank);
	} else	{
	candidates <- data.frame(taxon_no=as.numeric(),taxon=as.character(),taxon_rank=as.character());
	prob_children <- taxon;
	i <- 0;
	while (i < length(prob_children))	{
		i <- i+1;
		taxon_children_info <- updated_taxonomy[updated_taxonomy$parent_name %in% prob_children[i],];
		taxon_children_info <- taxon_children_info[taxon_children_info$flags %in% c("","B"),]
		poss_candidates <- data.frame(taxon_no=as.numeric(taxon_children_info$taxon_no),
									  taxon=as.character(taxon_children_info$taxon_name),
									  taxon_rank=taxon_children_info$taxon_rank);
		candidates <- rbind(candidates,poss_candidates[poss_candidates$taxon_rank %in% colnames(pbdb_finds),]);
		prob_children <- c(prob_children,poss_candidates$taxon[!poss_candidates$taxon_rank %in% colnames(pbdb_finds)]);
		}
	}

# look for taxa that were classified in this taxon but that have been removed
candidates$rank_no <- match(candidates$taxon_rank,taxonomic_rank);
candidates <- candidates[order(-candidates$rank_no),];
missing_taxonomy <- pbdb_taxonomy[1,];
cc <- 0;
while (cc < nrow(candidates))	{
	cc <- cc+1;
	taxon_col <- match(candidates$taxon_rank[cc],colnames(pbdb_taxonomy));
	taxon_info <- pbdb_taxonomy[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc],]
	taxon_no_list <- pbdb_taxonomy$taxon_no[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc]][!pbdb_taxonomy$taxon_no[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc]] %in% reloaded_taxonomy$taxon_no];
	if (length(taxon_no_list)>0) {
		if (candidates$rank_no[cc]>4) {
			print(paste("Editting taxa once or currently placed in the",candidates$taxon[cc]));
			} else	{
			print(paste("Editting",candidates$taxon[cc]));
			}
		missing_taxon_info <- accersi_taxonomic_data_for_list_of_taxon_nos(taxon_no_list,);
		missing_taxon_info[,!colnames(missing_taxon_info) %in% colnames(pbdb_taxonomy)] <- NULL;
		missing_taxonomy <- rbind(missing_taxonomy,missing_taxon_info)
		}
#	missing_taxa_nos <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$class=="Trilobita"][!pbdb_taxonomy$taxon_no[pbdb_taxonomy$class=="Trilobita"] %in% reloaded_taxonomy$taxon_no];
	}
if (nrow(missing_taxonomy)>1)	{
	missing_taxonomy <- missing_taxonomy[2:nrow(missing_taxonomy),];
	updated_taxonomy <- rbind(updated_taxonomy,missing_taxonomy);
	}

updated_taxonomy <- updated_taxonomy[order(updated_taxonomy$taxon_no,updated_taxonomy$orig_no),];
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];

#pbdb_taxonomy$taxon_no[hist(pbdb_taxonomy$taxon_no,breaks=c(0,pbdb_taxonomy$taxon_no),plot=F)$counts==2]
#pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no[hist(pbdb_taxonomy$taxon_no,breaks=c(0,pbdb_taxonomy$taxon_no),plot=F)$counts==2],]
#delete <- (1:nrow(pbdb_taxonomy))[pbdb_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no[hist(pbdb_taxonomy$taxon_no,breaks=c(0,pbdb_taxonomy$taxon_no),plot=F)$counts==2]][c(1,3,5)];
#delete <- delete[!is.na(delete)]
#pbdb_taxonomy <- pbdb_taxonomy[!(1:nrow(pbdb_taxonomy)) %in% delete,];
#updated_taxonomy <- updated_taxonomy[match(unique(updated_taxonomy$taxon_no),updated_taxonomy$taxon_no),];
#length(updated_taxonomy$taxon_no)
#length(unique(updated_taxonomy$taxon_no));
#pbdb_taxonomy <- rbind(pbdb_taxonomy,updated_taxonomy[!updated_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no,])
#pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no),];
taxon_no_counts <- hist(pbdb_taxonomy$taxon_no,breaks=0:max(pbdb_taxonomy$taxon_no),plot=FALSE)$counts;
names(taxon_no_counts) <- 1:max(pbdb_taxonomy$taxon_no);
taxon_no_counts <- taxon_no_counts[taxon_no_counts>1];
tnc <- 0;
while (tnc < length(taxon_no_counts))	{
	tnc <- tnc+1;
	txn <- as.numeric(names(taxon_no_counts)[tnc]);
	keeper <- match(max(pbdb_taxonomy$modified[pbdb_taxonomy$taxon_no %in% txn]),pbdb_taxonomy$modified[pbdb_taxonomy$taxon_no %in% txn]);
	txn_rows <- (1:nrow(pbdb_taxonomy))[pbdb_taxonomy$taxon_no %in% txn];
	txn_rows <- txn_rows[!(1:length(txn_rows)) %in% keeper];
	pbdb_taxonomy <- pbdb_taxonomy[!(1:nrow(pbdb_taxonomy)) %in% txn_rows,];
	}

pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% updated_taxonomy$taxon_no,] <- updated_taxonomy;
pbdb_taxonomy <- rbind(pbdb_taxonomy,new_taxonomy);
pbdb_taxonomy <- unique(pbdb_taxonomy);
rep_taxon_nos <- hist(match(pbdb_taxonomy$taxon_no,unique(pbdb_taxonomy$taxon_no)),breaks=0:max(pbdb_taxonomy$taxon_no),plot=F)$counts;
names(rep_taxon_nos) <- 1:length(rep_taxon_nos);
rep_taxon_nos <- rep_taxon_nos[rep_taxon_nos>1];
rt <- 0;
while (rt < length(rep_taxon_nos))	{
	rt <- rt+1;
	doppel <- as.numeric(names(rep_taxon_nos[rt]));
	# fill in the rest later!
	}
updated_genera <- reloaded_taxonomy$taxon_name[reloaded_taxonomy$taxon_rank %in% c("genus","subgenus")];
updated_genera_nos <- reloaded_taxonomy$taxon_no[reloaded_taxonomy$taxon_rank %in% c("genus","subgenus")];

print("Updating occurrence data");
# changed 15-12-2022 to change all finds
# updated 22-04-2023 to update subgenera first
sgfinds <- (1:nrow(pbdb_finds))[pbdb_finds$subgenus_no %in% pbdb_taxonomy$taxon_no];
pbdb_finds$genus[sgfinds] <- pbdb_taxonomy$genus[match(pbdb_finds$subgenus_no[sgfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$genus_no[sgfinds] <- pbdb_taxonomy$genus_no[match(pbdb_finds$subgenus_no[sgfinds],pbdb_taxonomy$taxon_no)];

gfinds <- (1:nrow(pbdb_finds))[pbdb_finds$genus_no %in% pbdb_taxonomy$taxon_no];
pbdb_finds$family[gfinds] <- pbdb_taxonomy$family[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$family_no[gfinds] <- pbdb_taxonomy$family_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order_no[gfinds] <- pbdb_taxonomy$order_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order[gfinds] <- pbdb_taxonomy$order[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$class_no[gfinds] <- pbdb_taxonomy$class_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$class[gfinds] <- pbdb_taxonomy$class[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$phylum_no[gfinds] <- pbdb_taxonomy$phylum_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$phylum[gfinds] <- pbdb_taxonomy$phylum[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$genus[gfinds] <- pbdb_taxonomy$genus[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$genus_no[gfinds] <- pbdb_taxonomy$genus_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];

pbdb_taxonomy <- clean_the_bastards(pbdb_taxonomy);
pbdb_finds <- clean_the_bastards(pbdb_finds);

pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
pbdb_data_list$pbdb_finds <- pbdb_finds;

print("Saving updated PBDB Data.");
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
return(pbdb_data_list);
}

update_taxonomic_data_in_pbdb_RData_one_taxon_only <- function(taxon,pbdb_data_list,pbdb_taxonomy_edits)	{
print(paste("Getting taxonomic information for",taxon));
reloaded_taxonomy <- update_taxonomy(taxon,inc_children=F);
reloaded_taxonomy <- reloaded_taxonomy[order(reloaded_taxonomy$taxon_no,reloaded_taxonomy$orig_no),];
pbdb_taxonomy_edits <- pbdb_taxonomy_edits[order(pbdb_taxonomy_edits$taxon_no,pbdb_taxonomy_edits$orig_no),];

reloaded_taxonomy[reloaded_taxonomy$taxon_no %in% pbdb_taxonomy_edits$taxon_no,] <- pbdb_taxonomy_edits[pbdb_taxonomy_edits$taxon_no %in% reloaded_taxonomy$taxon_no,];
reloaded_taxonomy <- reloaded_taxonomy[order(reloaded_taxonomy$taxon_no,reloaded_taxonomy$orig_no),];
reloaded_taxonomy <- clean_the_bastards(pbdb_data=reloaded_taxonomy);

pbdb_data_list$pbdb_taxonomy <- rbind(pbdb_data_list$pbdb_taxonomy,pbdb_taxonomy_edits[!pbdb_taxonomy_edits$taxon_no %in% pbdb_data_list$pbdb_taxonomy$taxon_no,]);
#pbdb_data_list$pbdb_taxonomy[pbdb_data_list$pbdb_taxonomy$taxon_name=="Erniettamorpha",] <- pbdb_taxonomy_edits[pbdb_taxonomy_edits$taxon_name=="Erniettamorpha",];
pbdb_taxonomy <- clean_the_bastards(pbdb_data_list$pbdb_taxonomy);
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];
updated_taxonomy <- reloaded_taxonomy[reloaded_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no,];
new_taxonomy <- reloaded_taxonomy[!reloaded_taxonomy$taxon_no %in% pbdb_taxonomy$taxon_no,];

pbdb_finds <- pbdb_data_list$pbdb_finds;
print(paste("Getting taxonomic information for taxa formerly classified in",taxon));
this_taxon_info <- reloaded_taxonomy[reloaded_taxonomy$taxon_name %in% taxon,];
if (nrow(this_taxon_info)>0)	this_taxon_info <- this_taxon_info[this_taxon_info$flags %in% c("","B"),];
if (nrow(this_taxon_info)==1)	{
	taxon_rank <- this_taxon_info$taxon_rank;
	} else	{
	if (nrow(this_taxon_info)==0)	{
		this_taxon_info <- updated_taxonomy[updated_taxonomy$taxon_name %in% taxon,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$flags %in% c("","B","BV"),];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$n_occs>0,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$orig_no==this_taxon_info$taxon_no,];
		if (nrow(this_taxon_info)>1)	this_taxon_info <- this_taxon_info[this_taxon_info$modified==max(this_taxon_info$modified),];
		} else	{
		return(print("Pete, we got problems!"));
		}
	}
if (taxon_rank %in% colnames(pbdb_finds))	{
	candidates <- data.frame(taxon_no=this_taxon_info$taxon_no,taxon=taxon,taxon_rank=taxon_rank);
	} else	{
	candidates <- data.frame(taxon_no=as.numeric(),taxon=as.character(),taxon_rank=as.character());
	prob_children <- taxon;
	i <- 0;
	while (i < length(prob_children))	{
		i <- i+1;
		taxon_children_info <- updated_taxonomy[updated_taxonomy$parent_name %in% prob_children[i],];
		taxon_children_info <- taxon_children_info[taxon_children_info$flags %in% c("","B"),]
		poss_candidates <- data.frame(taxon_no=as.numeric(taxon_children_info$taxon_no),
									  taxon=as.character(taxon_children_info$taxon_name),
									  taxon_rank=taxon_children_info$taxon_rank);
		candidates <- rbind(candidates,poss_candidates[poss_candidates$taxon_rank %in% colnames(pbdb_finds),]);
		prob_children <- c(prob_children,poss_candidates$taxon[!poss_candidates$taxon_rank %in% colnames(pbdb_finds)]);
		}
	}

# look for taxa that were classified in this taxon but that have been removed
candidates$rank_no <- match(candidates$taxon_rank,taxonomic_rank);
candidates <- candidates[order(-candidates$rank_no),];
missing_taxonomy <- pbdb_taxonomy[1,];
cc <- 0;
while (cc < nrow(candidates))	{
	cc <- cc+1;
	taxon_col <- match(candidates$taxon_rank[cc],colnames(pbdb_taxonomy));
	taxon_info <- pbdb_taxonomy[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc],]
	taxon_no_list <- pbdb_taxonomy$taxon_no[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc]][!pbdb_taxonomy$taxon_no[pbdb_taxonomy[,taxon_col] %in% candidates$taxon[cc]] %in% reloaded_taxonomy$taxon_no];
	if (length(taxon_no_list)>0) {
		if (candidates$rank_no[cc]>4) {
			print(paste("Editting taxa once or currently placed in the",candidates$taxon[cc]));
			} else	{
			print(paste("Editting",candidates$taxon[cc]));
			}
		missing_taxon_info <- accersi_taxonomic_data_for_list_of_taxon_nos(taxon_no_list,);
		missing_taxon_info[,!colnames(missing_taxon_info) %in% colnames(pbdb_taxonomy)] <- NULL;
		missing_taxonomy <- rbind(missing_taxonomy,missing_taxon_info)
		}
#	missing_taxa_nos <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$class=="Trilobita"][!pbdb_taxonomy$taxon_no[pbdb_taxonomy$class=="Trilobita"] %in% reloaded_taxonomy$taxon_no];
	}
if (nrow(missing_taxonomy)>1)	{
	missing_taxonomy <- missing_taxonomy[2:nrow(missing_taxonomy),];
	updated_taxonomy <- rbind(updated_taxonomy,missing_taxonomy);
	}

updated_taxonomy <- updated_taxonomy[order(updated_taxonomy$taxon_no,updated_taxonomy$orig_no),];
pbdb_taxonomy <- pbdb_taxonomy[order(pbdb_taxonomy$taxon_no,pbdb_taxonomy$orig_no),];

pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% updated_taxonomy$taxon_no,] <- updated_taxonomy;
pbdb_taxonomy <- rbind(pbdb_taxonomy,new_taxonomy);
pbdb_taxonomy <- unique(pbdb_taxonomy)
updated_genera <- reloaded_taxonomy$taxon_name[reloaded_taxonomy$taxon_rank %in% c("genus","subgenus")];
updated_genera_nos <- reloaded_taxonomy$taxon_no[reloaded_taxonomy$taxon_rank %in% c("genus","subgenus")];

print("Updating occurrence data");
# changed 15-12-2022 to change all finds
gfinds <- (1:nrow(pbdb_finds))[pbdb_finds$genus_no %in% pbdb_taxonomy$taxon_no];
pbdb_finds$family[gfinds] <- pbdb_taxonomy$family[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$family_no[gfinds] <- pbdb_taxonomy$family_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order_no[gfinds] <- pbdb_taxonomy$order_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$order[gfinds] <- pbdb_taxonomy$order[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$class_no[gfinds] <- pbdb_taxonomy$class_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$class[gfinds] <- pbdb_taxonomy$class[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$phylum_no[gfinds] <- pbdb_taxonomy$phylum_no[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];
pbdb_finds$phylum[gfinds] <- pbdb_taxonomy$phylum[match(pbdb_finds$genus_no[gfinds],pbdb_taxonomy$taxon_no)];

pbdb_taxonomy <- clean_the_bastards(pbdb_taxonomy);
pbdb_finds <- clean_the_bastards(pbdb_finds);

pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
pbdb_data_list$pbdb_finds <- pbdb_finds;

save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
return(pbdb_data_list);
}

fix_subgenera_parents_in_taxonomy_database_only <- function(pbdb_taxonomy)	{
subgenera <- pbdb_taxonomy[pbdb_taxonomy$accepted_rank == "subgenus",];
subgenera[1,]
parent_genus <- base::t(sapply(subgenera$accepted_name,divido_subgenus_names_from_genus_names))[,1];
parent_genus_no <- pbdb_taxonomy$taxon_no[match(parent_genus,pbdb_taxonomy$taxon_name)]

}

update_pbdb_taxonomy_RData <- function(pbdb_data_list,paleodb_fixes)	{
#start_time <- date();
pbdb_taxonomy <- put_pbdb_dataframes_into_proper_type(pbdb_data_list$pbdb_taxonomy);
# get taxonomic updates ####
print("Getting taxonomic opinions & data entered/modified since the last download");
#taxa_modified_after <- as.Date(max(pbdb_taxonomy$modified[!is.na(pbdb_taxonomy$modified)]))-7;
taxa_modified_after <- as.Date(max(pbdb_taxonomy$modified[nrow(pbdb_taxonomy)]))-2;
#taxa_modified_after <- "2010-04-30";
new_taxa <- update_taxonomy(taxa_modified_after = taxa_modified_after,inc_children=T);
new_taxa <- put_pbdb_dataframes_into_proper_type(new_taxa);
new_taxa <- clean_pbdb_fields(new_taxa);

new_opinions <- update_opinions(ops_modified_after = taxa_modified_after);
#new_opinions <- rbind(new_opinions,update_opinions("Ichnofossils",ops_modified_after = taxa_modified_after));
new_opinions <- put_pbdb_dataframes_into_proper_type(new_opinions);
new_opinions <- clean_pbdb_fields(new_opinions);

# get older opinions for taxa with new opinions
opinionated <- gsub(" ","%20",unique(new_opinions$taxon_name));
print("Collecting older opinions for taxa with new opinions");
old_opinions <- pbdb_opinions;
old_opinions <- old_opinions[old_opinions$opinion_no==0,];
for (i in 1:length(opinionated))
	old_opinions <- rbind(old_opinions,update_opinions(taxon=opinionated[i],inc_children = F));
#old_opinions <- base::t(pbapply::pbsapply(opinionated,update_opinions,inc_children = F));
# eliminate duplicate opinions;
opinion_nos <- pbdb_opinions$opinion_no;
opinion_no_counts <- hist(opinion_nos,breaks=0:max(opinion_nos),plot=F)$counts;
names(opinion_no_counts) <- 1:max(opinion_nos);
dup_opinions <- opinion_no_counts[opinion_no_counts>1];
dp <- 0;
while (dp < length(dup_opinions))	{
	dp <- dp+1;
	fix_me <- as.numeric(names(dup_opinions)[dp]);
	prob_rows <- (1:nrow(pbdb_opinions))[pbdb_opinions$opinion_no %in% fix_me];
	last_modified <- max(pbdb_opinions$modified[pbdb_opinions$opinion_no %in% fix_me]);
	keep_me <- prob_rows[pbdb_opinions$modified[pbdb_opinions$opinion_no %in% fix_me] %in% last_modified]
	kill_me <- prob_rows[!prob_rows %in% keep_me[1]];
	pbdb_opinions <- pbdb_opinions[(1:nrow(pbdb_opinions)) %in% kill_me,];
	}

for (i in 1:nrow(old_opinions))	{
	fix_me <- match(old_opinions$opinion_no[i],pbdb_opinions$opinion_no);
	if (!is.na(fix_me))
		pbdb_opinions[fix_me,] <- old_opinions[i,];
	}
#pbdb_opinions[pbdb_opinions$opinion_no %in% old_opinions$opinion_no,] <- nrow(old_opinions[old_opinions$opinion_no %in% pbdb_opinions$opinion_no,]);

print("Collecting updated taxonomic information.");
taxa_modified_after <- taxa_modified_after-5;
#taxa_modified_after <- "2010-04-30";
new_taxa <- update_taxonomy(taxa_modified_after = taxa_modified_after,inc_children = T);
#new_taxa <- rbind(new_taxa,update_taxonomy(taxon = "Ichnofossils",taxa_modified_after = taxa_modified_after,inc_children = T));
new_taxa <- put_pbdb_dataframes_into_proper_type(new_taxa);
new_taxa <- clean_pbdb_fields(new_taxa);
#new_taxa[new_taxa$taxon_name=="Agnostida",]
new_opinions_spc <- new_opinions[new_opinions$taxon_rank %in% c("species","subspecies"),];
new_opinions_spc <- new_opinions_spc[new_opinions_spc$opinion_type=="class",];
new_taxa <- new_taxa[order(new_taxa$taxon_no),];
redone_taxa <- new_taxa[new_taxa$taxon_no %in% pbdb_taxonomy$taxon_no,];
pbdb_taxonomy[pbdb_taxonomy$taxon_no %in% new_taxa$taxon_no,] <- redone_taxa;
totally_new_taxa <- new_taxa[!new_taxa$taxon_no %in% pbdb_taxonomy$taxon_no,];
pbdb_taxonomy <- rbind(pbdb_taxonomy,totally_new_taxa);

# put a bow on it and gift it to yourself ####
pbdb_data_list$pbdb_taxonomy <- pbdb_taxonomy;
#pbdb_data_list$pbdb_opinions <- pbdb_opinions;
save(pbdb_data_list,file=paste("~/Documents/R_Projects/Data_for_R/Paleobiology_Database.RData",sep=""));
#save(pbdb_data_list,file=paste("~/Documents/R_Projects/Marine_Tetrapod_Diversification/data/Paleobiology_Database.RData",sep=""));

pbdb_data_list_smol$pbdb_taxonomy <- pbdb_taxonomy;
save(pbdb_data_list_smol,file=paste(data_for_r,"/PBDB_Data_Smøl.RData",sep=""));
return(pbdb_data_list);
}


fixing_crap <- function()	{
test <- hist(pbdb_taxonomy$taxon_no,breaks=0:max(pbdb_taxonomy$taxon_no),plot=F)$counts
names(test) <- 1:max(pbdb_taxonomy$taxon_no);
test[test==2]
pbdb_taxonomy[pbdb_taxonomy$taxon_no==326445 ,]

test <- hist(pbdb_taxonomy_edits$taxon_no,breaks=0:max(pbdb_taxonomy_edits$taxon_no),plot=F)$counts
names(test) <- 1:max(pbdb_taxonomy_edits$taxon_no);
test[test==2]

pbdb_sites[pbdb_sites$collection_no %in% 218987,];
}
#pbdb_taxonomy$taxon_environment[pbdb_taxonomy$taxon_name=="Gigantopteridales"] <- "terrestrial"

kluge_update <- function(tn,pbdb_taxonomy,updated_taxonomy)	{
nt <- match(tn,pbdb_taxonomy$taxon_no)
#rbind(updated_taxonomy[i,],pbdb_taxonomy[nt,])
pbdb_taxonomy[nt,] <- updated_taxonomy[tn,];
return(pbdb_taxonomy);
}

{}
