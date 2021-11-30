#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal
# scribere: write
#install.packages("devtools", dependencies=TRUE);
library(paleobioDB);    #install.packages("paleobioDB", dependencies=TRUE);
library(stringr);       #install.packages("stringr", dependencies=TRUE)
#library(lettercase);	#devtools::install_github('decisionpatterns/lettercase')
library(Rcpp);			#devtools::install_github("https://github.com/RcppCore/Rcpp/")
library(rvest);
library(dplyr);
library(prodlim);		#install.packages("prodlim", dependencies=TRUE);
library(xml2);			#install.packages("xml2", dependencies=TRUE);

classic_articulate_brachiopod_taxa <- c("Chileata","Kutorginata","Obolellata","Rhynchonellata","Strophomenata");
classic_inarticulate_brachiopod_taxa <- c("Craniata","Lingulata","Paterinata","Tegulella");
gastropod_control_taxa <- c("Gastropoda","Helcionelloida","Paragastropoda","Polyplacophora","Tergomya");
bivalve_control_taxa <- c("Bivalvia","Rostroconchia","Scaphopoda");
uncertains <- c("?","cf.","aff.","ex_gr.");
taxon_qualifiers <- c("\\? ","cf. ","aff. ","ex_gr. ","n. sp. ","n. gen. "," informal");
nomens <- c("nomen dubium","nomen nudum","nomen oblitum","nomen vanum");
hell_no <- F;

carbonate_environments <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
siliciclastic_environments <- c("marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
marine_environments <- c("marine indet.",carbonate_environments,siliciclastic_environments);

sedimentary_rocks <- c("arenarie","areniscas","argillaceous","argilliti","ashes","ash","calcaerous","calcaire","carbonate","chalk","cherty","cherts","chert","clay","claystone","claystones","conglomerates","conglomerate","coquina","coquinas","diatomite","dolomites","dolomite","dolostones","dolostone","flags","glauconites","glauconite","glauconitics","glauconitic","gres","grauwacke","greywacke","greywackes","grits","grit","kalk","kalkmergel","limestone","limestones","limeston","limstone","ls.","ls","lst","lst.","lutitas","marlstones","marlstone","marl","marls","marly","marne","micrites","micrite","mergel","mudstones","mudstone","ooid","ooids","oolitic","pebbles","pebble","phosphatics","phosphatic","phosphorite","phosphorites","platy","qzt.","quartzite","quartzites","quarziti","reef","sands","sand","sandstone","sandstones","shales","schichten","schistes","shale","shaly","siltstones","siltstone","slates","slate","tillite","tillites","tuff","tuffs","volcanic","volcanics");
missing_taxon_assignment <- c("NP","NO","NC","NF","NG","");
missing_data_assignment <- c("NP","NO","NC","NF","NG","","coordinates not computable using this model");
paleodb_numeric_fields <- c("no","ma","size","occs","colls","geoplate","plate","subset","lat","lng","paleolat","paleolng","bin","max_ma","min_ma","n_colls","n_occs","gplate_no","splate_no","pubyr");
same_differences <- c("misspelling of","obsolete variant of","recombined as");
#taxonomic_rank <- c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily","family","superfamily","infraorder","suborder");
taxonomic_rank <- c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily","family","superfamily","infraorder","suborder","order","superorder","infraclass","subclass","class","superclass","subphylum","phylum","superphylum","subkingdom","kingdom","superkingdom","unranked clade","informal");
taxonomic_field <- c(taxonomic_rank,c("subspecies_no","species_no","subgenus_no","genus_no","subtribe_no","tribe_no","subfamily_no","family_no","superfamily_no","infraorder_no","suborder_no","order_no","superorder_no","infraclass_no","subclass_no","class_no","superclass_no","subphylum_no","phylum_no","superphylum_no","subkingdom_no","kingdom_no","superkingdom_no","unranked clade_no","informal_no"));
plant_phyla <- c("Aneurophytophyta","Angiospermae","Antherocerotophyta","Anthophyta","Archaeopteridophyta","Botryopteridiophyta","Bryophyta","Chlorophyta","Cladoxylophyta","Coniferophyta","Cycadeoideophyta","Cycadophyta","Cyanobacteria","Equisetophyta","Filicophyta","Ginkgophyta","Gnetophyta","Gymnospermae","Gymnospermophyta","Isoetophyta","Langiophytophyta","Lycophyta","Lycopodiophyta","Lycopodophyta","Magnoliophyta","Marattiophyta","Marchantiophyta","Moresnetiophyta","Noeggerathiophyta","Ophioglossophyta","Peltaspermophyta","Pinophyta","Polypodiophyta","Progymnospermophyta","Psilophytophyta","Pteridophyta","Pteridospermophyta","Rhodophyta","Rhyniophyta","Spermatophyta","Sphenophyllophyta","Sphenophyta","Thallophyta","Tracheophyta","Trimerophytophyta","Zosterophyllophyta");
protist_phyla <- c("Actinopoda","Bacillariophyta","Cyanobacteria","Foraminifera","Heterokontophyta","Ochrophyta","Prasinophyta","Radiolaria","Sarcodina");
obsolete <- c("subjective synonym of","objective synonym of","replaced by");
publication_years <- 1700:2100;
not_real_rock_names <- c("lower","middle","upper",0:100,as.character(as.roman(1:99)));
opinion_basis_rank <- c("stated with evidence","stated without evidence","","implied","second hand");
dummy_finds <- data.frame(occurrence_no=as.numeric(),record_type=as.character(),reid_no=as.numeric(),
						  flags=as.character(),
						  collection_no=as.numeric(),
						  identified_name=as.character(),identified_rank=as.character(),identified_no=as.numeric(),
						  difference=as.character(),
						  accepted_name=as.character(),accepted_rank=as.character(),accepted_no=as.numeric(),
						  early_interval=as.character(),late_interval=as.character(),max_ma=as.numeric(),min_ma=as.numeric(),
						  ref_author=as.character(),ref_pubyr=as.numeric(),reference_no=as.numeric(),
						  phylum=as.character(),phylum_no=as.numeric(),class=as.character(),class_no=as.numeric(),order=as.character(),order_no=as.numeric(),family=as.character(),family_no=as.numeric(),genus=as.character(),genus_no=as.numeric(),subgenus_no=as.numeric(),
						  occurrence_comments=as.character(),
						  authorizer=as.character(),enterer=as.character(),modifier=as.character(),
						  abund_value=as.character(),abund_unit=as.character(),
						  created=as.character(),modified=as.character(),
						  stringsAsFactors = F);
dummy_finds2 <- data.frame(occurrence_no=as.numeric(),record_type=as.character(),reid_no=as.numeric(),
						  flags=as.character(),
						  collection_no=as.numeric(),
						  identified_name=as.character(),identified_rank=as.character(),identified_no=as.numeric(),
						  difference=as.character(),
						  accepted_name=as.character(),accepted_rank=as.character(),accepted_no=as.numeric(),
						  early_interval=as.character(),late_interval=as.character(),max_ma=as.numeric(),min_ma=as.numeric(),
						  ref_author=as.character(),ref_pubyr=as.numeric(),reference_no=as.numeric(),
						  phylum=as.character(),phylum_no=as.numeric(),class=as.character(),class_no=as.numeric(),order=as.character(),order_no=as.numeric(),family=as.character(),family_no=as.numeric(),genus=as.character(),genus_no=as.numeric(),subgenus_no=as.numeric(),
						  occurrence_comments=as.character(),
						  authorizer=as.character(),enterer=as.character(),modifier=as.character(),
						  abund_value=as.character(),abund_unit=as.character(),
						  created=as.character(),modified=as.character(),accepted_name_orig=as.character(),
						  stringsAsFactors = F);

					##### ROUTINES TO DOWNLOAD DATA USING API #######
list_to_dataframe_for_pbdb_data <- function(pbdb_list)	{
ttl_entries <- 0;
for (nn in 1:nrow(pbdb_list))	ttl_entries <- ttl_entries+length(pbdb_list$created[[nn]]);
if (ttl_entries==0)	for (nn in 1:nrow(pbdb_list))	ttl_entries <- ttl_entries+length(pbdb_list$record_type[[nn]]);
unlist <- data.frame(array("",dim=c(ttl_entries,as.numeric(ncol(pbdb_list)))));
colnames(unlist) <- colnames(pbdb_list);
ttl_entries <- 0;
for (nn in 1:nrow(pbdb_list))	{
	if (length(pbdb_list$record_type[[nn]])>0)	{
		aa <- ttl_entries+1;
		zz <- ttl_entries+length(pbdb_list$record_type[[nn]]);
		for (cn in 1:ncol(pbdb_list))	{
			if (colnames(pbdb_list)[cn] %in% c("created","modified"))	{
				unlist[aa:zz,cn] <- pbdb_list[,cn][[nn]];
#				unlist[aa:zz,cn][unlist[aa:zz,cn] %in% "0000-00-00 00:00:00"] <- date();
				options(warn = -1);
				if (!is.na(as.numeric(unlist[aa,cn])))	{
					unlist[aa:zz,cn] <- as.character(as.Date.numeric(as.numeric(unlist[aa:zz,cn]),origin="1970-01-01"));
					} else if (sum(unlist[aa:zz,cn] %in% "0000-00-00 00:00:00")==0)	{
					unlist[aa:zz,cn] <- as.Date.character(pbdb_list[,cn][[nn]])
#					unlist[aa:zz,cn] <- as.character(as.Date.numeric(as.numeric(unlist[aa:zz,cn]),origin="1970-01-01"));
					}
				options(warn = 1);
				} else	{
				unlist[aa:zz,cn] <- pbdb_list[,cn][[nn]];
				}
			}
		ttl_entries <- ttl_entries+length(pbdb_list$record_type[[nn]]);
		}
	}
unlist <- put_pbdb_dataframes_into_proper_type(pbdb_data=unlist);
return(unlist);
}

put_pbdb_dataframes_into_proper_type <- function(pbdb_data)	{
# routine to make sure that characters are characters & numbers are numbers and that factors burn in hell.
# written 2020-05-19
options(warn=-1);
for (cn in 1:ncol(pbdb_data))	{
	if (is.list(pbdb_data[,cn]))	pbdb_data[,cn] <- simplify2array(pbdb_data[,cn]);
	old_info <- pbdb_data[,cn];
	cnm <- strsplit(x=colnames(pbdb_data)[cn],split="_")[[1]];
	if((sum(cnm %in% paleodb_numeric_fields) > 0) && colnames(pbdb_data)[cn]!="collection_size")	{
		old_info[old_info %in% missing_data_assignment] <- 0;
		pbdb_data[,cn] <- as.numeric(old_info)
		pbdb_data[is.na(pbdb_data[,cn]),cn] <- 0;
		} else	{
		pbdb_data[,cn] <- as.character(old_info);
		pbdb_data[is.na(pbdb_data[,cn]),cn] <- "";
		}
	}
if (!is.null(pbdb_data$created))	{
	pbdb_date <- pbdb_data$created;
	ndates1 <- (1:length(pbdb_date))[!is.na(as.numeric(pbdb_date))];
	ndates2 <- (1:length(pbdb_date))[is.na(as.numeric(pbdb_date))];
	pbdb_date[ndates1] <- as.character(as.Date.numeric(as.numeric(pbdb_date[ndates1]),origin="1970-01-01"));
	pbdb_date[ndates2] <- as.Date.character(as.Date.character(pbdb_data$created[ndates2]));

	pbdb_date <- pbdb_data$modified;
	ndates1 <- (1:length(pbdb_date))[!is.na(as.numeric(pbdb_date))];
	ndates2 <- (1:length(pbdb_date))[is.na(as.numeric(pbdb_date))];
	ndates3 <- ndates2[pbdb_date[ndates2]!="0000-00-00 00:00:00"];
	pbdb_date[ndates2][pbdb_date[ndates2]=="0000-00-00 00:00:00"] <- pbdb_data$created[ndates2[pbdb_date[ndates2]=="0000-00-00 00:00:00"]];
	pbdb_date[ndates3] <- as.Date.character(as.Date.character(pbdb_data$modified[ndates3]));

#		as.numeric(pbdb_date[!is.na(as.numeric(pbdb_date))]);
#	pbdb_date2 <- as.character(pbdb_date[is.na(as.numeric(pbdb_date))]);
#	if (!is.na(as.numeric(pbdb_date)))	{
#		pbdb_data$created <- as.Date.numeric(pbdb_date,origin="1970-01-01");
#		} else	{
#		pbdb_data$created <- as.Date.character(pbdb_data$created);
#		}
#	pbdb_date <- pbdb_data$modified;
#	if (!is.na(as.numeric(pbdb_date[1])))	{
#		pbdb_data$modified <- as.Date.numeric(as.numeric(pbdb_date),origin="1970-01-01");
#		} else	{
#		pbdb_data$modified <- as.Date.character(pbdb_data$modified);
#		}
	}
options(warn=1);
return(pbdb_data);
}

# modified 2021-05-10
convert_pbdb_date <- function(pbdb_date)	{
if (gsub("\\/","-",as.character(pbdb_date))!=as.character(pbdb_date))	{
	entry_date <- strsplit(pbdb_date," ")[[1]][1];
	entry_time <- strsplit(pbdb_date," ")[[1]][2];
	entry_dates <- as.numeric(strsplit(entry_date,"\\/")[[1]]);
	if (entry_dates[3] %in% 80:99)	{
		entry_dates[3] <- 1900+entry_dates[3];
		} else	{
		entry_dates[3] <- 2000+entry_dates[3];
		}
	if(length(strsplit(entry_time,":")[[1]])==2) entry_time <- paste(entry_time,":00",sep="");
	entry_times <- strsplit(entry_time,":")[[1]];
	if (is.na(entry_times[1]))	{
		return(paste(entry_dates[3],entry_dates[1],entry_dates[2],sep="-"));
		} else	{
		return(ISOdate(entry_dates[3],entry_dates[1],entry_dates[2],entry_times[1],entry_times[2],entry_times[3]));
		}
	} else if (gsub("-"," ",as.character(pbdb_date))!=as.character(pbdb_date))	{
	entry_dates <- as.numeric(strsplit(strsplit(as.character(pbdb_date)," ")[[1]][1],"-")[[1]]);
	if (length(strsplit(as.character(pbdb_date)," ")[[1]])==2)	{
		entry_time <- strsplit(as.character(pbdb_date)," ")[[1]][2];
		entry_times <- strsplit(entry_time,":")[[1]];
		return(ISOdate(entry_dates[3],entry_dates[2],entry_dates[1],entry_times[1],entry_times[2],entry_times[3]));
		} else	{
		return(ISOdate(entry_dates[1],entry_dates[2],entry_dates[3]));
		}
	} else if (!is.na(as.numeric(pbdb_date)))	{
	return(as.Date(as.numeric(pbdb_date),origin="1970-01-01"));
	} else {
	print(pbdb_date);
	}
}

# get occurrences for a taxon from some span of time & environment
# modified 2020-03-02
# get occurrences for a taxon from some span of time & environment
# modified 2020-03-02
# modified 2020-05-05
# modified 2021-07-09
accersi_occurrence_data <- function(taxa,onset="Archean",end="Holocene",basic_environments="terr,marine,unknown",species_only=T,clean_entered_taxa=T,directory="",save_files=T,output_type=".csv") {
# Arguments:
# 	taxa: proper taxonomic name
# 	onset: onset geological interval from which you want new records
# 	end: end geological interval from which you want new records
# 	basic_environments: environments to download ("terr,marine,unknown" for "terrestrial, marine, unknown")
# 	species_only: if true, then eliminate Genus sp. identifications.
#	clean_entered_taxa: remove tags from entered taxa. (These always are removed from accepted names)
#	directory: where to send output files
#	save_files: if true, then output final results
# 	file_format: the end tag on the output files: '.xls' for Excel, '.txt', '.csv" for comma-delimited
taxa <- paste(taxa, collapse = ",");
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
taxa <- gsub(" ","%20",taxa);
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
all_finds <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#fetch <- RCurl::getURL(http);
fetched <- gsub("\"","",simplify2array(strsplit(RCurl::getURL(http),"\r\n"))[,1]);
if (!is.na(match("THIS REQUEST RETURNED NO RECORDS",fetched)))	{
	return(dummy_finds);
	} else	{
#	all_finds <- utils::read.csv(text = fetch, header = FALSE, stringsAsFactors=hell_no);
	if (all_finds[1,1]=="Warning:")	{
		# this will happen only if something goes badly wrong.
#		if (sum(strsplit(all_finds[1,2]," ")[[1]] %in% c("did","not","match","any","name","in","the","taxonomy","table"))>=length(c("did","not","match","any","name","in","the","taxonomy","table")))	{
#			unentered <- T;
#			} else	{
#			entered <- F;
#			}
		cc <- match("occurrence_no",all_finds[,1])
		kluge_mc_kluge_face <- character();
		for (mm in 1:ncol(all_finds))
			kluge_mc_kluge_face <- c(kluge_mc_kluge_face,as.character(all_finds[cc,mm]));
		colnames(all_finds) <- kluge_mc_kluge_face;
		all_finds <- all_finds[(cc+1):nrow(all_finds),];
		all_finds$accepted_name[all_finds$accepted_name==""] <- all_finds$identified_name[all_finds$accepted_name==""];
		all_finds$accepted_rank[all_finds$accepted_rank==""] <- all_finds$identified_rank[all_finds$accepted_rank==""];
		taxon_name=all_finds$accepted_name[all_finds$genus==""];
		all_finds$genus[all_finds$genus==""] <- sapply(taxon_name,divido_genus_names_from_species_names)
		} #else	{
#		all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
#		}
	all_finds <- expello_na_from_matrix(data=all_finds,replacement = "");
	if (species_only)	{
#		xxx <- (1:nrow(all_finds))[all_finds$identified_rank %in% c("species","subspecies")]
		desired_finds <- all_finds[all_finds$identified_rank %in% c("species","subspecies"),];
		desired_finds <- subset(desired_finds,desired_finds$genus_no>0);
#		desired_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
		}	else	{
		desired_finds <- all_finds;
		}
	if (nrow(desired_finds)==0)	{
		return(desired_finds);
		} else	{
		noccr <- nrow(desired_finds);
		print("Cleaning misentered taxonomic uncertainties.")
		taxon_names <- desired_finds$identified_name[desired_finds$identified_rank %in% c("species","subspecies")];
		for (tn in 1:length(taxon_names))	taxon_names[tn] <- repair_misentered_uncertain_species(taxon_names[tn]);
		desired_finds$identified_name[desired_finds$identified_rank %in% c("species","subspecies")] <- unlist(pbapply::pbsapply(taxon_names,repair_misentered_uncertain_species));
		taxon_name <- desired_finds$identified_name;
		print("Flagging uncertain taxonomic assignments.");
		taxon_names <- desired_finds$identified_name;
		desired_finds$flags <- pbapply::pbsapply(taxon_names,identify_taxonomic_uncertainty);
#		taxon_name <- desired_finds$identified_name;
#		print("Separating uncertain species assignments.")
#		flags1 <- pbapply::pbsapply(taxon_name,revelare_uncertain_species_assignments);
#		flags3 <- flags2 <- rep("",noccr);
#		taxon_name <- desired_finds$identified_name[desired_finds$accepted_rank %in% c("genus","subgenus")];
#		if (length(taxon_name)>0)	{
#			flags2[desired_finds$accepted_rank %in% c("genus","subgenus")] <- sapply(taxon_name,revelare_uncertain_genus_assignments);
#			double <- (1:noccr)[flags1!=""][(1:noccr)[flags1!=""] %in% (1:noccr)[flags2!=""]];
#			flags3[flags1!=""] <- flags1[flags1!=""];
#			flags3[flags2!=""] <- flags2[flags2!=""];
#			flags3[double] <- paste(unique(flags2[flags2!=""]),unique(flags1[flags1!=""]),sep=", ");
#			desired_finds$flags <- simplify2array(flags3);
#			} else	{
#			desired_finds$flags <- flags1;
#			}

		# use flags field to note uncertain genus or species assignments.
#		desired_finds$flags <- sapply(taxon_name,identify_taxonomic_uncertainty);
#		for (tn in 1:nrow(desired_finds))	{
#			desired_finds$flags[tn] <- identify_taxonomic_uncertainty(taxon_name=desired_finds$identified_name[tn]);
#			}
		# some uncertain genus assignments will be clarified by updated generic assignments
#		uncertain_genera <- (1:noccr)[desired_finds$flags=="uncertain genus"];
#		recertain_genera <- uncertain_genera[desired_finds$difference[uncertain_genera]=="recombined as"]
#		desired_finds$flags[recertain_genera] <- "";

		if (clean_entered_taxa)	{
			print("Cleaning entered names...")
			cleaned_names <- pbapply::pbsapply(as.character(desired_finds$identified_name),mundify_taxon_names);
			desired_finds$identified_name <- cleaned_names;
			}	# removes tags such as "cf.", "?", "n. sp." from entered names

		if (species_only)	{
			entered_species <- sort(c((1:noccr)[desired_finds$accepted_rank=="species"],(1:noccr)[desired_finds$accepted_rank=="subspecies"]));
			unentered_species <- (1:noccr)[!(1:noccr) %in% entered_species];
			desired_finds$accepted_name[unentered_species] <- desired_finds$identified_name[unentered_species];
			noccr <- nrow(desired_finds);
			}
		print("Cleaning accepted names...")
		cleaned_names <- pbapply::pbsapply(as.character(desired_finds$accepted_name),mundify_taxon_names);
		desired_finds$accepted_name <- cleaned_names;

#		desired_finds$genus==""
#		taxon_names <- desired_finds$accepted_name;
#		sapply(taxon_names,divido_genus_names_from_species_names)
		# make sure that type subgenera are consistently Lophospira (Lophospira)
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		print("Separating subgenera.")
		genus_subgenus <- base::t(pbapply::pbsapply(genus_name,divido_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		s_g <- nrow(subgenera);
		if (is.null(s_g))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				s_g <- 1;
				} else	{
				subgenera <- matrix(data=0,nrow=0,ncol=2);
				s_g <- 0;
				}
			}
		type_subgenus <- subgenera[subgenera[,2]==subgenera[,1],];
		# first do non-type subgenera (e.g., Lophospira (Ruedemannia))
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(data=0,nrow=0,ncol=2);
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (!is.na(match(nontype_subgenus[sg,2],genus_subgenus[,1])))	{
	#			doubly_ranked <- genus_name[match(nontype_subgenus[sg,2],genus_subgenus[,1])]
	#			print(nontype_subgenus[sg,2]);
				doubly_ranked <- paste(nontype_subgenus[sg,1]," (",nontype_subgenus[sg,2],")",sep="");
				subgenus_finds <- (1:noccr)[desired_finds$genus==doubly_ranked];
				desired_finds$genus[subgenus_finds] <- nontype_subgenus[sg,2];
				nn <- match(nontype_subgenus[sg,2],genus_name);
				if (is.na(nn))	{
					# if name never is used alone, then add it to the genus lists
					nn <- match(nontype_subgenus[sg,2],genus_subgenus[,1]);
					genus_name <- insert_cell_into_vector_x(x=genus_name,new_value=as.character(genus_subgenus[nn,1]),cell_no=nn);
					genus_subgenus <- insert_row_into_matrix_x(x=genus_subgenus,new_row=c(genus_subgenus[nn,1],""),row_no=nn);
					}
				}
			}

		# now do type subgenera (e.g., Lophospira (Lophospira))
		t_s_g <- nrow(type_subgenus);
		if (is.null(t_s_g))	{
			if (length(type_subgenus)==2)	{
				t_s_g <- 1;
				type_subgenus <- matrix(data=type_subgenus,nrow=1,ncol=2)
				} else	{
				t_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < t_s_g)	{
			sg <- sg+1;
			subgenus_finds <- (1:noccr)[desired_finds$genus==type_subgenus[sg,1]];
			desired_finds$genus[subgenus_finds] <- paste(type_subgenus[sg,1]," (",type_subgenus[sg,2],")",sep="");
			}

		# now, make sure that subgenera are not listed in two different genera
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		genus_subgenus <- base::t(sapply(genus_name,divido_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		if (is.null(nrow(subgenera)))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				} else	{
				subgenera <- matrix(data=0,nrow=1,ncol=2);
				}
			}
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		if (is.null(nrow(nontype_subgenus)))	{
			if (length(nontype_subgenus)==2)	{
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(0,0,2);
				}
			}
		if (nrow(nontype_subgenus)>0)
			nontype_subgenus <- nontype_subgenus[nontype_subgenus[,2]!="",];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2)
				} else	{
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (sum(nontype_subgenus[sg,2]==nontype_subgenus[,2])>1)	{
				sgs <- (1:nt_s_g)[nontype_subgenus[,2] %in% nontype_subgenus[sg,2]]
				ssg <- 1;
				senior_entry <- paste(nontype_subgenus[sgs[1],1]," (",nontype_subgenus[sgs[1],2],")",sep="");
				while (ssg < sum(nontype_subgenus[sg,2]==nontype_subgenus[,2]))	{
					ssg <- ssg+1;
					double_entry <- paste(nontype_subgenus[sgs[ssg],1]," (",nontype_subgenus[sgs[ssg],2],")",sep="");
					desired_finds$genus[(1:noccr)[desired_finds$genus==double_entry]] <- senior_entry;
					}
				}
			}

		taxon_name <- desired_finds$accepted_name;
		species_epithet <- sapply(taxon_name,divido_species_epithets);
		desired_finds$accepted_name_orig <- desired_finds$accepted_name;
		desired_finds$accepted_name <- paste(desired_finds$genus,species_epithet);

		desired_finds <- expello_na_from_matrix(desired_finds,replacement="");

		if (save_files)	{
			taxa <- gsub(",","+",taxa);
			if (onset!=end)	{
				timespan <- paste(onset,"-",end,sep="");
				}	else	timespan <- onset;
			output <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="");
			if (directory!="")
				output <- paste(directory,output,sep="");
			if (output_type==".csv")	{
				output <- gsub("TRUE","",output);
				write.csv(desired_finds,file=output,row.names = FALSE,fileEncoding = "UTF-8");
				}	else	{
				output <- gsub("TRUE","",output);
				write.table(desired_finds,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
				}
			}
		}
	for (cn in 1:ncol(desired_finds))	{
		old_info <- desired_finds[,cn];
		cnm <- strsplit(x=colnames(desired_finds)[cn],split="_")[[1]];
		if(cnm[length(cnm)] %in% paleodb_numeric_fields)	{
			old_info[old_info %in% missing_data_assignment] <- 0;
			desired_finds[,cn] <- as.numeric(old_info);
			} else	{
			desired_finds[,cn] <- as.character(old_info);
			}
		}
	return(desired_finds);
	}
}

# get occurrences modified after a certain date. If the date is blank, then it will just be the last 5 years.
update_occurrence_data <- function(taxa="Life",onset="Archean",end="Holocene",basic_environments="terr,marine,unknown",occs_modified_after="1900-01-01",species_only=T,clean_entered_taxa=T,directory="",save_files=T,output_type=".csv"){
if (occs_modified_after=="1900-01-01")
	occs_modified_after <- strsplit(as.character(Sys.time()-5*366*24*60*60),split=" ")[[1]][1];
taxa <- paste(taxa, collapse = ",");
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
taxa <- gsub(" ","%20",taxa);
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&occs_modified_after=",occs_modified_after,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
all_finds <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
fetched <- gsub("\"","",simplify2array(strsplit(RCurl::getURL(http),"\r\n"))[,1]);
if (!is.na(match("THIS REQUEST RETURNED NO RECORDS",fetched)))	{
	return(dummy_finds);
	} else	{
#	all_finds <- utils::read.csv(text = fetch, header = FALSE, stringsAsFactors=hell_no);
	if (all_finds[1,1]=="Warning:")	{
		# this will happen only if something goes badly wrong.
#		if (sum(strsplit(all_finds[1,2]," ")[[1]] %in% c("did","not","match","any","name","in","the","taxonomy","table"))>=length(c("did","not","match","any","name","in","the","taxonomy","table")))	{
#			unentered <- T;
#			} else	{
#			entered <- F;
#			}
		cc <- match("occurrence_no",all_finds[,1])
		kluge_mc_kluge_face <- character();
		for (mm in 1:ncol(all_finds))
			kluge_mc_kluge_face <- c(kluge_mc_kluge_face,as.character(all_finds[cc,mm]));
		colnames(all_finds) <- kluge_mc_kluge_face;
		all_finds <- all_finds[(cc+1):nrow(all_finds),];
		all_finds$accepted_name[all_finds$accepted_name==""] <- all_finds$identified_name[all_finds$accepted_name==""];
		all_finds$accepted_rank[all_finds$accepted_rank==""] <- all_finds$identified_rank[all_finds$accepted_rank==""];
		taxon_name=all_finds$accepted_name[all_finds$genus==""];
		all_finds$genus[all_finds$genus==""] <- sapply(taxon_name,divido_genus_names_from_species_names)
		} #else	{
#		all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
#		}
	all_finds <- expello_na_from_matrix(data=all_finds,replacement = "");
	if (species_only)	{
		xxx <- (1:nrow(all_finds))[all_finds$identified_rank %in% c("species","subspecies")]
		desired_finds <- all_finds[xxx,];
		desired_finds <- subset(desired_finds,desired_finds$genus_no>0);
#		desired_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
		}	else	{
		desired_finds <- all_finds;
		}
	if (nrow(desired_finds)==0)	{
		return(desired_finds);
		} else	{
		noccr <- nrow(desired_finds);
		print("Cleaning misentered taxonomic uncertainties.")
		taxon_names <- desired_finds$identified_name[desired_finds$identified_rank %in% c("species","subspecies")];
		desired_finds$identified_name[desired_finds$identified_rank %in% c("species","subspecies")] <- unlist(pbapply::pbsapply(taxon_names,repair_misentered_uncertain_species));
		taxon_name <- desired_finds$identified_name;
		print("Flagging uncertain taxonomic assignments.");
		taxon_names <- desired_finds$identified_name;
		desired_finds$flags <- pbapply::pbsapply(taxon_names,identify_taxonomic_uncertainty);

#		flags1 <- pbapply::pbsapply(taxon_name,revelare_uncertain_species_assignments);
#		flags3 <- flags2 <- rep("",noccr);
#		taxon_name <- desired_finds$identified_name[desired_finds$accepted_rank %in% c("genus","subgenus")];
#		if (length(taxon_name)>0)	{
#			flags2[desired_finds$accepted_rank %in% c("genus","subgenus")] <- sapply(taxon_name,revelare_uncertain_genus_assignments);
#			double <- (1:noccr)[flags1!=""][(1:noccr)[flags1!=""] %in% (1:noccr)[flags2!=""]];
#			flags3[flags1!=""] <- flags1[flags1!=""];
#			flags3[flags2!=""] <- flags2[flags2!=""];
#			flags3[double] <- paste(unique(flags2[flags2!=""]),unique(flags1[flags1!=""]),sep=", ");
#			desired_finds$flags <- simplify2array(flags3);
#			} else	{
#			desired_finds$flags <- flags1;
#			}

		# use flags field to note uncertain genus or species assignments.
#		desired_finds$flags <- sapply(taxon_name,identify_taxonomic_uncertainty);
#		for (tn in 1:nrow(desired_finds))	{
#			desired_finds$flags[tn] <- identify_taxonomic_uncertainty(taxon_name=desired_finds$identified_name[tn]);
#			}
		# some uncertain genus assignments will be clarified by updated generic assignments
#		uncertain_genera <- (1:noccr)[desired_finds$flags=="uncertain genus"];
#		recertain_genera <- uncertain_genera[desired_finds$difference[uncertain_genera]=="recombined as"]
#		desired_finds$flags[recertain_genera] <- "";

		if (clean_entered_taxa)	{
			print("Cleaning entered names...")
			cleaned_names <- pbapply::pbsapply(as.character(desired_finds$identified_name),mundify_taxon_names);
			desired_finds$identified_name <- cleaned_names;
			}	# removes tags such as "cf.", "?", "n. sp." from entered names

		if (species_only)	{
			entered_species <- sort(c((1:noccr)[desired_finds$accepted_rank=="species"],(1:noccr)[desired_finds$accepted_rank=="subspecies"]));
			unentered_species <- (1:noccr)[!(1:noccr) %in% entered_species];
			desired_finds$accepted_name[unentered_species] <- desired_finds$identified_name[unentered_species];
			noccr <- nrow(desired_finds);
			}
		print("Cleaning accepted names...")
		cleaned_names <- pbapply::pbsapply(as.character(desired_finds$accepted_name),mundify_taxon_names);
		desired_finds$accepted_name <- cleaned_names;

#		desired_finds$genus==""
#		taxon_names <- desired_finds$accepted_name;
#		sapply(taxon_names,divido_genus_names_from_species_names)
		# make sure that type subgenera are consistently Lophospira (Lophospira)
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		print("Separating subgenera.")
		genus_subgenus <- base::t(pbapply::pbsapply(genus_name,divido_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		s_g <- nrow(subgenera);
		if (is.null(s_g))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				s_g <- 1;
				} else	{
				subgenera <- matrix(data=0,nrow=0,ncol=2);
				s_g <- 0;
				}
			}
		type_subgenus <- subgenera[subgenera[,2]==subgenera[,1],];
		# first do non-type subgenera (e.g., Lophospira (Ruedemannia))
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(data=0,nrow=0,ncol=2);
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (!is.na(match(nontype_subgenus[sg,2],genus_subgenus[,1])))	{
	#			doubly_ranked <- genus_name[match(nontype_subgenus[sg,2],genus_subgenus[,1])]
	#			print(nontype_subgenus[sg,2]);
				doubly_ranked <- paste(nontype_subgenus[sg,1]," (",nontype_subgenus[sg,2],")",sep="");
				subgenus_finds <- (1:noccr)[desired_finds$genus==doubly_ranked];
				desired_finds$genus[subgenus_finds] <- nontype_subgenus[sg,2];
				nn <- match(nontype_subgenus[sg,2],genus_name);
				if (is.na(nn))	{
					# if name never is used alone, then add it to the genus lists
					nn <- match(nontype_subgenus[sg,2],genus_subgenus[,1]);
					genus_name <- insert_cell_into_vector_x(x=genus_name,new_value=as.character(genus_subgenus[nn,1]),cell_no=nn);
					genus_subgenus <- insert_row_into_matrix_x(x=genus_subgenus,new_row=c(genus_subgenus[nn,1],""),row_no=nn);
					}
				}
			}

		# now do type subgenera (e.g., Lophospira (Lophospira))
		t_s_g <- nrow(type_subgenus);
		if (is.null(t_s_g))	{
			if (length(type_subgenus)==2)	{
				t_s_g <- 1;
				type_subgenus <- matrix(data=type_subgenus,nrow=1,ncol=2)
				} else	{
				t_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < t_s_g)	{
			sg <- sg+1;
			subgenus_finds <- (1:noccr)[desired_finds$genus==type_subgenus[sg,1]];
			desired_finds$genus[subgenus_finds] <- paste(type_subgenus[sg,1]," (",type_subgenus[sg,2],")",sep="");
			}

		# now, make sure that subgenera are not listed in two different genera
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		genus_subgenus <- base::t(sapply(genus_name,divido_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		if (is.null(nrow(subgenera)))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				} else	{
				subgenera <- matrix(data=0,nrow=1,ncol=2);
				}
			}
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		if (is.null(nrow(nontype_subgenus)))	{
			if (length(nontype_subgenus)==2)	{
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(0,0,2);
				}
			}
		if (nrow(nontype_subgenus)>0)
			nontype_subgenus <- nontype_subgenus[nontype_subgenus[,2]!="",];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2)
				} else	{
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (sum(nontype_subgenus[sg,2]==nontype_subgenus[,2])>1)	{
				sgs <- (1:nt_s_g)[nontype_subgenus[,2] %in% nontype_subgenus[sg,2]]
				ssg <- 1;
				senior_entry <- paste(nontype_subgenus[sgs[1],1]," (",nontype_subgenus[sgs[1],2],")",sep="");
				while (ssg < sum(nontype_subgenus[sg,2]==nontype_subgenus[,2]))	{
					ssg <- ssg+1;
					double_entry <- paste(nontype_subgenus[sgs[ssg],1]," (",nontype_subgenus[sgs[ssg],2],")",sep="");
					desired_finds$genus[(1:noccr)[desired_finds$genus==double_entry]] <- senior_entry;
					}
				}
			}

		taxon_name <- desired_finds$accepted_name;
		species_epithet <- sapply(taxon_name,divido_species_epithets);
		desired_finds$accepted_name_orig <- desired_finds$accepted_name;
		desired_finds$accepted_name <- paste(desired_finds$genus,species_epithet);

		desired_finds <- expello_na_from_matrix(desired_finds,replacement="");

		if (save_files)	{
			taxa <- gsub(",","+",taxa);
			if (onset!=end)	{
				timespan <- paste(onset,"-",end,sep="");
				}	else	timespan <- onset;
			output <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="");
			if (directory!="")
				output <- paste(directory,output,sep="");
			if (output_type==".csv")	{
				output <- gsub("TRUE","",output);
				write.csv(desired_finds,file=output,row.names = FALSE,fileEncoding = "UTF-8");
				}	else	{
				output <- gsub("TRUE","",output);
				write.table(desired_finds,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
				}
			}
		}
	for (cn in 1:ncol(desired_finds))	{
		old_info <- desired_finds[,cn];
		cnm <- strsplit(x=colnames(desired_finds)[cn],split="_")[[1]];
		if(cnm[length(cnm)] %in% paleodb_numeric_fields)	{
			old_info[old_info %in% missing_data_assignment] <- 0;
			desired_finds[,cn] <- as.numeric(old_info);
			} else	{
			desired_finds[,cn] <- as.character(old_info);
			}
		}
	return(desired_finds);
	}
}

# get occurrences from a single locality (coll_id = PaleoDB collection_no)
accersi_occurrences_from_one_paleodb_collection <- function(coll_id)	{
# coll_id: PaleoDB collection_no
# coll_id <- 12408;
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?coll_id=",coll_id,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
fetch <- RCurl::getURL(http);
coll_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding = "UTF-8");
taxon_names <- coll_finds$identified_name[coll_finds$identified_rank %in% c("species","subspecies")];
#taxon_names_fixed <- unlist(pbapply::pbsapply(taxon_names,repair_misentered_uncertain_species));
#coll_finds$identified_name[coll_finds$identified_rank %in% c("species","subspecies")] <- taxon_names_fixed;
coll_finds$identified_name[coll_finds$identified_rank %in% c("species","subspecies")] <- unlist(pbapply::pbsapply(taxon_names,repair_misentered_uncertain_species));
taxon_names <- coll_finds$identified_name;
coll_finds$flags <- pbapply::pbsapply(taxon_names,identify_taxonomic_uncertainty);

return(coll_finds);
}

accersi_occurrences_from_file_of_taxa <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- read.table(file=taxonfilename,stringsAsFactors = TRUE,sep="\t")[,1]
taxon_list <- sapply(taxon_name,mundify_taxon_names)
occompendium <- accersi_occurrences_for_list_of_taxa(taxon_list,analysis_name,output_file)
return(occompendium)
}

accersi_occurrences_for_one_taxon <- function(taxon)	{
# taxon: taxon name (can be any rank)
# lump_subgenera: if TRUE and taxon is a genus, then subgeneneric occurrences are included;
#	if false, then subgenus occurrences are culled.
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
all_finds <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#fetch <- RCurl::getURL(http);
fetched <- gsub("\"","",simplify2array(strsplit(RCurl::getURL(http),"\r\n"))[,1]);
if (!is.na(match("THIS REQUEST RETURNED NO RECORDS",fetched)))	{
	return(dummy_finds2);
	} else	{
#	all_finds <- utils::read.csv(text = fetch, header = FALSE, stringsAsFactors=hell_no);
	if (all_finds[1,1]=="Warning:")	{
		cc <- match("occurrence_no",all_finds[,1])
		kluge_mc_kluge_face <- character();
		for (mm in 1:ncol(all_finds))	kluge_mc_kluge_face <- c(kluge_mc_kluge_face,as.character(all_finds[cc,mm]));
		colnames(all_finds) <- kluge_mc_kluge_face;
		all_finds <- all_finds[(cc+1):nrow(all_finds),];
		all_finds$accepted_name[all_finds$accepted_name==""] <- all_finds$identified_name[all_finds$accepted_name==""];
		all_finds$accepted_rank[all_finds$accepted_rank==""] <- all_finds$identified_rank[all_finds$accepted_rank==""];
		taxon_name=all_finds$accepted_name[all_finds$genus==""];
		all_finds$genus[all_finds$genus==""] <- sapply(taxon_name,divido_genus_names_from_species_names)
		} #else	{
	all_finds <- expello_na_from_matrix(data=all_finds,replacement = "");
	#http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxon,"&show=full",sep="");
	#taxon_finds <- accersi_occurrence_data(taxa=taxon,species_only=species_only);
	taxon_name <- all_finds$accepted_name;
	species_epithet <- sapply(taxon_name,divido_species_epithets);
	all_finds$accepted_name_orig <- all_finds$accepted_name;
	all_finds$accepted_name <- paste(all_finds$genus,species_epithet);
	return(all_finds);
	}
}

accersi_occurrences_for_one_taxon_no <- function(taxon_no)	{
# taxon: taxon name (can be any rank)
# lump_subgenera: if TRUE and taxon is a genus, then subgeneneric occurrences are included;
#	if false, then subgenus occurrences are culled.
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?taxon_no=",taxon_no,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
all_finds <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#fetch <- RCurl::getURL(http);
fetched <- gsub("\"","",simplify2array(strsplit(RCurl::getURL(http),"\r\n"))[,1]);
if (!is.na(match("THIS REQUEST RETURNED NO RECORDS",fetched)))	{
	return(dummy_finds2);
	} else	{
#	all_finds <- utils::read.csv(text = fetch, header = FALSE, stringsAsFactors=hell_no);
	if (all_finds[1,1]=="Warning:")	{
		cc <- match("occurrence_no",all_finds[,1])
		kluge_mc_kluge_face <- character();
		for (mm in 1:ncol(all_finds))	kluge_mc_kluge_face <- c(kluge_mc_kluge_face,as.character(all_finds[cc,mm]));
		colnames(all_finds) <- kluge_mc_kluge_face;
		all_finds <- all_finds[(cc+1):nrow(all_finds),];
		all_finds$accepted_name[all_finds$accepted_name==""] <- all_finds$identified_name[all_finds$accepted_name==""];
		all_finds$accepted_rank[all_finds$accepted_rank==""] <- all_finds$identified_rank[all_finds$accepted_rank==""];
		taxon_name=all_finds$accepted_name[all_finds$genus==""];
		all_finds$genus[all_finds$genus==""] <- sapply(taxon_name,divido_genus_names_from_species_names)
		} #else	{
	all_finds <- expello_na_from_matrix(data=all_finds,replacement = "");
	#http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxon,"&show=full",sep="");
	#taxon_finds <- accersi_occurrence_data(taxa=taxon,species_only=species_only);
	taxon_name <- all_finds$accepted_name;
	species_epithet <- sapply(taxon_name,divido_species_epithets);
	all_finds$accepted_name_orig <- all_finds$accepted_name;
	all_finds$accepted_name <- paste(all_finds$genus,species_epithet);
	return(all_finds);
	}
}

accersi_occurrences_for_list_of_taxa <- function(taxon_list,lump_subgenera=F,species_only=T,paleogeography="gplates")	{
# remove any funny symbols from taxon names
taxon_list <- sapply(taxon_list,mundify_taxon_names);
ntaxa <- length(taxon_list);

#occurrences_compendium_list <- sapply(taxa,accersi_occurrences_for_one_taxon,species_only);
#occurrences_compendium_list <- base::t(sapply(taxa,accersi_occurrence_data,species_only=species_only,save_files=F));
occurrences_compendium <- c();
#tx <- 0;
for (tx in 1:ntaxa)	{
#	taxon_finds <- occurrences_compendium_list[[tx]];
#	tx <- tx+1;
#	taxa <- taxon_list[tx];
	taxon_finds <- accersi_occurrence_data(taxa=taxon_list[tx],species_only=species_only,save_files=F);
	if (!is.null(taxon_finds) && (nrow(taxon_finds)==0 && species_only))
		# we found nothing; let's see if we can find non-species finds
		taxon_finds <- accersi_occurrence_data(taxa=taxon_list[tx],species_only=F,save_files=F);
	if (!lump_subgenera && !is.null(nrow(taxon_finds)))	{
		taxon_info <- accersi_taxonomic_data_for_one_taxon(taxon=taxon_list[tx]);
		if (nrow(taxon_info)==1 || taxon_info[2,1]!="THIS REQUEST RETURNED NO RECORDS")	{
			backup_taxon_finds <- taxon_finds;
			if (taxon_info$taxon_rank=="genus")	{
				this_genus <- (1:nrow(taxon_finds))[taxon_finds$genus %in% c(taxon_list[tx],paste(taxon_list[tx]," (",taxon_list[tx],")",sep=""))];
				taxon_finds <- taxon_finds[this_genus,];
				# this kluge protects against wonky cases where the PaleoDB has conflicted information about genus/subgenus status
				if (nrow(taxon_finds)==0)	{
					taxon_finds <- backup_taxon_finds;
					genus_name <- taxon_finds$genus
					genus_subgenus <- sapply(genus_name,divido_subgenus_names_from_genus_names);
					taxon_finds$genus <- genus_subgenus[2,];
					this_genus <- (1:nrow(taxon_finds))[taxon_finds$genus %in% c(taxon_list[tx],paste(taxon_list[tx]," (",taxon_list[tx],")",sep=""))];
					taxon_finds <- taxon_finds[this_genus,];
					}
				}
			}
		}
	if (!is.null(nrow(taxon_finds)))	{
		taxon <- rep(taxon_list[tx],nrow(taxon_finds));
		taxon_finds <- cbind(taxon=as.character(taxon),taxon_finds);
		if (tx==1)	{
			occurrences_compendium <- taxon_finds;
			} else	{
			occurrences_compendium <- rbind(occurrences_compendium,taxon_finds);
			}
		}
#	print(dim(taxon_finds));
	}
print("finished initial run of taxa")
if (!is.null(occurrences_compendium))	{
	coll_id <- collection_no <- sort(unique(occurrences_compendium$collection_no));
	print("Getting collection data...")
	pbdb_list <- data.frame(base::t(pbapply::pbsapply(coll_id,accersi_data_for_one_collection)));
	collection_compendium <- list_to_dataframe_for_pbdb_data(pbdb_list);
#	c <- 0;
#	collection_compendium <- c();
#	while (c < length(collection_no))	{
#		c <- c+1;
#	for (c in 1:length(collection_no))	{
#		if (is.null(collection_compendium))	{
#	if (c==1)	{
#			collection_compendium <- accersi_single_locality_info(collection_no=collection_no[c],paleogeography=paleogeography);
#			} else	{
#			collection_compendium <- rbind(collection_compendium,accersi_single_locality_info(collection_no=collection_no[c],paleogeography=paleogeography));
#			}
#		}

	named_rock_units <- collection_compendium$formation;
	print("cleaning formations...");
	collection_compendium$formation <- pbapply::pbsapply(named_rock_units,mundify_rock_unit_names);
	named_rock_units <- collection_compendium$member;
	print("cleaning members...");
	collection_compendium$member <- pbapply::pbsapply(named_rock_units,mundify_rock_unit_names);
	named_rock_units <- collection_compendium$stratgroup;
	print("cleaning groups...");
	collection_compendium$stratgroup <- pbapply::pbsapply(named_rock_units,mundify_rock_unit_names);
	zone <- collection_compendium$zone[collection_compendium$zone!=""];
	print("cleaning zones...");
	collection_compendium$zone[collection_compendium$zone!=""] <- pbapply::pbsapply(zone,mundus_zone);
	web_text <- collection_compendium$collection_name;
	print("cleaning collection names...");
	collection_compendium$collection_name <- pbapply::pbsapply(web_text,mundify_web_text_dull);
	print("cleaning stratigraphic comments...");
	web_text <- collection_compendium$stratcomments;
	collection_compendium$stratcomments <- pbapply::pbsapply(web_text,mundify_web_text_dull);

	output <- list(collection_compendium,occurrences_compendium);
	}	else	{
	output <- list("","");
	}
names(output) <- c("collection_compendium","occurrences_compendium");
return(output);
}

update_occurrence_from_occurrence_no <- function(occurrence_no)	{
# remove any funny symbols from taxon names
http <- paste("https://paleobiodb.org/data1.2/occs/single.csv?id=",occurrence_no,"&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
updated_find <- read.csv(http,header=T,stringsAsFactors = hell_no);
return(updated_find);
}

accersi_occurrences_for_list_of_taxa_old <- function(taxon_list,lump_subgenera=F,species_only=T)	{
# remove any funny symbols from taxon names
taxon_list <- sapply(taxon_list,mundify_taxon_names);
ntaxa <- length(taxon_list);

#occurrences_compendium_list <- sapply(taxa,accersi_occurrences_for_one_taxon,species_only);
#occurrences_compendium_list <- base::t(sapply(taxa,accersi_occurrence_data,species_only=species_only,save_files=F));
occurrences_compendium <- c();
for (tx in 1:ntaxa)	{
#	taxon_finds <- occurrences_compendium_list[[tx]];
#	tx <- tx+1;
	taxa <- taxon_list[tx];
	taxon_finds <- accersi_occurrence_data(taxa,species_only=species_only,save_files=F);
	if (!is.null(taxon_finds) && (nrow(taxon_finds)==0 && species_only))
		taxon_finds <- accersi_occurrence_data(taxa,species_only=F,save_files=F);
	if (!lump_subgenera && !is.null(nrow(taxon_finds)))	{
		taxon_info <- accersi_taxonomic_data_for_one_taxon(taxon=taxa);
		if (taxon_info$taxon_rank=="genus")	{
			this_genus <- (1:nrow(taxon_finds))[taxon_finds$genus %in% c(taxa,paste(taxa," (",taxa,")",sep=""))];
			taxon_finds <- taxon_finds[this_genus,];
			}
		}
	if (!is.null(nrow(taxon_finds)))	{
		taxon <- rep(taxa,nrow(taxon_finds));
		taxon_finds <- cbind(taxon=as.character(taxon),taxon_finds);
		if (tx==1)	{
			occurrences_compendium <- taxon_finds;
			} else	{
			occurrences_compendium <- rbind(occurrences_compendium,taxon_finds);
			}
		}
	}

collection_no <- sort(unique(occurrences_compendium$collection_no));
for (c in 1:length(collection_no))	{
	if (c==1)	{
		collection_compendium <- accersi_single_locality_info(collection_no=collection_no[c]);
		} else	{
		collection_compendium <- rbind(collection_compendium,accersi_single_locality_info(collection_no=collection_no[c]));
		}
	}

named_rock_units <- collection_compendium$formation;
collection_compendium$formation <- sapply(named_rock_units,mundify_rock_unit_names);
named_rock_units <- collection_compendium$member;
collection_compendium$member <- sapply(named_rock_units,mundify_rock_unit_names);
named_rock_units <- collection_compendium$stratgroup;
collection_compendium$stratgroup <- sapply(named_rock_units,mundify_rock_unit_names);
zone <- collection_compendium$zone[collection_compendium$zone!=""];
collection_compendium$zone[collection_compendium$zone!=""] <- sapply(zone,mundus_zone);
web_text <- collection_compendium$collection_name;
collection_compendium$collection_name <- sapply(web_text,mundify_web_text_boring);
web_text <- collection_compendium$stratcomments;
collection_compendium$stratcomments <- sapply(web_text,mundify_web_text_boring);

output <- list(collection_compendium,occurrences_compendium);
names(output) <- c("collection_compendium","occurrences_compendium");
return(output);
}

accersi_compendium_from_file_of_taxa <- function(taxonfilename,analysis_name,output_file=TRUE)	{
taxon_name <- read.table(file=taxonfilename,stringsAsFactors = TRUE)[,1]
taxon_list <- sapply(taxon_name,mundify_taxon_names)
compendium <- accersi_compendium_for_list_of_taxa(taxon_list,analysis_name,output_file)
return(compendium)
}

accersi_compendium_for_list_of_taxa <- function(taxon_list,analysis_name,output_file=TRUE)	{
# remove any funny symbols from taxon names
taxon_name <- taxon_list
taxon_list <- sapply(taxon_name,mundify_taxon_names)
ntaxa <- length(taxon_list)

compendium <- data.frame()
#entries <- array(0,ntaxa)
tx <- 1
while (tx <= ntaxa)	{
#for (tx in 1:ntaxa)	{
	print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	http1 <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon_list[tx],"&rank=genus,subgenus&private&show=attr,app",sep="")
	accio <- RCurl::getURL(http1)
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))

	## If the taxon is absent, then something weird will happen: kill it with fire.
	if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
		taxon_info <- compendium[1,]
		for(i in 1:length(taxon_info))	taxon_info[i] <-"?"
		tn <- match("taxon_name",colnames(compendium))
		taxon_info[tn] <- taxon_list[tx]
		no <- match("n_occs",colnames(compendium))
		taxon_info[no] <- 0
		}

	# if multiple taxa are returned, then winnow it down to one
	if (nrow(taxon_info)>1)	{
		tt <- match(taxon_list[tx],taxon_info$taxon_name)
		if (!is.na(tt))	{
			taxon_info <- taxon_info[tt,]
			} else	{
			taxon_info <- taxon_info[1,]
			tn <- match("taxon_name",colnames(compendium))
			taxon_info[tn] <- taxon_list[tx]
			}
		}

	## clear NAs
	if (is.na(taxon_info$difference))	taxon_info$difference <- ""
	taxon_info[(1:length(taxon_info))[is.na(taxon_info)]] <- ""
#	entries[tx] <- nrow(taxon_info)
#	if (entries[tx]>1)	print(paste(tx," has problems!!!"))
	# if the taxon is in the PaleoDB, then get the oldest & youngest occurrences
	oldest_rock <- youngest_rock <- "?"
	if (max(taxon_info$n_occs)>0)	{
		http2 <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon_list[tx],"&show=coll,strat,refattr&limit=all",sep = "")
		accio <- RCurl::getURL(http2)
		all_finds <- utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
		all_finds <- subset(all_finds,all_finds$identified_rank=="species")
		ttl_finds <- nrow(all_finds)
		if (ttl_finds>0)	{
			named_rock_unit <- all_finds$formation
			all_finds$formation <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
			named_rock_unit <- all_finds$member
			all_finds$member <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
			### get oldest rocks
			oldest <- (1:ttl_finds)[all_finds$max_ma==max(all_finds$max_ma)]
			oldest_seds <- c()
			for (rr in 1:length(oldest))	{
				ol <- oldest[rr]
				if (all_finds$member[ol]=="")	{
					zr <- as.character(all_finds$formation[ol])
					}	else	{
					zr <- paste(as.character(all_finds$formation[ol])," (",as.character(all_finds$member[ol]),")",sep="")
					}
				if (zr!="")	oldest_seds <- c(oldest_seds,zr)
				}
			if (length(oldest_seds)==0)	oldest_seds <- ""
			oldest_seds <- sort(unique(oldest_seds))
			oldest_rock <- c()
			if (length(oldest_seds)>0)	{
				rr <- 1
				while (rr<length(oldest_seds))	{
					oldest_rock <- paste(oldest_rock,oldest_seds[rr],"; ",sep="")
					rr <- rr+1
					}
				oldest_rock <- paste(oldest_rock,oldest_seds[rr],sep="")
				}	else	oldest_rock <- ""

			### get youngest rocks
			youngest <- (1:ttl_finds)[all_finds$min_ma==min(all_finds$min_ma)]
			youngest_seds <- c()
			for (rr in 1:length(youngest))	{
				yn <- youngest[rr]
				if (all_finds$member[yn]=="")	{
					zr <- as.character(all_finds$formation[yn])
					}	else	{
					zr <- paste(as.character(all_finds$formation[yn])," (",as.character(all_finds$member[yn]),")",sep="")
					}
				if (zr!="")	youngest_seds <- c(youngest_seds,zr)
				}
			if (length(youngest_seds)==0)	youngest_seds <- ""
			youngest_seds <- sort(unique(youngest_seds))
			youngest_rock <- c()
			if (length(youngest_seds)>0)	{
				rr <- 1
				while (rr<length(youngest_seds))	{
					youngest_rock <- paste(youngest_rock,youngest_seds[rr],"; ",sep="")
					rr <- rr+1
					}
				youngest_rock <- paste(youngest_rock,youngest_seds[rr],sep="")
				}	else	youngest_rock <- ""
			}
		} else	{
		a <- match("firstapp_max_ma",colnames(taxon_info))
		z <- match("late_interval",colnames(taxon_info))
		taxon_info[a:z] <- "?"
		}

	gee <- base::t(data.frame(c(oldest_rock,youngest_rock)))
	colnames(gee) <- c("oldest_rocks","youngest_rocks")
#	gee_science <- rbind(gee_science,c(oldest_rock,youngest_rock))
	compendium <- rbind(compendium, merge(taxon_info,gee))
	tx <- tx+1
	}

if (output_file)	{
	output_file_name <- paste(analysis_name,"_Compendium.xls",sep="")
	write.table(compendium,output_file_name,row.names = FALSE,col.names = TRUE,sep="\t")
	}
return(compendium)
}

#https://paleobiodb.org/data1.2/occs/list.txt?base_name=Bellimurina%20(Cyphomenoidea)%20wisgoriensis&show=refattr,classext,rem,entname&limit=all
accersi_paleodb_collection_numbers_for_taxa <- function(taxon,species_only=TRUE,directory="",save_files=TRUE,output_type=".csv") {
taxa <- paste(taxon, collapse = ",");
taxa <- gsub(" ","%20",taxa);
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&show=refattr,classext,rem,entname&limit=all",sep = "");
fetch <- RCurl::getURL(http);
all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
if (species_only)	{
	desired_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"));
	}	else	{
	desired_finds <- all_finds;
	}
#cleaned_names <- sapply(as.character(desired_finds$identified_name),mundify_occurrence_identifications)
cleaned_names <- sapply(as.character(desired_finds$identified_name),mundify_taxon_names)
desired_finds$identified_name <- cleaned_names;
cleaned_names <- sapply(as.character(desired_finds$accepted_name),mundify_taxon_names);
desired_finds$accepted_name <- cleaned_names;

if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="");
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(desired_finds,file=output,row.names = FALSE,fileEncoding = "UTF-8");
		}	else	{
		write.table(desired_finds,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}
return(desired_finds)
}

accersi_abundance_data <- function(taxa,onset="Cambrian",end="Holocene",realm="marine",save_collections=TRUE,directory="",save_files=TRUE,output_type=".txt") {
taxa <- paste(taxa, collapse = ",")
realm <- paste(realm, collapse=",")
if (save_collections)	{
	http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&envtype=",realm,"&interval=",onset,",",end,"&show=abund,pres,coll,coords,loc,paleoloc,strat,stratext,lith,lithext,env,geo,ref,ent,entname,crmod&limit=all",sep = "")
	}	else	{
	http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&envtype=",realm,"&interval=",onset,",",end,"&show=abund,coll&limit=all",sep = "")
	}
fetch <- RCurl::getURL(http)
all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
abundance_finds <- subset(all_finds,all_finds$abund_value!="")
abundance_finds <- subset(abundance_finds,abundance_finds$abund_unit!="category")
abundance_finds <- subset(abundance_finds,abundance_finds$abund_unit!="rank")
abundance_finds <- subset(abundance_finds,abundance_finds$abund_unit!="elements")
abundance_finds <- subset(abundance_finds,abundance_finds$abund_unit!="fragments")
cleaned_names <- sapply(as.character(abundance_finds$identified_name),mundify_taxon_names)
abundance_finds$identified_name <- cleaned_names

if (save_collections)	{
	coll_info <- abundance_finds
	coll_info$occurrence_no <- NULL
	coll_info$record_type <- NULL
	coll_info$reid_no <- NULL
	coll_info$flags <- NULL
	coll_info$identified_name <- NULL
	coll_info$identified_no <- NULL
	coll_info$difference <- NULL
	coll_info$accepted_name <- NULL
	coll_info$accepted_no <- NULL
	coll_info$abund_value <- NULL
	coll_info$abund_unit <- NULL
	coll_info$identified_rank <- NULL
	coll_info$accepted_rank <- NULL
	coll_no <- unique(sort(coll_info$collection_no))
	noccr <- nrow(coll_info)
	coll_unq <- (1:noccr)[match(coll_no,coll_info$collection_no)]
	coll_info <- coll_info[coll_unq,]
	named_rock_unit <- coll_info$formation
	xxx <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
	coll_info$formation <- xxx
	named_rock_unit <- coll_info$stratgroup
	xxx <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
	coll_info$stratgroup <- xxx
	named_rock_unit <- coll_info$member
	xxx <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
	coll_info$member <- xxx
	}
if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Abundances",output_type,sep="")
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(abundance_finds,file=output,row.names = FALSE,fileEncoding = "UTF-8");
		}	else	{
		write.table(abundance_finds,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}

if (save_collections)	{
	aa <- list(abundance_finds,coll_info)
	names(aa) <- c("Abundance_Data","Collection_Info")
	return(aa)
	} else	{
	return(species_finds)
	}
}

# get collections from a certain span of time that include a particular taxon or set of taxa
# get collections from a certain span of time that include a particular taxon or set of taxa
# modified 2020-02-17
# modified 2020-05-18
accersi_collection_data <- function(taxa="Life",onset="Archean",end="Cenozoic",basic_environments="terr,marine,unknown",paleogeography="gplates",standardize_members=F,directory="",save_files=T,species_only=F,output_type=".csv") {
if (length(taxa)>1)	taxa <- paste(taxa, collapse = ",");
taxa <- gsub(" ","%20",taxa);
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
#http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr,entname,lith,env,crmod",sep="");
collections <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
if (!is.null(collections$formation.1))	{
	collections$formation.1 <- NULL;
	collections$member.1 <- NULL;
	collections$stratgroup.1 <- NULL;
	}
if (!is.null(collections$lithology1.1))	{
	collections$lithdescript.1 <- collections$lithology1.1 <- collections$lithadj1.1 <- collections$lithification1.1 <- collections$minor_lithology1.1 <- collections$fossilsfrom1.1 <- NULL;
	collections$lithology2.1 <- collections$lithadj2.1 <- collections$lithification2.1 <- collections$minor_lithology2.1 <- collections$fossilsfrom2.1 <- NULL;
	}
#fetch <- RCurl::getURL(http);
#collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
#	http <- paste("https://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,",&show=full,etbasis,strat,lith,env,timebins,timecompare,ref,ent,entname,crmod",sep="")
if (species_only)	{
	http <- paste("https://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=full",sep="");
	species_finds <- read.csv(http, header = TRUE, stringsAsFactors=hell_no);
	species_finds <- subset(species_finds,species_finds$identified_rank %in% c("species","subspecies"));
	unique_collections <- unique(species_finds$collection_no);
	save_colls <- match(unique_collections,collections$collection_no);
	save_colls <- save_colls[!is.na(save_colls)];
	collections <- collections[save_colls,];
	}
ttl_coll <- nrow(collections);

collections <- put_pbdb_dataframes_into_proper_type(pbdb_data=collections);
geographic_orphans <- collections$collection_no[is.na(collections$paleolat)];
geographic_orphans <- sort(unique(c(geographic_orphans,collections$collection_no[collections$geoplate==0])));
#coll_no <- coll_no[1:5]
if (length(geographic_orphans)>0) {
	 coll_no <- geographic_orphans;
	 if (paleogeography=="gplates")	{
		print("Using Scotese's model to get paleogeographic data for collections gplates cannot.")
	 	model <- "scotese";
	 	} else	{
		print("Using the gplates model to get paleogeographic data for collections Scotese's model cannot.")
	 	model <- "gplates";
		}
	revised_paleogeography <- data.frame(base::t(pbapply::pbsapply(coll_no,paleogeographic_orphanarium,model)));
	for (cc in 1:ncol(revised_paleogeography))	revised_paleogeography[,cc] <- as.vector(unlist(revised_paleogeography[,cc]));
	revised_paleogeography <- put_pbdb_dataframes_into_proper_type(revised_paleogeography);
	leelas <- match(revised_paleogeography$collection_no,collections$collection_no);
	mutants <- match(colnames(revised_paleogeography),colnames(collections));
	collections[leelas,mutants] <- revised_paleogeography;
	}
# clean up rock unit names
#clean_groups <- mundify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- collections$stratgroup;
print("Cleaning Groups.")
clean_groups <- pbapply::pbsapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
collections$stratgroup <- clean_groups;
named_rock_unit <- collections$formation;
print("Cleaning Formations.")
clean_formations <- pbapply::pbsapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
collections$formation <- clean_formations;
print("Cleaning Members.")
named_rock_unit <- collections$member;
clean_members <- pbapply::pbsapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
collections$member <- clean_members;
collections$collection_subset <- expello_na_from_vector(collections$collection_subset,"")
# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(collections$ref_pubyr[vote_formation]) > max(collections$ref_pubyr[vote_member]))	{
					### elevate member to formation in appropriate collections
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}	else if (max(collections$ref_pubyr[vote_formation]) < max(collections$ref_pubyr[vote_member]))	{
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
					} else	{
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}
				}
			}
		}
	}

## clean zones of question marks, aff.s, etc.
collections_w_zones <- (1:ttl_coll)[collections$zone!=""];
cwz <- length(collections_w_zones);
collections$zone <- as.character(collections$zone);
#zones <- as.character(collections$zone[collections_w_zones]);
#zz <- c();
#for (z in 1:length(zones))	{
#	print(c(zones[z],mundus_zone(zone=zones[z])))
#	zz <- c(zz,mundus_zone(zones[z]));
#	}
print("Cleaning Zones")
zone <- as.character(collections$zone[collections_w_zones]);
collections$zone[collections_w_zones] <- pbapply::pbsapply(zone,mundus_zone);
#sort(unique(collections$zone[collections_w_zones]))

#collections <- expello_na_from_matrix(collections,replacement = "");
collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Collections",output_type,sep="")
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(collections,file=output,row.names = FALSE,fileEncoding = "UTF-8");
		}	else	{
		write.table(collections,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}

#collections <- put_pbdb_dataframes_into_proper_type(collections);
#for (cn in 1:ncol(collections))	{
#cn <- 0;
#while (cn<ncol(collections))	{
#	cn <- cn+1;
#	print(cn);
#	old_info <- collections[,cn];
#	cnm <- strsplit(x=colnames(collections)[cn],split="_")[[1]]
#	if(cnm[length(cnm)] %in% paleodb_numeric_fields && colnames(collections)[cn]!="collection_size")	{
#		old_info[old_info %in% missing_data_assignment] <- 0;
#		old_info[old_info=="coordinates not computable using this model"] <- 0;
#		collections[,cn] <- as.numeric(old_info)
#		} else	{
#		collections[,cn] <- as.character(old_info);
#		}
#	}

return(collections)
}

#colls_modified_after=2021-01-01; basic_environments=c("marine","unknown");
#taxa="Life";onset="Archean";end="Cenozoic";basic_environments="terr,marine,unknown";colls_modified_after="1900-01-01";paleogeography="gplates";standardize_members=F;directory="";save_files=T;species_only=F;output_type=".csv";
update_collection_data <- function(taxa="Life",onset="Archean",end="Cenozoic",basic_environments="terr,marine,unknown",colls_modified_after="1900-01-01",paleogeography="gplates",standardize_members=F,directory="",save_files=T,species_only=F,output_type=".csv")	{
taxa <- paste(taxa, collapse = ",");
taxa <- gsub(" ","%20",taxa);
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
#if (colls_modified_after=="1900-01-01" || colls_modified_after==as.Date("0000-00-00"))
if (colls_modified_after=="1900-01-01")
	colls_modified_after <- strsplit(as.character(Sys.time()-5*366*24*60*60),split=" ")[[1]][1];
http <- paste("https://www.paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&colls_modified_after=",colls_modified_after,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
collections_update <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
if (!is.null(collections_update$formation.1))	{
	collections_update$formation.1 <- NULL;
	collections_update$member.1 <- NULL;
	collections_update$stratgroup.1 <- NULL;
	}
if (!is.null(collections_update$lithology1.1))	{
	collections_update$lithdescript.1 <- collections_update$lithology1.1 <- collections_update$lithadj1.1 <- collections_update$lithification1.1 <- collections_update$minor_lithology1.1 <- collections_update$fossilsfrom1.1 <- NULL;
	collections_update$lithology2.1 <- collections_update$lithadj2.1 <- collections_update$lithification2.1 <- collections_update$minor_lithology2.1 <- collections_update$fossilsfrom2.1 <- NULL;
	}
collections_update <- put_pbdb_dataframes_into_proper_type(pbdb_data=collections_update);
geographic_orphans <- unique(sort(c(collections_update$collection_no[is.na(collections_update$paleolat)],collections_update$collection_no[collections_update$geoplate==0])));
#coll_no <- coll_no[1:5]
if (length(geographic_orphans)>0) {
	 coll_no <- geographic_orphans;
	 if (paleogeography=="gplates")	{
		print("Using Scotese's model to get paleogeographic data for collections gplates cannot.")
	 	model <- "scotese";
	 	} else	{
		print("Using the gplates model to get paleogeographic data for collections Scotese's model cannot.")
	 	model <- "gplates";
		}
	revised_paleogeography <- data.frame(base::t(pbapply::pbsapply(coll_no,paleogeographic_orphanarium,model)));
	for (cc in 1:ncol(revised_paleogeography))	revised_paleogeography[,cc] <- as.vector(unlist(revised_paleogeography[,cc]));
	revised_paleogeography <- put_pbdb_dataframes_into_proper_type(revised_paleogeography);
	leelas <- match(revised_paleogeography$collection_no,collections_update$collection_no);
	mutants <- match(colnames(revised_paleogeography),colnames(collections_update));
	collections_update[leelas,mutants] <- revised_paleogeography;
	}

ttl_coll <- nrow(collections_update);
collections_w_zones <- (1:ttl_coll)[collections_update$zone!=""];
cwz <- length(collections_w_zones);
collections_update$zone <- as.character(collections_update$zone);
print("Cleaning Zones")
zone <- as.character(collections_update$zone[collections_w_zones]);
collections_update$zone[collections_w_zones] <- pbapply::pbsapply(zone,mundus_zone);
#sort(unique(collections$zone[collections_w_zones]))

#collections <- expello_na_from_matrix(collections,replacement = "");
collections_update$late_interval[collections_update$late_interval==""] <- collections_update$early_interval[collections_update$late_interval==""];
if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Collections",output_type,sep="")
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(collections_update,file=output,row.names = FALSE,fileEncoding = "UTF-8");
		}	else	{
		write.table(collections_update,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}

return(collections_update);
}

# get collections from a certain span of time that include a particular taxon or set of taxa
# coll_id <- collection_nos[1]
accersi_data_for_one_collection <- function(coll_id) {
#http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?coll_id=",coll_id,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
fetch <- RCurl::getURL(http);
collection <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");

# clean up rock unit names
#clean_groups <- mundify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
collection <- expello_na_from_vector(collection);

if (collection$stratgroup!="")
	collection$stratgroup <- mundify_rock_unit_names(named_rock_unit=as.character(collection$stratgroup),dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
if (collection$formation!="")
	collection$formation <- mundify_rock_unit_names(named_rock_unit=as.character(collection$formation),dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
if (collection$member!="")
	collection$member <- mundify_rock_unit_names(named_rock_unit=as.character(collection$member),dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);

## clean zones of question marks, aff.s, etc.
if (collection$zone!="")
	collection$zone <- mundus_zone(zone=as.character(collection$zone));

if (!is.null(collection$formation.1))	{
	collection$formation.1 <- NULL;
	collection$member.1 <- NULL;
	collection$stratgroup.1 <- NULL;
	}

if (!is.null(collection$lithology1.1))	{
	collection$lithdescript.1 <- collection$lithology1.1 <- collection$lithadj1.1 <- collection$lithification1.1 <- collection$minor_lithology1.1 <- collection$fossilsfrom1.1 <- NULL;
	collection$lithology2.1 <- collection$lithadj2.1 <- collection$lithification2.1 <- collection$minor_lithology2.1 <- collection$fossilsfrom2.1 <- NULL;
	}

collection <- put_pbdb_dataframes_into_proper_type(pbdb_data=collection);
if (is.na(collection$paleolat))	{
	scotesed <- paleogeographic_orphanarium(coll_no=coll_id,model="scotese");
	collection[,match(colnames(scotesed),colnames(collection))] <- scotesed;
	}

return(collection)
}

# get collections data for one locality
# modified 2020-05-02
# modified 2020-05-19
accersi_single_locality_info <- function(collection_no,paleogeography="gplates")	{
collection_no <- as.character(collection_no);
stupid_exponents <- strsplit(collection_no,"")[[1]];
if (sum(stupid_exponents=="+")==1 && sum(stupid_exponents=="e")==1)	{
	doh <- match("+",stupid_exponents);
	expn <- as.numeric(paste(stupid_exponents[(doh+1):length(stupid_exponents)],collapse=""));
	stupid_exponents2 <-paste(stupid_exponents[1:(doh-2)],collapse="");
	for (ee in 1:expn)	{
		stupid_exponents2 <- c(stupid_exponents2,"0");
		}
	collection_no <- paste(stupid_exponents2,collapse="");
	}

#httpC <- paste("https://paleobiodb.org/data1.2/colls/list.csv?id=",collection_no,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
#http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
httpC <- paste("https://paleobiodb.org/data1.2/colls/list.csv?id=",collection_no,"&pgm=",paleogeography,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,secref,refattr,ent,entname,crmod",sep="");
coll_info <- read.csv(httpC, header = TRUE, stringsAsFactors=hell_no);
options(warn=0);
coll_info <- put_pbdb_dataframes_into_proper_type(coll_info);
coll_info <- clear_na_from_vector(coll_info,"")
options(warn=1);
return(coll_info);
}

#collection_nos <-
accersi_collection_data_for_list_of_collection_nos <- function(collection_nos)	{
coll_id <- collection_nos;
pbdb_site_data <- data.frame(simplify2array(base::t(pbapply::pbsapply(coll_id,accersi_data_for_one_collection))));
for (nc in 1:ncol(pbdb_site_data))	pbdb_site_data[,nc] <- as.vector(unlist(pbdb_site_data[,nc]));
if (!is.null(pbdb_site_data$formation.1))
	pbdb_site_data$formation.1 <- pbdb_site_data$stratgroup.1 <- pbdb_site_data$member.1 <- NULL;
if (!is.null(pbdb_site_data$lithdescript.1))
	pbdb_site_data$lithdescript.1 <- pbdb_site_data$lithology1.1 <- pbdb_site_data$lithadj1.1 <- pbdb_site_data$lithification1.1 <- pbdb_site_data$minor_lithology1.1 <- pbdb_site_data$fossilsfrom1.1 <- pbdb_site_data$lithology2.1 <- pbdb_site_data$lithadj2.1 <- pbdb_site_data$lithification2.1 <- pbdb_site_data$minor_lithology2.1 <- pbdb_site_data$fossilsfrom2.1 <- NULL;
pbdb_site_data$created <- as.Date.numeric(as.numeric(pbdb_site_data$created),origin = "1970-01-01");
pbdb_site_data$modified <- as.Date.numeric(as.numeric(pbdb_site_data$modified),origin = "1970-01-01");
pbdb_site_data <- put_pbdb_dataframes_into_proper_type(pbdb_data=pbdb_site_data);
return(pbdb_site_data);
}

# get paleogeographic data from the other model (scotese if gplates, gplates if scotese)
paleogeographic_orphanarium <- function(coll_no,model="scotese")	{
http <- paste("https://paleobiodb.org/data1.2/colls/single.csv?id=",coll_no,"&pgm=",model,"&show=loc,paleoloc",sep="");
site_reboot <- read.csv(http, header = TRUE, stringsAsFactors=hell_no);
keepers <- c("collection_no","paleomodel","paleolng","paleolat","geoplate")
site_reboot <- site_reboot[,match(keepers,colnames(site_reboot))]
#colnames(site_reboot)
return(site_reboot);
}

#rock_unit <- "Ain El Guettar";
#rock_unit <- "Oum Ed Dhiab";
#rock_unit <- "Cape Phillips";
#rock_unit <- "Alum";
#rock_unit <- "Star Peak";
#rock_unit <- "A Congeria";
#rock_unit <- "Antimonio";
#rock_unit <- "Agreda";
#rock_unit <- "facies à Congeria";
#rock_unit <- "Honda";
#rock_unit <- "Upper Nawata";
accersi_rock_unit_data <- function(rock_unit,ma_lb=4000,ma_ub=0)	{
#rock_unit: name of rock unit
#ma_lb: oldest rocks to consider. (Here for homonyms)
#ma_ub: oldest rocks to consider. (Here for homonyms)
rock_unit_url <- gsub(" ","%20",rock_unit);
http <- paste("https://paleobiodb.org/data1.2/strata/list.csv?name=",rock_unit_url,"&show=gplates,splates",sep="");
rock_info <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
if (mundify_rock_unit_names(rock_unit,delete_rock_type = T)!=rock_unit)	{
	rock_unit_url <- gsub(" ","%20",mundify_rock_unit_names(rock_unit,delete_rock_type = T));
	http <- paste("https://paleobiodb.org/data1.2/strata/list.csv?name=",rock_unit_url,"&show=gplates,splates",sep="");
	rock_info <- unique(rbind(rock_info,read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8")));
	}
if (mundify_rock_unit_names(rock_unit,delete_informal=T)!=rock_unit)	{
	rock_unit_url <- gsub(" ","%20",mundify_rock_unit_names(rock_unit,delete_informal = T));
	http <- paste("https://paleobiodb.org/data1.2/strata/list.csv?name=",rock_unit_url,"&show=gplates,splates",sep="");
	rock_info <- unique(rbind(rock_info,read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8")));
	}
if (mundify_rock_unit_names(rock_unit,delete_rock_type=T,delete_informal=T)!=rock_unit)	{
	rock_unit_url <- gsub(" ","%20",mundify_rock_unit_names(rock_unit,delete_rock_type=T,delete_informal = T));
	http <- paste("https://paleobiodb.org/data1.2/strata/list.csv?name=",rock_unit_url,"&show=gplates,splates",sep="");
	rock_info <- unique(rbind(rock_info,read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8")));
	}
if (gsub(" ","-",mundify_rock_unit_names(rock_unit,delete_rock_type=T,delete_informal=T))!=mundify_rock_unit_names(rock_unit,delete_rock_type=T,delete_informal=T))	{
	rock_unit_url <- gsub(" ","-",mundify_rock_unit_names(rock_unit,delete_rock_type=T,delete_informal=T));
	http <- paste("https://paleobiodb.org/data1.2/strata/list.csv?name=",rock_unit_url,"&show=gplates,splates",sep="");
	rock_info <- unique(rbind(rock_info,read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8")));
	}
rock_info <- rock_info[rock_info$record_type!="THIS REQUEST RETURNED NO RECORDS",];
# if there are lumped rock units, then separate them.
# delete redundancies. Assume lone entries are the better ones
if (nrow(rock_info)>0)	{
	rock_info <- clear_na_from_matrix(rock_info,"");
	rock_info$group_alt <- rock_info$member_alt <- rock_info$formation_alt <- rep("",nrow(rock_info));

	rock_name <- rock_info$formation[rock_info$formation!=""];
	if (length(rock_name)>0)	{
		dr <- t(sapply(rock_name,separate_two_rock_names))
		rock_info$formation_alt[rock_info$formation!=""] <- unlist(dr[,2]);
		rock_info$formation[rock_info$formation!=""] <- unlist(dr[,1]);
		added_rock <- rock_info[rock_info$formation_alt!="",];
		added_rock$formation <- added_rock$formation_alt;
		rock_info <- rbind(rock_info,added_rock);
		}

	rock_name <- rock_info$member[rock_info$member!=""];
	if (length(rock_name)>0)	{
		dr <- t(sapply(rock_name,separate_two_rock_names));
		rock_info$member_alt[rock_info$member!=""] <- unlist(dr[,2]);
		rock_info$member[rock_info$member!=""] <- unlist(dr[,1]);
		added_rock <- rock_info[rock_info$member_alt!="",];
		added_rock$member <- added_rock$member_alt;
		rock_info <- rbind(rock_info,added_rock);
		}

	full_names <- rock_info$formation;
	full_names[rock_info$member!=""] <- paste(full_names[rock_info$member!=""]," (",rock_info$member[rock_info$member!=""],")",sep="");
	unique_full_names <- unique(full_names);
	if (length(unique_full_names)>0)
		rock_info <- rock_info[unique(match(full_names,unique_full_names)),];
	rock_info$member_alt <- rock_info$formation_alt <- rep("",nrow(rock_info));

	rock_name <- rock_info$group[rock_info$group!=""];
	if (length(rock_name)>0)	{
		dr <- t(sapply(rock_name,separate_two_rock_names));
		rock_info$group_alt[rock_info$group!=""] <- unlist(dr[,2]);
		rock_info$group[rock_info$group!=""] <- unlist(dr[,1]);
		added_rock <- rock_info[rock_info$group!="",];
		added_rock$group <- added_rock$group;
		rock_info <- rbind(rock_info,added_rock);
		}
	rock_info$group_alt <- rock_info$member_alt <- rock_info$formation_alt <- NULL;
	}

if (nrow(rock_info)>1)	{
	rock_info <- subset(rock_info,rock_info$min_ma<=ma_lb);
	rock_info <- subset(rock_info,rock_info$max_ma>=ma_ub);
	for (rr in 1:nrow(rock_info))	{
		plates <- as.numeric(strsplit(as.character(rock_info$gplate_no[rr]),",")[[1]]);
		plates <- plates[plates!=0];
		if (length(plates)==0)	plates <- 0;
		rock_info$gplate_no[rr] <- as.numeric(plates[1]);
		plates <- as.numeric(strsplit(as.character(rock_info$splate_no[rr]),",")[[1]]);
		plates <- plates[plates!=0];
		if (length(plates)==0)	plates <- 0;
		rock_info$splate_no[rr] <- plates[1];
		}
	rock_info <- put_pbdb_dataframes_into_proper_type(rock_info);

	named_rock_unit <- rock_info$formation[rock_info$formation!=""];
	if (length(named_rock_unit)>0)	rock_info$formation[rock_info$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names);
	named_rock_unit <- rock_info$member[rock_info$member!=""];
	if (length(named_rock_unit)>0)	rock_info$member[rock_info$member!=""] <- sapply(named_rock_unit,mundify_rock_unit_names);
	named_rock_unit <- rock_info$group[rock_info$group!=""];
	if (length(named_rock_unit)>0)	rock_info$group[rock_info$group!=""] <- sapply(named_rock_unit,mundify_rock_unit_names);
	strata_rank_cells <- which(rock_info==rock_unit,arr.ind = T);
	strata_rank_cols <- sort(unique(strata_rank_cells[,2]));
	strata_rank <- colnames(rock_info)[sort(unique(which(rock_info==rock_unit,arr.ind = T)[,2]))];
	if (length(strata_rank)==1 && nrow(rock_info)>1)	{
		if (strata_rank=="formation")	{
#			this_unit_info <- subset(rock_info,rock_info$formation==rock_unit);
			if (nrow(rock_info)>1)	{
				this_unit <- rock_info[1,];
				this_unit$member <- ""
				this_unit$max_ma <- max(rock_info$max_ma);
				this_unit$min_ma <- min(rock_info$min_ma);
				this_unit$n_colls <- sum(rock_info$n_colls);
				this_unit$n_occs <- sum(rock_info$n_occs);
				lithology <- c();
				for (rr in 1:nrow(rock_info))	{
					lithologies <- strsplit(gsub("/",",",rock_info$lithology[rr]),",")[[1]]
					lithology <- unique(c(lithology,lithologies));
					}
				lithology <- lithology[lithology!="not reported"];
				if (length(lithology)>1)	this_unit$lithology <- paste(lithology,collapse=",");
				}
			} else if (strata_rank=="member")	{
			### rewrite this to do it by each formation + member combo
			formations <- unique(rock_info$formation);
			this_unit <- rock_info[1:length(formations),];
			this_unit$formation <- formations;
			this_unit$member <- rep(rock_unit,length(formations));
			for (ff in 1:length(formations))	{
				this_unit_info <- subset(rock_info,rock_info$formation==formations[ff]);
				group <- this_unit_info$group[this_unit_info$group!=""]
				if (length(group)==0)	{
					this_unit$group[ff] <- "";
					} else	{
					this_unit$group[ff] <- group[1];
					}
				this_unit$max_ma[ff] <- max(this_unit_info$max_ma);
				this_unit$min_ma[ff] <- min(this_unit_info$min_ma);
				this_unit$n_colls[ff] <- sum(this_unit_info$n_colls);
				this_unit$n_occs[ff] <- sum(this_unit_info$n_occs);
				this_unit$gplate_no[ff] <- this_unit_info$gplate_no[1];
				this_unit$splate_no[ff] <- this_unit_info$splate_no[1];
				lithology <- c();
				for (rr in 1:nrow(this_unit_info))	{
					lithologies <- strsplit(gsub("/",",",this_unit_info$lithology[rr]),",")[[1]]
					lithology <- unique(c(lithology,lithologies));
					}
				lithology <- lithology[lithology!="not reported"];
				if (length(lithology)>1)	this_unit$lithology[ff] <- paste(lithology,collapse=",");
				}
			} else	{
			this_unit <- rock_info[1,];
			this_unit$formation <- this_unit$member <- "";
			this_unit$max_ma <- max(rock_info$max_ma);
			this_unit$min_ma <- min(rock_info$min_ma);
			this_unit$n_colls <- sum(rock_info$n_colls);
			this_unit$n_occs <- sum(rock_info$n_occs);
			lithology <- c();
			for (rr in 1:nrow(rock_info))	{
				lithologies <- strsplit(rock_info$lithology[rr],",")[[1]]
				lithology <- unique(c(lithology,lithologies));
				}
			lithology <- lithology[lithology!="not reported"];
			if (length(lithology)>1)	this_unit$lithology <- paste(lithology,collapse=",");
			}
		} else if (length(strata_rank)>1)	{
		this_unit_info <- rock_info[unique(strata_rank_cells[,1]),];
		this_unit <- this_unit_info[1,];
		this_unit$max_ma <- max(this_unit_info$max_ma);
		this_unit$min_ma <- min(this_unit_info$min_ma);
		this_unit$n_colls <- sum(this_unit_info$n_colls);
		this_unit$n_occs <- sum(this_unit_info$n_occs);
		lithology <- c();
		for (rr in 1:nrow(this_unit_info))	{
			lithologies <- strsplit(this_unit_info$lithology[rr],",")[[1]]
			lithology <- unique(c(lithology,lithologies));
			}
		lithology <- lithology[lithology!="not reported"];
		if (length(lithology)>1)	this_unit$lithology <- paste(lithology,collapse=",");
		} else	{
		this_unit <- rock_info
		}
	} else	{
	this_unit <- rock_info;
	}
return(this_unit);
}

# get information about rock units within a certain interval of time.
accersi_rock_unit_data_for_taxa <- function(taxa,onset="Cambrian",end="Holocene",standardize_members=TRUE,directory="",save_files=TRUE,output_type=".csv") {
taxa <- paste(taxa, collapse = ",");
http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr",sep="")
fetch <- RCurl::getURL(http);
collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
ttl_coll <- nrow(collections);

# clean up rock unit names
#clean_groups <- mundify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- collections$stratgroup
clean_groups <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
collections$stratgroup <- clean_groups
named_rock_unit <- collections$formation;
clean_formations <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
collections$formation <- clean_formations
named_rock_unit <- collections$member
clean_members <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
collections$member <- clean_members

# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(collections$ref_pubyr[vote_formation]) > max(collections$ref_pubyr[vote_member]))	{
					### elevate member to formation in appropriate collections
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}	else if (max(collections$ref_pubyr[vote_formation]) < max(collections$ref_pubyr[vote_member]))	{
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
					} else	{
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}
				}
			}
		}
	}
keep <- match(c("formation","member","stratgroup","zone","early_interval","late_interval","max_ma","min_ma"),colnames(collections))
remove <- (1:ncol(collections))[!(1:ncol(collections)) %in% keep]
rock_info <- collections;
rock_info <- rock_info[,-remove];
rock_info <- rock_info[,order(keep)];
rock_info <- unique(rock_info);

if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Rock_Units",output_type,sep="")
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(rock_info,file=output,row.names = FALSE,fileEncoding = "UTF-8");
		}	else	{
		write.table(rock_info,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}
return(rock_info)
}

# summarize rock unit data from a set of paleodb collections
accersi_rock_unit_data_from_paleodb_collections <- function(paleodb_collections,standardize_members=TRUE) {
ttl_coll <- nrow(paleodb_collections);

# clean up rock unit names
#clean_groups <- mundify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- paleodb_collections$stratgroup;
clean_groups <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
paleodb_collections$stratgroup <- clean_groups;
named_rock_unit <- paleodb_collections$formation;
clean_formations <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
paleodb_collections$formation <- clean_formations;
named_rock_unit <- paleodb_collections$member;
clean_members <- sapply(as.character(named_rock_unit),mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
paleodb_collections$member <- clean_members;

# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(paleodb_collections$ref_pubyr[vote_formation]) > max(paleodb_collections$ref_pubyr[vote_member]))	{;
					### elevate member to formation in appropriate collections
					paleodb_collections$formation[vote_member] <- member_or_formation[c];
					paleodb_collections$member[vote_member] <- "";
					}	else if (max(paleodb_collections$ref_pubyr[vote_formation]) < max(paleodb_collections$ref_pubyr[vote_member]))	{;
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(paleodb_collections$ref_pubyr[vote_member]),paleodb_collections$ref_pubyr[vote_member])];
					paleodb_collections$formation[vote_formation] <- paleodb_collections$formation[latest_opinion];
					paleodb_collections$member[vote_formation] <- member_or_formation[c];
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(paleodb_collections$ref_pubyr[vote_member]),paleodb_collections$ref_pubyr[vote_member])];
					paleodb_collections$formation[vote_formation] <- paleodb_collections$formation[latest_opinion];
					paleodb_collections$member[vote_formation] <- member_or_formation[c];
					} else	{
					paleodb_collections$formation[vote_member] <- member_or_formation[c];
					paleodb_collections$member[vote_member] <- "";
					}
				}
			}
		}
	}

keep <- match(c("formation","member","stratgroup","zone","early_interval","late_interval","max_ma","min_ma"),colnames(paleodb_collections));
if (!is.null(paleodb_collections$ma_lb))	{
	keep <- c(keep,match(c("ma_lb","ma_ub"),colnames(paleodb_collections)));
	}
if (!is.null(paleodb_collections$interval_lb))	{
	keep <- c(keep,match(c("interval_lb","interval_ub"),colnames(paleodb_collections)));
	}
if (!is.null(paleodb_collections$direct_ma))	{
	keep <- c(keep,match(c("direct_ma","direct_ma_error"),colnames(paleodb_collections)));
	}
remove <- (1:ncol(paleodb_collections))[!(1:ncol(paleodb_collections)) %in% keep];
rock_info <- paleodb_collections;
rock_info <- rock_info[,-remove];
rock_info <- rock_info[,order(keep)];
#ncolls <- nrow(rock_info);
#rrr <- (1:ncolls)[rock_info$formation %in% "Shuayb"]
#rock_info[rrr,]
rock_info <- unique(rock_info);

return(rock_info)
}

# routine to organize Paleobiologyy Database rock units to facilitate matching to external rock unit database
accersi_rock_units_by_time_and_general_environment <- function(first_interval,last_interval,environment,delete_rock_type=TRUE)	{
httpR <- paste("https://www.paleobiodb.org/data1.2/occs/strata.csv?base_name=Metazoa&interval=",first_interval,",",last_interval,"&envtype=",tolower(environment),"&show=coords,gplates,splates",sep="")
#httpR <- paste("https://www.paleobiodb.org/data1.2/occs/strata.csv?base_name=Metazoa&interval="Ordovician",Llandovery&envtype=marine&show=coords,gplates,splates"
accio <- RCurl::getURL(httpR)
rock_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
print("data downloaded")

formation_identified <- rock_info$formation <- as.character(rock_info$formation)
id_forms <- (1:nrow(rock_info))[rock_info$formation!=""]
named_rock_unit <- rock_info$formation[id_forms]
taxi <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
rock_info$formation[id_forms] <- taxi
print("formations cleaned")

member_identified <- rock_info$member <- as.character(rock_info$member)
id_membs <- (1:nrow(rock_info))[rock_info$member!=""]
named_rock_unit <- rock_info$member[id_membs]
taxi <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
rock_info$member[id_membs] <- taxi
print("members cleaned")

rock_info$group <- as.character(rock_info$group)
id_groups <- (1:nrow(rock_info))[rock_info$group!=""]
named_rock_unit <- rock_info$group[id_groups]
taxi <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
rock_info$group[id_groups] <- taxi
print("groups cleaned")

id_both <- id_forms[id_forms %in% id_membs]
rock_units <- rock_info$formation
rock_units[id_both] <- paste(rock_info$formation[id_both]," (",rock_info$member[id_both],")",sep="")
members_only <- id_membs[!id_membs %in% id_forms]
rock_units[members_only] <- rock_info$member[members_only]
group_only <- id_groups[!id_groups %in% sort(unique(c(id_forms,id_membs)))]
rock_units[group_only] <- paste(rock_info$group[group_only],"Unnamed")
rock_info <- cbind(rock_info,rock_units,formation_identified,member_identified)
unique_rocks <- unique(rock_units)
while (!is.na(match("",unique_rocks)))	{
	xxx <- match("",unique_rocks)
	unique_rocks[xxx:(length(unique_rocks)-1)] <- unique_rocks[(xxx+1):length(unique_rocks)]
	}
keepers <- match(unique_rocks,rock_units)
print("keepers found")
# reduce rock_info only to those records worth keeping
rock_info <- rock_info[keepers,]
# sort rock_info by full names of rock units
rock_info <- rock_info[order(rock_info$rock_units),]
# A kluge for the Kazahikstan bug
rock_info$cc_list[(1:nrow(rock_info))[rock_info$cc_list %in% ""]] <- NA
rock_info$cc_list <- expello_na_from_vector(as.character(rock_info$cc_list),"KZ")
print("KZ added");

if (delete_rock_type)	{
	id_forms <- (1:nrow(rock_info))[rock_info$formation!=""]
	cleaned_formations <- rock_info$formation
	named_rock_unit <- cleaned_formations[id_forms]
	taxi <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	cleaned_formations[id_forms] <- taxi
	print("formation rocks deleted")

	id_membs <- (1:nrow(rock_info))[rock_info$member!=""]
	cleaned_members <- rock_info$member
	named_rock_unit <- cleaned_members[id_membs]
	taxi <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	cleaned_members[id_membs] <- taxi
	print("member rocks deleted")

	id_both <- id_forms[id_forms %in% id_membs]
	rock_units_names_only <- cleaned_formations
	rock_units_names_only[id_both] <- paste(cleaned_formations[id_both]," (",cleaned_members[id_both],")",sep="")
	members_only <- id_membs[!id_membs %in% id_forms]
	rock_units_names_only[members_only] <- cleaned_members[members_only]
	rock_info  <- cbind(rock_info,rock_units_names_only)
	print("rock-free full names constructed")
	}

if (first_interval!=last_interval)	{
	output_file <- paste(environment,"_Rock_Units_",first_interval,"_to_",last_interval,".txt",sep="")
	} else	{
	output_file <- paste(environment,"_Rock_Units_",first_interval,".txt",sep="")
	}
rock_info$record_type <- NULL;
write.table(rock_info,output_file,col.names = TRUE,row.names = FALSE,sep="\t");
return(rock_info);
}

accersi_unique_rock_units_for_list_of_taxa <- function(taxon_list,analysis_name,output_file=TRUE)	{
# remove any funny symbols from taxon names
taxon_name <- taxon_list
taxon_list <- sapply(taxon_name,mundify_taxon_names)
ntaxa <- length(taxon_name)

rockpendium <- data.frame(row.names=FALSE,stringsAsFactors = FALSE)	# occurrences
rock_data <- data.frame(row.names=FALSE,stringsAsFactors = FALSE)
#while (tx <= ntaxa)	{
for (tx in 1:ntaxa)	{
	#print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	taxon <- gsub(" ","%20",taxon_list[tx])
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=attr,app",sep="")
	## If the taxon is absent, then something weird will happen: kill it with fire.
	accio <- RCurl::getURL(httpT)
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
	if (!(ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS")))	{
		http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxon,"&show=full",sep="")
		accio <- RCurl::getURL(http)
		taxon_finds <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
		if (nrow(taxon_finds)==1 && taxon_finds$collection_no=="THIS REQUEST RETURNED NO RECORDS")	{
			n_occs <- 0
			} else	{
			n_occs <- nrow(taxon_finds)
			taxon_finds <- expello_na_from_matrix(taxon_finds, "")
			named_rock_unit <- taxon_finds$formation
			taxon_finds$formation <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
			named_rock_unit <- taxon_finds$member
			taxon_finds$member <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
			rocks <- unique(cbind(as.character(taxon_finds$formation),as.character(taxon_finds$member)))
			rocks_full_info <- cbind(as.character(taxon_finds$formation),as.character(taxon_finds$member),as.character(taxon_finds$early_interval),as.character(taxon_finds$late_interval),as.character(taxon_finds$zone))
			colnames(rocks_full_info) <- c("formation","member","early_interval","late_interval","zone")
			rocks_full_info <- data.frame(rocks_full_info)
			rock_data <- rbind(rock_data,rocks_full_info)
			full_names_unique <- paste(rocks[,1],rocks[,2])
			full_names_all <- paste(taxon_finds$formation,taxon_finds$member)
			rock_finds <- coll_nos <- c()
			stages <- zonations <- c()
			for (i in 1:length(full_names_unique))	{
				rock_finds <- c(rock_finds,sum(full_names_all %in% full_names_unique[i]))
				colls <- taxon_finds$collection_no[(1:n_occs)[full_names_all %in% full_names_unique[i]]]
				xxx <- rocks_full_info[match(colls,taxon_finds$collection_no),]
				intervals <- c(as.character(xxx$early_interval),as.character(xxx$late_interval))
				intervals <- unique(intervals[intervals!=""])
				stage <- intervals[1]
				j <- 1
				while (j < length(intervals))	{
					j <- j+1
					stage <- paste(stage,", ",intervals[j],sep="")
					}
				zones <- sort(unique(as.character(xxx$zone)),decreasing = TRUE)
				zonation <- zones[1]
				j <- 1
				while (j < length(zones))	{
					j <- j+1
					if (zones[j]!="")
					zonation <- paste(zonation,", ",zones[j],sep="")
					}
				cn <- c()
				for (j in 1:length(colls))	{
					if (j < length(colls))	{
						cn <- paste(cn,colls[j],", ",sep="")
						} else	{
						cn <- paste(cn,colls[j],sep="")
						}
					}
				stages <- c(stages,stage)
				zonations <- c(zonations,zonation)
				coll_nos <- c(coll_nos,cn)
				}
			rocks <- cbind(rocks,rock_finds,stages,zonations,coll_nos)
			if (length(full_names_unique)>1)	{
				rocks <- rocks[order(rocks[,2]),]
				rocks <- rocks[order(rocks[,1]),]
				}
			if (rocks[1,1]=="" && rocks[1,2]=="")	rocks[1,1] <- "unreported"
			rocks <- cbind(rep(taxon_list[tx],length(full_names_unique)),rocks)
			rockpendium <- rbind(rockpendium,rocks)
			}
#		ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],taxon_info$n_occs[1]))
		}
#	taxon_list[tx]
#	tx <- tx+1
	}

colnames(rockpendium) <- c("taxon","formation","member","finds","stages","zones","collections")
rownames(rockpendium) <- NULL
rock_data <- unique(rock_data)
whole_name <- vector(length=nrow(rock_data))
form_only <- (1:nrow(rock_data))[rock_data$member==""]
form_memb <- (1:nrow(rock_data))[rock_data$member!=""]
whole_name[form_only] <- as.character(rock_data$formation[form_only])
whole_name[form_memb] <- paste(as.character(rock_data$formation[form_memb]),as.character(rock_data$member[form_memb]))
rock_data <- cbind(whole_name,rock_data)
rock_data <- rock_data[order(whole_name),]
rock_data <- subset(rock_data,rock_data$whole_name!="")
whole_name <- unique(rock_data$whole_name)
rock_data$early_interval <- as.character(rock_data$early_interval)
rock_data$late_interval <- as.character(rock_data$late_interval)
rock_data$zone <- as.character(rock_data$zone)
rock_data$whole_name <- as.character(rock_data$whole_name)
for (i in 1:length(whole_name))	{
	if (sum(rock_data$whole_name==whole_name[i])>1)	{
		mc <- (1:nrow(rock_data))[rock_data$whole_name %in% whole_name[i]]
#		rd <- subset(rock_data,rock_data$whole_name==as.character(whole_name[i]))
		rd <- rock_data[mc,]
#		rd$zone <- as.character(rd$zone)
		nn <- sum(rock_data$whole_name==whole_name[i])
		if (length(unique(rd$early_interval))>1)	{
			el <- sort(unique(rd$early_interval),decreasing = FALSE)
			j <- 1
			while (el[j]=="")	el <- el[2:length(el)]
			ne <- el[1]
			while (j<length(el))	{
				j <- j+1
				ne <- paste(ne,", ",as.character(el[j]),sep="")
				}
			for (j in 1:nn)	rd$early_interval[j] <- ne
			}
		if (length(unique(rd$late_interval))>1)	{
			ll <- sort(unique(rd$late_interval),decreasing = FALSE)
			j <- 1
			while (ll[j]=="")	ll <- ll[2:length(ll)]
			nl <- ll[1]
			while (j<length(ll))	{
				j <- j+1
				nl <- paste(nl,", ",as.character(ll[j]),sep="")
				}
			for (j in 1:nn)	rd$late_interval[j] <- nl
			}
		if (length(unique(rd$zone))>1)	{
			zz <- sort(unique(rd$zone),decreasing = FALSE)
			j <- 1
			while (zz[j]=="")	zz <- zz[2:length(zz)]
			nz <- zz[1]
			while (j<length(zz))	{
				j <- j+1
				nz <- paste(nz,", ",as.character(zz[j]),sep="")
				}
			for (j in 1:nn)	rd$zone[j] <- nz
			}
		rock_data[mc,] <- rd
		}
	}
rock_data <- unique(rock_data)
#unique_rocks <- rockpendium
#unique_rocks$taxon <- NULL
#unique_rocks$finds <- NULL
#unique_rocks$collections <- NULL
#rownames(unique_rocks) <- NULL
#unique_rocks <- unique(unique_rocks)
#whole_name <- c()
#for (u in 1:nrow(unique_rocks))	{
#	if (unique_rocks$member[u]=="")	{
#		whole_name <- c(whole_name,as.character(unique_rocks$formation[u]))
#		}	else	{
#		whole_name <- c(whole_name,paste(unique_rocks$formation[u]," (",unique_rocks$member[u],")",sep=""))
#		}
#	}
#unique_rocks <- cbind(unique_rocks,whole_name)
#unique_rocks <- unique_rocks[order(unique_rocks$whole_name),]

if (output_file)	{
	output_file1 <- paste(analysis_name,"_Rock_Units_per_Taxon.xls",sep="")
	write.table(rockpendium,output_file1,col.names=TRUE,row.names=FALSE,sep="\t")
	output_file2 <- paste(analysis_name,"_Unique_Rock_Units.xls",sep="")
	write.table(rock_data,output_file2,col.names=TRUE,row.names=FALSE,sep="\t")
	}

pendia <- list(rockpendium,rock_data)
names(pendia) <- c("Taxon_Info","Rocks_Units")
return(pendia)
}

accersi_unique_rock_units_for_higher_taxon <- function(higher_taxon,onset="Cambrian",end="Holocene",standardize_members=TRUE,directory="",save_files=TRUE,species_only=FALSE,output_type=".txt")	{
# remove any funny symbols from taxon names
#taxon_name <- taxon_list
#taxon_list <- sapply(taxon_name,mundify_taxon_names)
#ntaxa <- length(taxon_name)
http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",higher_taxon,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr",sep="")
fetch <- RCurl::getURL(http)
collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
#	http <- paste("https://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,",&show=full,etbasis,strat,lith,env,timebins,timecompare,ref,ent,entname,crmod",sep="")
if (species_only)	{
	http <- paste("https://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",higher_taxon,"&interval=",onset,",",end,",&show=full",sep="")
	fetch <- RCurl::getURL(http)
	species_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
	species_finds <- subset(species_finds,species_finds$identified_rank=="species")
	unique_collections <- unique(species_finds$collection_no)
	save_colls <- match(unique_collections,collections$collection_no)
	collections <- collections[save_colls,]
	}
ttl_coll <- nrow(collections)
named_rock_unit <- as.character(collections$formation)
formation_cleaned <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
named_rock_unit <- as.character(collections$member)
member_cleaned <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);

with_members <- (1:ttl_coll)[member_cleaned!=""]
without_members <- (1:ttl_coll)[member_cleaned==""]
with_formations <- (1:ttl_coll)[formation_cleaned!=""]
without_formations <- (1:ttl_coll)[formation_cleaned==""]
rock_units <- array("",dim=ttl_coll)
xxx <- with_formations[with_formations %in% without_members]
rock_units[xxx] <- formation_cleaned[xxx]
xxx <- with_formations[with_formations %in% with_members]
rock_units[xxx] <- paste(formation_cleaned[xxx]," (",member_cleaned[xxx],")",sep="")
xxx <- with_members[with_members %in% without_formations]
if (length(xxx)>0)
	rock_units[xxx] <- member_cleaned[xxx]
xxx <- without_formations[without_formations %in% without_members]
yyy <- xxx[collections$stratgroup[xxx]!=""]
for (c in 1:length(yyy))	{
	if (collections$late_interval[yyy[c]]=="")	{
		rock_units[yyy[c]] <- paste(collections$stratgroup[yyy[c]]," Group ",collections$early_interval[yyy[c]],sep="")
		} else	{
		rock_units[yyy[c]] <- paste(collections$stratgroup[yyy[c]]," Group ",collections$early_interval[yyy[c]],"-",collections$late_interval[yyy[c]],sep="")
		}
	}

rock_data <- data.frame(row.names=FALSE,stringsAsFactors = FALSE)
rock_data <- cbind(rock_data,rock_units,formation_cleaned,member_cleaned,collections$zone,collections$early_interval,collections$late_interval,collections$formation,collections$member)
colnames(rock_data) <- c("rock_unit","formation","member","zone","early_interval","late_interval","formation_identified","member_identified")
#rock_data[1,]
rockless <- (1:ttl_coll)[rock_data$rock_unit==""]
if (length(rockless)>0)	{
	lost_collections <- collections$collection_no[rockless]
	} else	{
	lost_collections <- c()
	}
rock_data <- subset(rock_data,rock_data$rock_unit!="")
sorted_rocks <- order(rock_data$rock_unit)
rock_data <- rock_data[sorted_rocks,]
rock_data <- unique(rock_data)
unique_rocks <- unique(rock_data$rock_unit)
urns <- match(unique_rocks,rock_data$rock_unit)

composite_rocks <- c()
localities <- c()
for (ur in 1:length(unique_rocks))	{
	rock_dummy <- array("",dim=ncol(rock_data))
	unique_rock_data <- subset(rock_data,rock_data$rock_unit==unique_rocks[ur])
#	urd <- data.frame(unique_rock_data[1,])
	rock_dummy[1] <- as.character(unique_rock_data$rock_unit[1])
	rock_dummy[2] <- as.character(unique_rock_data$formation[1])
	rock_dummy[3] <- as.character(unique_rock_data$member[1])
	if (nrow(unique_rock_data)>1)	{
		zones <- sort(unique(unique_rock_data$zone))
		if (zones[1]=="" || zones[1]==" ")	zones <- zones[2:length(zones)]
		if (length(zones)==1)	{
			rock_dummy[4] <- as.character(zones)
			} else if (length(zones)>1)	{
			rock_dummy[4] <- as.character(zones[1])
			for (z in 2:length(zones))	rock_dummy[4] <- paste(rock_dummy[4],", ",as.character(zones[z]),sep="")
			}
		lbs <- sort(unique(as.character(unique_rock_data$early_interval)))
		if (lbs[1]=="" || lbs[1]==" ")	lbs <- lbs[2:length(lbs)]
		if (length(lbs)==1)	{
			rock_dummy[5] <- as.character(lbs[1])
			} else if (length(lbs)>1)	{
			rock_dummy[5] <- as.character(lbs[1])
			for (z in 2:length(lbs))	rock_dummy[5] <- paste(rock_dummy[5],", ",as.character(lbs[z]),sep="")
			}
		ubs <- sort(unique(as.character(unique_rock_data$late_interval)))
		if (ubs[1]=="" || ubs[1]==" ")	{
			if (length(ubs)>1)	{
				ubs <- ubs[2:length(ubs)]
				} else	{
				ubs <- c()
				}
			}
		if (length(ubs)==1)	{
			rock_dummy[6] <- as.character(ubs[1])
			} else if (length(ubs)>1)	{
			rock_dummy[6] <- as.character(ubs[1])
			for (z in 2:length(ubs))	rock_dummy[6] <- paste(rock_dummy[6],", ",as.character(ubs[z]),sep="")
			}
		alt_form <- sort(unique(as.character(unique_rock_data$formation_identified)))
		if (!is.na(match(rock_dummy[2],alt_form)))	{
			alt_form[match(rock_dummy[2],alt_form)] <- ""
			alt_form <- sort(alt_form)
			if (length(alt_form)>1)	{
				alt_form <- alt_form[2:length(alt_form)]
				} else	{
				alt_form <- c()
				}
			}
		if (length(alt_form)==1)	{
			rock_dummy[7] <- as.character(alt_form[1])
			} else if (length(alt_form)>1)	{
			rock_dummy[7] <- as.character(alt_form[1])
			for (z in 2:length(alt_form))	rock_dummy[7] <- paste(rock_dummy[7],", ",as.character(alt_form[z]),sep="")
			}
		alt_memb <- sort(unique(as.character(unique_rock_data$member_identified)))
		if (!is.na(match(rock_dummy[3],alt_memb)))	{
			if (length(alt_memb)>1)	{
				alt_memb[match(rock_dummy[3],alt_memb)] <- ""
				alt_memb <- sort(alt_memb)
				alt_memb <- alt_memb[2:length(alt_memb)]
				} else	{
				alt_memb <- c()
				}
			}
		if (length(alt_memb)==1)	{
			rock_dummy[8] <- as.character(alt_memb[1])
			} else if (length(alt_memb)>1)	{
			rock_dummy[8] <- as.character(alt_memb[1])
			for (z in 2:length(alt_memb))	rock_dummy[8] <- paste(rock_dummy[8],", ",as.character(alt_memb[z]),sep="")
			} else	{
			rock_dummy[8] <- ""
			}
#		rock_dummy <- data.frame(rock_dummy,stringsAsFactors = TRUE)
#		rockpendium[ur,] <- as.character(rock_dummy)
		} else	{
		rock_dummy[4] <- as.character(unique_rock_data$zone[1])
		rock_dummy[5] <- as.character(unique_rock_data$early_interval[1])
		rock_dummy[6] <- as.character(unique_rock_data$late_interval[1])
		if (as.character(unique_rock_data$formation[1])!=as.character(unique_rock_data$formation_identified[1]))
			rock_dummy[7] <- as.character(unique_rock_data$formation_identified[1])
		if (as.character(unique_rock_data$member[1])!=as.character(unique_rock_data$member_identified[1]))
			rock_dummy[8] <- as.character(unique_rock_data$member_identified[1])
		}
	urc <- collections$collection_no[(1:ttl_coll)[rock_units %in% unique_rock_data$rock_unit]]
	if (length(urc)==1)	{
		localities <- c(localities,urc)
		} else if (length(urc)>1)	{
		loc <- urc[1]
		for (i in 2:length(urc))	loc <- paste(loc,", ",urc[i],sep="")
		localities <- c(localities,loc)
		}

	composite_rocks <- rbind(composite_rocks,rock_dummy)
	}
rownames(composite_rocks) <- rep("",nrow(composite_rocks))
rockpendium <- data.frame(composite_rocks,stringsAsFactors = TRUE)
colnames(rockpendium) <- colnames(rock_data)
rockpendium <- cbind(rockpendium,localities)
rownames(rockpendium) <- NULL

if (save_files)	{
	if (onset==end)	{
		output_file1 <- paste(higher_taxon,"_",onset,"_Unique_Rock_Units",output_type,sep="")
		} else{
		output_file1 <- paste(higher_taxon,"_",onset,"_to_",end,"_Unique_Rock_Units",output_type,sep="")
		}
	if (directory!="")
		output_file1 <- paste(directory,output_file1,sep="");
	write.table(rockpendium,output_file1,col.names=TRUE,row.names=FALSE,sep="\t")
	if (length(lost_collections)>0)	{
		if (onset==end)	{
			output_file2 <- paste(higher_taxon,"_",onset,"_Rockless_Collections",output_type,sep="")
			} else{
			output_file2 <- paste(higher_taxon,"_",onset,"_to_",end,"_Rockless_Collections",output_type,sep="")
			}
		if (directory!="")
			output_file2 <- paste(directory,output_file2,sep="");
		write.table(lost_collections,output_file2,col.names=FALSE,row.names=FALSE,sep="\t")
		}
	}

pendia <- list(rockpendium,lost_collections)
names(pendia) <- c("Rocks_Units","Uncertain_Collections")
return(pendia)
}

#rock_name <- "Bredasdorp & Varswater";
#rock_name <- "Bredasdorp/Varswater";
#rock_name <- "Bredasdorp and Varswater";
#rock_name <- "Bredasdorp Shale and Sandstone"
#rock_name <- c("Bredasdorp & Varswater","Bredasdorp/Varswater","Bredasdorp and Varswater","Bredasdorp Shale and Sandstone")
#ddd <- t(sapply(rock_name,separate_two_rock_names))
separate_two_rock_names <- function(rock_name)	{
rock_name <- gsub("/"," & ",rock_name);
rock_name <- gsub(" - "," & ",rock_name);
rock_name <- gsub(" or "," & ",rock_name);
if (gsub(" and "," & ",rock_name)!=rock_name)	{
	rock_molecularized <- strsplit(rock_name," ")[[1]];
	lithotypes <- (1:length(rock_molecularized))[tolower(rock_molecularized) %in% sedimentary_rocks];
	und <- match("and",tolower(rock_molecularized));
	if (length(lithotypes)>1 && und==(max(lithotypes)-1))	{

		} else	{
		rock_name <- gsub(" and "," & ",rock_name);
		}
	}
while (gsub(" &","&",rock_name)!=rock_name)
	rock_name <- gsub(" &","&",rock_name);
while (gsub("& ","&",rock_name)!=rock_name)
	rock_name <- gsub("& ","&",rock_name);
rock_molecularized <- strsplit(rock_name,"&")[[1]];
rock1 <- rock_molecularized[1];
if (length(rock_molecularized)>1)	{
	rock2 <- rock_molecularized[2];
	} else	{
	rock2 <- "";
	}
return(data.frame(rock1=as.character(rock1),rock2=as.character(rock2),stringsAsFactors = F));
}

# nnn <- paleodb_fixes$collection_no[paleodb_fixes$formation_alt!=""]
# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy.  Note that it gets only occurrences with specimen or individual occurrences
accersi_paleodb_abundance_data <- function(taxon,oldest,youngest,get_localities,file_format)	{

# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once

if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else {
	print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")}
	if(youngest %in% strat_names$nam) {
		end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
		} else {
		print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
		}

	if (oldest!=youngest)	{
		time_range <- paste(oldest,"-",youngest,sep="")
		} else {
		time_range <- oldest
		}

	httpO <- paste("https://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
	occurrences <- read.table(httpO, sep=',', header=T)
	abundances <- subset(occurrences,occurrences$abund_unit=="specimens")
	abundances <- rbind(abundances,subset(occurrences,occurrences$abund_unit=="individuals"))

	abundances$taxon_name <- mundify_taxon_names(abundances$taxon_name)
	filename <- paste(taxon,time_range,"Abundance_Records",sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
	} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
	}

if (get_localities==TRUE)  {
	localities <- unique(abundances$collection_no)
	coll <- length(localities)
	httpC <- paste("https://paleobiodb.org/data1.1/colls/list.txt?id=",localities[1],"&show=loc,time",sep="")
	collections <- read.table(httpC, sep=',', header=T)
	for (c in 2:coll) {
		httpC <- paste("https://paleobiodb.org/data1.1/colls/list.txt?id=",localities[c],"&show=loc,time",sep="")
	collections <- rbind(collections,read.table(httpC, sep=',', header=T))
	}

	collections$formation <- mundify_rock_unit_names(collections$formation,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	collections$member <- mundify_rock_unit_names(collections$member,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	collections$late_age <- expello_na_from_vector(collections$late_age, 0)
	collections$paleomodel <- expello_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2 <- expello_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3 <- expello_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4 <- expello_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2 <- expello_na_from_vector(collections$paleolng2, "")
	collections$paleolng3 <- expello_na_from_vector(collections$paleolng3, "")
	collections$paleolng4 <- expello_na_from_vector(collections$paleolng4, "")
	collections$paleolat2 <- expello_na_from_vector(collections$paleolat2, "")
	collections$paleolat3 <- expello_na_from_vector(collections$paleolat3, "")
	collections$paleolat4 <- expello_na_from_vector(collections$paleolat4, "")
	collections$geoplate2 <- expello_na_from_vector(collections$geoplate2, "")
	collections$regionalsection <- expello_na_from_vector(collections$regionalsection, "")
	collections$regionalbed <- expello_na_from_vector(collections$regionalbed, "")
	collections$regionalorder <- expello_na_from_vector(collections$regionalorder, "")
	collections$collection_subset <- expello_na_from_vector(collections$collection_subset, "")
		filename <- paste(taxon,time_range,"Collections",sep="_")
		filename <- paste(filename,file_format,sep="")
		if (file_format!=".csv")	{
			write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
			write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country <- function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,country,file_format)	{

	# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once

if(oldest %in% strat_names$nam) {
	start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
	} else {
	print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")
	}

if (youngest %in% strat_names$nam) {
	end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
	} else {
	print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
	}

if (oldest!=youngest)	{
	time_range <- paste(oldest,"-",youngest,sep="")
	} else {
	time_range <- oldest
	}

if (get_records==TRUE)	{
	if (get_abundances=="TRUE") {
		httpO <- paste("https://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
		} else {
		httpO <- paste("https://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="")
		}
	oc <- read.table(httpO, sep=',', header=T)
	occurrences <- subset(oc,oc$cc==country)
		#finds <- dim(occurrences)[1]
		# get occurrence data: if this can be modified to make it species-only, then do that.
	occurrences$taxon_name <- mundify_taxon_names(occurrences$taxon_name)
	filename <- paste(taxon,time_range,"Records_from",country,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}

if (get_localities==TRUE)	{
	httpC <- paste("https://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
	coll <- read.table(httpC, sep=',', header=T)
	collections <- subset(coll, coll$cc==country)
	collections$formation <- mundify_rock_unit_names(collections$formation,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	collections$member <- mundify_rock_unit_names(collections$member,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE);
	collections$late_interval <- expello_na_from_vector(collections$late_interval,"")
	collections$paleomodel <- expello_na_from_vector(collections$paleomodel, "")
	collections$paleomodel2 <- expello_na_from_vector(collections$paleomodel2, "")
	collections$paleomodel3 <- expello_na_from_vector(collections$paleomodel3, "")
	collections$paleomodel4 <- expello_na_from_vector(collections$paleomodel4, "")
	collections$paleolng2 <- expello_na_from_vector(collections$paleolng2, "")
	collections$paleolng3 <- expello_na_from_vector(collections$paleolng3, "")
	collections$paleolng4 <- expello_na_from_vector(collections$paleolng4, "")
	collections$paleolat2 <- expello_na_from_vector(collections$paleolat2, "")
	collections$paleolat3 <- expello_na_from_vector(collections$paleolat3, "")
	collections$paleolat4 <- expello_na_from_vector(collections$paleolat4, "")
	collections$geoplate2 <- expello_na_from_vector(collections$geoplate2, "")
	collections$regionalsection <- expello_na_from_vector(collections$regionalsection, "")
	collections$regionalbed <- expello_na_from_vector(collections$regionalbed, "")
	collections$regionalorder <- expello_na_from_vector(collections$regionalorder, "")
	collections$collection_subset <- expello_na_from_vector(collections$collection_subset, "")
	filename <- paste(taxon,time_range,"Collections_from",country,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}

if (get_taxonomy==TRUE)	{
	httpT <- paste("https://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
	taxonomy <- read.table(httpT, sep=',', header=T)
	filename <- paste(taxon,time_range,"Taxonomy",sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
		} else if (file_format==".csv")	{
		write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
		}
	}
}

# this gets occurrence, taxonomic and/or locality data for a taxon.  Enter 'TRUE' to get records, localities and/or taxonomy. This retrieves occurrences only from one country.
paleodb_data_download_by_country_two <- function(taxon,oldest,youngest,get_records,get_localities,get_abundances,get_taxonomy,environment,file_format)	{

if (environment=="terrestrial") {
  environs <- c("terrestrial indet.", "fluvial indet.", "alluvial fan", "channel lag", "coarse channel fill", "fine channel fill", "channel", "wet floodplain", "dry floodplain", "&quot;floodplain&quot;", "crevasse splay", "levee", "mire/swamp", "fluvial-lacustrine indet.", "delta plain", "fluvial-deltaic indet.", "lacustrine - large", "lacustrine - small", "pond", "crater lake", "lacustrine delta plain", "lacustrine interdistributary bay", "lacustrine delta front", "lacustrine prodelta", "lacustrine deltaic indet.", "lacustrine indet.", "dune", "interdune", "loess", "eolian indet.", "cave", "fissure fill", "sinkhole", "karst indet.", "tar", "mire/swamp", "spring", "glacial")
  } else if (environment=="reef") {
  environs <- c("reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef")
  } else if (environment=="carbonate")  {
  environs <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
  } else if (environment=="siliciclastic")  {
  environs <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  } else if (environment=="marine") {
  env1 <- c("carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)")
  env2 <- c("marine indet.","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","submarine fan","basinal (siliciclastic)","deep-water indet.")
  environs <- c(env1,env2)
  }

# get temporal information for search
strat_names <- pbdb_intervals(limit="all") #read all intervals at once

if(oldest %in% strat_names$nam) {
  start_ma <- strat_names$eag[match(oldest, strat_names$nam)]
  } else {
  print("Error! ",oldest," is not a term in our chronostratigraphic lexicon!")
  }

if (youngest %in% strat_names$nam) {
  end_ma <- strat_names$lag[match(youngest, strat_names$nam)]
  } else {
  print("Error! ",youngest," is not a term in our chronostratigraphic lexicon!")
  }

if (oldest!=youngest)	{
  time_range <- paste(oldest,"-",youngest,sep="")
  } else {
  time_range <- oldest
  }

if (get_records==TRUE)	{
  if (get_abundances=="TRUE") {
    httpO <- paste("https://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=abund,time,loc,crmod&limit=all",sep="")
    } else {
    httpO <- paste("https://paleobiodb.org/data1.1/occs/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,geo,crmod&limit=all",sep="")
    }
  oc <- read.table(httpO, sep=',', header=T)
  oc$environment <- expello_na_from_vector(oc$environment, "")
  occurrences <- subset(oc,oc$environment==environs)
  occurrences$reid_no <- expello_na_from_vector(occurrences$reid_no,"")
  occurrences$superceded <- expello_na_from_vector(occurrences$superceded,"")
  #finds <- dim(occurrences)[1]
  # get occurrence data: if this can be modified to make it species-only, then do that.
  occurrences$taxon_name <- mundify_taxon_names(occurrences$taxon_name)
  filename <- paste(taxon,time_range,"Records_from",environment,sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(occurrences,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(occurrences,filename,sep=",",eol="\n",row.names=FALSE)
    }
  }


if (get_localities==TRUE)	{
  httpC <- paste("https://paleobiodb.org/data1.1/colls/list.txt?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,geo,lithext,loc,paleoloc,crmod&limit=all",sep="")
  coll <- read.table(httpC, sep=',', header=T)
  collections <- subset(coll, coll$cc==country)
  collections$formation <- mundify_rock_unit_names(collections$formation,dehyphenate=TRUE,delete_rock_type=TRUE)
  collections$member <- mundify_rock_unit_names(collections$member,dehyphenate=TRUE,delete_rock_type=FALSE)
  collections$late_interval <- expello_na_from_vector(collections$late_interval,"")
  collections$paleomodel <- expello_na_from_vector(collections$paleomodel, "")
  collections$paleomodel2 <- expello_na_from_vector(collections$paleomodel2, "")
  collections$paleomodel3 <- expello_na_from_vector(collections$paleomodel3, "")
  collections$paleomodel4 <- expello_na_from_vector(collections$paleomodel4, "")
  collections$paleolng2 <- expello_na_from_vector(collections$paleolng2, "")
  collections$paleolng3 <- expello_na_from_vector(collections$paleolng3, "")
  collections$paleolng4 <- expello_na_from_vector(collections$paleolng4, "")
  collections$paleolat2 <- expello_na_from_vector(collections$paleolat2, "")
  collections$paleolat3 <- expello_na_from_vector(collections$paleolat3, "")
  collections$paleolat4 <- expello_na_from_vector(collections$paleolat4, "")
  collections$geoplate2 <- expello_na_from_vector(collections$geoplate2, "")
  collections$regionalsection <- expello_na_from_vector(collections$regionalsection, "")
  collections$regionalbed <- expello_na_from_vector(collections$regionalbed, "")
  collections$regionalorder <- expello_na_from_vector(collections$regionalorder, "")
  collections$collection_subset <- expello_na_from_vector(collections$collection_subset, "")
  filename <- paste(taxon,time_range,"Collections_from",country,sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(collections,filename,sep="\t",eol="\n",row.names=FALSE)
    } else if (file_format==".csv")	{
    write.table(collections,filename,sep=",",eol="\n",row.names=FALSE)
    }
}

if (get_taxonomy==TRUE)	{
  httpT <- paste("https://paleobiodb.org/data1.1/taxa/list.txt?name=",taxon,"&rel=all_children&show=attr&limit=99999&rank=genus,subgenus,species",sep="")
  taxonomy <- read.table(httpT, sep=',', header=T)
  filename <- paste(taxon,time_range,"Taxonomy",sep="_")
  filename <- paste(filename,file_format,sep="")
  if (file_format!=".csv")	{
    write.table(taxonomy,filename,sep="\t",eol="\n",row.names=FALSE)
	  } else if (file_format==".csv")	{
    write.table(taxonomy,filename,sep=",",eol="\n",row.names=FALSE)
	  }
	}
#occs.marine <- occs.raw %>% filter(cx_int_no %in% stage$interval_no) %>%filter(environment %in% marine_env)
}

accersi_information_for_list_of_collections <- function(collections_file,analysis_name="Add_These_Localities",save_file=TRUE)	{
collection_list <- read.table(collections_file,sep="")[,1]

collections_info <- c()
for (ll in 1:length(collection_list))	{
	collections_info <- rbind(collections_info,accersi_single_locality_info(locality_no=collection_list[ll]))
	}
if (save_file)	{
	output_file <- paste(analysis_name,".xls",sep="")
	write.table(collections_info,file=output_file,sep="\t",col.names = TRUE,row.names = FALSE)
	}
return(collections_info)
}

accersi_single_locality_info_old <- function(collection_no)	{
httpC <- paste("https://paleobiodb.org/data1.2/colls/list.csv?id=",collection_no,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,refattr,ent,entname,crmod",sep="");
accio <- RCurl::getURL(httpC)
coll_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
coll_info <- expello_na_from_vector(coll_info,"")
return(coll_info)
}

accersi_occurrences_from_occurrence_number <- function(occurrences_file,analysis_name="Add_These_Occurrences",save_file=TRUE)	{
occurrences_list <- read.table(occurrences_file,sep="")[,1]
occurrences_info <- c()
ll <- 0
#for (ll in 1:length(occurrences_list))	{
while (ll<=length(occurrences_list))	{
	ll <- ll+1
	occ_result <- accersi_single_occurrence_record(find_no=occurrences_list[ll])
#	occ_result[5]
	if (occ_result[1]!="No Record")	occurrences_info <- rbind(occurrences_info,occ_result)
	}
if (save_file)	{
	output_file <- paste(analysis_name,".xls",sep="")
	write.table(occurrences_info,file=output_file,sep="\t",col.names = TRUE,row.names = FALSE)
	}
return(occurrences_info)
}

accersi_single_occurrence_record <- function(find_no)	{
httpO <- paste("https://paleobiodb.org/data1.2/occs/list.csv?occ_id=",find_no,"&show=class,resgroup,ref,ent,entname,crmod",sep="")
accio <- RCurl::getURL(httpO);
j <- strsplit(accio,split="",fixed=TRUE)[[1]];
if ((j[2]=="o" && j[3]=="c") && (j[443]!="N" && j[444]!="O")) 	{
	coll_finds <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
	coll_finds <- expello_na_from_vector(coll_finds,"")
	} else {
	coll_finds <- array("No Record",dim=c(1))
	}
return(coll_finds)
}

								### Ecological Data ###
accersi_ecologic_guild_data <- function(taxon,onset="Cambrian",end="Holocene",taxon_level="genus,subgenus",directory="",save_files=TRUE,output_type=".csv")	{
taxon <- gsub(" ","%20",taxon);
taxon <- paste(taxon, collapse = ",");
taxon_level <- paste(taxon_level,collapse=",");
#http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=refattr,classext,rem,entname&limit=all",sep = "");
http <- paste("https://www.paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&rank=",taxon_level,"&interval=",onset,",",end,"&show=attr,app,ecospace,etbasis",sep="");
fetch <- RCurl::getURL(http);
taxon_guilds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
taxon_guilds <- expello_na_from_matrix(taxon_guilds,"");

taxon_level <- strsplit(taxon_level,split=",")[[1]];
ntaxa <- nrow(taxon_guilds)
nn <- (1:ntaxa)[!taxon_guilds$accepted_rank %in% taxon_level]
#taxon_guilds[nn[1:2],];
#dim(taxon_guilds);
#dim(taxon_guilds[taxon_guilds$accepted_rank %in% taxon_level,])
if (save_files)	{
	if (onset!=end)	{
		output <- paste(directory,onset,"-",end,"_",taxon,"_Basic_Ecological_Data",output_type,sep="");
		} else	{
		output <- paste(directory,end,"_",taxon,"_Basic_Ecological_Data",output_type,sep="");
		}
	if (output_type==".csv")	{
		write.csv(taxon_guilds,output,row.names = F,fileEncoding = "UTF-8");
		} else	{
		write.table(taxon_guilds,output,row.names = F,sep="\t");
		}
	}

return(taxon_guilds);
}

								### Taphonomic Data ###
accersi_taphonomic_data <- function(taxon,onset="Cambrian",end="Holocene",taxon_level="genus,subgenus",directory="",save_files=TRUE,output_type=".csv")	{
http <- paste("https://www.paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&rank=",taxon_level,"&show=attr,app,ttaph,etbasis",sep="")
fetch <- RCurl::getURL(http);
taphonomic_grades <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
taphonomic_grades <- expello_na_from_matrix(taphonomic_grades,"");
if (save_files)	{
	if (onset!=end)	{
		output <- paste(directory,onset,"-",end,"_",taxon,"_Basic_Taphonomic_Data",output_type,sep="");
		} else	{
		output <- paste(directory,end,"_",taxon,"_Basic_Taphonomic_Data",output_type,sep="");
		}
	if (output_type==".csv")	{
		write.csv(taphonomic_grades,output,row.names = F,fileEncoding = "UTF-8");
		} else	{
		write.table(taphonomic_grades,output,row.names = F,sep="\t");
		}
	}
return(taphonomic_grades);
}

								### Measurements ###
accersi_measurements <- function(taxon,onset="Cambrian",end="Holocene",directory="",save_files=TRUE,output_type=".csv")	{
http <- paste("https://www.paleobiodb.org/data1.2/specs/measurements.csv?base_name=",taxon,"&interval=",onset,",",end,"&show=full",sep="");
fetch <- RCurl::getURL(http);
measurements <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
measurements <- expello_na_from_matrix(measurements,"");
if (save_files)	{
	if (onset!=end)	{
		output <- paste(directory,onset,"-",end,"_",taxon,"_Measurements",output_type,sep="");
		} else	{
		output <- paste(directory,end,"_",taxon,"_Measurements",output_type,sep="");
		}
	if (output_type==".csv")	{
		write.csv(measurements,output,row.names = F,fileEncoding = "UTF-8");
		} else	{
		write.table(measurements,output,row.names = F,sep="\t");
		}
	}
return(measurements);
}

#paleodb_finds <- clade_finds
add_subgenus_names_to_paleodb_finds <- function(paleodb_finds)	{
# routine to provide a separate column giving subgenus names to paleodb occurence records
genus_name <- paleodb_finds$genus[paleodb_finds$genus!=""];
genus_subgenus <- pbapply::pbsapply(genus_name,divido_subgenus_names_from_genus_names);
paleodb_finds$genus[paleodb_finds$genus!=""] <- as.character(genus_subgenus[1,]);
subgenus <- rep("",nrow(paleodb_finds));
subgenus[paleodb_finds$genus!=""] <- as.character(genus_subgenus[2,]);
subgenus[subgenus==""] <- paleodb_finds$genus[subgenus==""];
if (is.na(match("subgenus",colnames(paleodb_finds))))	{
	paleodb_finds <- tibble::add_column(paleodb_finds, subgenus=as.character(subgenus), .after = match("genus_no",colnames(paleodb_finds)));
	} else	{
	paleodb_finds$subgenus <- as.character(subgenus);
	}
return(paleodb_finds);
}

#pbdb_finds <- clade_finds_orig;
add_subgenera_to_pbdb_finds <- function(pbdb_finds,both_subgenus_names=T)	{
# if strip_genera=T, then Loxoplocus (Lophospira) -> Loxoplocus
# if strip_genera=F, then Loxoplocus (Lophospira) -> Loxoplocus (Lophospira)
genus <- sort(unique(as.character(pbdb_finds$genus)));
genus_subgenus <- t(sapply(genus,divido_subgenus_names_from_genus_names));
if (both_subgenus_names) genus_subgenus[genus_subgenus[,2]!="",] <- rownames(genus_subgenus[genus_subgenus[,2]!="",]);
subgenus <- genus_subgenus[match(pbdb_finds$genus,rownames(genus_subgenus)),2]
if (is.na(match("subgenus",colnames(pbdb_finds))))	{
	pbdb_finds <- tibble::add_column(pbdb_finds, subgenus=as.character(subgenus), .after = match("genus_no",colnames(pbdb_finds)));
	} else	{
	pbdb_finds$subgenus <- as.character(subgenus);
	}
pbdb_finds$genus <- genus_subgenus[match(pbdb_finds$genus,genus),1];

#pbdb_finds <- tibble::add_column(pbdb_finds, subgenus, .after = match("genus_no",colnames(pbdb_finds)));
return(pbdb_finds);
}

# remove specified taxon_qualifiers
# (defaults to list above; to be used to save particular ones if needs be)
unqualify_identified_taxon <- function(identified_taxon,taxon_qualifiers=taxon_qualifiers)	{
for (tq in 1:length(taxon_qualifiers))
	identified_taxon <- gsub(taxon_qualifiers[tq],"",identified_taxon);
return(identified_taxon);
}

					##### TAXONOMIC DOWNLOAD ROUTINES #######
# taxon <- "Diplograptina";
# given a list of SPECIES (or subspecies), get the taxonomic data that the PaleoDB has (if any)
accersi_taxonomic_data_for_one_taxon <- function(taxon,inc_children=F,inc_parents=F,oldest="",youngest="")	{
# taxon: taxon name;
# inc_subtaxa: if T, then taxonomic data for all constituent taxa are included; if F, then only for this taxon.
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
if (oldest!="" && youngest=="")	{
	youngest <- "Cenozoic";
	} else if (youngest!="" && oldest=="")	{
	oldest <- "Hadean";
	}
if (oldest!="")	{
	interval <- paste("&interval=",oldest,",",youngest,sep="");
	} else	{
	interval <- "";
	}
if (inc_children)	{
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&variant=all&rel=all_children",interval,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr",sep="");
	} else if (inc_parents)	{
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&variant=all&rel=all_parents",interval,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr",sep="");
	} else	{
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&show=attr,parent,class,classext,refattr,crmod",sep="");
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,                   interval,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,ref,refattr",sep="");
	}
#             https://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Cypraeidae&rank=species,subspecies&op_type=all
#fetch <- RCurl::getURL(httpT);
#taxon_information <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
taxon_information <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");

if (!is.na(match(orig_taxon,taxon_information$taxon_name)) && !inc_children)	{
	taxon_information <- subset(taxon_information,taxon_information$taxon_name==orig_taxon);
	if (nrow(subset(taxon_information,taxon_information$flags=="B"))==1)
		taxon_information <- subset(taxon_information,taxon_information$flags=="B");
	if (nrow(taxon_information)>1)	{
		if (nrow(subset(taxon_information,taxon_information$difference==""))==1)	{
			taxon_information <- subset(taxon_information,taxon_information$difference=="");
			}
		}
	} #else if (!is.na(match(orig_taxon,taxon_information)))
taxon_information <- expello_na_from_matrix(data=taxon_information);
return(taxon_information);
}

accersi_taxonomic_data_for_one_taxon_no <- function(taxon_no,inc_children=F,inc_parents=F,oldest="",youngest="")	{
# taxon: taxon name;
# inc_subtaxa: if T, then taxonomic data for all constituent taxa are included; if F, then only for this taxon.
if (oldest!="" && youngest=="")	{
	youngest <- "Cenozoic";
	} else if (youngest!="" && oldest=="")	{
	oldest <- "Hadean";
	}
if (oldest!="")	{
	interval <- paste("&interval=",oldest,",",youngest,sep="");
	} else	{
	interval <- "";
	}
if (inc_children)	{
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?id=txn:",taxon_no,"&variant=all&rel=all_children",interval,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,pres,ref",sep="");
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?id=txn:",taxon_no,"&rel=all_children&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,pres,ref,refattr",sep="");
	} else if (inc_parents)	{
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?id=txn:",taxon_no,"&rel=all_parents&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,pres,ref,refattr",sep="");
	} else	{
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?match_name=",taxon,"&show=attr,parent,class,classext,refattr,crmod",sep="");
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/single.csv?id=txn:",taxon_no,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,pres,ref,refattr",sep="");
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/single.txt?id=txn:",taxon_no,interval,"&show=attr,common",sep="");
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?id=txn:=",taxon_no,interval,"&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,pres,ref,refattr",sep="");
	}
#             https://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Cypraeidae&rank=species,subspecies&op_type=all
#fetch <- RCurl::getURL(httpT);
#taxon_information <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
#options(warn=-1);
#test <- readline(httpT);
#options(warn=-1);
#if (is.matrix(test))	{
	taxon_information <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#	} else	{
#	recast <- as.numeric(simplify2array(strsplit(test[length(test)],"\"")));
#	recast <- recast[!is.na(recast)];
#	taxon_no <- recast[1];
#	taxon_information <- accersi_taxonomic_data_for_one_taxon_no(taxon_no,inc_children,inc_parents,oldest,youngest)
#	}
#if (!is.na(match(orig_taxon,taxon_information$taxon_name)) && !inc_children)	{
#	taxon_information <- subset(taxon_information,taxon_information$taxon_name==orig_taxon);
#	if (nrow(subset(taxon_information,taxon_information$flags=="B"))==1)
#		taxon_information <- subset(taxon_information,taxon_information$flags=="B");
#	if (nrow(taxon_information)>1)	{
#		if (nrow(subset(taxon_information,taxon_information$difference==""))==1)	{
#			taxon_information <- subset(taxon_information,taxon_information$difference=="");
#			}
#		}
#	} #else if (!is.na(match(orig_taxon,taxon_information)))
taxon_information <- expello_na_from_matrix(data=taxon_information);
return(taxon_information);
}

are_you_my_mommy <- function(taxon_no)	{
httpT <- paste("https://paleobiodb.org/data1.2/taxa/single.csv?id=txn:",taxon_no,"&show=parent,immparent",sep="");
taxon_information <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
taxon_information <- put_pbdb_dataframes_into_proper_type(taxon_information);
return(taxon_information);
}

accersi_taxonomic_data_for_list_of_species <- function(species_list) {
taxon <- species_list <- unique(species_list[species_list!=""]);
species <- pbapply::pbsapply(taxon,mundify_taxon_names);
ttl_spc <- length(species_list);
#species <- c();
#for (g in 1:length(species_list))	{
#	species_list[g] <- gsub("lower ","",species_list[g]);
#	species_list[g] <- gsub("middle ","",species_list[g]);
#	species_list[g] <- gsub("upper ","",species_list[g]);
#	species <- c(species,transmogrify_full_zone_names_to_species_names_only(species_list[g]));
#	}
#species_list <- unique(species_list);
#species_list <- species_list[order(species)];
#species <- species[order(species)];
#species_list <- species_list[species!=""];
#species <- species[species!=""];
entered_species <- unentered_species <- data.frame();
for (g in 1:ttl_spc)	{
	taxon <- species_list[g];
	taxon_information <- accersi_taxonomic_data_for_one_taxon(taxon);
	if (sum(taxon_information[,1] %in% "THIS REQUEST RETURNED NO RECORDS")>0)	{
		unentered_species <- rbind(unentered_species,c(species_list[g],species[g]));
		} else	{
		if (length(entered_species)==0)	{
			entered_species <- taxon_information;
			} else	{
			entered_species <- rbind(entered_species,taxon_information);
			}
		}
	}

if (nrow(unentered_species)>0)
	colnames(unentered_species) <- c("species","species_epithet");
output <- list(entered_species,unentered_species);
names(output) <- c("entered_species","unentered_species");
return(output);
}

# updated 2021-02-04: better homonym routine
accersi_taxonomic_data_for_list_of_taxa <- function(taxon_list) {
taxon_list <- unique(pbapply::pbsapply(taxon_list,mundify_taxon_names));

entered_taxon <- unentered_taxon <- homonyms <- c();
#for (g in 1:length(taxon_list))	{
g <- 0;
while (g < length(taxon_list))	{
	g <- g+1;
	taxon <- gsub(" ","%20",taxon_list[g]);
	taxon_information <- accersi_taxonomic_data_for_one_taxon(taxon);
	taxon_information <- taxon_information[!taxon_information$difference %in% nomens,];
	# routine for taxa somehow entered twice #
	if (nrow(taxon_information)>1)	{
		ttl_names <- nrow(taxon_information);
		attribution <- gsub("\\(","",taxon_information$taxon_attr);
		attribution <- gsub("\\)","",attribution);
		name_attribution <- unique(cbind(taxon_information$taxon_name,attribution));
		if (nrow(name_attribution)<nrow(taxon_information))	{
			keeper1 <- (1:nrow(taxon_information))[taxon_information$taxon_name %in% name_attribution[,1]];
			keeper2 <- (1:nrow(taxon_information))[attribution %in% name_attribution[,2]];
			keeper <- keeper1[keeper1 %in% keeper2];
			if (sum(taxon_information$n_occs[keeper]>0)>0 && sum(taxon_information$n_occs[keeper]<length(keeper)))	{
				cutter <- keeper[taxon_information$n_occs[keeper]==0];
				taxon_information <- taxon_information[(1:ttl_names)[!(1:ttl_names) %in% cutter],];
				} else if (sum(sum(taxon_information$n_occs[keeper]>0)==0))	{
				cutter <- keeper[2:length(keeper)];
				taxon_information <- taxon_information[(1:ttl_names)[!(1:ttl_names) %in% cutter],];
				}
			}
		}
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,"&variant=all&rel=all_children&show=ref,refattr",sep="");
#	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,"&variant=all&rel=all_children&show=attr,common,app,parent,immparent,size,class,classext,subcounts,ecospace,taphonomy,etbasis,pres,img,ref,refattr,crmod",sep="");
#             https://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Cypraeidae&rank=species,subspecies&op_type=all
#	fetch <- RCurl::getURL(httpT);
#	taxon_information <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
	relv_info <- (1:nrow(taxon_information))[taxon_information$taxon_name %in% taxon_list[g]];
	if (ncol(taxon_information)>2)	{
		if (is.null(entered_taxon) && length(relv_info)==1)	{
			entered_taxon <- taxon_information[relv_info,];
			} else if (length(relv_info)>1)	{
			homonyms <- c(homonyms,g);
			} else	{
			if (is.na(relv_info[1]))	{
				if (nrow(taxon_information)==1)	relv_info <- 1;
				}
			entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
			}
		}	else	{
		taxon <- gsub("%20"," ",taxon);
		unentered_taxon <- rbind(unentered_taxon,taxon);
		taxon_split <- simplify2array(strsplit(taxon," "));
		if (length(taxon_split)>1)	{
			genus_name <- divido_genus_names_from_species_names(taxon);
			subgenus_test <- divido_subgenus_names_from_genus_names(genus_name);
			if (is.na(match(subgenus_test[1],taxon_list)))
				taxon_list <- c(taxon_list,subgenus_test[1]);
			if (is.na(match(genus_name,taxon_list)) && subgenus_test[2]!="" && subgenus_test[1]!=subgenus_test[2])
				taxon_list <- c(taxon_list,genus_name);
			}
		}
#	print(unentered_taxon[g]);
	}

# if there are two legit names, then first see if we can
#	match the family, order or class of the other taxa.
# if not, then see if we can fit it in the right rank
hm <- 0;
while (hm < length(homonyms))	{
	hm <- hm+1;
	g <- homonyms[hm];
	taxon <- gsub(" ","%20",taxon_list[g]);
	taxon_information <- accersi_taxonomic_data_for_one_taxon(taxon);
	taxon_information <- taxon_information[!taxon_information$difference %in% nomens,];
	taxon_information <- taxon_information[taxon_information$taxon_name==taxon_list[g],];
	synos <- nrow(taxon_information);
	relv_info <- (1:synos)[taxon_information$family %in% entered_taxon$family];
	if (length(relv_info)==1)	{
		entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
		} else	{
		relv_info <- (1:synos)[taxon_information$order %in% entered_taxon$order];
		if (length(relv_info)==1)	{
			entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
			} else	{
			relv_info <- (1:synos)[taxon_information$class %in% entered_taxon$class];
			if (length(relv_info)==1)	{
				entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
				} else	{
				relv_info <- c();
				}
			}
		}
	if (length(relv_info)==0)	{
		relv_info <- (1:synos)[taxon_information$taxon_rank %in% entered_taxon$taxon_rank];
		if (length(relv_info)==1)	{
			entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
			} else	{
			taxon <- gsub("%20"," ",taxon);
			unentered_taxon <- rbind(unentered_taxon,taxon);
			}
		}
	}

output <- list(entered_taxon,unentered_taxon);
names(output) <- c("entered_taxa","unentered_taxa");
return(output);
}

# updated 2021-02-04: better homonym routine
# updated 2021-04-15: added sapply version.
# updated 2021-04-20:
accersi_taxonomic_data_for_list_of_taxa_prototype <- function(taxon_list) {
print("Cleaning taxon names.")
taxon <- taxon_list <- unique(pbapply::pbsapply(taxon_list,mundify_taxon_names));
print("Downloading taxonomic information from the PBDB")
taxon_info <- data.frame(simplify2array(base::t(pbapply::pbsapply(taxon[1:4],accersi_taxonomic_data_for_one_taxon))));
#for (cl in 1:ncol(taxon_info))	taxon_info[,cl] <- simplify2array(taxon_info[,cl]);

entered_taxon <- unentered_taxon <- c();
homonyms <- c();
for (g in 1:length(taxon_info))	{
	if (length(taxon_info[[g]])==2)	{
		unentered_taxon <- c(unentered_taxon,taxon_list[g]);
		} else	{
		if (length(entered_taxon)==0)	{
			entered_taxon <- taxon_info[[g]];
			} else	{
			entered_taxon <- rbind(entered_taxon,taxon_info[[g]]);
			}
		if (nrow(taxon_info[[g]])>1)	homonyms <- c(homonyms,taxon_list[g]);
		}
	}
entered_taxon$taxon_attr <- gsub("\\(","",entered_taxon$taxon_attr);
entered_taxon$taxon_attr <- gsub("\\)","",entered_taxon$taxon_attr);
entered_taxon <- put_pbdb_dataframes_into_proper_type(pbdb_data=entered_taxon);
# if there are two legit names, then first see if we can
#	match the family, order or class of the other taxa.
# if not, then see if we can fit it in the right rank
hm <- 0;
while (hm < length(homonyms))	{
	hm <- hm+1;
	g <- homonyms[hm];
	taxon <- gsub(" ","%20",taxon_list[g]);
	entered_taxon <- accersi_taxonomic_data_for_one_taxon(taxon);
	entered_taxon <- entered_taxon[!entered_taxon$difference %in% nomens,];
	entered_taxon <- entered_taxon[entered_taxon$taxon_name==taxon_list[g],];
	synos <- nrow(entered_taxon);
	relv_info <- (1:synos)[entered_taxon$family %in% entered_taxon$family];
	if (length(relv_info)==1)	{
		entered_taxon <- rbind(entered_taxon,entered_taxon[relv_info,]);
		} else	{
		relv_info <- (1:synos)[entered_taxon$order %in% entered_taxon$order];
		if (length(relv_info)==1)	{
			entered_taxon <- rbind(entered_taxon,entered_taxon[relv_info,]);
			} else	{
			relv_info <- (1:synos)[entered_taxon$class %in% entered_taxon$class];
			if (length(relv_info)==1)	{
				entered_taxon <- rbind(entered_taxon,entered_taxon[relv_info,]);
				} else	{
				relv_info <- c();
				}
			}
		}
	if (length(relv_info)==0)	{
		relv_info <- (1:synos)[entered_taxon$taxon_rank %in% entered_taxon$taxon_rank];
		if (length(relv_info)==1)	{
			entered_taxon <- rbind(entered_taxon,entered_taxon[relv_info,]);
			} else	{
			taxon <- gsub("%20"," ",taxon);
			unentered_taxon <- rbind(unentered_taxon,taxon);
			}
		}
	}

output <- list(entered_taxon,unentered_taxon);
names(output) <- c("entered_taxa","unentered_taxa");
return(output);
}

accersi_taxonomic_data_for_list_of_taxon_nos <- function(taxon_no_list) {
#SEQ  <- seq(1,100);
#pb   <- txtProgressBar(1, 100, style=3);
#TIME <- Sys.time();
#for(i in SEQ)	{
#	}
ntaxa <- length(taxon_no_list);
entered_taxon <- unentered_taxon <- c();
# there really shouldn't be an unentered taxa given numbers: but who knows what whackiness people might enter....
#pts <- round(ntaxa*(1:100)/100,0);
for (g in 1:length(taxon_no_list))	{
#	if (g %in% pts)	{
#		Sys.sleep(0.02);
#		setTxtProgressBar(pb, i)
#		}
	taxon_information <- accersi_taxonomic_data_for_one_taxon_no(taxon_no=taxon_no_list[g]);
	entered_no <- data.frame(entered_no=as.numeric(taxon_no_list[g]),stringsAsFactors = F);
	taxon_information <- cbind(entered_no,taxon_information);
	if (g==1)	{
#		entered_taxon <- taxon_information[relv_info,];
		entered_taxon <- taxon_information;
		} else if (ncol(taxon_information)>2)	{
#		entered_taxon <- rbind(entered_taxon,taxon_information[relv_info,]);
		entered_taxon <- rbind(entered_taxon,taxon_information);
		}
	}
#output <- list(entered_taxon,unentered_taxon);
#names(output) <- c("entered_taxa","unentered_taxa");
entered_taxon <- entered_taxon[!is.na(entered_taxon$entered_no),];
entered_taxon <- put_pbdb_dataframes_into_proper_type(pbdb_data=entered_taxon);
#entered_no <- data.frame(entered_no=as.numeric(taxon_no_list));
#entered_taxon <- cbind(entered_no,entered_taxon);
return(entered_taxon);
}

# get the taxonomic opinions that the PaleoDB has (if any) for a taxon given its name
accersi_taxonomic_opinions_for_one_taxon <- function(taxon,exact_match=T) {
# taxon: name of taxon
# exact_match: if true, then get only this taxon; if false, then get all constituent taxa, too.
taxon <- gsub(" ","%20",taxon);
taxon_opinions <- c();
if (exact_match)	{
	httpTO <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTO <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}
#accio <- RCurl::getURL(httpTO);
taxon_opinions <- read.csv(httpTO,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
return(taxon_opinions);
}

# get the taxonomic opinions that the PaleoDB has (if any) for a taxon given its name
accersi_taxonomic_opinions_for_one_taxon_no <- function(taxon_no,exact_match=T) {
# taxon_no: paleodb taxon number for a taxon
# exact_match: if true, then get only this taxon; if false, then get all constituent taxa, too.
taxon_opinions <- c();
if (exact_match)	{
	httpTO <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?taxon_id=",taxon_no,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTO <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_id=",taxon_no,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}
#accio <- RCurl::getURL(httpTO);
#https://www.paleobiodb.org/data1.2/opinions/list.json?author=Osborn
taxon_opinions <- read.csv(httpTO,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
return(taxon_opinions);
}

# get the latest authoritative opinion for a taxon based on PaleoDB number
accersi_latest_opinion_for_taxon_no <- function(taxon_no,exact_match=T)	{
# add something to make sure that the most recent assignment is still a valid higher taxon
	#print(taxon);
opinions <- accersi_taxonomic_opinions_for_one_taxon_no(taxon_no);
if (ncol(opinions)==2)	{
	xxx <- accersi_taxonomic_opinions_for_one_taxon("Lophospira perangulata");
	xxx[1,] <- "";
	xxx <- xxx[1,];
	return(xxx);
	} else	{
	opinion_ranks <- match(opinions$basis,opinion_basis_rank);
	opinions <- subset(opinions,opinion_ranks==min(opinion_ranks))
	opinions <- subset(opinions,opinions$pubyr==max(opinions$pubyr));
	if (nrow(opinions)>1)	{
		if (!is.na(match("class",opinions$opinion_type)))	{
			opinions <- opinions[match("class",opinions$opinion_type),];
			} else if (length(unique(opinions$parent_no))==1)	{
			opinions <- opinions[1,];
			} else	{
			possible_parents <- accersi_taxonomic_data_for_list_of_taxon_nos(unique(opinions$parent_no));
			taxon_ranks <- match(possible_parents$taxon_rank,accersi_taxonomic_ranks())
			opinions <- opinions[match(min(taxon_ranks),taxon_ranks),];
			}
		}
#	if (sum(opinions$basis=="stated with evidence")>0)	opinions <- subset(opinions,opinions$basis=="stated with evidence");
#	opinions <- subset(opinions,opinions$status=="belongs to");
	return(opinions);
	}
}

# get opinion by PBDB opinion number #
# get the taxonomic opinions that the PaleoDB has (if any) for a taxon given its name
#opinion_no <- 856773;
accersi_taxonomic_opinion_no <- function(opinion_no) {
# taxon_no: paleodb taxon number for a taxon
# exact_match: if true, then get only this taxon; if false, then get all constituent taxa, too.
#/data1.2/opinions/single.json?id=1000&show=entname
httpTO <- paste("https://www.paleobiodb.org/data1.2/opinions/single.txt?id=",opinion_no,"&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
httpTO <- paste("https://www.paleobiodb.org/data1.2/opinions/single.txt?id=",opinion_no,"&show=basis,ref,refattr,ent,entname,crmod",sep="");
#accio <- RCurl::getURL(httpTO);
#https://www.paleobiodb.org/data1.2/opinions/list.json?author=Osborn
taxon_opinions <- read.csv(httpTO,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
return(taxon_opinions);
}

# get the latest authoritative opinion for a taxon based on PaleoDB number
entangle_disconnected_genus_and_subgenus <- function(taxon_nos)	{
# add something to make sure that the most recent assignment is still a valid higher taxon
	#print(taxon);
opinions <- accersi_taxonomic_opinions_for_one_taxon_no(taxon_no=taxon_nos[1]);
for (tn in 2:length(taxon_nos))
	opinions <- rbind(opinions,accersi_taxonomic_opinions_for_one_taxon_no(taxon_no=taxon_nos[tn]));

opinion_ranks <- match(opinions$basis,opinion_basis_rank);
opinions <- subset(opinions,opinion_ranks==min(opinion_ranks))
opinions <- subset(opinions,opinions$pubyr==max(opinions$pubyr));
return(opinions);
}

# get the latest authoritative opinion for a list of taxon numbers
accersi_latest_opinion_for_taxon_no_list <- function(taxon_nos)	{
# add something to make sure that the most recent assignment is still a valid higher taxon
	#print(taxon);
ntaxa <- length(taxon_nos);
for (n in 1:ntaxa)	{
	if (n==1)	{
		opinions <- accersi_latest_opinion_for_taxon_no(taxon_no=taxon_nos[n]);
		} else	{
		opinions <- rbind(opinions,accersi_latest_opinion_for_taxon_no(taxon_no=taxon_nos[n]));
		}
	}
return(opinions);
}

# get the latest authoritative opinion for a taxon based on name
accersi_latest_opinion_for_taxon <- function(taxon,exact_match=T)	{
# add something to make sure that the most recent assignment is still a valid higher taxon
#print(taxon);
opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon,exact_match);
if (ncol(opinions)==2)	{
	xxx <- accersi_taxonomic_opinions_for_one_taxon("Lophospira perangulata");
	xxx[1,] <- "";
	xxx <- xxx[1,];
	return(xxx);
	} else	{
	if (sum(opinions$basis=="stated with evidence")>0)	opinions <- subset(opinions,opinions$basis=="stated with evidence");
	best_and_latest <- match(max(opinions$pubyr),opinions$pubyr);
	return(opinions[best_and_latest,]);
	}
}

# get a list of all the names by which a species is known
revelare_species_aliases <- function(species_name)	{
#opinions <- accersi_taxonomic_opinions_for_one_species(species_name);
opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon=species_name);
if (is.null(opinions))	{
	all_names <- species_name;
	} else	{
	all_names <- unique(c(as.character(opinions$taxon_name),as.character(opinions$child_name)));
	all_names <- all_names[!is.na(all_names)];
	all_names <- all_names[all_names!=""];
	}
return(all_names)
}

# generate a dictionary with the most recent opinions for species in PaleoDB;
accersi_species_dictionary <- function(paleodb_finds)	{
taxon_no <- entered_species_nos <- unique(paleodb_finds$accepted_no[!paleodb_finds$accepted_no %in% c(paleodb_finds$genus_no,paleodb_finds$subgenus_no)]);
entered_species_nos <- entered_species_nos[order(entered_species)];

species_dictionary <- c();
for (tl in 1:length(entered_species))
	species_dictionary <- rbind(species_dictionary,accersi_latest_opinion_for_taxon_no(taxon_no=entered_species_nos[tl]));

species_dictionary <- expello_na_from_matrix(species_dictionary,"");
species_dictionary$child_name[species_dictionary$child_name==""] <- species_dictionary$taxon_name[species_dictionary$child_name==""];
species_dictionary <- species_dictionary[species_dictionary$opinion_no!="",];
keepers <- match(c("opinion_no","taxon_rank","taxon_name","orig_no","child_name","child_spelling_no","status","parent_name","parent_no","parent_spelling_no","spelling_reason","basis","author","pubyr","ref_author","ref_pubyr","reference_no","primary_reference"),colnames(species_dictionary))
return(species_dictionary[,keepers]);
}

# generate a thesaurus that gives different binomial combinations for species in PaleoDB;
accersi_species_thesaurus <- function(paleodb_finds)	{
taxon_no <- entered_species_nos <- unique(paleodb_finds$accepted_no[!paleodb_finds$accepted_no %in% c(paleodb_finds$genus_no,paleodb_finds$subgenus_no)]);
entered_species_nos <- entered_species_nos[order(entered_species)];

species_thesaurus <- c();
for (tl in 1:length(entered_species))	{
	everyone_has_one <- accersi_taxonomic_opinions_for_one_taxon_no(taxon_no=entered_species_nos[tl],exact_match = T);
	current_genus <- divido_genus_names_from_species_names(everyone_has_one$taxon_name[1]);
	check_genus <- accersi_latest_opinion_for_taxon(current_genus,exact_match = T);
	if (check_genus$status != "belongs to")
		current_genus <- check_genus$parent_name;
	everyone_has_one$child_name <- expello_na_from_vector(everyone_has_one$child_name,replacement="");
	this_spc <- rbind(c(everyone_has_one$taxon_name[1],everyone_has_one$orig_no[1],everyone_has_one$taxon_name[1],everyone_has_one$orig_no[1],current_genus),
		cbind(everyone_has_one$child_name,everyone_has_one$child_spelling_no,rep(everyone_has_one$taxon_name[1],nrow(everyone_has_one)),rep(everyone_has_one$orig_no[1],nrow(everyone_has_one)),rep(current_genus,nrow(everyone_has_one))));
	this_spc <- this_spc[match(unique(this_spc[,1]),this_spc[,1]),]
	keepers <- (1:nrow(this_spc))[this_spc[,1]!=""];
	if (tl==1)	{
		species_thesaurus <- data.frame(species_name=as.character(this_spc[keepers,1]),name_no=as.numeric(this_spc[keepers,2]),species=as.character(this_spc[keepers,3]),taxon_no=as.numeric(this_spc[keepers,4]),genus=as.character(this_spc[keepers,5]));
		} else	{
		this_spc <- data.frame(species_name=as.character(this_spc[keepers,1]),name_no=as.numeric(this_spc[keepers,2]),species=as.character(this_spc[keepers,3]),taxon_no=as.numeric(this_spc[keepers,4]),genus=as.character(this_spc[keepers,5]));
		species_thesaurus <- rbind(species_thesaurus,this_spc);
		}
	}
return(species_thesaurus);
}

accersi_taxonomic_ranks <- function()	{
taxonomic_rank <- "subspecies"
taxonomic_rank <- c(taxonomic_rank,"species")
taxonomic_rank <- c(taxonomic_rank,"subgenus")
taxonomic_rank <- c(taxonomic_rank,"genus")
taxonomic_rank <- c(taxonomic_rank,"subtribe")
taxonomic_rank <- c(taxonomic_rank,"tribe")
taxonomic_rank <- c(taxonomic_rank,"subfamily")
taxonomic_rank <- c(taxonomic_rank,"family")
taxonomic_rank <- c(taxonomic_rank,"superfamily")
taxonomic_rank <- c(taxonomic_rank,"infraorder")
taxonomic_rank <- c(taxonomic_rank,"suborder")
taxonomic_rank <- c(taxonomic_rank,"order")
taxonomic_rank <- c(taxonomic_rank,"superorder")
taxonomic_rank <- c(taxonomic_rank,"infraclass")
taxonomic_rank <- c(taxonomic_rank,"subclass")
taxonomic_rank <- c(taxonomic_rank,"class")
taxonomic_rank <- c(taxonomic_rank,"superclass")
taxonomic_rank <- c(taxonomic_rank,"subphylum")
taxonomic_rank <- c(taxonomic_rank,"phylum")
taxonomic_rank <- c(taxonomic_rank,"superphylum")
taxonomic_rank <- c(taxonomic_rank,"subkingdom")
taxonomic_rank <- c(taxonomic_rank,"kingdom")
taxonomic_rank <- c(taxonomic_rank,"superkingdom")

return(taxonomic_rank)
}

accersi_parents_and_parent_ranks <- function(taxa,taxon_rank)	{
ntaxa <- length(taxa)
if (taxon_rank=="" || taxon_rank=="all"|| taxon_rank=="All")	{
	accioranks <- ""
	} else if (taxon_rank=="species")	{
	accioranks <- "&rank=species";
	} else if (taxon_rank=="genus")	{
	accioranks <- "&rank=genus,subgenus";
	} else if (taxon_rank=="family")	{
	accioranks <- "&rank=family";
	}
find_taxa <- paste(taxa, collapse = ",");
http <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",find_taxa,accioranks,"&op_type=class&private",sep="");
fetch <- RCurl::getURL(http);
opinions <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");

parents <- as.character(opinions$parent_name[match(taxa,opinions$taxon_name)]);
unique_parents <- sort(unique(parents));

if (taxon_rank=="" || taxon_rank=="all" || taxon_rank=="All")	{
	accioranks2 <- "&rank=genus,subgenus";
	} else if (taxon_rank=="species")	{
	accioranks2 <- "&rank=genus,subgenus";
	} else if (taxon_rank=="genus")	{
	accioranks2 <- "&rank=family";
	} else if (taxon_rank=="family")	{
	accioranks2 <- "&rank=min_family";
	}
find_taxa2 <- paste(unique_parents, collapse = ",");
http2 <- paste("https://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",find_taxa2,accioranks2,"&op_type=class&private",sep="");
fetch2 <- RCurl::getURL(http2);
opinions2 <- utils::read.csv(text = fetch2, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8");
used_opinions2 <- opinions2[match(unique_parents,opinions2$taxon_name),];
parent_ranks <- as.character(used_opinions2$taxon_rank[match(parents,used_opinions2$taxon_name)]);
return(data.frame(cbind(taxa,parents,parent_ranks)));
}

accersi_classification_tree <- function(taxon,high_rank,low_rank,file_format=".csv")	{
taxon_ranks <- accersi_taxonomic_ranks();
start <- match(high_rank,taxon_ranks);
end <- match(low_rank,taxon_ranks);
tranks <- taxon_ranks[start];
for (t in (start-1):end)	tranks <- paste(tranks,",",taxon_ranks[t],sep="");
http <- paste("https://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,sep="");
http <- paste(http,"&rank=",tranks,sep="");
#if (oldest!="")		http <- paste(http,"&interval=",oldest,",",youngest)
http <- paste(http,"&show=attr,parent",sep="")
taxonomy <- read.table(http, sep=',', header=T)
genera <- subset(taxonomy,taxonomy$taxon_rank=="genus")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
#word(subgenera$taxon_name[1],2)
subgenus_names <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenus_genera <- word(subgenera$taxon_name,1)
combined <- cbind(subgenus_genera,subgenus_names)
genus_names <- as.character(genera$taxon_name)
uniq_subgenus_genera <- unique(subgenus_genera)
gsg <- length(uniq_subgenus_genera)
np <- prob <- 0
prob_gen <- ""
for (g in 1:gsg)	{
	n <- match(uniq_subgenus_genera[g],genus_names)
	if (is.na(n)) {
		prob <- prob+1
		mtch <- 0
		if (prob==1)	{
			prob_gen <- uniq_subgenus_genera[g]
			} else prob_gen <- c(prob_gen,uniq_subgenus_genera[g])
		}	else {
		np <- np+1
		mtch <- genera$taxon_no[n]
		}
	if (g==1)	{
		matchgen_no <- mtch
		} else matchgen_no <- c(matchgen_no,mtch)
#	print(c(g,mtch,uniq_subgenus_genera[g]))
	}
if (prob>0)	{
	for (p in 1:prob)	{
		orig_mtch <- match(prob_gen[p],combined[,2])
		if (is.na(orig_mtch))	{
			xxx <- expello_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			line1 <- paste(combined[yyy,2]," is assigned to: ",combined[yyy,1]," (",combined[yyy,2],")",sep="")
			line2 <- paste("   however, ",combined[yyy,1]," is not a currently accepted genus name or spelling.",sep="")
			print(line1)
			print(line2)
			}	else {
			cur_comb <- paste(combined[orig_mtch,1]," (",combined[orig_mtch,2],")",sep="")
			sgn <- match(cur_comb,subgenera$taxon_name)
			author <- str_replace_all(as.character(subgenera$taxon_attr[sgn]), "[[:punct:]]", "")
			if (is.na(author))	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else if (author=="")	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else {
				line1 <- paste(combined[orig_mtch,2],author)
				}
			line1 <- paste(line1," currently is ",combined[orig_mtch,1]," (",combined[orig_mtch,2],"), yet we also have:",sep="")
			print(line1)
			xxx <- expello_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			if (length(yyy)==1)	{
				dy <- yyy
				yyy <- vector(length=1)
				yyy[1] <- dy
				}
			for (i in 1:length(yyy))	{
				line2 <- paste("     ",combined[yyy[i],1]," (",combined[yyy[i],2],")",sep="")
				print(line2)
				}
			}
		}
	}	else print("All genera with subgenera are last ranked as genera.")
}

accersi_classification <- function(taxon,ranks)	{
#https://www.paleobiodb.org/data1.2/taxa/list.csv?base_name=Lophospira&rank=genus,subgenus&private&show=class
http1 <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&rank=",ranks,"&show=class",sep="")
classification <- read.csv(http1,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#accio <- RCurl::getURL(http1)
#classification <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
return(classification)
}

revelare_inconsistent_genus_and_subgenus_ranks <- function(taxon,file_format)	{
#taxon: a higher taxon that includes genera and subgenera
#file_format:
#	csv: comma-delimited
#	txt: tab-delimited text
#	tab: tab-delimited text
#	xls: tab-delimited text that usually will be opened directly by Excel
tranks <- "genus,subgenus"
http <- paste("https://paleobiodb.org/data1.2/taxa/list.txt?base_name=",taxon,sep="")
http <- paste(http,"&rank=",tranks,sep="")
http <- paste(http,"&show=attr,parent",sep="")
taxonomy <- read.table(http, sep=',', header=T)
genera <- subset(taxonomy,taxonomy$taxon_rank=="genus")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
subgenus_names <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenus_genera <- word(subgenera$taxon_name,1)
combined <- cbind(subgenus_genera,subgenus_names)
genus_names <- as.character(genera$taxon_name)
uniq_subgenus_genera <- unique(subgenus_genera)
gsg <- length(uniq_subgenus_genera)
np <- prob <- 0
prob_gen <- ""
for (g in 1:gsg)	{
	n <- match(uniq_subgenus_genera[g],genus_names)
	if (is.na(n)) {
		prob <- prob+1
		mtch <- 0
		if (prob==1)	{
			prob_gen <- uniq_subgenus_genera[g]
			} else prob_gen <- c(prob_gen,uniq_subgenus_genera[g])
		}	else {
		np <- np+1
		mtch <- genera$taxon_no[n]
		}
	if (g==1)	{
		matchgen_no <- mtch
		} else matchgen_no <- c(matchgen_no,mtch)
	}
if (prob>0)	{
	output <- vector(length=1)
	outlines <- 1
	for (p in 1:prob)	{
		orig_mtch <- match(prob_gen[p],combined[,2])
		if (is.na(orig_mtch))	{
			xxx <- expello_na_from_vector(match(combined[,1],prob_gen[p]),0)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			sg <- paste(combined[yyy,1]," (",combined[yyy,2],")",sep="")
			nn <- match(sg,taxonomy$taxon_name)
			new_parent <- taxonomy$parent_name[nn]
			line1 <- paste(combined[yyy,2]," is assigned to: ",combined[yyy,1]," (",combined[yyy,2],")",sep="")
			line2 <- paste("however, ",combined[yyy,1]," is either synonymized into or corrected as ",new_parent,".",sep="")
			print(line1)
			if (outlines==1)	{
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
				output <- c(output,line1)
				}
			outlines <- outlines+1
			dum <- paste("     ",line2,sep="")
			print(dum)
			dum <- paste("",line2,sep="\t")
			output <- c(output,dum)
			}	else {
			cur_comb <- paste(combined[orig_mtch,1]," (",combined[orig_mtch,2],")",sep="")
			sgn <- match(cur_comb,subgenera$taxon_name)
			author <- str_replace_all(as.character(subgenera$taxon_attr[sgn]), "[[:punct:]]", "")
			if (is.na(author))	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else if (author=="")	{
				line1 <- paste(combined[orig_mtch,2],"(author unknown)")
				} else {
				line1 <- paste(combined[orig_mtch,2],author)
				}
			line1 <- paste(line1," currently is ",combined[orig_mtch,1]," (",combined[orig_mtch,2],"), yet we also have:",sep="")
			print(line1)
			if (outlines==1)	{
#				output[outlines,1] <- "st"
				output[outlines] <- line1
				}	else {
				outlines <- outlines+1
#				dum[1] <- "st"
#				dum[2] <- line1
				output <- c(output,line1)
				}
			xxx <- expello_na_from_vector(match(combined[,1],prob_gen[p]),0)
			ttl <- sum(xxx)
			for (i in 1:length(xxx))	xxx[i] <- i*(xxx[i])
			yyy <- subset(xxx,xxx>0)
			if (length(yyy)==1)	{
				dy <- yyy
				yyy <- vector(length=1)
				yyy[1] <- dy
				}
			for (j in 1:length(yyy))	{
				line2 <- paste(combined[yyy[j],1]," (",combined[yyy[j],2],")",sep="")
				dum <- paste("     ",line2,sep="")
				print(dum)
				dum <- paste("",line2,sep="\t")
				outlines <- outlines+1
				output <- c(output,dum)
				}
			}
		} # end problem cases
	outfile <- paste(taxon,"_genus_and_subgenus_conflicts.",file_format,sep="")
	if (file_format=="csv")	{
		flsp <- ","
		} else flsp <- "\t"
	write(output,outfile);
	}	else print("All genera with subgenera are last ranked as genera.")
}

revelare_taxonomy_for_one_taxon <- function(taxon,settle=F)	{
dud <- F;
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=full,immparent,classext,attr,app",sep="");
## If the taxon is absent, then something weird will happen: kill it with fire.
taxon_info <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#accio <- RCurl::getURL(httpT);
#taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"));
if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS" || nrow(which(taxon_info=="THIS REQUEST RETURNED NO RECORDS",arr.ind = T))>0))	{
#	ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0));
	dud <- T;
	# if unknown species, then get information about the genus
	if (settle && length(strsplit(orig_taxon,split=" ")[[1]])>1)	{
		taxon <- strsplit(orig_taxon,split=" ")[[1]][1];
		httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=full,immparent,classext,attr,app",sep="");
		taxon_info <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
		taxon_info <- taxon_info[match(taxon,taxon_info$taxon_name),];	# get rid of species & subgenera
		if (nrow(which(taxon_info=="THIS REQUEST RETURNED NO RECORDS",arr.ind = T))==0)
			dud <- F;
		}
	if (dud)	{
		httpTT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=Lophospira&show=full,immparent,classext,attr,app",sep="");
		accioTT <- RCurl::getURL(httpTT);
		dummy_info <- data.frame(utils::read.csv(text = accioTT, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"));
#	output_labels <- c("orig_no","taxon_no","record_type","flags","taxon_rank","taxon_name","taxon_attr","difference","accepted_no","accepted_rank","accepted_name","parent_no","reference_no","is_extant","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma","early_interval","late_interval");
		output_labels <- colnames(dummy_info);
		taxon_info <- data.frame(matrix("",1,length(output_labels)),stringsAsFactors = FALSE);
		colnames(taxon_info) <- output_labels;
		taxon_info$taxon_name <- orig_taxon;
		taxon_info$taxon_attr <- taxon_info$accepted_name <- "?";
		taxon_info$n_occs <- 0;
		}
	}
if (nrow(taxon_info)>1)	{
	orig_info <- taxon_info;
	taxon_info <- subset(taxon_info, taxon_info$taxon_name==orig_taxon);
	if (nrow(taxon_info)==0)	{
		taxon_info <- orig_info;
		opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon=orig_taxon);
		ttl_opinions <- nrow(opinions);
		orig_name <- (1:ttl_opinions)[opinions$child_name %in% orig_taxon];
		taxon_info <- taxon_info[match(opinions$taxon_name[orig_name],taxon_info$taxon_name),];
#		taxon_info$n_occs <- 0;
		}
	}
return(taxon_info);
}

#taxon <- "Mattheviidae"; max_rank <- "class";
accersi_parents_up_to_rank_X_for_one_taxon <- function(taxon,max_rank="phylum")	{
accersi_taxonomic_data_for_one_taxon()
}

#taxon_no <- 76763; max_rank <- "class";
accersi_parents_up_to_rank_X_for_one_taxon_no <- function(taxon_no,max_rank="phylum")	{
taxon_info <- accersi_taxonomic_data_for_one_taxon_no(taxon_no);
parent_info <- data.frame(taxon_no=as.numeric(taxon_info$parent_no),
						  taxon=as.character(taxon_info$parent_name),
						  taxon_rank=as.character(taxon_info$taxon_rank),
						  taxon_rank_no=as.numeric(match(taxon_info$taxon_rank,taxonomic_rank)),stringsAsFactors = F);
max_rank_no <- match(max_rank,taxonomic_rank);

while(max(parent_info$taxon_rank_no)<max_rank_no)	{
	taxon_no <- parent_info$taxon_no[nrow(parent_info)];
	taxon_info <- accersi_taxonomic_data_for_one_taxon_no(taxon_no);
	parent_info$taxon_rank[nrow(parent_info)] <- taxon_info$taxon_rank;
	parent_info$taxon_rank_no[nrow(parent_info)] <- match(taxon_info$taxon_rank,taxonomic_rank);
	if (parent_info$taxon_rank_no[nrow(parent_info)]>match("superkingdom",taxonomic_rank))
		parent_info$taxon_rank_no[nrow(parent_info)] <- -parent_info$taxon_rank_no[nrow(parent_info)];
	new_info <- data.frame(taxon_no=as.numeric(taxon_info$parent_no),
						   taxon=as.character(taxon_info$parent_name),
						   taxon_rank=as.character(taxon_info$taxon_rank),
						   taxon_rank_no=as.numeric(match(taxon_info$taxon_rank,taxonomic_rank)),stringsAsFactors = F);
	parent_info <- rbind(parent_info,new_info);
	}
return(parent_info[1:(1+sum(parent_info$taxon_rank_no<max_rank_no)),]);
}

revelare_parent_taxon_for_one_taxon <- function(taxon)	{
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=immparent",sep="");
## If the taxon is absent, then something weird will happen: kill it with fire.
taxon_info <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#accio <- RCurl::getURL(httpT);
#taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"));
if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
#	ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0));
	output_labels <- c("orig_no","taxon_no","record_type","flags","taxon_rank","taxon_name","taxon_attr","difference","accepted_no","accepted_rank","accepted_name","parent_no","reference_no","is_extant","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma","early_interval","late_interval");
	taxon_info <- data.frame(matrix("",1,length(output_labels)),stringsAsFactors = FALSE);
	colnames(taxon_info) <- output_labels;
	taxon_info$taxon_name <- orig_taxon;
	taxon_info$taxon_attr <- taxon_info$accepted_name <- "?";
	taxon_info$n_occs <- 0;
	}
if (nrow(taxon_info)>1)	{
	orig_info <- taxon_info;
	taxon_info <- subset(taxon_info, taxon_info$taxon_name==orig_taxon);
	if (nrow(taxon_info)==0)	{
		taxon_info <- orig_info;
#		opinions <- accersi_taxonomic_opinions_for_one_species(species_name=orig_taxon);
		opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon=orig_taxon);
		ttl_opinions <- nrow(opinions);
		orig_name <- (1:ttl_opinions)[opinions$child_name %in% orig_taxon];
		taxon_info <- taxon_info[match(opinions$taxon_name[orig_name],taxon_info$taxon_name),];
		taxon_info$n_occs <- 0;
		}
	}
return(taxon_info$parent_name);
}

revelare_parent_taxon_no_for_one_taxon <- function(taxon)	{
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=immparent",sep="");
## If the taxon is absent, then something weird will happen: kill it with fire.
taxon_info <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#accio <- RCurl::getURL(httpT);
#taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"));
if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
#	ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0));
	output_labels <- c("orig_no","taxon_no","record_type","flags","taxon_rank","taxon_name","taxon_attr","difference","accepted_no","accepted_rank","accepted_name","parent_no","reference_no","is_extant","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma","early_interval","late_interval");
	taxon_info <- data.frame(matrix("",1,length(output_labels)),stringsAsFactors = FALSE);
	colnames(taxon_info) <- output_labels;
	taxon_info$taxon_name <- orig_taxon;
	taxon_info$taxon_attr <- taxon_info$accepted_name <- "?";
	taxon_info$n_occs <- 0;
	}
if (nrow(taxon_info)>1)	{
	orig_info <- taxon_info;
	taxon_info <- subset(taxon_info, taxon_info$taxon_name==orig_taxon);
	if (nrow(taxon_info)==0)	{
		taxon_info <- orig_info;
#		opinions <- accersi_taxonomic_opinions_for_one_species(species_name=orig_taxon);
		opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon=orig_taxon);
		ttl_opinions <- nrow(opinions);
		orig_name <- (1:ttl_opinions)[opinions$child_name %in% orig_taxon];
		taxon_info <- taxon_info[match(opinions$taxon_name[orig_name],taxon_info$taxon_name),];
		taxon_info$n_occs <- 0;
		}
	}
return(taxon_info$parent_no);
}

revelare_taxon_rank_for_one_taxon <- function(taxon)	{
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=full",sep="");
## If the taxon is absent, then something weird will happen: kill it with fire.
taxon_info <- read.csv(httpT,header=T,stringsAsFactors=hell_no,encoding="UTF-8");
#accio <- RCurl::getURL(httpT);
#taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"));
if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
#	ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0));
	output_labels <- c("orig_no","taxon_no","record_type","flags","taxon_rank","taxon_name","taxon_attr","difference","accepted_no","accepted_rank","accepted_name","parent_no","reference_no","is_extant","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma","early_interval","late_interval");
	taxon_info <- data.frame(matrix("",1,length(output_labels)),stringsAsFactors = FALSE);
	colnames(taxon_info) <- output_labels;
	taxon_info$taxon_name <- orig_taxon;
	taxon_info$taxon_attr <- taxon_info$accepted_name <- "?";
	taxon_info$n_occs <- 0;
	}
if (nrow(taxon_info)>1)	{
	orig_info <- taxon_info;
	taxon_info <- subset(taxon_info, taxon_info$taxon_name==orig_taxon);
	if (nrow(taxon_info)==0)	{
		taxon_info <- orig_info;
#		opinions <- accersi_taxonomic_opinions_for_one_species(species_name=orig_taxon);
		opinions <- accersi_taxonomic_opinions_for_one_taxon(taxon=orig_taxon);
		ttl_opinions <- nrow(opinions);
		orig_name <- (1:ttl_opinions)[opinions$child_name %in% orig_taxon];
		taxon_info <- taxon_info[match(opinions$taxon_name[orig_name],taxon_info$taxon_name),];
		taxon_info$n_occs <- 0;
		}
	}
return(taxon_info$taxon_rank);
}

accersi_occurrences_for_list_of_taxa_old <- function(taxon_list,analysis_name="",output_file=TRUE)	{
# remove any funny symbols from taxon names
#taxon_name <- taxon_list
taxon_list <- sapply(taxon_list,mundify_taxon_names);
ntaxa <- length(taxon_list)

occompendium <- data.frame(row.names=FALSE)	# occurrences
collpendium <- data.frame(row.names=FALSE)		# collections
ttl_finds <- data.frame(row.names=FALSE)
#entries <- array(0,ntaxa)
#tx <- 1
#while (tx <= ntaxa)	{
for (tx in 1:ntaxa)	{
#	print(paste(taxon_list[tx],", ",tx," of ",ntaxa,sep=""))
	#paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon_list[tx],"&show=coll,strat,refattr&limit=all",sep = "")	accio <- RCurl::getURL(http1)
	taxon <- gsub(" ","%20",taxon_list[tx]);
	httpT <- paste("https://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=attr,app",sep="")
	## If the taxon is absent, then something weird will happen: kill it with fire.
	accio <- RCurl::getURL(httpT);
	taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
	if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS"))	{
		ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0))
		}	else	{
		http <- paste("https://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxon,"&show=full",sep="")
		accio <- RCurl::getURL(http)
		taxon_finds <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8"))
		if (nrow(taxon_finds)==1 && taxon_finds$collection_no=="THIS REQUEST RETURNED NO RECORDS")	{
			n_occs <- 0
			} else	{
			n_occs <- nrow(taxon_finds)
			if (n_occs==1)	{
				taxon_finds <- expello_na_from_vector(taxon_finds,"")
				} else	{
				taxon_finds <- expello_na_from_matrix(taxon_finds, "")
				}
			collpendium <- rbind(collpendium,taxon_finds)
			occompendium <- rbind(occompendium,cbind(taxon_finds$collection_no,rep(taxon_list[tx],n_occs)))
			}
		ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],taxon_info$n_occs[1]))
		}
#	tx <- tx+1
	}

# reduce collections to just unique ones & sort them
novel_coll <- sort(unique(collpendium$collection_no))
xx <- (1:nrow(collpendium))[match(novel_coll,collpendium$collection_no)]
collpendium <- collpendium[xx,]
# put column names on the occurrences & taxon finds
colnames(occompendium) <- c("collection_no","taxon")
colnames(ttl_finds) <- c("taxon","occurrences")

if (output_file)	{
	output_file1 <- paste(analysis_name,"_Collections.xls",sep="")
	output_file2 <- paste(analysis_name,"_Occurrences.xls",sep="")
	output_file3 <- paste(analysis_name,"_Finds_per_Taxon.xls",sep="")
	write.table(collpendium,output_file1,col.names=TRUE,row.names=FALSE,sep="\t")
	write.table(occompendium,output_file2,col.names=TRUE,row.names=FALSE,sep="\t")
	write.table(ttl_finds,output_file3,col.names=TRUE,row.names=FALSE,sep="\t")
	}
pendia <- list(occompendium,collpendium,ttl_finds)
names(pendia) <- c("Occurrences","Collections","Sampling")
return(pendia)
}

accersi_oldest_members_of_a_taxon_ideal <- function(taxon)	{
eras <- c("Cryogenian","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Paleogene","Neogene")
finds <- 0
stg <- 1
while (finds==0 && stg<length(eras))	{
	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.tsv?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],sep="")
#	occurrences <- read.table(httpO, sep='\t', header=T)
#	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],"&show=stratext",sep="")
	occurrences <- read.table(httpO, sep=',', header=T)
	if (occurrences[1,1]=="THIS REQUEST RETURNED NO RECORDS")	{
		stg <- stg+1
		}	else {
		finds <- dim(occurrences)[1]
		httpO <- paste("https://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&interval=",eras[stg],",",eras[stg],"&show=stratext",sep="")
		occurrences <- read.table(httpO, sep=',', header=T)
#		for (i in 1:dim(occurrences)[1])	for (j in 1:dim(occurrences[2]))	if (is.na(data[i,j]))	occurrences[i,j] <- ""
		}
	}
elder <- order(occurrences$max_ma,decreasing=TRUE)
filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
write.table(occurrences[elder,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
}

accersi_oldest_species_of_a_taxon_ideal <- function(taxon,contrain_oldest="Archaean",give_or_take=20,output_file=TRUE)	{
intervals <- c("Archaean","Proterozoic","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Cenozoic")
good <- FALSE
i <- match(contrain_oldest,intervals)
#"https://paleobiodb.org/data1.2/occs/list.txt?base_name=Rostroconchia&interval=Ediacaran&show=stratext"
while (!good)	{
	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",intervals[i],"&show=stratext",sep = "")
	fetch <- RCurl::getURL(httpO)
	all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
	if(nrow(all_finds)==1 && all_finds[1]=="THIS REQUEST RETURNED NO RECORDS")	{
		i <- i+1
		}	else	{
		species_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
		if (nrow(species_finds)>0)	{
			good <- TRUE
			} else	i <- i+1
		}
	}
maxma <- max(species_finds$max_ma)
elder <- subset(species_finds,species_finds$max_ma>=(maxma-(maxma/give_or_take)))
minma <- min(species_finds$max_ma)
if (abs(maxma-minma)<give_or_take && intervals[i]!="Cenozoic")	{
	i<-i+1
	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",intervals[i],"&show=stratext",sep = "")
	fetch <- RCurl::getURL(httpO)
	new_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
	new_spc_finds <- subset(new_finds,new_finds$identified_rank=="species")
	elder <- rbind(elder,subset(new_spc_finds,new_spc_finds$max_ma>=(maxma-(maxma/give_or_take))))
	}
eldest <- order(elder$max_ma,decreasing=TRUE)
if (output_file)	{
	filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
	write.table(elder[eldest,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
	}
return(elder[eldest,])
}

accersi_youngest_species_of_a_taxon_ideal <- function(taxon,contrain_youngest="Cenozoic",give_or_take=20,output_file=TRUE)	{
intervals <- c("Archaean","Proterozoic","Ediacaran","Cambrian","Ordovician","Silurian","Devonian","Carboniferous","Permian","Triassic","Jurassic","Cretaceous","Cenozoic")
good <- FALSE
i <- match(contrain_oldest,intervals)
#"https://paleobiodb.org/data1.2/occs/list.txt?base_name=Rostroconchia&interval=Ediacaran&show=stratext"
while (!good)	{
	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",intervals[i],"&show=stratext",sep = "")
	fetch <- RCurl::getURL(httpO)
	all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
	if(nrow(all_finds)==1 && all_finds[1]=="THIS REQUEST RETURNED NO RECORDS")	{
		i <- i+1
		}	else	{
		species_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
		if (nrow(species_finds)>0)	{
			good <- TRUE
			} else	i <- i+1
		}
	}
maxma <- max(species_finds$max_ma)
elder <- subset(species_finds,species_finds$max_ma>=(maxma-(maxma/give_or_take)))
minma <- min(species_finds$max_ma)
if (abs(maxma-minma)<give_or_take && intervals[i]!="Cenozoic")	{
	i<-i+1
	httpO <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",intervals[i],"&show=stratext",sep = "")
	fetch <- RCurl::getURL(httpO)
	new_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=FALSE,encoding="UTF-8")
	new_spc_finds <- subset(new_finds,new_finds$identified_rank=="species")
	elder <- rbind(elder,subset(new_spc_finds,new_spc_finds$max_ma>=(maxma-(maxma/give_or_take))))
	}
eldest <- order(elder$max_ma,decreasing=TRUE)
if (output_file)	{
	filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
	write.table(elder[eldest,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
	}
return(elder[eldest,])
}

accersi_oldest_species_of_a_taxon <- function(taxon,give_or_take=20)	{
httpO <- paste("https://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxon,"&show=stratext",sep="")
occurrences <- read.table(httpO, sep=',', header=T)
occurrences <- subset(occurrences,occurrences$identified_rank=="species")
maxma <- max(occurrences$max_ma)
elder <- subset(occurrences,occurrences$max_ma>=(maxma-(maxma/give_or_take)))
eldest <- order(elder$max_ma,decreasing=TRUE)
filename1 <- paste(taxon,"Oldest_Records.csv",sep="_")
write.table(elder[eldest,],filename1,sep=",",eol="\n",row.names=FALSE, col.names=TRUE)
return(elder[eldest,])
}

# get common higher taxon #
accersi_taxonomic_root <- function(taxon_list)	{
initial_taxonomic_data <- accersi_taxonomic_data_for_list_of_taxa(taxon_list);
parent_nos <-  unique(initial_taxonomic_data$entered_taxa$parent_no);
parent_info <- accersi_taxonomic_data_for_list_of_taxon_nos(parent_nos);
# start here!
parents <- data.frame(orig_no=as.numeric(parent_info$orig_no),
					  taxon_no=as.numeric(parent_info$taxon_no),
					  taxon=as.character(parent_info$taxon_name),
					  rank=as.character(parent_info$taxon_rank),
					  rank_no=as.numeric(match(parent_info$taxon_rank,taxonomic_rank)),
					  parent_no=as.numeric(parent_info$parent_no),
					  stringsAsFactors = hell_no);
unfound_parent_nos <- parents$parent_no[!parents$parent_no %in% c(parents$orig_no,parents$taxon_no)];
while (length(unfound_parent_nos)>1)	{
	grandparent_info <- accersi_taxonomic_data_for_list_of_taxon_nos(unfound_parent_nos);
	grandparents <- data.frame(orig_no=as.numeric(grandparent_info$orig_no),
							   taxon_no=as.numeric(grandparent_info$taxon_no),
							   taxon=as.character(grandparent_info$taxon_name),
							   rank=as.character(grandparent_info$taxon_rank),
							   rank_no=as.numeric(match(grandparent_info$taxon_rank,taxonomic_rank)),
							   parent_no=as.numeric(grandparent_info$parent_no),
							   stringsAsFactors = hell_no);
	parents <- unique(rbind(parents,grandparents));
	unfound_parent_nos <- parents$parent_no[!parents$parent_no %in% c(parents$orig_no,parents$taxon_no)];
	}
nested <- (1:nrow(parents))[!is.na(match(parents$parent_no,parents$orig_no))];
parents$parent_no[nested] <- parents$taxon_no[match(parents$parent_no[nested],parents$orig_no)];
root_no <- parents$taxon_no[parents$parent_no==unfound_parent_nos]
while(sum(parents$parent_no==root_no)==1)	{
	root_no <- parents$taxon_no[parents$parent_no==root_no];
	}
root <- parents$taxon[match(root_no,parents$taxon_no)];
return(root);
}

#grandparents$parent_no <- parents$parent_no[match(grandparents$taxon,parents$taxon)];
# take taxon without parent in parents
# count down until you get 2+ daughter taxa
#sum(grandparents$taxon_no[1] %in% parents$parent_no)


#root <- parents$taxon[match(grandparents$parent_no[grandparents$parent_no %in% parents$orig_no],parents$orig_no)];
#root_no <- parents$taxon_no[match(grandparents$parent_no[grandparents$parent_no %in% parents$orig_no],parents$orig_no)];



#root <- parents$taxon[parents$parent_no %in% unfound_parent_nos[1]];
#if (is.na(root))
#	root <- parents$taxon[parents$orig_no %in% unfound_parent_nos[1]];

accersi_taxonomic_hierarchy <- function(taxon_list)	{
initial_taxonomic_data <- accersi_taxonomic_data_for_list_of_taxa(taxon_list);
orig_taxa <- data.frame(orig_no=as.numeric(initial_taxonomic_data$entered_taxa$orig_no),
						taxon_no=as.numeric(initial_taxonomic_data$entered_taxa$taxon_no),
						taxon=as.character(initial_taxonomic_data$entered_taxa$taxon_name),
						rank=as.character(initial_taxonomic_data$entered_taxa$taxon_rank),
						rank_no=as.numeric(match(initial_taxonomic_data$entered_taxa$taxon_rank,taxonomic_rank)),
						parent_no=as.numeric(initial_taxonomic_data$entered_taxa$parent_no),
						parent=as.character(initial_taxonomic_data$entered_taxa$taxon_name[match(initial_taxonomic_data$entered_taxa$parent_no,initial_taxonomic_data$entered_taxa$orig_no)]),
						stringsAsFactors = hell_no);
parent_nos <-  unique(initial_taxonomic_data$entered_taxa$parent_no);
parent_info <- accersi_taxonomic_data_for_list_of_taxon_nos(parent_nos);
# start here!
parents <- data.frame(orig_no=as.numeric(parent_info$orig_no),
					  taxon_no=as.numeric(parent_info$taxon_no),
					  taxon=as.character(parent_info$taxon_name),
					  rank=as.character(parent_info$taxon_rank),
					  rank_no=as.numeric(match(parent_info$taxon_rank,taxonomic_rank)),
					  parent_no=as.numeric(parent_info$parent_no),
					  parent=as.character(parent_info$taxon_name[match(parent_info$parent_no,parent_info$orig_no)]),
					  stringsAsFactors = hell_no);
orig_taxa$parent <- parents$taxon[match(orig_taxa$parent_no,parents$orig_no)]
parents$parent[is.na(parents$parent)] <- parent_info$taxon_name[match(parent_info$parent_no[is.na(parents$parent)],parent_info$taxon_no)];
parents$parent[is.na(parents$parent)] <- "";
unfound_parent_nos <- unique(parents$parent_no[!parents$parent_no %in% c(parents$orig_no,parents$taxon_no)]);
orig_taxa <- rbind(orig_taxa,parents);
while (length(unfound_parent_nos)>1)	{
	grandparent_info <- accersi_taxonomic_data_for_list_of_taxon_nos(unfound_parent_nos);
	grandparents <- data.frame(orig_no=as.numeric(grandparent_info$orig_no),
							   taxon_no=as.numeric(grandparent_info$taxon_no),
							   taxon=as.character(grandparent_info$taxon_name),
							   rank=as.character(grandparent_info$taxon_rank),
							   rank_no=as.numeric(match(grandparent_info$taxon_rank,taxonomic_rank)),
							   parent_no=as.numeric(grandparent_info$parent_no),
							   parent=as.character(grandparent_info$taxon_name[match(grandparent_info$parent_no,grandparent_info$orig_no)]),
							   stringsAsFactors = hell_no);
	grandparents$parent[is.na(grandparents$parent)] <- "";
	orig_taxa <- unique(rbind(orig_taxa,grandparents));
	orig_taxa$parent <- orig_taxa$taxon[match(orig_taxa$parent_no,orig_taxa$orig_no)];
	orig_taxa$parent[is.na(orig_taxa$parent)] <- "";
#	match(parents$parent_no,parents$taxon_no)
	unfound_parent_nos <- unique(orig_taxa$parent_no[!orig_taxa$parent_no %in% c(orig_taxa$orig_no,orig_taxa$taxon_no)]);
	}
min_rank <- min(orig_taxa$rank_no);
orig_taxa$parent_rank_no <- as.numeric(orig_taxa$rank_no[match(orig_taxa$parent,orig_taxa$taxon)]);
orig_taxa$parent_rank_no[is.na(orig_taxa$parent_rank_no)] <- 26;
all_ranks <- sort(unique(orig_taxa$rank_no));

hierarchy <- subset(orig_taxa,orig_taxa$rank_no==min_rank);
hierarchy <- hierarchy[order(-hierarchy$parent_rank_no,hierarchy$parent_no),];
remainder <- subset(orig_taxa,!orig_taxa$orig_no %in% hierarchy$orig_no);
unranked_taxa <- unique(hierarchy$parent_no[nrow(hierarchy):1])
unranked_taxa <- unranked_taxa[!unranked_taxa %in% hierarchy$orig_no];
unranked_taxa <- unranked_taxa[!unranked_taxa %in% hierarchy$taxon_no];
#for (ur in 1:length(unranked_taxa))	{
# something goes NA crazy here.
while (ur < length(unranked_taxa))	{
	ur <- ur+1;
	hi_tx_no <- unranked_taxa[ur];
	higher_taxon <- orig_taxa[match(hi_tx_no,orig_taxa$orig_no),]
	subtaxa <- hierarchy$taxon_no[hierarchy$parent_no==hi_tx_no];
	st <- 1;
	while (st <= length(subtaxa))	{
		subtaxa <- c(subtaxa,hierarchy$taxon_no[hierarchy$parent_no==subtaxa[st]]);
		st <- st+1;
		}
	# find where this taxa first starts
	within <- sort(match(subtaxa,hierarchy$taxon_no));
	rn <- match(hi_tx_no,hierarchy$parent_no);
	below <- (rn:nrow(hierarchy))[!hierarchy$taxon_no[rn:nrow(hierarchy)] %in% subtaxa]
	if (rn>1)	{
		above <- 1:(rn-1);
		} else	{
		above <- c();
		}
	hierarchy <- rbind(hierarchy[above,],higher_taxon,hierarchy[within,],hierarchy[below,]);
	if (!higher_taxon$parent_no %in% unranked_taxa)	unranked_taxa <- c(unranked_taxa,higher_taxon$parent_no);
	unranked_taxa <- unranked_taxa[!is.na(unranked_taxa)];
	}

#for (ar in 1:length(all_ranks))	{
#	this_rank <- subset(orig_taxa,orig_taxa$rank_no==min_rank);
#	this_rank[order(-this_rank$parent_rank_no,this_rank$parent_no),]
#	}

#root <- orig_taxa$taxon[orig_taxa$parent_no %in% unfound_parent_nos[1]];
#if (is.na(root))
#	root <- orig_taxa$taxon[orig_taxa$orig_no %in% unfound_parent_nos[1]];
return(root);
}

# routine to return gendern neutral version of species name: taxon_name="Bonnia parvula"
spay_or_neuter_your_genus_species_combo <- function(taxon_name)	{
genus_name <- divido_genus_names_from_species_names(taxon_name);
epithet <- divido_species_epithets(taxon_name);
return(paste(genus_name,spay_or_neuter_your_species(epithet)));
}

# routine to return gendern neutral version of species name
spay_or_neuter_your_species <- function(species_epithet)	{
# standardize gender of species epithets
#If (Right (species_epithet; 2)="us";Left(species_epithet; Length (species_epithet )-2)  &  "a";
#	If ( Right (species_epithet; 2)="um"; Left(species_epithet; Length (species_epithet)-2)  &  "a";
#		 If ( Right (species_epithet; 2)="is"; Left(species_epithet; Length (species_epithet)-2)  &  "e";
#		 	 If ( Right (species_epithet; 2)="ei"; Left(species_epithet; Length (species_epithet)-2)  &  "ii"; species_epithet))))j <- simplify2array(strsplit(taxon_name," ")[[1]]);
nnames <- length(strsplit(species_epithet,split=" ")[[1]]);
if (nnames==2)	{
	species_epithet <- simplify2array(strsplit(species_epithet,split=" ")[[1]]);
	return(paste(sapply(species_epithet,spay_or_neuter_your_species),collapse=" "));
	} else if (nnames==1)	{
	spc_ep <- simplify2array(strsplit(species_epithet,split="")[[1]]);
	t_l <- length(spc_ep);
	p_l <- t_l-1;
	pp_l <- p_l-1;
	if (paste(spc_ep[p_l:t_l],collapse="")=="us")	{
		spc_ep[t_l] <- "m";
		} else if (spc_ep[t_l]=="a")	{
		spc_ep[t_l] <- "u";
		spc_ep <- c(spc_ep,"m")
		} else if (paste(spc_ep[p_l:t_l],collapse="")=="is")	{
		spc_ep[p_l] <- "e";
		spc_ep <- spc_ep[1:p_l];
		} else if (paste(spc_ep[p_l:t_l],collapse="")=="ii" || paste(spc_ep[p_l:t_l],collapse="")=="ei")	{
		spc_ep[p_l] <- "e";
		spc_ep <- spc_ep[1:p_l];
		} else if (spc_ep[t_l]=="i")	{
		spc_ep[t_l] <- "e";
		}
	return(paste(spc_ep,collapse=""));
	} else	{
	return("");
	}
}

#taxon_name <- "Redlichia chinensis - Kootenia gimmelfarbi" taxon_name <- "lower Fungochitina spinfera"
# taxon_name <- "Eospirifer"
# routine to separate species name from genus or genus (subgenus) name
divido_genus_names_from_species_names <- function(taxon_name)	{
#print(taxon_name);
k <- stringr::str_split(taxon_name,pattern="<")[[1]];
j <- stringr::str_split(k[1:(length(k)-1)],pattern=" ")[[1]];
#j <- stringr::str_split(taxon_name,pattern=" ")[[1]];
zone_detritus <- c("basal","lowermost","lower","middle","upper","uppermost","top");
j <- j[!j %in% zone_detritus];
j <- j[j!=""];
j <- j[!j %in% c("aff.","cf.","informal")]
genus_name <- "";
if (length(j)==2)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)==4)	{
	genus_name <- paste(j[1:2],collapse=" ");
	} else if (length(j)==3)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)>1)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)==1)	{
		jj <- strsplit(j,split="",fixed=TRUE)[[1]];
		if (jj[1]==toupper(jj[1]))
			genus_name <- j;
	}
return(genus_name)
}

# routine to separate subgenus & genus names from genus (subgenus) name
divido_subgenus_names_from_genus_names <- function(genus_name)	{
#j <- stringr::str_split(genus_name,pattern=" ")[[1]];
j <- strsplit(genus_name,split=" ")[[1]];
name_1 <- j[1];
if (length(j)==1)	{
	name_2 <- "";
	} else {
	jj <- strsplit(genus_name,split="")[[1]];
	if (sum(jj=="(")==1)	{
		name_2 <- paste(strsplit(j[2],"")[[1]][2:(length(strsplit(j[2],"")[[1]])-1)],collapse="");
		} else	{
		name_2 <- "";
		}
	}
return(c(name_1,name_2));
}

accersi_genus_subgenus_combinations <- function(genus_species_combo)	{
print(paste(genus_species_combo,sep=" "));
genus_name <- genus_species_combo[1];
species_name <- genus_species_combo[2];
if (species_name!="")	{
	gsg <- strsplit(genus_name,split=" ")[[1]];
	if (length(gsg)==2)	{
		g_s_g <- strsplit(gsg[2],split="")[[1]];
		if (g_s_g[1]=="(")	{
			gsg[2] <- paste(g_s_g[2:(length(g_s_g)-1)],collapse="");
			} else	{
			gsg[2] <- gsg[1];
			}
		} else {
		gsg <- c(gsg,gsg);
		}
	} else	{
	gsg <- c("","");
	}
#output <- matrix("",nrow=1,ncol=2);
#output[1,1] <- paste(gsg[1],species_name);
#output[1,2] <- paste(gsg[2],species_name);
genus_species <- paste(gsg[1],species_name);
subgenus_species <- paste(gsg[2],species_name);
#return(output)
return(data.frame(genus_species=as.character(genus_species),subgenus_species=as.character(subgenus_species)));
}

accersi_genus_subgenus_combinations_for_zones <- function(zone)	{
zone <- mundus_zone(zone);
genus_name <- transmogrify_full_zone_names_to_genus_names_only(zone);
output <- c("","");
if (length(genus_name)>0)	{
	species_epithet <- transmogrify_full_zone_names_to_species_names_only(zone);
	if (length(species_epithet)>0)	{
		gsg <- strsplit(genus_name,split=" ")[[1]];
		if (length(gsg)==2)	{
			g_s_g <- strsplit(gsg[2],split="")[[1]];
			if (g_s_g[1]=="(")	{
				gsg[2] <- paste(g_s_g[2:(length(g_s_g)-1)],collapse="");
				output <- c(paste(gsg[1],species_epithet),paste(gsg[2],species_epithet));
				} else	{
				output <- c(zone,zone);
				}
			} else if (length(gsg)==1)	{
			output <- c(zone,zone);
			}
		}
	}
return(output);
}

# for zones called things like "Zone G1" or "Ashgill Shelly Zone 1"
### returns vector that turns zone="Zone G-1 (Hintzeia celsaora)" into:
###		zone_name[1] = "Zone G1"
###		zone_name[2] = "G1"
aparecium_nontaxon_zone <- function(zone)	{
this_zone <- stringr::str_split(zone,pattern=" ")[[1]];
zone_title <- match("zone",tolower(this_zone));
zone_name <- array("",dim=2);
names(zone_name) <- c("non_taxon_zone","non_taxon_zone_label");
#zone_name <- "";
if (!is.na(zone_title))	{
	if (zone_title < length(this_zone))	{
		zone_name[1] <- this_zone[1];
		for (zt in 2:(zone_title+1))	{
			zone_name[1] <- paste(zone_name[1],this_zone[zt],sep=" ");
			}
		zone_name[1] <- gsub("-", "",zone_name[1]);
		zone_name[1] <- gsub("—", "",zone_name[1]);
		zone_name[2] <- this_zone[zone_title+1];
		zone_name[2] <- gsub("-", "",zone_name[2]);
		zone_name[2] <- gsub("—", "",zone_name[2]);
		}
	}
#zone_name <- data_frame("F")
return(zone_name);
}

# remove duplicate occurrences due to synonomies and/or lumping collections
divido_nonsynonymous_occurrences <- function(pbdb_finds)	{
reduced <- data.frame(collection_no=as.numeric(pbdb_finds$collection_no),accepted_name=as.character(pbdb_finds$accepted_name),stringsAsFactors = hell_no);
rownames(reduced) <- pbdb_finds$occurrence_no;
senior_only <- unique(reduced);
return(as.numeric(rownames(senior_only)));
}

# remove duplicate occurrences due to synonomies and/or lumping collections
remove_duplicate_occurrences_paleodb <- function(occurrences)	{
ttl_finds_cnd <- ttl_finds <- nrow(occurrences);
species <- sort(unique(occurrences$accepted_name));
otu <- length(species)

reduced_finds <- data.frame(collection_no=as.numeric(occurrences$collection_no),
							accepted_name=as.character(occurrences$accepted_name));
reduced_finds <- unique(reduced_finds);
keepers <- as.numeric(rownames(reduced_finds));
return(occurrences[keepers,]);
}

# update species names from newer opinions
transmogrify_accepted_species_name <- function(identified_name,accepted_genus)	{
#transmogrify_accepted_species_name <- function(name_to_fix)	{
#print(name_to_fix[1:2])
#accepted_genus <- name_to_fix[1];
species_epithet <- divido_species_epithets(identified_name);
return(paste(accepted_genus,species_epithet))
}

# routine to find and eliminate informal species designations that the PaleoDB considers species IDs
expello_indeterminate_species <- function(pbdb_finds)	{
taxon_name <- sort(unique(pbdb_finds$accepted_name));
species_epithets <- sapply(taxon_name,divido_species_epithets);
indet_species <- c("sp.");
indet_species <- c(indet_species,paste("sp.",LETTERS));
indet_species <- c(indet_species,paste("sp.",letters));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp. nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("indet.",i));
indet_species <- c(indet_species,"sensu");
species_epithets <- sapply(species_epithets,accersi_embedded_informal_names);
taxon_names <- taxon_name[!species_epithets %in% indet_species];
echino_species <- taxon_names <- unique(taxon_names);
taxon_names <- sapply(echino_species,echinoscrub);
taxon_names <- taxon_names[taxon_names!=""];
pbdb_finds <- pbdb_finds[pbdb_finds$accepted_name %in% taxon_names,];
return(pbdb_finds);
}

# get those hard to find informals!
accersi_embedded_informal_names <- function(species_epithets)	{
for (se in 1:length(species_epithets))	{
	if (length(strsplit(species_epithets[se],split=" ")[[1]])>1)	{
		if (sum(strsplit(species_epithets[se],split=" ")[[1]] %in% c(0:9,letters,LETTERS,"genus","sp.","indet."))>0)
			species_epithets[se] <- "sp. 1";
		}
	}
return(species_epithets)
}

### routines to standardize taxon names taxon_name="Sulcatospira? praecursor";
mundify_taxon_names <- function(taxon_name,keep_uncertainty=F)	{
# taxon_name: string giving species name
taxon_name <- gsub(" n\\. sp\\.\\?","",taxon_name);
taxon_name <- gsub(" n\\. sp\\.","",taxon_name);
taxon_name <- gsub("n\\. gen\\. ","",taxon_name);
taxon_name <- gsub(" n\\. subgen\\.","",taxon_name);
if (!keep_uncertainty)	{
	taxon_name <- gsub(" cf\\.","",taxon_name);
	taxon_name <- gsub(" aff\\.","",taxon_name);
	taxon_name <- gsub(" ex gr\\.","",taxon_name);
	taxon_name <- gsub("cf\\. ","",taxon_name);
	taxon_name <- gsub("aff\\. ","",taxon_name);
	taxon_name <- gsub(" ex gr\\.","",taxon_name);
	taxon_name_broken <- stringr::str_split(taxon_name,pattern=" ")[[1]];
	if (length(taxon_name_broken)>1 && (taxon_name_broken[2]=="gr." && taxon_name_broken[1]==taxon_name_broken[3]))
		taxon_name <- paste(taxon_name_broken[3:length(taxon_name_broken)],collapse=" ");
	taxon_name <- gsub(" gr\\.","",taxon_name);
	}
taxon_name <- gsub(" informal","",taxon_name);
taxon_name <- gsub("sensu lato","",taxon_name);
taxon_name <- gsub(" sensu lato","",taxon_name);
taxon_name <- gsub("\"", "",taxon_name);
taxon_name <- gsub(" \\?" ,"",taxon_name);
taxon_name <- gsub("\\? " ," ",taxon_name);
taxon_name <- gsub("\\?" ,"",taxon_name);
taxon_name <- gsub("  " ," ",taxon_name);
molecularized_name <- strsplit(taxon_name,split="")[[1]];
pbdb_notes <- (1:length(molecularized_name))[molecularized_name %in% c("<",">")];
if (length(pbdb_notes)>0)	{
	molecularized_name <- molecularized_name[!(1:length(molecularized_name)) %in% pbdb_notes[1]:pbdb_notes[2]];
	while (molecularized_name[length(molecularized_name)]==" ")
		molecularized_name <- molecularized_name[1:(length(molecularized_name)-1)];
	while (molecularized_name[1]==" ")
		molecularized_name <- molecularized_name[2:length(molecularized_name)];
	taxon_name <- paste(molecularized_name,collapse="");
	}
taxon_name <- paste(strsplit(taxon_name,split=" ")[[1]][!strsplit(taxon_name,split=" ")[[1]] %in% paste(LETTERS,".",sep="")],collapse=" ");
#paste(LETTERS,".",sep="");
return(taxon_name)
}

reveal_uncertain_taxonomic_assignments_for_collection <- function(desired_finds)	{
noccr <- nrow(desired_finds);
taxon_name <- desired_finds$identified_name;
print("Separating uncertain species assignments.")
flags1 <- pbapply::pbsapply(taxon_name,revelare_uncertain_species_assignments);
flags3 <- flags2 <- rep("",noccr);
taxon_name <- desired_finds$identified_name[desired_finds$accepted_rank %in% c("genus","subgenus")];
if (length(taxon_name)>0)	{
	flags2[desired_finds$accepted_rank %in% c("genus","subgenus")] <- sapply(taxon_name,revelare_uncertain_genus_assignments);
	double <- (1:noccr)[flags1!=""][(1:noccr)[flags1!=""] %in% (1:noccr)[flags2!=""]];
	flags3[flags1!=""] <- flags1[flags1!=""];
	flags3[flags2!=""] <- flags2[flags2!=""];
	flags3[double] <- paste(unique(flags2[flags2!=""]),unique(flags1[flags1!=""]),sep=", ");
	desired_finds$flags <- simplify2array(flags3);
	} else	{
	desired_finds$flags <- flags1;
	}
return(desired_finds);
}

### routine to identify uncertain taxon assignment type
revelare_uncertain_species_assignments <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- "";
cleaned_name <- mundify_taxon_names(taxon_name);
taxon_name <- gsub("n. gen. ","",taxon_name);
taxon_name <- gsub("n. sp. ","",taxon_name);
taxon_name <- gsub("n. subgen. ","",taxon_name);
species_epithet <- mundify_taxon_names(taxon_name=divido_species_epithets(cleaned_name));
genus_name <- divido_genus_names_from_species_names(taxon_name);
subgenus_name <- divido_subgenus_names_from_genus_names(genus_name)[2];
if (subgenus_name!="")
	subgenus_name <- paste("(",subgenus_name,")",sep="");
modified_taxon_name <- gsub(" ex gr\\."," ex_gr\\.",taxon_name);
taxon_components <- simplify2array(strsplit(modified_taxon_name," ")[[1]]);
t_c <- 1:length(taxon_components);
if (!is.na(match("?",taxon_components)) && max(t_c[taxon_components %in% "?"])==length(taxon_components))	{
	flags <- "uncertain species";
	} else if (sum(taxon_components %in% uncertains) > 0)	{
	btc <- t_c[taxon_components %in% uncertains];
	if (!is.na(match(subgenus_name,taxon_components)) && (match(subgenus_name,taxon_components)+1) %in% btc)	{
		flags <- "uncertain species";
		} else if (subgenus_name=="" && (match(genus_name,taxon_components)+1) %in% btc)	{
		flags <- "uncertain species";
		}
	}
return(flags);
}

# remember to use only ids of species or subspecies!!!! taxon_names <- taxon_names[40]
repair_misentered_uncertain_species <- function(taxon_name)	{
taxon_name_split <- strsplit(taxon_name," ")[[1]]
if (taxon_name_split[1] %in% c("cf.","aff."))	{
#	if (is.subgenus(paste(taxon_name_split[2:length(taxon_name_split)],collapse=" ")))	{
	if (gsub("\\(","",taxon_name)!=taxon_name)	{
		if (taxon_name_split[4] %in% c("n.","?","cf.","aff."))	{
			return(taxon_name);
			} else	{
			new_order <- c(2,3,1,4:length(taxon_name_split));
			return(paste(taxon_name_split[new_order],collapse=" "));
			}
		} else	{
		if (taxon_name_split[3] %in% c("n.","?","cf.","aff."))	{
			return(taxon_name);
			} else	{
			new_order <- c(2,1,3:length(taxon_name_split));
			return(paste(taxon_name_split[new_order],collapse = " "));
			}
#			paste(taxon_name_split[2],taxon_name_split[1],taxon_name_split[3:length(taxon_name_split)],collapse=" ");
		}
	} else	{
	return(paste(taxon_name_split,collapse=" "));
	}
}

# find taxon labels that are informal names
revelare_informal_taxa <- function(taxon_name)	{
molecularized_name <- strsplit(taxon_name,split=" ")[[1]];
atomized_name <- strsplit(taxon_name,split="")[[1]];

if ((sum(atomized_name %in% ".")+sum(atomized_name %in% as.character(0:9)))>0)	{
	return(T);
	} else if (sum(molecularized_name %in% c("sp.","indet.","informal"))>0 || length(molecularized_name)>4 || (length(molecularized_name)==4 && strsplit(molecularized_name[2],split="")[[1]][1]!="("))	{
	return(T);
	} else	{
	return(F);
	}
}

revelare_uncertain_genus_assignments <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- "";
cleaned_name <- mundify_taxon_names(taxon_name);
taxon_name <- gsub("n. gen. ","",taxon_name);
taxon_name <- gsub("n. sp. ","",taxon_name);
taxon_name <- gsub("n. subgen. ","",taxon_name);
#species_epithet <- mundify_taxon_names(taxon_name=divido_species_epithets(cleaned_name));
whole_genus_name <- divido_genus_names_from_species_names(taxon_name=cleaned_name);
gen_subgen <- divido_subgenus_names_from_genus_names(genus_name=whole_genus_name);
genus_name <- mundify_taxon_names(taxon_name=gen_subgen[1]);
genus_name_check <- divido_subgenus_names_from_genus_names(genus_name=divido_genus_names_from_species_names(taxon_name))[1];
subgenus_name <- mundify_taxon_names(taxon_name=gen_subgen[2]);
if (subgenus_name!="")
	subgenus_name <- paste("(",subgenus_name,")",sep="");
modified_taxon_name <- gsub(" ex gr\\."," ex_gr\\.",taxon_name);
taxon_components <- simplify2array(strsplit(modified_taxon_name," ")[[1]]);

#if (genus_name != genus_name_check || taxon_components[2]=="?" || sum(taxon_components[1] %in% uncertains)==1)	{
if (genus_name != genus_name_check || taxon_components[1] %in% uncertains)	{
	genus_name <- genus_name_check
	flags <- "uncertain genus";
	} else if (sum(taxon_components %in% "?)")==1)	{
	flags <- "uncertain genus";
	}
return(flags);
}

### routine to identify uncertain taxon assignment type; taxon_name <- "Troostella ? sp."
identify_taxonomic_uncertainty <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- revelare_uncertain_genus_assignments(taxon_name);
flags <- c(flags,revelare_uncertain_species_assignments(taxon_name));
#paste(flags,collapse=", ");
if (sum(flags!="")==2) {
	return(paste(flags,collapse=", "));
	} else if (sum(flags=="")==2)	{
	return("");
	} else	{
	return(flags[flags!=""]);
	}
}

# routine to separate species name from whole name
divido_species_epithets <- function(taxon_name)	{
#split species name (or subspecies names) off of genus species combo
#print(taxon_name);
j <- simplify2array(strsplit(taxon_name," ")[[1]]);
j <- j[j!=""];
species_name <- "";
wordy <- length(j);
if (wordy==2)	{
	if (strsplit(j[2],split="")[[1]][1]!="(")
		species_name <- j[2];
	} else if (wordy > 1) {
		if (simplify2array(strsplit(j[2],"")[[1]])[1]=="(")	{
		# subgenus
		species_name <- paste(j[3:wordy],collapse=" ");
		} else	{
		species_name <- paste(j[2:wordy],collapse=" ");
		}
	} else if (wordy==1)	{
		if(strsplit(j,split="")[[1]][1] == tolower(strsplit(j,split="")[[1]][1]))
			species_name <- j;
	}
return(species_name);
}

lump_subspecies_to_species <- function(species_name)	{
species_names <- strsplit(species_name," ")[[1]];
return(paste(species_names[1:(length(species_names)-1)],collapse=" "));
}

# species_name <- "Eostropheodonta chilcaensis parvula";
elevate_subspecies_to_species <- function(species_name)	{
species_names <- strsplit(species_name," ")[[1]];
ttl_names <- length(species_names);
return(paste(species_names[(1:length(species_names))[!(1:length(species_names)) %in% (ttl_names-1)]],collapse=" "));
}

#species_name <- "Liobolina (Guilinaspis) intermedia" smithi"
elevate_subspecies_to_species_whole <- function(species_name)	{
species_names <- strsplit(species_name," ")[[1]];
ttl_names <- length(species_names);
if (ttl_names==4)	{
	return(paste(species_names[c(1:2,4)],collapse=" "))
	} else if (ttl_names==3)	{
	second_name <- strsplit(species_names[2],"")[[1]];
	if (second_name[1]=="(")	{
		return(species_name);
		}	else	{
		return(paste(species_names[c(1,3)],collapse=" "));
		}
	} else	{
	return(species_name);
	}
#return(paste(species_names[(1:length(species_names))[!(1:length(species_names)) %in% (ttl_names-1)]],collapse=" "));
}

#species_name <- "Guilinaspis subcylindrica subcylindrica";
#species_name <- "Guilinaspis subcylindrica brevis"
elevate_type_subspecies_to_species <- function(species_name)	{
species_names <- strsplit(species_name," ")[[1]];
ttl_names <- length(species_names);
if (ttl_names>2 && species_names[ttl_names]==species_names[ttl_names-1])
	species_name <- paste(species_names[1:(ttl_names-1)],collapse=" ")
return(species_name);
}

# taxon_name <- "Eostropheodonta (Eostropheodonta)"
is.subgenus <- function(taxon_name)	{
if (length(simplify2array(strsplit(taxon_name," ")[[1]]))==2)	{
	second_name <- simplify2array(strsplit(simplify2array(strsplit(taxon_name," ")[[1]])[2],"")[[1]]);
	if (second_name[1]=="\\(" || second_name[1]=="(")	{
		return(T);
		} else	{
		return(F);
		}
	} else	{
	return(F);
	}
}

is.eponymous.subgenus <- function(taxon_name)	{
if (is.subgenus(taxon_name))	{
	names <- strsplit(taxon_name," ")[[1]];
	names[2] <- gsub("\\(","",names[2]);
	names[2] <- gsub("\\)","",names[2]);
	if (names[1]==names[2])	{
		return(T);
		} else	{
		return(F);
		}
	} else	{
	return(F);
	}
}

is.species <- function(taxon_name)	{
if (length(simplify2array(strsplit(taxon_name," ")[[1]]))==2)	{
	second_name <- simplify2array(strsplit(simplify2array(strsplit(taxon_name," ")[[1]])[2],"")[[1]]);
	if (second_name[1]=="\\(" || second_name[1]=="(")	{
		return(F);
		} else	{
		return(T);
		}
	} else	{
	return(F);
	}
}

# taxon_name <- "Eostropheodonta (Eostropheodonta) chilcaensis chilcaensis"
is.subspecies <- function(taxon_name)	{
if (length(simplify2array(strsplit(taxon_name," ")[[1]]))>2 && length(simplify2array(strsplit(taxon_name," ")[[1]]))<5)	{
	if (length(simplify2array(strsplit(taxon_name," ")[[1]]))==4)	{
		second_name <- simplify2array(strsplit(simplify2array(strsplit(taxon_name," ")[[1]])[2],"")[[1]]);
		if (second_name[1]=="\\(" || second_name[1]=="(")	{
			taxon_name <- gsub(" \\(","_(",taxon_name);
			}
		}
	second_name <- simplify2array(strsplit(simplify2array(strsplit(taxon_name," ")[[1]])[2],"")[[1]]);
	third_name <- simplify2array(strsplit(simplify2array(strsplit(taxon_name," ")[[1]])[3],"")[[1]]);
	if (second_name[1] %in% letters && third_name[1] %in% letters)	{
		return(T);
		} else	{
		return(F);
		}
	} else	{
	return(F);
	}
}

is.eponymous.subspecies	<- function(taxon_name)	{
if (is.subspecies(taxon_name))	{
	all_names <- strsplit(taxon_name," ")[[1]];
	ttl_names <- length(all_names);
	if (all_names[ttl_names]==all_names[ttl_names-1])	{
		return(T);
		} else	{
		return(F);
		}
	} else	{
	return(F);
	}
}

# pentameral funkiness begone!
echinoscrub <- function(echino_species)	{
echinobabble <- c("columnals","debris","stem","stems","holdfast","holdfasts","miscellanea","miscellaneus","ossicle","ossicles","plate","plates");
j <- simplify2array(strsplit(echino_species," ")[[1]]);
if (sum(tolower(j) %in% echinobabble)>0)	{
	return("");
	} else	{
	return(echino_species);
	}
}

# get summaries of higher taxa from pbdb finds
suprageneric_taxon_summaries <- function(higher_taxa,valid_genera,pbdb_finds)	{
suprageneric_summaries <- data.frame(superfamily=as.character(),
									 family=as.character(),
									 subfamily=as.character(),
									 families=as.numeric(),genera=as.numeric(),species=as.numeric(),
									 n_occs=as.numeric(),stringsAsFactors = hell_no);
dummy <- data.frame(superfamily=as.character(""),family=as.character(""),
					subfamily=as.character(""),families=as.numeric(0),
					genera=as.numeric(0),species=as.numeric(0),
					n_occs=as.numeric(0),stringsAsFactors = hell_no);
for (sfm in 1:length(higher_taxa))	{
	this_superfamily <- subset(valid_genera,valid_genera$superfamily==higher_taxa[sfm]);
	this_superfamily_finds <- subset(pbdb_finds,pbdb_finds$parent %in% valid_genera$taxon_name);
	these_families <- unique(this_superfamily$family);
	these_families <- these_families[these_families!="NO_FAMILY_SPECIFIED"];
	dummy$superfamily <- higher_taxa[sfm];
	dummy$family <- dummy$subfamily <- "";
	dummy$families <- length(unique(this_superfamily_finds$family));
	dummy$genera <- length(unique(this_superfamily_finds$parent));
	dummy$species <- length(unique(this_superfamily_finds$accepted_name));
	dummy$n_occs <- nrow(this_superfamily_finds);
	suprageneric_summaries <- rbind(suprageneric_summaries,dummy);
	for (tf in 1:length(these_families))	{
		dummy$family <- these_families[tf];
		dummy$subfamily <- "";
		this_family <- subset(this_superfamily,this_superfamily$family==these_families[tf]);
		this_family_finds <- subset(pbdb_finds,pbdb_finds$parent %in% this_family$taxon_name);
		these_subfamilies <- sort(unique(this_family$subfamily[this_family$subfamily!=""]));
		dummy$families <- 0;
		dummy$genera <- length(unique(this_family_finds$parent));
		dummy$species <- length(unique(this_family_finds$accepted_name));
		dummy$n_occs <- nrow(this_family_finds);
		suprageneric_summaries <- rbind(suprageneric_summaries,dummy);
		ts <- 0;
		while (ts < length(these_subfamilies))	{
			ts <- ts+1;
			dummy$subfamily <- these_subfamilies[ts];
			this_subfamily <- subset(this_family,this_family$parent_name==these_subfamilies[ts]);
			this_subfamily_finds <- subset(pbdb_finds,pbdb_finds$parent %in% this_subfamily$taxon_name);
			dummy$genera <- length(unique(this_subfamily_finds$parent));
			dummy$species <- length(unique(this_subfamily_finds$accepted_name));
			dummy$n_occs <- nrow(this_subfamily_finds);
			suprageneric_summaries <- rbind(suprageneric_summaries,dummy);
			}
		}
	}
return(suprageneric_summaries);
}

					##### REFERENCE DOWNLOADING ROUTINES #######
# look up references for data type
accersi_references_ris <- function(taxon="Life",ref_types=c("auth","class","ops","occs","specs","colls"),output_type="ris")	{

# taxon: taxonomic group
# ref_types: occs / colls / ops [opinions] / class [classification]
biblio <- ref_types[1];
if (length(ref_types)>1)	{
	for (b in 2: length(ref_types))	biblio <- paste(biblio,",",ref_types[b],sep="");
	}
taxa <- taxon[1];
if (length(taxon)>1)	{
	for (s in 2:length(taxon))	taxa <- paste(taxa,",",taxon[s],sep="");
	}
#http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=refattr,entname&limit=all",sep = "")
#http <- paste("https://www.paleobiodb.org/data1.2/taxa/refs.ris?datainfo&rowcount&base_name=",taxa,"&select=",biblio,sep="");
http <- paste("https://www.paleobiodb.org/data1.2/taxa/refs.ris?base_name=",taxa,"&select=",biblio,sep="");
fetch <- RCurl::getURL(http);
web_text <- gsub("///","",fetch);
fetch <- sapply(web_text,mundify_web_text_pizzaz);
write(fetch,file="pbdb_data.ris");
ris <- simplify2array(strsplit(fetch,split="\r\n")[1]);
write.table(x=ris,file="pbdb_data.ris",quote=FALSE,row.names = FALSE,col.names = FALSE);

#ris <- simplify2array(read.table(http, sep='\n', header=F));
#ris <- read.table(text = fetch, header = FALSE, stringsAsFactors=FALSE, sep="");
#ris <- gsub("///","",ris);
#ris <- scan(file="pbdb_data.ris",what=character(),sep="\n");
x <- c("reference_no","author1init","author1last","author2init","author2last","otherauthors","pubyr","reftitle","pubtitle","editors","pubvol","pubno","firstpage","lastpage","publication_type","basis","language","doi","comments");
#references <- data.frame(array("",dim=c(0,length(x))));
new_ref <- data.frame(array("",dim=c(1,length(x))),stringsAsFactors = TRUE);
colnames(new_ref) <- x;
pubtypes <- data.frame(c("abstract","book","book chapter","general","journal article","thesis","unpublished"));
rownames(pubtypes) <- c("ABST","BOOK","CHAP","GEN","JOUR","THES","UNPD");
ttl_refs <- 0;
references <- c();
i <- 0;
#for (i in 1:length(ris))  {
while (i<length(ris))  {
	i <- i+1;
	if (i%%500==0)	print(paste("Line",i))
#	new_ref;
	ris[i] <- mundify_web_text_pizzaz(web_text = ris[i]);
	j <- strsplit(ris[i],split="",fixed=TRUE)[[1]];
	while(j[length(j)]==" ")	j <- j[1:(length(j)-1)]
	word1 <- paste(j[1],j[2],sep="");
	word2 <- c();
	m <- 2+match("-",j);
	#TY	publication_type
	if (word1=="TY")	{
#		references <- rbind(references,new_ref);
		if (ttl_refs > 0)	references <- rbind(references,new_ref);
		ttl_refs <- ttl_refs+1;
		for (a in 1:length(x))	new_ref[a] <- "";
		ttl_authors <- 0;
		for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
		rtype <- match(word2,rownames(pubtypes));
		new_ref$publication_type <- as.character(pubtypes[rtype,]);
	#M3	thesis_type
		} else if (word1=="M3")	{
		for (k in m:length(j))	new_ref$publication_type <- paste(new_ref$publication_type,j[k],sep="");
	#ID	reference_no
		} else if (word1=="ID")	{
		m <- 1+match(":",j);
		for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
		new_ref$reference_no <- word2;
	#AU	author
		} else if (word1=="AU")	{
		ttl_authors <- ttl_authors+1;
		k <- 2+match("-",j);
		word2 <- c();
		if (ttl_authors==1)	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			new_ref$author1last <- word2;
			m <- k+1;
			word2 <- c();
			if (m<length(j)){
				for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
				new_ref$author1init <- word2;
				}
			} else 		if (ttl_authors==2)	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			new_ref$author2last <- word2;
			m <- k+1;
			word2 <- c();
			if (m<length(j))	{
				for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
				new_ref$author2init <- word2;
				}
			}	else	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			m<- k+1;
			word3 <- c();
			for (k in m:length(j))	word3 <- paste(word3,j[k],sep="");
			if (ttl_authors==3)	{
				new_ref$otherauthors <- paste(word3,word2);
				}	else	{
				new_ref$otherauthors <- paste(new_ref$otherauthors,", ",word3," ",word2,sep="");
				}
			}
	#A2	editors
		} else if (word1=="A2")	{
		for (k in m:length(j))	new_ref$editors <- paste(new_ref$editors,j[k],sep="");
	#PY	pubyr
		} else if (word1=="PY")	{
		for (k in m:(m+3))	new_ref$pubyr <- paste(new_ref$pubyr,j[k],sep="");
	#TI	reftitle
		} else if (word1=="TI")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$reftitle <- paste(new_ref$reftitle,j[k],sep="");
	#T2	pubtitle
		} else if (word1=="T2")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubtitle <- paste(new_ref$pubtitle,j[k],sep="");
	#VL	pubvol
		} else if (word1=="VL")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubvol <- paste(new_ref$pubvol,j[k],sep="");
	#IS	pubno
		} else if (word1=="IS")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubno <- paste(new_ref$pubno,j[k],sep="");
	#SP	pages
		} else if (word1=="SP")	{
		n <- match("-",j[m:length(j)]);
		if (!is.na(n))	{
			for (k in m:(m+n-2))		new_ref$firstpage <- paste(new_ref$firstpage,j[k],sep="");
			for (k in (m+n):length(j))	new_ref$lastpage <- paste(new_ref$lastpage,j[k],sep="");
			}	else	{
			if (m<=length(j))
				for (k in m:length(j))		new_ref$firstpage <- new_ref$lastpage <- paste(new_ref$firstpage,j[k],sep="");
			}
	#LA	language
		} else if (word1=="LA")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$language <- paste(new_ref$language,j[k],sep="");
	#DO	doi
		} else if (word1=="DO")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$doi <- paste(new_ref$doi,j[k],sep="");
	#N1	comments
		} else if (word1=="N1")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$comments <- paste(new_ref$comments,j[k],sep="");
		}
#	new_ref;
	}
#references <- rbind(references,new_ref);
write.table(x=references[order(as.numeric(references$reference_no)),],file="PaleoDB_Refs.txt",quote=FALSE,row.names = FALSE,col.names = TRUE,sep="\t");
output <- list(ris,references);
names(output) <- c("ris","references_table");
return(output);
}

accersi_references_csv <- function(taxon="Life",onset="Hadean",end="Phanerozoic",ref_types=c("auth","class","ops","occs","specs","colls"))	{
http <- paste("https://www.paleobiodb.org/data1.2/occs/refs.csv?base_name=",taxon,"&interval=",onset,",",end,"&select=",paste(ref_types,collapse=","),"&private&show=crmod",sep="");
library <- read.csv(http,header=T,stringsAsFactors=hell_no,encoding="UTF-8");

fetch <- RCurl::getURL(http);
web_text <- gsub("///","",fetch);
fetch <- sapply(web_text,mundify_web_text_pizzaz);
all_refs <- utils::read.csv(text = fetch, header = T, stringsAsFactors=hell_no,encoding="UTF-8");
#write(fetch,file="pbdb_data.ris");
#ris <- simplify2array(strsplit(fetch,split="\r\n")[1]);
write.csv(fetch,file="pbdb_data.ris",quote=FALSE,row.names = FALSE,col.names = FALSE);

#ris <- simplify2array(read.table(http, sep='\n', header=F));
#ris <- read.table(text = fetch, header = FALSE, stringsAsFactors=FALSE, sep="");
#ris <- gsub("///","",ris);
#ris <- scan(file="pbdb_data.ris",what=character(),sep="\n");
x <- c("reference_no","author1init","author1last","author2init","author2last","otherauthors","pubyr","reftitle","pubtitle","editors","pubvol","pubno","firstpage","lastpage","publication_type","basis","language","doi","comments");
#references <- data.frame(array("",dim=c(0,length(x))));
new_ref <- data.frame(array("",dim=c(1,length(x))),stringsAsFactors = TRUE);
colnames(new_ref) <- x;
pubtypes <- data.frame(c("abstract","book","book chapter","general","journal article","thesis","unpublished"));
rownames(pubtypes) <- c("ABST","BOOK","CHAP","GEN","JOUR","THES","UNPD");
ttl_refs <- 0;
references <- c();
i <- 0;
#for (i in 1:length(ris))  {
while (i<length(ris))  {
	i <- i+1;
	if (i%%500==0)	print(paste("Line",i))
#	new_ref;
	ris[i] <- mundify_web_text_pizzaz(web_text = ris[i]);
	j <- strsplit(ris[i],split="",fixed=TRUE)[[1]];
	while(j[length(j)]==" ")	j <- j[1:(length(j)-1)]
	word1 <- paste(j[1],j[2],sep="");
	word2 <- c();
	m <- 2+match("-",j);
	#TY	publication_type
	if (word1=="TY")	{
#		references <- rbind(references,new_ref);
		if (ttl_refs > 0)	references <- rbind(references,new_ref);
		ttl_refs <- ttl_refs+1;
		for (a in 1:length(x))	new_ref[a] <- "";
		ttl_authors <- 0;
		for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
		rtype <- match(word2,rownames(pubtypes));
		new_ref$publication_type <- as.character(pubtypes[rtype,]);
	#M3	thesis_type
		} else if (word1=="M3")	{
		for (k in m:length(j))	new_ref$publication_type <- paste(new_ref$publication_type,j[k],sep="");
	#ID	reference_no
		} else if (word1=="ID")	{
		m <- 1+match(":",j);
		for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
		new_ref$reference_no <- word2;
	#AU	author
		} else if (word1=="AU")	{
		ttl_authors <- ttl_authors+1;
		k <- 2+match("-",j);
		word2 <- c();
		if (ttl_authors==1)	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			new_ref$author1last <- word2;
			m <- k+1;
			word2 <- c();
			if (m<length(j)){
				for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
				new_ref$author1init <- word2;
				}
			} else 		if (ttl_authors==2)	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			new_ref$author2last <- word2;
			m <- k+1;
			word2 <- c();
			if (m<length(j))	{
				for (k in m:length(j))	word2 <- paste(word2,j[k],sep="");
				new_ref$author2init <- word2;
				}
			}	else	{
			while (j[k]!="," && k<=length(j))	{
				word2 <- paste(word2,j[k],sep="");
				k <- k+1;
				}
			m<- k+1;
			word3 <- c();
			for (k in m:length(j))	word3 <- paste(word3,j[k],sep="");
			if (ttl_authors==3)	{
				new_ref$otherauthors <- paste(word3,word2);
				}	else	{
				new_ref$otherauthors <- paste(new_ref$otherauthors,", ",word3," ",word2,sep="");
				}
			}
	#A2	editors
		} else if (word1=="A2")	{
		for (k in m:length(j))	new_ref$editors <- paste(new_ref$editors,j[k],sep="");
	#PY	pubyr
		} else if (word1=="PY")	{
		for (k in m:(m+3))	new_ref$pubyr <- paste(new_ref$pubyr,j[k],sep="");
	#TI	reftitle
		} else if (word1=="TI")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$reftitle <- paste(new_ref$reftitle,j[k],sep="");
	#T2	pubtitle
		} else if (word1=="T2")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubtitle <- paste(new_ref$pubtitle,j[k],sep="");
	#VL	pubvol
		} else if (word1=="VL")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubvol <- paste(new_ref$pubvol,j[k],sep="");
	#IS	pubno
		} else if (word1=="IS")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$pubno <- paste(new_ref$pubno,j[k],sep="");
	#SP	pages
		} else if (word1=="SP")	{
		n <- match("-",j[m:length(j)]);
		if (!is.na(n))	{
			for (k in m:(m+n-2))		new_ref$firstpage <- paste(new_ref$firstpage,j[k],sep="");
			for (k in (m+n):length(j))	new_ref$lastpage <- paste(new_ref$lastpage,j[k],sep="");
			}	else	{
			if (m<=length(j))
				for (k in m:length(j))		new_ref$firstpage <- new_ref$lastpage <- paste(new_ref$firstpage,j[k],sep="");
			}
	#LA	language
		} else if (word1=="LA")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$language <- paste(new_ref$language,j[k],sep="");
	#DO	doi
		} else if (word1=="DO")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$doi <- paste(new_ref$doi,j[k],sep="");
	#N1	comments
		} else if (word1=="N1")	{
		if (m<=length(j))	for (k in m:length(j))	new_ref$comments <- paste(new_ref$comments,j[k],sep="");
		}
#	new_ref;
	}
#references <- rbind(references,new_ref);
write.table(x=references[order(as.numeric(references$reference_no)),],file="PaleoDB_Refs.txt",quote=FALSE,row.names = FALSE,col.names = TRUE,sep="\t");
output <- list(ris,references);
names(output) <- c("ris","references_table");
return(output);
}

accersi_reference_by_number <- function(reference_no)	{
http <- paste("https://www.paleobiodb.org/data1.2/refs/single.csv?id=",reference_no,"&show=both,crmod",sep="");

reference <- read.csv(http,header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
reference <- clear_na_from_matrix(reference,"");
nxt_ref <- reference_no;
while (reference$created=="0000-00-00 00:00:00")	{
	nxt_ref <- nxt_ref+1;
	http <- paste("https://www.paleobiodb.org/data1.2/refs/single.csv?id=",nxt_ref,"&show=both,crmod",sep="");
	kluge <- read.csv(http,header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
	if (nrow(kluge)==1)	reference$created <- kluge$created
	}
if (reference$modified=="0000-00-00 00:00:00")	reference$modified <- reference$created;

reference <- put_pbdb_dataframes_into_proper_type(reference);
reference <- clear_na_from_matrix(reference,0);
return(reference);
}

accersi_references_by_numbers <- function(reference_nos)	{
#http <- paste("https://www.paleobiodb.org/data1.2/refs/list.txt?id=",reference_no,"&show=both",sep="");
#reference_no <- paste(reference_nos,collapse=",");
reference_no <- reference_nos;
pbdb_references <- data.frame(t(pbapply::pbsapply(reference_no,accersi_reference_by_number)),stringsAsFactors = hell_no);
for (cc in 1:ncol(pbdb_references))	pbdb_references[,cc] <- as.vector(unlist(pbdb_references[,cc]));
#http <- paste("https://www.paleobiodb.org/data1.2/refs/single.txt?id=",reference_no,"&show=both",sep="");
#pbdb_references <- read.csv(http, header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
pbdb_references <- clear_na_from_matrix(pbdb_references,"");
pbdb_references <- put_pbdb_dataframes_into_proper_type(pbdb_data=pbdb_references);
pbdb_references <- clear_na_from_matrix(pbdb_references,0);
return(pbdb_references);
}

accersi_references_by_author <- function(author)	{
author <- gsub(" ","%20",author);
http <- paste("https://www.paleobiodb.org/data1.2/refs/list.csv?ref_author=",author,"&show=both",sep="");
pbdb_references <- read.csv(http, header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
if (pbdb_references$reference_no[1]=="THIS REQUEST RETURNED NO RECORDS")
	pbdb_references <- pbdb_references[integer(),];
pbdb_references <- clear_na_from_matrix(pbdb_references,"");
pbdb_references <- put_pbdb_dataframes_into_proper_type(pbdb_references);
pbdb_references <- clear_na_from_matrix(pbdb_references,0);
return(pbdb_references);
}

accersi_references_for_taxon <- function(taxon)	{
#http <- paste("https://www.paleobiodb.org/data1.2/refs/list.txt?id=",reference_no,"&show=both",sep="");
#http <- paste("https://www.paleobiodb.org/data1.2/refs/single.txt?id=",reference_no,"&show=both",sep="");
#http <- paste("https://paleobiodb.org/data1.2/taxa/refs.csv?datainfo&rowcount&base_name=",taxon,"&select=auth,class,ops,occs,specs,colls",sep="")
#http <- paste("https://paleobiodb.org/data1.2/taxa/refs.csv?datainfo&rowcount&base_name=",taxon,"&select=auth,class,ops,occs,specs,colls",sep="");
taxon <- paste(taxon, collapse = ",");
taxon <- gsub(" ","%20",taxon);
http <- paste("https://paleobiodb.org/data1.2/taxa/refs.csv?base_name=",taxon,"&select=auth,class,ops,occs,specs,colls",sep="");
#fetch <- RCurl::getURL(http);
#web_text <- gsub("///","",fetch);
#fetch <- sapply(web_text,mundify_web_text_pizzaz);
pbdb_references <- read.csv(http, header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
pbdb_references <- clear_na_from_matrix(pbdb_references,"");
pbdb_references <- put_pbdb_dataframes_into_proper_type(pbdb_references);
pbdb_references <- clear_na_from_matrix(pbdb_references,0);
pbdb_references <- pbdb_references[order(pbdb_references$reference_no),];
#reference <- read.table(http, header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
return(pbdb_references);
}

#reference <- references[370,]
turn_reference_to_citation <- function(reference,style="et al.")	{
# if style = "Full" then all authors listed regardless of number
# if style = "et al." then 3+ authors get author1last et al.
reference$author1last <- gsub(", Jr.","",reference$author1last);
reference$author1last <- gsub(", III","",reference$author1last);
reference$author1last <- gsub(",Jr.","",reference$author1last);
reference$author2last <- gsub(", Jr.","",reference$author2last);
reference$author2last <- gsub(", III","",reference$author2last);
reference$author2last <- gsub(",Jr.","",reference$author2last);

if (reference$otherauthors!="" && tolower(style)=="full")	{
	reference$otherauthors <- gsub(", Jr.","",reference$otherauthors);
	reference$otherauthors <- gsub(", III","",reference$otherauthors);
	reference$otherauthors <- gsub(",Jr.","",reference$otherauthors);
	white_walkers <- strsplit(reference$otherauthors,",")[[1]]
	white_walkers <- white_walkers[!white_walkers %in% c("Jr"," Jr","III"," III")];
	author_list <- c(reference$author1last,reference$author2last);
	for (ww in 1:length(white_walkers))	{
		this_author <- strsplit(white_walkers[ww],". ",)[[1]];
		if (length(this_author)==1)	this_author <- strsplit(white_walkers[ww],"\\.",)[[1]];
		this_author_l <- this_author[length(this_author)];
		this_author_l <- strsplit(this_author_l,"")[[1]];
		while (this_author_l[1]==" " && length(this_author)>1)	this_author_l <- this_author_l[2:length(this_author_l)];
		while (this_author_l[length(this_author_l)]==" " && length(this_author)>1)	this_author_l <- this_author_l[1:(length(this_author_l)-1)];
		this_author_l <- paste(this_author_l,collapse="");
		author_list <- c(author_list,this_author_l);
		}
	author_list <- author_list[!author_list %in% c(""," ")];
	ttl_authors <- length(author_list);
	if (ttl_authors==2)	{
		citation <- paste(author_list,collapse=" and ");
		} else	{
		citation <- c();
		for (al in 1:(length(author_list)-2))	{
			citation <- paste(citation,paste(author_list[al],", ",sep=""),sep="")
			}
		citation <- paste(citation,paste(author_list[(ttl_authors-1):ttl_authors],collapse=" and "),sep="");
		}
	} else if (reference$otherauthors!="")	{
	citation <- paste(reference$author1last,"et al.",reference$pubyr);
	} else if (reference$author2last!="")	{
	citation <- paste(reference$author1last,"and",reference$author2last,reference$pubyr);
	} else	{
	citation <- paste(reference$author1last,reference$pubyr);
	}
return(citation);
}

accersi_all_pbdb_references <- function(first_ref=1,track_progress=T)	{
# first_ref: where to start. If you already have the first 5000 & don't care about
#	edits, then start with first_ref=5001.
duds <- 0;
reference_no <- first_ref;
references <- c();
while (duds < 10)	{
	if (track_progress && reference_no%%100==0)	print(paste(reference_no,date()));
	http <- paste("https://www.paleobiodb.org/data1.2/refs/single.csv?id=",reference_no,"&show=both,crmod",sep="");
	reference <- read.csv(http,header = TRUE, stringsAsFactors=FALSE,fileEncoding = "UTF-8");
	if (nrow(reference)==0)	{
		duds <- duds+1;
		} else if (is.null(references))	{
		references <- reference;
		duds <- 0;
		} else	{
		references <- rbind(references,reference);
		duds <- 0;
		}
	reference_no <- reference_no+1;
	}
references <- clear_na_from_matrix(references,"");
references <- put_pbdb_dataframes_into_proper_type(references);
references <- clear_na_from_matrix(references,0);
return(references);
}

#earliest_date <- "2005-01-01"
#latest_date <- "2010-01-01"
#latest_date <- "2010-01-01"
update_pbdb_references <- function(latest_date="",earliest_date="1998-11-01")	{
if (latest_date=="")	latest_date <- strsplit(as.character(Sys.time())," ")[[1]][1];
http <- paste("https://paleobiodb.org/data1.2/taxa/refs.csv?base_name=Life&refs_created_before=",latest_date,"&refs_modified_after=",earliest_date,"&select=auth,class,ops,occs,specs,colls&private&show=both,crmod",sep="");
#http <- paste("https://paleobiodb.org/data1.2/taxa/refs.csv?base_name=Life&refs_modified_after=",latest_date,"&select=auth,class,ops,occs,specs,colls&private&show=both,crmod",sep="");
new_references <- read.csv(http,header=T,stringsAsFactors=F,fileEncoding = "UTF-8");
new_references <- new_references[order(new_references$reference_no),];
new_references <- put_pbdb_dataframes_into_proper_type(new_references);
return(new_references);
}

					##### GENERAL DATA SCRUBBING #######
# collections <- ml_optimized_collections;
# collections$collection_subset <- as.numeric(collections$collection_subset)
# collections$collection_subset[is.na(collections$collection_subset)] <- 0;
lump_subset_sites <- function(collections,no_jack=F)	{
n_colls <- nrow(collections);
subsetted_coll_nos <- sort(as.numeric(unique(collections$collection_subset[collections$collection_subset>0])));
relv_subsetted <- subsetted_coll_nos[!is.na(match(subsetted_coll_nos,collections$collection_no))];
other_subsetted <- subsetted_coll_nos[is.na(match(subsetted_coll_nos,collections$collection_no))];
for (os in 1:length(other_subsetted))	{
	first_example <- collections$collection_no[match(other_subsetted[os],collections$collection_subset)];
	collections$collection_subset[(1:n_colls)[collections$collection_subset==other_subsetted[os]]] <- first_example
	}
collections$collection_subset[collections$collection_subset==collections$collection_no] <- 0;
subsetted_coll_nos <- sort(as.numeric(unique(collections$collection_subset[collections$collection_subset>0])));

if (no_jack)	{
	jack_coll_nos <- subsetted_coll_nos[collections$authorizer[match(subsetted_coll_nos,collections$collection_no)]=="J. Sepkoski"];
	jack_coll_nos <- sort(c(jack_coll_nos,subsetted_coll_nos[collections$authorizer[match(subsetted_coll_nos,collections$collection_no)]=="A. Miller" & collections$enterer[match(subsetted_coll_nos,collections$collection_no)]=="P. Novack-Gottshall"]));

	collections$collection_subset[collections$collection_subset %in% jack_coll_nos] <- 0;
	subsetted_coll_nos <- sort(as.numeric(unique(collections$collection_subset[collections$collection_subset>0])));
	}
rosetta <- collections[,colnames(collections) %in% c("collection_no","collection_subset")];
rosetta <- subset(rosetta,rosetta$collection_subset>0);
return(rosetta);
}

mundify_web_text_boring <- function(web_text)	{
web_text <- gsub("–","-",web_text);
web_text <- gsub("‑","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("‚Äì","–",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","“",web_text);
web_text <- gsub("‚Äù","”",web_text);
web_text <- gsub("‚Äò","‘",web_text);
web_text <- gsub("‚Äô","’",web_text);
web_text <- gsub("â€™","’",web_text);
web_text <- gsub("√°","a",web_text);
web_text <- gsub("√©","e",web_text);
#web_text <- gsub("ö","&ouml;",web_text);
#web_text <- gsub("ü","&uuml;",web_text);
web_text <- gsub("&#321;","L",web_text);
web_text <- gsub("&#305;","i",web_text);
web_text <- gsub("≈ç","o",web_text);
web_text <- transmogrify_diacritics(funky_text=web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("√∫","u",web_text);
web_text <- gsub("√Ω","y",web_text);
web_text <- gsub("\t","     ",web_text);
return(web_text);
}

mundify_web_text_pizzaz <- function(web_text)	{
web_text <- gsub("–","-",web_text);
web_text <- gsub("‑","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("‚Äì","–",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("&#8722;","-",web_text);
web_text <- gsub("‚Äú","“",web_text);
web_text <- gsub("â€œ","“",web_text);
web_text <- gsub("‚Äù","”",web_text);
web_text <- gsub("â€<9d>","”",web_text);
web_text <- gsub("‚Äò","‘",web_text);
web_text <- gsub("‚Äô","’",web_text);
web_text <- gsub("‚Äô","'",web_text);
web_text <- gsub("&#8209;","-",web_text);
web_text <- gsub(" &#730;","˚",web_text);
web_text <- gsub("&#730;","˚",web_text);
web_text <- gsub("Ã","à",web_text);
web_text <- gsub("√Ñ","Ä",web_text);
web_text <- gsub("√Ö","Å",web_text);
web_text <- gsub("√°","á",web_text);
web_text <- gsub("Ã¡","á",web_text);
web_text <- gsub("·","á",web_text);
web_text <- gsub("√•","å",web_text);
web_text <- gsub("√£","ã",web_text);
web_text <- gsub("√§","ä",web_text);
web_text <- gsub("√†","à",web_text);
web_text <- gsub("&#259;","ă",web_text);
web_text <- gsub("Ã¦","æ",web_text);
web_text <- gsub("Œ±","α",web_text);
web_text <- gsub("√ß","ç",web_text);
web_text <- gsub("ƒç","ç",web_text);
web_text <- gsub("&#268;","Č",web_text);
web_text <- gsub("Ã‡","Ç",web_text);
web_text <- gsub("&#269;","č",web_text);
web_text <- gsub("&#263;","ć",web_text);
web_text <- gsub("eÃà","ë",web_text);
web_text <- gsub("eÃÄ","è",web_text);
web_text <- gsub("√©","é",web_text);
web_text <- gsub("eÌ<81>","é",web_text);
web_text <- gsub("e&#769;","é",web_text);
web_text <- gsub("√™","ê",web_text);
web_text <- gsub("eÌ€","è",web_text);
web_text <- gsub("Ã¨","è",web_text);
web_text <- gsub("√®","è",web_text);
web_text <- gsub("√©","é",web_text);
web_text <- gsub("Ã©","é",web_text);
web_text <- gsub("Ãª","ê",web_text);
web_text <- gsub("√©","é",web_text);
web_text <- gsub("ƒó","ė",web_text);
web_text <- gsub("&#283;","ě",web_text);
web_text <- gsub("ƒ±","ı",web_text);
web_text <- gsub("√≠","í",web_text);
web_text <- gsub("√Æ","î",web_text);
web_text <- gsub("&#321;","Ł",web_text);
web_text <- gsub("&#322;","ł",web_text);
web_text <- gsub("√±","ñ",web_text);
web_text <- gsub("&#324;","ń",web_text);
web_text <- gsub("√ñ","Ö",web_text);
web_text <- gsub("√∏","o",web_text);
web_text <- gsub("Ã³","ó",web_text);
web_text <- gsub("√∂","ö",web_text);
web_text <- gsub("Ã¶","ö",web_text);
web_text <- gsub("√µ","õ",web_text);
web_text <- gsub("ˆ","ö",web_text);
web_text <- gsub("&#337;","ő",web_text);
web_text <- gsub("√ï","Õ",web_text);
web_text <- gsub("√¥","ô",web_text);
web_text <- gsub("Ã¶","ö",web_text);
web_text <- gsub("Ã¸","ø",web_text);
web_text <- gsub("&#345;","ř",web_text);
web_text <- gsub("≈†","Š",web_text);
web_text <- gsub("≈°","š",web_text);
web_text <- gsub("Å¡","š",web_text);
web_text <- gsub("&#351;","ş",web_text);
web_text <- gsub("√ú","Ü",web_text);
web_text <- gsub("√º","ü",web_text);
web_text <- gsub("Ã¼","ü",web_text);
web_text <- gsub("≈´","ū",web_text);
web_text <- gsub("√Ω","ý",web_text);
web_text <- gsub("Ã½","ý",web_text);
#web_text <- gsub("Ã�Å¾","íž",web_text);
web_text <- gsub("Å¾","ž",web_text);
web_text <- gsub("\t","     ",web_text);
return(web_text);
}

# Change ü to ue, etc.
#funky_text <- "Gromczakiewicz-Łomnicka"; transmogrify_diacritics(funky_text);
transmogrify_diacritics <- function(funky_text)	{
funky_text <- gsub("̂","",funky_text);
funky_text <- gsub("̊","",funky_text);
j <- strsplit(as.character(funky_text),split="",fixed=TRUE)[[1]];
eek <- c("à","á","â","ã","ā","ă","ȧ","ä","ả","å","ǎ","ȁ","ȃ","ą","ạ","ḁ","ẚ","ầ","ấ","ẫ","ẩ","ằ","ắ","ẵ","ẳ","ǡ","ǟ","ǻ","ậ","ặ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "a";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "A";
eek <- c("æ","ǽ","ǣ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ae";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Ae";
eek <- c("ć","ĉ","ċ","č","ƈ","ç","ḉ","ȼ","ć");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "c";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "C";
eek <- c("ḋ","ɗ","ḍ","ḏ","ḑ","ḓ","ď","đ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "d";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "D";
eek <- c("è","é","ê","ẽ","ē","ĕ","ė","ë","ẻ","ě","ȅ","ȇ","ẹ","ȩ","ę","ḙ","ḛ","ề","ế","ễ","ể","ḕ","ḗ","ệ","ḝ","é");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "e";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "E";
eek <- c("ǵ","ĝ","ḡ","ğ","ġ","ǧ","ɠ","ģ","ǥ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "g";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "G";
eek <- c("ĥ","ḣ","ȟ","ḥ","ḩ","ḫ","ẖ","ħ","ⱨ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "h";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "H";
eek <- c("ì","í","î","ĩ","ī","ĭ","ı","ï","ỉ","ǐ","ị","į","ȉ","ȋ","ḭ","ɨ","ḯ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "i";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "I";
eek <- c("ĵ","ǰ","ȷ","ɉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "j";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "J";
eek <- c("ḱ","ǩ","ḵ","ƙ","ḳ","ĸ","ⱪ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "k";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "K";
eek <- c("ĺ","ľ","ŀ","ł","ƚ","ł");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "l";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "L";
eek <- c("ḿ","ṁ","ṃ","ɱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "m";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "M";
eek <- c("ǹ","ń","ñ","ṅ","ň");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "n";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "N";
eek <- c("ò","ó","ô","õ","ō","ŏ","ȯ","ö","ỏ","ő","ǒ","ô","ȍ","ȏ","ơ","ǫ","ọ","ɵ","ø");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "o";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "O";
eek <- c("ṕ","ṗ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "p";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "P";
eek <- c("ŕ","ṙ","ř","ȑ","ȓ","ṛ","ŗ","ṟ","ṝ","ɍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "r";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "R";
eek <- c("ś","ŝ","ṡ","š","ṣ","ș","ş","ȿ","ṥ","ṧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "s";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "S";
eek <- c("ß");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ss";
eek <- c("ṫ","ẗ","ť","ƫ","ṭ","ț","ţ","ṱ","ṯ","ŧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "t";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "T";
eek <- c("ù","ú","û","ũ","ū","ŭ","ü","ủ","ű","ǔ","ȕ","ȗ","ụ","ṳ","ų","ṷ","ṵ","ů");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "u";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "U";
eek <- c("ṽ","ṿ","ⱱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "v";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "V";
eek <- c("ẁ","ẃ","ŵ","ẇ","ẅ","ẘ","ẉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "w";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "W";
eek <- c("ẋ","ẍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "x";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "X";
eek <- c("ỳ","ý","ŷ","ȳ","ẏ","ÿ","ỷ","ẙ","ɏ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "y";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Y";
eek <- c("ź","ẑ","ż","ž","ȥ","ẓ","ẕ","ƶ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "z";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Z";
j <- j[j!="̈"];
j <- j[j!="̀"];
j <- j[j!="̌"];
return(paste(j,collapse=""));
}

mundify_web_text <- function(web_text)	{
web_text <- gsub("̂","",web_text);
web_text <- gsub("­","",web_text);
web_text <- gsub("–","-",web_text);
web_text <- gsub("‐","-",web_text);
web_text <- gsub("‑","-",web_text);
web_text <- gsub("−","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("ʹ","\'",web_text);
web_text <- gsub("\u009d","",web_text);
web_text <- gsub("\035","",web_text);
web_text <- gsub("\x8","",web_text);
web_text <- gsub("‚Äì","-",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","\"",web_text);
web_text <- gsub("‚Äù","\"",web_text);
web_text <- gsub("‚Äò","\'",web_text);
web_text <- gsub("‚Äô","\'",web_text);
web_text <- gsub("â€™","\'",web_text);
web_text <- gsub("é","e",web_text);
web_text <- gsub("√©","e",web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("ý","y",web_text);
web_text <- gsub("√Ω","y",web_text);
return(web_text);
}

mundify_web_text_dull <- function(web_text)	{
web_text <- gsub("­","",web_text);
web_text <- gsub("–","-",web_text);
web_text <- gsub("‐","-",web_text);
web_text <- gsub("‑","-",web_text);
web_text <- gsub("−","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("ʹ","\'",web_text);
web_text <- gsub("‚Äì","-",web_text); # stri_trans_general("‚Äì","latin-ascii")
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","\"",web_text);
web_text <- gsub("‚Äù","\"",web_text);
web_text <- gsub("‚Äò","\'",web_text);
web_text <- gsub("‚Äô","\'",web_text);
web_text <- gsub("â€™","\'",web_text);
web_text <- gsub("√°","a",web_text);
web_text <- gsub("√©","e",web_text);
web_text <- gsub("&#321;","L",web_text);
web_text <- gsub("≈ç","o",web_text);
web_text <- transmogrify_diacritics(funky_text=web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("√∫","u",web_text);
web_text <- gsub("√Ω","y",web_text);
web_text <- gsub("\t","     ",web_text);
return(web_text);
}

# routine to cleanup rock unit names; named_rock_unit <- "Wynniatt (2 And 3)"; dehyphenate=T;delete_rock_type=T;delete_informal=T
#named_rock_unit <- "'Hailey Shales'"
#mundify_rock_unit_names(named_rock_unit,T,T,T);
# dehyphenate=T; delete_rock_type=T; delete_informal=T;
mundify_rock_unit_names <- function(named_rock_unit,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)	{
# was: clean_rock_unit_names
# named_rock_unit: string giving the name of a formation, member or group
# delete_rock_type: if true, the "Burgess Shale" becomes "Burgess" This is here because workers are
#	inconsistent about including rock-types in formation names
if (is.na(named_rock_unit))	named_rock_unit <- "";

named_rock_unit <- gsub("\u009d","",named_rock_unit);
named_rock_unit <- gsub("\035","",named_rock_unit);
named_rock_unit <- gsub("\x8","",named_rock_unit);
named_rock_unit <- gsub("\x96","n",named_rock_unit);
named_rock_unit <- transmogrify_to_title_case(named_rock_unit);

if (named_rock_unit=="")	{
	return(named_rock_unit);
	} else	{
	nru <- named_rock_unit;		# for debugging
#	named_rock_unit <- str_lowercase(named_rock_unit);		# function no longer works
#	named_rock_unit <- str_ucfirst(named_rock_unit);		# function no longer works
	named_rock_unit <- gsub(" \\(\\?\\)","",named_rock_unit);
	named_rock_unit <- gsub("\\(\\?\\) ","",named_rock_unit);
	named_rock_unit <- gsub("\\(\\?\\)","",named_rock_unit);
	named_rock_unit <- gsub("\\?","",named_rock_unit);
	named_rock_unit <- gsub("\"", "",named_rock_unit);
	named_rock_unit <- gsub("'s","s",named_rock_unit)
	named_rock_unit <- gsub("'s ", "s ",named_rock_unit);
	named_rock_unit <- gsub("’s ", "s ",named_rock_unit);
	named_rock_unit <- gsub("’", "\\'",named_rock_unit);
	named_rock_unit <- gsub("“", "",named_rock_unit);
	named_rock_unit <- gsub("”", "",named_rock_unit);
	named_rock_unit <- gsub("‘", "",named_rock_unit);
	named_rock_unit <- gsub("’", "",named_rock_unit);
	named_rock_unit <- gsub("\'", "",named_rock_unit);
	named_rock_unit <- gsub(" - ","-",named_rock_unit);
	named_rock_unit <- gsub(" — ","-",named_rock_unit);
	named_rock_unit <- gsub("Ste.-","Ste. ",named_rock_unit);
	named_rock_unit <- gsub("Ste-","Ste. ",named_rock_unit);
	named_rock_unit <- gsub("St.-","St. ",named_rock_unit);
	named_rock_unit <- gsub("St-","St. ",named_rock_unit);
	named_rock_unit <- gsub("#","",named_rock_unit);
	named_rock_unit <- gsub("lower part","",named_rock_unit);
	named_rock_unit <- gsub("middle part","",named_rock_unit);
	named_rock_unit <- gsub("upper part","",named_rock_unit);
	named_rock_unit <- gsub("-bearing","",named_rock_unit);
	named_rock_unit <- gsub("-horizon","",named_rock_unit);
	named_rock_unit <- gsub("Ã„","A",named_rock_unit);
	named_rock_unit <- gsub("ã„","a",named_rock_unit);
	named_rock_unit <- gsub("Á","A",named_rock_unit);
	named_rock_unit <- gsub("á","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¡","a",named_rock_unit);
	named_rock_unit <- gsub("ã¡","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¤","a",named_rock_unit);
	named_rock_unit <- gsub("ã¤","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¥","a",named_rock_unit);
	named_rock_unit <- gsub("ã¥","a",named_rock_unit);
	named_rock_unit <- gsub("√†","a",named_rock_unit);
	named_rock_unit <- gsub("√§","a",named_rock_unit);
	named_rock_unit <- gsub("Ã§","c",named_rock_unit);
	named_rock_unit <- gsub("ã§","c",named_rock_unit);
	named_rock_unit <- gsub("Ã©","e",named_rock_unit);
	named_rock_unit <- gsub("ã©","e",named_rock_unit);
	named_rock_unit <- gsub("Ã¨","e",named_rock_unit);
	named_rock_unit <- gsub("Ã±","n",named_rock_unit);
	named_rock_unit <- gsub("Ã–","O",named_rock_unit);
	named_rock_unit <- gsub("Ã¸","o",named_rock_unit);
	named_rock_unit <- gsub("√≥","o",named_rock_unit);
	named_rock_unit <- gsub("√∂","o",named_rock_unit);
	named_rock_unit <- gsub("Ã¶","o",named_rock_unit);
	named_rock_unit <- gsub("Ãµ","o",named_rock_unit);
	named_rock_unit <- gsub("Ã´","o",named_rock_unit);
	named_rock_unit <- gsub("ã´","o",named_rock_unit);
	named_rock_unit <- gsub("Ã¼","u",named_rock_unit);
	named_rock_unit <- gsub("ã¼","u",named_rock_unit);
	named_rock_unit <- gsub("≈´","u",named_rock_unit);
	named_rock_unit <- gsub("√Æ","i",named_rock_unit);
	named_rock_unit <- gsub("√æ","i",named_rock_unit);
	named_rock_unit <- gsub("≈†","S",named_rock_unit);
	named_rock_unit <- gsub("≈°","s",named_rock_unit);
	named_rock_unit <- gsub("‚Äì","-",named_rock_unit);
	named_rock_unit <- gsub("√•","a",named_rock_unit);
	named_rock_unit <- gsub("&#367;","u",named_rock_unit);
	named_rock_unit <- gsub("&#945;","α",named_rock_unit);

	named_rock_unit <- gsub(" Part Of "," ",named_rock_unit);
	named_rock_unit <- transmogrify_diacritics(named_rock_unit);

#	named_rock_unit_2 <- gsub("-"," ",named_rock_unit);
	n_r_u <- strsplit(named_rock_unit," ")[[1]];
	rock_names <- length(n_r_u);
	if (n_r_u[rock_names]=="Part" && rock_names>1)	{
		rock_names <- rock_names-1;
		n_r_u <- n_r_u[1:rock_names]
		}

	thesaurus <- cbind(c("ft","ft.","mt","mt.","ste.","ste","st","st.","ls","lst","limeston","limstone","limestonee","qzt.","sh","claystones","limestones","slates","shales","siltstones","sandstones"),
					   c("Fort","Fort","Mountain","Mountain","Sainte","Sainte","Saint","Saint","Limestone","Limestone","Limestone","Limestone","Limestone","Quartzite","Shale","Claystone","Limestone","Slate","Shale","Siltstone","Sandstone"));
	th_length <- nrow(thesaurus);

	edit <- (1:rock_names)[tolower(n_r_u) %in% thesaurus[,1]];
	if (length(edit)>0)	{
		n_r_u[edit] <- thesaurus[match(tolower(n_r_u[edit]), thesaurus[,1]),2]
		if (n_r_u[1]=="Mountain")	n_r_u[1] <- "Mount";
		}

	bad_words <- c("basal","bed","beds","between","biofacies","biozone","contact","couches","cyclothem","cycle","facies","facies","fazies","fm.","fm","formacion","formation","horizons","horizon","layer","level","lagoonal","member","mbr","mb","mb.","miembro","niveau","portion","section","series","shelly","stage","standard","suite","subst.","subunit","subsuite","subzone","unit","units","unknown","unnamed","zone","bottom","top","(lower)","(middle)","(upper)","(bottom)","(top)","regional","seam","thin-bedded","undifferentiated");
	uncensored <- (1:rock_names)[!tolower(n_r_u) %in% bad_words];

	named_rock_unit <- paste(n_r_u[uncensored],collapse = " ");

	if (dehyphenate)	{
		named_rock_unit <- gsub("-"," ",named_rock_unit);
		named_rock_unit <- gsub("—"," ",named_rock_unit);
		}

	if (delete_rock_type)	{
		named_rock_unit_2 <- gsub("-"," ",named_rock_unit);
		named_rock_unit_2 <- gsub("—"," ",named_rock_unit_2);
		named_rock_unit_2 <- gsub("/"," ",named_rock_unit_2);
		named_rock_unit_2 <- gsub("&"," ",named_rock_unit_2);
		n_r_u <- strsplit(named_rock_unit_2," ")[[1]];
		rock_names <- length(n_r_u);
#		bad_words_2 <- c("argillaceous","ashes","ash","calcaerous","calcaire","carbonate","chalk","cherts","chert","clay","claystone","claystones","conglomerates","conglomerate","coquina","coquinas","dolomites","dolomite","dolostones","dolostone","flags","glauconites","glauconite","glauconitics","glauconitic","gres","grauwacke","greywacke","greywackes","grits","grit","kalk","limestone","limestones","limeston","limstone","ls.","ls","lst","lst.","marlstones","marlstone","marl","marls","marly","micrites","micrite","mudstones","mudstone","ooid","ooids","phosphatics","phosphatic","phosphorite","phosphorites","qzt.","quartzite","quartzites","sandstone","sandstones","shales","schichten","schistes","shale","shaly","siltstones","siltstone","tillite","tillites","tuff","tuffs","volcanic","volcanics");
#		sedimentary_rocks <- sedimentary_rocks;
		uncensored <- (1:rock_names)[!tolower(n_r_u) %in% sedimentary_rocks];
#		if (length(uncensored) < (rock_names-1))	{
		if (length(uncensored) < rock_names)	{
			censored <- (1:rock_names)[tolower(n_r_u) %in% sedimentary_rocks];
			accessories <- unique(sort(c("and","et","y","de","di","del","du","de","a","della","del","y","di")));
			# get the "and" & "of the" words out if they link rock types to proper names
			poss_censored <- (1:rock_names)[tolower(n_r_u) %in% accessories];
			if (length(poss_censored)>0 && length(censored)>0)
				censored <- sort(c(censored,poss_censored[poss_censored>=min(censored-1) & poss_censored<=max(censored+1)]));
			uncensored <- uncensored[!uncensored %in% censored];
#			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "and"])>0)	{
#				et <- uncensored[tolower(n_r_u[uncensored]) %in% "and"];
#				if (et > min(censored) && et < max(censored))
#					uncensored <- (1:length(uncensored))[tolower(n_r_u[uncensored])!="and"];
#				}
#			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "et"])>0)	{
#				et <- uncensored[tolower(n_r_u[uncensored]) %in% "et"];
#				if (et > min(censored) && et < max(censored))
#					uncensored <- (1:length(uncensored))[tolower(n_r_u[uncensored])!="et"];
#				}
			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "or"])>0)	{
				et <- uncensored[tolower(n_r_u[uncensored]) %in% "or"];
				if (et > min(censored) && et < max(censored))
					uncensored <- (1:length(uncensored))[n_r_u[uncensored]!="or"];
				}
#			if (length(censored)==1 && censored<rock_names)	{
#				french <- c("du","de","a","della","del","y","di");
#				if (!is.na(match(tolower(n_r_u[censored+1]),french)))	{
#					uncensored <- uncensored[uncensored!=censored+1];
#					}
#				}
			}
		named_rock_unit <- paste(n_r_u[uncensored],collapse = " ");
#		print(named_rock_unit);
		}

	if (delete_informal)	{
		informals <- c("basal","base","inferieur","lower","lowermost","lowest","middle","upper","uppermost","superieur","informal","tongue","equivalent");
		n_r_u <- strsplit(named_rock_unit," ")[[1]];
		rock_names <- length(n_r_u);
		# roman numeral madness
		n_r_u <- n_r_u[!tolower(n_r_u) %in% tolower(as.character(as.roman(1:50)))];
		n_r_u <- n_r_u[!tolower(n_r_u) %in% tolower(paste(as.roman(1:50),"a",sep=""))];
		n_r_u <- n_r_u[!tolower(n_r_u) %in% tolower(paste(as.roman(1:50),"b",sep=""))];
		n_r_u <- n_r_u[!tolower(n_r_u) %in% tolower(paste(as.roman(1:50),"c",sep=""))];
		n_r_u <- n_r_u[!tolower(n_r_u) %in% tolower(paste(as.roman(1:50),"d",sep=""))];
		rock_names <- length(n_r_u);
		formals <- (1:rock_names)[!tolower(n_r_u) %in% informals];
		named_rock_unit <- paste(n_r_u[formals],collapse = " ");
		named_rock_unit <- gsub("1", "",named_rock_unit);
		named_rock_unit <- gsub("2", "",named_rock_unit);
		named_rock_unit <- gsub("3", "",named_rock_unit);
		named_rock_unit <- gsub("4", "",named_rock_unit);
		named_rock_unit <- gsub("5", "",named_rock_unit);
		named_rock_unit <- gsub("6", "",named_rock_unit);
		named_rock_unit <- gsub("7", "",named_rock_unit);
		named_rock_unit <- gsub("8", "",named_rock_unit);
		named_rock_unit <- gsub("9", "",named_rock_unit);
		named_rock_unit <- gsub("0", "",named_rock_unit);

		if (sum(letters %in% tolower(named_rock_unit))==1)
			named_rock_unit <- "";
		}
	named_rock_unit <- gsub(" \\( And \\)","",named_rock_unit);
	named_rock_unit <- gsub(" \\( and \\)","",named_rock_unit);
#	print(named_rock_unit);
	if (named_rock_unit=="")	{
		return(named_rock_unit);
		} else	{
		named_rock_unit <- gsub(" )",")",named_rock_unit);
		j <- strsplit(as.character(named_rock_unit),split="",fixed=TRUE)[[1]];
		if (j[1]=="\\'" || j[1]=="'")
			j <- j[2:length(j)];
		if (j[length(j)]=="\\'" || j[length(j)]=="'")
			j <- j[1:(length(j)-1)];
		if (length(j)>0 && j[length(j)]==" ")	{
			if (length(j)>1)	{
				while (j[length(j)]==" " && length(j)>1)	j <- j[1:(length(j)-1)];
#				rnr <- c();
#				for (i in 1:length(j))	{
#					rnr<- paste(rnr,j[i],sep="");
#					}
#				named_rock_unit <- rnr;
				named_rock_unit <- paste(j,collapse="");
				} else {
				named_rock_unit <- "";
				}
			}
		if (length(j)>0 && j[1]==" ")	{
			while (j[1]==" " && length(j)>1)	j <- j[2:length(j)];
#			rnr <- c();
#			for (i in 1:length(j))	{
#				rnr<- paste(rnr,j[i],sep="");
#				}
			named_rock_unit <- paste(j,collapse="");
#			named_rock_unit <- rnr;
			}
		named_rock_unit <- gsub("  "," ",named_rock_unit);
		detritus <- c(",","/","."," and ","and"," or ","or"," ");
		if (sum(detritus==tolower(named_rock_unit))==length((named_rock_unit)))
			named_rock_unit <- "";
#		print(named_rock_unit);

#		print(c(nru,named_rock_unit))		# for debugging
		return(named_rock_unit);
		}
	}
#return(named_rock_unit);
}

expello_formation_name_from_member_name <- function(rock_unit_name)	{
rock_unit_name <- gsub("\\)","",rock_unit_name);
formation <- simplify2vector(strsplit(rock_unit_name,split = " \\("))[1];
member <- simplify2vector(strsplit(rock_unit_name,split = " \\("))[2];
member <- gsub(formation,"",member);
member <- simplify2array(strsplit(member,split = " ")[1]);
return(paste(paste(formation,collapse=" ")," (",paste(member,collapse=" "),")",sep=""));
}

#named_rock_unit <- "Dunedin Pine Point"
# sometimes 2+ rock units are entered; this will split them
#split_lumped_pbdb_rock_units(named_rock_unit)
split_lumped_pbdb_rock_units <- function(named_rock_unit)	{
orig_name <- named_rock_unit;
named_rock_unit <- gsub(" And ","&",named_rock_unit);
named_rock_unit <- gsub(" and ","&",named_rock_unit);
named_rock_unit <- gsub(" Or ","&",named_rock_unit);
named_rock_unit <- gsub(" or ","&",named_rock_unit);
named_rock_unit <- gsub(" / ","&",named_rock_unit);
named_rock_unit <- gsub("/","&",named_rock_unit);
named_rock_unit <- gsub(",","&",named_rock_unit);
named_rock_unit <- gsub(" &","&",named_rock_unit);
named_rock_unit <- gsub("& ","&",named_rock_unit);
#splitters <- c("-","—","/","and","or",",");
#splitters <- c("/","and","or",",","&");
n_r_u <- as.vector(unlist(strsplit(named_rock_unit,"&")));
#n_r_u <- strsplit(named_rock_unit," ")[[1]];
#n_r_u <- strsplit(named_rock_unit,splitters)[[1]];
rock_names <- length(n_r_u);
if (rock_names==1)	{
	return(c(orig_name,""));
	} else	{
	return(n_r_u[1:2]);
	}
}

divido_lumped_rock_units <- function(named_rock_unit)	{
named_rock_unit <- gsub("-"," - ",named_rock_unit);
named_rock_unit <- gsub("—"," — ",named_rock_unit);
named_rock_unit <- gsub("/"," / ",named_rock_unit);
named_rock_unit <- gsub("&"," and ",named_rock_unit);
named_rock_unit <- gsub(","," , ",named_rock_unit);
named_rock_unit <- gsub("  "," ",named_rock_unit);
#splitters <- c("-","—","/","and","or",",");
splitters <- c("/","and","or",",");
n_r_u <- strsplit(named_rock_unit," ")[[1]];
#n_r_u <- strsplit(named_rock_unit,splitters)[[1]];
rock_names <- length(n_r_u);
splits <- (1:rock_names)[tolower(n_r_u) %in% splitters]
if (length(splits)>0)	{
	splits <- c(0,splits,rock_names+1);
	multi_rock_units <- c()
	for (i in 1:(length(splits)-1))	{
		multi_rock_units <- c(multi_rock_units,paste(n_r_u[(splits[i]+1):(splits[i+1]-1)],collapse=" "));
		}
	return(multi_rock_units);
	}	else	{
	return(c(named_rock_unit,""));
	}
}

deformalize_rock_unit_names <- function(rock_unit_name,informals=c("lower","middle","upper","uppermost","Lower","Middle","Upper","Uppermost"))	{
split_formation <- stringr::str_split_fixed(rock_unit_name,pattern=" ",n=2);
if (!is.na(match(tolower(split_formation[[1]]),informals)))	{
	if (length(split_formation[1,])==2)	{
		rock_unit_name <- split_formation[[2]];
		} else	{
		rock_unit_name <- ""
		}
	}
return(rock_unit_name)
}

transmogrify_informal_formation_addendum_to_member <- function(rock_unit_name,informals=c("lower","middle","upper","Lower","Middle","Upper"))	{
split_formation <- stringr::str_split_fixed(rock_unit_name,pattern=" ",n=2);
if (!is.na(match(tolower(split_formation[[1]]),informals)))	{
	xxx <- strsplit(split_formation[[1]],"")[[1]]
	while (xxx[1]=="")	xxx <- xxx[2:length(xxx)]
	xxx[1] <- toupper(xxx[1]);
	informal <- c();
	for (i in 1:length(xxx))	informal <- paste(informal,xxx[i],sep="");
	rock_unit_name <- paste(split_formation[[2]]," (",informal,")",sep="");
	}
return(rock_unit_name)
}

### routines for zones
# clean taxon entries of ?, aff., etc. zone <- zones[z]; zone <- "Bl3-4"
# editted 2021-02-12: routine to clear out standard incorrect entries fixed
# editted 2021-06-14: fixed cases where things like just NP entered.
mundus_zone <- function(zone,dbug=FALSE)	{
if (dbug)	write.csv(zone,"Test_Zone.csv",row.names = F);
zone <- mundify_web_text_boring(web_text = zone);
zone <- gsub("<sub>","",zone);
zone <- gsub("<-sub>","",zone);

zone <- gsub("- -"," - ",zone);
zone <- gsub(", ","-",zone);
zone <- gsub("—","-",zone);
zone <- gsub("-","-",zone);
zone <- gsub("–","-",zone);
zone <- gsub("-"," - ",zone);
zone <- gsub("  "," ",zone);
zone <- gsub("  "," ",zone);
zone_molecularized <- strsplit(zone," ")[[1]];
zone_molecularized <- zone_molecularized[zone_molecularized!=""];
while (zone_molecularized[length(zone_molecularized)]=="-")	zone_molecularized <- zone_molecularized[1:(length(zone_molecularized)-1)];
if (sum(zone_molecularized %in% c("NP","CC"))>0)	{
	if (length(strsplit(zone,"/")[[1]])>length(zone))	{
		nps <- strsplit(zone,"/")[[1]];
		for (nn in 1:length(nps))	{
			if(gsub(" ","",nps[nn])!=nps[nn])	{
				nps[nn] <- gsub(" ","",nps[nn]);
				} else	{
				if (sum(zone_molecularized %in% "NP")>0)	{
					nps[nn] <- paste("NP",nps[nn],sep="");
					} else if (sum(zone_molecularized %in% "CC")>0)	{
					nps[nn] <- paste("CC",nps[nn],sep="");
					}
				}
			}
		zone_molecularized <- paste(nps,collapse=" + ");
		} else	{
		nps <- (1:length(zone_molecularized))[zone_molecularized %in% c("NP","CC")];
		for (nn in length(nps):1)	{
			npn <- nps[nn]+1;
			if (npn < length(nps) && as.numeric(zone_molecularized[npn])>0 && as.numeric(zone_molecularized[npn])<30)	{
				zone_molecularized[nps[nn]] <- paste(zone_molecularized[nps[nn]:npn],collapse = "");
				zone_molecularized <- zone_molecularized[!(1:length(zone_molecularized)) %in% npn]
				}
			}
		}
	}

zone <- paste(zone_molecularized, collapse = " ");
#Cf5 + Cf6
zone <- gsub("  "," ",zone);
zone <- gsub("  "," ",zone);
#zone <- gsub("\xd3","",zone) √
#zone <- gsub("\xd2","",zone)
zone <- gsub("≈í¬±","α",zone);
zone <- gsub("√ü","ß",zone);
zone <- gsub("‚àö√º","ß",zone);
zone <- gsub("≈í‚â•","γ",zone);
zone <- gsub("≈í¬•","δ",zone);
zone <- gsub(" (\\?)","",zone);
zone <- gsub("(\\?)","",zone);
zone <- gsub("\\? ","",zone);
zone <- gsub(" \\?"," ",zone);
zone <- gsub("\\?","",zone);
zone <- gsub("—","-",zone);
zone <- gsub("-","-",zone);
zone <- gsub("–","-",zone);
zone <- gsub("/" ,"-",zone);
zone <- gsub("‚Äò,Äò","\"",zone);
zone <- gsub("“","\"",zone);
zone <- gsub("”","\"",zone);
zone <- gsub("‘","\"",zone);
zone <- gsub("’","\"",zone);
zone <- gsub("zone of ","",zone);
zone <- gsub("Zone of ","",zone);
zone <- gsub("lower part "," ",zone);
zone <- gsub("middle part "," ",zone);
zone <- gsub("upper part "," ",zone);
zone <- gsub(" lower part"," ",zone);
zone <- gsub(" middle part"," ",zone);
zone <- gsub(" upper part"," ",zone);
zone <- gsub("prob. ","",zone);
zone <- gsub("probably ","",zone);
zone <- gsub("MP","MP ",zone);
zone <- gsub("MP  ","MP ",zone);
zone <- gsub("MN","MN ",zone);
zone <- gsub("MN  ","MN ",zone);
zone <- gsub(" \\(Graptolite\\)" ,"",zone);
zone <- gsub(" \\(Conodont\\)" ,"",zone);
zone <- gsub(" \\(Trilobite\\)" ,"",zone);
zone <- gsub("graptolites" ,"",zone);
zone <- gsub("topmost" ,"",zone);
zone <- gsub("Topmost" ,"",zone);
zone <- gsub("lower - upper ","",zone);
zone <- gsub(" of the ","",zone);
for (py in 1:length(publication_years))
	zone <- gsub(publication_years[py],"",zone);
zone <- gsub("\\(\\)","",zone);
zone <- gsub("\\( \\)","",zone);

if (zone != "")	{
	zone_detritus <- c("assemblage","asssemblage","biozone","zone","zones","subzone","subzones","level","levels","bed","beds","layer","fauna","interval","local","total","range","ammonite","goniatite","trilobite","conodont","coral","graptolite","reference","base","max:","min:","close","between","of","the","spore","part","s.l.","s.l.\\)","(s.l.\\)","concurrent","LVF","lvf");
	zone_molecularized <- strsplit(zone," ")[[1]];
	zone_molecularized[zone_molecularized %in% c("or","and","to","through")] <- "-";
	molecules <- (1:length(zone_molecularized))[!tolower(zone_molecularized) %in% zone_detritus];
	if (length(molecules)>0)	{
#		zone <- paste(zone_molecularized[molecules],collapse="");
		zone_molecularized <- zone_molecularized[molecules];
		wrong_informals <- c("earliest","early","late","latest");
		informals <- c("lowermost","lower","upper","uppermost","middle");
		if (length(molecules[tolower(zone_molecularized) %in% wrong_informals])>0)	{
			fix_these <- molecules[tolower(zone_molecularized) %in% wrong_informals];
			zone_molecularized[fix_these] <- informals[match(tolower(zone_molecularized[fix_these]),wrong_informals)];
			}
		subsz <- molecules[tolower(zone_molecularized) %in% informals];
		if (length(subsz)>0)	{
			zone_molecularized[subsz] <- tolower(zone_molecularized[subsz]);
		# if zone is entered as Polygnathus costatus upper, rewrite to upper Polygnathus costatus
		# if zone is entered as Polygnathus costatus upper - Polygnathus ensensis lower, rewrite to
		#		upper Polygnathus costatus - lower Polygnathus ensensis
			if (sum(zone_molecularized=="-")==0)	{
				zone_molecularized <- c(zone_molecularized[subsz],zone_molecularized[!molecules %in% subsz]);
				} else	{
				separators <- c(0,molecules[zone_molecularized=="-"],1+max(length(zone_molecularized)));
				for (sz in 1:length(subsz))	{
					sz_sep <- max(separators[separators<subsz[sz]]);
					if (!zone_molecularized[sz_sep+1] %in% informals)	{
						sz_sep_n <- separators[match(sz_sep,separators)+1];
						rearranged <- zone_molecularized[(sz_sep+1):(sz_sep_n-1)]
						informal_term <- rearranged[rearranged %in% informals];
						zone_molecularized[(sz_sep+1):(sz_sep_n-1)] <- rearranged <- c(informal_term,rearranged[!rearranged %in% informals]);
						}
					}
				}
			}
		molecules <- 1:length(zone_molecularized);
		if (zone_molecularized[molecules[length(molecules)]]=="-")
			molecules <- molecules[1:(length(molecules)-1)];
		zone <- paste(zone_molecularized[molecules], collapse = " ");
		if (zone_molecularized[1]=="MN" || zone_molecularized[1]=="MP")	{
			breakpts <- (1:length(zone_molecularized))[zone_molecularized %in% "-"];
			if (length(breakpts)>0)	{
				for (bp in length(breakpts):1)	{
					if (zone_molecularized[breakpts[bp]+1] != zone_molecularized[1])	{
						zone_molecularized <- c(zone_molecularized[1:breakpts[bp]],zone_molecularized[1],zone_molecularized[(breakpts[bp]+1):length(zone_molecularized)]);
						}
					}
				zone <- paste(zone_molecularized,collapse=" ");
				}
			}
		} else	{
		zone <- "";
		}
	}
zone_atomized <- strsplit(zone,"")[[1]];
while (zone_atomized[1]==" " && length(zone_atomized)>1)	zone_atomized <- zone_atomized[2:length(zone_atomized)];
zone <- paste(zone_atomized,collapse="");
zone <- gsub(" - 0" ,"-0",zone);
zone <- gsub(" - 1" ,"-1",zone);
zone <- gsub(" - 2" ,"-2",zone);
zone <- gsub(" - 3" ,"-3",zone);
zone <- gsub(" - 4" ,"-4",zone);
zone <- gsub(" - 5" ,"-5",zone);
zone <- gsub(" - 6" ,"-6",zone);
zone <- gsub(" - 7" ,"-7",zone);
zone <- gsub(" - 8" ,"-8",zone);
zone <- gsub(" - 9" ,"-9",zone);
zone <- zone[!zone %in% 1:9];
if (length(zone)==0)	zone <- "";
return(zone);
}

# Separate "Redlichia chinensis - Kootenia gimmelfarbi" into "Redlichia chinensis" & "Kootenia gimmelfarbi"
# editted 2020-03-04
divido_zone <- function(zone)	{
ddd <- strsplit(zone,split="")[[1]];
ddd[ddd=="+"] <- "&";
zone <- paste(ddd,collapse="")
zone <- gsub("&","-",zone);
zone <- gsub("/" ,"-",zone);
zone <- gsub("\\+","-",zone);			# added 2019-12-10
zone <- gsub(" + ","-",zone);
zone <- gsub("–","-",zone);
zone <- gsub(" -" ,"-",zone);
zone <- gsub("- " ,"-",zone);
multizones <- strsplit(zone,split="-")[[1]];
return(multizones);
}

# standardize "upper / late and "lower / early" -> "upper" and "early"
mundus_stage <- function(stage,dbug=FALSE)	{
if (dbug)	print(stage);
stage <- gsub(" / ","/",stage);
stage <- gsub("Early/Lower","Early",stage);
stage <- gsub("Lower","Early",stage);
stage <- gsub("Late/Upper","Late",stage);
stage <- gsub("Upper","Late",stage);
stage <- transmogrify_diacritics(funky_text=stage);
return(stage);
}

# turn Rossodus manitouensis zone to manitouensis zone
transmogrify_full_zone_names_to_species_names_only <- function(zone)	{
zone <- mundus_zone(zone);
multizones <- divido_zone(zone);
#for (z in 1:length(multizones))	{
z <- 0;
poss_species <- array("",dim=length(multizones));
qualifier <- c("lowermost","lower","middle","upper","uppermost");
while (z < length(multizones))	{
	z <- z+1
	this_zone <- strsplit(multizones[z],split=" ")[[1]];
	zed <- names <- length(this_zone);
	v <- "";
	while (tolower(this_zone[zed])==this_zone[zed] && zed >= 1)	{
		this_zone[zed] <- gsub(")","",this_zone[zed]);
		if (!is.na(match(this_zone[1],qualifier)))	{
			qlf <- match(this_zone[1],qualifier);
			if (v=="") {
				poss_species[z] <- paste(qualifier[qlf]," ",this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
				} else	{
				poss_species[z] <- paste(qualifier[qlf],this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
				}
			} else	{
			poss_species[z] <- paste(this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
			}
		zed <- zed-1;
		v <- " ";
		}
	}
zone_species <- "";
if (sum(poss_species!="") > 0)	{
	to_link <- (1:length(poss_species))[!poss_species %in% ""];	# do this for "Agnostus smithi - Nevadella - Redlichia goofyi" zones
	zz <- 1;
	zone_species <- poss_species[to_link[zz]];					# it will now return "smithi-goofyi"
	while (zz < length(to_link))	{
		zz <- zz+1;
		zone_species <- paste(zone_species,"-",poss_species[to_link[zz]],sep="");
		}
	}
return(zone_species)
}

# turn Rossodus manitouensis zone to Rossodus zone
transmogrify_full_zone_names_to_genus_names_only <- function(zone)	{
zone <- mundus_zone(zone);
multizones <- divido_zone(zone);
if (length(multizones)>0)	{
	taxon_name <- multizones;
	genus_name <- sapply(taxon_name,divido_genus_names_from_species_names);
	return(paste(genus_name,collapse="-"));
	} else	{
	return("");
	}
}

# there are some cases where the PaleoDB does not provide paleogeographic info but Fossil Works does.
reparo_paleodb_paleogeography_with_fossilworks_data <- function(paleodb_collections,fossil_works_geography)	{
ncolls <- nrow(paleodb_collections);
l_p_1 <- (1:ncolls)[is.na(paleodb_collections$paleolat)];
l_p_2 <- (1:ncolls)[paleodb_collections$paleolat==""];
lost_paleogeography <- sort(unique(c(l_p_1,l_p_2)));
retrieved_paleogeography <- match(paleodb_collections$collection_no[lost_paleogeography],fossil_works_geography$collection_no);
lost_paleogeography <- lost_paleogeography[!is.na(retrieved_paleogeography)];
retrieved_paleogeography <- retrieved_paleogeography[!is.na(retrieved_paleogeography)];
still_broke <- sort(unique(c(l_p_1,l_p_2)))[!sort(unique(c(l_p_1,l_p_2))) %in% lost_paleogeography];
paleodb_collections$paleolat[still_broke] <- paleodb_collections$paleolng[still_broke] <- paleodb_collections$geoplate[still_broke] <- 0;
if (length(retrieved_paleogeography)>0)	{
	paleodb_collections$geoplate[lost_paleogeography] <- as.numeric(fossil_works_geography$plate[retrieved_paleogeography]);
	paleodb_collections$paleolat[lost_paleogeography] <- as.numeric(fossil_works_geography$paleolatdec[retrieved_paleogeography]);
	paleodb_collections$paleolng[lost_paleogeography] <- as.numeric(fossil_works_geography$paleolngdec[retrieved_paleogeography]);
	paleodb_collections$geoplate <- expello_na_from_vector(data=paleodb_collections$geoplate,replacement = 0);
	paleodb_collections$paleolat <- expello_na_from_vector(data=paleodb_collections$paleolat,replacement = 0);
	paleodb_collections$paleolng <- expello_na_from_vector(data=paleodb_collections$paleolng,replacement = 0);

	lost_geoplate <- (1:ncolls)[paleodb_collections$geoplate==0];
	}

return(paleodb_collections);
}

# this provides edits for corrections of rock units that cannot currently be edited online
reparo_unedittable_paleodb_rock_identification <- function(paleodb_collections,paleodb_rock_reidentifications)	{
ncolls <- nrow(paleodb_collections);
ttl_re_id <- c();
for (rr in 1:nrow(paleodb_rock_reidentifications))	{
	m1 <- (1:ncolls)[paleodb_collections$formation %in% paleodb_rock_reidentifications$Formation[rr]];
	m2 <- (1:ncolls)[paleodb_collections$zone %in% paleodb_rock_reidentifications$Member[rr]];
	m3 <- (1:ncolls)[paleodb_collections$zone %in% paleodb_rock_reidentifications$Zone[rr]];
	re_id <- m1[m1 %in% m2];
	re_id <- re_id[re_id %in% m3];
	paleodb_collections$formation[re_id] <- as.character(paleodb_rock_reidentifications$New_Formation[rr]);
	paleodb_collections$member[re_id] <- as.character(paleodb_rock_reidentifications$New_Member[rr]);
	paleodb_collections$zone[re_id] <- as.character(paleodb_rock_reidentifications$New_Zone[rr]);
	ttl_re_id <- c(ttl_re_id,re_id);
	}
return(paleodb_collections);
}

# edits for paleodb collections that cannot be editted online
reparo_unedittable_paleodb_collections <- function(paleodb_collections,paleodb_collection_edits)	{
editted_fields <- colnames(paleodb_collection_edits);

editted_colls <- match(paleodb_collections$collection_no,as.numeric(paleodb_collection_edits$collection_no));
editted_colls <- unique(editted_colls[!is.na(editted_colls)]);
colls_to_edit <- match(as.numeric(paleodb_collection_edits$collection_no),paleodb_collections$collection_no);
colls_to_edit <- unique(colls_to_edit[!is.na(colls_to_edit)]);
if (length(editted_colls)!=length(colls_to_edit))	{
	paleodb_collections$collection_no[colls_to_edit] %in% as.numeric(paleodb_collection_edits$collection_no[editted_colls])
	as.numeric(paleodb_collection_edits$collection_no[editted_colls]) %in% paleodb_collections$collection_no[colls_to_edit]
	}

for (ef in 1:length(editted_fields))	{
	if (editted_fields[ef]!="collection_no")	{
		field_to_edit <- match(editted_fields[ef],colnames(paleodb_collections));
		if (is.na(field_to_edit))	{
			relv_col_no <- match(editted_fields[ef],colnames(paleodb_collection_edits));
			dummy <- array("",dim=c(nrow(paleodb_collections),1));
			colnames(dummy) <- editted_fields[ef];
			dummy[colls_to_edit,1] <- paleodb_collection_edits[editted_colls,ef];
			paleodb_collections <- cbind(paleodb_collections,dummy);
			} else	{
#			length(paleodb_collections[colls_to_edit,field_to_edit]) <- length(paleodb_collection_edits[editted_colls,ef]);
			paleodb_collections[colls_to_edit,field_to_edit] <- paleodb_collection_edits[editted_colls,ef];
			}
#		cbind(paleodb_collections[colls_to_edit,field_to_edit],paleodb_collection_edits[editted_colls,2]);
		}
	}
return(paleodb_collections);
}

clean_pbdb_abundances <- function(pbdb_counts,mean_ranges=T)	{
options(warn=-1);
good_entries <- as.numeric(pbdb_counts$abund_value);
nfinds <- nrow(pbdb_counts);
bad_counts <- (1:nfinds)[is.na(good_entries)];
for (bc in 1:length(bad_counts))	{
	n <- bad_counts[bc];
	count_atomized <- strsplit(pbdb_counts$abund_value[n],"")[[1]];
	if (sum(!is.na(as.numeric(count_atomized)))>0)	{
		spcm <- gsub("\\(","",pbdb_counts$abund_value[n]);
		spcm <- gsub("\\)","",spcm);
		spcm <- gsub("\\?","",spcm);
		spcm <- gsub("~","",spcm);
		spcms <- strsplit(spcm,"-")[[1]];
		if (length(spcms)==2 && sum(is.na(as.numeric(spcms)))==0 && mean_ranges)	{
			spcms <- as.numeric(spcms);
			spcm <- round(exp(mean(log(spcms))),0);
			}
		}
	pbdb_counts$abund_value[n] <- spcm;
	}
pbdb_counts$abund_value <- as.numeric(pbdb_counts$abund_value);
options(warn=0);
return(pbdb_counts[!is.na(pbdb_counts$abund_value),]);
}

					##### ROUTINES TO ORGANIZE PALEODB DATA WITH EXTERNAL DATABASE ######
# add something so that if a range of zones is given and some exceed rocks' zones...
#	then discard those zones.
# routine to convert intervals to some standardized chronostratigraphic scale (e.g., international unts)
reset_paleodb_intervals_to_desired_time_scale <- function(collections,finest_chronostrat,time_scale)	{
ncolls <- nrow(collections);
collections$early_interval <- as.character(collections$early_interval);
collections$late_interval <- as.character(collections$late_interval);
single_bins <- sort(c((1:ncolls)[is.na(collections$late_interval)],(1:ncolls)[collections$late_interval==""]));
collections$late_interval[single_bins] <- collections$early_interval[single_bins]
#for (sb in 1:length(single_bins))
#	collections$late_interval[single_bins[sb]] <- collections$early_interval[single_bins[sb]];
#collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
problem_early_intervals <- (1:ncolls)[is.na(match(collections$early_interval,finest_chronostrat$interval))];
problem_late_intervals <- (1:ncolls)[is.na(match(collections$late_interval,finest_chronostrat$interval))];
pei <- 0;
while (pei < length(problem_early_intervals))	{
	pei <- pei+1;
	coll_no <- problem_early_intervals[pei];
	int_no <- match(collections$early_interval[coll_no],time_scale$interval);
	collections$early_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_lb[int_no]<=finest_chronostrat$ma_lb))];
	if (is.na(match(collections$early_interval[coll_no],finest_chronostrat$interval)))
		collections$early_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$max_ma[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat)
	}

pli <- 0;
while (pli < length(problem_late_intervals))	{
	pli <- pli+1;
	coll_no <- problem_late_intervals[pli];
	int_no <- match(collections$late_interval[coll_no],time_scale$interval);

	collections$late_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_ub[int_no]<finest_chronostrat$ma_lb))];
	if (is.na(collections$late_interval[coll_no]))
		collections$late_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$min_ma[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
	}
return(collections);
}

# refine PaleoDB dates given finer time scale information than the PaleoDB uses.
# redone 2020-05-22 to fix incorrect reversals of early & late intervals by PaleoDB
# redone 2021-02-12 to fix make sure that it recognizes zones used as intervals
redate_paleodb_collections_with_time_scale <- function(paleodb_collections,time_scale,zone_database,stratchron="International")	{
# paleodb_collections: dataframe where:
#	paleodb_collections$early_interval gives oldest possible interval
#	paleodb_collections$late_interval gives youngest possible interval (this is filled if blank)
#	paleodb_collections$max_ma gives early_interal onset (initially from PaleoDB, then from Gradstein)
#	paleodb_collections$min_ma gives late_interal end (initially from PaleoDB, then from Gradstein)
# time_scale: dataframe where:
#	time_scale$interval gives interval name
#	time_scale$ma_lb gives interval onset (given Gradstein et al. 2012)
#	time_scale$ma_ub gives interval end (given Gradstein et al. 2012)
# zone_database: dataframe where:
#	we search for the entered interval names in one of several versions of the zone name
#	zone_database$interval_lb gives onset interval
#	zone_database$interval_ub gives end interval
#	zone_database$ma_lb gives interval onset (given Gradstein et al. 2012)
#	zone_database$ma_ub gives interval end (given Gradstein et al. 2012)

if (length(unique(time_scale$scale))==1)	stratchron <- unique(time_scale$scale)[1];
ncolls <- nrow(paleodb_collections);
#e_a_l <- (1:ncolls)[paleodb_collections$late_interval!=""];
no_late <- (1:ncolls)[paleodb_collections$late_interval==""];
paleodb_collections$late_interval[no_late] <- paleodb_collections$early_interval[no_late];

# there are some zones used as stages: look for these!
paleodb_collections$late_interval[paleodb_collections$late_interval=="unnamed Pridoli stage"] <- "Pridoli";
pbdb_intervals <- unique(c(paleodb_collections$early_interval,paleodb_collections$late_interval[paleodb_collections$late_interval!=""]));
missing_pbdb_intervals <- pbdb_intervals[is.na(match(pbdb_intervals,time_scale$interval))];
standard_time_scale <- subset(time_scale,time_scale$scale==stratchron);
still_missing <- c();
mi <- 0;
while (mi < length(missing_pbdb_intervals))	{
	mi <- mi+1;
	poss_zone <- unique(which(zone_database==tolower(missing_pbdb_intervals[mi]),arr.ind = T)[,1]);
	if (length(poss_zone)>0)	{
		full_zone <- unique(zone_database$zone_sr[poss_zone]);
		if (length(full_zone)==1)	{
			poss_sites <- subset(paleodb_collections,(paleodb_collections$early_interval==missing_pbdb_intervals[mi] | paleodb_collections$late_interval==missing_pbdb_intervals[mi]));
			poss_int <- subset(standard_time_scale,standard_time_scale$ma_ub<zone_database$ma_lb[poss_zone[1]]&standard_time_scale$ma_lb>zone_database$ma_ub[poss_zone[1]]);
			best_int <- poss_int[match(min(poss_int$ma_lb-poss_int$ma_ub),poss_int$ma_lb-poss_int$ma_ub),];
			base_offset <- poss_int$ma_lb-zone_database$ma_lb[poss_zone[1]];
			paleodb_collections$early_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- poss_int$interval[match(min(base_offset[base_offset>=0]),base_offset)];
			top_offset <- poss_int$ma_ub-zone_database$ma_ub[poss_zone[1]];
			paleodb_collections$late_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- poss_int$interval[match(max(top_offset[top_offset<=0]),top_offset)];
#			poss_sites <- subset(paleodb_collections,(paleodb_collections$early_interval==missing_pbdb_intervals[mi] | paleodb_collections$late_interval==missing_pbdb_intervals[mi]));
			paleodb_collections$zone[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- full_zone[1];
#			site_info$early_interval[site_info$early_interval==missing_pbdb_intervals[mi]] <- best_int$interval;
#			site_info$late_interval[site_info$late_interval==missing_pbdb_intervals[mi]] <- best_int$interval;
			} else if (length(full_zone)>1)	{
			poss_sites <- subset(paleodb_collections,(paleodb_collections$early_interval==missing_pbdb_intervals[mi] | paleodb_collections$late_interval==missing_pbdb_intervals[mi]));
			poss_zone <- match(full_zone,zone_database$zone);
			for (ps in 1:length(poss_sites))	{
				lb1 <- poss_sites$max_ma[ps];
				ub1 <- poss_sites$min_ma[ps];
				disjuncts <- c();
				for (pz in 1:length(poss_zone))	{
					lb2 <- zone_database$ma_lb[poss_zone[pz]];
					ub2 <- zone_database$ma_ub[poss_zone[pz]];
					#overlaps <- rbind(overlaps,accersi_temporal_overlap(lb1,ub1,lb2,ub2));
					disjuncts <- c(disjuncts,min(abs(lb2-ub1),abs(lb1-ub2)));
					}
				best_zone <- poss_zone[match(min(disjuncts),disjuncts)];
				poss_int <- subset(standard_time_scale,standard_time_scale$ma_ub<zone_database$ma_lb[best_zone[1]]&standard_time_scale$ma_lb>zone_database$ma_ub[best_zone[1]]);
				base_offset <- poss_int$ma_lb-zone_database$ma_lb[best_zone[1]];
				paleodb_collections$early_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- poss_int$interval[match(min(base_offset[base_offset>=0]),base_offset)];
				top_offset <- poss_int$ma_ub-zone_database$ma_ub[best_zone[1]];
				paleodb_collections$late_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- poss_int$interval[match(max(top_offset[top_offset<=0]),top_offset)];
				paleodb_collections$zone[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- zone_database$zone_sr[best_zone[1]];
				}
			}
		} else if (missing_pbdb_intervals[mi] %in% standard_time_scale$st)	{
		poss_sites <- subset(paleodb_collections,(paleodb_collections$early_interval==missing_pbdb_intervals[mi] | paleodb_collections$late_interval==missing_pbdb_intervals[mi]));
		poss_int <- subset(standard_time_scale,standard_time_scale$st %in% missing_pbdb_intervals[mi]);
		best_int <- poss_int[match(min(poss_int$ma_lb-poss_int$ma_ub),poss_int$ma_lb-poss_int$ma_ub),];
#		base_offset <- poss_int$ma_lb-zone_database$ma_lb[poss_zone[1]];
		paleodb_collections$early_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- best_int$interval;
		paleodb_collections$late_interval[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- best_int$interval;
#			poss_sites <- subset(paleodb_collections,(paleodb_collections$early_interval==missing_pbdb_intervals[mi] | paleodb_collections$late_interval==missing_pbdb_intervals[mi]));
#		paleodb_collections$zone[match(poss_sites$collection_no,paleodb_collections$collection_no)] <- full_zone[1];
#			site_info$early_interval[site_info$early_interval==missing_pbdb_intervals[mi]] <- best_int$interval;
#			site_info$late_interval[site_info$late_interval==missing_pbdb_intervals[mi]] <- best_int$interval;
		} else	{
		still_missing <- c(still_missing,missing_pbdb_intervals[mi]);
		}
#	still_missing <- c(still_missing,missing_pbdb_intervals[mi]);
	}
#match("unnamed Pridoli stage",paleodb_collections$zone)
# 2020-05-22: the PaleoDB sometimes incorrectly flips early & late intervals
diff_late <- (1:ncolls)[paleodb_collections$late_interval!=paleodb_collections$early_interval];
effed <- diff_late[time_scale$ma_ub[match(paleodb_collections$late_interval[diff_late],time_scale$interval)]>=time_scale$ma_lb[match(paleodb_collections$early_interval[diff_late],time_scale$interval)]];
#paleodb_collections$early_interval[diff_late][is.na(match(paleodb_collections$early_interval[diff_late],time_scale$interval))]
if (length(effed)>0)	{
	redone_intervals <- cbind(paleodb_collections$late_interval[effed],paleodb_collections$early_interval[effed]);
	paleodb_collections$early_interval[effed] <- redone_intervals[,1];
	paleodb_collections$late_interval[effed] <- redone_intervals[,2];
	}

early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#paste(paleodb_collections$collection_no[paleodb_collections$early_interval %in% "Freboldi"],collapse=",");
if (sum(is.na(early_intervals))>0 || sum(is.na(late_intervals))>0)	{
	poss_zones <- sort(unique(c(paleodb_collections$early_interval[is.na(early_intervals)],paleodb_collections$late_intervals[is.na(late_intervals)])));
	trouble <- trouble_pz <- c();
	for (pz in 1:length(poss_zones))	{
#		pz <- pz+1;
		xxx <- unique(c(which(zone_database==poss_zones[pz],arr.ind = T)[,1],which(zone_database==tolower(poss_zones[pz]),arr.ind = T)[,1]));
		this_trouble <- c();
		if (length(xxx)==1)	{
			paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx];
			paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx];
			} else if (length(xxx)>1)	{
			if(length(unique(zone_database$interval_lb[xxx]))==1)	{
				temp_zones <- paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]]
				temp_zones[temp_zones==""] <- zone_database$zone_sr[xxx[1]];
				paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]] <- temp_zones;
				paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx[1]];
				} else	{
#				print(paste("eff me down",pz));
				#zone_database$zone[xxx];
				this_trouble <- c(this_trouble,xxx);
				}
			if(length(unique(zone_database$interval_ub[xxx]))==1)	{
				paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx[1]];
#				paleodb_collections$min_ma[paleodb_collections$late_interval %in% poss_zones[pz]];
				} else	{
#				print(paste("eff me up",pz));
				#print(zone_database[xxx,]);
				this_trouble <- c(this_trouble,xxx);
				}
			}
		trouble<- c(trouble,(unique(this_trouble)));
		trouble_pz <- c(trouble_pz,rep(pz,length(unique(this_trouble))));
		}
#	zone_database[trouble,];
	}
early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#old_maxes <- collections$max_ma
#old_mins <- collections$min_ma
# sum(is.na(late_intervals))
#(1:nrow(paleodb_collections))[is.na(late_intervals)]
#paleodb_collections$early_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
#paleodb_collections$late_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
paleodb_collections$max_ma[!is.na(early_intervals)] <- time_scale$ma_lb[early_intervals[!is.na(early_intervals)]];
paleodb_collections$min_ma[!is.na(late_intervals)] <- time_scale$ma_ub[late_intervals[!is.na(late_intervals)]];
misentered <- (1:ncolls)[paleodb_collections$max_ma==paleodb_collections$min_ma];
if (length(misentered)>0)	{
	dummy <- paleodb_collections$early_interval[misentered];
	paleodb_collections$early_interval[misentered] <- paleodb_collections$late_interval[misentered];
	paleodb_collections$late_interval[misentered] <- dummy;
	paleodb_collections$max_ma[misentered] <- time_scale$ma_lb[match(paleodb_collections$early_interval[misentered],time_scale$interval)];
	paleodb_collections$min_ma[misentered] <- time_scale$ma_ub[match(paleodb_collections$late_interval[misentered],time_scale$interval)];
	}
return(paleodb_collections);
}

# refine PaleoDB dates given the zone information
# redone 2021-03-08 to just use zone_database
redate_paleodb_collections_with_zones <- function(paleodb_collections,zone_database,time_scale,emend_paleodb=T)	{
# paleodb_collections: dataframe of collections data downloaded from PaleoDB
# zone_matches: vector giving zones (as numbers) for collections with ';' separating different zones
# zone_database: dataframe with
#	dataframe$ma_lb: onset of zone (with 485.4 = 485.4 million years ago)
#	dataframe$ma_ub: end of zone (with 443.4 = 443.4 million years ago)
# chronostrat: dataframe with chronostratigraphic information.
#	note: it helps to make this include onlyy the finest intervals
# emend_paleodb: if T, then max_ma, min_ma, early_interval & late_interval are edited
#	if F, then ma_lb, ma_ub, interval_lb & interval_ub are appended to paleodb_collections

ncolls <- nrow(paleodb_collections);
coll_w_zones <- (1:ncolls)[paleodb_collections$zone!=""];
zone <- paleodb_collections$zone[coll_w_zones];
paleodb_collections$zone[coll_w_zones] <- pbapply::pbsapply(zone,mundus_zone);
#for (cwz in 1:length(coll_w_zones))	test <- mundus_zone(zone=coll_w_zones[cwz]);
coll_w_zones <- (1:ncolls)[paleodb_collections$zone!=""];
unique_zones <- sort(unique(paleodb_collections$zone[coll_w_zones]));
for (uz in 1:length(unique_zones))	{
	poss_zones <- sort(unique(which(zone_database==unique_zones[uz],arr.ind=T)[,1]));
	if (length(poss_zones)==0)	{
		multizones <- gsub(" -","-",unique_zones[uz]);
		multizones <- gsub("- ","-",multizones);
		multizones <- strsplit(multizones,"-")[[1]];
		if (length(multizones)>1)	{
			poss_zones <- c();
			for (mz in 1:length(multizones))	poss_zones <- c(poss_zones,sort(unique(which(zone_database==multizones[mz],arr.ind=T)[,1])));
			}
		}
	relv_collections <- subset(paleodb_collections,paleodb_collections$zone==unique_zones[uz]);
	rc <- 0;

	while (rc < nrow(relv_collections) && length(poss_zones)>0)	{
		rc <- rc+1;
		lb <- c(relv_collections$max_ma[rc],zone_database$ma_lb[poss_zones]);
		ub <- c(relv_collections$min_ma[rc],zone_database$ma_ub[poss_zones]);
		overlaps <- accersi_temporal_overlap_multiple_cases(lb,ub)[(2:length(lb))-1,];
		if (max(overlaps[,1]-overlaps[,2])>0)	{
			zn <- match(max(overlaps[,1]-overlaps[,2]),overlaps[,1]-overlaps[,2]);
			relv_collections$max_ma[rc] <- min(relv_collections$max_ma[rc],zone_database$ma_lb[poss_zones[zn]]);
			relv_collections$min_ma[rc] <- max(relv_collections$min_ma[rc],zone_database$ma_ub[poss_zones[zn]]);
			if (relv_collections$max_ma[rc]==relv_collections$min_ma[rc])	{
				relv_collections$max_ma[rc] <- zone_database$ma_lb[poss_zones[zn]];
				relv_collections$min_ma[rc] <- zone_database$ma_ub[poss_zones[zn]];
				}
			}
		}
	paleodb_collections$max_ma[match(relv_collections$collection_no,paleodb_collections$collection_no)] <- relv_collections$max_ma;
	paleodb_collections$min_ma[match(relv_collections$collection_no,paleodb_collections$collection_no)] <- relv_collections$min_ma;
	}

return(paleodb_collections);
}

redate_paleodb_collections_with_zone_matches <- function(paleodb_collections,zone_matches,zone_database,time_scale,emend_paleodb=TRUE)	{
# paleodb_collections: dataframe of collections data downloaded from PaleoDB
# zone_matches: vector giving zones (as numbers) for collections with ';' separating different zones
# zone_database: dataframe with
#	dataframe$ma_lb: onset of zone (with 485.4 = 485.4 million years ago)
#	dataframe$ma_ub: end of zone (with 443.4 = 443.4 million years ago)
# chronostrat: dataframe with chronostratigraphic information.
#	note: it helps to make this include onlyy the finest intervals
# emend_paleodb: if T, then max_ma, min_ma, early_interval & late_interval are edited
#	if F, then ma_lb, ma_ub, interval_lb & interval_ub are appended to paleodb_collections
chronostrat_units <- unique(c(unique(zone_database$interval_lb),unique(zone_database$interval_ub),unique(paleodb_collections$early_interval),unique(paleodb_collections$late_interval)));
chronostrat_units <- chronostrat_units[chronostrat_units!=""];
time_scale <- subset(time_scale,time_scale$ma_ub<1.25*max(paleodb_collections$max_ma));
if (length(unique(time_scale$scale[match(chronostrat_units,time_scale$interval)]))==1)	{
	chronostrat <- accersi_hierarchical_timescale(chronostrat_units,time_scale,regional_scale=time_scale$scale[match(chronostrat_units[1],time_scale$interval)]);
	finest_chronostrat <- chronostrat[chronostrat$bin_first==chronostrat$bin_last,];
	} else	{
	chronostrat <- accersi_hierarchical_timescale(chronostrat_units,time_scale);
	chronostrat_b <- accersi_hierarchical_timescale(chronostrat_units,time_scale=subset(time_scale,time_scale$scale=="International"));
	finest_chronostrat <- chronostrat_b[chronostrat_b$bin_first==chronostrat_b$bin_last,];
	}
finest_chronostrat$ma_lb <- 0.001*round(finest_chronostrat$ma_lb/0.001,0);
finest_chronostrat$ma_ub <- 0.001*round(finest_chronostrat$ma_ub/0.001,0);

ncolls <- nrow(paleodb_collections);
colls_w_zones <- (1:ncolls)[zone_matches!=""];
c_w_z <- length(colls_w_zones);
ma_lb <- 0.001*round(paleodb_collections$max_ma/0.001,0);
ma_ub <- 0.001*round(paleodb_collections$min_ma/0.001,0);
interval_lb <- paleodb_collections$early_interval;
redate_these <- (1:ncolls)[is.na(match(interval_lb,chronostrat$interval))];
new_intervals <- chronostrat$interval[match(paleodb_collections$max_ma[redate_these],chronostrat$ma_lb)];
interval_lb[redate_these[!is.na(new_intervals)]] <- new_intervals[!is.na(new_intervals)];
# cases where max_ma is within interval
if (sum(is.na(new_intervals))>0)	{
	re_redate_these <- redate_these[is.na(new_intervals)];
	for (rdt in 1:length(re_redate_these))	{
		bin <- sum(paleodb_collections$max_ma[re_redate_these[rdt]]<=chronostrat$ma_lb)
		interval_lb[redate_these[match(re_redate_these[rdt],redate_these)]] <- chronostrat$interval[bin];
		}
	}
# fill in blank late intervals with early_interval
paleodb_collections$late_interval[paleodb_collections$late_interval==""] <- paleodb_collections$early_interval[paleodb_collections$late_interval==""];
interval_ub <- paleodb_collections$late_interval;
redate_these <- (1:ncolls)[is.na(match(interval_ub,chronostrat$interval))];
new_intervals <- chronostrat$interval[match(paleodb_collections$min_ma[redate_these],chronostrat$ma_ub)];
interval_ub[redate_these[!is.na(new_intervals)]] <- new_intervals[!is.na(new_intervals)];
# cases where min_ma is within interval
if (sum(is.na(new_intervals))>0)	{
	re_redate_these <- redate_these[is.na(new_intervals)];
	for (rdt in 1:length(re_redate_these))	{
		bin <- sum(paleodb_collections$min_ma[re_redate_these[rdt]]<=chronostrat$ma_ub)
		interval_ub[redate_these[match(re_redate_these[rdt],redate_these)]] <- chronostrat$interval[bin];
		}
	}

# make sure to retain original interval_ub if new one is adjacent
#badness <- c();
#ages <- cbind(ma_lb,ma_ub);
for (cz in 1:c_w_z)	{
	coll_no <- colls_w_zones[cz];
	zone_nos <- as.numeric(strsplit(zone_matches[coll_no],";")[[1]]);
	entered_zone <- paleodb_collections$zone[coll_no];
	relv_zones <- zone_database[zone_nos,];
	zone_types <- unique(relv_zones$zone_type);

	group_zones <- relv_zones[1:length(zone_types),];
#	group_zones$zone <- group_zones$zone_sr <- ""
	for (zt in 1:length(zone_types))	{
		group_zones$zone_type[zt] <- zone_types[zt];
		type_zones <- subset(relv_zones,relv_zones$zone_type==zone_types[zt]);
		onset <- match(max(type_zones$ma_lb),type_zones$ma_lb);
		end <- match(min(type_zones$ma_ub),type_zones$ma_ub);
		group_zones$ma_lb[zt] <- type_zones$ma_lb[onset];
		group_zones$interval_lb[zt] <- type_zones$interval_lb[onset];
		group_zones$ma_ub[zt] <- type_zones$ma_ub[end];
		group_zones$interval_ub[zt] <- type_zones$interval_ub[end];
		}
#	zone_overlap <- accersi_minimum_joint_ranges(lb=c(ma_lb[coll_no],group_zones$ma_lb),ub=c(ma_ub[coll_no],group_zones$ma_ub));
	# get overlap between zones and original range
	for (gz in 1:nrow(group_zones))	{
		if (gz==1)	{
			zone_overlap <- accersi_temporal_overlap(lb1=ma_lb[coll_no],ub1=ma_ub[coll_no],lb2=group_zones$ma_lb[gz],ub2=group_zones$ma_ub[gz]);
			} else	{
			zone_overlap <- rbind(zone_overlap,accersi_temporal_overlap(lb1=ma_lb[coll_no],ub1=ma_ub[coll_no],lb2=group_zones$ma_lb[gz],ub2=group_zones$ma_ub[gz]));
			}
		}
	if (sum(zone_overlap$ma_lb>0)>0)	{
		zone_overlap <- subset(zone_overlap,zone_overlap$ma_lb>0);
		# if zones are disjunct within alloted range, then the latest ma_lb will be younger than the oldest ma_ub
		ma_lb[coll_no] <- round(min(zone_overlap$ma_lb[zone_overlap$ma_lb>max(zone_overlap$ma_ub)]),3);
		# if zones are disjunct within alloted range, then the oldest ma_ub will be older than the youngest ma_lb
		ma_ub[coll_no] <- round(max(zone_overlap$ma_ub[zone_overlap$ma_ub<min(zone_overlap$ma_lb)]),3);
		interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
		interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
		} else if (sum(zone_overlap$ma_lb==0)==nrow(zone_overlap) && entered_zone!="")	{
#		print(paste("totally effed",cz,coll_no));
		if (nrow(group_zones)==1)	{
			ma_lb[coll_no] <- round(group_zones$ma_lb[1],3);
			ma_ub[coll_no] <- round(group_zones$ma_ub[1],3);
#			interval_lb[coll_no] <- group_zones$interval_lb[1];
#			interval_ub[coll_no] <- group_zones$interval_ub[1];
			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
			} else	{
#			joint_ranges <- accersi_minimum_joint_ranges(lb=group_zones$ma_lb,ub=group_zones$ma_ub);
			ma_lb[coll_no] <- round(min(group_zones$ma_lb[group_zones$ma_lb>max(group_zones$ma_ub)]),3);
			ma_ub[coll_no] <- round(max(group_zones$ma_ub[group_zones$ma_ub<min(group_zones$ma_lb)]),3);
			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
#			interval_lb[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=ma_lb[coll_no])];
#			interval_ub[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>ma_ub[coll_no])];
			}
#		} else	{
		# in principle, we reach this only if we found 2+ zone taxa in the occcurrences AND if
		#	1+ of those zone taxa has a range within the time entered into the PaleoDB
#		ok_zones <- zone_overlap[zone_overlap$ma_lb!=0,];
#		if (sum(ok_zones$ma_lb>max(ok_zones$ma_ub))>0)	{
#			ma_lb[coll_no] <- round(min(ok_zones$ma_lb[ok_zones$ma_lb>max(ok_zones$ma_ub)]),3);
#			ma_ub[coll_no] <- round(max(ok_zones$ma_ub[ok_zones$ma_ub<min(ok_zones$ma_lb)]),3);
#			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
#			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
#			interval_lb[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=ma_lb[coll_no])];
#			interval_ub[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>ma_ub[coll_no])];
#			} else	{
#			badness <- c(badness,paste("partiallyy effed",cz,coll_no));
#			}
		}
	}

if (emend_paleodb)	{
	paleodb_collections$max_ma <- ma_lb;
	paleodb_collections$min_ma <- ma_ub;
	paleodb_collections$early_interval <- interval_lb;
	paleodb_collections$late_interval <- interval_ub;
	} else	{
	paleodb_collections <- cbind(paleodb_collections,ma_lb=as.numeric(ma_lb),ma_ub=as.numeric(ma_ub),interval_lb=as.character(interval_lb),interval_ub=as.character(interval_ub));
	}
return(paleodb_collections);
}

refine_paleodb_collection_dates_with_zone_data_only <- function(paleodb_collections,paleodb_finds,zone_database,time_scale,hierarchical_chronostrat,finest_chronostrat,examine_finds=T,temporal_precision=0.05)	{
# paleodb_collections: collections downloaded from Paleobiology Database
# paleodb_finds: occurrences downloaded from Paleobiology Database
zone <- unique(paleodb_collections$zone[paleodb_collections$zone!=""]);
zone_ids <- sapply(zone,match_one_collection_zone_to_zone_database,zone_database);
zone_ids <- zone_ids[zone_ids!=""];
zone_ma_lbs <- zone_ma_ubs <- c();
for (zi in 1:length(zone_ids))	{
	zone_nos <- as.numeric(strsplit(zone_ids[zi],";")[[1]]);
	zone_ma_lbs <- c(zone_ma_lbs,zone_database$ma_lb[zone_nos]);
	zone_ma_ubs <- c(zone_ma_ubs,zone_database$ma_ub[zone_nos]);
	}
max_ma <- max(c(paleodb_collections$max_ma),zone_ma_lbs);
min_ma <- min(c(paleodb_collections$min_ma),zone_ma_ubs);

zone_database <- subset(zone_database,zone_database$ma_lb>min_ma);
zone_database <- subset(zone_database,zone_database$ma_ub<max_ma);
zone_database$ma_lb <- temporal_precision*round(zone_database$ma_lb/temporal_precision,0);
zone_database$ma_ub <- temporal_precision*round(zone_database$ma_ub/temporal_precision,0);
ma <- zone_database$ma_lb;
new_interval_lb <- sapply(ma,reassign_intervals_to_uniform_scale,uniform_time_scale=finest_chronostrat,onset=T);
zone_database$interval_lb[new_interval_lb!=""] <- new_interval_lb[new_interval_lb!=""];
ma <- zone_database$ma_ub;
new_interval_ub <- sapply(ma,reassign_intervals_to_uniform_scale,uniform_time_scale=finest_chronostrat,onset=F);
zone_database$interval_ub[new_interval_ub!=""] <- new_interval_lb[new_interval_ub!=""];
if(examine_finds)	{
	zone_info <- match_paleodb_collections_to_possible_zones(paleodb_collections,zone_database,paleodb_finds);
	} else	{
	zone_info <- match_collections_zones_to_zone_database(paleodb_collections,zone_database);
	}

# use zone matches for still more exact dates
zone_matches <- zone_info$zone_matches;
refined_collections <- redate_paleodb_collections_with_zones(paleodb_collections,zone_matches,zone_database,time_scale,emend_paleodb=F);
return(refined_collections);
}

# download occurrences from a collection & look for zone taxa
revelare_zone_taxa_in_one_paleodb_collection <- function(coll_id,zone_database)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
this_coll_finds <- accersi_occurrences_from_one_paleodb_collection(coll_id);
# prepare to separate out just species occurrences
this_coll_finds$identified_rank[this_coll_finds$identified_rank=="subspecies"] <- "species"
this_coll_finds$accepted_rank[this_coll_finds$accepted_rank=="subspecies"] <- "species"
this_coll_finds <- subset(this_coll_finds,this_coll_finds$identified_rank=="species");
if (nrow(this_coll_finds) > 0)	{
	taxon_name <- this_coll_finds$identified_name;
	# replace identified name with species name if we have taxonomic data
	taxon_name[this_coll_finds$accepted_rank=="species"] <- this_coll_finds$accepted_name[this_coll_finds$accepted_rank=="species"];
	taxon_name <- sapply(taxon_name,mundify_taxon_names);
	poss_zones <- c();
	for (tn in 1:length(taxon_name))
		poss_zones <- c(poss_zones,match_one_collection_zone_to_zone_database(zone=taxon_name[tn],zone_database));
	poss_zones <- unique(poss_zones)
	return(paste(poss_zones[poss_zones!=""],collapse=";"));
	} else	{
	return("");
	}
}

# look for zone taxa in already downloaded PaleoDB collections
revelare_zone_taxa_in_paleodb_finds <- function(coll_id,paleodb_finds,zone_database)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
this_coll_finds <- subset(paleodb_finds,paleodb_finds$collection_no==coll_id);
# prepare to separate out just species occurrences
this_coll_finds$identified_rank[this_coll_finds$identified_rank=="subspecies"] <- "species"
this_coll_finds$accepted_rank[this_coll_finds$accepted_rank=="subspecies"] <- "species"
this_coll_finds <- subset(this_coll_finds,this_coll_finds$identified_rank=="species");
this_coll_finds <- this_coll_finds[!this_coll_finds$flags %in% c("uncertain species","uncertain genus, uncertain species"),];
if (nrow(this_coll_finds) > 0)	{
	taxon_name <- this_coll_finds$identified_name;
	# replace identified name with species name if we have taxonomic data
	taxon_name[this_coll_finds$accepted_rank=="species"] <- this_coll_finds$accepted_name[this_coll_finds$accepted_rank=="species"];
	taxon_name <- sapply(taxon_name,mundify_taxon_names);
	poss_zones <- c();
	for (tn in 1:length(taxon_name))
		poss_zones <- c(poss_zones,match_one_collection_zone_to_zone_database(zone=taxon_name[tn],zone_database));
	poss_zones <- unique(poss_zones);
	return(paste(poss_zones[poss_zones!=""],collapse=";"));
	} else	{
	return("");
	}
}

# compare zones from PaleoDB collections to a zone database
match_paleodb_collections_to_possible_zones <- function(paleodb_collections,zone_database,paleodb_finds)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
ncolls <- nrow(paleodb_collections);
paleodb_zone_info <- match_collections_zones_to_zone_database(collections=paleodb_collections,zone_database);
colls_w_zones <- (1:ncolls)[paleodb_collections$zone!=""];
zone_matches <- paleodb_zone_info$zone_matches;
#unmatched_colls <- paleodb_collections$collection_no[colls_w_zones[zone_matches[colls_w_zones]==""]];
unmatched_zone_info <- data.frame(coll_no = as.numeric(paleodb_collections$collection_no[colls_w_zones[zone_matches[colls_w_zones]==""]]),zone=as.character(paleodb_collections$zone[colls_w_zones[zone_matches[colls_w_zones]==""]]),stringsAsFactors = F);

# look for zone taxa in assemblages
colls_wo_zones <- (1:ncolls)[zone_matches==""];
cwoz <- length(colls_wo_zones);
coll_id <- paleodb_collections$collection_no[colls_wo_zones];
if (is.null(paleodb_finds) || paleodb_finds=="")	{
	zone_taxa_present <- sapply(coll_id,revelare_zone_taxa_in_one_paleodb_collection,zone_database);
	zone_matches[colls_wo_zones] <- zone_taxa_present;
	}	else	{
	coll_id <- coll_id[coll_id %in% unique(paleodb_finds$collection_no)];
	zone_taxa_present <- sapply(coll_id,revelare_zone_taxa_in_paleodb_finds,paleodb_finds,zone_database);
	colls_wo_zones <- colls_wo_zones[paleodb_collections$collection_no[colls_wo_zones] %in% coll_id];
	zone_matches[colls_wo_zones] <- zone_taxa_present;
	}

output <- list(zone_matches,unmatched_zone_info);
names(output) <- c("zone_matches","unmatched_zone_info");
return(output);
}

# match a particular zone to an external zone database
match_one_collection_zone_to_zone_database <- function(zone,zone_database)	{
#print(zone);
test <- which(zone_database==zone,arr.ind=T);

# if nothing, the check to see if it is 2+ zones
if (nrow(test)==0)	{
	zones <- divido_zone(zone);
	zones <- zones[zones!=""];
	if (length(zones)>1)	{
		test <- matrix("",nrow=0,ncol=2);
		for (zz in 1:length(zones))
			test <- rbind(test,which(zone_database==zones[zz],arr.ind=T));
		}
	}

# if nothing, then check to see if it is a mismatch because of subgenus assignment
if (nrow(test)==0)	{
	genus <- divido_genus_names_from_species_names(zone);
	species <- divido_species_epithets(taxon_name=zone);
	if(length(strsplit(genus,split=" ")[[1]])==2 && species!="")	{
		gs <- strsplit(genus,split=" ")[[1]];
		jj <- strsplit((gs)[2],split="")[[1]];
		if (jj[1]=="(")	{
			gs[2] <- paste(jj[2:(length(jj)-1)],collapse="");
			zone_a <- paste(gs[1],species);
			zone_b <- paste(gs[2],species);
			test <- matrix("",nrow=0,ncol=2);
			test <- rbind(test,which(zone_database==zone_a,arr.ind=T));
			test <- rbind(test,which(zone_database==zone_b,arr.ind=T));
			}
		}
	}

zone_rows <- unique(test[,1]);
zone_data_reduced <- zone_database[zone_rows,];
zone_data_reduced[,(1:ncol(zone_database))[!colnames(zone_database) %in% c("zone_sr","Regional_Scale","ma_lb","ma_ub")]] <- NULL;
rownames(zone_data_reduced) <- zone_rows;
zone_data_reduced <- unique(zone_data_reduced);
zdb <- paste(unique(rownames(zone_data_reduced)),collapse=";");
return(zdb);
}

# match zones from numerous collections to an external database
match_collections_zones_to_zone_database <- function(collections,zone_database)	{
# reduce zone database to feasible zones
#zone_database <- subset(zone_database,as.numeric(zone_database$ma_lb)<(max(collections$max_ma)+25));
#zone_database <- subset(zone_database,as.numeric(zone_database$ma_ub)>(min(collections$min_ma)-25));

if (is.na(match("zone_epithet",colnames(zone_database))))	{
	taxon_name <- zone <- as.character(zone_database$zone);
	zone_epithet <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
	zone <- as.character(zone_database$zone_sr);
	zone_epithet_sr <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
	zone_database <- cbind(zone_database,zone_epithet,zone_epithet_sr);
	}

ncolls <- nrow(collections);
colls_w_zones <- (1:ncolls)[collections$zone!=""];
zone <- collections$zone[colls_w_zones];
zone <- sapply(zone,mundus_zone);
zone_matches <- rep("",ncolls);
zone_matches[colls_w_zones] <- sapply(zone,match_one_collection_zone_to_zone_database,zone_database);
unmatched_zones <- sort(unique(zone[zone_matches[colls_w_zones]==""]));

output <- list(zone_matches,unmatched_zones);
names(output) <- c("zone_matches","unmatched_zones");
return(output);
}

# construct a thesaurus for zones that provides several ways to find a zone's name.
accersi_zone_thesaurus <- function(zone_data)	{
### construct a zone thesaurus.
zone_thesaurus <- cbind(zone_data$zone,zone_data$zone_sr,zone_data$zone_species,zone_data$zone_species_sr,zone_data$non_taxon_zone,zone_data$non_taxon_zone_label,zone_data$non_taxon_zone_sr,zone_data$non_taxon_zone_label_sr);
zone_thesaurus <- zone_data;
zone_thesaurus$zone_type <- zone_thesaurus$Regional_Scale <- zone_thesaurus$ma_lb <- zone_thesaurus$ma_ub <- zone_thesaurus$interval_lb <- zone_thesaurus$interval_ub <- NULL;
zone_thesaurus <- unique(zone_thesaurus);
return(zone_thesaurus);
}

accersi_rock_to_zone_data <- function(external_rock_to_zone_database)	{
if (external_rock_to_zone_database!="")	{
	file_type <- revelare_file_type(filename = external_rock_to_zone_database);
	if (file_type=="csv")	{
		rock_to_zone_data <- read.csv(file=external_rock_to_zone_database,header=TRUE);
		}	else	{
		rock_to_zone_data <- read.table(file=external_rock_to_zone_database,header=TRUE,sep="\t");
		}
	} else	{
	rock_to_zone_data <- rock_unit_data$rock_to_zone_database;
	}
#rock_to_zone_data <- read.table(file=external_rock_to_zone_database,header=TRUE,stringsAsFactors = TRUE,sep="\t");
rock_to_zone_data$zone <- zone <- as.character(rock_to_zone_data$zone);
rock_to_zone_data$zone <- sapply(zone,mundus_zone);
rock_to_zone_data$zone_sr <- zone <- as.character(rock_to_zone_data$zone_sr);
rock_to_zone_data$zone_sr <- sapply(zone,mundus_zone);
rzd <- nrow(rock_to_zone_data);

zone <- rock_to_zone_data$zone;
zone_species <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
Nontaxon_Zone <- raster::t(sapply(zone,aparecium_nontaxon_zone));
rownames(Nontaxon_Zone) <- NULL;
Nontaxon_Zone <- data.frame(Nontaxon_Zone);
zone <- rock_to_zone_data$zone_sr;
zone_species_sr <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
Nontaxon_zone_sr <- raster::t(sapply(zone,aparecium_nontaxon_zone));
rownames(Nontaxon_zone_sr) <- NULL;
colnames(Nontaxon_zone_sr) <- c("non_taxon_zone_sr","non_taxon_zone_label_sr");
rock_to_zone_data <- cbind(rock_to_zone_data,zone_species,zone_species_sr,Nontaxon_Zone,Nontaxon_zone_sr);
return(rock_to_zone_data);
}

#entered_formations <- collections$formation;
#entered_members <- collections$member;
#entered_onsets <- collections$max_ma;
#entered_rock_and_onsets <-cbind(entered_formations,entered_members,entered_onsets);
paleodb_editing_kluges <- function(entered_rock_and_onsets)	{
entered_formations <- as.character(entered_rock_and_onsets[,1]);
entered_members <- as.character(entered_rock_and_onsets[,2]);
entered_onsets <- as.numeric(as.character(entered_rock_and_onsets[,3]));

ncoll <- length(entered_formations);
nc <- (1:ncoll)[entered_formations=="Ross Point"][(1:ncoll)[entered_formations=="Ross Point"] %in% (1:ncoll)[entered_onsets<541]];
entered_formations[nc] <- "Ross Brook";

nc <- (1:ncoll)[entered_formations=="Synder"];
entered_formations[nc] <- "Snyder";

nc <- (1:ncoll)[entered_formations=="Konprusy"];
entered_formations[nc] <- "Koneprusy";

nc <- (1:ncoll)[entered_members=="Konprusy"];
entered_members[nc] <- "Koneprusy";

nc <- (1:ncoll)[entered_members=="KashongShale"];
entered_members[nc] <- "Kashong Shale";

nc <- (1:ncoll)[entered_formations=="Szydowek"];
entered_formations[nc] <- "Szydlowek";

nc <- (1:ncoll)[entered_members=="D"][(1:ncoll)[entered_members=="D"] %in% (1:ncoll)[entered_formations==""]]
entered_formations[nc] <- "Solvik";
entered_members[nc] <- "Myren";

return(cbind(entered_formations,entered_members));
}

# routine to download PaleoDB data and use external stratigraphic and zone databases to put
#	more exact dates on collections than PaleoDB provides
# study: name of study (anything you like)
# collections: downloaded PaleoDB data.  MUST include:
#	collection_no
#	formation, member and stratgroup
#	zone
# external_strat_database: name of an external database giving detailed information on rocks.
# MUST give:
##	rock_no: number of rock unit (specific formation + member combination)
##	rock_no_sr: number of the senior synonym (which usually is the same as rock_no)
##	formation_no: the number of the formation (to link different members of a formation)
##	formation: name of formation (might be blank if group is the only information)
##	member: name of member (usually blank if only formation or group information given)
##	full_name: formation_Name (member_Name) or formation_Name if no member.  (Blank if only group)
##	group: Stratigraphic group (often blank)
##	interval_lb: name of chronostratigraphc unit in which rock-unit starts (I prefer stage-slices)
##	interval_ub: name of chronostratigraphc unit in which rock-unit ends (I prefer stage-slices)
##  ma_lb: oldest possible age of rock unit given some chronostratigraphic scale (use 500 for 500 Ma);
##	ma_ub; oldest possible age of rock unit given some chronostratigraphic scale (use 495 for 495 Ma);
# external_zone_database: name of external database giving zones associated with rock units.
# collections <- strophomenoid_collections;
# edited 2020-03-05
# edited 2020-04-09
# edited 2020-04-12
# edited 2020-05-20
# edited 2020-05-22
# edited 2021-01-27
# edited 2021-03-01
refine_collection_dates_with_external_database <- function(study="",collections,rock_database,zone_database,rock_to_zone_database,time_scale,directory="",save_files=TRUE,output_type=".csv")	{
collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
n_rocks <- nrow(rock_database);
oldest_rocks <- max(collections$max_ma)+50;
youngest_rocks <- min(collections$min_ma)-50;
wagner_rocks <- subset(rock_database,rock_database$ma_lb>=oldest_rocks & rock_database$ma_ub<=youngest_rocks);
redone_collections <- redone_collections_db <- match_paleodb_collections_to_external_stratigraphic_database(collections,wagner_rocks=wagner_rocks);
n_colls <- nrow(redone_collections);
#lost_time <- (1:ncoll)[!redone_collections$early_interval %in% timey_wimey$interval]
#redone_collections$early_interval[lost_time]
#(1:ncolls)[redone_collections$ma_lb<=redone_collections$ma_ub]
#cbind(redone_collections$ma_lb[fff],redone_collections$ma_ub[fff])
#### But wait, there's more!  The PaleoDB has zone data for some localities.  This can further restrict dates
zone_thesaurus <- accersi_zone_thesaurus(zone_database);
zd <- nrow(zone_database);
zt <- nrow(zone_thesaurus);
rzd <- nrow(rock_to_zone_database);

relv_intervals <- sort(unique(c(unique(redone_collections$interval_lb),unique(redone_collections$interval_ub))));
ttt <- match(relv_intervals,time_scale$interval);
relv_time_scale <- time_scale[ttt[!is.na(ttt)],];
#relv_time_scale <- relv_time_scale[!is.na(relv_time_scale),];
relv_time_scale <- relv_time_scale[order(-abs(relv_time_scale$ma_lb),-abs(relv_time_scale$ma_ub)),];

collections_w_zones <- (1:n_colls)[redone_collections$zone!=""];
cwz <- length(collections_w_zones);
redone_collections$zone <- as.character(redone_collections$zone)
zone <- as.character(redone_collections$zone[collections_w_zones]);
redone_collections$zone[collections_w_zones] <- sapply(zone,mundus_zone);
zone_species <- as.character(redone_collections$zone);
non_taxon_zone <- non_taxon_zone_label <- array("",dim=nrow(redone_collections));
zone_species[collections_w_zones] <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
non_taxon_zone_info <- raster::t(sapply(zone,aparecium_nontaxon_zone));
rownames(non_taxon_zone_info) <- NULL;
non_taxon_zone[collections_w_zones] <- non_taxon_zone_info[,1];
non_taxon_zone_label[collections_w_zones] <- non_taxon_zone_info[,2];
#redone_collections <- cbind(redone_collections,zone_species,non_taxon_zone,non_taxon_zone_label);
redone_collections$zone_species <- as.character(zone_species);
redone_collections$non_taxon_zone <- as.character(non_taxon_zone);
redone_collections$non_taxon_zone_label <- as.character(non_taxon_zone_label);

collections_w_zones <- (1:n_colls)[redone_collections$zone!=""];
cwz <- length(collections_w_zones);
zones_matched <- c();

cz <- 0;
while (cz < cwz)	{
	### blows up  when a stage is  entered as the zone!!!
	cz <- cz+1;
	pbdbc <- collections_w_zones[cz];	#redone_collections$collection_no[paldbc]
	zone_paleodb <- mundus_zone(zone=as.character(redone_collections$zone[pbdbc]));
	if (!is.na(match(zone_paleodb,zone_database$zone)))	{
		zone_paleodb_sr <- zone_database$zone_sr[match(zone_paleodb,zone_database$zone)];
		}	else	{
		zone_paleodb_sr <- "flibberty_gibberty";
		}
	zone_spc_paleodb <- as.character(redone_collections$zone_species[pbdbc]);
	zone_nontax_paleodb <- as.character(redone_collections$non_taxon_zone[pbdbc]);
	zone_nontax_lab_paleodb <- as.character(redone_collections$non_taxon_zone_label[pbdbc]);
	# poss_matches gives possible rock matches
	poss_matches <- (1:rzd)[rock_to_zone_database$rock_no_sr %in% redone_collections$rock_no_sr[pbdbc]];
	wagner_rock_unit <- match(redone_collections$rock_no_sr[pbdbc],rock_database$rock_no);

	### if we cannot find formation (member) combination, then just look for formation
	if (length(poss_matches)==0 || (redone_collections$rock_no_sr[pbdbc]>0 && rock_database$member[wagner_rock_unit]==""))
		poss_matches <- (1:rzd)[rock_to_zone_database$formation_no %in% redone_collections$formation_no[pbdbc]];
	#rock_to_zone_database[poss_matches,]
	# now, try to find this zone in zones associated with rock-unit
	matches <- c();
#	relevant_zones <- array("",dim=c(0,5));
#	relevant_zones <- data.frame(relevant_zones);
	relevant_zones <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
								 ma_ub=as.numeric(),interval_lb=as.character(),
								 interval_ub=as.character(),stringsAsFactors=F);
	if (length(poss_matches)>0)	{
		# match listed zones to PaleoDB entered zones
		matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_paleodb];
		# match senior names of zones to PaleoDB entered zones
		matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb];
		matches3 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb_sr];
		# match PaleoDB zone to species name alone
		matches4 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
		# match PaleoDB zone species name only to species name alone
		# match PaleoDB zone species name only to senior species name alone
		if (length(zone_spc_paleodb)>0 && zone_spc_paleodb!="")	{
			matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_spc_paleodb];
			matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_spc_paleodb]
			} else	{
			matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
			matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_paleodb];
			}
		# match PaleoDB zone species name only to senior species name alone
		matches7 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_nontax_paleodb];
		# match PaleoDB zone species name only to senior species name alone
		if (length(zone_nontax_paleodb)>0 && zone_nontax_paleodb != "")	{
			matches8 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_nontax_paleodb];
			matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_nontax_paleodb];
			} else	{
			matches8 <- poss_matches[rock_to_zone_database$non_taxon_zone[poss_matches] %in% zone_paleodb];
			matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_paleodb];
			}
		# match PaleoDB zone name only to a zone label
		matches <- sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7,matches8,matches9)));
		if (length(matches)>0)	{
			these_zones <- data.frame(zone_sr=as.character(rock_to_zone_database$zone_sr[matches]),ma_lb=as.numeric(rock_to_zone_database$ma_lb[matches]),
								 ma_ub=as.numeric(rock_to_zone_database$ma_ub[matches]),interval_lb=as.character(rock_to_zone_database$interval_lb[matches]),
								 interval_ub=as.character(rock_to_zone_database$interval_ub[matches]),stringsAsFactors=F);
			relevant_zones <- rbind(relevant_zones,these_zones);
			}
		}

	## Look out for multiple zones: if that is the caes, then get ages from all of them.
	if (length(poss_matches)>0 && length(matches)==0)	{
#		print(paste("couldn't find",as.character(zone_spc_paleodb)));
		multizones <- divido_zone(zone=zone_paleodb);
		mz <- length(multizones);
		zz <- 0;
		while (zz<mz)	{
#		if (mz>1)	{
#			multimatches <- c();
#			for (zz in 1:mz)	{
				# match whole name of zones to PaleoDB entered zone
			zz <- zz+1;
			multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
				# match senior names of zones to PaleoDB entered zones
			if (length(multimatches)==0)	{
				multimatches <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% multizones[zz]];
				# match species name only to PaleoDB entered zone
				if (length(multimatches)==0)	{
					multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
					# match senior synonym species name only to PaleoDB entered zone
					if (length(multimatches)==0)	{
						multimatches <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% multizones[zz]];
						if (length(multimatches)==0)	{
							multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
							if (length(multimatches)==0)	{
								multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
								if (length(multimatches)==0)	{
									multimatches <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% multizones[zz]];
									}
								}
							}
						}
					}
				}
			if (length(multimatches)>0)	{
				these_zones <- data.frame(zone_sr=as.character(rock_to_zone_database$zone_sr[multimatches]),ma_lb=as.numeric(rock_to_zone_database$ma_lb[multimatches]),
										  ma_ub=as.numeric(rock_to_zone_database$ma_ub[multimatches]),interval_lb=as.character(rock_to_zone_database$interval_lb[multimatches]),
										  interval_ub=as.character(rock_to_zone_database$interval_ub[multimatches]),stringsAsFactors=F);
				relevant_zones <- rbind(relevant_zones,these_zones);
				} else if (length(multimatches)==0)	{
#				mm1 <- (1:zt)[zone_thesaurus$zone %in% multizones[zz]];
#				mm2 <- (1:zt)[zone_thesaurus$zone_species %in% multizones[zz]];
#				mm3 <- (1:zt)[zone_thesaurus$zone_species_sr %in% multizones[zz]];
#				mm4 <- (1:zt)[zone_thesaurus$non_taxon_zone %in% multizones[zz]];
#				mm5 <- (1:zt)[zone_thesaurus$non_taxon_zone_label %in% multizones[zz]];
#				mm6 <- (1:zt)[zone_thesaurus$non_taxon_zone_sr %in% multizones[zz]];
#				mm7 <- (1:zt)[zone_thesaurus$non_taxon_zone_label_sr %in% multizones[zz]];
				multimatches <- sort(unique(which(zone_thesaurus==multizones[zz],arr.ind = T)[,1]));
#				multimatches <- (1:zt)[zone_thesaurus$zone_sr[zone_thesaurus$zone_species %in% multizones[zz]] %in% rock_to_zone_data$zone[poss_matches]];
#				multimatches <- sort(unique(c(mm1,mm2,mm3,mm4,mm5,mm6,mm7)));
				if (length(multimatches)>0)	{
					zone_taxa <- sort(unique(zone_thesaurus$zone_sr[multimatches]));
					poss_zones_all <- (1:zd)[zone_database$zone_sr %in% zone_taxa];
					rel_z <- cbind(as.character(zone_database$zone_sr[poss_zones_all]),as.numeric(zone_database$ma_lb[poss_zones_all]),as.numeric(zone_database$ma_ub[poss_zones_all]),as.character(zone_database$interval_lb[poss_zones_all]),as.character(zone_database$interval_ub[poss_zones_all]));
					for (ztx  in 1:length(zone_taxa))	{
						poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
						zlb <- max(as.numeric(poss_zones[,2]));
						zub <- min(as.numeric(poss_zones[,3]));
						if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])	{
							##HERE!!
#							relevant_zones <- rbind(relevant_zones,poss_zones);
							these_zones <- data.frame(zone_sr=as.character(zone_database$zone_sr[multimatches]),ma_lb=as.numeric(zone_database$ma_lb[multimatches]),
													  ma_ub=as.numeric(zone_database$ma_ub[multimatches]),interval_lb=as.character(zone_database$interval_lb[multimatches]),
													  interval_ub=as.character(zone_database$interval_ub[multimatches]),stringsAsFactors=F);
							relevant_zones <- rbind(relevant_zones,these_zones);
							}
						}
					}
				}
			if (length(multimatches)>0)	matches <- c(matches,multimatches);
			}
		}
	### routine to find if a senior synonym of the reported zone is in the database
#		if (length(matches)==0 && length(strsplit(zone_paleodb," ")[[1]])==1)	{
	## Look out for zones that are not "Amorphognathus ordovicicus" but "Zone D"
	if (length(poss_matches)>0 && length(matches)==0)	{
		aa <- (1:zd)[zone_database$zone %in% zone_paleodb];
		if (length(aa)==0)	{
			bb <- (1:zd)[zone_database$zone_species %in% zone_paleodb];
			} else	{
			bb <- c();
			}
		cc <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_paleodb];
		dd <- (1:zd)[zone_database$non_taxon_zone_label_sr %in% zone_paleodb];
		ee <- sort(unique(c(aa,bb,cc,dd)));
		if (length(ee)>0)	{
			try_this_zone <- unique(zone_database$zone_sr[ee]);
			try_this_zone_spc <- transmogrify_full_zone_names_to_species_names_only(try_this_zone);
			## now, repeat the matching game!
			matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% try_this_zone];
			# match senior names of zones to PaleoDB entered zones
			matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% try_this_zone];
			matches3 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% try_this_zone_spc];
			matches <- c(matches,sort(unique(c(matches1,matches2,matches3))));
			if (length(matches)==0)	{
				these_zones <- data.frame(zone_sr=as.character(zone_database$zone_sr[ee]),ma_lb=as.numeric(zone_database$ma_lb[ee]),
										  ma_ub=as.numeric(zone_database$ma_ub[ee]),interval_lb=as.character(zone_database$interval_lb[ee]),
										  interval_ub=as.character(zone_database$interval_ub[ee]),stringsAsFactors=F);
				relevant_zones <- rbind(relevant_zones,these_zones);
#				relevant_zones <- rbind(relevant_zones,cbind(zone_database$zone_sr[ee],zone_database$ma_lb[ee],zone_database$ma_ub[ee],as.character(zone_database$interval_lb[ee]),as.character(zone_database$interval_ub[ee])));
				matches <- length(zone_database$zone_sr[ee])
				}
			}
		}

	if (length(poss_matches)==0)	{
		multizones <- divido_zone(zone=zone_paleodb);
		if (length(multizones)>1)
			zone_paleodb <- multizones;
		matches1 <- (1:zd)[zone_database$zone %in% zone_paleodb];
		if (length(matches1)==0)	{
			matches2 <- (1:zd)[zone_database$zone_species %in% zone_spc_paleodb];
			matches2 <- matches2[!zone_database$zone_species[matches2]==""];
			} else	{
			matches2 <- c();
			}
		if (zone_nontax_paleodb!="")	{
			matches3 <- (1:zd)[zone_database$non_taxon_zone %in% zone_nontax_paleodb];
			} else	{
			matches3 <- c();
			}
		if (zone_nontax_lab_paleodb!="")	{
			matches4 <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_nontax_lab_paleodb];
			} else	{
			matches4 <- c();
			}
		matches <- sort(unique(c(matches1,matches2,matches3,matches4)));
		if (length(matches)==0)	{
			multizones <- divido_zone(zone_paleodb);
			mz <- length(multizones);
			if (mz>1)	{
				rel_z <- c();
	#			multimatches <- c();
				for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
					matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
					matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
					# match species names only of zones to PaleoDB entered zones
					matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
					matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
					multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
					if (length(multimatches)>0)	{
						rel_z <- rbind(rel_z,cbind(as.character(zone_database$zone_sr[multimatches]),as.numeric(zone_database$ma_lb[multimatches]),as.numeric(zone_database$ma_ub[multimatches]),as.character(zone_database$interval_lb[multimatches]),as.character(zone_database$interval_ub[multimatches])));
						matches <- c(matches,multimatches);
						}
					}
			} else	{
				rel_z <- c();
				}
			} else	{
			rel_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[matches]),ma_lb=as.numeric(zone_database$ma_lb[matches]),ma_ub=as.numeric(zone_database$ma_ub[matches]),interval_lb=as.character(zone_database$interval_lb[matches]),interval_ub=as.character(zone_database$interval_ub[matches]),stringsAsFactors = F);
			}
		if (length(matches)>0)	{
			zone_taxa <- sort(unique(rel_z$zone_sr));
		# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
		#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
			if (length(zone_taxa)>=1)	{
#			possible_zones <- c();
				for (ztx  in 1:length(zone_taxa))	{
					poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
					zlb <- max(as.numeric(poss_zones[,2]));
					zub <- min(as.numeric(poss_zones[,3]));
					if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])
						relevant_zones <- rbind(relevant_zones,poss_zones);
					}
				}
			}
		}	# end case where no rock unit is reported

	# now redate the collection based on zone or zones
	if (length(relevant_zones$zone_sr)==0)	{
		match1 <- (1:zd)[zone_database$zone %in% redone_collections$zone[pbdbc]];
		match2 <- (1:zd)[zone_database$zone_sr %in% redone_collections$zone[pbdbc]];
		if (redone_collections$zone_species[pbdbc]!="")	{
			match3 <- (1:zd)[zone_database$zone_species %in% redone_collections$zone_species[pbdbc]];
			match4 <- (1:zd)[zone_database$zone_species_sr %in% redone_collections$zone_species[pbdbc]];
			}	else	{
			match3 <- c();
			match4 <- c();
#			match3 <- match3[!zone_data$zone_species[match3] %in% ""];
			}
		matches <- sort(unique(c(match1,match2,match3,match4)));
		if (length(matches)==0)	{
			multizones <- divido_zone(zone_paleodb);
			mz <- length(multizones);
			if (mz>1)	{
				rel_z <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
								ma_ub=as.numeric(),interval_lb=as.character(),
								interval_ub=as.character(),stringsAsFactors = F);
	#			multimatches <- c();
				for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
					matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
					matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
					# match species names only of zones to PaleoDB entered zones
					matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
					matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
					multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
					if (length(multimatches)>0)	{
						nwr_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[multimatches]),ma_lb=as.numeric(zone_database$ma_lb[multimatches]),
											ma_ub=as.numeric(zone_database$ma_ub[multimatches]),interval_lb=as.character(zone_database$interval_lb[multimatches]),
											interval_ub=as.character(zone_database$interval_ub[multimatches]),stringsAsFactors = F);

						rel_z <- rbind(rel_z,nwr_z);
						matches <- c(matches,multimatches);
						}
					}
				} else	{
				rel_z <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
								ma_ub=as.numeric(),interval_lb=as.character(),
								interval_ub=as.character(),stringsAsFactors = F);
				}
			} else	{
			rel_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[matches]),
								ma_lb=as.numeric(zone_database$ma_lb[matches]),
								ma_ub=as.numeric(zone_database$ma_ub[matches]),
								interval_lb=as.character(zone_database$interval_lb[matches]),
								interval_ub=as.character(zone_database$interval_ub[matches]),
								stringsAsFactors = F);
			rel_z <- unique(rel_z);
#			rel_z <- cbind(as.character(zone_database$zone_sr[matches]),as.numeric(zone_database$ma_lb[matches]),as.numeric(zone_database$ma_ub[matches]),as.character(zone_database$interval_lb[matches]),as.character(zone_database$interval_ub[matches]));
			}
		if (length(matches)>0)	{
			zone_taxa <- sort(unique(rel_z$zone_sr));
		# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
		#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
			if (length(zone_taxa)>=1)	{
#			possible_zones <- c();
				for (ztx  in 1:length(zone_taxa))	{
					poss_zones <- subset(rel_z,rel_z$zone_sr==zone_taxa[ztx]);
					zlb <- max(as.numeric(poss_zones$ma_lb));
					zub <- min(as.numeric(poss_zones$ma_ub));
					if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])
						relevant_zones <- rbind(relevant_zones,poss_zones);
					}
				}
			}
		}

	# now, make adjustments
#	colnames(relevant_zones) <- c("zone_sr","ma_lb","ma_ub","interval_lb","interval_ub");
	if (nrow(relevant_zones)>0)	{
		zones_matched <- c(zones_matched,pbdbc);
		# if zone starts after earliest possible age for formation, adjust oldest age & stage
		ma_lb <- as.numeric(as.character(relevant_zones$ma_lb));
		ma_mx <- max(ma_lb);
		ma_ub <- as.numeric(as.character(relevant_zones$ma_ub));
		ma_mn <- min(ma_ub);
#		lbz <- match(ma_mx,ma_lb);
#		accersi_temporal_overlap(lb1=100,lb2=50,ub1=75,ub2=40)
		new_dates <- accersi_temporal_overlap(lb1=as.numeric(redone_collections$ma_lb[pbdbc]),lb2=ma_mx,ub1=as.numeric(redone_collections$ma_ub[pbdbc]),ub2=ma_mn);
		if (new_dates[1]!=0)	{
			redone_collections$ma_lb[pbdbc] <- as.numeric(new_dates[1]);
			redone_collections$ma_ub[pbdbc] <- as.numeric(new_dates[2]);
			redone_collections$interval_lb[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>=round(as.numeric(new_dates[1]),2))]);
			redone_collections$interval_ub[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>round(as.numeric(new_dates[2]),2))]);
			} else	{
			if (as.numeric(redone_collections$ref_pubyr[pbdbc])>=1995)	{
				redone_collections$ma_lb[pbdbc] <- ma_mx;
				redone_collections$ma_ub[pbdbc] <- ma_mn;
				redone_collections$interval_lb[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>=round(ma_mx,2))]);
				redone_collections$interval_ub[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>round(ma_mn,2))]);
				}
			}
#		if (redone_collections$ma_lb[paldbc] > ma_mx || redone_collections$ma_lb[paldbc]==0)	{
#			redone_collections$ma_lb[paldbc] <- ma_mx;
#			redone_collections$interval_lb[paldbc] <- as.character(relevant_zones$interval_lb[lbz]);
#			}
		# if zone ends before latest possible age for formation, adjust youngest age & stage
#		ma_ub <- as.numeric(as.character(relevant_zones$ma_ub));
#		ma_mn <- min(ma_ub);
#		ubz <- match(ma_mn,ma_ub);
#		if (redone_collections$ma_ub[paldbc] < ma_mn)	{
#			redone_collections$ma_ub[paldbc] <- ma_mn;
#			redone_collections$interval_ub[paldbc] <- as.character(relevant_zones$interval_ub[ubz]);
#			}
		}
#	print(zones_matched);
	}

# check for poorly restricted rocks that had better ages from Paleodb
ma_max <- max(collections$max_ma)+25;
ma_min <- min(collections$min_ma)-25;
timey_wimey <- subset(time_scale,time_scale$ma_lb<=ma_max);
timey_wimey <- subset(timey_wimey,timey_wimey$ma_ub>=ma_min);
stage <- timey_wimey$interval;
timey_wimey$interval <- sapply(stage,mundus_stage);

disjunct <- c();
for (cn in 1:n_colls)	{
	# if there is overlap, then reduce
	if (redone_collections$max_ma[cn]>redone_collections$ma_ub[cn] && redone_collections$ma_lb[cn]>redone_collections$min_ma[cn])	{
		wibbly_wobbly <- accersi_temporal_overlap(lb1=redone_collections$max_ma[cn],ub1=redone_collections$min_ma[cn],lb2=redone_collections$ma_lb[cn],ub2=redone_collections$ma_ub[cn]);
		if (wibbly_wobbly$ma_lb==redone_collections$max_ma[cn])
			redone_collections$interval_lb[cn] <- redone_collections$early_interval[cn];
		if (wibbly_wobbly$ma_ub==redone_collections$min_ma[cn])
			redone_collections$interval_ub[cn] <- redone_collections$late_interval[cn];
		redone_collections$ma_lb[cn] <- wibbly_wobbly$ma_lb;
		redone_collections$ma_ub[cn] <- wibbly_wobbly$ma_ub;
		} else	{
		disjunct <- c(disjunct,cn);
		}
	}
#ddd <- (1:ncoll)[redone_collections$ma_lb<=redone_collections$ma_ub]
# cull out any collections using time units not in the external time scale
#stage_numb_1 <- match(redone_collections$interval_lb,timey_wimey$interval);
#stage_numb_2 <- match(redone_collections$interval_ub,timey_wimey$interval);
stage_numb_1 <- match(collections$early_interval,timey_wimey$interval);
stage_numb_2 <- match(collections$late_interval,timey_wimey$interval);
cut_these <- (1:n_colls)[is.na(stage_numb_1)];
#redone_collections[cut_these[1],]
#collections$early_interval[cut_these]
if (length(cut_these) > 0)	{
	rd_coll <- nrow(redone_collections);
	pass <- (1:rd_coll)[!is.na(stage_numb_1)];
	rdc <- redone_collections[pass,];
	stage_numb_1 <- stage_numb_1[pass];
	stage_numb_2 <- stage_numb_2[pass];
	# fixed 2021-01-27
	early_interval <- collections$early_interval[pass];
	# fixed 2021-01-27
	late_interval <- collections$late_interval[pass];
#	double_check_extdb <- double_check_extdb[pass];
	rd_coll <- nrow(rdc);
	} else	{
	rd_coll <- n_colls;
	}
cut_these <- (1:rd_coll)[is.na(stage_numb_2)];
#collections$late_interval[cut_these];
if (length(cut_these) > 0)	{
	pass <- (1:rd_coll)[!is.na(stage_numb_2)];
	rdc <- redone_collections[pass,];
	stage_numb_1 <- stage_numb_1[pass];
	stage_numb_2 <- stage_numb_2[pass];
	early_interval <- early_interval[pass];
	late_interval <- late_interval[pass];
#	double_check_extdb <- double_check_extdb[pass];
	rd_coll <- nrow(rdc);
	}

missing <- collections_w_zones[!collections_w_zones %in% zones_matched];
untraceable_zone_information <- cbind(redone_collections$collection_no[missing],redone_collections$formation[missing],redone_collections$member[missing],redone_collections$zone[missing]);
colnames(untraceable_zone_information) <- c("collection_no","formation","member","zone");
untraceable_zone_information <- data.frame(untraceable_zone_information);

age <- redone_collections$ma_lb;
redone_collections$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = relv_time_scale);
age <- redone_collections$ma_ub;
redone_collections$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = relv_time_scale);
# 2020-05-22: Cases where rock units have much wider chronostratigraphic ranges than the PaleoDB assignments.
#	Presumably, there was a reason for the more restricted age: ergo, use that instead;
#n_colls <- nrow(redone_collections);
#pbdb_ranges_are_subsets <- (1:n_colls)[redone_collections$max_ma<redone_collections$ma_lb][(1:n_colls)[redone_collections$max_ma<redone_collections$ma_lb] %in% (1:n_colls)[redone_collections$min_ma>redone_collections$ma_ub]];
#redone_collections[pbdb_ranges_are_subsets,]
#redone_collections[match(redone_collections$collection_no[pbdb_ranges_are_subsets],pbdb_collections$collection_no),]

# 2019-03-07: routine to get unidentified rock units limited....
#ranges_paledb <- redone_collections$max_ma-redone_collections$min_ma;
#ranges_redone <- redone_collections$ma_lb-redone_collections$ma_ub;
#max_ma_downloaded <- max(redone_collections$max_ma);
#min_ma_downloaded <- min(redone_collections$min_ma);
#too_old <- (1:ncoll)[redone_collections$ma_lb > max_ma_downloaded];
#redone_collections$ma_lb[too_old] <- max_ma_downloaded;
#redone_collections$interval_lb[match(max_ma_downloaded,redone_collections$max_ma)]
#redone_collections$interval_lb
#too_young <- (1:ncoll)[redone_collections$ma_ub < min_ma_downloaded];

if (save_files)	{
	output1 <- paste(study,"_Recalibrated_Collections",output_type,sep="");
	output2 <- paste(study,"_Collection_Zone_Data_Untraced",output_type,sep="");
	if (directory!="")	{
		output1 <- paste(directory,output1,sep="");
		output2 <- paste(directory,output2,sep="");
		}
	if (output_type==".csv")	{
		write.csv(redone_collections,output1,row.names=FALSE);
		write.csv(untraceable_zone_information,output2,row.names=FALSE);
		} else	{
		write.table(redone_collections,output1,row.names=FALSE,col.names=TRUE,sep="\t");
		write.table(untraceable_zone_information,output2,row.names=FALSE,col.names=TRUE,sep="\t");
		}
	}

output <- list(redone_collections,untraceable_zone_information);
names(output) <- c("Recalibrated_Collections","Untraced_Collection_Zones");

return(output);
}

# routine to organize external rock unit database to facilitate matching to Paleobiology Database rock units
organize_rock_data_for_paleodb_matching	<- function(external_strat_database)	{
if (external_strat_database!="")	{
	file_type <- revelare_file_type(filename = external_strat_database);
	if (file_type=="csv")	{
		wagner_rocks <- read.csv(file=external_strat_database,header=TRUE,stringsAsFactors=FALSE,encoding="UTF-8");
		}	else	{
		wagner_rocks <- read.table(file=rock_unit_database,header=TRUE,stringsAsFactors=FALSE,sep="\t");
		}
	wagner_rocks$rock_no_sr <- as.numeric(wagner_rocks$rock_no_sr);
	wagner_rocks$rock_no <- as.numeric(wagner_rocks$rock_no);
	wagner_rocks$formation_no <- as.numeric(wagner_rocks$formation_no);
	wagner_rocks$ma_lb <- as.numeric(wagner_rocks$ma_lb);
	wagner_rocks$ma_ub <- as.numeric(wagner_rocks$ma_ub);
	}	else	{
	wagner_rocks <- rock_unit_data$rock_unit_database;
	}
n_rocks <- nrow(wagner_rocks);
#### shift formation & rock unit entries for group-only records to "•"
wagner_rocks$formation_no <- expello_na_from_vector(wagner_rocks$formation_no,0)
groups_only <- subset(wagner_rocks,wagner_rocks$group!="");
groups_only <- subset(groups_only,groups_only$formation=="");
groups_only <- subset(groups_only,groups_only$member=="");
groups_only_nos <- as.numeric(rownames(groups_only));
wagner_rocks$formation_no[groups_only_nos] <- wagner_rocks$rock_no[groups_only_nos];
wagner_rocks$group <- as.character(wagner_rocks$group);
wagner_rocks$formation <- as.character(wagner_rocks$formation);
wagner_rocks$member <- as.character(wagner_rocks$member);
wagner_rocks$full_name <- as.character(wagner_rocks$full_name);
# added 2021-04-03
wagner_rocks$full_name[groups_only_nos] <- paste("[",wagner_rocks$group[groups_only_nos],"]",sep="");
wagner_rocks$formation[groups_only_nos] <- wagner_rocks$member[groups_only_nos] <- wagner_rocks$group[(1:n_rocks)[wagner_rocks$group==""]] <- "";
#wagner_rocks$formation[groups_only_nos] <- "•"
#wagner_rocks$full_name[groups_only_nos] <- "§";
#wagner_rocks$member[groups_only_nos] <- "¢";
#wagner_rocks$group[(1:n_rocks)[wagner_rocks$group==""]] <- "£";

named_rock_unit <- as.character(wagner_rocks$formation);
print("Cleaning Formation Names...");
formation_clean_basic <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
formation_clean_no_rock <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
formation_clean_no_rock_formal <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);
member_clean_basic <- as.character(wagner_rocks$member);
member_ided <- (1:n_rocks)[wagner_rocks$member!=""];
named_rock_unit <- as.character(member_clean_basic[member_ided]);
#xxx <- c();
#for (nru in 1:length(named_rock_unit))
#	xxx <- c(xxx,mundify_rock_unit_names(named_rock_unit = named_rock_unit[nru],dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE))
print("Cleaning Member Names...");
member_clean_basic[member_ided] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
member_clean_no_rock <- member_clean_no_rock_formal <- member_clean_basic;
member_clean_no_rock[member_ided] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
member_clean_no_rock_formal[member_ided] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

named_rock_unit <- wagner_rocks$group;
print("Cleaning Group Names...");
group_clean_basic <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
group_clean_no_rock <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
group_clean_no_rock_formal <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

### suppose we starte with formation = "Burgess Shale" & member = "Lower)
####  rock_unit_clean_basic: "Burgess Shale (Lower)"
####  rock_unit_clean_no_rock: "Burgess (Lower)"
####  rock_unit_clean_no_rock_formal: "Burgess"
rock_unit_clean_basic <- formation_clean_basic;
member_ided <- (1:n_rocks)[member_clean_basic!=""];
rock_unit_clean_basic[member_ided] <- paste(formation_clean_basic[member_ided]," (",member_clean_basic[member_ided],")",sep="");
rock_unit_clean_no_rock <- formation_clean_no_rock;
member_ided <- (1:n_rocks)[member_clean_no_rock!=""];
rock_unit_clean_no_rock[member_ided] <- paste(formation_clean_no_rock[member_ided]," (",member_clean_no_rock[member_ided],")",sep="");
rock_unit_clean_no_rock_formal <- formation_clean_no_rock_formal;
member_ided <- (1:n_rocks)[member_clean_no_rock_formal!=""];
rock_unit_clean_no_rock_formal[member_ided] <- paste(formation_clean_no_rock_formal[member_ided]," (",member_clean_no_rock_formal[member_ided],")",sep="");
rock_unit_senior <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr,wagner_rocks$rock_no)];
organized_rocks <- data.frame(cbind(rock_unit_senior,formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal),stringsAsFactors = FALSE);
wagner_rocks <- cbind(wagner_rocks,organized_rocks);
#wagner_rocks <- cbind(wagner_rocks,as.character(rock_unit_senior),formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal);
return(wagner_rocks);
}

# routine to organize Paleobiology Database rock units to facilitate matching to external rock unit database
# modified 2020-02-15: add groups that are not uniquely specified
# modified 2020-02-25: add formations that are not uniquely specified
# modified 2020-03-09
accersi_rock_unit_thesaurus_given_paleodb_collections	<- function(paleodb_rocks)	{
# standardize group usage (these often are not entered)
if (is.na(match("group",colnames(paleodb_rocks))) && !is.na(match("stratgroup",colnames(paleodb_rocks))))
	colnames(paleodb_rocks)[match("stratgroup",colnames(paleodb_rocks))] <- "group";
formation_names <- unique(paleodb_rocks$formation[paleodb_rocks$formation!=""]);
for (fn in 1:length(formation_names))	{
	named_rock_unit <- paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]];
	named_rock_unit <- unique(named_rock_unit[named_rock_unit!=""]);
	if (length(named_rock_unit)==1)	{
		paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]] <- named_rock_unit;
		} else if (length(named_rock_unit)>1)	{
		formation_group <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);
		examples <- c();
		groupies <- paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]];
		for (fg in 1:length(formation_group))
			examples <- c(examples,sum(groupies==formation_group[fg]));
		paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]] <- unique(named_rock_unit[named_rock_unit!=""])[match(max(examples),examples)]
		}
	}
paleodb_rocks <- unique(paleodb_rocks);

n_rocks <- nrow(paleodb_rocks);
rock_no <- 1:n_rocks;
formation_no <- match(paleodb_rocks$formation,formation_names);
formation_no[is.na(formation_no)] <- 0;

#formation_no[(1:n_rocks)[is.na(formation_no)]] <- rock_no[(1:n_rocks)[is.na(formation_no)]];

#### shift formation & rock unit entries for group-only records to "•"
has_formation <- (1:n_rocks)[paleodb_rocks$formation!=""];
has_member <- (1:n_rocks)[paleodb_rocks$member!=""];
has_group <- (1:n_rocks)[paleodb_rocks$group!=""];
group_only <- has_group[!has_group %in% c(has_formation,has_member)];
formation_only <- has_formation[!has_formation %in% has_member];
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];

paleodb_rocks$group <- as.character(paleodb_rocks$group);
paleodb_rocks$formation <- as.character(paleodb_rocks$formation);
paleodb_rocks$member <- as.character(paleodb_rocks$member);
full_name <- paleodb_rocks$formation;
full_name[formation_and_member] <- paste(full_name[formation_and_member]," (",paleodb_rocks$member[formation_and_member],")",sep="");
full_name[member_only] <- paleodb_rocks$member[member_only];
full_name[group_only] <- paleodb_rocks$group[group_only];
formation_no[group_only] <- rock_no[group_only];

formation_clean_basic <- named_rock_unit <- as.character(paleodb_rocks$formation);
formation_clean_basic[has_formation] <- sapply(named_rock_unit[has_formation],mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
formation_clean_no_rock_formal <- formation_clean_no_rock <- formation_clean_basic;
formation_clean_no_rock[has_formation] <- sapply(named_rock_unit[has_formation],mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
formation_clean_no_rock_formal[has_formation] <- sapply(named_rock_unit[has_formation],mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);
member_clean_basic <- as.character(paleodb_rocks$member);
named_rock_unit <- as.character(member_clean_basic[has_member]);
member_clean_basic[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
member_clean_no_rock <- member_clean_no_rock_formal <- member_clean_basic;
member_clean_no_rock[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
member_clean_no_rock_formal[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

group_clean_basic <- paleodb_rocks$group;
named_rock_unit <- paleodb_rocks$group[has_group];
group_clean_basic[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
group_clean_no_rock <- group_clean_no_rock_formal <- group_clean_basic;
group_clean_no_rock[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
group_clean_no_rock_formal[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

### suppose we start with formation = "Burgess Shale" & member = "Lower)
####  rock_unit_clean_basic: "Burgess Shale (Lower)"
####  rock_unit_clean_no_rock: "Burgess (Lower)"
####  rock_unit_clean_no_rock_formal: "Burgess"
rock_unit_clean_basic <- formation_clean_basic;
has_member <- (1:n_rocks)[member_clean_basic!=""];	# some members are erased
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_basic[formation_and_member] <- paste(formation_clean_basic[formation_and_member]," (",member_clean_basic[formation_and_member],")",sep="");
rock_unit_clean_basic[group_only] <- group_clean_basic[group_only];

rock_unit_clean_no_rock <- formation_clean_no_rock;
has_member <- (1:n_rocks)[member_clean_no_rock!=""];	# some members are erased
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_no_rock[formation_and_member] <- paste(formation_clean_no_rock[formation_and_member]," (",member_clean_no_rock[formation_and_member],")",sep="");
rock_unit_clean_no_rock[group_only] <- group_clean_no_rock[group_only];

rock_unit_clean_no_rock_formal <- formation_clean_no_rock_formal;
has_formation <- (1:n_rocks)[formation_clean_no_rock_formal!=""]; # some formations lost
has_member <- (1:n_rocks)[member_clean_no_rock_formal!=""]; 	# some members lost
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_no_rock_formal[formation_and_member] <- paste(formation_clean_no_rock_formal[formation_and_member]," (",member_clean_no_rock_formal[formation_and_member],")",sep="");

rock_no_sr <- rock_no;
m_a_f <- has_member[member_clean_no_rock_formal[has_member] %in% formation_clean_no_rock_formal];
members_as_formations <- member_clean_no_rock_formal[m_a_f];
rock_no_sr[m_a_f] <- match(members_as_formations,formation_clean_no_rock_formal);
g_a_f <- group_only[group_clean_no_rock[group_only] %in% formation_clean_no_rock];
groups_as_formations <- group_clean_no_rock[g_a_f];
rock_no_sr[g_a_f] <- match(groups_as_formations,formation_clean_no_rock);
aaa <- formation_no[match(rock_no_sr,rock_no)]
formation_no_sr <- formation_no[match(rock_no_sr,rock_no)]

#senior_rock_units <- unique(rock_unit_clean_no_rock[rock_unit_clean_no_rock!=""]);
#r_u_c_n_r <- length(rock_unit_clean_no_rock[rock_unit_clean_no_rock!=""]);
#if (length(senior_rock_units) < r_u_c_n_r)	{
#	for (sr in 1:length(length(senior_rock_units)))	{
#		this_senior <- (1:r_u_c_n_r)[senior_rock_units[sr]==rock_unit_clean_no_rock];
#		if (length(this_senior)>1)	print(sr)
#		}
#	}
rock_unit_senior <- full_name[match(rock_no_sr,rock_no)];
organized_rocks <- data.frame(rock_no_sr=as.numeric(rock_no_sr),rock_no=as.numeric(rock_no),formation_no=as.numeric(formation_no),formation_no_sr=as.numeric(formation_no_sr),
			formation=as.character(paleodb_rocks$formation),member=as.character(paleodb_rocks$member),
			full_name=as.character(full_name),group=as.character(paleodb_rocks$group),
			interval_lb=as.character(paleodb_rocks$interval_lb,interval_ub=as.character(paleodb_rocks$interval_ub)),
			ma_lb=as.numeric(paleodb_rocks$ma_lb),ma_ub=as.numeric(paleodb_rocks$ma_ub),
			rock_unit_senior=as.character(rock_unit_senior),formation_clean_basic=as.character(formation_clean_basic),
			formation_clean_no_rock=as.character(formation_clean_no_rock),formation_clean_no_rock_formal=as.character(formation_clean_no_rock_formal),
			member_clean_basic=as.character(member_clean_basic),member_clean_no_rock=as.character(member_clean_no_rock),
			member_clean_no_rock_formal=as.character(member_clean_no_rock_formal),group_clean_basic=as.character(group_clean_basic),
			group_clean_no_rock=as.character(group_clean_no_rock),group_clean_no_rock_formal=as.character(group_clean_no_rock_formal),
			rock_unit_clean_basic=as.character(rock_unit_clean_basic),rock_unit_clean_no_rock=as.character(rock_unit_clean_no_rock),
			rock_unit_clean_no_rock_formal=as.character(rock_unit_clean_no_rock_formal),
			stringsAsFactors=hell_no);
unique_rocks <- sort(unique(rock_no_sr));
for (ur in 1:length(unique_rocks))	{
#	this_rock <- subset(organized_rocks,organized_rocks$rock_no_sr==41);
	syns <- (1:nrow(organized_rocks))[organized_rocks$rock_no_sr==ur];
	if (length(syns)>1)	{
#		print(ur);
		oldest <- match(max(organized_rocks$ma_lb[syns]),organized_rocks$ma_lb[syns]);
		youngest <- match(min(organized_rocks$ma_ub[syns]),organized_rocks$ma_ub[syns]);
		organized_rocks$interval_lb[syns] <- as.character(organized_rocks$interval_lb[syns[oldest]]);
		organized_rocks$interval_ub[syns] <- as.character(organized_rocks$interval_ub[syns[oldest]]);
		organized_rocks$ma_lb[syns] <- as.numeric(organized_rocks$ma_lb[syns[oldest]]);
		organized_rocks$ma_ub[syns] <- as.numeric(organized_rocks$ma_ub[syns[oldest]]);

		if (length(unique(organized_rocks$group[syns]))>1)	{
			stratgroup <- unique(organized_rocks$group[syns]);
			stratgroup <- stratgroup[stratgroup!=""];
			if (length(stratgroup)==1)
				organized_rocks$group[syns] <- stratgroup;
			}
		}
	}

final_formations <- unique(organized_rocks$formation);
final_formations <- final_formations[final_formations!=""];
for (ff in 1:length(final_formations))	{
	this_form <- subset(organized_rocks,organized_rocks$formation==final_formations[ff]);
	if (sum(this_form$member=="")==0)	{
		dummy_rock <- this_form;
		dummy_rock$member <- dummy_rock$member_clean_basic <- dummy_rock$member_clean_no_rock <- dummy_rock$member_clean_no_rock_formal <- "";
		dummy_rock$full_name <- dummy_rock$formation;
		dummy_rock$rock_no_sr <- dummy_rock$formation_no_sr;
		dummy_rock$rock_no <- dummy_rock$formation_no;
		dummy_rock$rock_unit_clean_basic <- mundify_rock_unit_names(dummy_rock$formation);
		dummy_rock$rock_unit_clean_no_rock <- mundify_rock_unit_names(dummy_rock$formation,delete_rock_type = T);
		dummy_rock$rock_unit_clean_no_rock_formal <- mundify_rock_unit_names(dummy_rock$formation,delete_rock_type = T,delete_informal = T);
		organized_rocks <- rbind(organized_rocks,dummy_rock);
		organized_rocks <- organized_rocks[order(organized_rocks$formation_no_sr,organized_rocks$formation_no,organized_rocks$rock_no_sr,organized_rocks$rock_no),];
		}
	}

final_groups <- unique(organized_rocks$group);
final_groups <- final_groups[final_groups!=""];
for (fg in 1:length(final_groups))	{
	this_group <- subset(organized_rocks,organized_rocks$group==final_groups[fg]);
	this_group_a <- subset(this_group,this_group$formation=="");
	this_group_b <- subset(this_group_a,this_group_a$member=="");
	if (nrow(this_group_b)==0)	{
		new_group <- this_group[1,];
		new_group$rock_no <- new_group$rock_no_sr <- new_group$formation_no <- new_group$formation_no_sr <- nrow(this_group)+1;
#		new_group$full_name <- cleaned_rocks$paleodb_clean_group_basic[coll_no];
		new_group$full_name <- final_groups[fg];
		new_group$formation <- new_group$member <- new_group$formation_clean_basic <- new_group$formation_clean_no_rock <- new_group$formation_clean_no_rock_formal <- new_group$member_clean_basic <- new_group$member_clean_no_rock <- new_group$member_clean_no_rock_formal <- "";
		new_group$rock_unit_clean_basic <- mundify_rock_unit_names(new_group$group,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
		new_group$rock_unit_clean_no_rock <- mundify_rock_unit_names(new_group$group,dehyphenate=TRUE,delete_rock_type = T,delete_informal = FALSE);
		new_group$rock_unit_clean_no_rock_formal <- mundify_rock_unit_names(new_group$group,dehyphenate=TRUE,delete_rock_type = T,delete_informal = T);
		new_group$ma_lb <- max(abs(this_group$ma_lb));
		new_group$ma_ub <- min(abs(this_group$ma_ub));
		new_group$interval_lb <- this_group$interval_lb[match(new_group$ma_lb,this_group$ma_lb)];
		new_group$interval_ub <- this_group$interval_ub[match(new_group$ma_ub,this_group$ma_ub)];
		new_group$rock_no_sr <- new_group$rock_no <- new_group$formation_no <- new_group$formation_no_sr <- nrow(organized_rocks)+1;
		organized_rocks <- rbind(organized_rocks,new_group);
		}
	}

#organized_rocks <- data.frame(cbind(rock_unit_senior,formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal),stringsAsFactors=hell_no);
#paleodb_rocks <- cbind(paleodb_rocks,organized_rocks);
#wagner_rocks <- cbind(wagner_rocks,as.character(rock_unit_senior),formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal);
return(organized_rocks);
}

match_pbdb_collections_to_paleodb_rock_thesaurus <- function(collection_no,paleodb_collections,paleodb_rock_thesaurus)	{
ncolls <- nrow(paleodb_collections);
colls_w_formations <- (1:ncolls)[paleodb_collections$formation!=""];
colls_w_members <- (1:ncolls)[paleodb_collections$member!=""];
colls_w_formations_and_members <- colls_w_formations[colls_w_formations %in% colls_w_members];
colls_w_members_only <- colls_w_members[!colls_w_members %in% colls_w_formations];
colls_w_groups <- (1:ncolls)[paleodb_collections$stratgroup!=""];
colls_w_groups_only <- colls_w_groups[!colls_w_groups %in% c(colls_w_formations,colls_w_members_only)];
colls_w_rocks <- sort(unique(c(colls_w_formations,colls_w_members,colls_w_groups)));

paleodb_clean_formation_basic <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_basic[colls_w_formations] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_member_basic <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_basic[colls_w_members] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_rock_unit_basic <- paleodb_clean_formation_basic;
paleodb_clean_rock_unit_basic[colls_w_formations_and_members] <- paste(paleodb_clean_formation_basic[colls_w_formations_and_members]," (",paleodb_clean_member_basic[colls_w_formations_and_members],")",sep="");
paleodb_clean_rock_unit_basic[colls_w_members_only] <- paleodb_clean_member_basic[colls_w_members_only];

paleodb_clean_group_basic <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_basic[colls_w_groups] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_formation_no_rock <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_no_rock[colls_w_formations] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_member_no_rock <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_no_rock[colls_w_members] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_rock_unit_no_rock <- paleodb_clean_formation_no_rock;
qqq <- colls_w_formations_and_members[paleodb_clean_member_no_rock[colls_w_formations_and_members]!=""];
paleodb_clean_rock_unit_no_rock[qqq] <- paste(paleodb_clean_formation_no_rock[qqq]," (",paleodb_clean_member_no_rock[qqq],")",sep="");
paleodb_clean_rock_unit_no_rock[colls_w_members_only] <- paleodb_clean_member_no_rock[colls_w_members_only];

paleodb_clean_group_no_rock <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_no_rock[colls_w_groups] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

cleaned_rocks <- data.frame(paleodb_clean_formation_basic=as.character(paleodb_clean_formation_basic),paleodb_clean_member_basic=as.character(paleodb_clean_member_basic),paleodb_clean_group_basic=as.character(paleodb_clean_group_basic),paleodb_clean_formation_no_rock=as.character(paleodb_clean_formation_no_rock),paleodb_clean_member_no_rock=as.character(paleodb_clean_member_no_rock),paleodb_clean_group_no_rock=as.character(paleodb_clean_group_no_rock),paleodb_clean_rock_unit_basic=as.character(paleodb_clean_rock_unit_basic),paleodb_clean_rock_unit_no_rock=as.character(paleodb_clean_rock_unit_no_rock),stringsAsFactors=hell_no);
cleaned_rocks <- expello_na_from_matrix(data=cleaned_rocks,replacement="");

unique_rock_units <- sort(unique(paleodb_clean_rock_unit_no_rock[paleodb_clean_rock_unit_no_rock!=""]));
rock_nos <- vector(length=ncolls);
nrocks <- nrow(paleodb_rock_thesaurus);
no_finds <- c();
for (uru in 1:length(unique_rock_units))	{
	rn <- (1:nrocks)[paleodb_rock_thesaurus$rock_unit_clean_no_rock==unique_rock_units[uru]];
#	rn <- unique(which(paleodb_rock_thesaurus$rock_unit_clean_no_rock==paleodb_clean_rock_unit_no_rock[uru],arr.ind = T)[,1]);
	if (length(rn)==1)	{
		rock_nos[(1:ncolls)[paleodb_clean_rock_unit_no_rock==unique_rock_units[uru]]] <- rn;
		} else	{
		no_finds <- c(no_finds,uru);
		}
	}

unided_colls <- colls_w_rocks[rock_nos[colls_w_rocks]==0];
unided_colls_w_formations <- unided_colls[paleodb_clean_rock_unit_basic[unided_colls]!=""];
#for (ucwf in 1:length(unided_colls_w_formations))	{
#for (uc in 1:length(unided_colls))	{
uc <- 0;
while (uc < length(unided_colls))	{
	uc <- uc+1;
	coll_no <- unided_colls[uc];
	if (cleaned_rocks$paleodb_clean_member_basic[coll_no]!="")	{
		rn <- match(cleaned_rocks$paleodb_clean_member_basic[coll_no],paleodb_rock_thesaurus$formation_clean_basic);
		if (is.na(rn))
			rn <- match(cleaned_rocks$paleodb_clean_member_no_rock[coll_no],paleodb_rock_thesaurus$formation_clean_no_rock);
#		rns <- unique(which(paleodb_rock_thesaurus==cleaned_rocks$paleodb_clean_member_basic[coll_no],arr.ind = T)[,1]);
		if (!is.na(rn))	{
			# assign collections with same member to this rock number
			# clear them out of the que, too.
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_member_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_member_basic[coll_no]];
			rock_nos[assigned] <- rn;
			unided_colls <- unided_colls[!unided_colls %in% assigned];
			uc <- uc-1;
			}
		}
	if (rock_nos[coll_no]==0 && cleaned_rocks$paleodb_clean_group_basic[coll_no]!="")	{
		# case where group is id'ed
		rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
		if (!is.na(rn) && length(rn)>1)	{
			# if there is 1+ matches
			if (sum(paleodb_rock_thesaurus$formation[rn]=="")==0)	{
				new_group <- paleodb_rock_thesaurus[rn[1],];
				new_group$rock_no <- new_group$rock_no_sr <- new_group$formation_no <- new_group$formation_no_sr <- nrow(paleodb_rock_thesaurus)+1;
				new_group$full_name <- cleaned_rocks$paleodb_clean_group_basic[coll_no];
				new_group$formation <- new_group$member <- new_group$formation_clean_basic <- new_group$formation_clean_no_rock <- new_group$formation_clean_no_rock_formal <- new_group$member_clean_basic <- new_group$member_clean_no_rock <- new_group$member_clean_no_rock_formal <- "";
				new_group$ma_lb <- max(abs(paleodb_rock_thesaurus$ma_lb[rn]));
				new_group$ma_ub <- min(abs(paleodb_rock_thesaurus$ma_ub[rn]));
				new_group$interval_lb <- paleodb_rock_thesaurus$interval_lb[rn[[match(new_group$ma_lb,paleodb_rock_thesaurus$ma_lb[rn])]]];
				new_group$interval_ub <- paleodb_rock_thesaurus$interval_ub[rn[[match(new_group$ma_ub,paleodb_rock_thesaurus$ma_ub[rn])]]];
				paleodb_rock_thesaurus <- rbind(paleodb_rock_thesaurus,new_group);
				rn <- nrow(paleodb_rock_thesaurus);
				} else	{
				rn <- rn[paleodb_rock_thesaurus$formation[rn]==""];
				}
			}
		if (is.na(rn) || length(rn)==0)
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_rock_unit_basic[coll_no]];
		if (is.na(rn))
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_no_rock %in% cleaned_rocks$paleodb_clean_rock_unit_no_rock[coll_no]];
		if (!is.na(rn) && length(rn)==1)	{
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_group_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
			assigned <- assigned[paleodb_collections$formation[assigned]==""];
			if (length(assigned)>0)	{
				rock_nos[assigned] <- rn;
				unided_colls <- unided_colls[!unided_colls %in% assigned];
				uc <- uc-1;
				}
			}
		}
#	print(c(uc,length(unided_colls)));
	}
return(rock_nos)
}

# routine to take rock unit thesaurus and provide numbers for unique rock units
# modified 2020-02-15
match_paleodb_collections_to_paleodb_rock_thesaurus <- function(paleodb_collections,paleodb_rock_thesaurus)	{
ncolls <- nrow(paleodb_collections);
colls_w_formations <- (1:ncolls)[paleodb_collections$formation!=""];
colls_w_members <- (1:ncolls)[paleodb_collections$member!=""];
colls_w_formations_and_members <- colls_w_formations[colls_w_formations %in% colls_w_members];
colls_w_members_only <- colls_w_members[!colls_w_members %in% colls_w_formations];
colls_w_groups <- (1:ncolls)[paleodb_collections$stratgroup!=""];
colls_w_groups_only <- colls_w_groups[!colls_w_groups %in% c(colls_w_formations,colls_w_members_only)];
colls_w_rocks <- sort(unique(c(colls_w_formations,colls_w_members,colls_w_groups)));

paleodb_clean_formation_basic <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_basic[colls_w_formations] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_member_basic <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_basic[colls_w_members] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_rock_unit_basic <- paleodb_clean_formation_basic;
paleodb_clean_rock_unit_basic[colls_w_formations_and_members] <- paste(paleodb_clean_formation_basic[colls_w_formations_and_members]," (",paleodb_clean_member_basic[colls_w_formations_and_members],")",sep="");
paleodb_clean_rock_unit_basic[colls_w_members_only] <- paleodb_clean_member_basic[colls_w_members_only];

paleodb_clean_group_basic <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_basic[colls_w_groups] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T);

paleodb_clean_formation_no_rock <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_no_rock[colls_w_formations] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_member_no_rock <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_no_rock[colls_w_members] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_rock_unit_no_rock <- paleodb_clean_formation_no_rock;
qqq <- colls_w_formations_and_members[paleodb_clean_member_no_rock[colls_w_formations_and_members]!=""];
paleodb_clean_rock_unit_no_rock[qqq] <- paste(paleodb_clean_formation_no_rock[qqq]," (",paleodb_clean_member_no_rock[qqq],")",sep="");
paleodb_clean_rock_unit_no_rock[colls_w_members_only] <- paleodb_clean_member_no_rock[colls_w_members_only];

paleodb_clean_group_no_rock <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_no_rock[colls_w_groups] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

cleaned_rocks <- data.frame(paleodb_clean_formation_basic=as.character(paleodb_clean_formation_basic),paleodb_clean_member_basic=as.character(paleodb_clean_member_basic),paleodb_clean_group_basic=as.character(paleodb_clean_group_basic),paleodb_clean_formation_no_rock=as.character(paleodb_clean_formation_no_rock),paleodb_clean_member_no_rock=as.character(paleodb_clean_member_no_rock),paleodb_clean_group_no_rock=as.character(paleodb_clean_group_no_rock),paleodb_clean_rock_unit_basic=as.character(paleodb_clean_rock_unit_basic),paleodb_clean_rock_unit_no_rock=as.character(paleodb_clean_rock_unit_no_rock),stringsAsFactors=hell_no);
cleaned_rocks <- expello_na_from_matrix(data=cleaned_rocks,replacement="");

unique_rock_units <- sort(unique(paleodb_clean_rock_unit_no_rock[paleodb_clean_rock_unit_no_rock!=""]));
rock_nos <- vector(length=ncolls);
nrocks <- nrow(paleodb_rock_thesaurus);
no_finds <- c();
for (uru in 1:length(unique_rock_units))	{
	rn <- (1:nrocks)[paleodb_rock_thesaurus$rock_unit_clean_no_rock==unique_rock_units[uru]];
#	rn <- unique(which(paleodb_rock_thesaurus$rock_unit_clean_no_rock==paleodb_clean_rock_unit_no_rock[uru],arr.ind = T)[,1]);
	if (length(rn)==1)	{
		rock_nos[(1:ncolls)[paleodb_clean_rock_unit_no_rock==unique_rock_units[uru]]] <- rn;
		} else	{
		no_finds <- c(no_finds,uru);
		}
	}

unided_colls <- colls_w_rocks[rock_nos[colls_w_rocks]==0];
unided_colls_w_formations <- unided_colls[paleodb_clean_rock_unit_basic[unided_colls]!=""];
#for (ucwf in 1:length(unided_colls_w_formations))	{
#for (uc in 1:length(unided_colls))	{
uc <- 0;
while (uc < length(unided_colls))	{
	uc <- uc+1;
	coll_no <- unided_colls[uc];
	if (cleaned_rocks$paleodb_clean_member_basic[coll_no]!="")	{
		rn <- match(cleaned_rocks$paleodb_clean_member_basic[coll_no],paleodb_rock_thesaurus$formation_clean_basic);
		if (is.na(rn))
			rn <- match(cleaned_rocks$paleodb_clean_member_no_rock[coll_no],paleodb_rock_thesaurus$formation_clean_no_rock);
#		rns <- unique(which(paleodb_rock_thesaurus==cleaned_rocks$paleodb_clean_member_basic[coll_no],arr.ind = T)[,1]);
		if (!is.na(rn))	{
			# assign collections with same member to this rock number
			# clear them out of the que, too.
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_member_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_member_basic[coll_no]];
			rock_nos[assigned] <- rn;
			unided_colls <- unided_colls[!unided_colls %in% assigned];
			uc <- uc-1;
			}
		}
	if (rock_nos[coll_no]==0 && cleaned_rocks$paleodb_clean_group_basic[coll_no]!="")	{
		# case where group is id'ed
		rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
		if (!is.na(rn) && length(rn)>1)	{
			# if there is 1+ matches
			if (sum(paleodb_rock_thesaurus$formation[rn]=="")==0)	{
				new_group <- paleodb_rock_thesaurus[rn[1],];
				new_group$rock_no <- new_group$rock_no_sr <- new_group$formation_no <- new_group$formation_no_sr <- nrow(paleodb_rock_thesaurus)+1;
				new_group$full_name <- cleaned_rocks$paleodb_clean_group_basic[coll_no];
				new_group$formation <- new_group$member <- new_group$formation_clean_basic <- new_group$formation_clean_no_rock <- new_group$formation_clean_no_rock_formal <- new_group$member_clean_basic <- new_group$member_clean_no_rock <- new_group$member_clean_no_rock_formal <- "";
				new_group$ma_lb <- max(abs(paleodb_rock_thesaurus$ma_lb[rn]));
				new_group$ma_ub <- min(abs(paleodb_rock_thesaurus$ma_ub[rn]));
				new_group$interval_lb <- paleodb_rock_thesaurus$interval_lb[rn[[match(new_group$ma_lb,paleodb_rock_thesaurus$ma_lb[rn])]]];
				new_group$interval_ub <- paleodb_rock_thesaurus$interval_ub[rn[[match(new_group$ma_ub,paleodb_rock_thesaurus$ma_ub[rn])]]];
				paleodb_rock_thesaurus <- rbind(paleodb_rock_thesaurus,new_group);
				rn <- nrow(paleodb_rock_thesaurus);
				} else	{
				rn <- rn[paleodb_rock_thesaurus$formation[rn]==""];
				}
			}
		if (is.na(rn) || length(rn)==0)
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_rock_unit_basic[coll_no]];
		if (is.na(rn))
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_no_rock %in% cleaned_rocks$paleodb_clean_rock_unit_no_rock[coll_no]];
		if (!is.na(rn) && length(rn)==1)	{
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_group_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
			assigned <- assigned[paleodb_collections$formation[assigned]==""];
			if (length(assigned)>0)	{
				rock_nos[assigned] <- rn;
				unided_colls <- unided_colls[!unided_colls %in% assigned];
				uc <- uc-1;
				}
			}
		}
#	print(c(uc,length(unided_colls)));
	}
return(rock_nos)
}

# routine to put unique numbers on unnamed rock units unique to their geoplate & time
# modified 2020-03-04
# modified 2021-03-15: code was assuming formation_no_sr was present: maybe not!
name_unnamed_rock_units <- function(paleodb_collections,finest_chronostrat)	{
# paleodb_collections <- refined_collections;
intervals <- finest_chronostrat$interval;
ncolls <- nrow(paleodb_collections);
if (is.null(paleodb_collections$bin_lb))	{
	paleodb_collections$bin_lb <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb,finest_chronostrat$interval)];
	paleodb_collections$bin_ub <- finest_chronostrat$bin_last[match(paleodb_collections$interval_ub,finest_chronostrat$interval)];
	}
paleodb_collections$rock_no_sr[is.na(paleodb_collections$rock_no_sr)] <- 0;
paleodb_collections$formation_no[is.na(paleodb_collections$formation_no)] <- 0;
for (i in min(paleodb_collections$bin_lb):max(paleodb_collections$bin_ub))	{
	interval_collections <- subset(paleodb_collections,paleodb_collections$interval_lb==intervals[i]);
	interval_collections <- subset(interval_collections,interval_collections$interval_ub==intervals[i]);
	nloc <- nrow(interval_collections);
	if (nrow(interval_collections)>0)	{
		interval_collections_named_rocks <- subset(interval_collections,interval_collections$rock_no_sr>0);
		interval_collections_unnamed_rocks <- subset(interval_collections,interval_collections$rock_no_sr==0);
		# assumed that unnamed rocks on geoplates with no other rocks are unique rock units
		rockless_plates <- unique(as.numeric(interval_collections_unnamed_rocks$geoplate));
		rockless_plates <- rockless_plates[!is.na(rockless_plates)];
		rock_ided_plates <- unique(as.numeric(interval_collections_named_rocks$geoplate));
		rock_ided_plates <- rock_ided_plates[!is.na(rock_ided_plates)];
		plate_needs_rock <- rockless_plates[!rockless_plates %in% rock_ided_plates];
		if (length(plate_needs_rock)>0)	{
			# add inferred rock units to collections
			implied_rocks <- paste("Geoplate",plate_needs_rock,intervals[i]);
			coll_w_nameable_rocks <- (1:nloc)[interval_collections$geoplate %in% rockless_plates[!rockless_plates %in% rock_ided_plates]];
			interval_collections$formation[coll_w_nameable_rocks] <- implied_rocks[match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates])];
			mx_rock_no <- max(paleodb_collections$rock_no_sr);
			interval_collections$rock_no_sr[coll_w_nameable_rocks] <- mx_rock_no+match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates]);
			update_these <- match(interval_collections$collection_no[coll_w_nameable_rocks],paleodb_collections$collection_no);
			paleodb_collections$formation[update_these] <- interval_collections$formation[coll_w_nameable_rocks];
			paleodb_collections$formation_sr[update_these] <- interval_collections$formation_sr[coll_w_nameable_rocks];
			paleodb_collections$rock_no_sr[update_these] <- interval_collections$rock_no_sr[coll_w_nameable_rocks];
			}
		} else	{
		mx_rock_no <- 0;
		}
	}

if (sum(paleodb_collections$rock_no_sr)==0)	{
#if (sum(paleodb_collections$rock_no_sr==0)>0)	{
	geoplates <- sort(unique(paleodb_collections$geoplate[paleodb_collections$geoplate>0]));
	geoplates <- geoplates[!is.na(geoplates)];
	if (sum(geoplates==0)>0)	{
		bleh <- subset(paleodb_collections,paleodb_collections$geoplate==0);
		bleh$cc[bleh$cc==""] <- "KZ";
		unique_ccs <- unique(bleh$cc);
		for (uc in 1:length(unique_ccs))	{
			this_cc <- subset(paleodb_collections,paleodb_collections$cc==unique_ccs[uc]);
			this_cc_prob <- subset(this_cc,this_cc$geoplate==0);
			this_cc <- subset(this_cc,this_cc$geoplate>0);
			other_plates <- hist(this_cc$geoplate,breaks=0:1000,plot=F)$count;
			obs_plates <- (1:1000)[other_plates!=0];
			other_plates <- other_plates[other_plates>0];
			paleodb_collections$geoplate[match(this_cc_prob$collection_no,paleodb_collections$collection_no)] <- obs_plates[match(max(other_plates),other_plates)];
			}
		geoplates <- as.numeric(unique(paleodb_collections$geoplate));
		}
	gp <- 0;
	while (gp < length(geoplates))	{
		gp <- gp+1;
		this_plate <- (1:nrow(paleodb_collections))[as.numeric(paleodb_collections$geoplate)==geoplates[gp]];
		paleodb_collections$rock_no_sr[this_plate] <- paleodb_collections$rock_no[this_plate] <- gp+mx_rock_no;
		if (!is.null(paleodb_collections$formation_no_sr))	{
			paleodb_collections$formation_no_sr[this_plate] <- paleodb_collections$formation_no[this_plate] <- gp+mx_rock_no;
			} else	{
			paleodb_collections$formation_no[this_plate] <- gp+mx_rock_no;
			}
		}
	}

return(paleodb_collections);
}

# modified 2020-02-15
name_unnamed_rock_units_old <- function(paleodb_collections,finest_chronostrat)	{
# paleodb_collections <- refined_collections;
intervals <- finest_chronostrat$interval;
ncolls <- nrow(paleodb_collections);
if (is.null(paleodb_collections$bin_lb))	{
	paleodb_collections$bin_lb <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb,finest_chronostrat$interval)];
	paleodb_collections$bin_ub <- finest_chronostrat$bin_last[match(paleodb_collections$interval_ub,finest_chronostrat$interval)];
	}
for (i in min(paleodb_collections$bin_lb):max(paleodb_collections$bin_ub))	{
	interval_collections <- subset(paleodb_collections,paleodb_collections$interval_lb==intervals[i]);
	interval_collections <- subset(interval_collections,interval_collections$interval_ub==intervals[i]);
	nloc <- nrow(interval_collections);
	if (nrow(interval_collections)>0)	{
		interval_collections_named_rocks <- subset(interval_collections,interval_collections$rock_no_sr>0);
		interval_collections_unnamed_rocks <- subset(interval_collections,interval_collections$rock_no_sr==0);
		# assumed that unnamed rocks on geoplates with no other rocks are unique rock units
		rockless_plates <- unique(as.numeric(interval_collections_unnamed_rocks$geoplate));
		rockless_plates <- rockless_plates[!is.na(rockless_plates)];
		rock_ided_plates <- unique(as.numeric(interval_collections_named_rocks$geoplate));
		rock_ided_plates <- rock_ided_plates[!is.na(rock_ided_plates)];
		plate_needs_rock <- rockless_plates[!rockless_plates %in% rock_ided_plates];
		if (length(plate_needs_rock)>0)	{
			# add inferred rock units to collections
			implied_rocks <- paste("Geoplate",plate_needs_rock,intervals[i]);
			coll_w_nameable_rocks <- (1:nloc)[interval_collections$geoplate %in% rockless_plates[!rockless_plates %in% rock_ided_plates]];
			interval_collections$formation[coll_w_nameable_rocks] <- implied_rocks[match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates])];
			mx_rock_no <- max(paleodb_collections$rock_no_sr);
			interval_collections$rock_no_sr[coll_w_nameable_rocks] <- mx_rock_no+match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates]);
			update_these <- match(interval_collections$collection_no[coll_w_nameable_rocks],paleodb_collections$collection_no);
			paleodb_collections$formation[update_these] <- interval_collections$formation[coll_w_nameable_rocks];
			paleodb_collections$formation_sr[update_these] <- interval_collections$formation_sr[coll_w_nameable_rocks];
			paleodb_collections$rock_no_sr[update_these] <- interval_collections$rock_no_sr[coll_w_nameable_rocks];
			}
		} else	{
		mx_rock_no <- 0;
		}
	}

if (sum(paleodb_collections$rock_no_sr)==0)	{
	geoplates <- sort(unique(paleodb_collections$geoplate[paleodb_collections$geoplate>0]));
	geoplates <- geoplates[!is.na(geoplates)];
	if (sum(geoplates==0)>0)	{
		bleh <- subset(paleodb_collections,paleodb_collections$geoplate==0);
		bleh$cc[bleh$cc==""] <- "KZ";
		unique_ccs <- unique(bleh$cc);
		for (uc in 1:length(unique_ccs))	{
			this_cc <- subset(paleodb_collections,paleodb_collections$cc==unique_ccs[uc]);
			this_cc_prob <- subset(this_cc,this_cc$geoplate==0);
			this_cc <- subset(this_cc,this_cc$geoplate>0);
			other_plates <- hist(this_cc$geoplate,breaks=0:1000,plot=F)$count;
			obs_plates <- (1:1000)[other_plates!=0];
			other_plates <- other_plates[other_plates>0];
			paleodb_collections$geoplate[match(this_cc_prob$collection_no,paleodb_collections$collection_no)] <- obs_plates[match(max(other_plates),other_plates)];
			}
		geoplates <- as.numeric(unique(paleodb_collections$geoplate));
		}
	gp <- 0;
	while (gp < length(geoplates))	{
		gp <- gp+1;
		this_plate <- (1:nrow(paleodb_collections))[as.numeric(paleodb_collections$geoplate)==geoplates[gp]];
		paleodb_collections$rock_no_sr[this_plate] <- paleodb_collections$rock_no[this_plate] <- gp+mx_rock_no;
#		paleodb_collections$formation_no_sr[this_plate] <- paleodb_collections$formation_no[this_plate] <- gp+mx_rock_no;
		paleodb_collections$formation_no[this_plate] <- paleodb_collections$formation_no[this_plate] <- gp+mx_rock_no;
		}
	}

return(paleodb_collections);
}

organize_zone_database_for_paleodb_matches <- function(external_zone_database="")	{
if (external_zone_database!="")	{
	file_type <- revelare_file_type(filename = external_zone_database);
	if (file_type=="csv")	{
		zone_data <- read.csv(file=external_zone_database,stringsAsFactors=FALSE,header=TRUE);
		}	else	{
		zone_data <- read.table(file=external_zone_database,stringsAsFactors=FALSE,header=TRUE,sep="\t");
		}
	} else	{
	zone_data <- rock_unit_data$zone_database;
	}
zd <- nrow(zone_data);

#zone_data$zone <- zone <- as.character(zone_data$zone);
#zone_data$zone <- sapply(zone,mundus_zone);
#zone_data$zone_sr <- zone <- as.character(zone_data$zone_sr);
#zone_data$zone_sr <- sapply(zone,mundus_zone);

zone <- zone_data$zone;
print("Separating out species epithets...")
zone_species <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_species_names_only);
#for (zz in 1:length(zone))
#	transmogrify_full_zone_names_to_species_names_only(zone[zz]);
#
print("Deal with zones not named after taxa...")
Nontaxon_Zone <- raster::t(pbapply::pbsapply(zone,aparecium_nontaxon_zone));
rownames(Nontaxon_Zone) <- NULL;
Nontaxon_Zone <- data.frame(non_taxon_zone=as.character(Nontaxon_Zone[,1]),non_taxon_zone_label=as.character(Nontaxon_Zone[,2]),stringsAsFactors = F);

zone <- zone_data$zone_sr;
print("Separating out species epithets for senior names...")
zone_species_sr <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_species_names_only);

print("Deal with senior names of zones not named after taxa...")
Nontaxon_zone_sr <- raster::t(pbapply::pbsapply(zone,aparecium_nontaxon_zone));
rownames(Nontaxon_zone_sr) <- NULL;
colnames(Nontaxon_zone_sr) <- c("non_taxon_zone_sr","non_taxon_zone_label_sr");
Nontaxon_zone_sr <- data.frame(non_taxon_zone_sr=as.character(Nontaxon_zone_sr[,1]),non_taxon_zone_label_sr=as.character(Nontaxon_zone_sr[,2]),stringsAsFactors = F);

zone <- zone_data$zone;
print("Find different genus-species combinations for the same zone...")
new_zone_combos <- base::t(pbapply::pbsapply(zone,accersi_genus_subgenus_combinations_for_zones));
#new_zone_combos <- c();
#for (z in 1:length(zone_data$zone))	{
#	new_zone_combos <- rbind(new_zone_combos,accersi_genus_subgenus_combinations_for_zones(zone=zone_data$zone[z]));
#	}
#zone_species <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
#print("Get zone epithets again...")
print("Saving zones as just species names (when possible)...");
zone_epithet <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_species_names_only);
print("Saving zones as just genus names (when possible)...");
genus_name <- zone_genera <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_genus_names_only);
print("Treat subgenera as genera as an option...");
zone_genera_subgenera <- base::t(pbapply::pbsapply(genus_name,divido_subgenus_names_from_genus_names));
abbreviated_zone <- c();
for (zgs in 1:nrow(zone_genera_subgenera))	{
	if (zone_genera_subgenera[zgs,2]!="" && zone_epithet[zgs]!="")	{
		abbreviated_zone <- c(abbreviated_zone,paste(strsplit(zone_genera_subgenera[zgs,1],split = "")[[1]][1],". (",strsplit(zone_genera_subgenera[zgs,2],split = "")[[1]][1],"). ",zone_epithet[zgs],sep=""));
		} else	if (zone_epithet[zgs]!="")	{
		split_genus <- strsplit(zone_genera_subgenera[zgs,1],split="-")[[1]];
		if (length(split_genus)>1)	{
			split_species <- strsplit(zone_epithet[zgs],split="-")[[1]];
			if (length(split_species)==length(split_genus))	{
				multizone <- c();
				for (ss in 1:length(split_genus))	{
					if (split_species[ss]!="")	{
						multizone <- c(multizone,paste(strsplit(split_genus[ss],split = "")[[1]][1],". ",split_species[ss],sep=""));
						} else	{
						multizone <- c(multizone,paste(strsplit(split_genus[ss],split = "")[[1]][1],"."));
						}
					}
				abbreviated_zone <- c(abbreviated_zone,paste(multizone,collapse=" - "));
				} else	{
				abbreviated_zone <- c(abbreviated_zone,"•");
				}
			} else	{
			abbreviated_zone <- c(abbreviated_zone,paste(strsplit(zone_genera_subgenera[zgs,1],split = "")[[1]][1],". ",zone_epithet[zgs],sep=""));
			}
		} else	{
		abbreviated_zone <- c(abbreviated_zone,"•");
		}
	}
zone <- zone_data$zone_sr;
print("Finding senior synonyms of species used for zones");
zone_epithet_sr <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_species_names_only);
Nontaxon_Zone <- data.frame(Nontaxon_Zone);
new_zone_info <- data.frame(zone_species=as.character(zone_species),zone_species_sr=as.character(zone_species_sr),non_taxon_zone=as.character(Nontaxon_Zone$non_taxon_zone),non_taxon_zone_label=as.character(Nontaxon_Zone$non_taxon_zone_label),non_taxon_zone_sr=as.character(Nontaxon_zone_sr$non_taxon_zone_sr),non_taxon_zone_label_sr=as.character(Nontaxon_zone_sr$non_taxon_zone_label_sr),
							genus_species_combo=as.character(new_zone_combos[,1]),subgenus_species_combo=as.character(new_zone_combos[,2]),zone_epithet=as.character(zone_epithet),zone_epithet_sr=as.character(zone_epithet_sr),abbreviated_zone=as.character(abbreviated_zone),stringsAsFactors = F);
#colnames(new_zone_combos) <- c("Genus_species_combo","Subgenus_species_combo","zone_epithet","zone_epithet_sr","abbreviated_zone");
old_fields <- colnames(new_zone_info)[colnames(new_zone_info) %in% colnames(zone_data)];
new_fields <- colnames(new_zone_info)[!colnames(new_zone_info) %in% colnames(zone_data)];
if (length(old_fields)>0)	zone_data[,match(old_fields,colnames(zone_data))] <- new_zone_info[,match(old_fields,colnames(new_zone_info))];

if (length(new_fields)>0)	zone_data <- cbind(zone_data,new_zone_info[,match(new_fields,colnames(new_zone_info))]);

rownames(zone_data) <- NULL;
return(zone_data);
}

organize_paleodb_zones_for_database_matches <- function(paleodb_collections)	{
ncolls <- nrow(paleodb_collections);
coll_w_zones <- (1:ncolls)[paleodb_collections$zone!=""];

zone <- as.character(paleodb_collections$zone[coll_w_zones]);
paleodb_collections$zone[coll_w_zones] <- sapply(zone,mundus_zone);
zone_species <- array("",dim=ncolls);
zone_species[coll_w_zones] <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
paleodb_collections <- tibble::add_column(paleodb_collections,zone_species,.after="zone");
return(paleodb_collections);
}

organize_pbdb_rock_data_old <- function(paleodb_collections)	{
pbdb_rocks <- unique(data.frame(group=as.character(paleodb_collections$stratgroup),
								formation=as.character(paleodb_collections$formation),
								member=as.character(paleodb_collections$member),stringsAsFactors = F));
if (!is.null(paleodb_collections$formation_alt))	{
	dummy <- unique(data.frame(group=as.character(paleodb_collections$stratgroup_alt),
							   formation=as.character(paleodb_collections$formation_alt),
							   member=as.character(paleodb_collections$member_alt),stringsAsFactors = F));
	pbdb_rocks <- unique(rbind(pbdb_rocks,dummy));
	}

pbdb_rocks <- pbdb_rocks[order(pbdb_rocks$formation,pbdb_rocks$member,pbdb_rocks$group),];
#pbdb_rocks[1:4,]
nrocks <- nrow(pbdb_rocks);
gg <- (1:nrocks)[pbdb_rocks$group!=""]; ff <- (1:nrocks)[pbdb_rocks$formation!=""]; mm <- (1:nrocks)[pbdb_rocks$member!=""];
pbdb_rocks <- pbdb_rocks[sort(unique(c(gg,ff,mm))),];
nrocks <- nrow(pbdb_rocks);
gg <- (1:nrocks)[pbdb_rocks$group!=""]; ff <- (1:nrocks)[pbdb_rocks$formation!=""]; mm <- (1:nrocks)[pbdb_rocks$member!=""];
pbdb_rocks$full_name <- pbdb_rocks$formation;
pbdb_rocks$full_name[mm[mm %in% ff]] <- paste(pbdb_rocks$formation[mm[mm %in% ff]]," (",pbdb_rocks$member[mm[mm %in% ff]],")",sep="");
pbdb_rocks$full_name[mm[!mm %in% ff]] <- paste("(",pbdb_rocks$member[mm[!mm %in% ff]],")",sep="");
pbdb_rocks$full_name[gg[!gg %in% c(ff,mm)]] <- paste("[",pbdb_rocks$group[gg[!gg %in% c(ff,mm)]],"]",sep="");

nsites <- nrow(paleodb_collections);
paleodb_collections$full_name <- paleodb_collections$formation;
paleodb_collections$full_name[paleodb_collections$member!=""] <- paste(paleodb_collections$formation[paleodb_collections$member!=""]," (",paleodb_collections$member[paleodb_collections$member!=""],")",sep="");
members_only <- (1:nsites)[(paleodb_collections$member!="") & (paleodb_collections$formation=="")];
paleodb_collections$full_name[members_only] <- paste("(",paleodb_collections$member[members_only],")",sep="");
group_only <- (1:nsites)[paleodb_collections$stratgroup!="" & (paleodb_collections$member=="") & (paleodb_collections$formation=="")];
paleodb_collections$full_name[group_only] <- paste("[",paleodb_collections$stratgroup[group_only],"]",sep="");

unique_full_names <- unique(pbdb_rocks$full_name);
kill_me <- dups <- c();
for (ur in 1:length(unique_full_names))	{
	if (sum(pbdb_rocks$full_name %in% unique_full_names[ur])>1)	{
		this_case <- subset(pbdb_rocks,pbdb_rocks$full_name %in% unique_full_names[ur]);
		if (sum(this_case$group=="")>0)	{
			kill_me <- c(kill_me,(1:nrocks)[pbdb_rocks$full_name %in% unique_full_names[ur]][pbdb_rocks$group[pbdb_rocks$full_name %in% unique_full_names[ur]]==""]);
			} # end case where on rock isn't put into a rock uniit.
		if (sum(this_case$group!="")>1)	{
			this_case <- subset(this_case,this_case$group!="");
			these_rock_collections <- subset(paleodb_collections,paleodb_collections$full_name==unique_full_names[ur]);
			these_rock_collections <- subset(these_rock_collections,these_rock_collections$stratgroup!="");
			relv_groups <- unique(these_rock_collections$stratgroup);
			for (rg in 1:(length(relv_groups)-1))	{
				gr <- rg;
				lb_a <- these_rock_collections$max_ma[these_rock_collections$stratgroup==relv_groups[rg]];
				ub_a <- these_rock_collections$min_ma[these_rock_collections$stratgroup==relv_groups[rg]];
				while (gr < length(relv_groups))	{
					gr <- gr+1;
					lb_b <- these_rock_collections$max_ma[these_rock_collections$stratgroup==relv_groups[gr]];
					ub_b <- these_rock_collections$min_ma[these_rock_collections$stratgroup==relv_groups[gr]];
					if (do_two_ranges_overlap(lb_a,ub_a,lb_b,ub_b))	{
						# if they overlap, then eliminate the one last used longest ago
						if (max(these_rock_collections$ref_pubyr[these_rock_collections$stratgroup==relv_groups[rg]])>max(these_rock_collections$ref_pubyr[these_rock_collections$stratgroup==relv_groups[gr]]))	{
							kill_me <- unique(c(kill_me,match(rownames(this_case)[rg],rownames(pbdb_rocks))));
							} else	{
							kill_me <- unique(c(kill_me,match(rownames(this_case)[gr],rownames(pbdb_rocks))));
							}
						}
					}
				} # end looking at whether groups overlap & which is latest opinion
			# add something that lets disjunct ages & areas indicate homonyms;
			} # end case with 2+ groups
		} # end case of duplicate names;
	}
pbdb_rocks <- pbdb_rocks[!(1:nrocks) %in% kill_me,];
nrocks <- nrow(pbdb_rocks);
geoplate_list <- site_list <- list();
for (nr in 1:nrocks)	{
	site_list <- rlist::list.append(site_list,paleodb_collections$collection_no[paleodb_collections$full_name %in% pbdb_rocks$full_name[nr]]);
	geoplate_list <- rlist::list.append(geoplate_list,paleodb_collections$geoplate[paleodb_collections$full_name %in% pbdb_rocks$full_name[nr]]);
	}
names(site_list) <- names(geoplate_list) <- pbdb_rocks$full_name;
# "formation_clean_basic"	"formation_clean_no_rock"	"formation_clean_no_rock_formal"
pbdb_rocks$formation_clean_basic <- pbdb_rocks$formation_clean_no_rock <- pbdb_rocks$formation_clean_no_rock_formal <- pbdb_rocks$formation;
named_rock_unit <- pbdb_rocks$formation_clean_basic[pbdb_rocks$formation_clean_basic!=""];
pbdb_rocks$formation_clean_basic[pbdb_rocks$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$formation_clean_no_rock[pbdb_rocks$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$formation_clean_no_rock_formal[pbdb_rocks$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

# "member_clean_basic"      "member_clean_no_rock"      "member_clean_no_rock_formal"
pbdb_rocks$member_clean_basic <- pbdb_rocks$member_clean_no_rock <- pbdb_rocks$member_clean_no_rock_formal <- pbdb_rocks$member;
named_rock_unit <- pbdb_rocks$member_clean_basic[pbdb_rocks$member_clean_basic!=""];
pbdb_rocks$member_clean_basic[pbdb_rocks$member!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$member_clean_no_rock[pbdb_rocks$member!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$member_clean_no_rock_formal[pbdb_rocks$member!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

# "group_clean_basic"       "group_clean_no_rock"       "group_clean_no_rock_formal"
pbdb_rocks$group_clean_basic <- pbdb_rocks$group_clean_no_rock <- pbdb_rocks$group_clean_no_rock_formal <- pbdb_rocks$group;
named_rock_unit <- pbdb_rocks$group_clean_basic[pbdb_rocks$group_clean_basic!=""];
pbdb_rocks$group_clean_basic[pbdb_rocks$group!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$group_clean_no_rock[pbdb_rocks$group!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$group_clean_no_rock_formal[pbdb_rocks$group!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

# "rock_unit_clean_basic"   "rock_unit_clean_no_rock"   "rock_unit_clean_no_rock_formal"
pbdb_rocks$rock_unit_clean_basic <- pbdb_rocks$rock_unit_clean_no_rock <- pbdb_rocks$rock_unit_clean_no_rock_formal <- pbdb_rocks$full_name;
named_rock_unit <- pbdb_rocks$rock_unit_clean_basic[pbdb_rocks$rock_unit_clean_basic!=""];
pbdb_rocks$rock_unit_clean_basic[pbdb_rocks$full_name!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$rock_unit_clean_no_rock[pbdb_rocks$full_name!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$rock_unit_clean_no_rock_formal[pbdb_rocks$full_name!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

output <- list(pbdb_rocks,site_list,geoplate_list);
return(output);
}

#paleodb_zone_info <- match_collections_zones_to_zone_database(collections,zone_database);
#test <- subset(rocks,rocks$stratgroup=="Honda")
# construct a database of rock units from the PaleoDB
#	returns two tables: one for rock units, and another for zones to which rock units are assigned.
# modified 2020-02-15
# modified 2021-02-21
# construct a database of rock units from the PaleoDB
#	returns two tables: one for rock units, and another for zones to which rock units are assigned.
#### rewrite to use construct_stratigraphic_data_base_from_paleodb funciton below!!!
construct_stratigraphic_data_base_from_paleodb <- function(taxa,onset="Cambrian",end="Holocene",save_files=TRUE,directory="",output_type=".csv")	{
# taxa: taxon or list of taxa for which you want collection data
# onset: oldest rocks
# end: youngest rocks
# save_files: whether to save files
# directory: where to send files (if saved)
# output_type: ".csv" for ".csv"; ".txt" or ".tab" or ".xls" for tab-delimited
rocks <- accersi_rock_unit_data(taxa,onset,end,standardize_members = FALSE,directory,save_files=FALSE);
ncolls <- nrow(rocks);
rocks <- rocks[order(rocks$formation,rocks$member,rocks$zone,rocks$stratgroup),];
has_formation <- (1:ncolls)[rocks$formation!=""];
has_group <- (1:ncolls)[rocks$stratgroup!=""];
has_member <- (1:ncolls)[rocks$member!=""];
keepers <- sort(unique(c(has_formation,has_group,has_member)));
rocks <- rocks[keepers,];

ncolls <- nrow(rocks);
has_formation <- (1:ncolls)[rocks$formation!=""];
has_group <- (1:ncolls)[rocks$stratgroup!=""];
has_member <- (1:ncolls)[rocks$member!=""];
has_zone <- (1:ncolls)[rocks$zone!=""];
has_group_only <- has_group[!has_group %in% has_formation];
has_group_only <- has_group_only[!has_group_only %in% has_member];

formation_entered <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_entered[has_formation] <- sapply(named_rock_unit,mundify_rock_unit_names);
member_entered <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_entered[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names);
group_entered <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_entered[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names);

zone <- rocks$zone;
zones_cleaned <-  sapply(zone,mundus_zone);

formation_dehyph <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_dehyph[has_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);
member_dehyph <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_dehyph[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);
group_dehyph <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_dehyph[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);

formation_rockless <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless[has_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
member_rockless <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
group_rockless <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);

formation_rockless_formal <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless_formal[has_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
member_rockless_formal <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless_formal[has_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
group_rockless_formal <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless_formal[has_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

rock_base <- data.frame(formation=as.character(formation_entered),member=as.character(member_entered),group=as.character(group_entered),
			  early_interval=as.character(rocks$early_interval),late_interval=as.character(rocks$late_interval),zone = as.character(zones_cleaned),
			  max_ma=as.numeric(rocks$max_ma),min_ma=as.numeric(rocks$min_ma),
			  formation_dehyph=as.character(formation_dehyph),member_dehyph=as.character(member_dehyph),group_dehyph=as.character(group_dehyph),
			  formation_rockless=as.character(formation_rockless),member_rockless=as.character(member_rockless),group_rockless=as.character(group_rockless),
			  formation_rockless_formal=as.character(formation_rockless_formal),member_rockless_formal=as.character(member_rockless_formal),group_rockless_formal=as.character(group_rockless_formal),
			  stringsAsFactors = FALSE);
nrocks <- nrow(rock_base);
formation_names <- unique(formation_rockless[formation_rockless!=""]);

for (ff in 1:length(formation_names))	{
	ff_records <- subset(rock_base,rock_base$formation_rockless==formation_names[ff]);
	ff_r <- length(ff_records$group);
#	rb_r <- (1:nrocks)[rock_base$formation_rockless==formation_names[ff]];
	rb_r <- as.numeric(rownames(ff_records))
	if (sum(ff_records$group!="")>0 && sum(ff_records$group!="")<ff_r)	{
		ff_g <- (1:ff_r)[ff_records$group!=""][1];
		rock_base$group[rb_r] <- as.character(ff_records$group[ff_g]);
		rock_base$group_dehyph[rb_r] <- as.character(ff_records$group_dehyph[ff_g]);
		rock_base$group_rockless[rb_r] <- as.character(ff_records$group_rockless[ff_g]);
		rock_base$group_rockless_formal[rb_r] <- as.character(ff_records$group_rockless_formal[ff_g]);
		}
	formation_names_entered <- unique(ff_records$formation);
	if (length(formation_names_entered)==2)	{
#		print(ff);
		if (length(formation_names_entered[formation_names_entered!=formation_names[ff]])==2)	{
#			print(ff);
			# compare dehyphenated names to rockless names
			rock_base$formation[rb_r] <- sort(formation_names_entered[unique(formation_dehyph[rb_r])!=formation_names[ff]],decreasing = TRUE)[1];
#			sort(ff_records$formation[ff_records$formation_dehyph==ff_records$formation_rockless[1]],decreasing = TRUE)[1];
			} else	{
			rock_base$formation[rb_r] <- formation_names_entered[formation_names_entered!=formation_names[ff]];
			}
		} else if (length(formation_names_entered)>2)	{
		print(ff);
		}

	unique_members_entered <- sort(unique(ff_records$member[ff_records$member!=""]));
	unique_members_rockless <- sort(unique(ff_records$member_rockless[ff_records$member_rockless!=""]));

	if (length(unique_members_entered) > length(unique_members_rockless))	{
		for (mm in 1:length(unique_members_rockless))	{
			mm_records <- subset(ff_records,ff_records$member_rockless==unique_members_rockless[mm]);
			unique_members <- unique(mm_records$member);
#			if (length(unique(mm_records$member_dehyph))==1)	{
			if (length(unique_members)==2)	{
				rb_m <- rb_r[rock_base$member_rockless[rb_r]==unique_members_rockless[mm]];
#				rock_base$member[rb_m] <- unique(mm_records$member[mm_records$member!=unique_members_rockless[mm]]);
				if (length(unique_members[unique_members!=unique_members_rockless[mm]])==1)	{
					rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
					} else	{
					unique_members_dehyph <- unique(mm_records$member_dehyph);
					if (length(unique_members_dehyph[unique_members_dehyph!=unique_members_rockless[mm]])==1)	{
						rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
						}
					}
#				if (length(mm_records$member[mm_records$member!=unique_members_rockless[mm]])>1)
#					print(c(ff,mm));
#				} else if (length(unique_members)>2)	{
#				print(c(ff,mm));
				}
#			mm <- mm+1;
			}
		}
	}

rock_base <- unique(rock_base);

nrocks <- nrow(rock_base);
formation_member_names <- rock_base$formation_rockless_formal;
has_formation <- (1:nrocks)[rock_base$formation_rockless!=""];
formation_member_names[has_formation] <- rock_base$formation_rockless[has_formation];
w_members_rockless <- (1:nrocks)[rock_base$member_rockless!=""];
w_member_and_formation <- w_members_rockless[w_members_rockless %in% has_formation];
formation_member_names[w_member_and_formation] <- paste(as.character(rock_base$formation_rockless[w_member_and_formation])," (",as.character(rock_base$member_rockless[w_member_and_formation]),")",sep="");
rock_base <- cbind(rock_base,formation_member_names);

rocks_to_time_scale <- rocks_to_zones <- c();
#unique_rock_units <- unique(rock_base$formation[rock_base$formation!=""]);
formation_member_names <- formation_member_names[formation_member_names!=""];
unique_rock_units <- unique(sort(c(formation_names,formation_member_names)));

for (uru in 1:length(unique_rock_units))	{
	ru <- match(unique_rock_units[uru],rock_base$formation_member_names);
	if (is.na(ru))	{
		ff <- match(unique_rock_units[uru],rock_base$formation_rockless);
		frm <- rock_base$formation[ff];
		grp <- rock_base$group[ff];
		mmb <- "";
		} else	{
		grp <- rock_base$group[ru];
		frm <- rock_base$formation[ru];
		mmb <- rock_base$member[ru];
		}
	if (is.na(match(unique_rock_units[uru],formation_names)))	{
		ur_records <- subset(rock_base,rock_base$formation_member_names==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		} else	{
		ur_records <- subset(rock_base,rock_base$formation_rockless==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		}
	}

unique_groups <- unique(rock_base$group_rockless[has_group_only]);
gu <- 0;
while (gu < length(unique_groups))	{
	gu <- gu+1;
	grp <- unique_groups[gu];
	frm <- mmb <- "";
	gr_records <- subset(rock_base,rock_base$group_rockless==grp);
	if(gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]!="")	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		} else	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$early_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		}
	unique_zones <- sort(unique(gr_records$zone[gr_records$zone!=""]));
	nz <- length(unique_zones);
	if (nz>0)
		rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
	}

rock_database <- data.frame(formation=as.character(rocks_to_time_scale[,1]),
							member=as.character(rocks_to_time_scale[,2]),
							group=as.character(rocks_to_time_scale[,3]),
							formation_member=as.character(rocks_to_time_scale[,4]),
							interval_lb=as.character(rocks_to_time_scale[,7]),
							interval_ub=as.character(rocks_to_time_scale[,8]),
							ma_lb=as.numeric(rocks_to_time_scale[,5]),
							ma_ub=as.numeric(rocks_to_time_scale[,6]),
							stringsAsFactors = FALSE
							);
rocks_to_zones <- data.frame(formation=as.character(rocks_to_zones[,1]),member=as.character(rocks_to_zones[,2]),group=as.character(rocks_to_zones[,3]),zone=as.character(rocks_to_zones[,4]),stringsAsFactors = FALSE)
output <- list(rock_database,rocks_to_zones);
if (save_files==TRUE)	{
	if (onset!=end)	{
		file1 <- paste(directory,onset,"-",end,"_",taxa,"_Stratigraphic_Database",output_type,sep="");
		file2 <- paste(directory,onset,"-",end,"_",taxa,"_Rocks_to_Zones_Database",output_type,sep="");
		} else	{
		file1 <- paste(directory,end,"_",taxa,"_Stratigraphic_Database",output_type,sep="");
		file2 <- paste(directory,end,"_",taxa,"_Rocks_to_Zones_Database",output_type,sep="");
		}
	if (output_type==".csv" || output_type=="csv")	{
		write.csv(rock_database,file=file1,row.names=F,fileEncoding = "UTF-8");
		write.csv(rocks_to_time_scale,file=file2,row.names=F,fileEncoding = "UTF-8");
		} else	{
		write.table(rock_database,file=file1,sep="\t",row.names=F,col.name=T);
		write.table(rocks_to_time_scale,file=file2,sep="\t",row.names=F,col.name=T);
		}
	}
names(output) <- c("Rock_Database","Rocks_to_Zones");
return(output);
}

construct_stratigraphic_data_base_from_paleodb_collections <- function(paleodb_collections,zone_database="",time_scale)	{
# paleodb_collections: data.frame of collections data from PaleoDB
rocks <- accersi_rock_unit_data_from_paleodb_collections(paleodb_collections,standardize_members = T);
nrock_combos <- nrow(rocks);
rocks <- rocks[order(rocks$formation,rocks$member,rocks$zone,rocks$stratgroup,-rocks$max_ma),];
has_formation <- (1:nrock_combos)[rocks$formation!=""];
has_group <- (1:nrock_combos)[rocks$stratgroup!=""];
has_member <- (1:nrock_combos)[rocks$member!=""];
keepers <- sort(unique(c(has_formation,has_group,has_member)));	# keep only information with rock names
rocks <- rocks[keepers,];
rocks <- rocks[order(rocks$formation,rocks$member,rocks$stratgroup,-rocks$max_ma),];
# redo numbering without unnamed rocks
nrock_combos <- nrow(rocks);
has_formation <- (1:nrock_combos)[rocks$formation!=""];
has_group <- (1:nrock_combos)[rocks$stratgroup!=""];
has_member <- (1:nrock_combos)[rocks$member!=""];
has_zone <- (1:nrock_combos)[rocks$zone!=""];
has_group_only <- has_group[!has_group %in% has_formation];
has_group_only <- has_group_only[!has_group_only %in% has_member];

formation_entered <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
print("Cleaning Formation Names");
formation_entered[has_formation] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
member_entered <- rocks$member;
named_rock_unit <- rocks$member[has_member];
print("Cleaning Member Names");
member_entered[has_member] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);
group_entered <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
print("Cleaning Group Names");
group_entered[has_group] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names);

zones_cleaned <- zone <- rocks$zone;
print("Cleaning Zones");
zones_cleaned[zone!=""] <- pbapply::pbsapply(zone[zone!=""],mundus_zone);

formation_dehyph <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
print("More rock cleaning because cleanliness is apotheosis...");
formation_dehyph[has_formation] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);
member_dehyph <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_dehyph[has_member] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);
group_dehyph <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_dehyph[has_group] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE);

formation_rockless <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless[has_formation] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
member_rockless <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless[has_member] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
group_rockless <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless[has_group] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);

formation_rockless_formal <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless_formal[has_formation] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
member_rockless_formal <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless_formal[has_member] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
group_rockless_formal <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless_formal[has_group] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

# rock_base has all of the combinations of rock units, intervals and zones.
if (!is.null(rocks$interval_lb))	{
	rocks$early_interval <- rocks$interval_lb;
	rocks$late_interval <- rocks$interval_ub;
	}
if (!is.null(rocks$ma_lb))	{
	rocks$max_ma <- rocks$ma_lb;
	rocks$min_ma <- rocks$ma_ub;
	}
rock_base <- data.frame(formation=as.character(formation_entered),member=as.character(member_entered),group=as.character(group_entered),
			  early_interval=as.character(rocks$early_interval),late_interval=as.character(rocks$late_interval),zone = as.character(zones_cleaned),
			  max_ma=as.numeric(rocks$max_ma),min_ma=as.numeric(rocks$min_ma),
			  formation_dehyph=as.character(formation_dehyph),member_dehyph=as.character(member_dehyph),group_dehyph=as.character(group_dehyph),
			  formation_rockless=as.character(formation_rockless),member_rockless=as.character(member_rockless),group_rockless=as.character(group_rockless),
			  formation_rockless_formal=as.character(formation_rockless_formal),member_rockless_formal=as.character(member_rockless_formal),group_rockless_formal=as.character(group_rockless_formal),
			  stringsAsFactors = FALSE);
nrock_combos <- nrow(rock_base);
formation_names <- unique(formation_rockless[formation_rockless!=""]);

### START HERE!!!!
for (ff in 1:length(formation_names))	{
	ff_records <- subset(rock_base,rock_base$formation_rockless==formation_names[ff]);
	ff_r <- length(ff_records$group);
#	rb_r <- (1:nrocks)[rock_base$formation_rockless==formation_names[ff]];
	rb_r <- as.numeric(rownames(ff_records))
	if (sum(ff_records$group!="")>0 && sum(ff_records$group!="")<ff_r)	{
		ff_g <- (1:ff_r)[ff_records$group!=""][1];
		rock_base$group[rb_r] <- as.character(ff_records$group[ff_g]);
		rock_base$group_dehyph[rb_r] <- as.character(ff_records$group_dehyph[ff_g]);
		rock_base$group_rockless[rb_r] <- as.character(ff_records$group_rockless[ff_g]);
		rock_base$group_rockless_formal[rb_r] <- as.character(ff_records$group_rockless_formal[ff_g]);
		}
	formation_names_entered <- unique(ff_records$formation);
	if (length(formation_names_entered)==2)	{
#		print(ff);
		if (length(formation_names_entered[formation_names_entered!=formation_names[ff]])==2)	{
#			print(ff);
			# compare dehyphenated names to rockless names
			rock_base$formation[rb_r] <- sort(formation_names_entered[unique(formation_dehyph[rb_r])!=formation_names[ff]],decreasing = TRUE)[1];
#			sort(ff_records$formation[ff_records$formation_dehyph==ff_records$formation_rockless[1]],decreasing = TRUE)[1];
			} else	{
			rock_base$formation[rb_r] <- formation_names_entered[formation_names_entered!=formation_names[ff]];
			}
		} else if (length(formation_names_entered)>2)	{
		print(ff);
		}

	unique_members_entered <- sort(unique(ff_records$member[ff_records$member!=""]));
	unique_members_rockless <- sort(unique(ff_records$member_rockless[ff_records$member_rockless!=""]));

	if (length(unique_members_entered) > length(unique_members_rockless))	{
		for (mm in 1:length(unique_members_rockless))	{
			mm_records <- subset(ff_records,ff_records$member_rockless==unique_members_rockless[mm]);
			unique_members <- unique(mm_records$member);
#			if (length(unique(mm_records$member_dehyph))==1)	{
			if (length(unique_members)==2)	{
				rb_m <- rb_r[rock_base$member_rockless[rb_r]==unique_members_rockless[mm]];
#				rock_base$member[rb_m] <- unique(mm_records$member[mm_records$member!=unique_members_rockless[mm]]);
				if (length(unique_members[unique_members!=unique_members_rockless[mm]])==1)	{
					rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
					} else	{
					unique_members_dehyph <- unique(mm_records$member_dehyph);
					if (length(unique_members_dehyph[unique_members_dehyph!=unique_members_rockless[mm]])==1)	{
						rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
						}
					}
				}
			}
		}
	}

rock_base <- unique(rock_base);	# all remaining rock unit + stratigraphic assignment combinations

nrock_combos <- nrow(rock_base);
formation_member_names <- rock_base$formation_rockless_formal;
has_formation <- (1:nrock_combos)[rock_base$formation_rockless!=""];
formation_member_names[has_formation] <- rock_base$formation_rockless[has_formation];
w_members_rockless <- (1:nrock_combos)[rock_base$member_rockless!=""];
w_member_and_formation <- w_members_rockless[w_members_rockless %in% has_formation];
formation_member_names[w_member_and_formation] <- paste(as.character(rock_base$formation_rockless[w_member_and_formation])," (",as.character(rock_base$member_rockless[w_member_and_formation]),")",sep="");
rock_base <- cbind(rock_base,formation_member_names);

rocks_to_time_scale <- rocks_to_zones <- c();
#unique_rock_units <- unique(rock_base$formation[rock_base$formation!=""]);
formation_member_names <- formation_member_names[formation_member_names!=""];
unique_rock_units <- unique(sort(c(formation_names,formation_member_names)));

for (uru in 1:length(unique_rock_units))	{
	ru <- match(unique_rock_units[uru],rock_base$formation_member_names);
	if (is.na(ru))	{
		ff <- match(unique_rock_units[uru],rock_base$formation_rockless);
		frm <- rock_base$formation[ff];
		grp <- rock_base$group[ff];
		mmb <- "";
		} else	{
		grp <- rock_base$group[ru];
		frm <- rock_base$formation[ru];
		mmb <- rock_base$member[ru];
		}
	if (is.na(match(unique_rock_units[uru],formation_names)))	{
		# case where rock unit is not simply a formation
		ur_records <- subset(rock_base,rock_base$formation_member_names==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		} else	{
		# case where rock unit is simply a formation
		ur_records <- subset(rock_base,rock_base$formation_rockless==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			# add rock info to rocks_to_time_scale
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			# fill out table linking rocks to zones.
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		}
	}

unique_groups <- unique(rock_base$group_rockless[has_group_only]);
gu <- 0;
while (gu < length(unique_groups))	{
	gu <- gu+1;
	grp <- unique_groups[gu];
	frm <- mmb <- "";
	gr_records <- subset(rock_base,rock_base$group_rockless==grp);
	if(gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]!="")	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		} else	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$early_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		}
	unique_zones <- sort(unique(gr_records$zone[gr_records$zone!=""]));
	nz <- length(unique_zones);
	if (nz>0)
		rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
	}

rock_database <- data.frame(formation=as.character(rocks_to_time_scale[,1]),
							member=as.character(rocks_to_time_scale[,2]),
							group=as.character(rocks_to_time_scale[,3]),
							formation_member=as.character(rocks_to_time_scale[,4]),
							interval_lb=as.character(rocks_to_time_scale[,7]),
							interval_ub=as.character(rocks_to_time_scale[,8]),
							ma_lb=as.numeric(rocks_to_time_scale[,5]),
							ma_ub=as.numeric(rocks_to_time_scale[,6]),
							stringsAsFactors = FALSE);
rocks_to_zones <- data.frame(formation=as.character(rocks_to_zones[,1]),member=as.character(rocks_to_zones[,2]),group=as.character(rocks_to_zones[,3]),zone=as.character(rocks_to_zones[,4]),stringsAsFactors = FALSE);
nrocks <- nrow(rock_database);
r_t_z <- nrow(rocks_to_zones);

#if (!is.null(zone_database))	{
if (is.data.frame(zone_database))	{
	formations_to_zones <- (1:nrow(rocks_to_zones))[rocks_to_zones$formation!=""];
	members_to_zones <- (1:nrow(rocks_to_zones))[rocks_to_zones$member!=""];
	formations_and_members_to_zones <- members_to_zones[members_to_zones %in% formations_to_zones];
	members_only_to_zones <- members_to_zones[!members_to_zones %in% formations_to_zones];
	rock_units_w_zones <- rocks_to_zones$formation[formations_to_zones];
	rock_units_w_zones[formations_and_members_to_zones] <- paste(rocks_to_zones$formation[formations_and_members_to_zones]," (",rocks_to_zones$member[formations_and_members_to_zones],")",sep="");
	zoned_formations <- unique(rocks_to_zones$formation[rocks_to_zones$formation!=""]);
	zf <- 0;
	while (zf < length(zoned_formations))	{
		zf <- zf+1;
		this_form <- subset(rocks_to_zones,rocks_to_zones$formation==zoned_formations[zf]);
		# first collapse all of the zones into one string separated by '; '; then splitusing '; '; because multizone collections use '; '
		#		to separate zones, this will put all zones in a string 'a; b; c; d' which we then parts to c(a,b,c,d)
		form_zones <- unique(strsplit(paste(this_form$zone,collapse="; "),split="; ")[[1]]);
#		rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
		rn <- match(zoned_formations[zf],rock_database$formation_member);
		if (is.na(rn))	{
#			rocks_to_zones$formation[zf] <- mundify_rock_unit_names(named_rock_unit=rocks_to_zones$formation[zf],dehyphenate=T,delete_rock_type=T);
#			rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
			try_this <- mundify_rock_unit_names(named_rock_unit=zoned_formations[zf],dehyphenate=T,delete_rock_type=T);
			rn <- match(try_this,rock_database$formation_member);
			}
		if (is.na(rn))	{
#			rocks_to_zones$formation[zf] <- mundify_rock_unit_names(named_rock_unit=rocks_to_zones$formation[zf],dehyphenate=T,delete_rock_type=T,delete_informal = T);
#			rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
			try_this <- mundify_rock_unit_names(named_rock_unit=zoned_formations[zf],dehyphenate=T,delete_rock_type=T,delete_informal = T);
			rn <- match(try_this,rock_database$formation_member);
			}
		if (!is.na(rn))	{
			zone_nos <- c();
			for (zn in 1:length(form_zones))
				zone_nos <- c(zone_nos,unique(which(zone_database==form_zones[zn],arr.ind = T)[,1]));
			zone_nos <- sort(unique(zone_nos));
			if (length(zone_nos)>0)	{
				# if zones make rock younger than oldest possible age, then adjust accordingly
				this_rock_base <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				if (nrow(this_rock_base)==0)
					this_rock_base <- subset(rock_base,rock_base$formation_rockless==zoned_formations[zf]);
				if (nrow(this_rock_base)==0)
					this_rock_base <- subset(rock_base,rock_base$formation_rockless_formal==zoned_formations[zf]);
				zoneless_formation <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				zoneless_formation <- subset(zoneless_formation,zoneless_formation$zone=="");
				if (rock_database$ma_lb[rn] > max(zone_database$ma_lb[zone_nos]) && (nrow(zoneless_formation)==0 || rock_database$ma_lb[rn] > max(zoneless_formation$max_ma)))	{
					if (nrow(zoneless_formation)==0 || max(zone_database$ma_lb[zone_nos]) > max(zoneless_formation$max_ma))	{
						oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
						rock_database$ma_lb[rn] <- zone_database$ma_lb[oldest_zone];
						rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
						} else	{
						oldest_zone <- match(max(zoneless_formation$max_ma),zoneless_formation$max_ma);
						rock_database$ma_lb[rn] <- zoneless_formation$max_ma[oldest_zone];
						rock_database$interval_lb[rn] <- zoneless_formation$early_interval[oldest_zone];
						}
					}
				# if zones make rock older than youngest possible age, then adjust accordingly
				if (rock_database$ma_ub[rn] < min(zone_database$ma_ub[zone_nos]) && (nrow(zoneless_formation)==0 || rock_database$ma_ub[rn] < min(zoneless_formation$min_ma)))	{
					if (nrow(zoneless_formation)==0 || min(zone_database$ma_ub[zone_nos]) < min(zoneless_formation$min_ma))	{
						youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
						rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
						rock_database$interval_ub[rn] <- zone_database$interval_ub[youngest_zone];
						} else	{
						youngest_zone <- match(min(zoneless_formation$min_ma),zoneless_formation$min_ma);
						rock_database$ma_ub[rn] <- zoneless_formation$min_ma[youngest_zone];
						rock_database$interval_ub[rn] <- zoneless_formation$late_interval[youngest_zone];
						}
					}
				}
			formation_members <- unique(this_form$member[this_form$member!=""]);
			fm <- 0;
			while (fm < length(formation_members))	{
				fm <- fm+1;
				whole_name <- paste(zoned_formations[zf]," (",formation_members[fm],")",sep="");
				this_rock_base <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				this_member_base <- subset(this_rock_base,this_rock_base$member==formation_members[fm]);
				memb_zones <- unique(strsplit(paste(this_member_base$zone,collapse="; "),split="; ")[[1]]);
				zoneless_member <- subset(this_member_base,this_member_base$zone=="");
				rnm <- match(whole_name,rock_database$formation_member);
				if (is.na(rnm))	{
					whole_name <- paste(zoned_formations[zf]," (",mundify_rock_unit_names(named_rock_unit=formation_members[fm],dehyphenate=T,delete_rock_type=T),")",sep="");
					rnm <- match(whole_name,rock_database$formation_member);
					}
				if (is.na(rnm))	{
					whole_name <- paste(zoned_formations[zf]," (",mundify_rock_unit_names(named_rock_unit=formation_members[fm],dehyphenate=T,delete_rock_type=T,delete_informal=T),")",sep="");
					rnm <- match(whole_name,rock_database$formation_member);
					}
				if (!is.na(rnm))	{
					zone_nos <- c();
					memb_zones <- memb_zones[memb_zones!=""];
					zn <- 0;
					while (zn < length(memb_zones))	{
						zn <- zn+1;
						zone_nos <- c(zone_nos,unique(which(zone_database==memb_zones[zn],arr.ind = T)[,1]));
						}
					zone_nos <- sort(unique(zone_nos));
					if (length(zone_nos)>0)	{
						# if zones make rock younger than oldest possible age, then adjust accordingly
						if (rock_database$ma_lb[rnm] > max(zone_database$ma_lb[zone_nos]) && (nrow(zoneless_member)==0 || rock_database$ma_lb[rnm] > max(zoneless_member$max_ma)))	{
							if (nrow(zoneless_member)==0 || max(zone_database$ma_lb[zone_nos]) < max(zoneless_member$max_ma))	{
								oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
								rock_database$ma_lb[rnm] <- zone_database$ma_lb[oldest_zone];
								rock_database$interval_lb[rnm] <- zone_database$interval_lb[oldest_zone];
								} else	{
								oldest_zone <- match(max(zoneless_member$max_ma),zoneless_member$max_ma);
								rock_database$ma_lb[rnm] <- zoneless_member$max_ma[oldest_zone];
								rock_database$interval_lb[rnm] <- zoneless_member$early_interval[oldest_zone];
								}
							}
						# if zones make rock older than youngest possible age, then adjust accordingly
						if (rock_database$ma_ub[rnm] < min(zone_database$ma_ub[zone_nos]) && (nrow(zoneless_member)==0 || rock_database$ma_ub[rnm] < min(zoneless_member$min_ma)))	{
							if (nrow(zoneless_member)==0 || min(zone_database$ma_ub[zone_nos]) > min(zoneless_member$min_ma))	{
								youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
								rock_database$ma_ub[rnm] <- zone_database$ma_ub[youngest_zone];
								rock_database$interval_ub[rnm] <- zone_database$interval_ub[youngest_zone];
								} else	{
								youngest_zone <- match(min(zoneless_member$min_ma),zoneless_member$min_ma);
								rock_database$ma_ub[rnm] <- zoneless_member$min_ma[youngest_zone];
								rock_database$interval_ub[rnm] <- zoneless_member$late_interval[youngest_zone];
								}
							}
						}
					}
				}
			}
#		print(zf);
		}

	group_names <- unique(rock_database$group[rock_database$formation_member==""]);
#	group_names <- group_names[group_names!=""];
	gg <- 0;
	while (gg < length(group_names))	{
	#	ff_records <- subset(rock_base,rock_base$formation_rockless==formation_names[ff]);
		gg <- gg+1;
		gg_records <- unique(rbind(subset(rock_base,rock_base$group==group_names[gg]),subset(rock_base,rock_base$group_dehyph==group_names[gg]),subset(rock_base,rock_base$group_rockless==group_names[gg])));
		gb_r <- as.numeric(rownames(gg_records))
		gg_records <- subset(gg_records,gg_records$zone!="");
		gn <- (1:nrow(rock_base))[rock_base$group %in% group_names]
		group_zones <- unique(strsplit(paste(rock_base$zone[gn],collapse="; "),split="; ")[[1]]);
		group_zones <- group_zones[group_zones!=""];
		if (length(group_zones)>0)	{
			zone_nos <- c();
			for (zn in 1:length(group_zones))
				zone_nos <- c(zone_nos,unique(which(zone_database==group_zones[zn],arr.ind = T)[,1]));
			zone_nos <- sort(unique(zone_nos));
			if (length(zone_nos)>0)	{
				rn <- (1:nrow(rock_database))[rock_database$group %in% group_names[gg]];
				rn <- rn[rock_database$formation_member[rn]==""];
				if (!is.na(rn))	{
					# if zones make rock younger than oldest possible age, then adjust accordingly
					if (rock_database$ma_lb[rn] > max(zone_database$ma_lb[zone_nos]))	{
						oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
						# if this is a younger stratigraphic interval than the PaleoDB originally allowed, then change it IF there are not
						#	unzoned collections that make it older
						if (time_scale$ma_ub[match(rock_database$interval_lb[rn],time_scale$interval)] > zone_database$ma_lb[oldest_zone])	{
							other_colls <- (1:nrow(paleodb_collections))[paleodb_collections$stratgroup==group_names[gg]]
							zoneless <- other_colls[paleodb_collections$zone[other_colls]==""];
							if (length(zoneless)>0 && max(time_scale$ma_ub[match(paleodb_collections$early_interval[zoneless],time_scale$interval)]) < zone_database$ma_lb[oldest_zone])	{
								rock_database$ma_lb[rn] <- zone_database$ma_lb[oldest_zone];
								rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
								}
							}
						rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
						}
					# if zones make rock older than youngest possible age, then adjust accordingly
					if (rock_database$ma_ub[rn] < min(zone_database$ma_ub[zone_nos]))	{
						youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
						# if this is an older stratigraphic interval than the PaleoDB originally allowed, then change it IF there are not
						#	unzoned collections that make it younger
						rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
						if (time_scale$ma_lb[match(rock_database$interval_ub[rn],time_scale$interval)] < zone_database$ma_ub[youngest_zone])	{
							other_colls <- (1:nrow(paleodb_collections))[paleodb_collections$stratgroup==group_names[gg]]
							zoneless <- other_colls[paleodb_collections$zone[other_colls]==""];
							if (length(zoneless)>0 && min(time_scale$ma_lb[match(paleodb_collections$late_interval[zoneless],time_scale$interval)]) < zone_database$ma_lb[youngest_zone])	{
								rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
								rock_database$interval_ub[rn] <- zone_database$interval_ub[youngest_zone];
								}
							}
						}
					}
				}
			}
		}
	}
output <- list(rock_database,rocks_to_zones);
names(output) <- c("Rock_Database","Rocks_to_Zones");
return(output);
}

#mmmm <- subset(rocks_to_zones,rocks_to_zones$group=="Honda")
# routine to take rock unit thesaurus and provide numbers for unique rock units
# modified 2020-02-15
# modified 2020-06-08
# 2021-04-15: I had to abandon this because it stopped working for reasons unclear
match_paleodb_collections_to_external_stratigraphic_database <- function(collections,wagner_rocks)	{
n_rocks <- nrow(wagner_rocks);
clean_formations <- wagner_rocks$formation;
named_rock_unit <- wagner_rocks$formation[wagner_rocks$formation!=""];
clean_formations[wagner_rocks$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=FALSE,delete_informal=FALSE);
unique_formations <- sort(unique(clean_formations[clean_formations!=""]));
rock_units <- wagner_rocks$formation;
membered <- (1:n_rocks)[wagner_rocks$member!=""];
rock_units[membered] <- paste(wagner_rocks$formation[membered]," (",wagner_rocks$member[membered],")",sep="");

homonym_rocks <- subset(wagner_rocks,wagner_rocks$rock_no<1);
for (uf in 1:length(unique_formations))	{
	if (sum(wagner_rocks$full_name %in% unique_formations[uf])>1)	{
		xxx <- which(unique_formations[uf]==rock_units,arr.ind=TRUE);
		homonym_rocks <- rbind(homonym_rocks,wagner_rocks[xxx,]);
		}
	}

#rock_units[wagner_rocks$formation=="Haverford Mudstone"]
### Now, take PaleoDB collections and match them to external database
# separate collections with different types of stratigraphic information
n_coll <- nrow(collections);
paleodb_coll_w_formation <- (1:n_coll)[collections$formation!=""];
paleodb_coll_w_member <- (1:n_coll)[collections$member!=""];
paleodb_coll_w_group <- (1:n_coll)[collections$stratgroup!=""];
paleodb_coll_w_formation_and_member <- paleodb_coll_w_member[paleodb_coll_w_member %in% paleodb_coll_w_formation];
paleodb_coll_w_member_only <-  paleodb_coll_w_member[!paleodb_coll_w_member %in% paleodb_coll_w_formation];
paleodb_coll_w_formation_only <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_member];
paleodb_coll_w_group_only <- paleodb_coll_w_group[!paleodb_coll_w_group %in% c(paleodb_coll_w_formation,paleodb_coll_w_member)];
collections_w_rock_names <- sort(unique(c(paleodb_coll_w_formation,paleodb_coll_w_member,paleodb_coll_w_group)))
paleodb_coll_w_rock_names <- length(collections_w_rock_names);

paleodb_coll_entered_rock_unit <- array("",dim=n_coll);
paleodb_coll_entered_rock_unit[paleodb_coll_w_formation_only] <- collections$formation[paleodb_coll_w_formation_only];
paleodb_coll_entered_rock_unit[paleodb_coll_w_member_only] <- collections$member[paleodb_coll_w_member_only];
paleodb_coll_entered_rock_unit[paleodb_coll_w_formation_and_member] <- paste(as.character(collections$formation[paleodb_coll_w_formation_and_member])," (",as.character(collections$member[paleodb_coll_w_formation_and_member]),")",sep="");

### first pass: just check entered names (clean rock units basic)
paleodb_clean_formation_basic <- paleodb_clean_formation_no_rock <- paleodb_clean_formation_no_rock_formal <- collections$formation;
named_rock_unit <- paleodb_clean_formation_basic[paleodb_coll_w_formation];
paleodb_clean_formation_basic[paleodb_coll_w_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal=FALSE);
paleodb_clean_formation_no_rock[paleodb_coll_w_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=FALSE);
paleodb_clean_formation_no_rock_formal[paleodb_coll_w_formation] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

paleodb_clean_member_basic <- paleodb_clean_member_no_rock <- paleodb_clean_member_no_rock_formal <- collections$member;
named_rock_unit <- paleodb_clean_member_basic[paleodb_coll_w_member];
paleodb_clean_member_basic[paleodb_coll_w_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=FALSE,delete_informal=FALSE);
paleodb_clean_member_no_rock[paleodb_coll_w_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=FALSE);
paleodb_clean_member_no_rock_formal[paleodb_coll_w_member] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
# get lists of collections with members, members only, formations+members and formations only
#		after different cleaning
paleodb_coll_w_member_nr <- paleodb_coll_w_member[paleodb_clean_member_no_rock[paleodb_coll_w_member]!=""];
paleodb_coll_w_member_nr_f <- paleodb_coll_w_member[paleodb_clean_member_no_rock_formal[paleodb_coll_w_member]!=""];
paleodb_coll_w_member_only_nr <- paleodb_coll_w_member_nr[paleodb_clean_formation_no_rock[paleodb_coll_w_member_nr]==""];
paleodb_coll_w_member_only_nr_f <- paleodb_coll_w_member_nr_f[paleodb_clean_formation_no_rock_formal[paleodb_coll_w_member_nr_f]==""];
paleodb_coll_w_formation_and_member_nr <- paleodb_coll_w_member_nr[!paleodb_coll_w_member_nr %in% paleodb_coll_w_member_only_nr];
paleodb_coll_w_formation_and_member_nr_f <- paleodb_coll_w_member_nr_f[!paleodb_coll_w_member_nr_f %in% paleodb_coll_w_member_only_nr_f];
paleodb_coll_w_formation_only_nr <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_formation_and_member_nr];
paleodb_coll_w_formation_only_nr_f <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_formation_and_member_nr_f];

paleodb_coll_w_formation_or_member <- sort(unique(c(paleodb_coll_w_formation,paleodb_coll_w_member)));

paleodb_clean_rock_unit_basic <- paleodb_clean_formation_basic;
paleodb_clean_rock_unit_basic[paleodb_coll_w_member_only] <- paleodb_clean_member_basic[paleodb_coll_w_member_only];
paleodb_clean_rock_unit_basic[paleodb_coll_w_formation_and_member] <- paste(paleodb_clean_formation_basic[paleodb_coll_w_formation_and_member]," (",paleodb_clean_member_basic[paleodb_coll_w_formation_and_member],")",sep="");
paleodb_clean_rock_unit_no_rock <- paleodb_clean_formation_no_rock;
paleodb_clean_rock_unit_no_rock[paleodb_coll_w_member_only_nr] <- paleodb_clean_member_no_rock[paleodb_coll_w_member_only_nr];
paleodb_clean_rock_unit_no_rock[paleodb_coll_w_formation_and_member_nr] <- paste(paleodb_clean_formation_no_rock[paleodb_coll_w_formation_and_member_nr]," (",paleodb_clean_member_no_rock[paleodb_coll_w_formation_and_member_nr],")",sep="");
paleodb_clean_rock_unit_no_rock_formal <- paleodb_clean_formation_no_rock_formal;
paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_member_only_nr_f] <- paleodb_clean_member_no_rock_formal[paleodb_coll_w_member_only_nr_f];
paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f] <- paste(paleodb_clean_formation_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f]," (",paleodb_clean_member_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f],")",sep="");

paleodb_clean_group_basic <- paleodb_clean_group_no_rock <- paleodb_clean_group_no_rock_formal <- collections$stratgroup;
named_rock_unit <- paleodb_clean_group_basic[paleodb_coll_w_group];
paleodb_clean_group_basic[paleodb_coll_w_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal=FALSE);
paleodb_clean_group_no_rock[paleodb_coll_w_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=FALSE);
paleodb_clean_group_no_rock_formal[paleodb_coll_w_group] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

# if exact same name entered into formation & member, then make it just formation name.
same_name <- paleodb_coll_w_formation_and_member[paleodb_clean_formation_basic[paleodb_coll_w_formation_and_member]==paleodb_clean_member_basic[paleodb_coll_w_formation_and_member]];
paleodb_clean_member_no_rock_formal[same_name] <- paleodb_clean_member_no_rock[same_name] <- paleodb_clean_member_basic[same_name]<- "";
paleodb_clean_rock_unit_no_rock_formal[same_name] <- paleodb_clean_formation_no_rock_formal[same_name];
paleodb_clean_rock_unit_no_rock[same_name] <- paleodb_clean_formation_no_rock[same_name];
paleodb_clean_rock_unit_basic[same_name] <- paleodb_clean_formation_basic[same_name];

# set up new fields to be added
rock_unit_senior2 <- rock_unit_senior <- interval_lb <- interval_ub <- array("?",dim=n_coll);
ma_lb <- ma_ub <- rock_no <- rock_no_sr <- formation_no <- rock2_no <- rock2_no_sr <- formation2_no <- clean_match <- array(0,dim=n_coll);

# test 1: does name match any of the basic full names?
test_cols <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock"),colnames(wagner_rocks));
examined_matrix <- wagner_rocks[,test_cols];
for (i in 1:ncol(examined_matrix))	examined_matrix[,i] <- tolower(examined_matrix[,i]);

#find_me <- "Haverford";
find_me <- paleodb_clean_rock_unit_basic[paleodb_coll_w_formation_or_member];
external_rock_rows <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
find_me <- paleodb_clean_rock_unit_no_rock[paleodb_coll_w_formation_or_member];
e_r_r_2 <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
find_me <- paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_formation_or_member];
e_r_r_3 <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
external_rock_rows[external_rock_rows < 1] <- e_r_r_2[external_rock_rows < 1];
external_rock_rows[external_rock_rows < 1] <- e_r_r_3[external_rock_rows < 1];
# positive integers in testing are row numbers from stratigraphic database for matching rocks
# -1 means that no matches were made.
# -2 means that 2+ matches were made.
# paleodb_coll_w_formation_or_member[testing>0] gives row numbers of paleodb collections for which we've matched rocks
emended_rocks <- external_rock_rows[external_rock_rows>0];
updated_collections <- paleodb_coll_w_formation_or_member[external_rock_rows>0];
ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[external_rock_rows>0]];
ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[external_rock_rows>0]];
rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[external_rock_rows>0]];
rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[external_rock_rows>0]];
formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[external_rock_rows>0]];
rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[external_rock_rows>0]];
interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[external_rock_rows>0]];
interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[external_rock_rows>0]];
clean_match[updated_collections] <- 1;

# now, deal with stubborn members & formations....
peles <- sort(c(paleodb_coll_w_formation_only,paleodb_coll_w_member_only));	# one name units
need_info_still <- peles[ma_lb[peles]==0];
n_i_p <- length(need_info_still);
test_cols_p <- match(c("group","formation","member","group_clean_basic","formation_clean_basic","member_clean_basic","group_clean_no_rock","formation_clean_no_rock","member_clean_no_rock","group_clean_no_rock_formal","formation_clean_no_rock_formal","member_clean_no_rock_formal"),colnames(wagner_rocks));
examined_matrix_p <- wagner_rocks[,test_cols_p];
for (i in 1:length(test_cols_p))	examined_matrix_p[,i] <- tolower(examined_matrix_p[,i]);

test_cols_g <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),colnames(wagner_rocks));
examined_matrix_g <- wagner_rocks[,test_cols_g];
for (i in 1:length(test_cols_g))	examined_matrix_g[,i] <- tolower(examined_matrix_g[,i]);
#matrix_row_match("el paso",examined_matrix_p)

find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
external_rock_rows_p <- sapply(find_me,matrix_row_match,examined_matrix_p);
find_me <- tolower(paleodb_clean_rock_unit_no_rock[need_info_still]);
e_r_r_p_2 <- sapply(find_me,matrix_row_match,examined_matrix_p);
e_r_r_p_2[find_me==""] <- -1;
external_rock_rows_p[external_rock_rows_p<1] <- e_r_r_p_2[external_rock_rows_p<1];
find_me <- tolower(paleodb_clean_rock_unit_no_rock_formal[need_info_still]);
e_r_r_p_3 <- sapply(find_me,matrix_row_match,examined_matrix_p);
e_r_r_p_3[find_me==""] <- -1;
external_rock_rows_p[external_rock_rows_p<1] <- e_r_r_p_3[external_rock_rows_p<1];
#cbind(need_info_still,external_rock_rows_p)
# collections with members with 2+ ids in external database.
#still_need <- (1:length(need_info_member))[external_rock_rows_m==-2];
for_homer <- (1:n_i_p)[external_rock_rows_p==-2];
two_plus <- need_info_still[external_rock_rows_p==-2];
find_me_2 <- tolower(paleodb_clean_rock_unit_no_rock_formal[two_plus]);
two_plus_names <- tolower(unique(find_me_2));
# use try_again again!
tp <- 0;
while (tp < length(two_plus_names))	{
#for (tp in 1:length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	if (nrow(xxx)>0)	{
		yyy <- unique(xxx[,1]);
		if (length(unique(wagner_rocks$rock_no_sr[yyy]))==1)	{
			# gets cases where rock unit is a member in external data base
			# get needy collections with this as member
			tried_m <- (1:n_i_p)[tolower(paleodb_clean_member_no_rock_formal[need_info_still])==two_plus_names[tp]];
			tried_f <- (1:n_i_p)[tolower(paleodb_clean_formation_no_rock_formal[need_info_still])==two_plus_names[tp]];
			tried <- unique(c(tried_m,tried_f));
			external_rock_rows_p[tried] <- match(wagner_rocks$rock_no_sr[yyy[1]],wagner_rocks$rock_no)
			} else if (length(unique(wagner_rocks$formation_no[yyy]))==1)	{
			# get needy collections with this as formation
			# gets cases where rock unit is a formation in external data base
			tried <- (1:n_i_p)[tolower(paleodb_clean_formation_no_rock_formal[need_info_still])==two_plus_names[tp]];
			external_rock_rows_p[tried] <- match(wagner_rocks$formation_no[yyy[1]],wagner_rocks$rock_no);
			} else	{
			vvv <- which(examined_matrix_g==two_plus_names[tp],arr.ind=TRUE);
			uuu <- unique(vvv[,1]);
#			if (length(uuu[wagner_rocks$formation[uuu]=="•"])==1)	{
			if (length(uuu[wagner_rocks$formation[uuu]==""])==1)	{
				# gets cases where rock unit is a group in external data base
				tried <- (1:n_i_p)[find_me %in% two_plus_names[tp]];
#				tried <- (1:n_i_p)[tolower(paleodb_clean_group_no_rock_formal[need_info_still])==two_plus_names[tp]];
				external_rock_rows_p[tried] <- match(uuu[wagner_rocks$formation[uuu]==""],wagner_rocks$rock_no);
#				external_rock_rows_p[tried] <- match(uuu[wagner_rocks$formation[uuu]=="•"],wagner_rocks$rock_no);
				}
			}
		}
#	print(cbind(paleodb_clean_rock_unit_basic[two_plus],external_rock_rows_p[for_homer]))
	}

emended_rocks <- (1:n_i_p)[external_rock_rows_p>0];
emended_rocks <- emended_rocks[!is.na(emended_rocks)];	# KLUGE!!!!
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections) > 0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows_p[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows_p[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows_p[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows_p[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows_p[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows_p[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows_p[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows_p[emended_rocks]];
	clean_match[updated_collections] <- 1;
	}
### done with stubborn rocks

# look at things with just group ids
need_info_still <- collections_w_rock_names[ma_lb[collections_w_rock_names]==0];
need_info_still_f <- need_info_still[paleodb_clean_formation_basic[need_info_still]!=""];
need_info_still <- need_info_still[!need_info_still %in% need_info_still_f];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_group_no_rock_formal[need_info_still]);
external_rock_rows <- sapply(find_me,matrix_row_match,examined_matrix_p);
#cbind(paleodb_clean_group_no_rock_formal[need_info_still_g],external_rock_rows_g);
for_bart <- (1:n_i_s)[external_rock_rows==-2];
two_plus <- need_info_still[external_rock_rows==-2];
find_me <- tolower(paleodb_clean_group_no_rock_formal[two_plus]);
two_plus_names <- tolower(unique(find_me));
# use try_again again!
tp <- 0;
while (tp < length(two_plus_names))	{
#for (tp in 1:length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	yyy <- unique(xxx[,1]);
	vvv <- c();
#	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="• (¢)"])==1)	{
#		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="• (¢)")];
	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==""])==1)	{
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="")];
		}	else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])==1)	{
		# failing that, look to see if the group name is a stand-alone formation name
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])];
		} else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])>1)	{
		vvv <- yyy[wagner_rocks$member[yyy]==""];
		}
	# get all of the collections with this rock unit
	if (length(vvv)==1)	{
		tried <- (1:n_i_s)[tolower(paleodb_clean_group_no_rock_formal[need_info_still])==two_plus_names[tp]];
		external_rock_rows[tried] <- vvv[1];
		}
#	print(cbind(paleodb_clean_group_basic[two_plus],external_rock_rows_g[for_bart]));
	}

emended_rocks <- (1:n_i_s)[external_rock_rows>0];
emended_rocks <- emended_rocks[!is.na(emended_rocks)];		# KLUGE!!!!!
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections)>0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[emended_rocks]];
	clean_match[updated_collections] <- 1;
	}
# we are getting near the point of just burning it all and making it look like an electrical thing
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
external_rock_rows <- sapply(find_me,matrix_row_match,examined_matrix);
for_bender <- (1:n_i_s)[external_rock_rows==-2];
two_plus <- need_info_still[external_rock_rows==-2];
find_me <- tolower(paleodb_clean_rock_unit_basic[two_plus]);
two_plus_names <- tolower(unique(find_me));
tp <- 0;
#for (tp in 1:length(two_plus_names))	{
while (tp < length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	yyy <- unique(xxx[,1]);
	vvv <- c();
	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==""])==1)	{
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="")];
		}	else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])==1)	{
		# failing that, look to see if the group name is a stand-alone formation name
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])];
		} else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])>1)	{
		vvv <- yyy[wagner_rocks$member[yyy]==""];
		}
	# get all the collections with this rock unit
	if (length(vvv)==1)	{
		tried <- (1:n_i_s)[tolower(paleodb_clean_rock_unit_no_rock_formal[need_info_still])==two_plus_names[tp]];
		external_rock_rows[tried] <- vvv[1];
		}
	}

emended_rocks <- (1:n_i_s)[external_rock_rows>0];
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections)>0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[emended_rocks]];
	clean_match[updated_collections] <- 0;
	}

# I am getting stabby, so let's cut up some units with 2+ names.....
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#funky_text <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#find_me <- sapply(funky_text,transmogrify_diacritics);
cc <- 0;
#for (cc in 1:n_i_s)	{
while (cc < n_i_s)	{
	cc <- cc+1;
	named_rock_unit <- collections$formation[need_info_still[cc]];
	named_rock_unit <- divido_lumped_rock_units(named_rock_unit);
	if (length(named_rock_unit)>1)	{
		named_rock_unit <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
		find_me2 <- tolower(named_rock_unit);
		external_rock_rows <- sapply(find_me2,matrix_row_match,examined_matrix_p);
		if (sum(external_rock_rows==-2)>0)	{
			two_plus_names <- find_me2[external_rock_rows==-2];
			for (tp in 1:length(two_plus_names))	{
				xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
				yyy <- unique(xxx[,1]);
				zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$formation_no[yyy]];
				if (length(zzz)==0)
					zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$rock_no_sr[yyy]];
				if (length(zzz)==1)	{
					external_rock_rows[match(two_plus_names[tp],find_me2)] <- zzz[1];
					} else	{
					vvv <- yyy[wagner_rocks$formation[yyy]=="" && wagner_rocks$member[yyy]==""];
					if (length(vvv)==1)	{
						external_rock_rows[match(two_plus_names[tp],find_me2)] <- vvv[1];
						} else	{
						uuu <- yyy*(wagner_rocks$formation[yyy]==wagner_rocks$group[yyy]) * (wagner_rocks$member[yyy]=="");
						uuu <- uuu[uuu>1];
						if (length(uuu)==1)	{
							external_rock_rows[match(two_plus_names[tp],find_me2)] <- uuu[1];
							}
						}
					}
				}
			}
		if (sum(external_rock_rows>0)==1)	{
			# do a thing
			external_rock_rows <- external_rock_rows[external_rock_rows>0];
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]]
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]]
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]]
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[1]];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[1]];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows[1]];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows[1]];
			clean_match[need_info_still[cc]] <- 1;
			} else if (sum(external_rock_rows>0)>1)	{
			# do another thing
			external_rock_rows <- external_rock_rows[external_rock_rows>0];
			if (length(external_rock_rows)>2)	{
				e_r_r <- 1:length(external_rock_rows);
				aa <- match(max(wagner_rocks$ma_lb[external_rock_rows]),wagner_rocks$ma_lb[external_rock_rows]);
				r_e_r_r <- e_r_r[!e_r_r %in% aa];
				zz <- r_e_r_r[match(min(wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]]),wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]])];
				mm <- r_e_r_r[!r_e_r_r %in% zz];
				external_rock_rows <- external_rock_rows[c(aa,zz,mm)];
				}
			ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1:2]]);
			ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1:2]]);
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
			rock2_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[2]];
			rock2_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[2]];
			formation2_no[need_info_still[cc]] <-  wagner_rocks$formation_no[external_rock_rows[2]];
#			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[2]];
#			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[2]],wagner_rocks$rock_no)];
			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[2]];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[match(max(wagner_rocks$ma_lb[external_rock_rows[1:2]]),wagner_rocks$ma_lb[external_rock_rows[1:2]])]];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[match(min(wagner_rocks$ma_ub[external_rock_rows[1:2]]),wagner_rocks$ma_ub[external_rock_rows[1:2]])]];
			clean_match[need_info_still[cc]] <- 1;
			}
		} else	{
		named_rock_unit <- collections$member[need_info_still[cc]];
		named_rock_unit <- divido_lumped_rock_units(named_rock_unit);
		if (length(named_rock_unit)>1)	{
			named_rock_unit <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
			find_me2 <- tolower(named_rock_unit);
			external_rock_rows <- sapply(find_me2,matrix_row_match,examined_matrix_p);
			if (sum(external_rock_rows==-2)>0)	{
				two_plus_names <- find_me2[external_rock_rows==-2];
				for (tp in 1:length(two_plus_names))	{
					xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
					yyy <- unique(xxx[,1]);
					zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$formation_no[yyy]];
					if (length(zzz)==0)
						zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$rock_no_sr[yyy]];
					if (length(zzz)==1)	{
						external_rock_rows[match(two_plus_names[tp],find_me2)] <- zzz[1];
						} else	{
						vvv <- yyy[wagner_rocks$formation[yyy]=="" && wagner_rocks$member[yyy]==""];
						if (length(vvv)==1)	{
							external_rock_rows[match(two_plus_names[tp],find_me2)] <- vvv[1];
							} else	{
							uuu <- yyy*(wagner_rocks$formation[yyy]==wagner_rocks$group[yyy]) * (wagner_rocks$member[yyy]=="");
							uuu <- uuu[uuu>1];
							if (length(uuu)==1)	{
								external_rock_rows[match(two_plus_names[tp],find_me2)] <- uuu[1];
								}
							}
						}
					}
				}
			if (sum(external_rock_rows>0)==1)	{
				# do a thing
				external_rock_rows <- external_rock_rows[external_rock_rows>0];
				ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1]]);
				ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1]]);
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]]
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]]
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]]
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[1]];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[1]];
				clean_match[need_info_still[cc]] <- 1;
				} else if (sum(external_rock_rows>0)>1)	{
				# do another thing
				external_rock_rows <- external_rock_rows[external_rock_rows>0];
				if (length(external_rock_rows)>2)	{
					e_r_r <- 1:length(external_rock_rows);
					aa <- match(max(wagner_rocks$ma_lb[external_rock_rows]),wagner_rocks$ma_lb[external_rock_rows]);
					r_e_r_r <- e_r_r[!e_r_r %in% aa];
					zz <- r_e_r_r[match(min(wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]]),wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]])];
					mm <- r_e_r_r[!r_e_r_r %in% zz];
					external_rock_rows <- external_rock_rows[c(aa,zz,mm)];
					}
				ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1:2]]);
				ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1:2]]);
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
				rock2_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[2]];
				rock2_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[2]];
				formation2_no[need_info_still[cc]] <-  wagner_rocks$formation_no[external_rock_rows[2]];
#				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[2]];
#				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[2]],wagner_rocks$rock_no)];
				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[2]];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[match(max(wagner_rocks$ma_lb[external_rock_rows[1:2]]),wagner_rocks$ma_lb[external_rock_rows[1:2]])]];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[match(min(wagner_rocks$ma_ub[external_rock_rows[1:2]]),wagner_rocks$ma_ub[external_rock_rows[1:2]])]];
				clean_match[need_info_still[cc]] <- 1;
				}
			}
		}
	}
#cbind(need_info_still,collections$collection_no[need_info_still],paleodb_clean_rock_unit_basic[need_info_still],ma_lb[need_info_still])

#### last of it ???
pcwfom <- length(paleodb_coll_w_formation_or_member);
#(1:pcwfom)[is.na(ma_lb[paleodb_coll_w_formation_or_member])]
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#for (cc in 1:n_i_s)	{
cc <- 0;
while (cc<n_i_s)	{
	cc <- cc+1;
	if (paleodb_clean_formation_no_rock_formal[need_info_still[cc]]!="" && paleodb_clean_member_no_rock_formal[need_info_still[cc]]!="")	{
		find_me2 <- paleodb_clean_member_no_rock_formal[need_info_still[cc]];
		external_rock_rows_m <- sapply(tolower(find_me2),matrix_row_match,examined_matrix_p);
		find_me3 <- paleodb_clean_formation_no_rock_formal[need_info_still[cc]];
		external_rock_rows_f <- sapply(tolower(find_me3),matrix_row_match,examined_matrix_p);
		if (external_rock_rows_f>0 && external_rock_rows_m==-1)	{
#			print("One possible formation");
			external_rock_rows <- external_rock_rows_f;
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
			} else if (external_rock_rows_f==-2 && external_rock_rows_m==-1)	{
#			print("Multiple possible formations");
			xxx <- which(examined_matrix_p==tolower(find_me3),arr.ind=TRUE);
			yyy <- unique(xxx[,1]);
			yyy <- yyy[wagner_rocks$rock_no[yyy]==wagner_rocks$formation_no[yyy]];
			if (length(yyy)==0)
				yyy <- yyy[wagner_rocks$rock_no_sr[yyy]==wagner_rocks$formation_no[yyy]];
			if (length(yyy)>1)
				yyy <- yyy[tolower(wagner_rocks$formation[yyy])==tolower(find_me3)]
			if (length(yyy)==1)	{
				external_rock_rows <- yyy[1];
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
				ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
				ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
				}
			} else if (external_rock_rows_m>0)	{
#			print("One possible member");
			external_rock_rows <- external_rock_rows_m;
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
			} else if (external_rock_rows_m==-2)	{
#			print("Multiple possible members");
			xxx <- which(examined_matrix_p==tolower(find_me2),arr.ind=TRUE);
			yyy <- unique(xxx[,1]);
			yyy <- yyy[wagner_rocks$rock_no[yyy]==wagner_rocks$wagner_rocks$rock_no_sr[yyy]];
			if (length(yyy)==0)	{
				yyy <- unique(xxx[,1]);
				yyy <- yyy[wagner_rocks$rock_unit_clean_no_rock_formal[yyy]==find_me2];
				}
			if (length(yyy)>1)
				yyy <- yyy[tolower(wagner_rocks$member[yyy])==tolower(find_me2)]
			if (length(yyy)==1)	{
				external_rock_rows <- yyy[1];
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
				ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
				ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
				}
			}
#		print(cbind(find_me[1:cc],ma_lb[need_info_still[1:cc]]));
		}
	}

# I've had it: I'm going Arya on any last undated collections!!!!
collections$max_ma <- as.numeric(collections$max_ma);
collections$min_ma <- as.numeric(collections$min_ma);
ncolls <- nrow(collections);
# I do not know why this is needed; a mystery NA appears if I just set the vectors equal
#ma_lb[ma_lb==0] <- as.numeric(collections$max_ma[ma_lb==0]);
#ma_ub[ma_ub==0] <- as.numeric(collections$min_ma[ma_ub==0]);
undated <- (1:ncolls)[ma_lb==0];
for (ud in 1:length(undated))
	ma_lb[undated[ud]] <- as.numeric(collections$max_ma[undated[ud]]);
undated <- (1:ncolls)[ma_ub==0];
for (ud in 1:length(undated))
	ma_ub[undated[ud]] <- as.numeric(collections$min_ma[undated[ud]]);

homonym_names <- unique(homonym_rocks$formation);
for (hr in 1:length(homonym_names))	{
	these_cases <- subset(homonym_rocks,homonym_rocks$formation==homonym_names[hr])
	cases <- unique(c((1:ncolls)[paleodb_clean_rock_unit_no_rock_formal %in% homonym_names[hr]],
			(1:ncolls)[paleodb_clean_rock_unit_no_rock %in% homonym_names[hr]],
			(1:ncolls)[paleodb_clean_rock_unit_basic %in% homonym_names[hr]]));
	if (length(unique(these_cases$geoplate))==nrow(these_cases))	{
		ddd <- match(as.numeric(collections$geoplate[cases]),these_cases$geoplate);
		if (sum(is.na(ddd))>0)	{
			eee <- ddd[is.na(ddd)];
			}
		cases <- cases[!is.na(ddd)];
		ddd <- ddd[!is.na(ddd)];
		rock_no[cases] <- these_cases$rock_no[ddd];
		rock_no_sr[cases] <- these_cases$rock_no_sr[ddd];
		formation_no[cases] <- these_cases$formation_no[ddd];
		interval_lb[cases] <- these_cases$interval_lb[ddd];
		interval_ub[cases] <- these_cases$interval_ub[ddd];
		ma_lb[cases[ma_lb[cases]>these_cases$ma_lb[ddd]]] <- these_cases$ma_lb[ddd[ma_lb[cases]>these_cases$ma_lb[ddd]]];
		ma_ub[cases[ma_ub[cases]<these_cases$ma_ub[ddd]]] <- these_cases$ma_ub[ddd[ma_ub[cases]<these_cases$ma_ub[ddd]]];
		} else	{
		### do something here, Pete!
		}
	}

collections$rock_no_sr <- as.numeric(rock_no_sr);
collections$rock_no <- as.numeric(rock_no);
collections$formation_no <- as.numeric(formation_no);
collections$rock_unit_senior <- as.character(rock_unit_senior);
collections$rock2_no_sr <- as.numeric(rock2_no_sr);
collections$rock2_no <- as.numeric(rock2_no);
collections$formation2_no <- as.numeric(formation2_no);
collections$rock_unit_senior2 <- as.character(rock_unit_senior2);
collections$ma_lb <- as.numeric(ma_lb);
collections$ma_ub <- as.numeric(ma_ub);
collections$interval_lb <- as.character(interval_lb);
collections$interval_ub <- as.character(interval_ub);
collections$clean_match <- as.numeric(clean_match);
#poss_subset <- subset(wagner_rocks,wagner_rocks$formation_no==collections$rock_no_sr[worsened_collections[1]]);

collections$interval_lb[collections$interval_lb=="?"] <- as.character(collections$early_interval[interval_lb=="?"]);
collections$interval_ub[collections$interval_ub=="?"] <- as.character(collections$late_interval[interval_ub=="?"]);

### find any dates that are much too broad
# um, finest_chronostrat isn't in this routine, pete....
#overbroad_collections <- subset(collections,collections$ma_lb>=collections$max_ma & collections$ma_ub<=collections$min_ma);
#age <- overbroad_collections$ma_lb <- overbroad_collections$max_ma;
#overbroad_collections$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale = finest_chronostrat);
#age <- overbroad_collections$ma_ub <- overbroad_collections$min_ma;
#overbroad_collections$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale = finest_chronostrat);
#collections[match(overbroad_collections$collection_no,collections$collection_no),] <- overbroad_collections;

#dim(overbroad_collections)
#coll_spans_new <- abs(collections$ma_lb-collections$ma_ub);
#coll_spans_old <- abs(collections$max_ma-collections$min_ma);
#coll_spans_dif <- coll_spans_new-coll_spans_old;
#hist(coll_spans_new)
#hist(coll_spans_old)
#ggg <- seq(-40,60,by=5);
#ddd <- hist(coll_spans_new-coll_spans_old,breaks=ggg)$counts;
#names(ddd) <- ggg[1:(length(ggg)-1)];
#worsened_collections <- (1:nrow(collections))[order(coll_spans_dif,decreasing=T)];

#redone_collections <- cbind(collections,rock_no_sr,rock_no,formation_no,rock_unit_senior,rock2_no_sr,rock2_no,formation2_no,rock_unit_senior2,ma_lb,ma_ub,interval_lb,interval_ub,clean_match);
#redone_collections <- collections;
return(collections);
}

revelare_undocumented_rock_units <- function(coll_data)	{
ncoll <- nrow(coll_data);
unidentified_rocks <- (1:ncoll)[coll_data$rock_unit_senior=="?"];
sneaky_formations <- unidentified_rocks[coll_data$formation[unidentified_rocks]!=""];
sneaky_members <- unidentified_rocks[coll_data$member[unidentified_rocks]!=""];
sneaky_rocks <- sort(unique(c(sneaky_formations,sneaky_members)));

undocumented_rocks <- data.frame(formation=as.character(coll_data$formation[sneaky_rocks]),member=as.character(coll_data$member[sneaky_rocks]),early_interval=as.character(coll_data$early_interval[sneaky_rocks]),late_interval=as.character(coll_data$late_interval[sneaky_rocks]),zone=as.character(coll_data$zone[sneaky_rocks]),stringsAsFactors = FALSE);
undocumented_rocks <- unique(undocumented_rocks);
rock_unit <- as.character(undocumented_rocks$formation);
has_member <- (1:nrow(undocumented_rocks))[undocumented_rocks$member!=""];
rock_unit[has_member] <- paste(undocumented_rocks$formation[has_member]," (",undocumented_rocks$member[has_member],")",sep="");
has_member_no_formation <- has_member[undocumented_rocks$formation[has_member]==""];
rock_unit[has_member_no_formation] <- as.character(undocumented_rocks$member[has_member_no_formation]);
undocumented_rocks <- cbind(rock_unit=as.character(rock_unit),undocumented_rocks);
rock_units <- sort(unique(rock_unit));

rock_units_all <- coll_data$rock_unit_senior;
rock_units_all[sneaky_formations] <- coll_data$formation[sneaky_formations]
rock_units_all[sneaky_members]  <- paste(coll_data$formation[sneaky_members]," (",coll_data$member[sneaky_members],")",sep="");
rock_units_all[sneaky_members[!sneaky_members %in% sneaky_formations]] <- coll_data$member[sneaky_members[!sneaky_members %in% sneaky_formations]];

#doh <- c();
#for (i in 1:length(rock_units))
#	doh <- c(doh,sum(rock_units[i]==undocumented_rocks$rock_unit))
undocumented_rock_summary <- data.frame(stringsAsFactors = FALSE);
#rr <- 0;
for (rr in 1:length(rock_units))	{
#	dumb_ass <- subset(undocumented_rocks,undocumented_rocks$rock_unit==rock_units[rr]);
#	relv_rocks <- as.character(relv_rocks);
#	rr <- rr+1;
	relv_rocks <- subset(undocumented_rocks,undocumented_rocks$rock_unit==rock_units[rr]);
	relv_coll <- (1:ncoll)[rock_units_all==as.character(relv_rocks$rock_unit[1])];
	relv_coll_no <- paste(coll_data$collection_no[relv_coll],collapse=";");
	relv_coll_no <- rep(relv_coll_no,nrow(relv_rocks));
	ttl_colls <- length(relv_coll);
	relv_rocks <- cbind(relv_rocks,ttl_colls,relv_coll_no);
	if (nrow(relv_rocks)==1)	{
#		if (rr==1)	{
#			undocumented_rock_summary <- relv_rocks;
#			} else	{
			undocumented_rock_summary <- rbind(undocumented_rock_summary,relv_rocks);
#			}
		} else	{
		earlies <- unique(relv_rocks$early_interval);
		earlies <- paste(earlies,collapse=", ");
		laters <- unique(relv_rocks$late_interval[relv_rocks$late_interval!=""]);
		laters <- paste(laters,collapse=", ");
		zones <- unique(relv_rocks$zone);
		zones <- zones[zones!=""];
		if (length(zones)>1)	{
			zones <- paste(zones,collapse=", ");
			} else if (length(zones)==0)	{
			zones <- "";
			}
#		relv_rocks$early_interval[1] <- as.character(relv_rocks$early_interval[1])
#		relv_rocks$early_interval[1] <- relv_rocks$early_interval[1];
		relv_rocks$early_interval[1] <- earlies;
		relv_rocks$late_interval[1] <- laters;
		relv_rocks$zone[1] <- zones;

#		if (rr==1)	{
#			undocumented_rock_summary <- relv_rocks[1,];
#			} else	{
			undocumented_rock_summary <- rbind(undocumented_rock_summary,relv_rocks[1,]);
#			}
		}
#	print(c(rr,rock_units[rr]));
	}
return(undocumented_rock_summary);
}

# tally numbers of rock units
# modified 2020-02-15
number_unique_rock_units <- function(paleodb_collections,zone_database,time_scale)	{
paleodb_rocks_info <- construct_stratigraphic_data_base_from_paleodb(paleodb_collections=paleodb_collections,zone_database,time_scale);
paleodb_rocks <- paleodb_rocks_info$Rock_Database;
paleodb_rocks_to_zones <- paleodb_rocks_info$Rocks_to_Zones;
paleodb_rock_thesaurus <- accersi_rock_unit_thesaurus_given_paleodb_collections(paleodb_rocks);
paleodb_collections <- expello_na_from_matrix(paleodb_collections,replacement="");
paleodb_collections$rock_no <- as.numeric(match_paleodb_collections_to_paleodb_rock_thesaurus(paleodb_collections=paleodb_collections,paleodb_rock_thesaurus));
paleodb_collections$rock_no_sr <- paleodb_rock_thesaurus$rock_no_sr[match(paleodb_collections$rock_no,paleodb_rock_thesaurus$rock_no)];
paleodb_collections$formation_no <- paleodb_rock_thesaurus$formation_no[match(paleodb_collections$rock_no,paleodb_rock_thesaurus$rock_no)];
paleodb_collections$formation_no_sr <- paleodb_rock_thesaurus$formation_no_sr[match(paleodb_collections$rock_no_sr,paleodb_rock_thesaurus$rock_no)];
paleodb_collections <- expello_na_from_matrix(paleodb_collections,replacement=0);
return(paleodb_collections);
}

					##### Redone routines after some of the above stopped working for unknown reasons ######
#wagner_rocks <- rock_database;	paleodb_collections <- trilo_sites_reset;
organize_pbdb_rock_data <- function(paleodb_collections,geosplit=F,max_gap=25)	{
# max_gap: biggest gap allowed for two sets of sites with same name on same continent to be considered the same rock unit
print("Allocating rock thesaurus data.frame.")
pbdb_rocks <- unique(data.frame(group=as.character(paleodb_collections$stratgroup),
								formation=as.character(paleodb_collections$formation),
								member=as.character(paleodb_collections$member),stringsAsFactors = F));
if (!is.null(paleodb_collections$formation_alt))	{
	dummy <- unique(data.frame(group=as.character(paleodb_collections$stratgroup_alt),
							   formation=as.character(paleodb_collections$formation_alt),
							   member=as.character(paleodb_collections$member_alt),stringsAsFactors = F));
	pbdb_rocks <- unique(rbind(pbdb_rocks,dummy));
	}

#print("Ordering Rock Units.")
pbdb_rocks <- pbdb_rocks[order(pbdb_rocks$formation,pbdb_rocks$member,pbdb_rocks$group),];
#pbdb_rocks[1:4,]
nrocks <- nrow(pbdb_rocks);
gg <- (1:nrocks)[pbdb_rocks$group!=""]; ff <- (1:nrocks)[pbdb_rocks$formation!=""]; mm <- (1:nrocks)[pbdb_rocks$member!=""];
pbdb_rocks <- pbdb_rocks[sort(unique(c(gg,ff,mm))),];
nrocks <- nrow(pbdb_rocks);
gg <- (1:nrocks)[pbdb_rocks$group!=""]; ff <- (1:nrocks)[pbdb_rocks$formation!=""]; mm <- (1:nrocks)[pbdb_rocks$member!=""];
pbdb_rocks$full_name <- pbdb_rocks$formation;
pbdb_rocks$full_name[mm[mm %in% ff]] <- paste(pbdb_rocks$formation[mm[mm %in% ff]]," (",pbdb_rocks$member[mm[mm %in% ff]],")",sep="");
pbdb_rocks$full_name[mm[!mm %in% ff]] <- paste("(",pbdb_rocks$member[mm[!mm %in% ff]],")",sep="");
pbdb_rocks$full_name[gg[!gg %in% c(ff,mm)]] <- paste("[",pbdb_rocks$group[gg[!gg %in% c(ff,mm)]],"]",sep="");

#print("Construct full names for rock units");
nsites <- nrow(paleodb_collections);
paleodb_collections$full_name <- paleodb_collections$formation;
paleodb_collections$full_name[paleodb_collections$member!=""] <- paste(paleodb_collections$formation[paleodb_collections$member!=""]," (",paleodb_collections$member[paleodb_collections$member!=""],")",sep="");
members_only <- (1:nsites)[(paleodb_collections$member!="") & (paleodb_collections$formation=="")];
paleodb_collections$full_name[members_only] <- paste("(",paleodb_collections$member[members_only],")",sep="");
group_only <- (1:nsites)[paleodb_collections$stratgroup!="" & (paleodb_collections$member=="") & (paleodb_collections$formation=="")];
paleodb_collections$full_name[group_only] <- paste("[",paleodb_collections$stratgroup[group_only],"]",sep="");

unique_full_names <- unique(pbdb_rocks$full_name);
kill_me <- dups <- c();
for (ur in 1:length(unique_full_names))	{
	if (sum(pbdb_rocks$full_name %in% unique_full_names[ur])>1)	{
		this_case <- subset(pbdb_rocks,pbdb_rocks$full_name %in% unique_full_names[ur]);
		if (sum(this_case$group=="")>0)	{
			kill_me <- c(kill_me,(1:nrocks)[pbdb_rocks$full_name %in% unique_full_names[ur]][pbdb_rocks$group[pbdb_rocks$full_name %in% unique_full_names[ur]]==""]);
			} # end case where on rock isn't put into a rock uniit.
		if (sum(this_case$group!="")>1)	{
			this_case <- subset(this_case,this_case$group!="");
			these_rock_collections <- subset(paleodb_collections,paleodb_collections$full_name==unique_full_names[ur]);
			these_rock_collections <- subset(these_rock_collections,these_rock_collections$stratgroup!="");
			relv_groups <- unique(these_rock_collections$stratgroup);
			for (rg in 1:(length(relv_groups)-1))	{
				gr <- rg;
				lb_a <- these_rock_collections$max_ma[these_rock_collections$stratgroup==relv_groups[rg]];
				ub_a <- these_rock_collections$min_ma[these_rock_collections$stratgroup==relv_groups[rg]];
				while (gr < length(relv_groups))	{
					gr <- gr+1;
					lb_b <- these_rock_collections$max_ma[these_rock_collections$stratgroup==relv_groups[gr]];
					ub_b <- these_rock_collections$min_ma[these_rock_collections$stratgroup==relv_groups[gr]];
					if (do_two_ranges_overlap(lb_a,ub_a,lb_b,ub_b))	{
						# if they overlap, then eliminate the one last used longest ago
						if (max(these_rock_collections$ref_pubyr[these_rock_collections$stratgroup==relv_groups[rg]])>max(these_rock_collections$ref_pubyr[these_rock_collections$stratgroup==relv_groups[gr]]))	{
							kill_me <- unique(c(kill_me,match(rownames(this_case)[rg],rownames(pbdb_rocks))));
							} else	{
							kill_me <- unique(c(kill_me,match(rownames(this_case)[gr],rownames(pbdb_rocks))));
							}
						}
					}
				} # end looking at whether groups overlap & which is latest opinion
			# add something that lets disjunct ages & areas indicate homonyms;
			} # end case with 2+ groups
		} # end case of duplicate names;
	}
pbdb_rocks <- pbdb_rocks[!(1:nrocks) %in% kill_me,];
nrocks <- nrow(pbdb_rocks);
geoplate_list <- site_list <- list();
print("Creating lists of localities and geoplates for each rock unit");
for (nr in 1:nrocks)	{
	site_list <- rlist::list.append(site_list,paleodb_collections$collection_no[paleodb_collections$full_name %in% pbdb_rocks$full_name[nr]]);
	geoplate_list <- rlist::list.append(geoplate_list,unique(paleodb_collections$geoplate[paleodb_collections$full_name %in% pbdb_rocks$full_name[nr]]));
	}
names(site_list) <- names(geoplate_list) <- pbdb_rocks$full_name;
# "formation_clean_basic"	"formation_clean_no_rock"	"formation_clean_no_rock_formal"

print("Cleaning formation names.");
pbdb_rocks$formation_clean_no_rock_formal <- pbdb_rocks$formation_clean_no_rock <- pbdb_rocks$formation_clean_basic <- pbdb_rocks$formation;
named_rock_unit <- pbdb_rocks$formation_clean_basic[pbdb_rocks$formation_clean_basic!=""];
#mundify_rock_unit_names(pbdb_rocks$formation[rock_no])
pbdb_rocks$formation_clean_basic[pbdb_rocks$formation!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$formation_clean_no_rock[pbdb_rocks$formation!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$formation_clean_no_rock_formal[pbdb_rocks$formation!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

# "member_clean_basic"      "member_clean_no_rock"      "member_clean_no_rock_formal"
print("Cleaning member names.");
pbdb_rocks$member_clean_no_rock_formal <- pbdb_rocks$member_clean_no_rock <- pbdb_rocks$member_clean_basic <- pbdb_rocks$member;
named_rock_unit <- pbdb_rocks$member_clean_basic[pbdb_rocks$member_clean_basic!=""];
pbdb_rocks$member_clean_basic[pbdb_rocks$member!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$member_clean_no_rock[pbdb_rocks$member!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$member_clean_no_rock_formal[pbdb_rocks$member!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);
pbdb_rocks$member_clean_basic[pbdb_rocks$member_clean_basic=="\\(\\)"] <- pbdb_rocks$member_clean_no_rock[pbdb_rocks$member_clean_no_rock=="()"] <- pbdb_rocks$member_clean_no_rock[pbdb_rocks$member_clean_no_rock_formal=="()"] <- "";

# "group_clean_basic"       "group_clean_no_rock"       "group_clean_no_rock_formal"
print("Cleaning group names.");
pbdb_rocks$group_clean_basic <- pbdb_rocks$group_clean_no_rock <- pbdb_rocks$group_clean_no_rock_formal <- pbdb_rocks$group;
named_rock_unit <- pbdb_rocks$group_clean_basic[pbdb_rocks$group_clean_basic!=""];
pbdb_rocks$group_clean_basic[pbdb_rocks$group!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$group_clean_no_rock[pbdb_rocks$group!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
pbdb_rocks$group_clean_no_rock_formal[pbdb_rocks$group!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

# "rock_unit_clean_basic"   "rock_unit_clean_no_rock"   "rock_unit_clean_no_rock_formal"
pbdb_rocks$rock_unit_clean_basic <- pbdb_rocks$formation_clean_basic;
fff <- (1:nrow(pbdb_rocks))[pbdb_rocks$formation_clean_basic!=""];
mmm <- (1:nrow(pbdb_rocks))[pbdb_rocks$member_clean_basic!=""];
ggg <- (1:nrow(pbdb_rocks))[pbdb_rocks$group_clean_basic!=""];
pbdb_rocks$rock_unit_clean_basic[mmm[mmm %in% fff]] <- paste(pbdb_rocks$formation_clean_basic[mmm[mmm %in% fff]]," (",pbdb_rocks$member_clean_basic[mmm[mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_basic[mmm[!mmm %in% fff]] <- paste("(",pbdb_rocks$member_clean_basic[mmm[!mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_basic[ggg[!ggg %in% c(fff,mmm)]] <- paste("[",pbdb_rocks$group_clean_basic[ggg[!ggg %in% c(fff,mmm)]],"]",sep="");
#pbdb_rocks$rock_unit_clean_basic <- pbdb_rocks$rock_unit_clean_no_rock <- pbdb_rocks$rock_unit_clean_no_rock_formal <- pbdb_rocks$full_name;
#named_rock_unit <- pbdb_rocks$rock_unit_clean_basic[pbdb_rocks$rock_unit_clean_basic!=""];
#pbdb_rocks$rock_unit_clean_basic[pbdb_rocks$full_name!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=F,delete_informal=F);
pbdb_rocks$rock_unit_clean_no_rock <- pbdb_rocks$formation_clean_no_rock;
fff <- (1:nrow(pbdb_rocks))[pbdb_rocks$formation_clean_no_rock!=""];
mmm <- (1:nrow(pbdb_rocks))[pbdb_rocks$member_clean_no_rock!=""];
ggg <- (1:nrow(pbdb_rocks))[pbdb_rocks$group_clean_no_rock!=""];
pbdb_rocks$rock_unit_clean_no_rock[mmm[mmm %in% fff]] <- paste(pbdb_rocks$formation_clean_no_rock[mmm[mmm %in% fff]]," (",pbdb_rocks$member_clean_no_rock[mmm[mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_no_rock[mmm[!mmm %in% fff]] <- paste("(",pbdb_rocks$member_clean_no_rock[mmm[!mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_no_rock[ggg[!ggg %in% c(fff,mmm)]] <- paste("[",pbdb_rocks$group_clean_no_rock[ggg[!ggg %in% c(fff,mmm)]],"]",sep="");

pbdb_rocks$rock_unit_clean_no_rock_formal <- pbdb_rocks$formation_clean_no_rock_formal;
fff <- (1:nrow(pbdb_rocks))[pbdb_rocks$formation_clean_no_rock_formal!=""];
mmm <- (1:nrow(pbdb_rocks))[pbdb_rocks$member_clean_no_rock_formal!=""];
ggg <- (1:nrow(pbdb_rocks))[pbdb_rocks$group_clean_no_rock_formal!=""];
pbdb_rocks$rock_unit_clean_no_rock_formal[mmm[mmm %in% fff]] <- paste(pbdb_rocks$formation_clean_no_rock_formal[mmm[mmm %in% fff]]," (",pbdb_rocks$member_clean_no_rock_formal[mmm[mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_no_rock_formal[mmm[!mmm %in% fff]] <- paste("(",pbdb_rocks$member_clean_no_rock_formal[mmm[!mmm %in% fff]],")",sep="");
pbdb_rocks$rock_unit_clean_no_rock_formal[ggg[!ggg %in% c(fff,mmm)]] <- paste("[",pbdb_rocks$group_clean_no_rock_formal[ggg[!ggg %in% c(fff,mmm)]],"]",sep="");
#pbdb_rocks$rock_unit_clean_no_rock_formal[pbdb_rocks$full_name!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
#pbdb_rocks$rock_unit_clean_no_rock[pbdb_rocks$full_name!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
#pbdb_rocks$rock_unit_clean_no_rock[pbdb_rocks$full_name!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
#pbdb_rocks$rock_unit_clean_no_rock_formal[pbdb_rocks$full_name!=""] <- pbapply::pbsapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=T);

pbdb_rocks$rock_unit_clean_no_rock <- gsub(" \\(\\)","",pbdb_rocks$rock_unit_clean_no_rock);
pbdb_rocks$rock_unit_clean_basic <- gsub(" \\(\\)","",pbdb_rocks$rock_unit_clean_basic);
pbdb_rocks$rock_unit_clean_no_rock <- gsub(" \\(\\)","",pbdb_rocks$rock_unit_clean_no_rock);
pbdb_rocks$rock_unit_clean_no_rock_formal <- gsub(" \\(\\)","",pbdb_rocks$rock_unit_clean_no_rock_formal);
pbdb_rocks$rock_unit_clean_basic <- gsub("\\(\\)","",pbdb_rocks$rock_unit_clean_basic);
pbdb_rocks$rock_unit_clean_no_rock <- gsub("\\(\\)","",pbdb_rocks$rock_unit_clean_no_rock);
pbdb_rocks$rock_unit_clean_no_rock_formal <- gsub("\\(\\)","",pbdb_rocks$rock_unit_clean_no_rock_formal);

## Use Geoplates to split possible homonyms.
if (geosplit)	{
nrocks <- nrow(pbdb_rocks);
def_hom <- poss_hom <- c();
max_ma <- min_ma <- length(nrocks);
for (nr in 1:nrocks)	{
	rock_sites <- subset(paleodb_collections,paleodb_collections$full_name==pbdb_rocks$full_name[nr]);
	if (nrow(rock_sites)==0 && !is.null(rock_sites$formation_alt))
		rock_sites <- subset(paleodb_collections,(paleodb_collections$formation_alt==pbdb_rocks$formation[nr] & paleodb_collections$member_alt==pbdb_rocks$member[nr]));
	max_ma[nr] <- max(rock_sites$max_ma);
	min_ma[nr] <- min(rock_sites$min_ma);
	if (length(geoplate_list[[nr]])>1)	{
		gplates <- sort(unique(geoplate_list[[nr]]));
		conts <- unique(floor(gplates/100));
		if (length(conts)==1)	{
			lbs <- ubs <- vector(length=length(gplates));
			names(lbs) <- names(ubs) <- gplates;
			for (gp in 1:length(gplates))	{
				plate_sites <- subset(rock_sites,rock_sites$geoplate==gplates[gp])
				lbs[gp]	<- max(plate_sites$max_ma);
				ubs[gp]	<- min(plate_sites$min_ma);
				}
			site_gaps <- find_stratigraphic_range_gaps(lbs,ubs);
			site_gaps[site_gaps<max_gap] <- 0;
			site_gaps <- site_gaps[site_gaps>0];
			if (length(site_gaps)>0)	{
				poss_hom <- c(poss_hom,nr);
				}
			} else if (length(conts)>1)	{
			site_conts <- floor(rock_sites$geoplate/100);
			site_ccs <- rock_sites$cc;
			unique_site_conts <- unique(site_conts);
			unique(site_ccs);
			for (sc in 1:length(site_conts))	{
				cont_sites <- subset(rock_sites,site_conts==unique_site_conts[sc]);
				if (sc==1)	{
					max_ma[nr] <- max(cont_sites$max_ma);
					min_ma[nr] <- min(cont_sites$min_ma);
					}
				}
			def_hom <- c(def_hom,nr);
			}
		}
	}
}

print("Naming three output data.frames");
output <- list(pbdb_rocks,site_list,geoplate_list);
names(output) <- c("pbdb_rocks","site_list","geoplate_list");
return(output);
}

# rock_no <- 512
# pbdb_rocks <- organized_pbdb_rocks$pbdb_rocks; wagner_rocks <- rock_database;
# pbdb_rocks <- tolower_dataframe(pbdb_rocks); wagner_rocks <- lower_wagner_rocks;
match_organized_pbdb_rock_to_external_database <- function(rock_no,pbdb_rocks,wagner_rocks,geoplate_list)	{
pbdb_unit_fields <- colnames(pbdb_rocks);
pbdb_rock_unit_fields <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock","rock_unit_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_member_fields <- match(c("member","member_clean_basic","member_clean_no_rock","member_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_formation_fields <- match(c("formation","formation_clean_basic","formation_clean_no_rock","formation_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_group_fields <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),pbdb_unit_fields);

wag_unit_fields <- colnames(wagner_rocks);
wag_rock_unit_fields <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock","rock_unit_clean_no_rock_formal"),wag_unit_fields);
wag_member_fields <- match(c("member","member_clean_basic","member_clean_no_rock","member_clean_no_rock_formal"),wag_unit_fields);
wag_formation_fields <- match(c("formation","formation_clean_basic","formation_clean_no_rock","formation_clean_no_rock_formal"),wag_unit_fields);
wag_group_fields <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),wag_unit_fields);

wrocks <- nrow(wagner_rocks);
#print(pbdb_rocks$full_name[rock_no]);
orig_rock_info <- pbdb_rocks[rock_no,];
poss_matches <- unique(which(wagner_rocks==pbdb_rocks$full_name[rock_no],arr.ind = T)[,1]);
if (pbdb_rocks$rock_unit_clean_basic[rock_no]!="")			poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_basic[rock_no],arr.ind = T)[,1])));
if (pbdb_rocks$rock_unit_clean_no_rock[rock_no]!="")		poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_no_rock[rock_no],arr.ind = T)[,1])));
if (pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no],arr.ind = T)[,1])));
if (length(poss_matches)==0)	{
	if (gsub("\\[","",pbdb_rocks$full_name[rock_no])!=pbdb_rocks$full_name[rock_no])	{
		# this is a group
		poss_matches <- unique(which(wagner_rocks==pbdb_rocks$group[rock_no],arr.ind = T)[,1]);
		poss_matches <- poss_matches[wagner_rocks$group[poss_matches]==pbdb_rocks$group[rock_no]];
		poss_matches <- poss_matches[wagner_rocks$formation[poss_matches]=="" & wagner_rocks$member[poss_matches]==""];
		if (length(poss_matches)==0)	{
			poss_matches <- unique(which(wagner_rocks==pbdb_rocks$group[rock_no],arr.ind = T)[,1]);
			poss_matches <- poss_matches[wagner_rocks$formation[poss_matches]==pbdb_rocks$group[rock_no]];
			poss_matches <- poss_matches[wagner_rocks$member[poss_matches]==""];
			}
		} else if (gsub("\\(","",pbdb_rocks$full_name[rock_no])!=pbdb_rocks$full_name[rock_no])	{
		# this is a member that is not entered independently
		if (pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	{
		# member is just upper or limestone or upper limestone
			if (pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_no_rock_formal[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member_clean_no_rock[rock_no]!="")			poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_no_rock[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member_clean_basic[rock_no]!="")				poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_basic[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member[rock_no]!="")							poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member[rock_no],arr.ind = T)[,1])));
			if (length(poss_matches)==0)	{
				poss_matches <- unique(which(wagner_rocks==pbdb_rocks$formation_clean_no_rock_formal[rock_no],arr.ind = T)[,1]);
				poss_matches <- poss_matches[wagner_rocks$member_clean_no_rock_formal[poss_matches]!=""];
				retain <- c();
				for (pm in 1:length(poss_matches))	{
					wag_version <- unlist(strsplit(wagner_rocks$member_clean_no_rock_formal[poss_matches[pm]]," "));
					pbdb_version <- unlist(strsplit(pbdb_rocks$member_clean_no_rock_formal[rock_no]," "));
					if (sum(wag_version %in% pbdb_version)>0)	retain <- c(retain,pm)
					}
				poss_matches <- poss_matches[retain];
				}
			} else	{
			if (pbdb_rocks$formation[rock_no]!="")						poss_matches <- unique(which(wagner_rocks==pbdb_rocks$formation[rock_no],arr.ind = T)[,1]);
			if (pbdb_rocks$formation_clean_basic[rock_no]!="")			poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_basic[rock_no],arr.ind = T)[,1]));
			if (pbdb_rocks$formation_clean_no_rock[rock_no]!="")		poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_no_rock[rock_no],arr.ind = T)[,1]));
			if (pbdb_rocks$formation_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_no_rock_formal[rock_no],arr.ind = T)[,1]));
			}
		} else	{
		hide <- wagner_rocks;
		go_seek <- pbdb_rocks[rock_no,pbdb_rock_unit_fields]
		poss_matches <- unique(as.vector(unlist(sapply(go_seek,find_row_with_value,hide))));
		}
	}

if (length(poss_matches)==0 && gsub("'","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no])	{
	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$full_name[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_basic[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_no_rock[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)>0)		pbdb_rocks[rock_no,] <- gsub("''","",pbdb_rocks[rock_no,]);
	}

if (length(poss_matches)==0)	{
	# case where an apostrophe is presenting the names from being recognized.
	this_rock <- unique(c(pbdb_rocks$full_name[rock_no],pbdb_rocks$rock_unit_clean_basic[rock_no],pbdb_rocks$rock_unit_clean_no_rock[rock_no],pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]));
	this_rock <- this_rock[this_rock!=""];
	if (length(this_rock)>0)
		poss_matches <- unique(c((1:wrocks)[gsub("'","",wagner_rocks$full_name) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_basic) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock_formal) %in% this_rock]));
	}

if (length(poss_matches)==0 && pbdb_rocks$member[rock_no]!="")	{
	# case where an apostrophe is presenting the names from being recognized.
	this_rock <- unique(c(pbdb_rocks$formation[rock_no],pbdb_rocks$formation_clean_basic[rock_no],pbdb_rocks$formation_clean_no_rock[rock_no],pbdb_rocks$formation_clean_no_rock_formal[rock_no]));
	this_rock <- this_rock[this_rock!=""];
	if (length(this_rock)>0)
		poss_matches <- unique(c((1:wrocks)[gsub("'","",wagner_rocks$full_name) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_basic) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock_formal) %in% this_rock]));
	}

# find groups!
if (length(poss_matches)==0 && (pbdb_rocks$formation[rock_no]=="" & pbdb_rocks$member[rock_no]==""))	{
	# case where an apostrophe is presenting the names from being recognized.
	this_rock <- unique(c(pbdb_rocks$group[rock_no],pbdb_rocks$group_clean_basic[rock_no],pbdb_rocks$group_clean_no_rock[rock_no],pbdb_rocks$group_clean_no_rock_formal[rock_no]));
	this_rock <- this_rock[this_rock!=""];
	if (length(this_rock)>0)	{
		poss_matches <- unique(c((1:wrocks)[gsub("'","",wagner_rocks$group) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$group_clean_basic) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$group_clean_no_rock) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$group_clean_no_rock_formal) %in% this_rock]));
		poss_matches <- poss_matches[wagner_rocks$formation[poss_matches]==""];
		} else	{
		# see if the group is stored as a formation
		poss_matches <- unique(c((1:wrocks)[gsub("'","",wagner_rocks$formation) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$formation_clean_basic) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$formation_clean_no_rock) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$formation_clean_no_rock_formal) %in% this_rock]));
		poss_matches <- poss_matches[wagner_rocks$member[poss_matches]==""];
		}
	}

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1)	{
#	rock_no_sr <- unique(wagner_rocks$rock_no_sr[poss_matches]);
#	poss_matches <- poss_matches[1];
	retain <- c();
	for (i in 1:length(poss_matches))	{
		if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$full_name[poss_matches[i]]))>0)
			retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_basic[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_no_rock[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_no_rock_formal[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}

#		testing <- which(wagner_rocks$full_name[poss_matches[i]]==pbdb_rocks[rock_no,pbdb_rocks[rock_no,]!=""],arr.ind=T);
#		if (nrow(testing)>0)	retain <- c(retain,i);
	if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))==1)	{
		poss_matches <- poss_matches[retain];
		} else if (length(retain)>0)	{
		# fix case where "formal" name used for member or formation, but not formation or member.
		# e.g., Zlichov (Chynice Limestone)
		if (gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no] && pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	{
			retain2 <- retain[wagner_rocks$member[poss_matches[retain]]!=""];
			if (length(retain2)>0)	retain <- retain2;
		# fix case where Alum Shale (Upper) is in database, but Alum Shale (Upper Shale) is not
			} else if ((gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no]) && pbdb_rocks$member_clean_no_rock[rock_no]!="") {
			retain2 <- retain[wagner_rocks$member_clean_no_rock[poss_matches[retain]]==pbdb_rocks$member_clean_no_rock[rock_no]];
			if (length(retain2)>0)	retain <- retain2;
			} else if (gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no] && pbdb_rocks$member_clean_no_rock_formal[rock_no]=="") {
			retain2 <- retain[wagner_rocks$member[poss_matches[retain]]==""];
			if (length(retain2)>0)	retain <- retain2;
			} else	{
			retain2 <- retain;
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && pbdb_rocks$full_name[rock_no]!=pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no])	{
		# case where "informal" name used for formation matches a formation (member) name
		# e.g., Upper Bringewood & Bringewood (Upper)
			molecularized_full <- strsplit(pbdb_rocks$full_name[rock_no]," ")[[1]];
			molecularized_formal <- strsplit(pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]," ")[[1]];
			possible_member <- paste(paste(molecularized_full[molecularized_full %in% molecularized_formal],collapse=" "),
									 " (",molecularized_full[!molecularized_full %in% molecularized_formal],")",sep="");
			retain2 <- unique(which(wagner_rocks[poss_matches,]==possible_member,arr.ind = T)[,1]);
			if (length(retain2)>0)	retain <- retain2;
#			poss_matches <- poss_matches[retain];
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (pbdb_rocks$formation[rock_no]!="" && pbdb_rocks$member[rock_no]=="")))	{
		# PBDB formation with member that gets erased
			if (sum(wagner_rocks$member[poss_matches[retain]]!="")>0)	{
				retain2 <- retain[wagner_rocks$member[poss_matches[retain]]!=""];
				} else	{
				retain2 <- retain;
				}
			if (length(retain2)>0 && length(retain2)<length(retain))	{
				retain <- retain2;
#				poss_matches <- poss_matches[retain];
				}
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (pbdb_rocks$formation[rock_no]!="" && pbdb_rocks$group[rock_no]==""))	{
		# PBDB formation that matches a group
			connects_f <- connects_g <- vector(length=length(retain));
			for (rr in 1:length(retain))	{
				if (pbdb_rocks$formation[rock_no]=="")
					connects_g[rr] <- sum(wagner_rocks[poss_matches[retain[rr]],wag_group_fields]==pbdb_rocks[rock_no,pbdb_group_fields]);
				connects_f[rr] <- sum(wagner_rocks[poss_matches[retain[rr]],wag_formation_fields]==pbdb_rocks[rock_no,pbdb_formation_fields]);
				}
			if (max(connects_f)>max(connects_g) && sum(connects_f %in% max(connects_f))==1)	retain <- retain[match(max(connects_f),connects_f)];
			if (max(connects_g)>max(connects_f) && sum(connects_g %in% max(connects_g))==1)	retain <- retain[match(max(connects_g),connects_g)];
#			if (length(retain)>0)	poss_matches <- poss_matches[retain];
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1)	{
			# case where a "cleaner" version of a name is an exact match for another rock unit.
			# example: Moscow (Lower Windom) will match both Moscow (Lower Windom) and Moscow (Windom)
			match_fits <- vector(length=length(retain));
			for (rr in 1:length(retain))
				match_fits[rr] <- sum(!is.na(match(pbdb_rocks[rock_no,pbdb_rock_unit_fields],wagner_rocks[poss_matches[retain[rr]],wag_rock_unit_fields])));
			if (max(match_fits)>0)	retain <- retain[match_fits %in% max(match_fits)];
			}
		poss_matches <- poss_matches[retain];
		} else	{
		#wagner_rocks[poss_matches,]
		}
	}

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1)
	poss_matches <- poss_matches[wagner_rocks$geoplate[poss_matches] %in% geoplate_list[[rock_no]]]

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1 && pbdb_rocks$group!="")
	poss_matches <- poss_matches[wagner_rocks$group[poss_matches] %in% pbdb_rocks$group[[rock_no]]]

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))==1)	{
	rock_no_sr <- wagner_rocks$rock_no_sr[poss_matches[1]];
	} else	{
	rock_no_sr <- 0;
	}
# why is this here?
#if (sum(pbdb_rocks[rock_no,]!=orig_rock_info)>0)	pbdb_rocks[rock_no,] <- orig_rock_info;
return(rock_no_sr);
}

possible_matches_for_organized_pbdb_rock_to_external_database <- function(rock_no,pbdb_rocks,wagner_rocks,geoplate_list)	{
pbdb_unit_fields <- colnames(pbdb_rocks);
pbdb_rock_unit_fields <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock","rock_unit_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_member_fields <- match(c("member","member_clean_basic","member_clean_no_rock","member_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_formation_fields <- match(c("formation","formation_clean_basic","formation_clean_no_rock","formation_clean_no_rock_formal"),pbdb_unit_fields);
pbdb_group_fields <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),pbdb_unit_fields);

wag_unit_fields <- colnames(wagner_rocks);
wag_rock_unit_fields <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock","rock_unit_clean_no_rock_formal"),wag_unit_fields);
wag_member_fields <- match(c("member","member_clean_basic","member_clean_no_rock","member_clean_no_rock_formal"),wag_unit_fields);
wag_formation_fields <- match(c("formation","formation_clean_basic","formation_clean_no_rock","formation_clean_no_rock_formal"),wag_unit_fields);
wag_group_fields <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),wag_unit_fields);

wrocks <- nrow(wagner_rocks);
#print(pbdb_rocks$full_name[rock_no]);
orig_rock_info <- pbdb_rocks[rock_no,];
poss_matches <- unique(which(wagner_rocks==pbdb_rocks$full_name[rock_no],arr.ind = T)[,1]);
if (pbdb_rocks$rock_unit_clean_basic[rock_no]!="")			poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_basic[rock_no],arr.ind = T)[,1])));
if (pbdb_rocks$rock_unit_clean_no_rock[rock_no]!="")		poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_no_rock[rock_no],arr.ind = T)[,1])));
if (pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no],arr.ind = T)[,1])));
if (length(poss_matches)==0)	{
	if (gsub("\\[","",pbdb_rocks$full_name[rock_no])!=pbdb_rocks$full_name[rock_no])	{
		# this is a group
		poss_matches <- unique(which(wagner_rocks==pbdb_rocks$group[rock_no],arr.ind = T)[,1]);
		poss_matches <- poss_matches[wagner_rocks$group[poss_matches]==pbdb_rocks$group[rock_no]];
		poss_matches <- poss_matches[wagner_rocks$formation[poss_matches]=="" & wagner_rocks$member[poss_matches]==""];
		if (length(poss_matches)==0)	{
			poss_matches <- unique(which(wagner_rocks==pbdb_rocks$group[rock_no],arr.ind = T)[,1]);
			poss_matches <- poss_matches[wagner_rocks$formation[poss_matches]==pbdb_rocks$group[rock_no]];
			poss_matches <- poss_matches[wagner_rocks$member[poss_matches]==""];
			}
		} else if (gsub("\\(","",pbdb_rocks$full_name[rock_no])!=pbdb_rocks$full_name[rock_no])	{
		# this is a member that is not entered independently
		if (pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	{
		# member is just upper or limestone or upper limestone
			if (pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_no_rock_formal[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member_clean_no_rock[rock_no]!="")			poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_no_rock[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member_clean_basic[rock_no]!="")				poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member_clean_basic[rock_no],arr.ind = T)[,1])));
			if (pbdb_rocks$member[rock_no]!="")							poss_matches <- unique(c(poss_matches,unique(which(wagner_rocks==pbdb_rocks$member[rock_no],arr.ind = T)[,1])));
			if (length(poss_matches)==0)	{
				poss_matches <- unique(which(wagner_rocks==pbdb_rocks$formation_clean_no_rock_formal[rock_no],arr.ind = T)[,1]);
				poss_matches <- poss_matches[wagner_rocks$member_clean_no_rock_formal[poss_matches]!=""];
				retain <- c();
				for (pm in 1:length(poss_matches))	{
					wag_version <- unlist(strsplit(wagner_rocks$member_clean_no_rock_formal[poss_matches[pm]]," "));
					pbdb_version <- unlist(strsplit(pbdb_rocks$member_clean_no_rock_formal[rock_no]," "));
					if (sum(wag_version %in% pbdb_version)>0)	retain <- c(retain,pm)
					}
				poss_matches <- poss_matches[retain];
				}
			} else	{
			if (pbdb_rocks$formation[rock_no]!="")						poss_matches <- unique(which(wagner_rocks==pbdb_rocks$formation[rock_no],arr.ind = T)[,1]);
			if (pbdb_rocks$formation_clean_basic[rock_no]!="")			poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_basic[rock_no],arr.ind = T)[,1]));
			if (pbdb_rocks$formation_clean_no_rock[rock_no]!="")		poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_no_rock[rock_no],arr.ind = T)[,1]));
			if (pbdb_rocks$formation_clean_no_rock_formal[rock_no]!="")	poss_matches <- unique(c(poss_matches,which(wagner_rocks==pbdb_rocks$formation_clean_no_rock_formal[rock_no],arr.ind = T)[,1]));
			}
		} else	{
		hide <- wagner_rocks;
		go_seek <- pbdb_rocks[rock_no,pbdb_rock_unit_fields]
		poss_matches <- unique(as.vector(unlist(sapply(go_seek,find_row_with_value,hide))));
		}
	}

if (length(poss_matches)==0 && gsub("'","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no])	{
	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$full_name[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_basic[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_no_rock[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)==0)	poss_matches <- unique(which(wagner_rocks==gsub("'","",pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]),arr.ind = T)[,1]);
	if (length(poss_matches)>0)		pbdb_rocks[rock_no,] <- gsub("''","",pbdb_rocks[rock_no,]);
	}

if (length(poss_matches)==0)	{
	# case where an apostrophe is presenting the names from being recognized.
	this_rock <- unique(c(pbdb_rocks$full_name[rock_no],pbdb_rocks$rock_unit_clean_basic[rock_no],pbdb_rocks$rock_unit_clean_no_rock[rock_no],pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]));
	this_rock <- this_rock[this_rock!=""];
	if (length(this_rock)>0)
		poss_matches <- unique(c((1:wrocks)[gsub("'","",wagner_rocks$full_name) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_basic) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock) %in% this_rock],
								 (1:wrocks)[gsub("'","",wagner_rocks$rock_unit_clean_no_rock_formal) %in% this_rock]));
	}

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1)	{
#	rock_no_sr <- unique(wagner_rocks$rock_no_sr[poss_matches]);
#	poss_matches <- poss_matches[1];
	retain <- c();
	for (i in 1:length(poss_matches))	{
		if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$full_name[poss_matches[i]]))>0)
			retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_basic[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_no_rock[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}
	if (length(retain)==0)	{
		for (i in 1:length(poss_matches))
			if (sum(tolower(pbdb_rocks[rock_no,]) %in% tolower(wagner_rocks$rock_unit_clean_no_rock_formal[poss_matches[i]]))>0)
				retain <- c(retain,i);
		}

#		testing <- which(wagner_rocks$full_name[poss_matches[i]]==pbdb_rocks[rock_no,pbdb_rocks[rock_no,]!=""],arr.ind=T);
#		if (nrow(testing)>0)	retain <- c(retain,i);
	if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))==1)	{
		poss_matches <- poss_matches[retain];
		} else if (length(retain)>0)	{
		# fix case where "formal" name used for member or formation, but not formation or member.
		# e.g., Zlichov (Chynice Limestone)
		if (gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no] && pbdb_rocks$member_clean_no_rock_formal[rock_no]!="")	{
			retain2 <- retain[wagner_rocks$member[poss_matches[retain]]!=""];
			if (length(retain2)>0)	retain <- retain2;
		# fix case where Alum Shale (Upper) is in database, but Alum Shale (Upper Shale) is not
			} else if ((gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no]) && pbdb_rocks$member_clean_no_rock[rock_no]!="") {
			retain2 <- retain[wagner_rocks$member_clean_no_rock[poss_matches[retain]]==pbdb_rocks$member_clean_no_rock[rock_no]];
			if (length(retain2)>0)	retain <- retain2;
			} else if (gsub("\\(","",pbdb_rocks$full_name[rock_no]) != pbdb_rocks$full_name[rock_no] && pbdb_rocks$member_clean_no_rock_formal[rock_no]=="") {
			retain2 <- retain[wagner_rocks$member[poss_matches[retain]]==""];
			if (length(retain2)>0)	retain <- retain2;
			} else	{
			retain2 <- retain;
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && pbdb_rocks$full_name[rock_no]!=pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no])	{
		# case where "informal" name used for formation matches a formation (member) name
		# e.g., Upper Bringewood & Bringewood (Upper)
			molecularized_full <- strsplit(pbdb_rocks$full_name[rock_no]," ")[[1]];
			molecularized_formal <- strsplit(pbdb_rocks$rock_unit_clean_no_rock_formal[rock_no]," ")[[1]];
			possible_member <- paste(paste(molecularized_full[molecularized_full %in% molecularized_formal],collapse=" "),
									 " (",molecularized_full[!molecularized_full %in% molecularized_formal],")",sep="");
			retain2 <- unique(which(wagner_rocks[poss_matches,]==possible_member,arr.ind = T)[,1]);
			if (length(retain2)>0)	retain <- retain2;
#			poss_matches <- poss_matches[retain];
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (pbdb_rocks$formation[rock_no]!="" && pbdb_rocks$member[rock_no]=="")))	{
		# PBDB formation with member that gets erased
			if (sum(wagner_rocks$member[poss_matches[retain]]!="")>0)	{
				retain2 <- retain[wagner_rocks$member[poss_matches[retain]]!=""];
				} else	{
				retain2 <- retain;
				}
			if (length(retain2)>0 && length(retain2)<length(retain))	{
				retain <- retain2;
#				poss_matches <- poss_matches[retain];
				}
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1 && (pbdb_rocks$formation[rock_no]!="" && pbdb_rocks$group[rock_no]==""))	{
		# PBDB formation that matches a group
			connects_f <- connects_g <- vector(length=length(retain));
			for (rr in 1:length(retain))	{
				if (pbdb_rocks$formation[rock_no]=="")
					connects_g[rr] <- sum(wagner_rocks[poss_matches[retain[rr]],wag_group_fields]==pbdb_rocks[rock_no,pbdb_group_fields]);
				connects_f[rr] <- sum(wagner_rocks[poss_matches[retain[rr]],wag_formation_fields]==pbdb_rocks[rock_no,pbdb_formation_fields]);
				}
			if (max(connects_f)>max(connects_g) && sum(connects_f %in% max(connects_f))==1)	retain <- retain[match(max(connects_f),connects_f)];
			if (max(connects_g)>max(connects_f) && sum(connects_g %in% max(connects_g))==1)	retain <- retain[match(max(connects_g),connects_g)];
#			if (length(retain)>0)	poss_matches <- poss_matches[retain];
			}
		if (length(unique(wagner_rocks$rock_no_sr[poss_matches[retain]]))>1)	{
			# case where a "cleaner" version of a name is an exact match for another rock unit.
			# example: Moscow (Lower Windom) will match both Moscow (Lower Windom) and Moscow (Windom)
			match_fits <- vector(length=length(retain));
			for (rr in 1:length(retain))
				match_fits[rr] <- sum(!is.na(match(pbdb_rocks[rock_no,pbdb_rock_unit_fields],wagner_rocks[poss_matches[retain[rr]],wag_rock_unit_fields])));
			if (max(match_fits)>0)	retain <- retain[match_fits %in% max(match_fits)];
			}
		poss_matches <- poss_matches[retain];
		} else	{
		#wagner_rocks[poss_matches,]
		}
	}

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1)
	poss_matches <- poss_matches[wagner_rocks$geoplate[poss_matches] %in% geoplate_list[[rock_no]]]

if (length(unique(wagner_rocks$rock_no_sr[poss_matches]))>1 && pbdb_rocks$group!="")
	poss_matches <- poss_matches[wagner_rocks$group[poss_matches] %in% pbdb_rocks$group[[rock_no]]]

return(poss_matches);
}

# rock_no <- 1425; rock_no <- 1440;
match_possible_homonym_rocks_to_external_database <- function(rock_no,pbdb_rocks,wagner_rocks,paleodb_collections,site_list,geoplate_list)	{
poss_matches <- possible_matches_for_organized_pbdb_rock_to_external_database(rock_no,pbdb_rocks,wagner_rocks,geoplate_list);
if (length(poss_matches)>0)	{
	nsites <- nrow(paleodb_collections);
	lbw <- wagner_rocks$ma_lb[poss_matches];
	ubw <- wagner_rocks$ma_ub[poss_matches];
	ttl_comps <- length(ubw);
	for (ss in 1:length(site_list[[rock_no]]))	{
		site_row <- (1:nsites)[paleodb_collections$collection_no==site_list[[rock_no]][ss]];
		lbs <- c(lbw,paleodb_collections$max_ma[site_row]);
		ubs <- c(ubw,paleodb_collections$min_ma[site_row]);
		overlaps <- accersi_temporal_overlap_matrix(lbs,ubs);
		# overlaps[3,] gives overlap of site with possible formations
		if (max(overlaps[ttl_comps+1,1:ttl_comps])>0)	{
			best_match <- poss_matches[match(max(overlaps[ttl_comps+1,]),overlaps[ttl_comps+1,])];
			if (paleodb_collections$rock_no_sr[site_row]==0)	{
				paleodb_collections$rock_no_sr[site_row] <- wagner_rocks$rock_no[best_match];
				paleodb_collections$formation_no[site_row] <- wagner_rocks$formation_no[best_match];
				paleodb_collections$rock_unit_senior[site_row] <- wagner_rocks$full_name[best_match];
				} else	{
				paleodb_collections$rock2_no_sr[site_row] <- wagner_rocks$rock_no[best_match];
				paleodb_collections$formation2_no[site_row] <- wagner_rocks$formation_no[best_match];
				paleodb_collections$rock_unit_senior2[site_row] <- wagner_rocks$full_name[best_match];
				}
			}
		}
	}
return(paleodb_collections);
}

refine_pbdb_collections_w_stratigraphic_database <- function(paleodb_collections,pbdb_rocks,wagner_rocks,site_list,finest_chronostrat,geoplate_list)	{
nsites <- nrow(paleodb_collections);
rock_no_sr <- formation_no <- rock2_no_sr <- formation2_no <- rep(0,nsites);
rock_unit_senior <- rock_unit_senior2 <- rep("",nsites);

nrocks <- length(site_list);
for (nr in 1:nrocks)	{
	if (pbdb_rocks$rock_no_sr[nr]!=0)	{
		extn_row <- match(pbdb_rocks$rock_no_sr[nr],wagner_rocks$rock_no);
		coll_nos <- site_list[[nr]];
		coll_rows <- match(coll_nos,paleodb_collections$collection_no);

		primary_finds <- coll_rows[rock_no_sr[coll_rows]==0];
		secondary_finds <- coll_rows[rock_no_sr[coll_rows]!=0];
		rock_no_sr[primary_finds] <- rock2_no_sr[secondary_finds] <- pbdb_rocks$rock_no_sr[nr];
		formation_no[primary_finds] <- formation2_no[secondary_finds] <- wagner_rocks$formation_no[extn_row];
		rock_unit_senior[primary_finds] <- rock_unit_senior2[secondary_finds] <- wagner_rocks$full_name[extn_row];
		}
	}
paleodb_collections$rock_no_sr <- as.numeric(rock_no_sr);
paleodb_collections$formation_no <- as.numeric(formation_no);
paleodb_collections$rock_unit_senior<- as.character(rock_unit_senior); #wagner_rocks$full_name[match(rock_no_sr[rock_no_sr!=0],wagner_rocks$rock_no)];
paleodb_collections$rock2_no_sr <- as.numeric(rock2_no_sr);
paleodb_collections$formation2_no <- as.numeric(formation2_no);
paleodb_collections$rock_unit_senior2 <- as.character(rock_unit_senior2);

missing_rock_nos <- (1:nrocks)[pbdb_rocks$rock_no_sr==0];
mr <- 0;
while (mr < length(missing_rock_nos))	{
	mr <- mr+1;
	paleodb_collections <- match_possible_homonym_rocks_to_external_database(rock_no=missing_rock_nos[mr],pbdb_rocks,wagner_rocks,paleodb_collections,site_list,geoplate_list);
	}

paleodb_collections$ma_lb <- paleodb_collections$max_ma;
paleodb_collections$ma_ub <- paleodb_collections$min_ma;

# NAs have crept in here somehow: fix that. 
nn <- match(paleodb_collections$rock_no_sr[paleodb_collections$rock_no_sr!=0],wagner_rocks$rock_no)
paleodb_collections$rock_no_sr[paleodb_collections$rock_no_sr!=0][is.na(nn)]
paleodb_collections$ma_lb[paleodb_collections$rock_no_sr!=0] <- wagner_rocks$ma_lb[match(paleodb_collections$rock_no_sr[paleodb_collections$rock_no_sr!=0],wagner_rocks$rock_no)];
paleodb_collections$ma_ub[paleodb_collections$rock_no_sr!=0] <- wagner_rocks$ma_ub[match(paleodb_collections$rock_no_sr[paleodb_collections$rock_no_sr!=0],wagner_rocks$rock_no)];
# original data might refine rock data.
ncoll <- nrow(paleodb_collections);
wonky_sites <- (1:ncoll)[paleodb_collections$max_ma<paleodb_collections$ma_lb];
paleodb_collections$ma_lb[paleodb_collections$max_ma<paleodb_collections$ma_lb] <- paleodb_collections$max_ma[paleodb_collections$max_ma<paleodb_collections$ma_lb];
paleodb_collections$ma_ub[paleodb_collections$min_ma>paleodb_collections$ma_ub] <- paleodb_collections$min_ma[paleodb_collections$min_ma>paleodb_collections$ma_ub];
# if original PBDB date was completely off, then these won't overlap. Assume rock age is correct
prob_sites <- (1:nsites)[paleodb_collections$ma_lb<=paleodb_collections$ma_ub];
paleodb_collections$ma_lb[prob_sites] <- wagner_rocks$ma_lb[match(paleodb_collections$rock_no_sr[prob_sites],wagner_rocks$rock_no)];
paleodb_collections$ma_ub[prob_sites] <- wagner_rocks$ma_ub[match(paleodb_collections$rock_no_sr[prob_sites],wagner_rocks$rock_no)];

print("Rebinning collections into intervals.")
age <- paleodb_collections$ma_lb;
paleodb_collections$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
age <- paleodb_collections$ma_ub;
paleodb_collections$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);

#paleodb_collections[match(site_list[[missing_rocks[2]]],paleodb_collections$collection_no),]
return(paleodb_collections)
}

further_refine_pbdb_collections_w_zones <- function(redone_collections,zone_database,rock_to_zone_database,finest_chronostrat)	{
nsites <- nrow(redone_collections);
zone_thesaurus <- accersi_zone_thesaurus(zone_database);
zd <- nrow(zone_database);
zt <- nrow(zone_thesaurus);
rzd <- nrow(rock_to_zone_database);

relv_intervals <- sort(unique(c(unique(redone_collections$interval_lb),unique(redone_collections$interval_ub))));
ttt <- match(relv_intervals,finest_chronostrat$interval);
ttt <- ttt[!is.na(ttt)];
relv_time_scale <- finest_chronostrat[min(ttt):max(ttt),];
#relv_time_scale <- relv_time_scale[!is.na(relv_time_scale),];
relv_time_scale <- relv_time_scale[order(-abs(relv_time_scale$ma_lb),-abs(relv_time_scale$ma_ub)),];

collections_w_zones <- (1:nsites)[redone_collections$zone!=""];
if (length(collections_w_zones)>0)	{
	cwz <- length(collections_w_zones);
	redone_collections$zone <- as.character(redone_collections$zone)
	zone <- as.character(redone_collections$zone[collections_w_zones]);
	print("Cleaning PBDB Zones");
#for (cz in 1:cwz)	redone_collections$zone[collections_w_zones[cz]] <- mundus_zone(zone[cz]);
	redone_collections$zone[collections_w_zones] <- pbapply::pbsapply(zone,mundus_zone);

	zone_species <- as.character(redone_collections$zone);
	non_taxon_zone <- non_taxon_zone_label <- array("",dim=nrow(redone_collections));
	print("Creating thesaurus of possible labels for zones")
	zone_species[collections_w_zones] <- pbapply::pbsapply(zone,transmogrify_full_zone_names_to_species_names_only);
	print("Looking for zones not named for taxa");
	non_taxon_zone_info <- raster::t(pbapply::pbsapply(zone,aparecium_nontaxon_zone));
	rownames(non_taxon_zone_info) <- NULL;
	non_taxon_zone[collections_w_zones] <- non_taxon_zone_info[,1];
	non_taxon_zone_label[collections_w_zones] <- non_taxon_zone_info[,2];
	#redone_collections <- cbind(redone_collections,zone_species,non_taxon_zone,non_taxon_zone_label);
	redone_collections$zone_species <- as.character(zone_species);
	redone_collections$non_taxon_zone <- as.character(non_taxon_zone);
	redone_collections$non_taxon_zone_label <- as.character(non_taxon_zone_label);

	collections_w_zones <- (1:nsites)[redone_collections$zone!=""];
	cwz <- length(collections_w_zones);
	zones_matched <- c();

	cz <- 0;
	while (cz < cwz)	{
		### blows up  when a stage is  entered as the zone!!!
		cz <- cz+1;
		pbdbc <- collections_w_zones[cz];	#redone_collections$collection_no[paldbc]
		zone_paleodb <- mundus_zone(zone=as.character(redone_collections$zone[pbdbc]));
		if (!is.na(match(zone_paleodb,zone_database$zone)))	{
			zone_paleodb_sr <- zone_database$zone_sr[match(zone_paleodb,zone_database$zone)];
			}	else	{
			zone_paleodb_sr <- "flibberty_gibberty";
			}
		zone_spc_paleodb <- as.character(redone_collections$zone_species[pbdbc]);
		zone_nontax_paleodb <- as.character(redone_collections$non_taxon_zone[pbdbc]);
		zone_nontax_lab_paleodb <- as.character(redone_collections$non_taxon_zone_label[pbdbc]);
		# poss_matches gives possible rock matches
		poss_matches <- (1:rzd)[rock_to_zone_database$rock_no_sr %in% redone_collections$rock_no_sr[pbdbc]];
		wagner_rock_unit <- match(redone_collections$rock_no_sr[pbdbc],rock_database$rock_no);

		### if we cannot find formation (member) combination, then just look for formation
		if (length(poss_matches)==0 || (redone_collections$rock_no_sr[pbdbc]>0 && rock_database$member[wagner_rock_unit]==""))
			poss_matches <- (1:rzd)[rock_to_zone_database$formation_no %in% redone_collections$formation_no[pbdbc]];
		#rock_to_zone_database[poss_matches,]
		# now, try to find this zone in zones associated with rock-unit
		matches <- c();
	#	relevant_zones <- array("",dim=c(0,5));
	#	relevant_zones <- data.frame(relevant_zones);
		relevant_zones <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
									 ma_ub=as.numeric(),interval_lb=as.character(),
									 interval_ub=as.character(),stringsAsFactors=F);
		if (length(poss_matches)>0)	{
		# match listed zones to PaleoDB entered zones
			matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_paleodb];
			# match senior names of zones to PaleoDB entered zones
			matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb];
			matches3 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb_sr];
			# match PaleoDB zone to species name alone
			matches4 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
			# match PaleoDB zone species name only to species name alone
			# match PaleoDB zone species name only to senior species name alone
			if (length(zone_spc_paleodb)>0 && zone_spc_paleodb!="")	{
				matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_spc_paleodb];
				matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_spc_paleodb]
				} else	{
				matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
				matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_paleodb];
				}
			# match PaleoDB zone species name only to senior species name alone
			matches7 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_nontax_paleodb];
			# match PaleoDB zone species name only to senior species name alone
			if (length(zone_nontax_paleodb)>0 && zone_nontax_paleodb != "")	{
				matches8 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_nontax_paleodb];
				matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_nontax_paleodb];
				} else	{
				matches8 <- poss_matches[rock_to_zone_database$non_taxon_zone[poss_matches] %in% zone_paleodb];
				matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_paleodb];
				}
			# match PaleoDB zone name only to a zone label
			matches <- sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7,matches8,matches9)));
			if (length(matches)>0)	{
				these_zones <- data.frame(zone_sr=as.character(rock_to_zone_database$zone_sr[matches]),ma_lb=as.numeric(rock_to_zone_database$ma_lb[matches]),
										  ma_ub=as.numeric(rock_to_zone_database$ma_ub[matches]),interval_lb=as.character(rock_to_zone_database$interval_lb[matches]),
										  interval_ub=as.character(rock_to_zone_database$interval_ub[matches]),stringsAsFactors=F);
				relevant_zones <- rbind(relevant_zones,these_zones);
				}
			}

		## Look out for multiple zones: if that is the caes, then get ages from all of them.
		if (length(poss_matches)>0 && length(matches)==0)	{
	#		print(paste("couldn't find",as.character(zone_spc_paleodb)));
			multizones <- divido_zone(zone=zone_paleodb);
			mz <- length(multizones);
			zz <- 0;
			while (zz<mz)	{
	#		if (mz>1)	{
	#			multimatches <- c();
	#			for (zz in 1:mz)	{
					# match whole name of zones to PaleoDB entered zone
				zz <- zz+1;
				multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
				if (length(multimatches)==0)	{
					multimatches <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% multizones[zz]];
					# match species name only to PaleoDB entered zone
					if (length(multimatches)==0)	{
						multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
						# match senior synonym species name only to PaleoDB entered zone
						if (length(multimatches)==0)	{
							multimatches <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% multizones[zz]];
							if (length(multimatches)==0)	{
								multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
								if (length(multimatches)==0)	{
									multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
									if (length(multimatches)==0)	{
										multimatches <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% multizones[zz]];
										}
									}
								}
							}
						}
					}
				if (length(multimatches)>0)	{
					these_zones <- data.frame(zone_sr=as.character(rock_to_zone_database$zone_sr[multimatches]),ma_lb=as.numeric(rock_to_zone_database$ma_lb[multimatches]),
											  ma_ub=as.numeric(rock_to_zone_database$ma_ub[multimatches]),interval_lb=as.character(rock_to_zone_database$interval_lb[multimatches]),
											  interval_ub=as.character(rock_to_zone_database$interval_ub[multimatches]),stringsAsFactors=F);
					relevant_zones <- rbind(relevant_zones,these_zones);
					} else if (length(multimatches)==0)	{
	#				mm1 <- (1:zt)[zone_thesaurus$zone %in% multizones[zz]];
	#				mm2 <- (1:zt)[zone_thesaurus$zone_species %in% multizones[zz]];
	#				mm3 <- (1:zt)[zone_thesaurus$zone_species_sr %in% multizones[zz]];
	#				mm4 <- (1:zt)[zone_thesaurus$non_taxon_zone %in% multizones[zz]];
	#				mm5 <- (1:zt)[zone_thesaurus$non_taxon_zone_label %in% multizones[zz]];
	#				mm6 <- (1:zt)[zone_thesaurus$non_taxon_zone_sr %in% multizones[zz]];
	#				mm7 <- (1:zt)[zone_thesaurus$non_taxon_zone_label_sr %in% multizones[zz]];
					multimatches <- sort(unique(which(zone_thesaurus==multizones[zz],arr.ind = T)[,1]));
	#				multimatches <- (1:zt)[zone_thesaurus$zone_sr[zone_thesaurus$zone_species %in% multizones[zz]] %in% rock_to_zone_data$zone[poss_matches]];
	#				multimatches <- sort(unique(c(mm1,mm2,mm3,mm4,mm5,mm6,mm7)));
					if (length(multimatches)>0)	{
						zone_taxa <- sort(unique(zone_thesaurus$zone_sr[multimatches]));
						poss_zones_all <- (1:zd)[zone_database$zone_sr %in% zone_taxa];
						rel_z <- cbind(as.character(zone_database$zone_sr[poss_zones_all]),as.numeric(zone_database$ma_lb[poss_zones_all]),as.numeric(zone_database$ma_ub[poss_zones_all]),as.character(zone_database$interval_lb[poss_zones_all]),as.character(zone_database$interval_ub[poss_zones_all]));
						for (ztx  in 1:length(zone_taxa))	{
							poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
							zlb <- max(as.numeric(poss_zones[,2]));
							zub <- min(as.numeric(poss_zones[,3]));
							if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])	{
								##HERE!!
	#							relevant_zones <- rbind(relevant_zones,poss_zones);
								these_zones <- data.frame(zone_sr=as.character(zone_database$zone_sr[multimatches]),ma_lb=as.numeric(zone_database$ma_lb[multimatches]),
														  ma_ub=as.numeric(zone_database$ma_ub[multimatches]),interval_lb=as.character(zone_database$interval_lb[multimatches]),
														  interval_ub=as.character(zone_database$interval_ub[multimatches]),stringsAsFactors=F);
								relevant_zones <- rbind(relevant_zones,these_zones);
								}
							}
						}
					}
				if (length(multimatches)>0)	matches <- c(matches,multimatches);
				}
			}
		### routine to find if a senior synonym of the reported zone is in the database
	#		if (length(matches)==0 && length(strsplit(zone_paleodb," ")[[1]])==1)	{
		## Look out for zones that are not "Amorphognathus ordovicicus" but "Zone D"
		if (length(poss_matches)>0 && length(matches)==0)	{
			aa <- (1:zd)[zone_database$zone %in% zone_paleodb];
			if (length(aa)==0)	{
				bb <- (1:zd)[zone_database$zone_species %in% zone_paleodb];
				} else	{
				bb <- c();
				}
			cc <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_paleodb];
			dd <- (1:zd)[zone_database$non_taxon_zone_label_sr %in% zone_paleodb];
			ee <- sort(unique(c(aa,bb,cc,dd)));
			if (length(ee)>0)	{
				try_this_zone <- unique(zone_database$zone_sr[ee]);
				try_this_zone_spc <- transmogrify_full_zone_names_to_species_names_only(try_this_zone);
				## now, repeat the matching game!
				matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% try_this_zone];
				# match senior names of zones to PaleoDB entered zones
				matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% try_this_zone];
				matches3 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% try_this_zone_spc];
				matches <- c(matches,sort(unique(c(matches1,matches2,matches3))));
				if (length(matches)==0)	{
					these_zones <- data.frame(zone_sr=as.character(zone_database$zone_sr[ee]),ma_lb=as.numeric(zone_database$ma_lb[ee]),
											  ma_ub=as.numeric(zone_database$ma_ub[ee]),interval_lb=as.character(zone_database$interval_lb[ee]),
											  interval_ub=as.character(zone_database$interval_ub[ee]),stringsAsFactors=F);
					relevant_zones <- rbind(relevant_zones,these_zones);
	#				relevant_zones <- rbind(relevant_zones,cbind(zone_database$zone_sr[ee],zone_database$ma_lb[ee],zone_database$ma_ub[ee],as.character(zone_database$interval_lb[ee]),as.character(zone_database$interval_ub[ee])));
					matches <- length(zone_database$zone_sr[ee])
					}
				}
			}

		if (length(poss_matches)==0)	{
			multizones <- divido_zone(zone=zone_paleodb);
			if (length(multizones)>1)
				zone_paleodb <- multizones;
			matches1 <- (1:zd)[zone_database$zone %in% zone_paleodb];
			if (length(matches1)==0)	{
				matches2 <- (1:zd)[zone_database$zone_species %in% zone_spc_paleodb];
				matches2 <- matches2[!zone_database$zone_species[matches2]==""];
				} else	{
				matches2 <- c();
				}
			if (zone_nontax_paleodb!="")	{
				matches3 <- (1:zd)[zone_database$non_taxon_zone %in% zone_nontax_paleodb];
				} else	{
				matches3 <- c();
				}
			if (zone_nontax_lab_paleodb!="")	{
				matches4 <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_nontax_lab_paleodb];
				} else	{
				matches4 <- c();
				}
			matches <- sort(unique(c(matches1,matches2,matches3,matches4)));
			if (length(matches)==0)	{
				multizones <- divido_zone(zone_paleodb);
				mz <- length(multizones);
				if (mz>1)	{
					rel_z <- c();
		#			multimatches <- c();
					for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
						matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
						matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
					# match species names only of zones to PaleoDB entered zones
						matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
						matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
						matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
						matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
						matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
						multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
						if (length(multimatches)>0)	{
							rel_z <- rbind(rel_z,cbind(as.character(zone_database$zone_sr[multimatches]),as.numeric(zone_database$ma_lb[multimatches]),as.numeric(zone_database$ma_ub[multimatches]),as.character(zone_database$interval_lb[multimatches]),as.character(zone_database$interval_ub[multimatches])));
							matches <- c(matches,multimatches);
							}
						}
				} else	{
					rel_z <- c();
					}
				} else	{
				rel_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[matches]),ma_lb=as.numeric(zone_database$ma_lb[matches]),ma_ub=as.numeric(zone_database$ma_ub[matches]),interval_lb=as.character(zone_database$interval_lb[matches]),interval_ub=as.character(zone_database$interval_ub[matches]),stringsAsFactors = F);
				}
			if (length(matches)>0)	{
				zone_taxa <- sort(unique(rel_z$zone_sr));
			# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
			#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
				if (length(zone_taxa)>=1)	{
	#			possible_zones <- c();
					for (ztx  in 1:length(zone_taxa))	{
						poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
						zlb <- max(as.numeric(poss_zones[,2]));
						zub <- min(as.numeric(poss_zones[,3]));
						if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])
							relevant_zones <- rbind(relevant_zones,poss_zones);
						}
					}
				}
			}	# end case where no rock unit is reported

		# now redate the collection based on zone or zones
		if (length(relevant_zones$zone_sr)==0)	{
			match1 <- (1:zd)[zone_database$zone %in% redone_collections$zone[pbdbc]];
			match2 <- (1:zd)[zone_database$zone_sr %in% redone_collections$zone[pbdbc]];
			if (redone_collections$zone_species[pbdbc]!="")	{
				match3 <- (1:zd)[zone_database$zone_species %in% redone_collections$zone_species[pbdbc]];
				match4 <- (1:zd)[zone_database$zone_species_sr %in% redone_collections$zone_species[pbdbc]];
				}	else	{
				match3 <- c();
				match4 <- c();
	#			match3 <- match3[!zone_data$zone_species[match3] %in% ""];
				}
			matches <- sort(unique(c(match1,match2,match3,match4)));
			if (length(matches)==0)	{
				multizones <- divido_zone(zone_paleodb);
				mz <- length(multizones);
				if (mz>1)	{
					rel_z <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
										ma_ub=as.numeric(),interval_lb=as.character(),
										interval_ub=as.character(),stringsAsFactors = F);
	#			multimatches <- c();
					for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
						matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
						# match senior names of zones to PaleoDB entered zones
						matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
						# match species names only of zones to PaleoDB entered zones
						matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
						matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
						matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
						matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
						matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
						multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
						if (length(multimatches)>0)	{
							nwr_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[multimatches]),ma_lb=as.numeric(zone_database$ma_lb[multimatches]),
												ma_ub=as.numeric(zone_database$ma_ub[multimatches]),interval_lb=as.character(zone_database$interval_lb[multimatches]),
												interval_ub=as.character(zone_database$interval_ub[multimatches]),stringsAsFactors = F);

							rel_z <- rbind(rel_z,nwr_z);
							matches <- c(matches,multimatches);
							}
						}
					} else	{
					rel_z <- data.frame(zone_sr=as.character(),ma_lb=as.numeric(),
										ma_ub=as.numeric(),interval_lb=as.character(),
										interval_ub=as.character(),stringsAsFactors = F);
					}
  } else	{
	rel_z <- data.frame(zone_sr=as.character(zone_database$zone_sr[matches]),
	          ma_lb=as.numeric(zone_database$ma_lb[matches]),
	          ma_ub=as.numeric(zone_database$ma_ub[matches]),
	          interval_lb=as.character(zone_database$interval_lb[matches]),
	          interval_ub=as.character(zone_database$interval_ub[matches]),
	          stringsAsFactors = F);
	rel_z <- unique(rel_z);
	#			rel_z <- cbind(as.character(zone_database$zone_sr[matches]),as.numeric(zone_database$ma_lb[matches]),as.numeric(zone_database$ma_ub[matches]),as.character(zone_database$interval_lb[matches]),as.character(zone_database$interval_ub[matches]));
	}
	if (length(matches)>0)	{
	zone_taxa <- sort(unique(rel_z$zone_sr));
	# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
	#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
	if (length(zone_taxa)>=1)	{
	#			possible_zones <- c();
	  for (ztx  in 1:length(zone_taxa))	{
	    poss_zones <- subset(rel_z,rel_z$zone_sr==zone_taxa[ztx]);
	    zlb <- max(as.numeric(poss_zones$ma_lb));
	    zub <- min(as.numeric(poss_zones$ma_ub));
	    if (zlb >= redone_collections$min_ma[pbdbc] && zub <= redone_collections$max_ma[pbdbc])
	      relevant_zones <- rbind(relevant_zones,poss_zones);
	    }
	  }
	}
	}

	# now, make adjustments
	#	colnames(relevant_zones) <- c("zone_sr","ma_lb","ma_ub","interval_lb","interval_ub");
	if (nrow(relevant_zones)>0)	{
	zones_matched <- c(zones_matched,pbdbc);
	# if zone starts after earliest possible age for formation, adjust oldest age & stage
	ma_lb <- as.numeric(as.character(relevant_zones$ma_lb));
	ma_mx <- max(ma_lb);
	ma_ub <- as.numeric(as.character(relevant_zones$ma_ub));
	ma_mn <- min(ma_ub);
	#		lbz <- match(ma_mx,ma_lb);
	#		accersi_temporal_overlap(lb1=100,lb2=50,ub1=75,ub2=40)
	new_dates <- accersi_temporal_overlap(lb1=as.numeric(redone_collections$ma_lb[pbdbc]),lb2=ma_mx,ub1=as.numeric(redone_collections$ma_ub[pbdbc]),ub2=ma_mn);
	if (new_dates[1]!=0)	{
	redone_collections$ma_lb[pbdbc] <- as.numeric(new_dates[1]);
	redone_collections$ma_ub[pbdbc] <- as.numeric(new_dates[2]);
	redone_collections$interval_lb[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>=round(as.numeric(new_dates[1]),2))]);
	redone_collections$interval_ub[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>round(as.numeric(new_dates[2]),2))]);
	} else	{
	if (as.numeric(redone_collections$ref_pubyr[pbdbc])>=1995)	{
	  redone_collections$ma_lb[pbdbc] <- ma_mx;
	  redone_collections$ma_ub[pbdbc] <- ma_mn;
	  redone_collections$interval_lb[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>=round(ma_mx,2))]);
	  redone_collections$interval_ub[pbdbc] <- as.character(relv_time_scale$interval[sum(round(relv_time_scale$ma_lb,2)>round(ma_mn,2))]);
	  }
	}
	#		if (redone_collections$ma_lb[paldbc] > ma_mx || redone_collections$ma_lb[paldbc]==0)	{
	#			redone_collections$ma_lb[paldbc] <- ma_mx;
	#			redone_collections$interval_lb[paldbc] <- as.character(relevant_zones$interval_lb[lbz]);
	#			}
	# if zone ends before latest possible age for formation, adjust youngest age & stage
	#		ma_ub <- as.numeric(as.character(relevant_zones$ma_ub));
	#		ma_mn <- min(ma_ub);
	#		ubz <- match(ma_mn,ma_ub);
	#		if (redone_collections$ma_ub[paldbc] < ma_mn)	{
	#			redone_collections$ma_ub[paldbc] <- ma_mn;
	#			redone_collections$interval_ub[paldbc] <- as.character(relevant_zones$interval_ub[ubz]);
	#			}
	}
	#	print(zones_matched);
	}

	print("Rebinning collections into intervals.")
	age <- redone_collections$ma_lb;
	redone_collections$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
	age <- redone_collections$ma_ub;
	redone_collections$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
	}
return(redone_collections);
}

# reduce external database to the time frame (give or take ma_fuzz) for the pbdb data
accersi_temporally_relevant_database <- function(database,paleodb_collections,finest_chronostrat,time_scale_fields=c("early_interval","late_interval"),ma_fuzz=25)	{
tfs <- match(time_scale_fields,colnames(paleodb_collections));
relv_slices <- unique(c(unique(paleodb_collections[,tfs[1]]),unique(paleodb_collections[,tfs[2]])));
aa <- max(abs(finest_chronostrat$ma_lb[match(relv_slices,finest_chronostrat$interval)]))+ma_fuzz;
zz <- min(abs(finest_chronostrat$ma_ub[match(relv_slices,finest_chronostrat$interval)]))-ma_fuzz;
papa_bear <- (1:nrow(database))[database$ma_lb > zz];
mama_bear <- (1:nrow(database))[database$ma_ub < aa];
baby_bear <- mama_bear[mama_bear %in% papa_bear];
return(database[baby_bear,]);
}

#paleodb_collections=clade_sites_reset
# routine_to_run_all_of_the_above
refine_pbdb_collections_w_external_databases <- function(paleodb_collections,rock_database,zone_database,rock_to_zone_database,finest_chronostrat)	{
organized_pbdb_rocks <- organize_pbdb_rock_data(paleodb_collections=paleodb_collections);

pbdb_rocks <- organized_pbdb_rocks$pbdb_rocks;
site_list <- organized_pbdb_rocks$site_list;
geoplate_list <- organized_pbdb_rocks$geoplate_list;
nrocks <- nrow(pbdb_rocks);
lower_pbdb_rocks <- tolower_dataframe(pbdb_rocks);
lower_wagner_rocks <- tolower_dataframe(rock_database);
rock_no <- 1:nrow(pbdb_rocks);
print("Matching rock units in PBDB collections to rock unit database")
rock_no_sr <- pbapply::pbsapply(rock_no,match_organized_pbdb_rock_to_external_database,pbdb_rocks=lower_pbdb_rocks,wagner_rocks=lower_wagner_rocks,geoplate_list=geoplate_list);
#rock_no_sr <- vector(length=nrocks)
#for (rr in 1:nrocks)	{
#	if (rr%%100==0)	print(rr/nrocks);
#	rock_no_sr[rr] <- match_organized_pbdb_rock_to_external_database(rock_no[rr],pbdb_rocks=lower_pbdb_rocks,wagner_rocks=lower_wagner_rocks,geoplate_list=geoplate_list);
#	}
missing_rocks <- pbdb_rocks$full_name[rock_no_sr==0];
missing_rock_nos <- match(missing_rocks,pbdb_rocks$full_name);
if (is.null(pbdb_rocks$rock_no_sr))	{
	pbdb_rocks <- cbind(rock_no_sr=as.numeric(rock_no_sr),pbdb_rocks);
	} else	{
	pbdb_rocks$rock_no_sr <- rock_no_sr;
	}
print("Adding rock unit numbers to PBDB collections & redating the collections as needed");
#dummy <- rock_database[is.na(rock_database$rock_unit_senior),];
redone_collections <- refine_pbdb_collections_w_stratigraphic_database(paleodb_collections=paleodb_collections,pbdb_rocks=pbdb_rocks,wagner_rocks=rock_database,site_list,finest_chronostrat,geoplate_list);

#if (is.null(redone_collections$ma_lb))	{
#	age <- redone_collections$ma_lb;
#	age <- redone_collections$ma_ub;
#	}
if (is.null(redone_collections$interval_lb))	{
	print("Adding new interval names to collections");
	redone_collections$interval_lb <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
	redone_collections$interval_ub <- pbapply::pbsapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
	}
print("Further redating PBDB collections given zone data");
refined_collections <- further_refine_pbdb_collections_w_zones(redone_collections,zone_database,rock_to_zone_database,finest_chronostrat);

output <- list(refined_collections,missing_rocks,pbdb_rocks,site_list);
names(output) <- c("refined_collections","missing_rocks","PBDB_rocks","PBDB_rock_site_list");
return(output)
}
