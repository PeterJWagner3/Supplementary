# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundo: clean
# percursant: scour
# revelare: reveal
source('~/Documents/R_Projects/Common_R_Source_Files/Wagner_kluges.r')  #
library(paleobioDB);	#install.packages("paleobioDB", dependencies=TRUE)
library(magrittr);		#install.packages("magrittr", dependencies=TRUE)
library(stringr);		#install.packages("stringr", dependencies=TRUE)

bad_taxa <- c("nomen dubium","nomen vanum","nomen nudum","nomen oblitum");
juniors <- c("subjective synonym of","objective synonym of","replaced by","invalid subgroup of");
#basic_rock_types <- c("dolostone","dolomite","limestone","sandstone","shale");
#not_real_rock_names <- tolower(c("basal","lowermost","lower","middle","upper","uppermost",0:100,as.character(as.roman(1:99)),basic_rock_types));
differences <- c("corrected to","invalid subgroup of","misspelling of, corrected to","misspelling of, invalid subgroup of","misspelling of, nomen dubium","misspelling of, nomen nudum","misspelling of, nomen vanum","misspelling of, objective synonym of","misspelling of, obsolete variant of","misspelling of, recombined as","misspelling of, replaced by","misspelling of, species not entered","misspelling of, subjective synonym of","misspelling of, subjective synonym of, species not entered","nomen dubium","nomen nudum","nomen oblitum","nomen vanum","objective synonym of","objective synonym of, species not entered","obsolete variant of","recombined as","replaced by","replaced by, species not entered","species not entered","subjective synonym of","subjective synonym of, species not entered");
unentered <- c("misspelling of, species not entered","misspelling of, subjective synonym of, species not entered","objective synonym of, species not entered","replaced by, species not entered","species not entered","subjective synonym of, species not entered");

hell_no <- F;

get_general_environments <- function()	{
# function to list particular environments that are either marine or terrestrial
marine <- c("marine indet.","Carbonate marine","carbonate indet.","peritidal","shallow subtidal indet.","open shallow subtidal","lagoonal/restricted shallow subtidal","sand shoal","reef, buildup or bioherm","perireef or subreef","intrashelf/intraplatform reef","platform/shelf-margin reef","slope/ramp reef","basin reef","deep subtidal ramp","deep subtidal shelf","deep subtidal indet.","offshore ramp","offshore shelf","offshore indet.","slope","basinal (carbonate)","basinal (siliceous)","Siliciclastic marine","marginal marine indet.","coastal indet.","estuary/bay","lagoonal","paralic indet.","delta plain","interdistributary bay","delta front","prodelta","deltaic indet.","foreshore","shoreface","transition zone/lower shoreface","offshore","coastal indet.","submarine fan","basinal (siliciclastic)","basinal (siliceous)","basinal (carbonate)","deep-water indet.")
terrestrial <- c("terrestrial indet.","fluvial indet.","alluvial fan","channel lag","coarse channel fill","fine channel fill","\"channel\"","channel","wet floodplain","dry floodplain","\"floodplain\"","floodplain","crevasse splay","levee","mire/swamp","fluvial-lacustrine indet.","delta plain","fluvial-deltaic indet.","lacustrine - large","lacustrine - small","pond","crater lake","lacustrine delta plain","lacustrine interdistributary bay","lacustrine delta front","lacustrine prodelta","lacustrine deltaic indet.","lacustrine indet.")
output <- list(marine,terrestrial)
names(output) <- c("marine","terrestrial")
return(output)
}

#library(BioPhysConnectoR)
# data curation scripts

# Taxonomic Curation ####
# Routine to curate species occurrences in the PaleoDB
# taxon: name of taxon
# oldest: e.g., Cambrian
# end: e.g., Permian
# start_date: oldest entered records (1998-01-01)
# end_date: most recent date (defaults to today)
# file_format: .csv gives csv file; everything else is tab-delimited
# output_files: if TRUE, then files are output; otherwise, only a list of data.frame tables is output;
# taxa <- c("Proetidae","Phillipsiidae");
# Routine to find taxa last assigned to a higher taxon that now is defunct
accersi_records_of_unentered_species <- function(species_finds)	{
# function to separate records of species that do not have authority data
types_of_diffs <- as.character(unique(species_finds$difference))
tods <- length(types_of_diffs)
unentered <- c()
for (t in 1:tods)	{
	issues <- simplify2array(strsplit(types_of_diffs[t],", "))
	if (!is.na(match("species not entered",issues)))
		unentered <- c(unentered,t)
	}
if (length(unentered)>0)	{
	unentered_species_finds <- subset(species_finds,species_finds$difference==types_of_diffs[unentered[1]])
	if (length(unentered)>1)	{
		for (u in 2:length(unentered))	{
			unentered_species_finds <- rbind(unentered_species_finds,subset(species_finds,species_finds$difference==types_of_diffs[unentered[u]]))	
			}
		}
	}
return(unentered_species_finds)
}

# get sources of species pbdb_finds: this can help get taxonomic information
accersi_sources_for_taxon_finds <- function(taxon,pbdb_finds)	{
# function to get original sources for occurrences of a particular taxon
sources <- unique(pbdb_finds$reference_no[pbdb_finds$identified_name==taxon])
output <- data.frame(identified_name=rep(taxon,length(sources)),reference_no=as.numeric(sources),citation=as.character(unique(paste(pbdb_finds$ref_author[pbdb_finds$identified_name %in% taxon],pbdb_finds$ref_pubyr[pbdb_finds$identified_name %in% taxon]))))
return(output)
}

# get species names 
accersi_species_epithets <- function(binomen)	{
# function to extract the species epithet of a binomial combination
parts <- simplify2array(strsplit(binomen," "))
return(parts[length(parts)])
}

# separate out taxa last assigned to currently invalid taxon
accersi_taxa_assigned_to_invalid_higher_taxon <- function(taxon,file_format)	{
httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxon,"&order=hierarchy&limit=all",sep="")
opinions <- read.table(httpO, sep=',', header=T)
#opinions[,13] <- clear_na_from_vector(opinions[,13],"")
#opinions[,13] <- as.numeric(opinions[,13])
taxa <- dim(opinions)[1]
statuses <- vector(length=max(opinions$orig_no))
parent_status <- vector(length=taxa)
parent <- vector(length=max(opinions$orig_no))
parent_no <- vector(length=max(opinions$orig_no))
authority <- vector(length=max(opinions$orig_no))
probs <- 0
for (t in 1:taxa) {
  tn <- as.numeric(opinions$child_spelling_no[t])
  statuses[tn] <- as.character(opinions$status[t])
  tn <- as.numeric(opinions$orig_no[t])
  statuses[tn] <- as.character(opinions$status[t])
  parent[tn] <- as.character(opinions$parent_name[t])
  parent_no[tn] <- as.numeric(opinions$parent_no[t])
  authority[tn] <- paste(as.character(opinions$author[t]),as.character(opinions$pubyr[t]),sep=" ")
  if (t>1 && opinions$status[t]=="belongs to") {
    ht <- as.numeric(opinions$parent_no[t])
    if (statuses[ht]=="belongs to") {
      parent_status[t] <- "No worries"
      } else if (statuses[ht]=="replaced by")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"has been replaced by",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="subjective synonym of")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"was subjectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="objective synonym of")  {
      pt <- as.numeric(opinions$parent_no[t])
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"was objectively synonymized with",parent[pt],"by",authority[ht],sep=" ")
      probs <- probs+1
      } else if (statuses[ht]=="nomen dubium" || (statuses[ht]=="nomen nudum" || statuses[ht]=="nomen oblitum")) {
        parent_status[t] <- paste(as.character(opinions$parent_name[t]),"is now a",statuses[ht],"according to",authority[ht],sep=" ")
        probs <- probs+1
      } else if (statuses[ht]=="invalid subgroup of") {
      parent_status[t] <- paste(as.character(opinions$parent_name[t]),"is now an invalid subgroup of",parent[pt],"according to",authority[ht],sep=" ")
      probs <- probs+1
      } else {
      parent_status[t] <- paste(as.character(opinions$parent_name[t])," is ",statuses[ht],sep="")
      probs <- probs+1
      }
    }
  }
output <- matrix(0,probs,3)
prc <- 0
for (t in 2:taxa) {
  if (opinions$status[t]=="belongs to" && parent_status[t]!="No worries"){
    prc <- prc+1
    output[prc,1] <- as.numeric(opinions$orig_no[t])
    output[prc,2] <- as.character(opinions$taxon_name[t])
    output[prc,3] <- as.character(parent_status[t])
#    print(c(as.numeric(opinions$orig_no[t]),as.character(opinions$taxon_name[t]),as.character(parent_status[t])))
    }
  }
colnames(output) <- c("orig_no","taxon_name","issue")
problem_children <- paste(taxon,"Problem_Children",sep="_")
problem_children <- paste(problem_children,file_format,sep="")

if (file_format!=".csv")  {
  write.table(output,problem_children,sep="\t",eol="\n",row.names=FALSE, col.names=TRUE)
  }
if (file_format==".csv")	{
  write.table(output,problem_children,sep=",",eol="\n",row.names=FALSE,col.names=TRUE)
  }

}

accersi_taxon_authors <- function(taxon_no)	{
# function to get authors of taxa with authority data
http <- paste("http://paleobiodb.org/data1.2/taxa/single.csv?id=txn:",taxon_no,"&show=attr",sep="")
GotURL <- RCurl::getURL(http)
taxonomy <- utils::read.csv(text = GotURL, header = TRUE, stringsAsFactors=TRUE)
return(as.character(taxonomy$taxon_attr))
}

curate_entered_species <- function(taxa,onset="Archean",end="Phanerozoic",start_date="1998-01-01",end_date="",file_format=".csv",output_files=T) {
# taxonomy curation scripts
# Arguments:
# 	taxa: proper taxonomic name
# 	onset: onset geological interval from which you want new records
# 	end: end geological interval from which you want new records
# 	start_date: the onset entered records that you want to include
# 	end_date: the most recent records that you want to include. This defaults to today!
# 	file_format: the end tag on the output files: '.xls' for Excel, '.txt', '.csv" for comma-delimited
# get temporal information for search
http_time <- "http://www.paleobiodb.org/data1.2/intervals/list.csv?scale=all";
strat_info <- read.csv(http_time, sep=',', header=T);
strat_names <- strat_info$interval_name;

# MAKE SURE WE RECOGNIZE THE OLDEST AND YOUNGEST AGES!!!
if(onset %in% strat_names) {
	start_ma <- strat_info$max_ma[match(onset, strat_names)];
	} else	{
	print("Error! ",onset," is not a term in our chronostratigraphic lexicon!")
	}

if(end %in% strat_names) {
	end_ma <- strat_info$min_ma[match(end, strat_names)];
	} else {
	print("Error! ",end," is not a term in our chronostratigraphic lexicon!")
	}

if (end_date=="")	end_date <- strsplit(as.character(Sys.time())," ")[[1]][1];
taxa <- paste(taxa,collapse=",");
http <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&occs_created_after=",start_date,"&occs_modified_before=",end_date,"&show=ref,subgenus",sep="")
raw_occurrences <- read.csv(http, header=T,stringsAsFactors = hell_no);
occurrences <- subset(raw_occurrences,raw_occurrences$identified_rank %in% c("species","subspecies"));
occurrences <- put_pbdb_dataframes_into_proper_type(occurrences);
occurrences$identified_name_original <- occurrences$identified_name;
occurrences <- expello_indeterminate_species(occurrences);
pbdb_finds <- nrow(occurrences);
# get occurrence data: if this can be modified to make it species-only, then do that.
taxon_name <- occurrences$identified_name;
occurrences$identified_name <- sapply(taxon_name,mundify_taxon_names);
occurrences <- clear_na_from_matrix(occurrences,"");

ownerless <- occurrences[occurrences$difference %in% unentered,];
#sort(unique(unentered_species_finds$difference))
## find all unentered species; look for "species not entered" as any part of the "difference" field
#for (i in 1:pbdb_finds)	{
#	wrds <- length(strsplit(as.character(occurrences$difference[i])," ")[[1]]);
#	if (wrds>3)	{
#		state <- strsplit(as.character(occurrences$difference[i])," ")[[1]];
#		if (state[wrds-2]=="species" && (state[wrds-1]=="not" && state[wrds]=="entered"))
#			occurrences$difference[i] <- as.character("species not entered");
#		}
#	}
#ownerless <- subset(occurrences,occurrences$difference=="species not entered")

dummy <- order(ownerless$identified_name);
author_needed <- ownerless[dummy,]	# occurrence information ONLY for species without author information
unentered_taxa <- unique(author_needed$identified_name);
unentered_taxa_ttl <- length(unique(author_needed$identified_name));
epithet <- unentered_taxa_finds <- unentered_taxa_main_ref <- vector(length=unentered_taxa_ttl);
for (t in 1:unentered_taxa_ttl)	{
	unentered_taxa_finds[t] <- length(subset(author_needed$identified_name,author_needed$identified_name==unentered_taxa[t]));
	epithet[t] <- divido_species_epithets(unentered_taxa[t])
	if (unentered_taxa_finds[t]==1)	{
		unentered_taxa_main_ref[t] <- as.character(subset(author_needed$primary_reference,author_needed$identified_name==unentered_taxa[t]))
		}	else {
		all_refs <- subset(author_needed$primary_reference,author_needed$identified_name==unentered_taxa[t])
		refs <- unique(all_refs)	# list of references for this species
		rr <- length(refs)		# number of references for this species
		if (rr==1)	{			# if only one, then it's easy
			unentered_taxa_main_ref[t] <- as.character(refs[1])
			} else {
			mr <- mxr <- 0	# counters to find the ref with most occurrences for this species
			for (r in 1:rr)	{
				if (length(subset(all_refs,all_refs==refs[r]))>mxr)	{
					mxr <- length(subset(all_refs,all_refs==refs[r]))
					mr <- r
					}
				}
			unentered_taxa_main_ref[t] <- as.character(refs[mr])	# tally the dominant reference
			}
		}
	}
web_text <- unentered_taxa_main_ref;
unentered_taxa_main_ref <- sapply(web_text,mundify_web_text_boring);

# put together the "needy species" for output.
#species_needy <- cbind(unentered_taxa,unentered_taxa_finds,unentered_taxa_main_ref);
species_needy <- data.frame(species=as.character(unentered_taxa),
							epithet=as.character(epithet),
							records=as.numeric(unentered_taxa_finds),
							primary_reference=as.character(unentered_taxa_main_ref),
							stringsAsFactors = hell_no);
#colnames(species_needy) <- c("Species","Records","Primary Reference");

# separate species with known authors  
#species_finds <- subset(occurrences,occurrences$identified_rank=="species");
ttl_finds <- nrow(occurrences);
# read through species.  Separate senior species from replaced ones.
#ownerless <- subset(species_finds,species_finds$difference=="species not entered");

# authored names that are supplanted
synonym_finds <- subset(occurrences,occurrences$difference %in% juniors);

# authored names that are eliminated
invalid_finds <- subset(occurrences,occurrences$difference %in% bad_taxa);

imperfect_occr_no <- sort(c(ownerless$occurrence_no,synonym_finds$occurrence_no,invalid_finds$occurrence_no));
authored_occr_no <- (1:ttl_finds)[!occurrences$occurrence_no %in% imperfect_occr_no]

valid_entered_species_finds <- occurrences[authored_occr_no,];
genera <- sort(unique(valid_entered_species_finds$genus));

author_known <- rbind(valid_entered_species_finds,synonym_finds);	# all entered species with known authorship
#author_known <- rbind(valid_entered_species_finds,synonym_finds,invalid_finds);	# all entered species with known authorship
#dalek <- (1:ncol(entered_taxa))[!colnames(entered_taxa) %in% c("accepted_name","accepted_no","identified_name","identified_no")]
#for (i in length(dalek):1)
#	author_known[,dalek[i]] <- NULL;
entered_senior <- length(unique(valid_entered_species_finds$accepted_no));
entered_junior <- entered_invalid <- length(unique(invalid_finds$identified_name));

entered_taxa <- data.frame(accepted_name=as.character(unique(as.character(valid_entered_species_finds$accepted_name))),
		accepted_no=as.numeric(unique(frak_it(valid_entered_species_finds$accepted_no))));
synonym_data <- unique(data.frame(identified_name=as.character(synonym_finds$identified_name),
								  identified_no=as.numeric(synonym_finds$identified_no),
								  stringsAsFactors = hell_no));
identified_name <- synonym_data$identified_name;
identified_no <- synonym_data$identified_no;
synonymized <- data.frame(accepted_name=as.character(identified_name),accepted_no=as.numeric(identified_no));
entered_taxa <- rbind(entered_taxa,synonymized);
entered_taxa_ttl <- nrow(entered_taxa);

### get information about constituent species
taxon_information <- accersi_taxonomic_data_for_one_taxon(taxon=taxa,inc_children=T);
taxon_information <- subset(taxon_information,taxon_information$taxon_rank %in% c("species","subspecies"));
#ww <- timestamp(prefix="",suffix = "",quiet=F);

#aa <- timestamp(prefix="",suffix = "",quiet=F);
httpO <- paste("http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=",taxa,"&rank=species,subspecies&op_type=all&show=ref,refattr",sep="");
taxon_opinions <- read.csv(httpO,header=T,stringsAsFactors = hell_no);
taxon_opinions_final <- subset(taxon_opinions,taxon_opinions$opinion_type=="class");

genera <- genera[(1:length(genera))[double_jeapordy(genera)]];
ngen <- length(genera);
dummy <- rep(0,ngen);
genus_information <- data.frame(genus=as.character(genera),total_species=as.numeric(dummy),entered_species=as.numeric(dummy),unentered_species=as.numeric(dummy),species_w_no_finds=as.numeric(dummy),stringsAsFactors = hell_no);
for (g in 1:ngen)	{
	if (length(strsplit(genera[g]," ")[[1]])==1)	{
		genus_name <- c(genera[g],paste(genera[g]," (",genera[g],")",sep=""));
		} else	{
		genus_name <- genera[g];
		}
	genus_finds_entered <- subset(valid_entered_species_finds,valid_entered_species_finds$genus %in% genus_name);
	genus_finds_unentered <- subset(ownerless,ownerless$genus %in% genus_name);
	
	genus_information$entered_species[g] <- length(unique(genus_finds_entered$accepted_name));
	genus_information$unentered_species[g] <- length(unique(genus_finds_unentered$accepted_name));
	genus_information$total_species[g] <- genus_information$entered_species[g] + genus_information$unentered_species[g];

	species_information <- accersi_taxonomic_data_for_one_taxon(taxon=genera[g],inc_children=T);
	species_information <- subset(species_information,species_information$taxon_rank %in% c("species","subspecies"));
	accepted_species <- unique(species_information$accepted_name);
	species_information <- subset(species_information,species_information$taxon_name %in% accepted_species);
	genus_information$species_w_no_finds[g] <- sum(species_information$n_occs==0);
	}

entered_taxa_finds <- entered_taxa_main_ref <- entered_taxa_last_opinion <- entered_taxa_opinions_ttl <- entered_taxa_first_opinion <- original_name <- original_author <- original_year <- vector(length=entered_taxa_ttl);

for (t in 1:entered_taxa_ttl)	{
	# find most common reference for occurrences
	species_name <- as.character(entered_taxa$accepted_name[t]);
	if (t<=entered_senior)	{
		xxx <- subset(author_known,author_known$accepted_name==species_name);
		entered_taxa_finds[t] <- length(unique(xxx$collection_no));
		} else {
		xxx <- subset(author_known,author_known$identified_name==species_name);
		entered_taxa_finds[t] <- length(unique(xxx$collection_no));
		}
	tx_refs <- length(unique(xxx$reference_no));
	if (tx_refs==1)	{
		entered_taxa_main_ref[t] <- as.character(xxx$primary_reference[1]);
		} else {
		all_refs <- as.character(unique(xxx$primary_reference));
		all_ref_nos <- unique(xxx$reference_no);
		ref_finds <- c()
		for (rr in 1:tx_refs)
			ref_finds <- c(ref_finds,sum(frak_it(xxx$reference_no)==frak_it(all_ref_nos[rr])));
		mxr <- match(max(ref_finds),ref_finds);
		entered_taxa_main_ref[t] <- as.character(unique(xxx$primary_reference))[mxr];
		}	# end case of multiple pbdb_finds
	# get taxonomic opinions
	paleodb_taxon_no <- entered_taxa$accepted_no[t];
	# get original taxa name & information
	tx <- match(paleodb_taxon_no,taxon_information$taxon_no);
	if (is.na(tx))
		which(taxon_information==paleodb_taxon_no,arr.ind = T);
	# if name has changed, then original number & current number will be different
	if (!is.na(tx))	{
		if (taxon_information$orig_no[tx]!=taxon_information$taxon_no[tx])	tx <- match(taxon_information$orig_no[tx],taxon_information$taxon_no)
		spc_opinions <- subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]);
		oldest_opinion <- match(min(spc_opinions$pubyr),spc_opinions$pubyr);
		full_name <- strsplit(as.character(species_name)," ")[[1]]
		epitaph <- full_name[length(full_name)];
		original_name[t] <- paste(as.character(spc_opinions$parent_name[oldest_opinion]),epitaph);
		original_author[t] <- paste(as.character(spc_opinions$author[oldest_opinion]),as.character(spc_opinions$pubyr[oldest_opinion]));
		original_year[t] <- frak_it(spc_opinions$pubyr[oldest_opinion]);
		if(t==1)	{
			entered_taxon_info <- taxon_information[tx,];
			colnames(entered_taxon_info) <- colnames(taxon_information);
			entered_taxon_opinions <- subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]);
			} else {
			entered_taxon_info <- rbind(entered_taxon_info,taxon_information[tx,]);
			entered_taxon_opinions <- rbind(entered_taxon_opinions,subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]));
			}	# concatenate new information onto growing entered_taxon_info matrix
		entered_taxa_opinions_ttl[t] <- nrow(subset(taxon_opinions,taxon_opinions$orig_no==taxon_information$orig_no[tx]));
		txo <- match(taxon_information$orig_no[tx],taxon_opinions_final$orig_no);
		entered_taxa_last_opinion[t] <- as.character(taxon_opinions_final$primary_reference[txo]);
		}
	}

web_text <- entered_taxa_main_ref;
entered_taxa_main_ref <- mundify_web_text_boring(web_text);
web_text <- entered_taxa_last_opinion;
entered_taxa_last_opinion <- mundify_web_text_boring(web_text);

#print(c(length(entered_taxon_info$orig_no),length(entered_taxon_info$taxon_name),length(entered_taxon_info$accepted_no),length(entered_taxon_info$accepted_name),length(entered_taxa_finds),length(entered_taxon_info$ref_author),length(entered_taxon_info$ref_pubyr),length(entered_taxa_opinions_ttl),length(entered_taxa_last_opinion)))
#species_entered <- cbind(entered_taxon_info$orig_no,entered_taxon_info$taxon_name,entered_taxon_info$accepted_no,entered_taxon_info$accepted_name,entered_taxa_finds,entered_taxon_info$ref_author,entered_taxon_info$ref_pubyr,entered_taxa_opinions_ttl,entered_taxa_last_opinion)
species_entered <- data.frame(orig_no=as.numeric(entered_taxon_info$orig_no),
							  original_name=as.character(original_name),
							  original_author=as.character(original_author),
							  original_pubyr=as.numeric(entered_taxon_info$ref_pubyr),
							  accepted_no=as.numeric(entered_taxon_info$accepted_no),
							  accepted_name=as.character(entered_taxon_info$accepted_name),
							  records=as.numeric(entered_taxa_finds),
							  total_opinions=as.numeric(entered_taxa_opinions_ttl),
							  current_opinion=as.character(entered_taxa_last_opinion),
							  stringsAsFactors = hell_no);
#colnames(species_entered) <- c("orig_no","taxon_name","accepted_no","accepted_name",
#	"records","original_authors","ref_pubyr","total_opinions","current_opinion")

# set up file name.  If just one time unit used, then use that name instead of onset and end time unit names
valid_finds <- subset(occurrences,!occurrences$occurrence_no %in% invalid_finds$occurrence_no);
# tab delimited or comma delimited
if (output_files)	{
	if (onset!=end)	{
		time_range <- paste(onset,"_to_",end,sep="")
		} else {
		time_range <- onset
		}
	filename1 <- paste(taxa,"Entered",time_range,start_date,"to",end_date,sep="_");
	filename2 <- paste(taxa,"Unentered",time_range,start_date,"to",end_date,sep="_");
	filename3 <- paste(taxa,"Entered_Opinions",time_range,start_date,"to",end_date,sep="_");
	filename4 <- paste(taxa,"Invalid_Taxa",time_range,start_date,"to",end_date,sep="_");
	filename5 <- paste(taxa,"Junior_Synonym_Finds",time_range,start_date,"to",end_date,sep="_");
	filename6 <- paste(taxa,"All_Valid_Species_Finds",time_range,start_date,"to",end_date,sep="_");
	filename7 <- paste(taxa,"Genus_Information",time_range,start_date,"to",end_date,sep="_");

	filename1 <- paste(filename1,file_format,sep="");
	filename2 <- paste(filename2,file_format,sep="");
	filename3 <- paste(filename3,file_format,sep="");
	filename4 <- paste(filename4,file_format,sep="");
	filename5 <- paste(filename5,file_format,sep="");
	filename6 <- paste(filename6,file_format,sep="");
	filename7 <- paste(filename7,file_format,sep="");
	if (file_format==".csv")	{
		write.csv(species_entered,filename1,row.names=F,fileEncoding = "UTF-8");
		write.csv(species_needy,filename2,row.names=F,fileEncoding = "UTF-8");
		write.csv(entered_taxon_opinions,filename3,row.names=F,fileEncoding = "UTF-8");
		write.csv(invalid_finds,filename4,row.names=F,fileEncoding = "UTF-8");
		write.csv(synonym_data,filename5,row.names=F,fileEncoding = "UTF-8");
		write.csv(valid_finds,filename6,row.names=F,fileEncoding = "UTF-8");
		write.csv(genus_information,filename7,row.names=F,fileEncoding = "UTF-8");
		} else if (output_files) 	{
		write.table(species_entered,filename1,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(species_needy,filename2,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(entered_taxon_opinions,filename3,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(invalid_finds,filename4,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(synonym_finds,filename5,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(valid_finds,filename6,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		write.table(genus_information,filename7,sep="\t",eol="\n",col.names=T,row.names=F,fileEncoding = "UTF-8");
		}
	}
output <- list(species_needy,species_entered,entered_taxon_opinions,synonym_finds,invalid_finds,valid_finds,genus_information);
names(output) <- c("unentered_species","entered_species","entered_opinions","synonym_finds","invalid_finds","valid_finds","genus_information");
return(output);
}

# Routine to find species with given name (e.g., "smithi") within a suprageneric taxon
find_species_name_within_higher_taxon <- function(clade,species,sub)	{
httpM <- paste("http://www.paleobiodb.org/data1.1/taxa/list.csv?name=",clade,"&rel=all_children&show=attr,app",sep="")
members <- read.table(httpM, sep=',', header=T)
member_species <- subset(members, members$rank=="species")
if (sub==TRUE)
	member_species <- rbind(member_species,subset(members, members$rank=="subspecies"))

sp_ct <- dim(member_species)[1]

matches <- 0
if (sp_ct>0){
  for (s in 1:sp_ct)	{
  	words <- length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  	species_name <- strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  	if (species_name==species)
  		matches <- matches+1
  	}

  if (matches>0)	{	
  	poss_species <- vector(length=matches)
  	m <- 0
  	for (s in 1:sp_ct)	{
  		words <- length(strsplit(as.character(member_species$taxon_name[s])," ")[[1]])
  		species_name <- strsplit(as.character(member_species$taxon_name[s])," ")[[1]][words]
  		if (species_name==species)	{
  			m <- m+1
  			poss_species[m] <- as.character(member_species$taxon_name[s])
  			}
  		}
  	} else {
  	  poss_species <- vector(length=1)
  	  poss_species[1] <- paste("no matches")
  	  
  	}
  } else {
    poss_species <- vector(length=1)
    poss_species[1] <- paste(clade," not in database")
  }

return(poss_species)
}

# clean occurrence ids
mundare_occurrence_identifications <- function(taxon_name)	{
# function to remove question marks, cf.s, affs., etc.
text <- substring(taxon_name, seq(1, nchar(taxon_name), 1), seq(1, nchar(taxon_name), 1))
ttlch <- nchar(taxon_name)
if (text[1]=="\"")	{
	taxon_name <- pracma::strcat(text[2:ttlch])
	text <- text[2:ttlch]
	}
ttlch <- nchar(taxon_name)
if (text[ttlch]=="\"")
	taxon_name <- pracma::strcat(text[1:(ttlch-1)])
taxon_name <- gsub(" n\\. sp\\.","",taxon_name)
taxon_name <- gsub("n\\. gen\\. ","",taxon_name)
taxon_name <- gsub(" n\\. subgen\\.","",taxon_name)
taxon_name <- gsub(" cf\\.","",taxon_name)
taxon_name <- gsub("cf\\. ","",taxon_name)
taxon_name <- gsub(" aff\\.","",taxon_name)
taxon_name <- gsub("aff\\. ","",taxon_name)
taxon_name <- gsub(" ex gr\\.","",taxon_name)
taxon_name <- gsub(" informal","",taxon_name)
taxon_name <- gsub(" sensu lato","",taxon_name)
taxon_name <- gsub("\" "," ",taxon_name)
taxon_name <- gsub(" \""," ",taxon_name)
taxon_name <- gsub("\" ","",taxon_name)
taxon_name <- gsub(" \\?" ,"",taxon_name)
taxon_name <- gsub("\\? " ,"",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
taxon_name <- gsub("  " ," ",taxon_name)
return(as.character(taxon_name))
}

# get information on species for curation;
paleodb_entered_species_for_curation <- function(taxa,onset="Cambrian",end="Holocene",oldest_entry,latest_entry=NA,save_files=TRUE,output_type=".txt") {
# Relevant functions
# function to get occurrence data for a group from a general interval
#	of time and entered between particular dates
# Will return lists, but also can output csv, txt, tab, etc., files directly
# Returns three basic results:
# 1. Genus species combinations with no authority data
# 2. Sources for occurrences of those species
# 3. Genus species combinations for species with authority data
timespan <- onset;
if (onset!=end)	timespan <- paste(timespan,"-",end,sep="");
if (is.na(latest_entry))	{
	now <- as.character(Sys.time())
	latest_entry <- strsplit(now," ")[[1]][1]
	}

all_finds <- accersi_occurrence_data(taxa=taxa,onset=onset,end=end,species_only=T,save_files = F);
species_finds <- all_finds[all_finds$identified_rank %in% c("species","subspecies"),];
cleaned_names <- pbapply::pbsapply(as.character(species_finds$identified_name),mundare_occurrence_identifications);
species_finds$identified_name <- cleaned_names;

# set aside pbdb_finds for species with no taxonomic information
unentered_species_finds <- accersi_records_of_unentered_species(species_finds);
unentered_species <- sort(unique(unentered_species_finds$identified_name));
u_s <- length(unentered_species);
# get the references generating the species lacking taxonomic information
#taxon <- unentered_species;
nfinds <- vector(length=u_s); names(nfinds) <- unentered_species;
print("Getting occurrence references.")
unentered_species_sources <- data.frame(identified_name=as.character(),reference_no=as.numeric(),citation=as.character());
for (us in 1:u_s)	{
	new_info <- accersi_sources_for_taxon_finds(taxon=unentered_species[us],pbdb_finds=all_finds);
	unentered_species_sources <- rbind(unentered_species_sources,new_info);
	nfinds <- sum(unentered_species_finds$identified %in% unentered_species[us]);
	}
print("Getting species epithets");
unentered_epithets <- pbapply::pbsapply(unentered_species,accersi_species_epithets)
unentered_species_info <- data.frame(identified_name=as.character(unentered_species),epithet=as.character(unentered_epithets),nfinds=as.numeric(nfinds));

#ttl_finds <- nrow(species_finds);
#xx <- (1:ttl_finds)[!species_finds$occurrence_no %in% unentered_species_finds$occurrence_no]
entered_species_finds <- species_finds[!species_finds$occurrence_no %in% unentered_species_finds$occurrence_no,];
redone <- entered_species_finds[order(entered_species_finds$accepted_name,entered_species_finds$identified_name),]
entered_species <- data.frame(accepted_no=as.numeric(redone$accepted_no),accepted_name=as.character(redone$accepted_name),identified_no=as.numeric(redone$identified_no),identified_name=as.character(redone$identified_name));
#entered_species <- 	unique(cbind(redone$accepted_no,as.character(redone$accepted_name),redone$identified_no,redone$identified_name))
#colnames(entered_species) <- c("accepted_no","accepted_name","identified_no","identified_name")
print("Getting authority information for entered species.")
entered_species$authority <- pbapply::pbsapply(as.numeric(entered_species[,3]),accersi_taxon_authors);
#unentered_species
if (save_files)	{
	output1 <- paste(timespan,"_",taxa,"_Unentered_Species",output_type,sep="")
	output2 <- paste(timespan,"_",taxa,"_Sources_for_Unentered_Species",output_type,sep="")
	output3 <- paste(timespan,"_",taxa,"_Entered_Species",output_type,sep="")
	output4 <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="")
	if (output_type==".csv")	{
		write.csv(unentered_species_info,output1,row.names = F);
		write.csv(unentered_species_sources,output2,row.names = F);
		write.csv(entered_species,output3,row.names = F);
		write.csv(species_finds,output4,row.names = F);
		} else if (output_type %in% c(".xls",".xlsx"))	{
		writexl::write_xlsx(unentered_species_info,output1);
		writexl::write_xlsx(unentered_species_sources,output2);
		writexl::write_xlsx(entered_species,output3);
		writexl::write_xlsx(species_finds,output4);
		} else	{
		write.table(unentered_species_info,file=output1,sep = "\t",row.names = F)
		write.table(unentered_species_sources,file=output2,sep = sepr,row.names = F)
		write.table(entered_species,file=output3,sep = sepr,row.names = F)
		write.table(species_finds,file=output4,sep = sepr,row.names = F)
		}
	}

output <- list(unentered_species,unentered_species_sources,entered_species)
names(output) <- c("Unentered_Species","Unentered_Species_Sources","Entered_Species")
return(output)
}

revelare_youngest_daughters <- function(taxon,earliest="Proterozoic",latest="Holocene",output_results=T,file_format="csv")	{
#http <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&taxon_reso=species&interval=",earliest,",",latest,sep="")
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",earliest,",",latest,"&envtype=terr,marine,unknown&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
occurrences <- read.csv(http, header=T,stringsAsFactors = F,fileEncoding = "UTF-8");
occurrences <- occurrences[occurrences$identified_rank %in% c("species","subspecies"),];
lb_youngest <- min(occurrences$max_ma);
ub_youngest <- min(occurrences$min_ma);
#records <- nrow(occurrences);
arwens <- occurrences[occurrences$min_ma<lb_youngest,];
if (output_results) {
	brats <- paste("Oldest_",taxon,"_species.csv",sep="")
	write.csv(arwens,brats,row.names=F,fileEncoding = "UTF-8");
	}
return(finwes);
}

# routine to find the oldest species in a taxon
#lumos_vetustissima_filiam <- function(taxon,earliest,latest)	{
revelare_oldest_daughters <- function(taxon,earliest="Proterozoic",latest="Holocene",output_results=T,file_format="csv")	{
#http <- paste("http://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&taxon_reso=species&interval=",earliest,",",latest,sep="")
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxon,"&interval=",earliest,",",latest,"&envtype=terr,marine,unknown&show=refattr,classext,rem,entname,abund,crmod&limit=all",sep = "");
occurrences <- read.csv(http, header=T,stringsAsFactors = F,fileEncoding = "UTF-8");
occurrences <- occurrences[occurrences$identified_rank %in% c("species","subspecies"),];
lb_oldest <- max(occurrences$max_ma);
ub_oldest <- max(occurrences$min_ma);
#records <- nrow(occurrences);
finwes <- occurrences[occurrences$max_ma>ub_oldest,];
if (output_results) {
	methuselahs <- paste("Oldest_",taxon,"_species.csv",sep="")
	write.csv(finwes,methuselahs,row.names=F,fileEncoding = "UTF-8");
	}
return(finwes);
}


# Taxonomic Stratigraphic Range Curation ####
# Routine to find gaps in occurrences for taxa
#revelare_foraminibus_in_petra <- function(taxa,onset,end,signif_gap,file_format=".csv")	{
accersi_collections_with_taxon <- function(taxon,pbdb_finds)	{
# function to get pbdb_sites containing a particular taxon
ttl_finds <- dim(pbdb_finds)[1]
return((1:ttl_finds)[pbdb_finds$identified_name %in% taxon])
}

# find genera assigned to a higher taxon that have no occurrences
accersi_genera_in_absentio <- function(higher_taxon,file_format="csv")	{
#higher_taxon: a higher taxon that includes genera and subgenera
#tranks: lower taxon ranks to include; write 'genus,subgenus' to get genera & subgenera
#file_format:
#	csv: comma-delimited
#	txt: tab-delimited text
#	tab: tab-delimited text
#	xls: tab-delimited text that usually will be opened directly by Excel
#tranks <- "genus,subgenus"
#http <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",higher_taxon,sep="")
#http <- paste(http,"&rank=",tranks,sep="")
#http <- paste(http,"&show=attr,parent",sep="")
if (length(higher_taxon)>1)	higher_taxon <- paste(higher_taxon,collapse=",");
http <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",higher_taxon,"&rank=genus,subgenus&show=attr,parent",sep="");
taxonomy <- read.csv(http,header=T,stringsAsFactors = F,fileEncoding = "UTF-8");
absent <- subset(taxonomy,taxonomy$n_occs==0)
invalids <- subset(absent,absent$difference!="")
valids <- subset(absent,absent$difference=="")
subgenera <- subset(taxonomy,taxonomy$taxon_rank=="subgenus")
subgenus_name <- str_replace_all(word(subgenera$taxon_name,2), "[[:punct:]]", "")
subgenera <- cbind(subgenera,subgenus_name);
cases <- nrow(valids);
ex <- 0;
i <- 1;
# look for subgenera such as Cardium (Cardium) where Cardium has occurrences
while (i<=cases)	{
	if (as.character(valids$taxon_rank[i])=="subgenus")	{
		sgn <- str_replace_all(word(valids$taxon_name[i],2), "[[:punct:]]", "")
		gen <- word(valids$taxon_name[i],1)
		gen_no <- match(sgn,taxonomy$taxon_name)
		if (is.na(gen_no)==FALSE)	{
			doh <- paste(i-1,cases-1,ex+1,as.character(valids$taxon_name[i]),sep=" ")
			if (i==1)	{
				valids <- valids[2:cases,]
				}	else if (i<cases) {
				valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
				} else valids <- valids[1:(i-1),]
			i <- i-1
			cases <- cases-1
			ex <- ex+1
			}	
		# end case where we have a subgenus
		} else {
		sgn_no <- match(valids$taxon_name[i],subgenera$subgenus_name)
		if (is.na(sgn_no)==FALSE)	{
			if (subgenera$n_occs[sgn_no]>0)	{
				if (i==1)	{
					valids <- valids[2:cases,]
					}	else if (i<cases) {
					valids <- rbind(valids[1:(i-1),],valids[(i+1):cases,])
					} else valids <- valids[1:(i-1),]
				i <- i-1
				cases <- cases-1
				ex <- ex+1
				}
			} 	
		}	# end case where we have a genus	
	i <- i+1
	}

if (cases>0)	{
	filename <- paste(higher_taxon,"_Missing_from_PBDB.",file_format,sep="");
	write.csv(valids,filename,row.names=FALSE,fileEncoding = "UTF-8");
	}	else print("?!?!  All of the genera have occurrences!!! High Five!!!!");
}

count_collections_with_taxon <- function(taxon,pbdb_finds)	{
# function to get pbdb_sites containing a particular taxon
return(length((1:dim(pbdb_finds)[1])[pbdb_finds$identified_name %in% taxon]))
}

# find genus entered as both genus AND subgenus
double_jeapordy <- function(genera)	{
legit <- rep(T,length(genera));
for (g in 1:length(genera))	{
	if (length(strsplit(genera[g]," ")[[1]])==1)	{
		automatic_subgenus <- paste(genera[g]," (",genera[g],")",sep="");
		if (!is.na(match(automatic_subgenus,genera)))	{
			legit[match(automatic_subgenus,genera)] <- F;
			}
		}
	}
return(legit);
}

reveal_stratigrapic_gaps_in_taxon_ranges <- function(taxa,onset,end,signif_gap,file_format=".csv",output_file=T)	{
# get temporal information for search
# get temporal information for search
http_time <- "http://www.paleobiodb.org/data1.2/intervals/list.csv?scale=1";
strat_info <- read.table(http_time, sep=',', header=T);
strat_names <- strat_info$interval_name;

# MAKE SURE WE RECOGNIZE THE OLDEST AND YOUNGEST AGES!!!
if(onset %in% strat_names) {
	start_ma <- strat_info$max_ma[match(onset, strat_names)];
	} else	{ print("Error! ",onset," is not a term in our chronostratigraphic lexicon!");
	}

if(end %in% strat_names) {
	end_ma <- strat_info$min_ma[match(end, strat_names)];
	} else {
	print("Error! ",end," is not a term in our chronostratigraphic lexicon!")
	}

# get genus names
httpG <- paste("http://paleobiodb.org/data1.1/taxa/list.csv?name=",taxon,"&rel=all_children&show=attr&rank=genus,subgenus&limit=9999",sep="")
genera <- read.csv(httpG, header=T,stringsAsFactors = hell_no);
httpF <- paste("http://paleobiodb.org/data1.1/occs/list.csv?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,loc,crmod&limit=all",sep="");
occurrences <- read.csv(httpF, header=T,stringsAsFactors = hell_no);
	
# remove n. sp., aff., etc.
taxon_name <- occurrences$taxon_name;
occurrences$taxon_name <- sapply(taxon_name,mundify_taxon_names);
taxon_name <- occurrences$matched_name;
occurrences$matched_name <- sapply(taxon_name,mundify_taxon_names);
#occurrences$late_age <- clear_na_from_vector(occurrences$late_age,0)
if (length(is.na(occurrences$late_age))>0)
	occurrences$late_age[is.na(occurrences$late_age)] <- 0;
occurrences$late_age <- as.numeric(occurrences$late_age);
# sometimes taxon_name comes up NA; replace it with matched_name
occurrences <- expello_na_from_matrix(occurrences);
mid_age <- (occurrences$early_age+occurrences$late_age)/2;
occurrences <- cbind(occurrences,mid_age=as.numeric(mid_age));
pbdb_finds <- length(occurrences$matched_name);
found_genus <- vector(length=pbdb_finds);

# separate genus name
for (i in 1:pbdb_finds)	found_genus[i] <- strsplit(as.character(occurrences$matched_name[i])," ")[[1]][1];

occurrences$found_genus <- as.character(found_genus);
sortdata <- c("found_genus","mid_age");
# flip-flop age so that we can get oldest members first
occurrences$mid_age <- -1*occurrences$mid_age;
occurrences <- occurrences[do.call("order",occurrences[sortdata]),];
occurrences$mid_age <- -1*occurrences$mid_age;
	
#k <- 0;
gc <- gaps <- c();
for (i in 1:(pbdb_finds-1))  {
	j <- i+1;
#	k <- k+1;
	if ((occurrences$taxon_rank[i]=="genus" || occurrences$taxon_rank[i]=="subgenus") || occurrences$taxon_rank[i]=="species")	
		if (as.character(occurrences$found_genus[i])==as.character(occurrences$found_genus[j]) && (occurrences$late_age[i]-occurrences$early_age[j])>signif_gap)	{
			gc <- rbind(gc,c(i,j));
			gaps <- c(gaps,occurrences$late_age[i]-occurrences$early_age[j]);
			}
	}
genus_gaps <- data.frame(genus=as.character(occurrences$found_genus[gc[,1]]),
						 taxon_older=as.character(occurrences$taxon_name[gc[,1]]),
						 collection_no_older=as.numeric(occurrences$collection_no[gc[,1]]),
						 early_interval_older=as.character(occurrences$early_interval[gc[,1]]),
						 late_interval_older=as.character(occurrences$late_interval[gc[,1]]),
						 early_age_older=as.numeric(occurrences$early_age[gc[,1]]),
						 late_age_older=as.numeric(occurrences$late_age[gc[,1]]),
						 taxon_younger=as.character(occurrences$taxon_name[gc[,2]]),
						 collection_no_younger=as.numeric(occurrences$collection_no[gc[,2]]),
						 early_interval_younger=as.character(occurrences$early_interval[gc[,2]]),
						 late_interval_younger=as.character(occurrences$late_interval[gc[,2]]),
						 early_age_younger=as.numeric(occurrences$early_age[gc[,2]]),
						 late_age_younger=as.numeric(occurrences$late_age[gc[,2]]),
						 gap_myr=as.numeric(gaps),stringsAsFactors = hell_no);

if (output_file)	{
	if (onset!=end)	{
		time_range <- paste(onset,"-",end,sep="")
		} else {
		time_range <- onset	
		}
	filename <- paste("Suspicious_Gaps",onset,time_range,sep="_")
	filename <- paste(filename,file_format,sep="")
	if (file_format!=".csv")	{
		write.csv(genus_gaps,filename,row.names=FALSE)
		}	else	{
		write.table(genus_gaps,filename1,sep=",",eol="\n",row.names=FALSE)
		}
	}
return(genus_gaps)
}

# get distribution of gaps within taxon ranges
revelare_minimum_gaps_among_ranges <- function(strat_ranges)	{
if (!is.data.frame(strat_ranges))
	strat_ranges <- data.frame(strat_ranges);
colnames(strat_ranges) <- c("FAD","LAD");
nfinds <- nrow(strat_ranges);
rownames(strat_ranges) <- 1:nfinds;
#gap_before <- gap_after <- array(0,dim=nfinds);
gaps <- array(0,dim=nfinds-1);
for (n in 1:(nfinds-1))
	if (strat_ranges$FAD[n+1]>strat_ranges$LAD[n])
		gaps[n] <- strat_ranges$FAD[n+1]-strat_ranges$LAD[n];
#min_gap <- sum(gap_pre) + sum(gap_pst);
#min_gap <- sum(gap_before);
return(gaps);
}

# Rock Unit Curation Routines ####
# Collect information about rock units to make stratigraphy consistent
accersi_latest_rock_unit_assignment <- function(formations,members,ref_pubyr,reference_no)	{
# function to get the latest formation/member combination for rocks with
#	multiple such combinations and/or rankings in the PaleoDB.
member_or_formation <- members[(1:length(members))[members %in% formations]]
member_or_formation <- member_or_formation[member_or_formation!=""]
ttl_coll <- length(members)
last_authors <- vector(length=ttl_coll)
if (length(member_or_formation)>0)	{
	for (c in 1:length(member_or_formation))	{
		# Use latest opinion.  If latest opinion is "member," then reassign
		#	all pbdb_sites to the latest formation/member combo
		# If the latest opinion is tied, then go with majority rule.  If that
		#	is tied, too, then just make the damned thing a formation....
		vote_formation <- (1:ttl_coll)[formations %in% member_or_formation[c]]
		vote_member <- (1:ttl_coll)[members %in% member_or_formation[c]]
		if (length(vote_member)>0 && length(vote_formation)>0)	{
			if (max(ref_pubyr[vote_formation]) > max(ref_pubyr[vote_member]) || (max(ref_pubyr[vote_formation]) == max(ref_pubyr[vote_member]) && length(vote_formation) > length(vote_member)))	{
				### elevate member to formation in appropriate pbdb_sites
				formations[vote_member] <- member_or_formation[c]
				members[vote_member] <- ""
				latest_opinion <- vote_formation[match(max(ref_pubyr[vote_formation]),ref_pubyr[vote_formation])]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				}	else if ((max(ref_pubyr[vote_formation]) < max(ref_pubyr[vote_member])) || (max(ref_pubyr[vote_formation]) == max(ref_pubyr[vote_member]) && length(vote_formation) < length(vote_member)))	{
				latest_opinion <- vote_member[match(max(ref_pubyr[vote_member]),ref_pubyr[vote_member])]
				formations[vote_formation] <- formations[latest_opinion]
				members[vote_formation] <- member_or_formation[c]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				} else if (sum(vote_formation==vote_member)==length(vote_member))	{
				members[vote_member] <- ""
				} else	{
				latest_opinion <- vote_member[match(max(ref_pubyr[vote_member]),ref_pubyr[vote_member])]
				formations[vote_formation] <- formations[latest_opinion]
				members[vote_formation] <- member_or_formation[c]
				last_authors[c(vote_formation,vote_member)] <- reference_no[latest_opinion]
				}
			}
#		c <- c+1
		}
	}
output <- cbind(formations,members)
colnames(output) <- c("Formation_Latest","Member_Latest")
return(output)
}

confused_rock_assignments <- function(formations,members,reference_no,ref_pubyr,primary_reference)	{
# function to separate out rock units that seem to have both formation and
#	member status at different localities.  This returns the different references
#	expressing or using the opinions in order for someone to decide what is
#	the best combination to use now.  
member_or_formation <- members[(1:length(members))[members %in% formations]];
member_or_formation <- member_or_formation[member_or_formation!=""];
ttl_coll <- length(members);
dazed_and <- c();
dazed_and <- data.frame(formation=as.character(),member=as.character(),ref_pubyr=as.numeric(),reference_no=as.numeric(),primary_reference=as.character());
if (length(member_or_formation)>0)	{
	confused <- sort(unique(member_or_formation))
	confused <- confused[nchar(confused)>1]
	for (c in 1:length(confused))	{
		# Use latest opinion.  If latest opinion is "member," then reassign
		#	all pbdb_sites to the latest formation/member combo
		# If the latest opinion is tied, then go with majority rule.  If that
		#	is tied, too, then just make the damned thing a formation....
		vote_formation <- (1:ttl_coll)[formations %in% confused[c]]
		vote_member <- (1:ttl_coll)[members %in% confused[c]]
		rocks <- c(vote_formation,vote_member)
		all_refs <- reference_no[rocks]
		uniq_refs <-unique(all_refs)
		relv_colls <- rocks[match(uniq_refs,all_refs)]
		dazed_and <- rbind(dazed_and,cbind(formations=as.character(formations[relv_colls]),members=as.character(members[relv_colls]),ref_pubyr=as.numeric(ref_pubyr[relv_colls]),reference_no=as.numeric(reference_no[relv_colls]),primary_reference=as.character(primary_reference[relv_colls])))
		}
	}
#colnames(dazed_and) <- c("formation","member","ref_pubyr","reference_no","primary_reference");
return(dazed_and);
}

look_up_pbdb_formations <- function(taxon,onset,end,file_format=".csv")	{
	# get temporal information for search
http_time <- "http://www.paleobiodb.org/data1.2/intervals/list.csv?scale=all";
strat_info <- read.table(http_time, sep=',', header=T);
strat_names <- strat_info$interval_name;
	
# MAKE SURE WE RECOGNIZE THE OLDEST AND YOUNGEST AGES!!!
if(onset %in% strat_names) {
	start_ma <- strat_info$max_ma[match(onset, strat_names)];
	} else	{
	print("Error! ",onset," is not a term in our chronostratigraphic lexicon!")
	}

if(end %in% strat_names) {
	end_ma <- strat_info$min_ma[match(end, strat_names)];
	} else {
	print("Error! ",end," is not a term in our chronostratigraphic lexicon!")
	}

taxa <- paste(taxa, collapse = ",");
http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr,entname",sep="")
fetch <- RCurl::getURL(http)
pbdb_sites <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=TRUE);
#	
#	http <- paste("http://paleobiodb.org/data1.1/colls/list.csv?base_name=",taxon,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
#http <- paste("http://paleobiodb.org/data1.1/colls/list.csv?",request,"&max_ma=",start_ma,"&min_ma=",end_ma,"&show=time,stratext,crmod&limit=all",sep="")
#pbdb_sites  <- read.table(http, sep=',', header=T)
pbdb_sites <- clear_na_from_matrix(pbdb_sites, "")
ncoll <- nrow(pbdb_sites);
mid_ma <- (pbdb_sites$max_ma+pbdb_sites$min_ma)/2
pbdb_sites <- cbind(pbdb_sites,mid_ma)
named_rock_unit <- as.character(pbdb_sites$formation);
pbdb_sites$formation <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=FALSE);
named_rock_unit <- as.character(pbdb_sites$member);
pbdb_sites$member <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=FALSE);

zone <- pbdb_sites$zone
pbdb_sites$zone <- sapply(zone,turgio_zone);

# remove quotation debris
for (i in 1:ncol(pbdb_sites))	{
	pbdb_sites[,i] <- gsub("â€œ","",pbdb_sites[,i])
	pbdb_sites[,i] <- gsub("â€\u009d","",pbdb_sites[,i])
	}
	
rock_unit <- pbdb_sites$formation;
w_memb <- (1:ncoll)[pbdb_sites$member!=""];
rock_unit[w_memb] <- paste(pbdb_sites$formation[w_memb]," (",pbdb_sites$member[w_memb],")",sep="");
pbdb_sites <- cbind(pbdb_sites,rock_unit);
sortdata <- c("rock_unit","mid_ma","early_interval","zone")
pbdb_sites <- pbdb_sites[do.call("order",pbdb_sites[sortdata]),];	# cool trick for sorting data!

formation_member <- unique(cbind(as.character(pbdb_sites$formation[pbdb_sites$rock_unit!=""]),as.character(pbdb_sites$member[pbdb_sites$rock_unit!=""])));
reduced_collections <- data.frame(formation=as.character(pbdb_sites$formation),member=as.character(pbdb_sites$member),rock_unit=as.character(pbdb_sites$rock_unit),early_interval=as.character(pbdb_sites$early_interval),late_interval=as.character(pbdb_sites$late_interval),zone=as.character(pbdb_sites$zone),max_ma=as.numeric(pbdb_sites$max_ma),min_ma=as.numeric(pbdb_sites$min_ma),mid_ma=as.numeric(pbdb_sites$mid_ma));
#match("Tropic",orphans$member)
orphans <- subset(pbdb_sites,pbdb_sites$rock_unit=="");
rocky_collections <- subset(reduced_collections,reduced_collections$rock_unit!="");
named_rocks <- unique(rocky_collections$rock_unit);
entered_formations <- sort(unique(rocky_collections$formation[rocky_collections$formation!=""]));
entered_members <- sort(unique(rocky_collections$member[rocky_collections$member!=""]));

confused_rocks <- entered_formations[entered_formations %in% entered_members];

ttl_ints <- other_intervals <- early_interval <- late_interval <- zones <- max_ma <- min_ma <- mid_ma <- c();
for (ru in 1:length(named_rocks))	{
	this_rock <- subset(rocky_collections,rocky_collections$rock_unit==named_rocks[ru]);
	this_rock_intervals <- unique(c(as.character(this_rock$early_interval),as.character(this_rock$late_interval[this_rock$late_interval!=""])));
	ttl_ints <- c(ttl_ints,length(this_rock_intervals));
#	strat_no <- match(this_rock_intervals,strat_info$interval_name);
#	all_intervals <- strat_info[strat_no,];
	mm <- match(max(this_rock$max_ma),this_rock$max_ma);
	max_ma <- c(max_ma,this_rock$max_ma[mm]);
	early_interval <- c(early_interval,as.character(this_rock$early_interval[mm]));
	nn <- match(min(this_rock$min_ma),this_rock$min_ma);
	min_ma <- c(min_ma,this_rock$min_ma[nn]);
	if (this_rock$late_interval[nn]!="")	{
		late_interval <- c(late_interval,as.character(this_rock$late_interval[nn]));
		} else	{
		late_interval <- c(late_interval,as.character(this_rock$early_interval[nn]));
		}
	rock_zone <- sort(unique(as.character(this_rock$zone[this_rock$zone!=""])));
	if (length(rock_zone)==1)	{
		zones <- c(zones,rock_zone)
		} else if (length(rock_zone)>1)	{
		zones <- c(zones,paste(rock_zone, collapse = ", "));
		} else	{
		zones <- c(zones,"");
		}
	other_ints <- this_rock_intervals[!this_rock_intervals %in% c(early_interval[ru],late_interval[ru])];
	if (length(other_ints)==1)	{
		other_intervals <- c(other_intervals,other_ints);
		} else if (length(other_ints)>1)	{
		other_intervals <- c(other_intervals,paste(other_ints, collapse = ", "));
		} else	{
		other_intervals <- c(other_intervals,"");
		}
	}
#xxx <- cbind(frak_it(ttl_ints),as.character(other_intervals));
rock_unit_information <- data.frame(formation=as.character(formation_member[,1]),
member=as.character(formation_member[,2]),
early_interval=as.character(early_interval),late_interval=as.character(late_interval),
zones=as.character(zones),max_ma=as.numeric(max_ma),min_ma=as.numeric(min_ma),
mid_ma=(as.numeric(max_ma)+as.numeric(min_ma))/2);

if (onset!=end)	{
	time_range <- paste(onset,"-",end,sep="")
	} else {
	time_range <- onset
	}

filename1 <- paste("Rock_Unit_Information",taxon,time_range,sep="_")
filename1 <- paste(filename1,file_format,sep="");
if (file_format!=".csv")	
	write.table(rock_unit_information,filename1,sep="\t",eol="\n",row.names=FALSE);
if (file_format==".csv")	
	write.csv(rock_unit_information,filename1,row.names=FALSE);

filename2 <- paste("Formationless_Collections",taxon,time_range,sep="_");
filename2 <- paste(filename2,file_format,sep="");
if (file_format!=".csv")	
	write.table(orphans,filename2,sep="\t",eol="\n",row.names=FALSE);
if (file_format==".csv")	
	write.csv(orphans,filename2,row.names=FALSE);
output <- list(rock_unit_information,orphans,confused_rocks);
names(output) <- c("rock_unit_information","collections_sans_formations","member_or_formation");
return(output);
}

mundare_rock_unit_names <- function(named_rock_unit,delete_rock_type=TRUE)	{
# This removes quotes, accented letters (or their corruptions), lithology addenda
# 	etc. from the rock unit name.
ttlch <- nchar(named_rock_unit)
if (ttlch>0 && named_rock_unit != " ")	{
	text <- substring(named_rock_unit, seq(1, nchar(named_rock_unit), 1), seq(1, nchar(named_rock_unit), 1))
	while (text[1]==" " && ttlch>1)	{
		named_rock_unit <- pracma::strcat(text[2:ttlch])
		text <- text[2:ttlch]
		ttlch <- ttlch-1
		}
	if (text[1]=="\"")	{
		named_rock_unit <- pracma::strcat(text[2:ttlch])
		text <- text[2:ttlch]
		}
	ttlch <- nchar(named_rock_unit)
	if (text[ttlch]=="\"")
		named_rock_unit <- pracma::strcat(text[1:(ttlch-1)])
	named_rock_unit <- gsub("\\?","",named_rock_unit)
	named_rock_unit <- gsub("\"", "",named_rock_unit)
	named_rock_unit <- gsub("Ft. \\.","Fort ",named_rock_unit)
	named_rock_unit <- gsub(" Fm\\.", "",named_rock_unit)
	named_rock_unit <- gsub(" Fm", "",named_rock_unit)
	named_rock_unit <- gsub(" Formation", "",named_rock_unit)
	named_rock_unit <- gsub(" formation", "",named_rock_unit)
	named_rock_unit <- gsub(" Member", "",named_rock_unit)
	named_rock_unit <- gsub(" member", "",named_rock_unit)
	named_rock_unit <- gsub(" Mbr ", "",named_rock_unit)
	named_rock_unit <- gsub(" Mbr", "",named_rock_unit)
	named_rock_unit <- gsub("Á","A",named_rock_unit)
	named_rock_unit <- gsub("Ä","A",named_rock_unit)
	named_rock_unit <- gsub("ä","a",named_rock_unit)
	named_rock_unit <- gsub("á","a",named_rock_unit)
	named_rock_unit <- gsub("ç","c",named_rock_unit)
	named_rock_unit <- gsub("é","e",named_rock_unit)
	named_rock_unit <- gsub("è","e",named_rock_unit)
	named_rock_unit <- gsub("ñ","n",named_rock_unit)
	named_rock_unit <- gsub("Ö","O",named_rock_unit)
	named_rock_unit <- gsub("Ø","O",named_rock_unit)
	named_rock_unit <- gsub("ø","o",named_rock_unit)
	named_rock_unit <- gsub("ó","o",named_rock_unit)
	named_rock_unit <- gsub("ö","o",named_rock_unit)
	named_rock_unit <- gsub("õ","o",named_rock_unit)
	named_rock_unit <- gsub("ü","u",named_rock_unit)
	#named_rock_unit <- gsub("’","\\'",,named_rock_unit)
	named_rock_unit <- gsub("Ã„","A",named_rock_unit)
	named_rock_unit <- gsub("Á","A",named_rock_unit)
	named_rock_unit <- gsub("Ã¥","a",named_rock_unit)
	named_rock_unit <- gsub("Ã¡","a",named_rock_unit)
	named_rock_unit <- gsub("Ã¤","a",named_rock_unit)
	named_rock_unit <- gsub("Ã¥","a",named_rock_unit)
	named_rock_unit <- gsub("Ã§","c",named_rock_unit)
	named_rock_unit <- gsub("Ã©","e",named_rock_unit)
	named_rock_unit <- gsub("Ã¨","e",named_rock_unit)
	named_rock_unit <- gsub("Ã±","n",named_rock_unit)
	named_rock_unit <- gsub("Ã–","O",named_rock_unit)
	named_rock_unit <- gsub("Ã¸","o",named_rock_unit)
	named_rock_unit <- gsub("Ã¶","o",named_rock_unit)
	named_rock_unit <- gsub("Ãµ","o",named_rock_unit)
	named_rock_unit <- gsub("Ã´","o",named_rock_unit)
	named_rock_unit <- gsub("Ã¼","u",named_rock_unit)
	#named_rock_unit <- gsub("â€™","\\'",,named_rock_unit)
	
	if (delete_rock_type)	{
		named_rock_unit <- gsub(" Ls", "",named_rock_unit)
		named_rock_unit <- gsub(" Lst", "",named_rock_unit)
		named_rock_unit <- gsub(" Ls\\.", "",named_rock_unit)
		named_rock_unit <- gsub(" Lst\\.", "",named_rock_unit)
		named_rock_unit <- gsub(" Qzt\\.", "",named_rock_unit)
		named_rock_unit <- gsub(" Quartzite", "",named_rock_unit)
		named_rock_unit <- gsub(" limestone", "",named_rock_unit)
		named_rock_unit <- gsub(" Limestone", "",named_rock_unit)
		named_rock_unit <- gsub(" Limeston", "",named_rock_unit)
		named_rock_unit <- gsub(" Dolomite", "",named_rock_unit)
		named_rock_unit <- gsub(" Sandstone", "",named_rock_unit)
		named_rock_unit <- gsub(" Shale", "",named_rock_unit)
		}
	return(named_rock_unit)
	}	else	{
	named_rock_unit <- ""
	return(named_rock_unit)
	}
}

paleodb_rock_units_for_curation <- function(onset="Cambrian", end="Holocene", taxa=c("Animalia","Protista","Plantae"),system="All",delete_rock_type=TRUE,save_files=TRUE,output_type=".txt")	{
# function for general curation of rock units in PaleoDB.
# Routine returns (and will output as files) three summaries:
#	1.  All stages (local & international) to which a rock unit is assigned
#	2.  All zones (local & international) to which a rock unit is assigned
#	3.  All combinations for rock units alternativel classified as formations
#		or members, or different formation assignments for members
# Input: onset & end for oldest & latest rocks
# taxa: leave this set in order to get basically everything
# system: "all" "marine" or "terrestrial
# delete rock type: if TRUE< then "BLANK Limestone" or "BLAH Sandstone" has
#	lithology removed from name; this make Antelope Valley and Antelope Valley
#	Limestone the same rock unit
# save_files: if true, then output is saved by function with name based on
#	time span & system.
# output_type: ending (e.g., .txt or .csv) for file.
print(c(onset,end,taxa,system,delete_rock_type,save_files,output_type));
#taxa <- paste(taxa, collapse = ","); print("Line 1147");
if (tolower(system)!="all")	{
	pbdb_sites <- accersi_collection_data(taxa=taxa,onset=onset,end=end,basic_environments=system,save_files=F);
	} else	{
	pbdb_sites <- accersi_collection_data(taxa=taxa,onset=onset,end=end,save_files=F);
	}
ttl_coll <- nrow(pbdb_sites);

formation_members <- c();
formation_members_sep <- matrix(0,ttl_coll,2);
for (c in 1:ttl_coll)	{
	formation_members_sep[c,] <- c(as.character(pbdb_sites$formation[c]),as.character(pbdb_sites$member[c]))
	if (pbdb_sites$formation[c]!="" && pbdb_sites$member[c]!="")	{
		formation_members <- c(formation_members,as.character(paste(pbdb_sites$formation[c]," [",pbdb_sites$member[c],"]",sep="")))
		} else if (pbdb_sites$formation[c]!="")	{
		formation_members <- c(formation_members,as.character(pbdb_sites$formation[c]))
		} else if (pbdb_sites$member[c]!="")	{
		formation_members <- c(formation_members,as.character(paste("[",pbdb_sites$member[c],"]",sep="")))
		}	else if (pbdb_sites$stratgroup[c]!="")	{
		ad_hoc <- paste(as.character(pbdb_sites$stratgroup[c])," Group Unnamed Formation",sep="")
		formation_members <- c(formation_members,ad_hoc)
		formation_members_sep[c,1] <- as.character(ad_hoc)
		}	else	{
		formation_members <- c(formation_members,"")
		}
#	c <- c+1
	}
pbdb_sites <- cbind(pbdb_sites,formation_members);
rock_order <- order(formation_members);
formation_members <- sort(unique(formation_members));
formation_members_sep <- unique(formation_members_sep[rock_order,]);
rock_units <- length(formation_members);

rock_unit_zones <- rock_unit_stages <- c();
for (r in 1:rock_units)	{
	if (formation_members[r]!="" && formation_members[r]!=" ")	{
		## for some reason, subset command is malfunctioning here....
		cases <- (1:ttl_coll)[pbdb_sites$formation_members %in% formation_members[r]];
		this_rock <- pbdb_sites[cases,];
		zones <- as.character(unique(this_rock$zone[this_rock$zone!=""]));
		if (length(zones)>0)	{
			zn <- length(zones);
			zones <- sort(zones);
			rock_unit_zones <- rbind(rock_unit_zones,cbind(rep(formation_members_sep[r,1],zn),rep(formation_members_sep[r,2],zn),zones));
			}
		stages <- sort(unique(c(as.character(this_rock$early_interval[this_rock$early_interval!=""]),as.character(this_rock$late_interval[this_rock$late_interval!=""]))));
		st <- length(stages);
		if (st>0)	rock_unit_stages <- rbind(rock_unit_stages,cbind(rep(formation_members_sep[r,1],st),rep(formation_members_sep[r,2],st),stages));
		}
#	r <- r+1
	}
colnames(rock_unit_stages) <- c("Formation","Member","Stage")
colnames(rock_unit_zones) <- c("Formation","Member","Zone")
rock_unit_stages <- data.frame(rock_unit_stages);
rock_unit_zones <- data.frame(rock_unit_zones);
### prepare all combinations of formations & members, zones & strat_units

#accersi_latest_rock_unit_assignment(formations=relevant_collections$formation,members=relevant_collections$member,ref_pubyr=relevant_collections$ref_pubyr,reference_no=relevant_collections$reference_no)
confused_rocks <- confused_rock_assignments(formations=pbdb_sites$formation,members=pbdb_sites$member,reference_no=pbdb_sites$reference_no,ref_pubyr=pbdb_sites$ref_pubyr,primary_reference=pbdb_sites$primary_reference);

if (save_files)	{
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output1 <- paste(timespan,"_Stage_Assignments_for_",system,"_Strata",output_type,sep="")
	output2 <- paste(timespan,"_Zone_Assignments_for_",system,"_Strata",output_type,sep="")
	output3 <- paste(timespan,"_",system,"_Strata_Entered_as_Formations_and_Members",output_type,sep="")
	output4 <- paste(timespan,"_",system,"_Relevant_Collections",output_type,sep="")
	if (output_type==".csv")	{
		write.csv(rock_unit_stages,output1,row.names = F);
		write.csv(rock_unit_zones,output2,row.names = F);
		write.csv(confused_rocks,output3,row.names = F);
		write.csv(pbdb_sites,output4,row.names = F);
		} else if (output_type %in% c(".xls",".xlsx"))	{
		writexl::write_xlsx(rock_unit_stages,output1);
		writexl::write_xlsx(rock_unit_zones,output2);
		writexl::write_xlsx(confused_rocks,output3);
		writexl::write_xlsx(species_finds,output4);
		} else	{
		write.table(rock_unit_stages,file=output1,sep = "\t",row.names = FALSE,col.names=TRUE);
		write.table(rock_unit_zones,file=output2,sep = "\t",row.names = FALSE,col.names=TRUE);
		write.table(confused_rocks,file=output3,sep = "\t",row.names = FALSE,col.names=TRUE);
		write.table(pbdb_sites,file=output4,sep = "\t",row.names = FALSE,col.names=TRUE);
		}
	}
output <- list(rock_unit_stages,rock_unit_zones,confused_rocks);
names(output) <- c("Stage_Assignments","Zone Assignments","Inconsistent_Ranks");
#print("That's all folks....")
return(output);
}

paleodb_vett_formations <- function(taxon,onset="Archean",end="Phanerozoic",basic_environments="marine,unknown,terr",directory="",file_format=".csv")	{
	# get temporal information for search
http_time <- "http://www.paleobiodb.org/data1.2/intervals/list.csv?scale=all";
strat_info <- read.table(http_time, sep=',', header=T);
strat_names <- strat_info$interval_name;
	
# MAKE SURE WE RECOGNIZE THE OLDEST AND YOUNGEST AGES!!!
if(onset %in% strat_names) {
	start_ma <- strat_info$max_ma[match(onset, strat_names)];
	} else	{
	print("Error! ",onset," is not a term in our chronostratigraphic lexicon!")
	}

if(end %in% strat_names) {
	end_ma <- strat_info$min_ma[match(end, strat_names)];
	} else {
	print("Error! ",end," is not a term in our chronostratigraphic lexicon!")
	}

#taxa <- paste(taxa, collapse = ",");
pbdb_sites <- accersi_collection_data(taxa=taxa,onset=onset,end=end,basic_environments=basic_environments,species_only=F,save_files=F);

ncoll <- nrow(pbdb_sites);
pbdb_sites$mid_ma <- as.numeric((pbdb_sites$max_ma+pbdb_sites$min_ma)/2);
#named_rock_unit <- as.character(pbdb_sites$formation[pbdb_sites$formation!=""]);
#pbdb_sites$formation[pbdb_sites$formation!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
#named_rock_unit <- as.character(pbdb_sites$member[pbdb_sites$member!=""]);
#pbdb_sites$member[pbdb_sites$member!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);
#named_rock_unit <- as.character(pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""]);
#pbdb_sites$stratgroup[pbdb_sites$stratgroup!=""] <- sapply(named_rock_unit,mundify_rock_unit_names,dehyphenate=T,delete_rock_type=T,delete_informal=F);

#zone <- pbdb_sites$zone[pbdb_sites$zone!=""];
#pbdb_sites$zone[pbdb_sites$zone!=""] <- sapply(zone,turgio_zone);

# remove quotation debris
for (i in 1:ncol(pbdb_sites))	{
	pbdb_sites[,i] <- gsub("â€œ","",pbdb_sites[,i])
	pbdb_sites[,i] <- gsub("â€\u009d","",pbdb_sites[,i])
	}
	
collections_w_group <- (1:ncoll)[pbdb_sites$stratgroup!=""];
collections_w_formations <- (1:ncoll)[pbdb_sites$formation!=""];
collections_w_members <- (1:ncoll)[pbdb_sites$member!=""];
collections_w_any_rock <- unique(sort(c(collections_w_group,collections_w_formations,collections_w_members)));
collections_w_form_and_memb <- collections_w_members[collections_w_members %in% collections_w_formations];
members_only <- collections_w_members[!collections_w_members %in% collections_w_formations];
group_only <- collections_w_group[!collections_w_group %in% c(collections_w_formations,collections_w_members)];

pbdb_sites$rock_unit <- pbdb_sites$formation;
pbdb_sites$rock_unit[collections_w_form_and_memb] <- paste(pbdb_sites$formation[collections_w_form_and_memb]," (",pbdb_sites$member[collections_w_form_and_memb],")",sep="");
#pbdb_sites$rock_unit[memb_wo_form[!memb_wo_form %in% w_group]] <- pbdb_sites$member[memb_wo_form[!memb_wo_form %in% w_group]];
pbdb_sites$rock_unit[members_only] <- pbdb_sites$member[members_only];
pbdb_sites$rock_unit[group_only] <- paste("[",pbdb_sites$stratgroup[group_only],"]",sep="");

sortdata <- c("rock_unit","mid_ma","early_interval","zone")
pbdb_sites <- pbdb_sites[do.call("order",pbdb_sites[sortdata]),];	# cool trick for sorting data!

formation_member <- unique(data.frame(formation=as.character(pbdb_sites$formation[pbdb_sites$rock_unit!=""]),
									  member=as.character(pbdb_sites$member[pbdb_sites$rock_unit!=""]),
									  stringsAsFactors = F));
any_rock2 <- (1:ncoll)[pbdb_sites$rock_unit!=""];
reduced_collections <- unique(data.frame(formation=as.character(pbdb_sites$formation[any_rock2]),
										 member=as.character(pbdb_sites$member[any_rock2]),
										 rock_unit=as.character(pbdb_sites$rock_unit[any_rock2]),
										 group=as.character(pbdb_sites$stratgroup[any_rock2]),
										 early_interval=as.character(pbdb_sites$early_interval[any_rock2]),
										 late_interval=as.character(pbdb_sites$late_interval[any_rock2]),
										 zone=as.character(pbdb_sites$zone[any_rock2]),
										 max_ma=as.numeric(pbdb_sites$max_ma[any_rock2]),
										 min_ma=as.numeric(pbdb_sites$min_ma[any_rock2]),
										 mid_ma=as.numeric(pbdb_sites$mid_ma[any_rock2]),
										 stringsAsFactors = hell_no));
formations <- unique(reduced_collections$formation[reduced_collections$formation!=""]);
rcoll <- nrow(reduced_collections);
for (ff in 1:length(formations))	{
	this_rock <- subset(reduced_collections,reduced_collections$formation==formations[ff]);
	if (nrow(this_rock)>1 && sum(this_rock$group!="")>0)	{
		diff_groups <- unique(this_rock$group[this_rock$group!=""]);
		if (length(diff_groups)==1)	{
			this_rock$group <- diff_groups;
			update <- (1:rcoll)[reduced_collections$formation %in% formations[ff]];
			reduced_collections$group[update] <- diff_groups;
			} else	{
			max_matches <- 0;
			for (dg in 1:length(diff_groups))	{
				if (max_matches < sum(this_rock$group==diff_groups[dg]))	{
					max_matches <- sum(this_rock$group==diff_groups[dg]);
					best_bet <- diff_groups[dg];
					} else if (max_matches == sum(this_rock$group==diff_groups[dg]))	{
					best_bet <- c(best_bet,diff_groups[dg]);
					}
				}
			update <- (1:rcoll)[reduced_collections$formation %in% formations[ff]];
			if (length(best_bet)==1)	{
				reduced_collections$group[update] <- best_bet;
				} else {
				latest_pub <- 1000;
				for (bb in 1:length(best_bet))	{
					this_group_and_rock <- subset(pbdb_sites,pbdb_sites$formation==formations[ff]);
					this_group_and_rock <- subset(this_group_and_rock,this_group_and_rock$stratgroup==best_bet[bb]);
					if (latest_pub < max(as.numeric(this_group_and_rock$ref_pubyr)))	{
						latest_pub < max(as.numeric(this_group_and_rock$ref_pubyr))
						best_bet_2 <- best_bet[bb];
						}
					}
				reduced_collections$group[update] <- best_bet_2[1];
				}
			}
		}
	}
group_formation_member <- unique(data.frame(formation=as.character(reduced_collections$formation),
											member=as.character(reduced_collections$member),
											group=as.character(reduced_collections$group),
											rock_unit=as.character(reduced_collections$rock_unit),
											stringsAsFactors = hell_no));

orphans <- subset(pbdb_sites,pbdb_sites$rock_unit=="");
rocky_collections <- subset(reduced_collections,reduced_collections$rock_unit!="");
named_rocks <- unique(rocky_collections$rock_unit);
entered_formations <- sort(unique(rocky_collections$formation[rocky_collections$formation!=""]));
entered_members <- sort(unique(rocky_collections$member[rocky_collections$member!=""]));

confused_rocks <- entered_formations[entered_formations %in% entered_members];
confused_rocks <- confused_rocks[!confused_rocks %in% not_real_rock_names];

ttl_ints <- other_intervals <- early_interval <- late_interval <- zones <- max_ma <- min_ma <- mid_ma <- c();
for (ru in 1:length(named_rocks))	{
	this_rock <- subset(rocky_collections,rocky_collections$rock_unit==named_rocks[ru]);
	this_rock_intervals <- unique(c(as.character(this_rock$early_interval),as.character(this_rock$late_interval[this_rock$late_interval!=""])));
	ttl_ints <- c(ttl_ints,length(this_rock_intervals));

	mm <- match(max(this_rock$max_ma),this_rock$max_ma);
	max_ma <- c(max_ma,this_rock$max_ma[mm]);
	early_interval <- c(early_interval,as.character(this_rock$early_interval[mm]));
	nn <- match(min(this_rock$min_ma),this_rock$min_ma);
	min_ma <- c(min_ma,this_rock$min_ma[nn]);
	if (this_rock$late_interval[nn]!="")	{
		late_interval <- c(late_interval,as.character(this_rock$late_interval[nn]));
		} else	{
		late_interval <- c(late_interval,as.character(this_rock$early_interval[nn]));
		}
	rock_zone <- sort(unique(as.character(this_rock$zone[this_rock$zone!=""])));
	if (length(rock_zone)==1)	{
		zones <- c(zones,rock_zone)
		} else if (length(rock_zone)>1)	{
		zones <- c(zones,paste(rock_zone, collapse = ", "));
		} else	{
		zones <- c(zones,"");
		}
	other_ints <- this_rock_intervals[!this_rock_intervals %in% c(early_interval[ru],late_interval[ru])];
	if (length(other_ints)==1)	{
		other_intervals <- c(other_intervals,other_ints);
		} else if (length(other_ints)>1)	{
		other_intervals <- c(other_intervals,paste(other_ints, collapse = ", "));
		} else	{
		other_intervals <- c(other_intervals,"");
		}
	}
#xxx <- cbind(frak_it(ttl_ints),as.character(other_intervals));
rock_unit_information <- data.frame(formation=as.character(reduced_collections$formation[match(named_rocks,reduced_collections$rock_unit)]),
									member=as.character(reduced_collections$member[match(named_rocks,reduced_collections$rock_unit)]),
									group=as.character(reduced_collections$group[match(named_rocks,reduced_collections$rock_unit)]),
									early_interval=as.character(early_interval),
									late_interval=as.character(late_interval),
									zones=as.character(zones),
									max_ma=as.numeric(max_ma),
									min_ma=as.numeric(min_ma),
									mid_ma=(as.numeric(max_ma)+as.numeric(min_ma))/2,
									stringsAsFactors = hell_no);

if (onset!=end)	{
	time_range <- paste(onset,"-",end,sep="")
	} else {
	time_range <- onset
	}

filename1 <- paste("Rock_Unit_Information",taxon,time_range,sep="_")
filename1 <- paste(filename1,file_format,sep="");
if (file_format!=".csv")	
	write.table(rock_unit_information,filename1,sep="\t",eol="\n",row.names=FALSE);
if (file_format==".csv")	
	write.csv(rock_unit_information,filename1,row.names=FALSE);

filename2 <- paste("Formationless_Collections",taxon,time_range,sep="_");
filename2 <- paste(filename2,file_format,sep="");
if (file_format!=".csv")	
	write.table(orphans,filename2,sep="\t",eol="\n",row.names=FALSE);
if (file_format==".csv")	
	write.csv(orphans,filename2,row.names=FALSE);
output <- list(rock_unit_information,orphans,confused_rocks);
names(output) <- c("rock_unit_information","collections_sans_formations","member_or_formation");
return(output);
}
