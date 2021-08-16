#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal
# tired of search & replace? Use Rename in Scope
hell_no <- FALSE;	# well, it's true.
dependent_directory <- "~/Documents/R_Projects/Dependent_Characters/"
setwd(dependent_directory);
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";					# directory to folder where you keep R-Data files
library('greekLetters');
greek_to_me <- c("Alpha","Beta","Gamma","Delta","Epsilon","Zeta","Eta","Theta","Iota","Kappa","Lambda","Mu","Nu","Xi","Omicron","Pi","Rho","Sigma","Tau","Upsilon","Phi","Chi","Psi","Omega");

#if(!require(base64enc)) install.packages("base64enc");	library(base64enc);
#if(!require(encode)) install.packages("encode"); 		library(encode);
source(paste(common_source_folder,'Data_Downloading_v4.r',sep=""))  #
source(paste(common_source_folder,'Disparity.r',sep=""))  #
source(paste(common_source_folder,'Nexus_File_Routines.r',sep=""))  #
source(paste(common_source_folder,'Wagner_kluges.r',sep=""))  #
source(paste(common_source_folder,"Wagner_Stats_and_Probability_101.r",sep=""))  #
load(paste(data_for_R_folder,'Character_Data.RData',sep=""));
INAP <- -22; UNKNOWN <- -11;
allstates <- c(0:9,letter_states,more_letter_states);

#### LOAD DATA ####
data_file <- file.choose();															# choose a nexus file
matrix_name <- strsplit(data_file,"/")[[1]][length(strsplit(data_file,"/")[[1]])];	# get file name alone
character_info <- accersi_data_from_nexus_file(data_file);							# get information about nexus file
chmatrix <- character_info$Matrix;													# get character data
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);															# total characters
nstates <- character_info$States;
chtypes <- character_info$State_Types;

# Get Dependent - Independent Pairs ####
dchar <- unique(which(chmatrix==INAP,arr.ind=T)[,2]);								# separate characters with inapplicables
independents <- (1:nchars);
# get the "parent" character upon which each dependent character depends
parent_chars <- pbapply::pbsapply(dchar,find_independent_character,independents,chmatrix,UNKNOWN,INAP);
#find_independent_character(51,independents,chmatrix,UNKNOWN,INAP);
# create data.frame with additive characters, with dependent & independent.
additives <- data.frame(dependents=as.numeric(dchar),independents=as.numeric(parent_chars));

# Get the Unique Character State combinations & recodings ####
unique_indies <- unique(additives$independents[additives$independents>0]);
unique_indies <- unique_indies[!unique_indies %in% additives$dependents];
redone <- list();
for (ui in 1:length(unique_indies))	{
	ind_char <- unique_indies[ui];
	dep_chars <- additives$dependents[additives$independents==ind_char];
	tag_alongs <- additives$dependents[additives$independents %in% dep_chars];
	dep_chars <- sort(c(dep_chars,tag_alongs));
	while (length(tag_alongs)>0)	{
		new_deps <- additives$dependents[additives$independents %in% tag_alongs];
		dep_chars <- sort(c(dep_chars,new_deps));
		tag_alongs <- additives$dependents[additives$independents %in% new_deps];
		}
	secondary_dependencies <- c(ind_char,additives$independents[match(dep_chars,additives$dependents)])
	combos <- unique(chmatrix[,c(ind_char,dep_chars)]);
	combos <- combos[!(rowMaxs(combos)==UNKNOWN & rowMins(combos)==UNKNOWN),];
	out_states <- unique(chmatrix[chmatrix[,dep_chars[1]]==INAP,ind_char]);
	chmatrix[chmatrix[,ind_char] %in% out_states,dep_chars] <- INAP;

	redone <- rlist::list.append(redone,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix,secondary_dependencies=secondary_dependencies,INAP=INAP,UNKNOWN=UNKNOWN,theoretical=T));
	}
names(redone) <- unique_indies;

chmatrix_recoded <- convert_character_matrix_to_character(chmatrix);
nchars <- ncol(chmatrix);
max_states <- max(character_info$States);
rstates <- nstates;
names(rstates) <- 1:nchars;
for (ui in length(unique_indies):1)	{
	rchar <- ncol(chmatrix_recoded);
	ind_char <- unique_indies[ui];
	dep_chars <- additives$dependents[additives$independents==ind_char];
	tag_alongs <- additives$dependents[additives$independents %in% dep_chars];
	dep_chars <- sort(c(dep_chars,tag_alongs));
	while (length(tag_alongs)>0)	{
		new_deps <- additives$dependents[additives$independents %in% tag_alongs];
		dep_chars <- sort(c(dep_chars,new_deps));
		tag_alongs <- additives$dependents[additives$independents %in% new_deps];
		}
	garbo <- (1:rchar)[!(1:rchar) %in% c((1:ind_char),dep_chars)];
	greta <- (1:rchar)[!(1:rchar) %in% c(garbo,ind_char,dep_chars)];
	chmatrix_recoded <- cbind(chmatrix_recoded[,greta],redone[[ui]]$new_character,chmatrix_recoded[,garbo]);
	max_states <- max(max_states,nrow(redone[[ui]]$Q));
	rstates <- c(rstates[greta],nrow(redone[[ui]]$Q),rstates[garbo]);
	}

colnames(chmatrix_recoded) <- names(rstates);
rchars <- ncol(chmatrix_recoded);

matrix_name_parts <- strsplit(matrix_name,"_")[[1]];
coauthor_split <- match("&",matrix_name_parts);
etal_split <- match("et",matrix_name_parts);
hyphen_split <- match("-",matrix_name_parts);
if (is.na(coauthor_split) & is.na(etal_split) & is.na(hyphen_split))	{
	new_file_name <- paste(c(matrix_name_parts[1:(length(matrix_name_parts)-2)],"Recoded",matrix_name_parts[(length(matrix_name_parts)-1):length(matrix_name_parts)]),collapse="_")
	revbayes_file_name <- paste(c(matrix_name_parts[1:(length(matrix_name_parts)-2)],"Dummy",matrix_name_parts[(length(matrix_name_parts)-1):length(matrix_name_parts)]),collapse="_");
	} else if (!is.na(coauthor_split))	{
	new_file_name <- paste(c(matrix_name_parts[1:(coauthor_split-2)],"Recoded",matrix_name_parts[(coauthor_split-1):length(matrix_name_parts)]),collapse="_");
	revbayes_file_name <- paste(c(matrix_name_parts[1:(coauthor_split-2)],"Dummy",matrix_name_parts[(coauthor_split-1):length(matrix_name_parts)]),collapse="_");
	} else if (!is.na(etal_split))	{
	new_file_name <- paste(c(matrix_name_parts[1:(etal_split-2)],"Recoded",matrix_name_parts[(etal_split-1):length(matrix_name_parts)]),collapse="_");
	revbayes_file_name <- paste(c(matrix_name_parts[1:(etal_split-2)],"Dummy",matrix_name_parts[(etal_split-1):length(matrix_name_parts)]),collapse="_");
	}
scribio_nexus_file_from_chmatrix_character(ch_matrix_ch=chmatrix_recoded,new_file_name,max_states);

# now, break up matrices for RevBayes
rbstates <- rstates[names(rstates)!=""];
rbtypes <- chtypes[as.numeric(names(rstates[names(rstates)!=""]))];
retained_chars <- as.numeric(names(rbstates));
char_type_combos <- char_types <- data.frame(character=as.numeric(retained_chars),
						 states=as.numeric(rbstates),
						 ordering=as.character(rbtypes),stringsAsFactors = F);
char_type_combos$character <- NULL;
char_type_combos <- unique(char_type_combos[order(char_type_combos[,1]),]);

for (nn in 1:nrow(char_type_combos))	{
	this_set <- subset(char_types,char_types$states==char_type_combos$states[nn]);
	this_set <- subset(this_set,this_set$ordering==char_type_combos$ordering[nn]);
	dummy_sub <- paste(char_type_combos$states[nn],"States",char_type_combos$ordering[nn],sep="_");
	new_file_name_a <- gsub("Dummy",dummy_sub,revbayes_file_name);
#	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch=ch_matrix_ch[,this_set$character],new_file_name=new_file_name_a,max_states=char_type_combos$states[nn]);
	if (length(this_set$character)==1)	{
		another_prop <- array(chmatrix[,this_set$character],dim=c(notu,1))
		rownames(another_prop) <- rownames(chmatrix);
		scribio_nexus_file_from_chmatrix(ch_matrix=another_prop,new_file_name=new_file_name_a);
		} else	{
		scribio_nexus_file_from_chmatrix(ch_matrix=chmatrix[,this_set$character],new_file_name=new_file_name_a);
		}
	}
# start here: break up

rescored <- (1:rchars)[names(rstates)==""];
for (nn in 1:length(rescored))	{
	new_file_name_a <- gsub("Dummy",paste("Rescored",greek_to_me[nn],sep="_"),revbayes_file_name);
#	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch=ch_matrix_ch[,this_set$character],new_file_name=new_file_name_a,max_states=char_type_combos$states[nn]);
	max_states <- nrow(redone[[nn]]$Q);
	another_prop <- array(chmatrix_recoded[,rescored[nn]],dim=c(notu,1));
	rownames(another_prop) <- rownames(chmatrix);
	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch = another_prop,new_file_name = new_file_name_a,max_states=max_states);
	}

{}
