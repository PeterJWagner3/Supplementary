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
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";					# directory to folder where you keep R-Data files
source(paste(common_source_folder,'Data_Downloading_v4.r',sep=""))  #
source(paste(common_source_folder,'Nexus_File_Routines.r',sep=""))  #
source(paste(common_source_folder,'Wagner_kluges.r',sep=""))  #
load(paste(data_for_R_folder,'Character_Data.RData',sep=""));
INAP <- -22; UNKNOWN <- -11;
brachiopod_matrices <- names(character_database$Lophophorates$Brachiopods);

bm <- 27;
matrix_name <- brachiopod_matrices[bm];
#character_info <- accersi_data_from_chosen_nexus_file(polymorphs = T,UNKNOWN=-11,INAP=-22);
#character_info <- accersi_data_from_RData(matrix_name=matrix_name,character_database);
chmatrix <- character_info$Matrix;
dependents <- unique(which(chmatrix==INAP,arr.ind=T)[,2]);
nchars <- ncol(chmatrix);
independents <- (1:nchars)[!(1:nchars) %in% dependents];
dchar <- dependents;
parent_chars <- sapply(dchar,find_independent_character,independents,chmatrix,UNKNOWN,INAP);
additives <- data.frame(dependents=as.numeric(dependents),independents=as.numeric(parent_chars));

unique_indies <- unique(additives$independents[additives$independents>0]);
redone <- list();
for (ui in 1:length(unique_indies))	{
	ind_char <- unique_indies[ui];
	dep_chars <- additives$dependents[additives$independents==ind_char];
	redone <- rlist::list.append(redone,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix,INAP,UNKNOWN));
	}




{}
