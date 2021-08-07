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
if(!require(base64enc)) install.packages("base64enc");	library(base64enc);
if(!require(encode)) install.packages("encode"); 		library(encode);
source(paste(common_source_folder,'Data_Downloading_v4.r',sep=""))  #
source(paste(common_source_folder,'Disparity.r',sep=""))  #
source(paste(common_source_folder,'Nexus_File_Routines.r',sep=""))  #
source(paste(common_source_folder,'Wagner_kluges.r',sep=""))  #
load(paste(data_for_R_folder,'Character_Data.RData',sep=""));
INAP <- -22; UNKNOWN <- -11;

#### LOAD DATA ####
data_file <- file.choose();															# choose a nexus file
matrix_name <- strsplit(data_file,"/")[[1]][length(strsplit(data_file,"/")[[1]])];	# get file name alone
character_info <- accersi_data_from_nexus_file(data_file);							# get information about nexus file
chmatrix <- character_info$Matrix;													# get character data
notu <- nrow(chmatrix);

# Get Dependent - Independent Pairs ####
nchars <- ncol(chmatrix);															# total characters
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
redone_1 <- redone_2 <- list();
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
	combos <- unique(chmatrix[,c(ind_char,dep_chars)]);
	out_states <- unique(chmatrix[chmatrix[,dep_chars[1]]==INAP,ind_char]);
	chmatrix[chmatrix[,ind_char] %in% out_states,dep_chars] <- INAP;

	redone_1 <- rlist::list.append(redone_1,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix,INAP,UNKNOWN,multichanges=F,theoretical=F,unknown_inap_sep=F));
	redone_2 <- rlist::list.append(redone_2,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix,INAP,UNKNOWN,multichanges=F,theoretical=T,unknown_inap_sep=F));
	}
names(redone_1) <- names(redone_2) <- unique_indies;

# compare transition matrices under both schemes
pi2_1 <- expm::expm(0.05*redone_1[[2]]$Q);
pi2_2 <- expm::expm(0.05*redone_2[[2]]$Q);

unadditives <- additives[additives$independents==0,];

if (nrow(unadditives)>0)	{
	new_indies <- array(-11,dim=c(notu,nrow(unadditives)));
	colnames(new_indies) <- unadditives$dependents;
	for (ua in 1:nrow(unadditives))	{
		new_indies[chmatrix[,unadditives$dependents[ua]]==INAP,ua] <- 0;
		new_indies[chmatrix[,unadditives$dependents[ua]]>=0,ua] <- 1;
		polys <- (1:notu)[chmatrix[,unadditives$dependents[ua]]<0];
		polys <- polys[!chmatrix[polys,unadditives$dependents[ua]] %in% c(INAP,UNKNOWN)];
		new_indies[polys,ua] <- 1;
		unadditives$independents[ua] <- ua;
		}
	}

ua <- 1;
while (ua < nrow(unadditives))	{
	inap_taxa <- (1:notu)[chmatrix[,unadditives$dependents[ua]]==INAP];
	au <- ua+1;
	while (au <= nrow(unadditives))	{
		inap_taxa_2 <- (1:notu)[chmatrix[,unadditives$dependents[au]]==INAP];
		if (sum(inap_taxa_2 %in% inap_taxa)==length(inap_taxa_2))	unadditives$independents[au] <- unadditives$independents[ua]
		au <- au+1;
		}
	ua <- ua+1;
	}

unique_indies2 <- unique(unadditives$independents);
for (ui in 1:length(unique_indies2))	{
	ind_char <- unique_indies2[ui];
	dep_chars <- unadditives$dependents[unadditives$independents==ind_char];
	chmatrix_d <- chmatrix;
	chmatrix_d[,ind_char] <- new_indies[,ui];

	redone_1 <- rlist::list.append(redone_1,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix=chmatrix_d,INAP,UNKNOWN,theoretical=F,unknown_inap_sep=F));
	redone_2 <- rlist::list.append(redone_2,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix=chmatrix_d,INAP,UNKNOWN,theoretical=T,unknown_inap_sep=F));
	}

names(redone_1) <- names(redone_2) <- c(names(redone_1)[names(redone_1)!=""],unadditives$dependents[match(unique_indies2,unadditives$independents)]);


{}
