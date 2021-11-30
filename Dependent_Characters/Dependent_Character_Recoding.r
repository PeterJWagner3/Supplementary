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
source(paste(common_source_folder,'disparity.r',sep=""))  #
source(paste(common_source_folder,'Nexus_File_Routines.r',sep=""))  #
source(paste(common_source_folder,'paleophylogeny_routines.r',sep=""))  #
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
notu <- nrow(chmatrix);				nchars <- ncol(chmatrix);															# total characters
nstates <- character_info$States;	chtypes <- character_info$State_Types;	nchars <- length(character_info$States);

# Get Dependent - Independent Pairs ####
dchar <- unique(which(chmatrix==INAP,arr.ind=T)[,2]);								# separate characters with inapplicables
independents <- (1:nchars);
# get the "parent" character upon which each dependent character depends
#for (dc in 1:length(dchar))	dd <- find_independent_character(dchar=dchar[dc],independents,chmatrix,UNKNOWN,INAP);
parent_chars <- pbapply::pbsapply(dchar,find_independent_character,independents,chmatrix,UNKNOWN,INAP);
# create data.frame with additive characters, with dependent & independent.
additives <- data.frame(dependents=as.numeric(dchar),independents=as.numeric(parent_chars));
#unique(chmatrix[,c(46:48)])
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
	char_dependencies <- c(ind_char,additives$independents[match(dep_chars,additives$dependents)])
	combos <- unique(chmatrix[,c(ind_char,dep_chars)]);
	combos <- combos[!(rowMaxs(combos)==UNKNOWN & rowMins(combos)==UNKNOWN),];
	out_states <- unique(chmatrix[chmatrix[,dep_chars[1]]==INAP,ind_char]);
	accidental_unknowns <- which(chmatrix[chmatrix[,ind_char] %in% out_states,dep_chars]==UNKNOWN,arr.ind = T);
	if (length(accidental_unknowns)>0)	{
		if (!is.matrix(accidental_unknowns))	accidental_unknowns <- array(accidental_unknowns,dim=c(1,2));
		for (au in 1:nrow(accidental_unknowns))	{
			fch <- c(ind_char,dep_chars)[accidental_unknowns[2]];
			chmatrix[accidental_unknowns[1],fch] <- INAP;
			}
		}
	
	redone <- rlist::list.append(redone,transmogrify_additive_dependents_to_multistate(ind_char,dep_chars,chmatrix,char_dependencies,INAP,UNKNOWN));
	}
names(redone) <- unique_indies;

tstates <- tacit_indies <- additives$dependents[additives$independents==0];
names(tstates) <- tacit_indies;
redone_2 <- list();
redone_all <- redone;
ti <- 0;
#for (ti in 1:length(tacit_indies))		{
while (ti < length(tacit_indies))	{
	ti <- ti+1;
	tch <- tacit_indies[ti];
	new_states <- orig_states <- chmatrix[,tch];
	polys <- (1:notu)[orig_states < 0 & !orig_states %in% c(UNKNOWN,INAP)];
	new_states[orig_states>=0] <- orig_states[orig_states>=0]+1;
	new_states[orig_states==INAP] <- 0;
	pl <- 0;
	while (pl < length(polys))	{
		pl <- pl+1
		new_states[polys[pl]] <- ravel_polymorph(unravel_polymorph_badass(orig_states[polys[pl]])+1);
#		print(ti);
		}
	tstates[ti] <- max(new_states+1);
	Q <- array(1,dim=c(tstates[ti],tstates[ti]));
	Q[1,] <- 1/nstates[tch];
	diag(Q) <- 0; 
	diag(Q) <- -(rowSums(Q));
	colnames(Q) <- rownames(Q) <- c("-",(1:nstates[tch])-1)
	x <- list(orig_states,Q,new_states)
	names(x) <- names(redone[[1]]);
	redone_2 <- rlist::list.append(redone_2,x);
	redone_all <- rlist::list.append(redone_all,x);
	}
names(redone_2) <- tacit_indies;
names(redone_all) <- c(names(redone),names(redone_2));
#order(as.numeric(names(redone_all)))
redone_all <- redone_all[order(as.numeric(names(redone_all)))];
#names(redone_all)
new_matrix <- array("",dim=c(notu,length(redone_all)));
for (ui in 1:ncol(new_matrix))	new_matrix[,ui] <- redone_all[[ui]]$new_character;
colnames(new_matrix) <- names(redone_all);
rownames(new_matrix) <- character_info$OTUs;

chmatrix_recoded <- convert_character_matrix_to_character(chmatrix);
replacements <- match(as.numeric(colnames(new_matrix)),as.numeric(colnames(chmatrix_recoded)));
orig_colnames <- colnames(chmatrix_recoded);
chmatrix_recoded[,replacements] <- new_matrix;
retentions <- (1:nchars)[!(1:nchars) %in% additives$dependents];
chmatrix_recoded <- chmatrix_recoded[,retentions];

rchars <- ncol(chmatrix_recoded);
names(nstates) <- 1:nchars;
rstates <- nstates[names(nstates) %in% as.numeric(colnames(chmatrix_recoded))];
additives$independents[additives$independents==0] <- additives$dependents[additives$independents==0];

for (ai in 1:length(redone_all))	{
	rs <- match(as.numeric(names(redone_all)[ai]),names(rstates));
	rstates[rs] <- nrow(redone_all[[ai]]$Q);
	}
max_states <- max(rstates);

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
#  ADJUST THIS to exclude altered hierarchical characters
rbstates <- rstates[names(rstates)!=""];
rbtypes <- chtypes[as.numeric(names(rstates[names(rstates)!=""]))];
retained_chars <- as.numeric(names(rbstates));
char_type_combos <- char_types <- data.frame(character=as.numeric(retained_chars),
						 states=as.numeric(rbstates),
						 ordering=as.character(rbtypes),stringsAsFactors = F);
char_type_combos$character <- NULL;
char_type_combos <- unique(char_type_combos[order(char_type_combos[,1]),]);

hierarchical <- as.numeric(names(redone_all));
for (nn in 1:nrow(char_type_combos))	{
	this_set <- subset(char_types,char_types$states==char_type_combos$states[nn]);
	this_set <- subset(this_set,this_set$ordering==char_type_combos$ordering[nn]);
	this_set <- this_set[!this_set$character %in% hierarchical,];
	if (nrow(this_set)>0)	{
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
	}
# start here: break up
lgm <- length(greek_to_me);
while (length(redone_all) > length(greek_to_me))	{
	andreikelo <- greek_to_me[(1+length(greek_to_me)-lgm):length(greek_to_me)];
	for (lg in 1:lgm)	andreikelo[lg] <- paste(andreikelo[lg],"_",greek_to_me[lg],sep="");
	greek_to_me <- c(greek_to_me,andreikelo);
	}
for (nn in 1:length(redone))	{
	new_file_name_a <- gsub("Dummy",paste("Rescored",greek_to_me[nn],sep="_"),revbayes_file_name);
#	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch=ch_matrix_ch[,this_set$character],new_file_name=new_file_name_a,max_states=char_type_combos$states[nn]);
	max_states <- nrow(redone[[nn]]$Q);
	nc <- match(as.numeric(names(redone)[nn]),as.numeric(colnames(chmatrix_recoded)));
	another_prop <- array(chmatrix_recoded[,nc],dim=c(notu,1));
	rownames(another_prop) <- rownames(chmatrix);
	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch = another_prop,new_file_name = new_file_name_a,max_states=max_states);
	}
redone_2_states <- c();
nnn <- 0;
while (nnn < length(redone_2))	{
	nnn <- nnn+1;
	redone_2_states <- c(redone_2_states,nrow(redone_2[[nnn]]$Q));
	}
# now, consolidate tacit hierarcicals 	

for (nnn in 1:length(redone_2))	{
	new_file_name_a <- gsub("Dummy",paste("Rescored",greek_to_me[nn+nnn],sep="_"),revbayes_file_name);
#	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch=ch_matrix_ch[,this_set$character],new_file_name=new_file_name_a,max_states=char_type_combos$states[nn]);
	max_states <- nrow(redone[[nn]]$Q);
	nc <- match(as.numeric(names(redone)[nn]),as.numeric(colnames(chmatrix_recoded)));
	another_prop <- array(chmatrix_recoded[,nc],dim=c(notu,1));
	rownames(another_prop) <- rownames(chmatrix);
	scribio_nexus_file_from_chmatrix_character(ch_matrix_ch = another_prop,new_file_name = new_file_name_a,max_states=max_states);
	}

for (nn in 1:length(redone))	{
	new_file_name_b <- gsub("Dummy",paste("Q_Matrix",greek_to_me[nn],sep="_"),revbayes_file_name);
	new_file_name_b <- gsub(".nex",".txt",new_file_name_b);
	write.table(redone[[nn]]$Q,new_file_name_b,sep="\t");
	}

for (nnn in 1:length(redone_2))	{
	new_file_name_b <- gsub("Dummy",paste("Q_Matrix",greek_to_me[nn+nnn],sep="_"),revbayes_file_name);
	new_file_name_b <- gsub(".nex",".txt",new_file_name_b);
	write.table(redone[[nn]]$Q,new_file_name_b,sep="\t");
	}

nchars <- ncol(chmatrix);
max_states <- max(character_info$States);
rstates <- nstates;
names(rstates) <- 1:nchars;
# go through all characters; do the routine below for those with coded independents;
#		just swap out character vector for those with tacit independents
testing <- c();
for (ui in length(unique_indies):1)	{
	testing <- cbind(redone_all[[ui]]$new_character,testing);
	colnames(testing) <- names(redone_all[length(unique_indies):ui]);
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
	new_charnames <- colnames(chmatrix_recoded)[c(greta,ind_char,garbo)];
	chmatrix_recoded <- cbind(chmatrix_recoded[,greta],redone_all[[ui]]$new_character,chmatrix_recoded[,garbo]);
	max_states <- max(max_states,nrow(redone_all[[ui]]$Q));
	rstates <- c(rstates[greta],nrow(redone_all[[ui]]$Q),rstates[garbo]);
	names(rstates) <- colnames(chmatrix_recoded) <- new_charnames;
	}

{}
