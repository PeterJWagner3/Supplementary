#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal

### udpated with routines written for first phylogenetics course at UNL
MAXNO <- 1.797693e+308;
newick_verbotten <- c(".","?","\"","\'");
letter_states <- LETTERS[!LETTERS %in% c("I","O")];
more_letter_states <- tolower(letter_states);
zzzz <- 0.25;
UNKNOWN <- -11;
INAP <- -22;
polymorphs=TRUE;

#### HOUSE-CLEANING ####
clear_matrix_na_with_another_cell_value <- function(data,j, k)	{
size <- dim(data)
for (i in 1:size[1])	{
	if (is.na(data[i,j]))	data[i,j] <- data[i,k]
	}
return(data)
}

count_characters_in_string <- function(string_to_count)	{
j <- strsplit(string_to_count,split="",fixed=TRUE)[[1]];
return(length(j));
}

mundify_nexus_text <- function(nexus_line)	{
#nexus_line <- gsub("\xd4","",nexus_line);
#nexus_line <- gsub("\xd5","",nexus_line);
#nexus_line <- gsub("\x87","a",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\xa3\xa4\xbc","n",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x93\xa0\xbc","ae",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x83\xa0\xbc","c",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x96\x8c\xbc","",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8c\x93\xa4\xbc","\'",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\xb3\xa0\xbc","",nexus_line);
nexus_line <- gsub("\'\'\'","\'\"",nexus_line);
nexus_line <- gsub("\'\' ","\" ",nexus_line);
#nexus_line <- gsub("\xac","",nexus_line);
#nexus_line <- gsub("\xa0","†",nexus_line);
#nexus_line <- gsub("\x80","",nexus_line);
#nexus_line <- gsub("\xd1","",nexus_line);
#nexus_line <- gsub("\xc9","?",nexus_line);
#nexus_line <- gsub("\xe1","a",nexus_line);
#nexus_line <- gsub("\xe9","e",nexus_line);
#nexus_line <- gsub("\x8e","e",nexus_line);
#nexus_line <- gsub("\x8f","e",nexus_line);
#nexus_line <- gsub("\x92","i",nexus_line);
#nexus_line <- gsub("\xbf","o",nexus_line);
#nexus_line <- gsub("\x9a","o",nexus_line);
#nexus_line <- gsub("\x97","o",nexus_line);
#nexus_line <- gsub("\xf6","ö",nexus_line);
#nexus_line <- gsub("\xfc","ue",nexus_line);
#nexus_line <- gsub("\xb0","˚",nexus_line);
#nexus_line <- gsub("\xba","˚",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8e\x93\xa4\xbc","o",nexus_line);
#nexus_line <- gsub("\x9f","ue",nexus_line);
#nexus_line <- gsub("\xd0","-",nexus_line);
#nexus_line <- gsub("\xd2","\"",nexus_line);
#nexus_line <- gsub("\xd3","\"",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x8d\x86\x90\xbc","\'",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x8d\x86\x8c\xbc","ƒ",nexus_line);
#nexus_line <- gsub("\xdf","ß",nexus_line);
#nexus_line <- gsub("\xa7","ß",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x8c\xa6\x88\xbc","≤",nexus_line);
#nexus_line <- gsub("\xb3","≥",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x8d\x96\x8c\xbc","≈",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x98\xa6\x98\xbc","˚",nexus_line);
#nexus_line <- gsub("\xb6","∂",nexus_line);
#nexus_line <- gsub("\xc6","∆",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x8d\xb6\x88\xbc","∑",nexus_line);
#nexus_line <- gsub("\xfc\xbe\x99\x86\x88\xbc","Ω",nexus_line);
#nexus_line <- gsub("\xa5"," ",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ä","A",nexus_line);
nexus_line <- gsub("ä","a",nexus_line);
nexus_line <- gsub("á","a",nexus_line);
nexus_line <- gsub("å","a",nexus_line);
nexus_line <- gsub("Ç","C",nexus_line);
nexus_line <- gsub("ç","c",nexus_line);
nexus_line <- gsub("č","c",nexus_line);
nexus_line <- gsub("é","e",nexus_line);
nexus_line <- gsub("è","e",nexus_line);
nexus_line <- gsub("ê","e",nexus_line);
nexus_line <- gsub("ė","e",nexus_line);
nexus_line <- gsub("î","i",nexus_line);
nexus_line <- gsub("Î","I",nexus_line);
nexus_line <- gsub("ñ","n",nexus_line);
nexus_line <- gsub("Ö","O",nexus_line);
nexus_line <- gsub("Ø","O",nexus_line);
nexus_line <- gsub("ø","o",nexus_line);
nexus_line <- gsub("ó","o",nexus_line);
nexus_line <- gsub("ö","o",nexus_line);
nexus_line <- gsub("õ","o",nexus_line);
nexus_line <- gsub("Š","S",nexus_line);
nexus_line <- gsub("š","s",nexus_line);
nexus_line <- gsub("ů","u",nexus_line);
nexus_line <- gsub("ü","u",nexus_line);
nexus_line <- gsub("’","’",nexus_line);
#nexus_line <- gsub("\x88","a",nexus_line);
nexus_line <- gsub("Ã„","A",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ã¡","a",nexus_line);
nexus_line <- gsub("Ã¤","a",nexus_line);
nexus_line <- gsub("Ã¥","a",nexus_line);
nexus_line <- gsub("Ã§","c",nexus_line);
nexus_line <- gsub("Ã©","e",nexus_line);
nexus_line <- gsub("Ã¨","e",nexus_line);
nexus_line <- gsub("Ã±","n",nexus_line);
nexus_line <- gsub("Ã–","O",nexus_line);
nexus_line <- gsub("Ã¸","o",nexus_line);
nexus_line <- gsub("Ã¶","o",nexus_line);
nexus_line <- gsub("Ãµ","o",nexus_line);
nexus_line <- gsub("Ã´","o",nexus_line);
nexus_line <- gsub("Ã¼","u",nexus_line);
nexus_line <- gsub("√Æ","i",nexus_line);
nexus_line <- gsub("≈†","S",nexus_line);
nexus_line <- gsub("≈°","s",nexus_line);
nexus_line <- gsub("√•","a",nexus_line);
nexus_line <- gsub("&#367;","u",nexus_line);
nexus_line <- gsub("&#945;","α",nexus_line);
return(nexus_line);
}

# Turn Lophospira sp. or Lophospira sp. A to Lophospira
reduce_genus_sp_to_genus <- function(taxon_name)	{
taxon_names <- strsplit(taxon_name," ")[[1]];
indet_species <- c("sp.");
indet_species <- c(indet_species,paste("sp.",LETTERS));
indet_species <- c(indet_species,paste("sp.",letters));
indet_species <- c(indet_species,"spp.");
indet_species <- c(indet_species,paste("spp.",LETTERS));
indet_species <- c(indet_species,paste("spp.",letters));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp. nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("indet.",i));

if (sum(taxon_names %in% indet_species)>0)
	taxon_names <- taxon_names[(1:length(taxon_names))[!taxon_names %in% indet_species]];
taxon_name <- paste(taxon_names,collapse=" ");
return(taxon_name);
}

accersi_study_name <- function(nexus_file_name)	{
n_f_n <- simplify2array(strsplit(nexus_file_name,"/"));
n_f_n <- n_f_n[length(n_f_n)];
n_f_n <- gsub("\\.nex","",n_f_n);
filename_parts <- simplify2array(strsplit(n_f_n,"_"));

if (sum(filename_parts %in% "&")==1)	{
	author_start <- match("&",filename_parts)-1;
	} else if (sum(filename_parts %in% c("et","al"))==2)	{
	author_start <- match("et",filename_parts)-1;
	} else	{
	author_start <- length(filename_parts)-2;
	}
study_group <- filename_parts[1:(author_start-1)];
citation <- paste("(",paste(filename_parts[!filename_parts %in% study_group],collapse=" "),")",sep="");
return(paste(study_group,citation));
}

name_length <- function(taxon_name)	{
return(length(strsplit(taxon_name,"")[[1]]))
}

#### ROUTINES TO READ CHARACTER MATRIX INFORMATION IN NEXUS FILE ####
# routine to read nexus file of Mesquite or Maclade format & return important infromation
accersi_data_from_nexus_file <- function(nexus_file_name, polymorphs=TRUE, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
output <- accersi_data_from_nexus_vector(nexus=nexus,polymorphs=polymorphs,UNKNOWN=UNKNOWN,INAP=INAP,rate_partitions=rate_partitions,trend_partitions=trend_partitions);
return(output)
}

# routine to prompt user for a nexus file of Mesquite or Maclade format & return important infromation
accersi_data_from_chosen_nexus_file <- function(polymorphs=TRUE, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
#print("Choose the nexus file you with to analyze: ");
print("Choose the nexus file that you wish to analyze: ");
flush.console();
Sys.sleep(zzzz);
nexus_file_name <- file.choose();
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
output <- accersi_data_from_nexus_vector(nexus=nexus,polymorphs=polymorphs,UNKNOWN=UNKNOWN,INAP=INAP,rate_partitions=rate_partitions,trend_partitions=trend_partitions);
return(output)
}

# routine to read nexus information Mesquite or Maclade format & return important infromation
accersi_data_from_nexus_vector <- function(nexus, polymorphs=TRUE, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
ml <- 0;
#i <- 1

for (i in 1:length(nexus))  {
	nexus[i] <- mundify_nexus_text(nexus_line = nexus[i]);
#	j <- simplify2array(strsplit(nexus[i],split="",fixed=TRUE))[,1];
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]];
#	strsplit(string=nexus[i],pattern="")
	if (length(j)>ml) ml <- length(j);
	}
ml <- ml+1;	# LENGTH OF LONGEST LINE

# file is now a vector of characters.  Turn it into a matrix with one char per cell
nexusfile <- matrix("\n",length(nexus),ml);
for (i in 1:length(nexus))  {
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]];
	for (k in 1:length(j))      nexusfile[i,k] <- j[k]
	if ((length(j)+2)<ml)
		for (k in (length(j)+2):ml) nexusfile[i,k] <- ""
	}

top <- match("matrix",tolower(nexus));
if (is.na(top)) top <- match("\tmatrix",tolower(nexus));
if (is.na(top))	{
	top <- 0;
	ln <- 1;		# this is the row with the word "Matrix": character data starts next.
	while (top==0)	{
		em_nexus <- gsub("\t","",nexus[ln]);
		nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
		if (!is.na(match("matrix",tolower(nexus_words))))	{
			top <- ln;
			} else	ln <- ln+1;
		}
	}
top <- top+1;	# this will give the first row of data
# skip the comment text denoting character numbers (if present)
while(nexusfile[top,1]=="[" || nexusfile[top,1]==" ") top <- top+1;

# get character & state informations
all_states <- c();
missing <- "?";
gap <- "-";
notu <- nchars <- strat <- range <- geog <- 0;
for (i in 2:top)  {
	while ((nexusfile[i,1]=="[" || nexusfile[i,1]=="\n") && i<top)	i <- i+1;
	em_nexus <- gsub("\t","",nexus[i]);
	em_nexus <- gsub("="," = ",em_nexus);
	em_nexus <- gsub(";"," ; ",em_nexus);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	nexus_words <- nexus_words[nexus_words!=""];
	if (!is.na(match("ntax",tolower(nexus_words))) || !is.na(match("ntaxa",tolower(nexus_words))))	{
		j <- 1+match("ntax",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("ntaxa",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		notu <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("nchar",tolower(nexus_words))) || !is.na(match("nchars",tolower(nexus_words))))	{
		j <- 1+match("nchar",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("nchars",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		nchars <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("gap",tolower(nexus_words))))	{
		if (nexus_words[match("gap",tolower(nexus_words))+1]=="=")	{
			j <- 1+match("gap",tolower(nexus_words));
			while(nexus_words[j]=="=")	j <- j+1;
			gap <- nexus_words[j];
			}
		}
	if (!is.na(match("missing",tolower(nexus_words))))	{
		if (nexus_words[match("missing",tolower(nexus_words))+1]=="=")	{
			j <- 1+match("missing",tolower(nexus_words));
			while(nexus_words[j]=="=")	j <- j+1;
			missing <- nexus_words[j];
			}
		}
	if (!is.na(match("symbols",tolower(nexus_words))))	{
		j <- match("symbols",tolower(nexus_words))+1;
		while (nexus_words[j] %in% c("=","\""))	j <- j+1;
		jj <- min(((j+1):length(nexus_words))[tolower(nexus_words)[((j+1):length(nexus_words))] %in% c("missing","gap",";")]-1);
#		jj <- j+match(";",nexus_words[(j+1):length(nexus_words)])-1;
		all_states <- gsub("\"","",nexus_words[j:jj]);
#		if (tolower(all_states) %in% "")
		}
	if (!is.na(match("fa",tolower(nexus_words))) || !is.na(match("fka",tolower(nexus_words))))	{
		nexus_words[tolower(nexus_words)=="fka"] <- "fa";
		strat <- as.numeric(nexus_words[match("fa",tolower(nexus_words))-1]);
		}
	if (!is.na(match("la",tolower(nexus_words))) || !is.na(match("lka",tolower(nexus_words))))	{
		nexus_words[tolower(nexus_words)=="lka"] <- "la";
		range <- as.numeric(nexus_words[match("la",tolower(nexus_words))-1]);
		}
	if (!is.na(match("geog",tolower(nexus_words))))	{
		geog <- as.numeric(nexus_words[match("geog",tolower(nexus_words))-1]);
		} else if (!is.na(match("geography",tolower(nexus_words))))	{
		geog <- as.numeric(nexus_words[match("geography",tolower(nexus_words))-1]);
		}
	}
if (is.null(all_states))	all_states <- 0:9;

extra <- 0;
if (strat>0)	{
	if (range>0)	{
		nchars <- nchars-2
		extra <- 2
		} else {
		nchars <- nchars-1
		extra <- 1
		}
	strat_ranges <- matrix(0,notu,2)
	}
if (geog>0)	{
	nchars <- nchars-1
	geography <- vector(length=notu)
	extra <- extra+1
	}

taxa <- vector(length=notu);
nstates <- array(0,dim=nchars);
chmatrix <- array(0,dim=c(notu,nchars));
tx <- 1;

# look for outgroup designation
exclude <- outgroup <- -1;
if (!is.na(match("BEGIN SETS;",nexus)))	{
	tx_pt <- match("BEGIN SETS;",nexus);	# look at taxon partitions
	look_for_outgroup <- TRUE;
	while (look_for_outgroup)	{
		tx_pt <- 1+tx_pt;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("-"," - ",yyy);
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub("  -  "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("outgroup",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("outgroup",tolower(xxx));
			while (xxx[jj]==":" || xxx[jj]=="=")	jj <- jj+1;
			outgroup <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					outgroup <- c(outgroup,((as.numeric(outgroup[length(outgroup)])+1):as.numeric(xxx[jj])));
					} else	{
					outgroup <- c(outgroup,xxx[jj]);
					}
				jj <- jj+1;
				}
			look_for_outgroup <- FALSE;
			} else	{
			if (tolower(nexus[tx_pt])=="end;" || tolower(nexus[tx_pt])=="\tend;")
				look_for_outgroup <- FALSE;
			}
		}

	# look for characters to exclude
	tx_pt <- match("BEGIN SETS;",nexus);
	xxx <- strsplit(paste(nexusfile[tx_pt-1,],collapse = "")," ");
	while(tolower(xxx[1])!="end")	{
		tx_pt <- tx_pt+1;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (length(xxx)==0 || is.na(xxx[1]))	xxx <- "";
#		if (!is.na(xxx) && !is.null(xxx) && xxx!="")	{
		if (xxx[1]=="charpartition")	{
			if (xxx[1]=="charpartition" && !is.na(match("exclude",tolower(xxx))))	{
				ttl_ln <- length(xxx);
				jj <- 1+match("exclude",tolower(xxx));
				while (xxx[jj]==":")	jj <- jj+1;
				exclude <- c();
				while (xxx[jj]!="," && jj<ttl_ln)	{
					if (xxx[jj]=="-")	{
						jj <- jj+1;
						exclude <- c(exclude,((as.numeric(exclude[length(exclude)])+1):as.numeric(xxx[jj])));
						} else	{
						exclude <- c(exclude,as.numeric(xxx[jj]));
						}
					jj <- jj+1;
					}
				}
			}
#		xxx[1];
#		tx_pt;
		}
	}

if (rate_partitions!="")	{
	ln <- (1:length(nexus))[gsub(rate_partitions,"",nexus)!=nexus];
	ln_min <- match("BEGIN SETS;",nexus);
	if (length(ln)>0)	{
		if (length(ln)>1)	ln <- min(ln[ln>ln_min])
		nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
		nexus[ln] <- gsub("  -"," -",nexus[ln]);
		nexus[ln] <- gsub("-  ","- ",nexus[ln]);
		nexus[ln] <- gsub("'","",nexus[ln]);
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		breakup_this_line <- gsub(",","",breakup_this_line);
		breakup_this_line <- gsub(";","",breakup_this_line);
		breakup_this_line <- breakup_this_line[breakup_this_line!=""];
		breakup_this_line <- breakup_this_line[match(rate_partitions,breakup_this_line):length(breakup_this_line)];
		kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
		partition_names <- breakup_this_line[kk-1];
		kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
		character_rate_partitions <- rep("",nchars);
		for (pn in 1:length(partition_names))	{
					ll <- kk[pn]+1;
					this_part <- as.numeric(breakup_this_line[ll]);
					ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
					if (pn < length(partition_names))	{
						break_cell <- kk[pn+1]-1;
						} else	{
						break_cell <- kk[pn+1];
						}
					while (ll<break_cell)	{
						if (breakup_this_line[ll]=="-")	{
							ll <- ll+1;
							this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
							} else	{
							this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
							}
						ll <- ll+1;
						}
					character_rate_partitions[this_part] <- partition_names[pn];
					}
		got_splits <- TRUE;
		} else	character_rate_partitions <- rep("imagine",nchars);
	} else	character_rate_partitions <- rep("imagine",nchars);

if (rate_partitions=="discard_this")	{
	ln <- (1:length(nexus))[gsub("Shell_Type_Partitions","",nexus)!=nexus];
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(rate_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(rate_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_rate_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_rate_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- TRUE;
			}
		}
	}

if (trend_partitions!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(trend_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(trend_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_trend_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_trend_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- TRUE;
			}
		}
	} else	character_trend_partitions <- rep("square",nchars);

state_orders <- rep("unordered",nchars);

if (!is.na(match("BEGIN ASSUMPTIONS;",nexus)))	{
	tx_pt <- 1+match("BEGIN ASSUMPTIONS;",nexus);	# look at taxon partitions
	while (tolower(nexus[tx_pt])!="end;")	{
#		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",nexus[tx_pt]);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("ord:",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("ord:",xxx);
			while (xxx[jj]==":")	jj <- jj+1;
			ordered <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					ordered <- c(ordered,((as.numeric(ordered[length(ordered)])+1):as.numeric(xxx[jj])));
					} else	{
					ordered <- c(ordered,as.numeric(xxx[jj]));
					}
				jj <- jj+1;
				}
			state_orders[ordered] <- "ordered";
			}
		tx_pt <- 1+tx_pt;
		}
	}
mxln <- length(nexusfile[top,]);
s <- top;
# te all of the taxon names
for (tx in 1:notu)	{
	# first, read taxon name
	#### look for quotations###
	s <- top+tx-1;
	endline <- match("\n",nexusfile[s,]);
	if (is.na(endline))	endline <- length(nexusfile[s,]);
	if (nexusfile[s,1]=="'" || nexusfile[s,2]=="'")	{
		jj <- ((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"]);
		i <- max((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"])
		taxa[tx] <- pracma::strcat(nexusfile[s,(jj[1]+1):(jj[2]-1)])
		i <- i+1
		while (nexusfile[s,i]==" " && i<ncol(nexusfile))	i <- i+1
		}	else	{
		i <- 1;
		if (nexusfile[s,1]!="\"")  {
			while (nexusfile[s,i]=="\t")	i <- i+1;
			taxa[tx] <- nexusfile[s,i]
			i <- i+1
			while (nexusfile[s,i]!=" " && nexusfile[s,i]!='\t' && i<ncol(nexusfile))	{
				if (nexusfile[s,i]!="_")	{
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
					} else {
					taxa[tx] <- paste0(taxa[tx]," ")
					}
				i <- i+1
				}
			}	else {
			taxa[tx] <- nexusfile[s,2];
			i <- 3;
#			while (nexusfile[s,i]!=" " && nexusfile[s,i+1]!=" " && i<ncol(nexusfile))	{
			while (nexusfile[s,i]!=" " && i<ncol(nexusfile))	{
				if (as.character(nexusfile[s,i])!="\"")
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
				i <- i+1;
				#print(taxa[tx]);
				}
			}
		# now, get to characters
		i <- (i:endline)[!nexusfile[s,i:endline] %in% c(" ","\t")][1];
#		while ((nexusfile[s,i]==" " || nexusfile[s,i]=="\t") && i<ncol(nexusfile))
#			i <- i+1
		}
	k <- i;
	if ((endline-k)==(nchars+extra))	{
		# true if there are no polymorphic characters for the taxon
		dummy <- nexusfile[s,k:(endline-1)];
		dummy[dummy==missing] <- UNKNOWN;
		dummy[dummy==gap] <- INAP;
		letterstate <- dummy[!dummy %in% c(UNKNOWN,INAP)];
		dummy[!dummy %in% c(UNKNOWN,INAP)] <- sapply(letterstate,switch_letter_state_to_numeric,all_states);
		chmatrix[tx,] <- as.numeric(dummy[1:nchars]);
		if (strat>0)	{
			strat_ranges[tx,1] <- strat_ranges[tx,2] <- as.numeric(dummy[strat])
			if (range>0)	strat_ranges[tx,2] <- as.numeric(dummy[range])
			}
		if (geog>0)	geography[tx]=as.numeric(nexusfile[geog,i])
		for (c in 1:nchars)	{
			if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
			}
		} else	{
#		for (c in 1:(nchars+extra))	{
		c <- 0;
		while (c < (nchars+extra))	{
			c <- c+1;
			#print(c);
			if (c<=nchars)	{
				if (nexusfile[s,i]=="(" || nexusfile[s,i]=="{")	{
					if (polymorphs==TRUE || polymorphs==1)	{
						# added 2020-11-28: sometimes polymorphics come in out-of-order
						riteparens <- (i:endline)[nexusfile[s,i:endline] %in% c(")","}")];
						ddd <- (i+1):(riteparens[1]-1);
						polysites <- ddd[!nexusfile[s,ddd] %in% c(",","&"," ")]
						polystates <- nexusfile[s,polysites];
						for (ps in 1:length(polystates))
							if (!polystates[ps] %in% 0:9)
								polystates[ps] <- switch_letter_state_to_numeric(polystates[ps],all_states=all_states);
						nexusfile[s,polysites] <- sort(as.numeric(polystates));
#						nexusfile[s,(i+1):(riteparens[1]-1)] <- sort(nexusfile[s,(i+1):(riteparens[1]-1)]);
						i <- i+1;
						w <- as.numeric(nexusfile[s,i])
						chmatrix[tx,c] <- -1*as.numeric(nexusfile[s,i])
						if ((1+w)>nstates[c])  nstates[c] <- 1+w;
						i <- i+1
						j <- 1
						while (nexusfile[s,i]!=")" && nexusfile[s,i]!="}" && i<ncol(nexusfile))	{
							if (nexusfile[s,i]!="," && nexusfile[s,i]!=" ")	{
								w <- as.numeric(nexusfile[s,i])
								if ((w+1)>nstates[c])	nstates[c] <- w+1
								chmatrix[tx,c] <- chmatrix[tx,c]-((10^j)*w)
								i <- i+1
								j <- j+1
								} else {
								i <- i+1
								}
							}
						}	else {
						chmatrix[tx,c] <- UNKNOWN;
						while (nexusfile[s,i]!=')' && nexusfile[s,i]!="}")	i <- i+1;
						}
					} else if (nexusfile[s,i]==missing)	{
					chmatrix[tx,c] <- UNKNOWN;
					}	else if (nexusfile[s,i]==gap)	{
					chmatrix[tx,c] <- INAP;
					} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")  {
					chmatrix[tx,c] <- switch_letter_state_to_numeric(nexusfile[s,i],all_states=all_states);
					}	else if (nexusfile[s,i]>="0" && nexusfile[s,i]<="9") {
					chmatrix[tx,c] <- as.numeric(nexusfile[s,i]);
					}
				if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1;
				if (i < (endline-1)) i <- i+1;
				}	else {
				if (c==strat)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,1]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,1]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					if (range==0)	strat_ranges[tx,2] <- strat_ranges[tx,1]
					i <- i+1
					} else if (c==range)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,2]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,2]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					i <- i+1
					} else if (c==geog)	{
						if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
							geography[tx]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
							geography[tx]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					}
				} # end non-morphological data
#			print(nexusfile[s,k:83]);
#			print(chmatrix[tx,])
			if (nexusfile[s,i+1]=="\n" || i==(mxln-1)) c <- nchars+extra;
			}
		}
	
	atomizer <- strsplit(taxa[tx],"")[[1]];
	while (atomizer[length(atomizer)]==" " && length(atomizer)>1)	atomizer <- atomizer[1:(length(atomizer)-1)];
	taxa[tx] <- paste(atomizer,collapse="");
	taxa[tx] <- gsub("__","_",taxa[tx]);
#	chmatrix[tx,];
#	tx <- tx+1;
#	s <- s+1
	}
	#x <- list(taxa,chmatrix,strat_ranges,geography)
	#return (list(taxa,chmatrix,strat_ranges,geography))

chmatrix <- mundify_character_matrix(chmatrix,minst=0,UNKNOWN,INAP);	# clean up coding
nstates <- count_states(chmatrix,UNKNOWN,INAP);

tree_found <- 0;
while (s<length(nexus) && tree_found==0)	{
	while (nexus[s]!= "BEGIN TREES; " && s<length(nexus))	s <- s+1;
	if (s<length(nexus))	{
		while (tree_found==0 && s<length(nexus))	{
			s <- s+1
			jj <- strsplit(nexus[s],split=c("\t"," "),fixed=TRUE)[[1]];
			jj <- paste(jj,collapse="")
			jj <- strsplit(jj,split=" ",fixed=TRUE)[[1]];
			if (sum(jj=="TREE")>0 || sum(jj=="tree")>0)	tree_found <- 1;
			}
#		s <- s+notu;
#		while (jj[i]=="")	jj[i] <- NULL;
#		while (j[1]=="\t")	j <- j[2:length(j)];
#		if (j[1]=="TRUE" && j[2]=="R" && j[3]=="E")	{
#			while (j!="(")	j <- j[2:length(j)];
#			}
		newick_string <- jj[length(jj)];
		newick_string <- fix_newick_ancestors(jj[length(jj)])
		tree <- read_newick_string(newick_string);
		tree_found <- 1
		s <- length(nexus);
		}
	}

row.names(chmatrix) <- taxa;

unscored_taxa <- c();
for (n in 1:notu)	{
	if (sum(chmatrix[n,]==UNKNOWN)==nchars)
		unscored_taxa <- c(unscored_taxa,n);
	}

if (nchars<10)	{
	colnames(chmatrix) <- 1:nchars;
	} else if (nchars<100)	{
	colnames(chmatrix) <- c(paste(0,(1:9),sep=""),10:nchars);
	} else if (nchars<1000)	{
	colnames(chmatrix) <- c(paste(00,(1:9),sep=""),paste(0,(10:99),sep=""),100:nchars);
	}
if (exclude[1]!=-1)	{
	keepers <- (1:nchars)[!(1:nchars) %in% exclude];
	chmatrix <- chmatrix[,keepers];
	nstates <- nstates[keepers];
	state_orders <- state_orders[keepers];
	character_rate_partitions <- character_rate_partitions[keepers];
	}

if (strat!=0 && geog!=0 && tree_found==1)  {
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,geography,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else if (strat!=0)  {
	if (geog!=0)  {
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,geography,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else if (tree_found!=0)	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (geog!=0)  {
	if (tree_found!=0)	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,geography,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,geography,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (tree_found!=0) {
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else	{
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	}

return(output)
}

#nexus_name <- "Anopliidae_Hanger_&_Strong_2000.nex";
accersi_matrix_location_from_RData <- function(matrix_name,character_database)	{
glossary <- character_database$Glossary;
nexus_names <- pbapply::pbsapply(glossary,accersi_nexus_name_in_RData);
return(match(matrix_name,nexus_names));
}

accersi_matrix_data_from_RData <- function(matrix_name,character_database)	{
nn <- accersi_matrix_location_from_RData(matrix_name,character_database);
pathway <- character_database$Glossary[nn];
pathways <- strsplit(pathway,"\\$")[[1]];
pp <- match(pathways[1],names(character_database));
cl <- match(pathways[2],names(character_database[[pp]]));
if (pathways[2]==matrix_name)	{
	nexus <- unlist(character_database[[pp]][cl]);
	return(nexus);
	} else	{
	od <- match(pathways[3],names(character_database[[pp]][[cl]]));
	if (pathways[3]==matrix_name)	{
		nexus <- unlist(character_database[[pp]][[cl]][od]);
		return(nexus);
		} else	{
		sf <- match(pathways[3],names(character_database[[pp]][[cl]][[od]]));
		nexus <- unlist(character_database[[pp]][[cl]][[od]][sf]);
		return(nexus);
		}
	}
}

accersi_nexus_name_in_RData <- function(glossary)	{
return(strsplit(glossary,"\\$")[[1]][length(strsplit(glossary,"\\$")[[1]])])
}

accersi_data_from_RData <- function(matrix_name, character_database, polymorphs=TRUE, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
nexus <- accersi_matrix_data_from_RData(matrix_name=matrix_name,character_database);
output <- accersi_data_from_nexus_vector(nexus=as.character(nexus), polymorphs=polymorphs, UNKNOWN=UNKNOWN, INAP=INAP, rate_partitions=rate_partitions, trend_partitions=trend_partitions);
return(output);
}

#### DEAL WITH TRICKY CHARACTERS ####
switch_letter_state_to_numeric <- function(state,all_states=c(0:9,LETTERS[!LETTERS %in% c("I","O")]))  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
# 2020-12-01: allows for i's and o's, but assumes that they are not there!
# -1 is for 0 to be zero
return(match(state,all_states)-1);
}

switch_letter_state_to_numeric_old <- function(state)  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
if (state > 9)	{
	state <- toupper(state)
	poss_letter_states <- toupper(letters[!letters %in% c("i","o")]);
	return(9+match(state,poss_letter_states));
	} else	{
	return(state);
	}
}

switch_numeric_state_to_letter <- function(state)  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
if (state > 9)	{
#	state <- toupper(state)
	poss_letter_states <- toupper(letters[!letters %in% c("i","o")]);
	return(poss_letter_states[state-9]);
	} else	{
	return(state);
	}
}

unravel_polymorph_badass <- function(poly,minst=0)	{
combo <- -1*poly;
state_test <- as.numeric(strsplit(x=as.character(combo),split="")[[1]])
if (state_test==sort(state_test,decreasing = TRUE) && length(unique(state_test))==length(state_test))	{
	sts <- 1+floor(log10(abs(combo)))
	polymorphics <- vector(length=sts)
	base <- 10^(sts-1)
	for (s in 1:sts)	{
		polymorphics[s] <- floor(abs(combo)/base)
		combo <- combo%%base
		base <- base/10
		}
	} else	{
	breakpt <- match(max(state_test),state_test);
	if (breakpt > 2)	{
		i <- 1;
		while (i < breakpt)	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else if (sum(state_test<minst)>0)	{
		# this should happen
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else if ((length(state_test) %% 2)==0)	{
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else	{
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		}
	}
return (polymorphics);
}

unravel_polymorph <- function(poly)	{
combo <- -1*poly
sts <- 1+floor(log10(abs(combo)))
polymorphics <- vector(length=sts)

base <- 10^(sts-1)
for (s in 1:sts)	{
	polymorphics[s] <- floor(abs(combo)/base)
	combo <- combo%%base
	base <- base/10
	}
return (polymorphics)
}

ravel_polymorph <- function(polystates)	{
polystates <- sort(polystates,decreasing = TRUE);
polym <- polystates[1];
for (st in 2:length(polystates))	polym <- (10*polym)+polystates[st]
return(-1*polym)
}

ravel_polymorph_for_file <- function(polystates)	{
polystates <- sort(polystates,decreasing = FALSE);
return(paste("(",paste(polystates,collapse=""),")",sep=""));
}

fix_states_unique_to_polymorphic_combinations <- function(chrmatrix,UNKNOWN=-11,INAP=-22)	{
nstates <- count_states(chrmatrix);
min_states <- minimum_state(chrmatrix);
max_states <- maximum_state(chrmatrix);
nchars <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
polychs <- c();
for (ch in 1:nchars)	{
	resolved_character <- chrmatrix[chrmatrix[,ch]>=0,ch];
	resolved_states <- sort(unique(resolved_character));
	if (length(resolved_states)<nstates[ch])	{
		unresolved_states <- (min_states[ch]:max_states[ch])[!(min_states[ch]:max_states[ch]) %in% resolved_states];
		poly_scores <- chrmatrix[,ch][(1:notu)[chrmatrix[,ch]<0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)]]
		poly_otus <- (1:notu)[chrmatrix[,ch]<0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)];
		potus <- length(poly_otus);
		for (po in 1:potus)	{
			if (sum(unresolved_states %in% unravel_polymorph(chrmatrix[poly_otus[po],ch]))==1)	{
				chrmatrix[poly_otus[po],ch] <- unresolved_states[unresolved_states %in% unravel_polymorph(chrmatrix[poly_otus[po],ch])];
				} else if (sum(unresolved_states %in% unravel_polymorph(chrmatrix[poly_otus[po],ch]))>1) {
				print(ch);
				}
			}
		}
	}
return(chrmatrix);
}

change_polymorphs_with_all_states_to_unknowns <- function(chrmatrix,UNKNOWN=-11,INAP=-22)	{
nstates <- count_states(chmatrix=chrmatrix,UNKNOWN=UNKNOWN,INAP=INAP);
min_states <- minimum_state(chrmatrix);
max_states <- maximum_state(chrmatrix);
nchars <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
examples <- c();
for (ch in 1:nchars)	{
	resolved_character <- chrmatrix[chrmatrix[,ch]>=0,ch];
	resolved_states <- sort(unique(resolved_character));
	poly_scores <- chrmatrix[,ch][(1:notu)[chrmatrix[,ch]<0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)]]
	poly_otus <- (1:notu)[chrmatrix[,ch]<0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)];
	potus <- length(poly_otus);
	po <- 0;
	while (po < potus)	{
		po <- 1+po;
		if (length(unravel_polymorph(chrmatrix[poly_otus[po],ch]))==nstates[ch])	{
			examples <- rbind(examples,c(poly_otus[po],ch));
			chrmatrix[poly_otus[po],ch] <- UNKNOWN;
			}
		}
	}
#chrmatrix[examples]
return(chrmatrix);
}

convert_binary_polymorphs_to_ordered_three_states <- function(chrmatrix,chrtypes,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
polymorphs <- find_polymorphic_characters(chrmatrix);
nstates <- count_states(chrmatrix);
binary <- (1:nch)[nstates==2];
binary_polymorphs <- intersect(polymorphs,binary);
b_p <- length(binary_polymorphs);
if (b_p > 0)	{
	nstates[binary_polymorphs] <- 3;
	chrtypes[binary_polymorphs] <- 1;
	for (pc in 1:b_p)	{
		chr <- binary_polymorphs[pc]
		chrmatrix[chrmatrix[,chr]==1,chr] <- 2;
		chrmatrix[chrmatrix[,chr]==-10,chr] <- 1;
		}
	}
output <- list(chrmatrix,chrtypes);
names(output) <- c("chrmatrix","chrtypes");
return(output);
}

convert_binary_polymorphs_to_unordered_three_states <- function(chrmatrix,chrtypes,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
polymorphs <- find_polymorphic_characters(chrmatrix);
nstates <- count_states(chrmatrix);
binary <- (1:nch)[nstates==2];
binary_polymorphs <- intersect(polymorphs,binary);
b_p <- length(binary_polymorphs);
if (b_p > 0)	{
	nstates[binary_polymorphs] <- 3;
	chrtypes[binary_polymorphs] <- 0;
	for (pc in 1:b_p)	{
		chr <- binary_polymorphs[pc]
		chrmatrix[chrmatrix[,chr]==1,chr] <- 2;
		chrmatrix[chrmatrix[,chr]==-10,chr] <- 1;
		}
	}
output <- list(chrmatrix,chrtypes);
names(output) <- c("chrmatrix","chrtypes");
return(output);
}

convert_binary_polymorphs_to_three_states <- function(chrmatrix,chrtypes,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
polymorphs <- find_polymorphic_characters(chrmatrix);
nstates <- count_states(chrmatrix);
binary <- (1:nch)[nstates==2];
binary_polymorphs <- intersect(polymorphs,binary);
b_p <- length(binary_polymorphs);
if (b_p > 0)	{
	nstates[binary_polymorphs] <- 3;
#	chrtypes[binary_polymorphs] <- 0;
	for (pc in 1:b_p)	{
		chr <- binary_polymorphs[pc]
		chrmatrix[chrmatrix[,chr]==1,chr] <- 2;
		chrmatrix[chrmatrix[,chr]==-10,chr] <- 1;
		}
	}
return(chrmatrix);
}

find_polymorphic_characters <- function(chrmatrix,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chrmatrix);
notu <- nrow(chrmatrix);
polymorphs <- c();
for (ch in 1:nch)
	#chrmatrix[,ch][chrmatrix[,ch] < 0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)]
	if (length((1:notu)[chrmatrix[,ch] < 0 & !chrmatrix[,ch] %in% c(UNKNOWN,INAP)])>0)
		polymorphs <- c(polymorphs,ch)
return(polymorphs);
}

#### SUMMARIZE CHARACTER DATA ####
# count taxa scored with something other than missing or inapplicable
count_scored_characters_per_otu <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
notu <- nrow(chmatrix);
scored <- vector(length=nchars);
#scored <- c();
for (s in 1:notu)
	scored[s] <- nchars - (sum(chmatrix[s,]==UNKNOWN)+sum(chmatrix[s,]==INAP));
#	scored <- c(scored,notu - (sum(chmatrix[s,]==UNKNOWN)+sum(chmatrix[s,]==INAP)));
return(scored);
}

# count missing and/or inapplicable per otu
count_scored_otus_per_character <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
if (is.matrix(chmatrix))	{
	nchars <- ncol(chmatrix);
	notu <- nrow(chmatrix);
	} else	{
	nchars <- 1;
	notu <- length(chmatrix);
	dummy <- array(0,dim=c(length(chmatrix),1));
	dummy[,1] <- chmatrix;
	chmatrix <- dummy;
	}
#scored <- vector(length=nchars)
scored <- c();
for (c in 1:nchars)
	scored <- c(scored,notu - (sum(chmatrix[,c]==UNKNOWN)+sum(chmatrix[,c]==INAP)));
return(scored);
}

# count missing and/or inapplicable per otu
count_scored_otus_per_character_state <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chmatrix);
#notu <- nrow(chmatrix);
#scored <- vector(length=nch)
scored <- array(0,dim=c(nch,max(chstates)))
for (c in 1:nch)	{
	for (st in 1:chstates[c])	{
		stt <- st-1
		scored[c,st] <- sum(chmatrix[,c]==stt)
		}
	}
return(scored);
}

# count characters with autapomorphic taxa
count_autapomorphic_characters <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
nchars <- ncol(chmatrix);
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- 0
for (c in 1:nchars)	{
	if(sum(notus_per_chstate[c,1:chstates[c]]==1)>0)	{
		autaps <- autaps+1;
		}
	}
return(autaps);
}

# count states coding only one taxon
count_autapomorphic_states <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
# count number of states that are autapomorphic
nchars <- ncol(chmatrix);
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- 0
for (ch in 1:nchars)
	autaps <- autaps+sum(notus_per_chstate[ch,1:chstates[ch]]==1);
return(autaps);
}

# routine to list all characters with at least one autapomorphic state
list_autapomorphic_characters <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
#scored <- vector(length=nch)
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- c()
for (ch in 1:nchars)	{
	if(sum(notus_per_chstate[ch,1:chstates[ch]]==1)>0)	{
		autaps <- c(autaps,ch);
		}
	}
return(autaps);
}

# routine to list all character states that are autapomorphic
list_autapomorphic_states <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
#scored <- vector(length=nch)
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- c()
for (c in 1:nchars)	{
	if(sum(notus_per_chstate[c,1:chstates[c]]==1)>0)	{
		autaps <- rbind(autaps,c(c,chstates[c]));
		}
	}
return(autaps);
}

# count the number of scorings that are polymorphic. There can be 1 per character per taxon
count_polymorphic_scorings <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
polymorphs <- vector(length=nchars);
for (ch in 1:nchars)	{
	char_states <- chmatrix[,ch];
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	polys <- char_states[char_states<0];
	polymorphs[ch] <- length(polys);
	}
return(polymorphs);
}

# count the number of characters with a polymorphic scoring
count_polymorphic_characters <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
polymorphs <- count_polymorphic_scorings(chmatrix,UNKNOWN,INAP);
return(sum(polymorphs>0));
}

# count the number of states per character that show polymorphism
count_polymorphic_states <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))
	chmatrix <- data.frame(ch=chmatrix);
nchars <- ncol(chmatrix);
polymorphic_states <- vector(length=nchars);
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	polys <- char_states[char_states<0];
	poly_states <- c();
	pp <- 0;
	while (pp < length(polys))	{
		pp <- pp+1;
		poly_states <- unique(c(poly_states,unravel_polymorph_badass(polys[pp])));
		}
	polymorphic_states[ch] <- length(polys);
	}
return(polymorphic_states);
}

# get number of states for each character
count_states <- function(chmatrix,UNKNOWN=-11,INAP=-22,include_polymorphs=TRUE)	{
  if (!is.matrix(chmatrix))	{
    chmatrix <- data.frame(ch=chmatrix);
  }
  nchars <- ncol(chmatrix);
  nstates <- c();
  for (ch in 1:nchars)	{
    char_states <- sort(unique(chmatrix[,ch]));
    char_states <- char_states[char_states!=UNKNOWN];
    char_states <- char_states[char_states!=INAP];
    unique_states <- char_states[char_states>=0]
    if (sum(char_states<0)>0 && include_polymorphs)	{
      polys <- char_states[char_states<0];
      polystates <- c();
      for (pp in polys)	{
        polystates <- c(polystates,unravel_polymorph(poly=pp));
      }
      unique_states <- sort(unique(c(unique_states,polystates)));
    }
    nstates <- c(nstates,length(unique_states));
  }	# pick up here!!!
  return(nstates);
}

# get number of states for each character
list_states <- function(chmatrix,UNKNOWN=-11,INAP=-22,include_polymorphs=TRUE)	{
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
state_list <- list();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	unique_states <- char_states[char_states>=0]
	if (sum(char_states<0)>0 && include_polymorphs)	{
		polys <- char_states[char_states<0];
		polystates <- c();
		for (pp in polys)	{
			polystates <- c(polystates,unravel_polymorph(poly=pp));
			}
		unique_states <- sort(unique(c(unique_states,polystates)));
		}
	state_list <- rlist::list.append(state_list,sort(unique_states));
	}	# pick up here!!!
return(state_list);
}

count_states_old <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
nstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	nstates <- c(nstates,length(char_states));
	}
return(nstates);
}

maximum_state <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
maxstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (length(char_states)==0)	{
		char_states <- 0;
		}	else if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	maxstates <- c(maxstates,max(char_states));
	}
return(maxstates);
}

minimum_state <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
minstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (length(char_states)==0)	{
		char_states <- 0;
		} else if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	minstates <- c(minstates,min(char_states));
	}
return(minstates);
}

#### FIND & RECODE DEPENDENT - INDEPENDENT SETS ####

# modified 2022-02-05
find_independent_character <- function(dchar,independents,chmatrix,UNKNOWN,INAP)	{
dstates <- sum(unique(chmatrix[,dchar])>=0)
pre_indies <- independents[independents < dchar];
pi <- p_inds <- length(pre_indies);
ind_char <- 0;
for (pi in p_inds:1)	{
	ic <- pre_indies[pi];
	unique_combos <- unique(chmatrix[,c(ic,dchar)]);
	unique_combos_2 <- unique_combos[!unique_combos[,1] %in% c(INAP,UNKNOWN),];
	# exclude -- or ?1 pairings
	u_c <- nrow(unique_combos);
	u_c2 <- nrow(unique_combos_2);
	if (is.null(u_c2))	if (length(unique_combos_2)>0)	{u_c2 <- 1;} else {u_c2 <- 0;} 
	missing_match <- (1:u_c)[unique_combos[,1] %in% UNKNOWN][!(1:u_c)[unique_combos[,1] %in% UNKNOWN] %in% (1:u_c)[unique_combos[,2] %in% UNKNOWN]];
	inap_match <- (1:u_c)[unique_combos[,1] %in% INAP][!(1:u_c)[unique_combos[,1] %in% INAP] %in% (1:u_c)[unique_combos[,2] %in% c(UNKNOWN,INAP)]]
#	sum((unique_combos[,1] %in% INAP) * (unique_combos[,2] %in% INAP));
#	double_inap <- double_inap + sum((unique_combos[,1] %in% UNKNOWN) * (!unique_combos[,2] %in% c(UNKNOWN,INAP)));
	if (u_c2>=dstates && (length(missing_match)>0 || length(inap_match)>0))	{
		scored_matches <- unique(unique_combos_2[,1])[!unique(unique_combos_2[,1]) %in% unique(unique_combos_2[!unique_combos_2[,2] %in% c(INAP,UNKNOWN),1])];
		if (length(scored_matches)==1 && sum(scored_matches >= 0))	{
			ind_char <- ic;
			pi <- 1;
			return(ind_char);
			}
		} else if (u_c2>1)	{
		unique_combos <- unique_combos_2;
		if (!is.matrix(unique_combos))	unique_combos <- array(unique_combos,dim=c(1,2));
		unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
		if (!is.matrix(unique_combos))	unique_combos <- array(unique_combos,dim=c(1,2));
		poly_combos <- subset(unique_combos,unique_combos[,1]<0);
		unique_combos <- unique_combos[unique_combos[,1] >= 0,];
		unique_combos_ap <- unique_combos[!unique_combos[,2] %in% INAP,];
		unique_combos_inap <- unique_combos[unique_combos[,2] %in% INAP,];
		if (!is.matrix(unique_combos_ap))	{
			if (length(unique_combos_ap)==2) {
				unique_combos_ap <- array(unique_combos_ap,dim=c(1,2));
				} else	{
				unique_combos_ap <- cbind(0:5,0:5);
				}
			}
		if (!is.matrix(unique_combos_inap))	{
			if (length(unique_combos_inap)==2) {
				unique_combos_inap <- array(unique_combos_inap,dim=c(1,2));
				} else	{
				unique_combos_inap <- cbind(0:5,0:5);
				}
			}
		if (nrow(unique_combos_ap)>0 && nrow(unique_combos_inap)>0)	{
			if (sum(unique_combos_ap[,1] %in% unique_combos_inap[,1])==0)	{
				ind_char <- ic;
				pi <- 1;
				return(ind_char);
				}
			}	# end case where ic might be the "parent" character
		} else if (unique_combos[1]>=0 && unique_combos[2]==INAP)	{
		ind_char <- ic;
		pi <- 1;
		return(ind_char);
		}
	}
return(ind_char);
}

#transmogrify_additive_dependents_to_multistate <- function(ind_char,dep_chars,chmatrix,INAP,UNKNOWN,multichanges=F,theoretical=F,unknown_inap_sep=F)	{
#transmogrify_additive_dependents_to_multistate           (ind_char,dep_chars,chmatrix,char_dependencies,INAP,UNKNOWN,theoretical=TRUE)
transmogrify_additive_dependents_to_multistate <- function(ind_char,dep_chars,chmatrix,char_dependencies,INAP=-22,UNKNOWN=-11)	{
# theortical: if TRUE, then Q matrix & recoding allows for all theoretically possible combinations
# unknown_inap_sep: code unknowns separately; this is a bad idea that I regret and that we'll cut
notu <- nrow(chmatrix);
all_chars <- c(ind_char,dep_chars);
if (is.null(names(char_dependencies)) || is.na(names(char_dependencies[1])))
	names(char_dependencies) <- all_chars;
#indies <- all_chars[match(unique(char_dependencies),all_chars)];
combos <- chmatrix[,all_chars];
#combos[1,2] <- -10;
combos <- combos[!(rowMaxs(combos)==UNKNOWN & rowMins(combos)==UNKNOWN),];
hchars <- ncol(combos);
hstates <- count_states(combos);
names(hstates) <- all_chars;
#combos <- combos[combos[,2]!=INAP,];
ucombos <- unique(combos);
polyms <- unique(which(ucombos < 0 & !ucombos %in% INAP,arr.ind = TRUE)[,1]);
ucombos <- ucombos[!(1:nrow(ucombos)) %in% polyms,];
all_poss_combos <- accersi_all_theoretical_character_state_combinations(obs_combinations=ucombos,char_dependencies = char_dependencies,INAP=INAP);
a_p_c <- nrow(all_poss_combos);
#combo_dependency <- accersi_key_character_for_hierarchical_combination(all_poss_combos,char_dependencies);
#keystates <- accersi_key_states_for_independent_character_in_hierarchical(all_poss_combos = all_poss_combos,char_dependencies = char_dependencies)
key_combos <- accersi_gateway_states_among_hierarchical_characters(all_poss_combos,char_dependencies,hstates=hstates)
ind_states <- sort(unique(ucombos[,1]));
if (ncol(key_combos)>1)	{
	Q <- construct_Q_matrix_for_divergent_hierarchical_multistate(all_poss_combos,char_dependencies,UNKNOWN=UNKNOWN,INAP=INAP);
	# start here!!!! (why (-21)?)
	new_multistate <- accersi_hierarchical_multistate(all_poss_combos,obs_combos=chmatrix[,all_chars]);
	} else if (length(all_poss_combos)>0 & ncol(key_combos)==1)	{
	# get rid of missing examples
	Q <- construct_Q_matrix_for_hierarchical_multistate(all_poss_combos,char_dependencies,UNKNOWN=UNKNOWN,INAP=INAP);
#	Q <- construct_Q_matrix_for_divergent_hierarchical_multistate(all_poss_combos)
	new_multistate <- accersi_hierarchical_multistate(all_poss_combos,obs_combos=chmatrix[,all_chars]);
	#cbind(new_multistate,chmatrix[,c(ind_char,dep_chars)])
	} else	{
	# only one dependent state exists, making this just a normal character
	combos <- chmatrix[,c(ind_char,dep_chars)];
	combos <- combos[!(rowMaxs(combos)==UNKNOWN & rowMins(combos)==UNKNOWN),];
	combos <- unique(combos);
	all_poss_combos <- combos;
	new_multistate <- chmatrix[,ind_char];
	k <- max(2,nrow(combos));
	Q <- construct_Q_matrix_unordered(k);
	}
output <- list(all_poss_combos,Q,new_multistate);
names(output) <- c("unique_combinations","Q","new_character");
return(output);
}
#obs_combinations <- ucombos; char_dependencies <- secondary_dependencies;

accersi_all_theoretical_character_state_combinations <- function(obs_combinations,char_dependencies=c(),INAP=-22,UNKNOWN=-11)	{
ucombos <- unique(obs_combinations);
ochars <- ncol(ucombos);				# no. original characters involved
all_chars <- as.numeric(colnames(ucombos));
	#combos <- combos[order(combos[,1],combos[,2],combos[,3]),];
all_poss_combos <- c();
for (nch in ncol(ucombos):1)	{
	if (is.null(all_poss_combos))	{
		all_poss_combos <- array(sort(unique(ucombos[,nch])),dim=c(length(unique(ucombos[,nch])),1));
		} else	{
		new_states <- sort(unique(ucombos[,nch]));
		prior_combos <- all_poss_combos;
		all_poss_combos <- cbind(rep(new_states[1],nrow(prior_combos)),prior_combos);
		for (ns in 2:length(new_states))
			all_poss_combos <- rbind(all_poss_combos,cbind(rep(new_states[ns],nrow(prior_combos)),prior_combos));
		}
	}
colnames(all_poss_combos) <- colnames(ucombos);
state_combos <- all_poss_combos[,1];
for (dc in 2:ncol(all_poss_combos))	{
	dummyc <- as.character(all_poss_combos[,dc]);
	dummyc <- gsub("-22","-",dummyc);
	for (sc in 1:length(state_combos))	state_combos[sc] <- paste(state_combos[sc],dummyc[sc],sep="");
	}
rownames(all_poss_combos) <- state_combos;
dependent_chars <- c();
if (length(char_dependencies)>0)	{
	ind_char_2 <- match(ind_char,all_chars);
	for (nch in 1:ochars)	if (sum(ucombos[,nch] %in% INAP)>0)	dependent_chars <- c(dependent_chars,nch);
	for (dc in 1:length(dependent_chars))	{
		a_p_c <- nrow(all_poss_combos);
		dch <- dependent_chars[dc];
		ich <- match(char_dependencies[dependent_chars[dc]],all_chars)
		keyst <- unique(ucombos[!ucombos[,c(ich,dch)][,2] %in% INAP,ich]);
		keyst <- keyst[keyst != INAP];
		keepers1 <- (1:a_p_c)[all_poss_combos[,ich] %in% keyst & !all_poss_combos[,dch] %in% INAP];
		keepers2 <- (1:a_p_c)[!all_poss_combos[,ich] %in% keyst & all_poss_combos[,dch] %in% INAP];
		keepers <- sort(c(keepers1,keepers2));
		all_poss_combos <- all_poss_combos[keepers,];
#		if (ich!=ind_char_2)	{
#			keepers <- keepers[!all_poss_combos[,dch] %in% INAP];
#			all_poss_combos <- all_poss_combos[keepers,];
#			}
		}
	}
return(all_poss_combos);
}

accersi_hierarchical_multistate <- function(all_poss_combos,obs_combos)	{
# condense hierarchical additive scheme into one multistate
new_multistate <- vector(length=notu);
for (nn in 1:notu)	new_multistate[nn] <- match_vector_to_matrix_row(test_vector=obs_combos[nn,],test_matrix=all_poss_combos)-1;
all_states <- c(0:9,letter_states,more_letter_states);
if (nrow(all_poss_combos)>10)	if (is.numeric(new_multistate))	new_multistate[!is.na(new_multistate)] <- all_states[1+new_multistate[!is.na(new_multistate)]];
hstates <- count_states(all_poss_combos);
# look for uncoded examples; these will be either polymorphs or unknowns
prob_child <- (1:notu)[is.na(new_multistate)];
#combos <- chmatrix[,c(ind_char,dep_chars)];
pc <- 0;
polymorphs <- unique(obs_combos[obs_combos<0])
polymorphs <- polymorphs[!polymorphs %in% c(UNKNOWN,INAP)];
while (pc < length(prob_child))	{
	pc <- pc+1;
	pcc <- prob_child[pc];
	#obs_combos[pcc,]
	if (obs_combos[pcc,1]==UNKNOWN)	{
		new_multistate[pcc] <- "?";
#		new_multistate[pcc] <- UNKNOWN;
		} else	{
				# make it polymorphic for all possible states
#		this_combo <- obs_combos[pcc,(0:dchars)+1];
		this_combo <- obs_combos[pcc,];
	#		doofi <- (1:length(this_combo))[this_combo==UNKNOWN];
			# this if/else probably is unneeded now!
		if (sum(this_combo %in% polymorphs)==0)	{
			# must be case where 1+ dependent is unknown
			set_chars <- (1:length(this_combo))[this_combo!=UNKNOWN];
			set_states <- this_combo[set_chars];
			poss_combos <- all_poss_combos;
			ss <- 0;
			while (ss < length(set_chars))	{
				ss <- ss+1;
				poss_combos <- subset(poss_combos,poss_combos[,set_chars[ss]]==set_states[ss]);
				}
			polymorph <- 0;
			if (nrow(all_poss_combos)>10)	polymorph <- c();
			for (cp in 1:nrow(poss_combos))	{
				this_state <- (row.match(poss_combos[cp,],as.data.frame(all_poss_combos))-1);
				if (nrow(all_poss_combos)>10)	{
					polymorph <- c(polymorph,all_states[this_state]);
	#				base64encode(this_state);
	#					for (i in 1:5)	print(base64encode(10^(i-1)));
					} else	{
					polymorph <- polymorph-((10^(cp-1))*this_state);
					}
				}
			} else	{
			set_chars <- (1:ochars)[!this_combo %in% c(UNKNOWN,polymorphs)];				unset_chars <- (1:ochars)[this_combo %in% c(UNKNOWN,polymorphs)];
			un <- 0;
			missing <- unset_chars[this_combo[unset_chars]==UNKNOWN];
			polys <- unset_chars[!this_combo[unset_chars] %in% UNKNOWN];
			while (un < length(missing))	{
				un <- un+1;
				mchar <- missing[un];
				rstates <- (1:hstates[mchar])-1;
				this_combo[mchar] <- 0;
				for (rs in 1:length(rstates))
					this_combo[mchar] <- this_combo[rs]-rstates[rs]*(10^(rs-1));
				}
	#			set_chars <- set_chars[set_chars>1];
			set_states <- this_combo[set_chars];
			poss_combos <- all_poss_combos;
	#		print(nrow(poss_combos));
			# reduce the possible combinations to those consistent
			sc <- length(set_chars);
	#		for (sc in 1:length(set_chars))	{
			while (sc > 1)	{
				poss_combos <- subset(poss_combos,poss_combos[,set_chars[sc]]==set_states[sc]);
	#			print(nrow(poss_combos));
				sc <- sc-1;
				}
			uc <- 0;
			while (uc < length(polys))	{
				uc <- uc+1;
				poss_combos <- poss_combos[poss_combos[,polys[uc]] %in% unravel_polymorph_badass(this_combo[polys[uc]]),];
				}
			pstates <- match(rownames(poss_combos),rownames(all_poss_combos));
			polymorph <- 0;
			if (nrow(all_poss_combos)>10)	polymorph <- c();
			for (cp in 1:nrow(poss_combos))	{
				this_state <- (row.match(poss_combos[cp,],as.data.frame(all_poss_combos))-1);
				if (nrow(all_poss_combos)>10)	{
					polymorph <- c(polymorph,all_states[this_state]);
	#				base64encode(this_state);
	#				for (i in 1:5)	print(base64encode(10^(i-1)));
					} else	{
					polymorph <- polymorph-((10^(cp-1))*this_state);
					}
				}
			}
	#	if (nrow(unq_combos)>10)	{
	#		if (is.numeric(new_multistate))
	#			new_multistate[!is.na(new_multistate)] <- all_states[1+new_multistate[!is.na(new_multistate)]]
		new_multistate[pcc] <- paste("{",paste(polymorph,collapse=""),"}",sep="");
		new_multistate[pcc] <- gsub("-","",new_multistate[pcc]);
	#			} else	{
	#			new_multistate[pcc] <- polymorph;
	#			}
		}
	}
return(new_multistate);
}

accersi_key_states_for_independent_character_in_hierarchical <- function(all_poss_combos,char_dependencies,INAP=-22,UNKNOWN=-11)	{
all_chars <- as.numeric(colnames(all_poss_combos));
indies <- all_chars[match(unique(char_dependencies),all_chars)];
dep_chars <- as.numeric(colnames(all_poss_combos))[as.numeric(colnames(all_poss_combos))!=char_dependencies];
ostates <- count_states(all_poss_combos);
keystates <- array(-1,dim=c(length(dep_chars),max(ostates[match(indies,all_chars)])));
#unkeystates <- sort(unique(ucombos[,1]));
for (dc in 1:length(dep_chars))	{
	kc <- match(char_dependencies[dc+1],all_chars)
	ks <- unique(combos[!combos[,dc+1] %in% c(INAP,UNKNOWN),kc]);
	ks <- ks[ks!=INAP];
	keystates[dc,1:length(ks)] <- ks;
#	unkeystates <- unkeystates[unkeystates!=ks];
	}
rownames(keystates) <- dep_chars; # rownames(keystates) <- c(dep_chars,49);
return(keystates);
}

accersi_key_character_for_hierarchical_combination <- function(all_poss_combos,char_dependencies,INAP=-22,UNKNOWN=-11)	{
all_chars <- as.numeric(colnames(all_poss_combos));
ochars <- ncol(all_poss_combos);
hchars <- nrow(all_poss_combos);
combo_dependency <- vector(length=hchars);
for (hc in 1:hchars)	combo_dependency[hc] <- max(0,char_dependencies[(1:ochars)[!all_poss_combos[hc,] %in% INAP]]);
return(combo_dependency);
}

accersi_gateway_states_among_hierarchical_characters <- function(all_poss_combos,char_dependencies,hstates)	{
all_chars <- as.numeric(colnames(all_poss_combos));
keystates <- accersi_key_states_for_independent_character_in_hierarchical(all_poss_combos = all_poss_combos,char_dependencies = char_dependencies)
key_combos <- array(-1,dim=c(length(unique(char_dependencies)),max(hstates[match(unique(char_dependencies),all_chars)])))
rownames(key_combos) <- unique(char_dependencies);
for (i in 1:nrow(keystates))	{
	j <- match(as.numeric(rownames(keystates))[i],all_chars);	# get dependency
	k <- match(char_dependencies[j],rownames(key_combos));		# match dependency
	h <- keystates[i,keystates[i,]>=0];
	h <- h[!h %in% key_combos[k,]];
	if (length(h)>0)	{
		g <- 1+sum(key_combos[k,]>-1);
		g <- g:(g+length(h)-1);
		key_combos[k,g] <- h;
		}
	}
mx_states <- as.array(colMaxs(key_combos));
new_dim <- orig_dim <- dim(key_combos);
new_dim[2] <- sum(mx_states>=0);
if (min(new_dim)==1)	{
	orig_rownames <- rownames(key_combos);
	key_combos <- array(key_combos[,1:new_dim[2]],dim=c(new_dim[1],new_dim[2]));
	rownames(key_combos) <- orig_rownames;
	colnames(key_combos) <- letters[1:new_dim[2]];
	}
return(key_combos);
}

#### MODIFY & WRITE NEXUS FILES ####
# routine to "clean" character matrix (e.g., remove gaps in coding, standarize minimum states, etc.)
mundify_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix);	# replaces spc to standardize coding.
nchars <- ncol(chmatrix);
min_states <- minimum_state(chmatrix,UNKNOWN=UNKNOWN,INAP=INAP);
max_states <- maximum_state(chmatrix,UNKNOWN=UNKNOWN,INAP=INAP);
for (ch in 1:nchars)	{
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP]);
	if (length(rem)>0)	{
		test <- chmatrix[-rem,ch]
		}	else test <- chmatrix[,ch];
	# check polymorphics for anything that needs to be changed
	if (length(rem) < notu)	{
		polys <- sum(test<0);	# taxa with polymorphic scores
		coded <- sort(unique(test[test>=0]));
		if (polys>0)	{
#			examps <- test[test<0]
			polycoded <- sort(unique(test[test<0]))
			for (i in 1:length(polycoded))	{
				polystates <- unravel_polymorph(polycoded[i])
				coded <- sort(unique(c(coded,polystates)))
#				if (min(polystates)<minstch)	minstch <- min(polystates)
				}
			} else	{
			polycoded <- c();
			}
		minstch <- min(coded);
	# eliminate gaps in states
		if (sum(!min(coded):max(coded) %in% coded)>0)	{
#			new_codes <- match(coded,coded)-(1-minst);
			new_codes <- match(coded,coded)-(1-minstch);
			for (st in 1:length(coded))	{
				if (coded[st]!=new_codes[st])	{
					rec <- (1:notu)[chmatrix[,ch]==coded[st]];
					chmatrix[rec,ch] <- new_codes[st];
					redo_poly <- 1;
					while (redo_poly <= length(polycoded))	{
						polystates <- unravel_polymorph(polycoded[redo_poly]);
						polystates[polystates==coded[st]] <- new_codes[st]
						newpolystates <- ravel_polymorph(polystates);
						testp <- (1:notu)[chmatrix[,ch]==polycoded[redo_poly]];
#						if (newpolystates != polycoded[redo_poly])	{
						chmatrix[testp,ch] <- newpolystates;
#							}
#						polycoded[redo_poly] <- chmatrix[,ch][chmatrix[,ch] %in% polycoded[redo_poly]] <- newpolystates
						polycoded[redo_poly] <- newpolystates;
						redo_poly <- redo_poly+1;
						}
					coded[st] <- new_codes[st];
					}
				}
			}
		# standardize minimum state
		# simple cheat: subtract 1111 to polymorphics
		if (minstch!=minst)	{
			adj <- minst-minstch;
			test2 <- (1:notu)[chmatrix[,ch]>=0];
			chmatrix[test2,ch] <- chmatrix[test2,ch]+adj;
			if (polys>0)	{
				examps2 <- polycoded;
				for (i in 1:length(polycoded))	{
					testp <- (1:notu)[chmatrix[,ch]==polycoded[i]];
					examps2[i] <- examps2[i]-(adj*floor(10^floor(log10(abs(polycoded[i])))*10/9));
					chmatrix[testp,ch] <- examps2[i];
					}
				polycoded <- examps2;
				} # end rescoring of polytomies
			} # end rescaling stats
		} # end case where rem < notu
	} # end search of characters;
return(chmatrix);
}

# routine to remove invariant and/or unscored characters from matrix
remove_invariants_from_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix)	# replaces spc to standardize coding.
ncharss <- ncol(chmatrix)
rem_char <- c()
for (ch in 1:ncharss)	{
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP])
	if (length(rem)>0)	{
		test <- chmatrix[-rem,ch]
		} else	test <- chmatrix[,ch]
	if (length(unique(test))<2)	rem_char <- c(rem_char,ch)
	}
return(chmatrix[,-rem_char])
}

# generate composite score from 2+ scored taxa
accersi_composite_scores <- function(mini_matrix,return_polymorph=TRUE,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(mini_matrix);
nchars <- ncol(mini_matrix);
composite_score <- c();
ch <- 1;
for (ch in 1:nchars)	{
	composite_states <- unique(mini_matrix[,ch]);
	if (length(composite_states)==1)	{
		composite_score <- c(composite_score,composite_states);
		} else	{
		composite_states <- composite_states[composite_states!=INAP];
		if (length(composite_states)>1)	{
			composite_states <- composite_states[composite_states!=UNKNOWN];
			if (length(composite_states)>1)	{
				polyscored <- composite_states[composite_states<0];
				if (length(polyscored)>0)	{
					accersi_states <- sapply(polyscored,unravel_polymorph);
					composite_states <- sort(unique(c(accersi_states,composite_states[composite_states>=0])));
					}
				composite_score <- c(composite_score,ravel_polymorph(composite_states));
				} else	{
				composite_score <- c(composite_score,composite_states);
				}
			} else	{
			composite_score <- c(composite_score,composite_states);
			}
		}
#	ch <- ch+1;
#	composite_score;
	}
#rbind(mini_matrix,composite_score)
return(composite_score);
}

scribio_nexus_file_from_chmatrix <- function(ch_matrix,new_file_name,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(ch_matrix);
taxon_names <- rownames(ch_matrix);
nchars <- ncol(ch_matrix);
nstates <- count_states(chmatrix = ch_matrix);

nexus_file_content <- c();
nexus_file_content <- rbind("#NEXUS","","BEGIN DATA;")
nexus_file_content <- rbind(nexus_file_content,paste("	DIMENSIONS  NTAX=",notu," NCHAR=",nchars,";",sep=""));
if (max(nstates)<10) {
	state_symbols <- " ";
	for (st in 1:max(nstates))
		state_symbols <- paste(state_symbols,st-1,sep=" ");
	} else	{
	mxl <- max(nstates)-10;
	letter_states <- LETTERS[!LETTERS %in% c("I","O")][1:mxl]
	all_states <- c(0:9,letter_states);
	state_symbols <- paste(all_states,collapse=" ");
	}
nexus_file_content <- rbind(nexus_file_content,paste("	FORMAT DATATYPE = STANDARD RESPECTCASE GAP = - MISSING = ? SYMBOLS = \"",state_symbols,"\";"));
nexus_file_content <- rbind(nexus_file_content,"	MATRIX");

string_to_count <- taxon_names;
name_lengths <- sapply(string_to_count,count_characters_in_string);
max_name_length <- max(name_lengths);
need_quotes <- c(".","(",")","[","]");
for (nn in 1:notu)	{
	test_name <- strsplit(taxon_names[nn],split="",fixed=TRUE)[[1]]
	if (sum(test_name %in% need_quotes)==0)	{
		taxon <- gsub(" ","_",taxon_names[nn]);
		} else	{
		taxon <- paste("\"",taxon_names[nn],"\"",sep="");
		name_lengths[nn] <- name_lengths[nn]+2;
		}
	this_line <- paste("\t",taxon,paste(rep(" ",(5+(max_name_length-name_lengths[nn]))),collapse=""),sep="");
	otu_code <- c();
	for (ch in 1:nchars)	{
		if (ch_matrix[nn,ch]>=0 && ch_matrix[nn,ch]<=9)	{
			otu_code <- paste(otu_code,ch_matrix[nn,ch],sep="");
			} else if (ch_matrix[nn,ch]>9)	{
			otu_code <- paste(otu_code,all_states[1+ch_matrix[nn,ch]],sep="");	# note: we need +1 because of state 0
			} else if (ch_matrix[nn,ch]==UNKNOWN)	{
			otu_code <- paste(otu_code,"?",sep="");
			} else if (ch_matrix[nn,ch]==INAP)	{
			otu_code <- paste(otu_code,"-",sep="");
			} else if (ch_matrix[nn,ch]<0)	{
			polystates <- strsplit(as.character(ch_matrix[nn,ch]),split="",fixed=TRUE)[[1]];
			polystates <- as.numeric(polystates[polystates!="-"]);
			otu_code <- paste(otu_code,ravel_polymorph_for_file(polystates),sep="");
			}
		}
	nexus_file_content <- rbind(nexus_file_content,paste(this_line,otu_code,sep=""));
	}

nexus_file_content <- rbind(nexus_file_content,";");
nexus_file_content <- rbind(nexus_file_content,"END;");
nexus_file_content <- rbind(nexus_file_content,"begin mrbayes;");
nexus_file_content <- rbind(nexus_file_content,"	set autoclose=yes nowarn=yes;");
nexus_file_content <- rbind(nexus_file_content,"	lset nst=6 rates=invgamma;");
nexus_file_content <- rbind(nexus_file_content,"	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); ");
nexus_file_content <- rbind(nexus_file_content,"	prset applyto=(all) ratepr=variable;");
nexus_file_content <- rbind(nexus_file_content,"	mcmcp ngen= 100000000 relburnin=yes burninfrac=0.25 printfreq=10000  samplefreq=10000 nchains=4 savebrlens=yes;");
nexus_file_content <- rbind(nexus_file_content,"	mcmc;");
nexus_file_content <- rbind(nexus_file_content,"	sumt;");
nexus_file_content <- rbind(nexus_file_content,"end;");
write(nexus_file_content,file=new_file_name);
}

# write nexus file from chmatrix that already is converted to character
scribio_nexus_file_from_chmatrix_character <- function(ch_matrix_ch,new_file_name,max_states,unknown="?",inap="-")	{
notu <- nrow(ch_matrix_ch);
taxon_names <- rownames(ch_matrix_ch);
nchars <- ncol(ch_matrix_ch);

nexus_file_content <- c();
nexus_file_content <- rbind("#NEXUS","","BEGIN DATA;")
nexus_file_content <- rbind(nexus_file_content,paste("	DIMENSIONS  NTAX=",notu," NCHAR=",nchars,";",sep=""));
if (max_states<10) {
	state_symbols <- " ";
	for (st in 1:max_states)
		state_symbols <- paste(state_symbols,st-1,sep=" ");
	} else	{
#	mxl <- max_states-10;
	all_states <- c(0:9,letter_states,more_letter_states);
	state_symbols <- paste(all_states[1:max_states],collapse=" ");
	}
nexus_file_content <- rbind(nexus_file_content,paste("	FORMAT DATATYPE = STANDARD RESPECTCASE GAP = ",inap," MISSING = ",unknown," SYMBOLS = \"",state_symbols,"\";"));
nexus_file_content <- rbind(nexus_file_content,"	MATRIX");

string_to_count <- taxon_names;
name_lengths <- sapply(string_to_count,count_characters_in_string);
max_name_length <- max(name_lengths);
need_quotes <- c(".","(",")","[","]");
for (nn in 1:notu)	{
	test_name <- strsplit(taxon_names[nn],split="",fixed=TRUE)[[1]]
	if (sum(test_name %in% need_quotes)==0)	{
		taxon <- gsub(" ","_",taxon_names[nn]);
		} else	{
		taxon <- paste("\"",taxon_names[nn],"\"",sep="");
		name_lengths[nn] <- name_lengths[nn]+2;
		}
	this_line <- paste("\t",taxon,paste(rep(" ",(5+(max_name_length-name_lengths[nn]))),collapse=""),sep="");
	this_line <- paste(this_line,paste(ch_matrix_ch[nn,],collapse=""),sep="");
	nexus_file_content <- rbind(nexus_file_content,this_line);
	}

nexus_file_content <- rbind(nexus_file_content,";");
nexus_file_content <- rbind(nexus_file_content,"END;");
nexus_file_content <- rbind(nexus_file_content,"begin mrbayes;");
nexus_file_content <- rbind(nexus_file_content,"	set autoclose=yes nowarn=yes;");
nexus_file_content <- rbind(nexus_file_content,"	lset nst=6 rates=invgamma;");
nexus_file_content <- rbind(nexus_file_content,"	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); ");
nexus_file_content <- rbind(nexus_file_content,"	prset applyto=(all) ratepr=variable;");
nexus_file_content <- rbind(nexus_file_content,"	mcmcp ngen= 100000000 relburnin=yes burninfrac=0.25 printfreq=10000  samplefreq=10000 nchains=4 savebrlens=yes;");
nexus_file_content <- rbind(nexus_file_content,"	mcmc;");
nexus_file_content <- rbind(nexus_file_content,"	sumt;");
nexus_file_content <- rbind(nexus_file_content,"end;");
write(nexus_file_content,file=new_file_name);
}

ravel_polymorph_for_file <- function(polystates)	{
polystates <- sort(polystates,decreasing = FALSE);
return(paste("(",paste(polystates,collapse=""),")",sep=""));
}

convert_character_matrix_to_character <- function(chmatrix,UNKNOWN=-11,INAP=-22,missing="?",gap="-")	{
notu <- nrow(chmatrix);
nchars <- ncol(chmatrix);
chmatrix_char <- chmatrix;
allstates <- c(0:9,letter_states,more_letter_states);
for (ch in 1:nchars)	{
	these_states <- chmatrix[,ch];
	coded_notu <- (1:notu)[these_states>=0];
	gap_notu <- (1:notu)[these_states==INAP];
	miss_notu <- (1:notu)[these_states==UNKNOWN];
	poly_notu <- (1:notu)[!(1:notu) %in% c(coded_notu,gap_notu,miss_notu)];
	chmatrix_char[coded_notu,ch] <- as.character(allstates[1+chmatrix[coded_notu,ch]]);
	chmatrix_char[miss_notu,ch] <- missing;
	chmatrix_char[gap_notu,ch] <- gap;
	pn <- 0;
	while (pn < length(poly_notu))	{
		pn <- pn+1;
		chmatrix_char[poly_notu[pn],ch] <- ravel_polymorph_for_file(allstates[1+sort(unravel_polymorph_badass(chmatrix[poly_notu[pn],ch]))]);
		}
	}
return(chmatrix_char);
}

#### READ NEWICK FILES ####
#### convert (1,(2,3)) to vector_tree = 4 5 5 -1 4
read_newick_tree_from_chosen_file <- function() {
newicktree_file <- file.choose();
newick_tree <- scan(file=newicktree_file,what=character(),sep="\n");
nexus_string <- strsplit(newick_tree,split="",fixed=TRUE)[[1]]
nodes <- 0
for (i in 1:length(nexus_string))		if (nexus_string[i]=="(")	nodes <- nodes+1
# get clades
clades <- vector(length=nodes)
for (c in 1:nodes)	clades[c] <- c
# get taxa
notu <- p <- 0
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		otu <- as.numeric(nexus_string[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (nexus_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (nexus_string[i]==",")	{
		sp <- p <- 0
		} else if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		sp <- as.numeric(nexus_string[i])+(sp*10)
		p <- p+1
		if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

read_newick_tree_from_file <- function(newicktree_file) {
newick_tree <- scan(file=newicktree_file,what=character(),sep="\n")
nexus_string <- strsplit(newick_tree,split="",fixed=TRUE)[[1]]
nodes <- 0
for (i in 1:length(nexus_string))		if (nexus_string[i]=="(")	nodes <- nodes+1
# get clades
clades <- 1:nodes;
#clades <- vector(length=nodes)
#for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0;
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		otu <- as.numeric(nexus_string[i])+(otu * 10)
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (nexus_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (nexus_string[i]==",")	{
		sp <- p <- 0
		} else if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		sp <- as.numeric(nexus_string[i])+(sp*10);
		p <- p+1
		if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c];
		}
	}

return(vector_tree)
}

#(1,(2,4,(9,(7,(10,((18,(26,(33,36)),(34,40)),(19,27))),(13,25,(14,16)))),(3,15,(5,(11,17,24,(21,29),(22,31)),(28,(20,(30,32))))),(6,(8,39,(12,23,35,37,38)))))
#### convert vector_tree = 4 5 5 -1 4 to (1,(2,3))
#### 	where number is the htu number of the clade to which a species or htu belong
#### does not work yet
write_newick_string_with_taxon_names_from_vector_tree <- function(vector_tree,otu_names) {
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
notu <- match(-1,vector_tree)-1;
nNodes <- length(vector_tree)-notu;
all_names <- c(otu_names,(1:nNodes)+notu);
newick_string <- paste(match(-1,vector_tree),";",sep="");
for (nd in 1:nNodes)  {
  this_node <- notu+nd;
  all_names[mat_tree[nd,mat_tree[nd,]>0]]
  newick_string <- gsub(this_node,paste("(",paste(all_names[mat_tree[nd,mat_tree[nd,]>0]],collapse=","),")",sep=""),newick_string);
#  newick_string <- gsub(this_node,paste("(",paste(mat_tree[nd,mat_tree[nd,]>0],collapse=","),")",sep=""),newick_string)
  }
return(newick_string);
}

write_newick_string_from_vector_tree <- function(vector_tree) {
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
notu <- match(-1,vector_tree)-1;
nNodes <- length(vector_tree)-notu;
newick_string <- paste(match(-1,vector_tree),";",sep="");
for (nd in 1:nNodes)  {
  this_node <- notu+nd;
  newick_string <- gsub(this_node,paste("(",paste(mat_tree[nd,mat_tree[nd,]>0],collapse=","),")",sep=""),newick_string)
  }
return(newick_string);
}

# newick_string <- newick_string_taxa_only
read_newick_string <- function(newick_string) {
nodes <- 0;
if (length(newick_string)==1)	newick_string <- strsplit(newick_string,split="",fixed=TRUE)[[1]];
nodes <- sum(newick_string=="(");
#for (i in 1:length(newick_string))		if (newick_string[i]=="(")	nodes <- nodes+1;
# get clades
#clades <- vector(length=nodes);
#for (c in 1:nodes)	clades[c] <- c;
clades <- 1:nodes;
# get taxa
newick_string_one <- paste(newick_string,collapse="");
newick_string_one <- gsub("\\(","•",newick_string_one);
newick_string_one <- gsub("\\)","•",newick_string_one);
newick_string_one <- gsub(",","•",newick_string_one);
newick_string_one <- gsub(";","•",newick_string_one);
newick_string_one <- newick_string_one[newick_string_one!="•"];
newick_string_one <- strsplit(newick_string_one,"•")[[1]];
newick_string_one <- as.numeric(newick_string_one[newick_string_one!=""]);
#while (newick_string_two!=newick_string_one)	{
#	newick_string_two <- newick_string_one;
#	newick_string_one <- gsub("••","•",newick_string_one);
#	}
#options(warn=0);
#newick_string_one <- as.numeric(strsplit(newick_string_one,"•")[[1]]);
#options(warn=1);
#newick_string_one <- newick_string_one[!is.na(newick_string_one)];
notu <- max(newick_string_one);
#notu <- p <- 0
#for (i in 1:length(newick_string))	{
#	if (newick_string[i]>="0" && newick_string[i]<="9")	{
#		otu <- as.numeric(newick_string[i])+(otu * (10^p))
#		p <- p+1
#		if (otu>notu)	notu <- otu
#		} else {
#		p <- otu <- 0
#		}
#	}
vector_tree <- vector(length=notu+nodes);
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(newick_string))	{
	if (newick_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (newick_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (newick_string[i]==",")	{
		sp <- p <- 0
		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		sp <- as.numeric(newick_string[i])+(sp*10)
		p <- p+1
		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

# written 2022-08-19
count_taxa_in_newick_string <- function(newick_string)	{
# newick_string should look like: "(12,(7,(8,1,3,(2,(5,4,9,(11,(10,6)))))));"
newick_string <- gsub(";","",newick_string);
atomized_newick <- strsplit(newick_string,"")[[1]];
this_clade <- paste(atomized_newick,collapse="");
#	print(this_clade)
this_clade <- gsub("\\(","",this_clade);
this_clade <- gsub(")","",this_clade);
these_prog <- sort(as.numeric(strsplit(this_clade,",")[[1]]));
return(length(these_prog));
}

# written for Cinctan project
transform_newick_string_to_venn_tree <- function(newick_string)	{
# newick_string should look like: "(12,(7,(8,1,3,(2,(5,4,9,(11,(10,6)))))));"
atomized_newick <- strsplit(newick_string,"")[[1]];
l_m <- length(atomized_newick);
clade_bounds_l <- (1:l_m)[atomized_newick %in% "("];					# clade bounds left
clade_bounds_e <- clade_bounds_r <- (1:l_m)[atomized_newick %in% ")"];	# clade bounds right
nNodes <- length(clade_bounds_l);	# number of clades;
names(clade_bounds_l) <- names(clade_bounds_r) <- names(clade_bounds_e) <- 1:nNodes;
# get the first possible right paren ending this clade;
for (nn in 1:nNodes)	clade_bounds_e[nn] <- sum(clade_bounds_r>clade_bounds_l[nn]);
clade_bounds_e_unq <- unique(clade_bounds_e);
clade_boundaries <- clade_boundaries_orig <- cbind(clade_bounds_l,clade_bounds_r,clade_bounds_e);
#1 <- nn <- 1
while (length(clade_bounds_e_unq)>0)	{
	this_group <- sum(clade_boundaries_orig[,3]==clade_bounds_e_unq[1]);
	if (length(clade_bounds_e_unq)>1)	{
		this_group_starts <- this_group-sum(clade_boundaries_orig[,2]<clade_boundaries_orig[this_group+1,1]);
		} else	{
		this_group_starts <- 0;
		}
	clade_boundaries_orig[1:this_group,2] <- sort(clade_boundaries_orig[1:this_group,2],decreasing=TRUE);
	clade_boundaries[clade_boundaries[,1] %in% clade_boundaries_orig[,1],] <- clade_boundaries_orig;
	clade_boundaries_orig <- clade_boundaries_orig[!(1:nrow(clade_boundaries_orig)) %in% ((this_group_starts+1):this_group),]
	if (nrow(clade_boundaries_orig)>0)	clade_boundaries_orig[,2] <- sort(clade_boundaries_orig[,2]);
	n <- 0;
	while (n < nrow(clade_boundaries_orig))	{
		n <- n+1;
		clade_boundaries_orig[n,3] <- sum(clade_boundaries_orig[,2]>clade_boundaries_orig[n,1]);
		}
	clade_bounds_e_unq <- unique(clade_boundaries_orig[,3]);
	}

notu <- count_taxa_in_newick_string(newick_string);
venn_tree <- array(0,dim=c(nNodes,notu));
for (nn in 1:nNodes)	{
	lp <- clade_boundaries[nn,1];
	rp <- clade_boundaries[nn,2];
	this_clade <- paste(atomized_newick[lp:rp],collapse="");
#	print(this_clade)
	this_clade <- gsub("\\(","",this_clade);
	this_clade <- gsub(")","",this_clade);
	these_prog <- sort(as.numeric(strsplit(this_clade,",")[[1]]));
	venn_tree[nn,1:length(these_prog)] <- these_prog;
	}
return(venn_tree);
}

# written for Cinctan project
accersi_clade_boundaries_from_newick_string <- function(newick_string)	{
atomized_newick <- strsplit(newick_string,"")[[1]];
l_m <- length(atomized_newick);
clade_bounds_l <- (1:l_m)[atomized_newick %in% "("];
clade_bounds_e <- clade_bounds_r <- (1:l_m)[atomized_newick %in% ")"];
nNodes <- length(clade_bounds_l);	# number of clades;
names(clade_bounds_l) <- names(clade_bounds_r) <- names(clade_bounds_e) <- 1:nNodes;
# get the first possible right paren ending this clade;
for (nn in 1:nNodes)	clade_bounds_e[nn] <- sum(clade_bounds_r>clade_bounds_l[nn])
clade_bounds_e_unq <- unique(clade_bounds_e);
clade_boundaries <- clade_boundaries_orig <- cbind(clade_bounds_l,clade_bounds_r,clade_bounds_e);
#1 <- nn <- 1
while (length(clade_bounds_e_unq)>0)	{
	this_group <- sum(clade_boundaries_orig[,3]==clade_bounds_e_unq[1]);
	if (length(clade_bounds_e_unq)>1)	{
		this_group_starts <- this_group-sum(clade_boundaries_orig[,2]<clade_boundaries_orig[this_group+1,1]);
		} else	{
		this_group_starts <- 0;
		}
	clade_boundaries_orig[1:this_group,2] <- sort(clade_boundaries_orig[1:this_group,2],decreasing=TRUE);
	clade_boundaries[clade_boundaries[,1] %in% clade_boundaries_orig[,1],] <- clade_boundaries_orig;
	clade_boundaries_orig <- clade_boundaries_orig[!(1:nrow(clade_boundaries_orig)) %in% ((this_group_starts+1):this_group),]
	if (nrow(clade_boundaries_orig)>0)	clade_boundaries_orig[,2] <- sort(clade_boundaries_orig[,2]);
	n <- 0;
	while (n < nrow(clade_boundaries_orig))	{
		n <- n+1;
		clade_boundaries_orig[n,3] <- sum(clade_boundaries_orig[,2]>clade_boundaries_orig[n,1]);
		}
	clade_bounds_e_unq <- unique(clade_boundaries_orig[,3]);
	}
clade_boundaries <- data.frame(lp=as.numeric(clade_boundaries[,1]),
							   rp=as.numeric(clade_boundaries[,2]))
return(clade_boundaries);
}

#### 	where number is the htu number of the clade to which a species or htu belong
# newick_string_ancestored <- newick_string_taxa_only_raw
# written for Cinctan project
# updated 2020-12-30 to allow outgroup to be ancestral
#newick_string_ancestored <- newick_string_full
#i <- match("Rozanovicystis_triangularis",otu_names)
identify_newick_ancestors <- function(newick_string_ancestored,otu_names)	{
notu <- length(otu_names);
ancestors <- c();
for (i in 1:notu)	{
	flibberty <- strsplit(newick_string_ancestored,otu_names[i])[[1]][1];
	jibberty <- strsplit(flibberty,"")[[1]];
	if (jibberty[length(jibberty)]==")")	ancestors <- c(ancestors,otu_names[i]);
	}
return(ancestors);
}

#newick_string_ancestored <- newick_string_full
find_newick_ancestors <- function(newick_string_ancestored)	{
atomic_ancestral <- strsplit(newick_string_ancestored,"")[[1]];
a_a <- length(atomic_ancestral);
l_paren <- (1:a_a)[atomic_ancestral=="("];
r_paren <- (1:a_a)[atomic_ancestral==")"];
sisters <- (1:a_a)[atomic_ancestral==","];
otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
notu <- 1;
breaks <- c();
for (i in 2:length(otu_nos))	{
	if ((otu_nos[i]-1)>otu_nos[i-1])	{
		notu <- notu+1;
		breaks <- c(breaks,otu_nos[i]-1);
		}
	}
sampled_ancestors <- array(0,dim=notu);
ancestral_starts <- 1+breaks[atomic_ancestral[breaks]==")"];
ab <- 0;
while (ab < length(ancestral_starts))	{
	ab <- ab+1;
	dd <- ancestral_starts[ab];
	st_hr <- match(dd,otu_nos);
	this_anc <- as.numeric(atomic_ancestral[otu_nos[st_hr]]);
	while(st_hr < length(otu_nos) && otu_nos[st_hr+1]==(otu_nos[st_hr]+1))	{
		st_hr <- st_hr+1;
		this_anc <- (10*this_anc)+as.numeric(atomic_ancestral[otu_nos[st_hr]]);
		}
	sampled_ancestors[this_anc] <- 1;
	}
return(sampled_ancestors);
}

# written for Cinctan project
#newick_string_ancestored <- newick_string_full
fix_newick_ancestors <- function(newick_string_ancestored)	{
atomic_ancestral <- strsplit(newick_string_ancestored,"")[[1]];
a_a <- length(atomic_ancestral);
l_paren <- (1:a_a)[atomic_ancestral=="("];
r_paren <- (1:a_a)[atomic_ancestral==")"];
sisters <- (1:a_a)[atomic_ancestral==","];
otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
for (rp in length(r_paren):1)	{
	if (!is.na(match(1,otu_nos-r_paren[rp])))	{
		an <- r_paren[rp];
		an_no <- c();
		while (atomic_ancestral[an+1] %in% as.character(0:9))	{
			an_no  <- c(an_no,an+1);
			an <- an+1;
			}
		atomic_ancestral <- c(atomic_ancestral[1:(r_paren[rp]-1)],
							  ",",
							  atomic_ancestral[an_no],
							  ")",
							  atomic_ancestral[(an+1):a_a]);
		a_a <- length(atomic_ancestral);
		l_paren <- (1:a_a)[atomic_ancestral=="("];
		r_paren <- (1:a_a)[atomic_ancestral==")"];
		sisters <- (1:a_a)[atomic_ancestral==","];
		otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
		otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
		}
	}
revised_newick_string <- paste(atomic_ancestral,collapse="");
return(revised_newick_string);
}

# written to deal with BEAST
#newick_string <- mctree2;
extract_otu_names_from_newick_string <- function(newick_string)	{
newick_string <- gsub("\\(","•",newick_string);
newick_string <- gsub("\\)","•",newick_string);
newick_string <- gsub(",","•",newick_string);
newick_string <- gsub(":","•",newick_string);
newick_string <- gsub("\\{","•",newick_string);
newick_string <- gsub("\\}","•",newick_string);
newick_string <- gsub("\\[","•",newick_string);
newick_string <- gsub("\\[","•",newick_string);
newick_string <- gsub("\\;","•",newick_string);
for (i in 0:9)	newick_string <- gsub(i,"",newick_string);
newick_string <- gsub("\\.","•",newick_string);
taxon_names <- strsplit(newick_string,"•")[[1]];
taxon_names <- taxon_names[taxon_names!=""];
taxon_names <- sort(taxon_names);
return(taxon_names);
}

#newick_string_bl <- mctree2;
read_newick_string_for_branch_lengths <- function(newick_string_bl,taxon_names=NULL)	{
if (is.null(taxon_names))
	taxon_names <- extract_otu_names_from_newick_string(newick_string_bl);
notu <- length(taxon_names);
nclades <- sum(strsplit(newick_string_bl,"")[[1]]=="(");
name_length <- sapply(taxon_names,count_characters_in_string);
name_length <- sort(name_length,decreasing = TRUE);
branch_length <- vector(length=(notu+nclades));

clade_names <- rep(paste("node_",paste(rep(0,ceiling(log10(nclades))),collapse=""),sep=""),nclades);
for (cl in 1:nclades)	{
	atomize <- strsplit(clade_names[cl],"")[[1]];
	latom <- length(atomize);
	dgts <- floor(log10(cl));
	atomize[(latom-dgts):latom] <- strsplit(as.character(cl),"")[[1]];
	clade_names[cl] <- paste(atomize,collapse="");
	}
names(branch_length) <- c(taxon_names,clade_names);

newick_string_bl_taxa_only <- newick_string_bl;
for (i in 0:9)	newick_string_bl_taxa_only <- gsub(i,"",newick_string_bl_taxa_only);
newick_string_bl_taxa_only <- gsub(":","",newick_string_bl_taxa_only);
newick_string_bl_taxa_only <- gsub("\\.","",newick_string_bl_taxa_only);
for (sp in 1:notu)	{
	spn <- match(names(name_length)[sp],taxon_names);
	newick_string_bl_taxa_only <- gsub(names(name_length)[sp],spn,newick_string_bl_taxa_only);
	newick_string_bl <- gsub(names(name_length)[sp],"•",newick_string_bl);
	newick_string_bl_split <- strsplit(newick_string_bl,"•")[[1]][2];
	atomized <- strsplit(newick_string_bl_split,"")[[1]];
	next_paren <- match(")",atomized);
	next_comma <- match(",",atomized);
	if (is.na(next_comma))	next_comma <- length(atomized);
	branch_length[spn] <- as.numeric(paste(atomized[2:min(next_comma-1,next_paren-1)],collapse=""));
	newick_string_bl <- gsub("•",spn,newick_string_bl);
#	newick_string_bl <- gsub("•",spn,newick_string_bl);
	}
vector_tree <- read_newick_string(newick_string_bl_taxa_only);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
atomized <- strsplit(newick_string_bl,"")[[1]];
clade_strs <- (1:length(atomized))[atomized %in% "("];
names(clade_strs) <- 1:nclades;
clade_ends <- (1:length(atomized))[atomized %in% ")"];

for (i in 1:length(clade_ends)){
	rp <- max(clade_strs[clade_strs<clade_ends[i]]);
	names(clade_ends)[i] <- names(clade_strs)[match(rp,clade_strs)];
	clade_strs <- clade_strs[clade_strs!=rp];
	}
colons <- (1:length(atomized))[atomized %in% ":"];
names(clade_ends) <- length(clade_ends):1;
clade_colons <- 1+clade_ends[(clade_ends+1) %in% colons];
for (i in 1:length(clade_ends))	{
	htu <- as.numeric(names(clade_colons)[i])+notu;
	brl <- "";
	j <- 1+clade_colons[i];
	while (atomized[j] %in% c(".",0:9))	{
		brl <- paste(brl,atomized[j],sep="");
		j <- j+1;
		}
	branch_length[htu] <- as.numeric(brl);
	}
return(branch_length);
}

#named_newick_string <- newick_string <- mctree2;
#newick_string <- newick_string_bl;
transmogrify_named_newick_to_numbered_newick <- function(named_newick_string,taxon_names=NULL)	{
if (is.null(taxon_names))
		taxon_names <- extract_otu_names_from_newick_string(named_newick_string);
notu <- length(taxon_names)
name_length <- sapply(taxon_names,count_characters_in_string);
name_length <- sort(name_length,decreasing = TRUE);
numbered_newick_string <- named_newick_string;
for (sp in 1:notu)	{
	txno <- match(names(name_length)[sp],taxon_names);
	numbered_newick_string <- gsub(taxon_names[txno],txno,numbered_newick_string);
	}
return(numbered_newick_string);
}

#newick_string <- mctree2
read_newick_string_with_branch_lengths <- function(newick_string,taxon_names=NULL) {
if (is.null(taxon_names) && gsub("u","",newick_string)!=newick_string)	{
	taxon_names <- extract_otu_names_from_newick_string(newick_string);
	notu <- length(taxon_names);
	numbered_newick_string <- transmogrify_named_newick_to_numbered_newick(newick_string)
	} else if (is.null(taxon_names))	{
	taxon_names <- 1:notu;
	numbered_newick_string <- newick_string;
	} else	{
	numbered_newick_string <- transmogrify_named_newick_to_numbered_newick(newick_string,taxon_names=taxon_names);
	}
nodes <- 0;
if (length(numbered_newick_string)==1)	numbered_newick_string <- strsplit(numbered_newick_string,split="",fixed=TRUE)[[1]];
nodes <- sum(numbered_newick_string=="(");
node_names <- rep(paste("node_",paste(rep(0,ceiling(log10(nodes))),collapse=""),sep=""),nodes);
for (cl in 1:nodes)	{
	atomize <- strsplit(node_names[cl],"")[[1]];
	latom <- length(atomize);
	dgts <- floor(log10(cl));
	atomize[(latom-dgts):latom] <- strsplit(as.character(cl),"")[[1]];
	node_names[cl] <- paste(atomize,collapse="");
	}
branch_length <- rep(0,length(c(taxon_names,node_names)));
names(branch_length) <- c(taxon_names,node_names);

clades <- rep(-1,nodes);
# get taxa
newick_string_one <- paste(numbered_newick_string,collapse="");
newick_string_one <- gsub("\\(","•",newick_string_one);
newick_string_one <- gsub("\\)","•",newick_string_one);
newick_string_one <- gsub(",","•",newick_string_one);
newick_string_one <- gsub(";","•",newick_string_one);
newick_string_one <- gsub(":","•",newick_string_one);
newick_string_one <- newick_string_one[newick_string_one!="•"];
newick_string_one <- strsplit(newick_string_one,"•")[[1]];
newick_string_one <- newick_string_one[newick_string_one!=""];
newick_string_one <- as.numeric(newick_string_one[newick_string_one!=""]);
notu <- max(newick_string_one[ceiling(newick_string_one)==newick_string_one])

htu_height <- branch_length <- vector_tree <- vector(length=notu+nodes);
names(htu_height) <- names(branch_length) <- names(vector_tree) <- c(taxon_names,node_names);
#for (c in 1:nodes)	clades[c] <- -1;
cl <- c <- 0;
i <- 1;
#for (i in 1:length(newick_string))	{
while (i <= length(numbered_newick_string)) {
	if (numbered_newick_string[i]=="(")	{
		sp <- p <- 0;
		cl <- cl+1;
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu;
			htu_height[notu+cl] <- c;
			} else vector_tree[notu+1] <- -1
		c <- c+1;
		clades[c] <- cl;
		} else if (numbered_newick_string[i]==")")	{
		last_taxon <- notu+clades[c];
		c <- c-1;
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==",")	{
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==":")	{
		brl <- "";
		while (numbered_newick_string[i+1] %in% c(".",0:9))	{
			brl <- paste(brl,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}
		branch_length[last_taxon] <- as.numeric(paste(brl,collapse=""));
#		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		} else if (numbered_newick_string[i] %in% 0:9)	{
#		sp <- as.numeric(newick_string[i])+(sp*10)
		sp <- numbered_newick_string[i];
		while (numbered_newick_string[i+1] %in% 0:9)	{
			sp <- paste(sp,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}		#p <- p+1;
		sp <- last_taxon <- as.numeric(sp);
		vector_tree[sp] <- notu+clades[c];
		htu_height[sp] <- c;
#		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c];
		}
#	print(paste("clade",max(clades[clades!=-1])))
#	print(paste(numbered_newick_string[1:i],collapse=""))
#	print(branch_length)
	i <- i+1;
	}
output <- list(vector_tree,branch_length,htu_height);
names(output) <- c("vector_tree","branch_lengths","htu_heights");
return(output)
}

#numbered_newick_string <- newick_string_taxa_only_atomized
accersi_ttu_heights <- function(numbered_newick_string)	{
if (length(numbered_newick_string)==1)	numbered_newick_string <- strsplit(numbered_newick_string,"")[[1]];
vector_tree <- read_newick_string(numbered_newick_string);
cl <- c <- 0;
i <- 1;
notu <- match(-1,vector_tree)-1;
clades <- length(unique(vector_tree[vector_tree>0]));
ttu_height <- vector(length=notu+clades);
#for (i in 1:length(newick_string))	{
while (i <= length(numbered_newick_string)) {
	if (numbered_newick_string[i]=="(")	{
		sp <- p <- 0;
		cl <- cl+1;
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu;
			ttu_height[notu+cl] <- c;
			} else vector_tree[notu+1] <- -1
		c <- c+1;
		clades[c] <- cl;
		} else if (numbered_newick_string[i]==")")	{
		last_taxon <- notu+clades[c];
		c <- c-1;
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==",")	{
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==":")	{
		brl <- "";
		while (numbered_newick_string[i+1] %in% c(".",0:9))	{
			brl <- paste(brl,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}
#		branch_length[last_taxon] <- as.numeric(paste(brl,collapse=""));
#		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		} else if (numbered_newick_string[i] %in% 0:9)	{
#		sp <- as.numeric(newick_string[i])+(sp*10)
		sp <- numbered_newick_string[i];
		while (numbered_newick_string[i+1] %in% 0:9)	{
			sp <- paste(sp,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}		#p <- p+1;
		sp <- last_taxon <- as.numeric(sp);
		vector_tree[sp] <- notu+clades[c];
		ttu_height[sp] <- c;
#		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c];
		}
#	print(paste("clade",max(clades[clades!=-1])))
#	print(paste(numbered_newick_string[1:i],collapse=""))
#	print(branch_length)
	i <- i+1;
	}
return(ttu_height);
}

# write this out! 2024-06-14
accersi_daughter_clades <- function(numbered_newick_string)	{
if (length(numbered_newick_string)==1)	numbered_newick_string <- strsplit(numbered_newick_string,"")[[1]];
vector_tree <- read_newick_string(numbered_newick_string);
cl <- c <- 0;
i <- 1;
notu <- match(-1,vector_tree)-1;
clades <- length(unique(vector_tree[vector_tree>0]));
ttu_height <- vector(length=notu+clades);
#for (i in 1:length(newick_string))	{
while (i <= length(numbered_newick_string)) {
	if (numbered_newick_string[i]=="(")	{
		sp <- p <- 0;
		cl <- cl+1;
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu;
			ttu_height[notu+cl] <- c;
			} else vector_tree[notu+1] <- -1
		c <- c+1;
		clades[c] <- cl;
		} else if (numbered_newick_string[i]==")")	{
		last_taxon <- notu+clades[c];
		c <- c-1;
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==",")	{
		sp <- p <- 0;
		} else if (numbered_newick_string[i]==":")	{
		brl <- "";
		while (numbered_newick_string[i+1] %in% c(".",0:9))	{
			brl <- paste(brl,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}
#		branch_length[last_taxon] <- as.numeric(paste(brl,collapse=""));
#		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		} else if (numbered_newick_string[i] %in% 0:9)	{
#		sp <- as.numeric(newick_string[i])+(sp*10)
		sp <- numbered_newick_string[i];
		while (numbered_newick_string[i+1] %in% 0:9)	{
			sp <- paste(sp,numbered_newick_string[i+1],sep="");
			i <- i+1;
			}		#p <- p+1;
		sp <- last_taxon <- as.numeric(sp);
		vector_tree[sp] <- notu+clades[c];
		ttu_height[sp] <- c;
#		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c];
		}
#	print(paste("clade",max(clades[clades!=-1])))
#	print(paste(numbered_newick_string[1:i],collapse=""))
#	print(branch_length)
	i <- i+1;
	}
return(ttu_height);
}
# written for Cinctan project
# heavily modified 2020-12
# modified 2024-06-13 to take into account having both "asper" and "asperiformis";
#read_newick_string_mcmc <- function(newick_string_full,otu_names) {
read_newick_string_mcmc <- function(newick_string_full,otu_names) {
otu_names_nex <- gsub(" ","_",otu_names);
otu_names_nex[gsub("\\(","",otu_names_nex)!=otu_names_nex] <- paste("\"",otu_names[gsub("\\(","",otu_names_nex)!=otu_names_nex],"\"",sep="");
otu_name_lengths <- sapply(otu_names_nex,name_length);
otu_name_lengths <- sort(otu_name_lengths,decreasing = TRUE);

simple_newick_string <- molecularize <- strsplit(newick_string_full,split="")[[1]];
left_brackets <- (1:length(molecularize))[molecularize %in% "["];
right_brackets <- (1:length(molecularize))[molecularize %in% "]"];
for (br in 1:length(left_brackets))	molecularize[left_brackets[br]:right_brackets[br]] <- "";
newick_string_taxa_only_rawwest <- newick_string_taxa_raw <- newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
notu <- length(otu_names);
branch_durations <- array(0,dim=(2*notu)-1);
for (nn in 1:notu)	{
	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	dd <- 1+match("•",dummy_newick);
	b_d <- dummy_newick[dd];
	dd <- dd+1;
	while (dummy_newick[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,dummy_newick[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[nn] <- as.numeric(b_d);
	}
for (i in 0:9)	newick_string_taxa_only <- gsub(i,"",newick_string_taxa_only);
newick_string_taxa_only <- gsub(":","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("-","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("\\.","",newick_string_taxa_only);
for (mm in 1:notu)	{
#mm <- 0;
#while (mm < notu)	{
#	mm <- mm+1;
	nn <- match(names(otu_name_lengths)[mm],otu_names_nex);
	newick_string_taxa_only_old <- newick_string_taxa_only;
	newick_string_taxa_only <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_only);
	newick_string_taxa_raw <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_raw);
#	print(newick_string_taxa_only);
	if (newick_string_taxa_only_old==newick_string_taxa_only)	{
		print(paste("problem taxon",nn));
		}
	}
#xx <- read_newick_string_with_branch_lengths(newick_string_taxa_raw);
ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
names(ancestral) <- otu_names_nex;
ancestorless_tree <- fix_newick_ancestors(newick_string_taxa_only);
newick_string_taxa_only_atomized <- strsplit(ancestorless_tree,"")[[1]];
ttu_heights <- accersi_ttu_heights(newick_string_taxa_only_atomized);
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa]!=";")	{
	newick_string_taxa_only_atomized <- c(newick_string_taxa_only_atomized,";");
	nstoa <- length(newick_string_taxa_only_atomized);
	}
newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
#ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
#names(ancestral) <- otu_names_nex;
newick_string_taxa_only_raw <- newick_string_taxa_only;
newick_string_taxa_only <- fix_newick_ancestors(newick_string_taxa_only);
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa-1]!=")")	{
	newick_string_taxa_only_atomized[nstoa] <- ")";
	newick_string_taxa_only_atomized <- c("(",newick_string_taxa_only_atomized,";");
	newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
	}

vector_tree <- read_newick_string(newick_string_taxa_only);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
newick_string_taxa_raw <- strsplit(newick_string_taxa_raw,split="")[[1]];
clade_ends <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ")"];
colons <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ":"];
names(clade_ends) <- length(clade_ends):1;
clade_colons <- 1+clade_ends[(clade_ends+1) %in% colons];
for (cc in 1:length(clade_colons))	{
	n_node <- notu+as.numeric(names(clade_colons)[cc]);
	dd <- clade_colons[cc]+1;
	b_d <- newick_string_taxa_raw[dd];
	dd <- dd+1;
	while (newick_string_taxa_raw[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,newick_string_taxa_raw[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[n_node] <- as.numeric(b_d);
	}
nNodes <- length(clade_ends);
newick_string_taxa_raw <- paste(newick_string_taxa_raw,collapse="");
if (nNodes<10)	{
	node_names <- paste("node_",1:nNodes,sep="")
	} else if (nNodes<100)	{
	node_names <- c(paste("node_0",1:9,sep=""),paste("node_",10:nNodes,sep=""));
	} else	{
	node_names <- c(paste("node_00",1:9,sep=""),paste("node_0",10:99,sep=""),paste("node_",100:nNodes,sep=""));
	}
names(branch_durations) <- c(otu_names_nex,node_names);
#rbind(xx$branch_lengths,branch_durations)

#vector_tree_raw <- read_newick_string(newick_string_taxa_only_raw);
#vector_tree <- read_newick_string(newick_string_taxa_only);
molecularize <- strsplit(newick_string_taxa_only,split="")[[1]];
if (match("(",molecularize)>1)	{
	mm <- match("(",molecularize);
	molecularize <- molecularize[mm:length(molecularize)];
	newick_string_taxa_only <- paste(molecularize,collapse="");
	}

venn_tree_newick <- transform_newick_string_to_venn_tree(newick_string = newick_string_taxa_only);

# have: newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,ancestral,branch_durations);
# need: clade_posteriors, prob_ancestor,hpd_age;
# uset mat_tree and newick_string_taxa_only to figure out which node is what number
newick_rem_info <- gsub("sampled_ancestor=","•",newick_string_full);
newick_rem_info <- gsub("height_95%_HPD","age_95%_HPD",newick_rem_info);
if (gsub("age_95%_HPD=","§",newick_rem_info)!=newick_rem_info)	{
	newick_rem_info <- gsub("age_95%_HPD=","§",newick_rem_info);
	} else	{
	newick_rem_info <- gsub("brlen_95%_HPD=","§",newick_rem_info);
	}
newick_rem_info <- gsub("posterior=","¶",newick_rem_info);
hpd_age <- data.frame(lb=as.numeric(rep(0,notu+nNodes)),ub=as.numeric(rep(0,notu+nNodes)));
molecularized <- strsplit(newick_rem_info,"")[[1]];# print("line 2275");
molecules <- length(molecularized);
panc_boundaries <- (1:molecules)[molecularized=="•"];
post_boundaries <- (1:molecules)[molecularized=="¶"];
hpd_age_boundaries <- (1:molecules)[molecularized=="§"];
taxon_boundaries <- array(0,dim=c(notu,2)); #print("line 2280");
for (mm in 1:notu)	{
	nn <- match(names(otu_name_lengths)[mm],otu_names_nex);
	taxon_dummy <- gsub(paste(otu_names_nex[nn],"\\[",sep=""),"£[",newick_rem_info);
#	if (length(strsplit(taxon_dummy,"£")[[1]])>2)	{
#		
#		}
	taxon_boundaries[nn,1] <- (1:length(strsplit(taxon_dummy,"")[[1]]))[strsplit(taxon_dummy,"")[[1]]=="£"];
	taxon_boundaries[nn,2] <- taxon_boundaries[nn,1]+length(strsplit(otu_names_nex[nn],"")[[1]])-1;
	}
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_rem_info);
colnames(taxon_boundaries) <- colnames(clade_boundaries);
tu_boundaries <- rbind(taxon_boundaries,clade_boundaries);
rownames(tu_boundaries) <- all_names <- c(otu_names_nex,node_names);
tu_boundaries <- tu_boundaries[order(tu_boundaries$rp),];# print("line 2290");
brackets_l <- (1:molecules)[molecularized=="{"];
brackets_r <- (1:molecules)[molecularized=="}"];
hp <- 0;
while (hp < length(hpd_age_boundaries))	{
	hp <- hp+1;
	tx <- sum(tu_boundaries$rp<hpd_age_boundaries[hp]);
	txn <- match(rownames(tu_boundaries)[tx],all_names);
	i <- brackets_l[1+sum(hpd_age_boundaries[hp]>brackets_l)]+1;
	j <- brackets_r[1+sum(hpd_age_boundaries[hp]>brackets_r)]-1;
	# put lower bound 
	hpd_age$lb[txn] <- max(as.numeric(strsplit(paste(molecularized[i:j],collapse=""),",")[[1]]));
	hpd_age$ub[txn] <- min(as.numeric(strsplit(paste(molecularized[i:j],collapse=""),",")[[1]]));
	}
rownames(hpd_age) <- all_names;
hpd_age[vector_tree[(1:notu)[ancestral==1]],] <- hpd_age[(1:notu)[ancestral==1],];

rownames(taxon_boundaries) <- otu_names_nex;#  print("line 2305");
taxon_boundaries <- taxon_boundaries[order(taxon_boundaries[,1]),];
brackets_l <- (1:molecules)[molecularized=="["];
brackets_r <- (1:molecules)[molecularized=="]"];
prob_ancestor <- array(0,dim=notu);
pb <- 0;  #print("line 2310");
while (pb < length(panc_boundaries))	{
#for (pb in 1:length(panc_boundaries))	{
	pb <- pb+1;
	tx <- sum(taxon_boundaries[,1]<panc_boundaries[pb]);
	txn <- match(rownames(taxon_boundaries)[tx],otu_names_nex);
	i <- panc_boundaries[pb]+1;
	pranc <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		pranc <- paste(pranc,molecularized[i],sep="");
		i <- i+1;
		}
	prob_ancestor[txn] <- as.numeric(pranc);
	}
names(prob_ancestor) <- otu_names_nex;
ancestral[prob_ancestor>=0.5] <- 1; # added 2023-09-25: we seem to miss some ancestors!

clade_posteriors <- array(0,dim=nNodes);#  print("line 2325");
rownames(clade_boundaries) <- names(clade_posteriors) <- node_names;
clade_boundaries <- clade_boundaries[order(clade_boundaries$rp),];
pp <- 0;
while (pp < length(post_boundaries))	{
	pp <- pp+1;
	cl <- sum(clade_boundaries$rp<post_boundaries[pp]);
	cld <- match(rownames(clade_boundaries)[cl],node_names);
	i <- post_boundaries[pp]+1;
	postp <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		postp <- paste(postp,molecularized[i],sep="");
		i <- i+1;
		}
	clade_posteriors[cld] <- as.numeric(postp);
#	molecularized[i:i+-1:10]
	}

hpd_age$md <- (hpd_age$lb+hpd_age$ub)/2;
for (mt in nrow(mat_tree):1)	{
	htu <- mt+notu;
	#branch_durations[mat_tree[mt,mat_tree[mt,]>0]]
	if (sum(hpd_age$md[mat_tree[mt,mat_tree[mt,]>0]]>hpd_age$md[htu])>0)	{
#		hpd_age$md[htu] <- min(hpd_age$md[mat_tree[mt,mat_tree[mt,]>0]])-branch_durations[htu]
		offset <- max(hpd_age$md[mat_tree[mt,mat_tree[mt,]>0]]-hpd_age$md[htu]);
		min(hpd_age$md[mat_tree[mt,mat_tree[mt,]>0]])-branch_durations[htu]
		hpd_age[htu,] <- hpd_age[htu,]+offset;
#		print(mt);
		}
	}
#sum(strsplit(newick_rem_info,"")[[1]]=="•");
#sum(strsplit(newick_rem_info,"")[[1]]=="§");
#sum(strsplit(newick_rem_info,"")[[1]]=="¶");
output <- list(newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,clade_posteriors,ancestral,prob_ancestor,hpd_age,branch_durations);
names(output) <- c("newick","newick_raw","vector_tree","clade_posteriors","ancestral","ancestored_nodes","hpd_age","branch_durations");
return(output);
}

read_newick_string_mcmc_numbered <- function(newick_string_full,otu_names) {
atomized <- strsplit(newick_string_full,"")[[1]];
atomized <- atomized[match("(",atomized):length(atomized)]
newick_string_full <- paste(atomized,collapse = "");

otu_names_nex <- gsub(" ","_",otu_names);
otu_names_nex[gsub("\\(","",otu_names_nex)!=otu_names_nex] <- paste("\"",otu_names[gsub("\\(","",otu_names_nex)!=otu_names_nex],"\"",sep="");

simple_newick_string <- molecularize <- strsplit(newick_string_full,split="")[[1]];
left_brackets <- (1:length(molecularize))[molecularize %in% "["];
right_brackets <- (1:length(molecularize))[molecularize %in% "]"];
for (br in 1:length(left_brackets))	molecularize[left_brackets[br]:right_brackets[br]] <- "";
newick_string_taxa_raw <- newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
notu <- length(otu_names);
branch_durations <- array(0,dim=(2*notu)-1);
node_names <- c();
for (i in 1:(notu-1))	node_names <- c(node_names,paste("node_",paste(rep(0,floor(log10(notu-1))-floor(log10(i+0.5))),collapse=""),i,sep=""));
names(branch_durations) <- c(otu_names_nex,node_names);
for (nn in 1:notu)	{
	try1 <- paste("\\(",nn,":",sep="");
	rep1 <- "(•";
	try2 <- paste(",",nn,":",sep="");
	rep2 <- ",•";
	dummy_newick <- gsub(try1,rep1,newick_string_taxa_only);
	dummy_newick <- gsub(try2,rep2,dummy_newick);
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	dd <- 1+match("•",dummy_newick);
	b_d <- dummy_newick[dd];
	dd <- dd+1;
	while (dummy_newick[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,dummy_newick[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[nn] <- as.numeric(b_d);
#	names(branch_durations)[nn] <- otu_names_nex[nn];
	}
#for (i in 0:9)	newick_string_taxa_only <- gsub(i,"",newick_string_taxa_only);
simple_newick_string <- paste(molecularize,collapse=""); 

key <- data.frame(variable=as.character(c("age","duration","posterior")),symbol=as.character(c("§","¢","¶")));
newick_string_denoted <- gsub("height_95%_HPD","age_95%_HPD",newick_string_full);
newick_string_denoted <- gsub("age_95%_HPD=",key$symbol[key$variable %in% c("age","height")],newick_string_denoted);
newick_string_denoted <- gsub("length_95%_HPD=",key$symbol[key$variable %in% c("duration","length")],newick_string_denoted);
newick_string_denoted <- gsub("posterior=",key$symbol[key$variable %in% "posterior"],newick_string_denoted);
# find where clade left & right parens are in character string.
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_string_denoted);
nNodes <- nrow(clade_boundaries);
otu_ends <- vector(length=notu);
options(warn=2);
for (nn in 1:notu)	{
	taxon_label <- paste(rep("†",1+floor(log10(nn))),collapse="")
	try1 <- paste("\\(",nn,"\\[",sep="");
	rep1 <- paste("(",taxon_label,"[",sep="");
	try2 <- paste(",",nn,"\\[",sep="");
	rep2 <- paste(",",taxon_label,"[",sep="");
	dummy_newick <- gsub(try1,rep1,newick_string_denoted);
	dummy_newick <- gsub(try2,rep2,dummy_newick);
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	otu_ends[nn] <- max((1:length(dummy_newick))[dummy_newick %in% "†"]);
	}
names(otu_ends) <- otu_names_nex;
#we now have boundaries for taxa & clades in simple_newick_string;
molecularize2 <- strsplit(newick_string_denoted,"")[[1]];
right_brackets <- (1:length(molecularize2))[molecularize2 %in% "]"];
left_brackets <- (1:length(molecularize2))[molecularize2 %in% "["];
colons <- (1:length(molecularize2))[molecularize2 %in% ":"];
commas <- (1:length(molecularize2))[molecularize2 %in% ","];
taxon_breaks <- sort(c(commas,left_brackets,clade_boundaries$rp));
for (cc in 2:nNodes)	{
#	htu <- notu+length(clade_starts[clade_starts<clade_ends[cc]]);
	htu <- cc+notu;
	this_colon <- colons[colons>clade_boundaries$rp[cc]][1]+1;
	this_taxon_break <- taxon_breaks[taxon_breaks>this_colon][1]-1;
	branch_durations[htu] <- as.numeric(paste(molecularize2[this_colon:this_taxon_break],collapse=""));
	}

# delete branch durations from molecularize;
molecularize <- strsplit(paste(molecularize,collapse=""),"")[[1]];
colons1 <- (1:length(molecularize))[molecularize %in% ":"];
commas1 <- (1:length(molecularize))[molecularize %in% ","];
clade_boundaries1 <- accersi_clade_boundaries_from_newick_string(paste(molecularize,collapse=""));
taxon_ends <- sort(c(commas1,clade_boundaries1$lp,clade_boundaries1$rp));
for (cc in 1:length(colons1))	{
	str <- colons1[cc];
	if (length(taxon_ends[taxon_ends>colons1[cc]])>0)	{
		end <- min(taxon_ends[taxon_ends>colons1[cc]])-1;
		} else	{
		end <- length(molecularize)-1;
		}
	molecularize[str:end] <- "";
	}

# get cladogram;
newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
vector_tree <- read_newick_string(newick_string_taxa_only);			# we have cladogram;
names(vector_tree) <- c(otu_names_nex,node_names);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);	# we have cladogram 2;

# get ancestors;
ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
print("naming ancestors");
names(ancestral) <- otu_names_nex;
if (sum(ancestral)==0)
	ancestral[branch_durations[1:notu]==0] <- 1;
ancestored_nodes <- unique(vector_tree[(1:notu)[ancestral==1]]);

# get hpd information
hpd_dur <- hpd_age <- data.frame(lb=as.numeric(rep(0,notu+nNodes)),ub=as.numeric(rep(0,notu+nNodes)));
rownames(hpd_age) <- rownames(hpd_dur) <- names(branch_durations);
#molecularize2 <- strsplit(newick_string_denoted,"")[[1]];
for (nn in 1:notu)	{
#	ancestral[nn]
	lb <- min(left_brackets[left_brackets>otu_ends[nn]]);
	rb <- min(right_brackets[right_brackets>lb]);
	taxon_info <- molecularize2[(lb-2):rb];
	age_break <- match(key$symbol[key$variable=="age" ],taxon_info);
	if (!is.na(age_break))	{
		relv_taxon_info <- taxon_info[age_break:length(taxon_info)];
		curlies <- (1:length(relv_taxon_info))[relv_taxon_info %in% c("{","}")][1:2];
		hpd_age[nn,] <- as.numeric(strsplit(paste(relv_taxon_info[(1+curlies[1]):(curlies[2]-1)],collapse=""),",")[[1]]);
		} else	{
		taxon_info <- paste(taxon_info,collapse="");
		taxon_info <- gsub("height=",key$symbol[key$variable=="age"],taxon_info);
		taxon_info <- gsub("age=",key$symbol[key$variable=="age"],taxon_info);
		taxon_info <- strsplit(taxon_info,"")[[1]];
		age_break <- match(key$symbol[key$variable=="age" ],taxon_info);
		relv_taxon_info <- taxon_info[(age_break+1):length(taxon_info)];
		curlies <- min((1:length(relv_taxon_info))[relv_taxon_info %in% c(",","]")])-1;
		x <- as.numeric(paste(relv_taxon_info[1:curlies],collapse=""));
		hpd_age[nn,] <- c(x,x);
		}
	
	duration_break <- match(key$symbol[key$variable=="duration" ],taxon_info);
	if (!is.na(duration_break))	{
		relv_taxon_info <- taxon_info[duration_break:length(taxon_info)];
		curlies <- (1:length(relv_taxon_info))[relv_taxon_info %in% c("{","}")][1:2];
		hpd_dur[nn,] <- as.numeric(strsplit(paste(relv_taxon_info[(1+curlies[1]):(curlies[2]-1)],collapse=""),",")[[1]]);
		}
	}
clade_posteriors <- vector(length=nNodes);
names(clade_posteriors) <- node_names;
#nNodes <- length(clade_posteriors);
for (nn in 1:nNodes)	{
	htu <- nn+notu;
	lb <- min(left_brackets[left_brackets>clade_boundaries$rp[nn]]);
	rb <- min(right_brackets[right_brackets>lb]);
	clade_info <- molecularize2[lb:rb];
	
	age_break <- match(key$symbol[key$variable=="age"],clade_info);
	if (!is.na(age_break))	{
		relv_clade_info <- clade_info[age_break:length(clade_info)];
		curlies <- (1:length(relv_clade_info))[relv_clade_info %in% c("{","}")][1:2];
		hpd_age[htu,] <- as.numeric(strsplit(paste(relv_clade_info[(1+curlies[1]):(curlies[2]-1)],collapse=""),",")[[1]]);
		} else	{
		clade_info <- paste(clade_info,collapse="");
		clade_info <- gsub("height=",key$symbol[key$variable=="age"],clade_info);
		clade_info <- gsub("age=",key$symbol[key$variable=="age"],clade_info);
		clade_info <- strsplit(clade_info,"")[[1]];
		age_break <- match(key$symbol[key$variable=="age"],clade_info);
		curlies <- (1:length(clade_info))[clade_info %in% c(",","]")];
		curlies <- min(curlies[curlies>age_break])
		x <- as.numeric(paste(clade_info[(age_break+1):(curlies-1)],collapse=""));
		hpd_age[htu,] <- c(x,x);
		}
	
	duration_break <- match(key$symbol[key$variable=="duration" ],clade_info);
	if (!is.na(duration_break))	{
		relv_clade_info <- clade_info[duration_break:length(clade_info)];
		curlies <- (1:length(relv_clade_info))[relv_clade_info %in% c("{","}")][1:2];
		hpd_dur[htu,] <- as.numeric(strsplit(paste(relv_clade_info[(1+curlies[1]):(curlies[2]-1)],collapse=""),",")[[1]]);
		} else	{
		clade_info <- paste(clade_info,collapse="");
		clade_info <- gsub("length=",key$symbol[key$variable=="duration"],clade_info);
		clade_info <- strsplit(clade_info,"")[[1]];
		duration_break <- match(key$symbol[key$variable=="duration"],clade_info);
		if (!is.na(duration_break))	{
			curlies <- (1:length(clade_info))[clade_info %in% c(",","]")];
			curlies <- min(curlies[curlies>duration_break]);
			x <- as.numeric(paste(clade_info[(duration_break+1):(curlies-1)],collapse=""));
			} else	{
			# added 2024-05-22 to handle FBDMM output
			x <- branch_durations[htu];
			}
			hpd_dur[htu,] <- c(x,x);
		}

	posterior_break <- match(key$symbol[key$variable=="posterior"],clade_info);
	relv_clade_info <- clade_info[posterior_break:length(clade_info)];
	comma <- min((1:length(relv_clade_info))[relv_clade_info %in% ","])-1;
#	comma <- min(comma[comma>posterior_break]);
	clade_posteriors[nn] <- as.numeric(paste(relv_clade_info[2:comma],collapse=""));
#	as.numeric(strsplit(paste(relv_clade_info[(1+curlies[1]):(curlies[2]-1)],collapse=""),",")[[1]]);
	}

output <- list(newick_string_taxa_only,vector_tree,clade_posteriors,ancestral,ancestored_nodes,hpd_age,hpd_dur,branch_durations);
names(output) <- c("newick","vector_tree","clade_posteriors","ancestral","ancestored_nodes","hpd_age","hpd_dur","branch_durations");
return(output);
}

read_newick_string_mcmc_beast <- function(newick_string_full,otu_names) {
otu_names_nex <- gsub(" ","_",otu_names);
otu_names_nex[gsub("\\(","",otu_names_nex)!=otu_names_nex] <- paste("\"",otu_names[gsub("\\(","",otu_names_nex)!=otu_names_nex],"\"",sep="");
if (gsub(otu_names_nex[1],"",newick_string_full)!=newick_string_full)	{
	named <- TRUE
	} else	named <- FALSE;
simple_newick_string <- molecularize <- strsplit(newick_string_full,split="")[[1]];
left_brackets <- (1:length(molecularize))[molecularize %in% "["];
right_brackets <- (1:length(molecularize))[molecularize %in% "]"];
for (br in 1:length(left_brackets))	molecularize[left_brackets[br]:right_brackets[br]] <- "";
newick_string_taxa_only_rawwest <- newick_string_taxa_raw <- newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
notu <- length(otu_names);
branch_durations <- array(0,dim=(2*notu)-1);
for (nn in 1:notu)	{
#	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),"\\[",sep=""),"•",newick_string_taxa_only);
#	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
	if (named)	{
		dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
		} else	{
		try1 <- paste("\\(",nn,"\\[&",sep="");
		try2 <- paste(",",nn,"\\[&",sep="");
		try3 <- paste("\\(",nn,sep="");
		try4 <- paste(",",nn,sep="");
#		dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
		dummy_newick <- gsub(paste(try1,":",sep=""),"•",newick_string_taxa_only);
		dummy_newick <- gsub(paste(try2,":",sep=""),"•",dummy_newick);
		dummy_newick <- gsub(paste(try3,":",sep=""),"•",dummy_newick);
		dummy_newick <- gsub(paste(try4,":",sep=""),"•",dummy_newick);
		}
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	dd <- 1+match("•",dummy_newick);
	b_d <- dummy_newick[dd];
	dd <- dd+1;
	while (dummy_newick[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,dummy_newick[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[nn] <- as.numeric(b_d);
	}
for (i in 0:9)	newick_string_taxa_only <- gsub(i,"",newick_string_taxa_only);
newick_string_taxa_only <- gsub(":","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("-","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("\\.","",newick_string_taxa_only);
for (nn in notu:1)	{
	newick_string_taxa_only <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_only);
	newick_string_taxa_raw <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_raw);
#	newick_string_taxa_only <- gsub(paste(otu_names_nex[nn],":",sep=""),as.character(nn),newick_string_taxa_only);
#	newick_string_taxa_raw <- gsub(paste(otu_names_nex[nn],":",sep=""),as.character(nn),newick_string_taxa_raw);
	}
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa]!=";")	{
	newick_string_taxa_only_atomized <- c(newick_string_taxa_only_atomized,";");
	nstoa <- length(newick_string_taxa_only_atomized);
	}
newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
#ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
ancestral <- rep(0,notu);
ancestral[branch_durations[1:notu]==0] <- 1;
names(ancestral) <- otu_names_nex;
newick_string_taxa_only_raw <- newick_string_taxa_only;
newick_string_taxa_only <- fix_newick_ancestors(newick_string_taxa_only);
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa-1]!=")")	{
	newick_string_taxa_only_atomized[nstoa] <- ")";
	newick_string_taxa_only_atomized <- c("(",newick_string_taxa_only_atomized,";");
	newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
	}
newick_string_taxa_only <- gsub(" ","",newick_string_taxa_only);
vector_tree <- read_newick_string(newick_string=newick_string_taxa_only);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
newick_string_taxa_raw <- strsplit(newick_string_taxa_raw,split="")[[1]];
clade_ends <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ")"];
colons <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ":"];
names(clade_ends) <- length(clade_ends):1;
clade_colons <- 1+clade_ends[(clade_ends+1) %in% colons];
for (cc in 1:length(clade_colons))	{
	n_node <- notu+as.numeric(names(clade_colons)[cc]);
	dd <- clade_colons[cc]+1;
	b_d <- newick_string_taxa_raw[dd];
	dd <- dd+1;
	while (newick_string_taxa_raw[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,newick_string_taxa_raw[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[n_node] <- as.numeric(b_d);
	}
nNodes <- length(clade_ends);
newick_string_taxa_raw <- paste(newick_string_taxa_raw,collapse="");
if (nNodes<10)	{
	node_names <- paste("node_",1:nNodes,sep="")
	} else if (nNodes<100)	{
	node_names <- c(paste("node_0",1:9,sep=""),paste("node_",10:nNodes,sep=""));
	} else	{
	node_names <- c(paste("node_00",1:9,sep=""),paste("node_0",10:99,sep=""),paste("node_",100:nNodes,sep=""));
	}
names(branch_durations) <- c(otu_names_nex,node_names);

#vector_tree_raw <- read_newick_string(newick_string_taxa_only_raw);
#vector_tree <- read_newick_string(newick_string_taxa_only);
molecularize <- strsplit(newick_string_taxa_only,split="")[[1]];
if (match("(",molecularize)>1)	{
	mm <- match("(",molecularize);
	molecularize <- molecularize[mm:length(molecularize)];
	newick_string_taxa_only <- paste(molecularize,collapse="");
	}

venn_tree_newick <- transform_newick_string_to_venn_tree(newick_string = newick_string_taxa_only);

# have: newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,ancestral,branch_durations);
# need: clade_posteriors, prob_ancestor,hpd_age;
# uset mat_tree and newick_string_taxa_only to figure out which node is what number
newick_rem_info <- gsub("sampled_ancestor=","•",newick_string_full);
newick_rem_info <- gsub("height_95%_HPD","age_95%_HPD",newick_rem_info);
#newick_rem_info <- gsub("length_95%_HPD","age_95%_HPD",newick_rem_info);
if (gsub("age_95%_HPD=","§",newick_rem_info)!=newick_rem_info)	{
	newick_rem_info <- gsub("age_95%_HPD=","§",newick_rem_info);
	} else	{
	newick_rem_info <- gsub("brlen_95%_HPD=","§",newick_rem_info);
	}
newick_rem_info <- gsub("posterior=","¶",newick_rem_info);
hpd_age <- data.frame(lb=as.numeric(rep(0,notu+nNodes)),ub=as.numeric(rep(0,notu+nNodes)));
molecularized <- strsplit(newick_rem_info,"")[[1]];# print("line 2275");
molecules <- length(molecularized);
panc_boundaries <- (1:molecules)[molecularized=="•"]; # where are ancestors
post_boundaries <- (1:molecules)[molecularized=="¶"]; # where are posterior probabilities
hpd_age_boundaries <- (1:molecules)[molecularized=="§"];	# where are branch lengths
taxon_boundaries <- array(0,dim=c(notu,2)); #print("line 2280");
for (nn in notu:1)	{
	taxon_dummy <- gsub(otu_names_nex[nn],"£",newick_rem_info);
	ndummy <- length(strsplit(taxon_dummy,"")[[1]]);
	taxon_boundaries[nn,1] <- (1:ndummy)[strsplit(taxon_dummy,"")[[1]]=="£"];
#	taxon_boundaries[nn,1] <- (1:length(strsplit(taxon_dummy,"")[[1]]))[strsplit(taxon_dummy,"")[[1]]=="£"];
	taxon_boundaries[nn,2] <- taxon_boundaries[nn,1]+length(strsplit(otu_names_nex[nn],"")[[1]])-1;
	}
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_rem_info);
colnames(taxon_boundaries) <- colnames(clade_boundaries);
tu_boundaries <- rbind(taxon_boundaries,clade_boundaries);
rownames(tu_boundaries) <- all_names <- c(otu_names_nex,node_names);
tu_boundaries <- tu_boundaries[order(tu_boundaries$rp),];# print("line 2290");
brackets_l <- (1:molecules)[molecularized=="{"];
brackets_r <- (1:molecules)[molecularized=="}"];
hp <- 0;
rownames(hpd_age) <- all_names;
while (hp < length(hpd_age_boundaries))	{
	hp <- hp+1;
	tx <- sum(tu_boundaries$rp<hpd_age_boundaries[hp]);
	txn <- match(rownames(tu_boundaries)[tx],all_names);
	i <- brackets_l[1+sum(hpd_age_boundaries[hp]>brackets_l)]+1;
	j <- brackets_r[1+sum(hpd_age_boundaries[hp]>brackets_r)]-1;
	hpd_age[txn,] <- as.numeric(strsplit(paste(molecularized[i:j],collapse=""),",")[[1]]);
	}
hpd_age[vector_tree[(1:notu)[ancestral==1]],] <- hpd_age[(1:notu)[ancestral==1],];

rownames(taxon_boundaries) <- otu_names_nex;#  print("line 2305");
taxon_boundaries <- taxon_boundaries[order(taxon_boundaries[,1]),];
brackets_l <- (1:molecules)[molecularized=="["];
brackets_r <- (1:molecules)[molecularized=="]"];
prob_ancestor <- array(0,dim=notu);
pb <- 0;  #print("line 2310");
while (pb < length(panc_boundaries))	{
#for (pb in 1:length(panc_boundaries))	{
	pb <- pb+1;
	tx <- sum(taxon_boundaries[,1]<panc_boundaries[pb]);
	txn <- match(rownames(taxon_boundaries)[tx],otu_names_nex);
	i <- panc_boundaries[pb]+1;
	pranc <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		pranc <- paste(pranc,molecularized[i],sep="");
		i <- i+1;
		}
	prob_ancestor[txn] <- as.numeric(pranc);
	}
names(prob_ancestor) <- otu_names_nex;
ancestral[prob_ancestor>=0.5] <- 1; # added 2023-09-25: we seem to miss some ancestors!

clade_posteriors <- array(0,dim=nNodes);#  print("line 2325");
rownames(clade_boundaries) <- names(clade_posteriors) <- node_names;
clade_boundaries <- clade_boundaries[order(clade_boundaries$rp),];
pp <- 0;
while (pp < length(post_boundaries))	{
	pp <- pp+1;
	cl <- sum(clade_boundaries$rp<post_boundaries[pp]);
	cld <- match(rownames(clade_boundaries)[cl],node_names);
	i <- post_boundaries[pp]+1;
	postp <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		postp <- paste(postp,molecularized[i],sep="");
		i <- i+1;
		}
	clade_posteriors[cld] <- as.numeric(postp);
#	molecularized[i:i+-1:10]
	}
#sum(strsplit(newick_rem_info,"")[[1]]=="•")
#sum(strsplit(newick_rem_info,"")[[1]]=="§")
#sum(strsplit(newick_rem_info,"")[[1]]=="¶")
output <- list(newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,clade_posteriors,ancestral,prob_ancestor,hpd_age,branch_durations);
names(output) <- c("newick_modified","newick","vector_tree","clade_posteriors","ancestral","prob_ancestor","hpd_age","branch_durations");
return(output);
}

# for (taxon1:branch_lenght)
read_newick_string_mcmc_simple <- function(newick_string_full,otu_names) {
otu_names_nex <- gsub(" ","_",otu_names);
otu_names_nex[gsub("\\(","",otu_names_nex)!=otu_names_nex] <- paste("\"",otu_names[gsub("\\(","",otu_names_nex)!=otu_names_nex],"\"",sep="");
if (gsub(otu_names_nex[1],"",newick_string_full)!=newick_string_full)	{
	named <- TRUE
	} else	named <- FALSE;
simple_newick_string <- molecularize <- strsplit(newick_string_full,split="")[[1]];
left_brackets <- (1:length(molecularize))[molecularize %in% "["];
right_brackets <- (1:length(molecularize))[molecularize %in% "]"];
for (br in 1:length(left_brackets))	molecularize[left_brackets[br]:right_brackets[br]] <- "";
newick_string_taxa_only_rawwest <- newick_string_taxa_raw <- newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
notu <- length(otu_names);
branch_durations <- array(0,dim=(2*notu)-1);
for (nn in 1:notu)	{
#	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),"\\[",sep=""),"•",newick_string_taxa_only);
#	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
	if (named)	{
		dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
		} else	{
		try1 <- paste("\\(",nn,"\\[&",sep="");
		try2 <- paste(",",nn,"\\[&",sep="");
		dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
		}
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	dd <- 1+match("•",dummy_newick);
	b_d <- dummy_newick[dd];
	dd <- dd+1;
	while (dummy_newick[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,dummy_newick[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[nn] <- as.numeric(b_d);
	}
for (i in 0:9)	newick_string_taxa_only <- gsub(i,"",newick_string_taxa_only);
newick_string_taxa_only <- gsub(":","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("-","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("\\.","",newick_string_taxa_only);
for (nn in notu:1)	{
	newick_string_taxa_only <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_only);
	newick_string_taxa_raw <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_raw);
#	newick_string_taxa_only <- gsub(paste(otu_names_nex[nn],":",sep=""),as.character(nn),newick_string_taxa_only);
#	newick_string_taxa_raw <- gsub(paste(otu_names_nex[nn],":",sep=""),as.character(nn),newick_string_taxa_raw);
	}
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa]!=";")	{
	newick_string_taxa_only_atomized <- c(newick_string_taxa_only_atomized,";");
	nstoa <- length(newick_string_taxa_only_atomized);
	}
newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
#ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
ancestral <- rep(0,notu);
ancestral[branch_durations[1:notu]==0] <- 1;
print("Naming ancestors");
names(ancestral) <- otu_names_nex;
newick_string_taxa_only_raw <- newick_string_taxa_only;
newick_string_taxa_only <- fix_newick_ancestors(newick_string_taxa_only);
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa-1]!=")")	{
	newick_string_taxa_only_atomized[nstoa] <- ")";
	newick_string_taxa_only_atomized <- c("(",newick_string_taxa_only_atomized,";");
	newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
	}
newick_string_taxa_only <- gsub(" ","",newick_string_taxa_only);
vector_tree <- read_newick_string(newick_string=newick_string_taxa_only);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
newick_string_taxa_raw <- strsplit(newick_string_taxa_raw,split="")[[1]];
clade_ends <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ")"];
colons <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ":"];
names(clade_ends) <- length(clade_ends):1;
clade_colons <- 1+clade_ends[(clade_ends+1) %in% colons];
for (cc in 1:length(clade_colons))	{
	n_node <- notu+as.numeric(names(clade_colons)[cc]);
	dd <- clade_colons[cc]+1;
	b_d <- newick_string_taxa_raw[dd];
	dd <- dd+1;
	while (newick_string_taxa_raw[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,newick_string_taxa_raw[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[n_node] <- as.numeric(b_d);
	}
nNodes <- length(clade_ends);
newick_string_taxa_raw <- paste(newick_string_taxa_raw,collapse="");
if (nNodes<10)	{
	node_names <- paste("node_",1:nNodes,sep="")
	} else if (nNodes<100)	{
	node_names <- c(paste("node_0",1:9,sep=""),paste("node_",10:nNodes,sep=""));
	} else	{
	node_names <- c(paste("node_00",1:9,sep=""),paste("node_0",10:99,sep=""),paste("node_",100:nNodes,sep=""));
	}
names(branch_durations) <- c(otu_names_nex,node_names);

#vector_tree_raw <- read_newick_string(newick_string_taxa_only_raw);
#vector_tree <- read_newick_string(newick_string_taxa_only);
molecularize <- strsplit(newick_string_taxa_only,split="")[[1]];
if (match("(",molecularize)>1)	{
	mm <- match("(",molecularize);
	molecularize <- molecularize[mm:length(molecularize)];
	newick_string_taxa_only <- paste(molecularize,collapse="");
	}

venn_tree_newick <- transform_newick_string_to_venn_tree(newick_string = newick_string_taxa_only);

# have: newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,ancestral,branch_durations);
# need: clade_posteriors, prob_ancestor,hpd_age;
# uset mat_tree and newick_string_taxa_only to figure out which node is what number
newick_rem_info <- gsub("sampled_ancestor=","•",newick_string_full);
newick_rem_info <- gsub("height_95%_HPD","age_95%_HPD",newick_rem_info);
#newick_rem_info <- gsub("length_95%_HPD","age_95%_HPD",newick_rem_info);
if (gsub("age_95%_HPD=","§",newick_rem_info)!=newick_rem_info)	{
	newick_rem_info <- gsub("age_95%_HPD=","§",newick_rem_info);
	} else	{
	newick_rem_info <- gsub("brlen_95%_HPD=","§",newick_rem_info);
	}
newick_rem_info <- gsub("posterior=","¶",newick_rem_info);
hpd_age <- data.frame(lb=as.numeric(rep(0,notu+nNodes)),ub=as.numeric(rep(0,notu+nNodes)));
molecularized <- strsplit(newick_rem_info,"")[[1]];# print("line 2275");
molecules <- length(molecularized);
panc_boundaries <- (1:molecules)[molecularized=="•"]; # where are ancestors
post_boundaries <- (1:molecules)[molecularized=="¶"]; # where are posterior probabilities
hpd_age_boundaries <- (1:molecules)[molecularized=="§"];	# where are branch lengths
taxon_boundaries <- array(0,dim=c(notu,2)); #print("line 2280");
for (nn in notu:1)	{
	taxon_dummy <- gsub(otu_names_nex[nn],"£",newick_rem_info);
	ndummy <- length(strsplit(taxon_dummy,"")[[1]]);
	taxon_boundaries[nn,1] <- (1:ndummy)[strsplit(taxon_dummy,"")[[1]]=="£"];
#	taxon_boundaries[nn,1] <- (1:length(strsplit(taxon_dummy,"")[[1]]))[strsplit(taxon_dummy,"")[[1]]=="£"];
	taxon_boundaries[nn,2] <- taxon_boundaries[nn,1]+length(strsplit(otu_names_nex[nn],"")[[1]])-1;
	}
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_rem_info);
colnames(taxon_boundaries) <- colnames(clade_boundaries);
tu_boundaries <- rbind(taxon_boundaries,clade_boundaries);
rownames(tu_boundaries) <- all_names <- c(otu_names_nex,node_names);
tu_boundaries <- tu_boundaries[order(tu_boundaries$rp),];# print("line 2290");
brackets_l <- (1:molecules)[molecularized=="{"];
brackets_r <- (1:molecules)[molecularized=="}"];
hp <- 0;
rownames(hpd_age) <- all_names;
while (hp < length(hpd_age_boundaries))	{
	hp <- hp+1;
	tx <- sum(tu_boundaries$rp<hpd_age_boundaries[hp]);
	txn <- match(rownames(tu_boundaries)[tx],all_names);
	i <- brackets_l[1+sum(hpd_age_boundaries[hp]>brackets_l)]+1;
	j <- brackets_r[1+sum(hpd_age_boundaries[hp]>brackets_r)]-1;
	hpd_age[txn,] <- as.numeric(strsplit(paste(molecularized[i:j],collapse=""),",")[[1]]);
	}
hpd_age[vector_tree[(1:notu)[ancestral==1]],] <- hpd_age[(1:notu)[ancestral==1],];

rownames(taxon_boundaries) <- otu_names_nex;#  print("line 2305");
taxon_boundaries <- taxon_boundaries[order(taxon_boundaries[,1]),];
brackets_l <- (1:molecules)[molecularized=="["];
brackets_r <- (1:molecules)[molecularized=="]"];
prob_ancestor <- array(0,dim=notu);
pb <- 0;  #print("line 2310");
while (pb < length(panc_boundaries))	{
#for (pb in 1:length(panc_boundaries))	{
	pb <- pb+1;
	tx <- sum(taxon_boundaries[,1]<panc_boundaries[pb]);
	txn <- match(rownames(taxon_boundaries)[tx],otu_names_nex);
	i <- panc_boundaries[pb]+1;
	pranc <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		pranc <- paste(pranc,molecularized[i],sep="");
		i <- i+1;
		}
	prob_ancestor[txn] <- as.numeric(pranc);
	}
names(prob_ancestor) <- otu_names_nex;
ancestral[prob_ancestor>=0.5] <- 1; # added 2023-09-25: we seem to miss some ancestors!

clade_posteriors <- array(0,dim=nNodes);#  print("line 2325");
rownames(clade_boundaries) <- names(clade_posteriors) <- node_names;
clade_boundaries <- clade_boundaries[order(clade_boundaries$rp),];
pp <- 0;
while (pp < length(post_boundaries))	{
	pp <- pp+1;
	cl <- sum(clade_boundaries$rp<post_boundaries[pp]);
	cld <- match(rownames(clade_boundaries)[cl],node_names);
	i <- post_boundaries[pp]+1;
	postp <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		postp <- paste(postp,molecularized[i],sep="");
		i <- i+1;
		}
	clade_posteriors[cld] <- as.numeric(postp);
#	molecularized[i:i+-1:10]
	}
#sum(strsplit(newick_rem_info,"")[[1]]=="•")
#sum(strsplit(newick_rem_info,"")[[1]]=="§")
#sum(strsplit(newick_rem_info,"")[[1]]=="¶")
output <- list(newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,clade_posteriors,ancestral,prob_ancestor,hpd_age,branch_durations);
names(output) <- c("newick_modified","newick","vector_tree","clade_posteriors","ancestral","prob_ancestor","hpd_age","branch_durations");
return(output);
}

# character_rates <- branch_rates[1,]
# branch_durations <- branch_durations[1,]
# taxon_info <- fossil_fas;
#(1:length(character_rates))[as.numeric(character_rates[1,]) %in% as.numeric(character_rates[1,42])]
accersi_rev_bayes_index_numbers_for_otus <- function(newick_string_full,otu_names)	{
notu <- length(otu_names);
ancestors <- identify_newick_ancestors(newick_string_full,otu_names);
output <- data.frame(taxon=as.character(otu_names),index=as.numeric(1:notu));
for (i in 1:notu)	{
	flibberty <- strsplit(newick_string_full,otu_names[i])[[1]][2]
	jibberty <- strsplit(flibberty,"\\]")[[1]][1];
	if (!otu_names[i] %in% ancestors)	{
#		flibberty <- strsplit(newick_string_full,otu_names[i])[[1]][2]
#		jibberty <- strsplit(flibberty,"\\]")[[1]][1];
		output$index[i] <- as.numeric(gsub("\\[&index=","",jibberty));
		} else	{
#		output$index[i] <- 0;
		output$index[i] <- -as.numeric(gsub("\\[&index=","",jibberty));
		}
	}
return(output);
}

reduce_indexed_newick_tree_to_only_indices <- function(newick_string_full,otu_names)	{
# first eliminate branch length information
#basics <- read_newick_string_mcmc(newick_string_full,otu_names);
#ancestors <- identify_newick_ancestors(newick_string_full,otu_names);
notu <- length(otu_names);
#matrix_tree <- transform_vector_tree_to_matrix_tree(read_newick_string(basics$newick_modified));
#otu_indices <- accersi_rev_bayes_index_numbers_for_otus(newick_string_full,otu_names);

molecularized_newick <- strsplit(newick_string_full,":")[[1]];
for (i in 2:length(molecularized_newick))	{
	atomized <- strsplit(molecularized_newick[i],"")[[1]];
	while (atomized[1] %in% c(0:9,"."))	atomized <- atomized[2:length(atomized)];
	molecularized_newick[i] <- paste(atomized,collapse="");
	}

reduced_newick <- paste(molecularized_newick,collapse="");
#sampled_ancestors <- identify_newick_ancestors(best_tree,otu_names);
for (i in 1:notu)	reduced_newick <- gsub(otu_names[i],"",reduced_newick);
return(reduced_newick);
}

read_newick_string_with_revbayes_indices <- function(newick_string_full,otu_names) {
basics <- read_newick_string_mcmc(newick_string_full,otu_names);
ttu <- max(basics$vector_tree);
notu <- length(otu_names);
htus <- (notu+1):ttu;
reduced_newick <- reduce_indexed_newick_tree_to_only_indices(newick_string_full,otu_names);
matrix_tree <- transform_vector_tree_to_matrix_tree(read_newick_string(basics$newick_modified));
nnodes <- nrow(matrix_tree);
otu_indices <- accersi_rev_bayes_index_numbers_for_otus(newick_string_full,otu_names);
#otu_index <- paste("\\[&index=",sort(otu_indices$index[otu_indices$index>0]),"\\]",sep="");
otu_index <- paste("\\[&index=",sort(otu_indices$index),"\\]",sep="");
node_index <- paste("\\[&index=",(notu+1):ttu,"\\]",sep="");
for (i in 1:length(otu_index))  reduced_newick <- gsub(otu_index[i],i,reduced_newick);
for (i in 1:nnodes)  {
	htu <- notu+i;
	reduced_newick <- gsub(node_index[i],htu,reduced_newick);
	}
reduced_newick_fixed <- fix_newick_ancestors(reduced_newick);
full_matrix_tree <- transform_vector_tree_to_matrix_tree(read_newick_string(reduced_newick_fixed))
full_venn_tree <- transform_matrix_tree_to_venn_tree(full_matrix_tree);
redone_matrix_tree <- node_labels <- c();
for (h in 1:nnodes)	{
	nn <- which(full_matrix_tree==htus[h],arr.ind=TRUE)[,1];
	redone_matrix_tree <- rbind(redone_matrix_tree,full_matrix_tree[nn,!full_matrix_tree[nn,] %in% htus[h]]);
	node_labels <- c(htus[h],node_labels);
	}
rownames(redone_matrix_tree) <- node_labels;
redone_matrix_tree[,1][redone_matrix_tree[,1]>ttu] <- redone_matrix_tree[,1][redone_matrix_tree[,1]>ttu]-nnodes;
redone_matrix_tree[,2][redone_matrix_tree[,2]>ttu] <- redone_matrix_tree[,2][redone_matrix_tree[,2]>ttu]-nnodes;

nodes <- 0;
if (length(reduced_newick)==1)	reduced_newick <- strsplit(reduced_newick,split="",fixed=TRUE)[[1]];
for (i in 1:length(reduced_newick))		if (reduced_newick[i]=="(")	nodes <- nodes+1;
# get clades
clades <- vector(length=nodes);
for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0
for (i in 1:length(reduced_newick))	{
	if (reduced_newick[i]>="0" && reduced_newick[i]<="9")	{
		otu <- as.numeric(reduced_newick[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(reduced_newick))	{
	if (reduced_newick[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (reduced_newick[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (reduced_newick[i]==",")	{
		sp <- p <- 0
		} else if (reduced_newick[i]>="0" && reduced_newick[i]<="9")	{
		sp <- as.numeric(reduced_newick[i])+(sp*10)
		p <- p+1
		if (reduced_newick[i+1]<"0" || reduced_newick[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}


why_must_life_be_hard <- function(newick_string_full,otu_names)	{
# first eliminate branch length information
basics <- read_newick_string_mcmc(newick_string_full,otu_names);
ancestors <- identify_newick_ancestors(newick_string_full,otu_names);
notu <- length(otu_names);
matrix_tree <- transform_vector_tree_to_matrix_tree(read_newick_string(basics$newick_modified));
otu_indices <- accersi_rev_bayes_index_numbers_for_otus(newick_string_full,otu_names);

molecularized_newick <- strsplit(newick_string_full,":")[[1]];
for (i in 2:length(molecularized_newick))	{
	atomized <- strsplit(molecularized_newick[i],"")[[1]];
	while (atomized[1] %in% c(0:9,"."))	atomized <- atomized[2:length(atomized)];
	molecularized_newick[i] <- paste(atomized,collapse="");
	}
notu <- length(otu_names);
reduced_newick <- paste(molecularized_newick,collapse="");
sampled_ancestors <- identify_newick_ancestors(best_tree,otu_names);
for (i in 1:notu)	reduced_newick <- gsub(otu_names[i],"",reduced_newick);
indices <- paste("\\[&index=",1:((2*notu)-1),"\\]",sep="");
node_index <- index_node <- vector(length=length(indices));
rownames(matrix_tree) <- -1:-nnodes;
for (i in 6:notu)	{
	notu_node <- which(matrix_tree==i,arr.ind=TRUE)[,1];
  this_index <- paste("\\[&index=",otu_indices$index[i],"\\]",sep="");
#  if (gsub(indices[i],"",reduced_newick)!=reduced_newick)	{
  if (gsub(this_index,"",reduced_newick)!=reduced_newick)	{
#		flibberty <- gsub(indices[i],"•",reduced_newick);
		flibberty <- gsub(this_index,"•",reduced_newick);
		jibberty <- strsplit(flibberty,"•")[[1]][2];
		blibberty <- strsplit(jibberty,"\\)")[[1]];
		blibberty <- blibberty[blibberty!=""];
		blabberty <- strsplit(blibberty[1],"]")[[1]][1];
		node_index[notu_node] <- rownames(matrix_tree)[notu_node] <- as.numeric(gsub("\\[&index=","",blabberty));
		}
	}
}

accersi_rev_bayes_index_numbers_for_nodes <- function(newick_string_full,otu_names)	{
basics <- read_newick_string_mcmc(newick_string_full,otu_names);
ancestors <- identify_newick_ancestors(newick_string_full,otu_names);
notu <- length(otu_names);
matrix_tree <- transform_vector_tree_to_matrix_tree(read_newick_string(basics$newick_modified));
nnodes <- nrow(matrix_tree);
otu_indices <- accersi_rev_bayes_index_numbers_for_otus(newick_string_full,otu_names);
reduced_newick <- newick_string_full;
for (i in 1:nrow(otu_indices))	{
	reduced_newick <- gsub(otu_indices$taxon[i],otu_indices$index[i],reduced_newick);
	dummy <- paste("\\[&index=",otu_indices$index[i],"\\]",sep="");
	reduced_newick <- gsub(dummy,"",reduced_newick);
	}

#molecularized_reduced_newick <- strsplit(reduced_newick,":")[[1]];
molecularized_reduced_newick <- strsplit(reduced_newick,"\\)")[[1]];
node_indices <- vector(length=nnodes);
nn <- nnodes;
node_names <- c();
for (j in 2:length(molecularized_reduced_newick))	{
	flibberty <- molecularized_reduced_newick[j];
	jibberty <- strsplit(flibberty,":")[[1]][1]
	this_index <- gsub("\\[&index=","",jibberty);
	node_indices[nn] <- as.numeric(gsub("\\]","",this_index));
	zeros <- floor(log10(nnodes))-floor(log10(nn));
	node_names <- c(paste("node_",rep(0,zeros),nn,sep=""),node_names);
	nn <- nn-1;
	}
#for (j in 1:length(molecularized_reduced_newick))	{
#	atomized_newick_part <- strsplit(molecularized_reduced_newick[j],"")[[1]];
#	while (atomized_newick_part[1] %in% c(0:9,".") && length(atomized_newick_part)>0)	atomized_newick_part <- atomized_newick_part[2:length(atomized_newick_part)];
#	molecularized_reduced_newick[j] <- paste(atomized_newick_part,collapse="");
#	}
#reduced_newick <- paste(molecularized_reduced_newick,collapse="");

output <- data.frame(ttu=as.numeric((notu+1):(notu+nnodes)),
										 taxon=as.character(node_names),
										 index=as.numeric(node_indices));
#molecularized_reduced_newick <- strsplit(reduced_newick,"\\)")[[1]];
#molecularized_reduced_newick <- molecularized_reduced_newick[!gsub("index","",molecularized_reduced_newick) %in% molecularized_reduced_newick];
#for (nn in 1:nnodes)	{
#	node_no <- nnodes+1-nn;
#	flibberty <- strsplit(reduced_newick,"\\)")[[1]][nn+1]
#	flibberty <- molecularized_reduced_newick[nn]
#	if (gsub("index","",flibberty)==flibberty)	{
#		flibberty <- strsplit(reduced_newick,"\\)")[[1]][i+1];
#		i <- i+1;
#		}
#	jibberty <- strsplit(flibberty,"\\]")[[1]][1];
#	output$index[node_no] <- as.numeric(gsub("\\[&index=","",jibberty));
#	nn <- nn+1;
#	print(output[node_no:nnodes,])
#	}
return(output);
}

accersi_rev_bayes_index_numbers_for_otus_and_nodes <- function(newick_string_full,otu_names)	{
basics <- read_newick_string_mcmc(newick_string_full,otu_names);
notu <- length(otu_names);
vector_tree <- read_newick_string(basics$newick_modified);

output1 <- accersi_rev_bayes_index_numbers_for_otus(newick_string_full,otu_names);
output1$index[output1$index==0] <- (1:notu)[!(1:notu) %in% output1$index[output1$index!=0]];
output2 <- accersi_rev_bayes_index_numbers_for_nodes(newick_string_full,otu_names);
output1$ttu <- 1:nrow(output1);
output1 <- output1[,match(colnames(output2),colnames(output1))];
output <- rbind(output1,output2);
return(output);
}

link_mcmc_newick_string_to_rates <- function(newick_string_full,character_rates,branch_durations,divergence_dates,taxon_info) {
notu <- nrow(taxon_info);
#colnames(character_rates) <- gsub("branch_rates","",colnames(character_rates))
redone_newick <- gsub("\\[&","•",newick_string_full);
redone_newick <- gsub("\\]:","•",redone_newick);
for (i in 1:notu)	redone_newick <- gsub(taxon_info$taxon[i],paste("•",taxon_info$taxon[i],sep=""),redone_newick);
redone_newick <- gsub("\\(","",redone_newick);
redone_newick <- gsub("\\)","",redone_newick);
redone_newick <- gsub(",","",redone_newick);
redone_newick <- gsub(";","",redone_newick);
redone_newick <- gsub("index=","",redone_newick);
redone_newick_molecularized <- strsplit(redone_newick,"•")[[1]];
redone_newick_molecularized <- redone_newick_molecularized[redone_newick_molecularized!=""];
r_n_m <- length(redone_newick_molecularized);
tree_rate_info <- data.frame(taxon=as.character(),index=as.numeric(),branch_duration=as.numeric(),character_rate=as.numeric());
i <- 0;
dummy <- data.frame(taxon=as.character("a"),index=as.numeric(0),branch_duration=as.numeric(0),character_rate=as.numeric(0));
ancestral_branches <- (1:notu)[branch_durations[1,]<=0];
while (i < r_n_m)	{
	i <- i+1;
#for (i in 1:r_n_m)	{
	if (redone_newick_molecularized[i] %in% taxon_info$taxon)	{
		dummy$taxon <- redone_newick_molecularized[i];
		} else if (dummy$index[1]==0 && (as.numeric(redone_newick_molecularized[i]))==round((as.numeric(redone_newick_molecularized[i])),0))	{
		dummy$index <- round(as.numeric(redone_newick_molecularized[i]),0);
		if (dummy$index<length(character_rates))
			dummy$character_rate <- as.numeric(character_rates[dummy$index]);
		if (dummy$index>notu)	{
			zeros <- floor(log10(notu))-floor(log10(dummy$index-notu));
			if (dummy$taxon=="a")	{
				dummy$taxon <- paste("node_",rep(0,zeros),dummy$index-notu,sep="");
				} else	{
				dummy <- rbind(dummy,dummy);
				dummy$taxon[2] <- paste("node_",rep(0,zeros),dummy$index[2]-notu,sep="");
				txn <- match(dummy$taxon[1],taxon_info$taxon);
				dummy$index[1] <- ancestral_branches[divergence_dates[ancestral_branches] %in% taxon_info$min_age[txn]][1];
				ancestral_branches <- ancestral_branches[!ancestral_branches %in% dummy$taxon[1]];
				}
			}
		} else	{
		dummy$branch_duration <- as.numeric(redone_newick_molecularized[i]);
		tree_rate_info <- rbind(tree_rate_info,dummy);
#		print(tree_rate_info);
		dummy <- data.frame(taxon=as.character("a"),index=as.numeric(0),branch_duration=as.numeric(0),character_rate=as.numeric(0));
		}
	}
#tree_rate_info[order(tree_rate_info$index),]
htu <- max(tree_rate_info$index);
branch_durations[branch_durations==0] <- -1;
rosetta <- match(round(tree_rate_info$branch_duration,4),round(branch_durations,4));
#tree_rate_info$branch_duration[is.na(rosetta)]
if (sum(is.na(rosetta))>0)	rosetta[is.na(rosetta)] <- match(round(tree_rate_info$branch_duration,3),round(branch_durations,3))[is.na(rosetta)];
if (sum(is.na(rosetta))>0)	rosetta[is.na(rosetta)] <- match(round(tree_rate_info$branch_duration,2),round(branch_durations,2))[is.na(rosetta)];
rosetta[tree_rate_info$branch_duration==0] <- htu;
tree_rate_info$onset_date <- tree_rate_info$divergence_dates <- rep(0,nrow(tree_rate_info));
tree_rate_info$divergence_dates[tree_rate_info$index<htu] <- as.numeric(divergence_dates[rosetta[rosetta<htu]]);
tree_rate_info$onset_date[tree_rate_info$index<htu] <- tree_rate_info$divergence_dates[tree_rate_info$index<htu]+tree_rate_info$branch_duration[tree_rate_info$index<htu];
tree_rate_info$mid_point <- (tree_rate_info$onset_date+tree_rate_info$divergence_dates)/2;

setup_tree <- newick_string_full;
for (i in 1:nrow(tree_rate_info))	setup_tree <- gsub(paste(":",tree_rate_info$branch_duration[i],sep=""),"",setup_tree);
setup_tree <- gsub(000,"",setup_tree);
setup_tree <- gsub(00,"",setup_tree);
for (i in 1:max(tree_rate_info$index))	setup_tree <- gsub(paste("\\[&index=",i,"\\]",sep=""),"",setup_tree);
setup_tree <- gsub("\\.","",setup_tree);
setup_tree <- gsub(";","",setup_tree);
i <- 0;
while (i < notu)	{
	i <- i+1;
	nn <- match(taxon_info$taxon[i],tree_rate_info$taxon);
	setup_tree <- gsub(taxon_info$taxon[i],tree_rate_info$index[nn],setup_tree);
#	print(setup_tree);
	}

#setup_tree_old <- setup_tree;
setup_tree <- fix_newick_ancestors(setup_tree)
vector_tree <- read_newick_string(setup_tree);

i <- 0;
dummy <- data.frame(taxon=as.character("a"),index=as.numeric(0),branch_duration=as.numeric(0),character_rate=as.numeric(0));
ancestral_branches <- (1:notu)[branch_durations[1,]<=0];
while (i < r_n_m)	{
	i <- i+1;
#for (i in 1:r_n_m)	{
	if (redone_newick_molecularized[i] %in% taxon_info$taxon)	{
		dummy$taxon <- redone_newick_molecularized[i];
		} else if (dummy$index[1]==0 && (as.numeric(redone_newick_molecularized[i]))==round((as.numeric(redone_newick_molecularized[i])),0))	{
		dummy$index <- round(as.numeric(redone_newick_molecularized[i]),0);
		if (dummy$index<length(character_rates))
			dummy$character_rate <- as.numeric(character_rates[dummy$index]);
		if (dummy$index>notu)	{
			zeros <- floor(log10(notu))-floor(log10(dummy$index-notu));
			if (dummy$taxon=="a")	{
				dummy$taxon <- paste("node_",rep(0,zeros),dummy$index-notu,sep="");
				} else	{
				dummy <- rbind(dummy,dummy);
				dummy$taxon[2] <- paste("node_",rep(0,zeros),dummy$index[2]-notu,sep="");
				txn <- match(dummy$taxon[1],taxon_info$taxon);
				dummy$index[1] <- ancestral_branches[divergence_dates[ancestral_branches] %in% taxon_info$min_age[txn]][1];
				ancestral_branches <- ancestral_branches[!ancestral_branches %in% dummy$taxon[1]];
				}
			}
		} else	{
		dummy$branch_duration <- as.numeric(redone_newick_molecularized[i]);
		tree_rate_info <- rbind(tree_rate_info,dummy);
#		print(tree_rate_info);
		dummy <- data.frame(taxon=as.character("a"),index=as.numeric(0),branch_duration=as.numeric(0),character_rate=as.numeric(0));
		}
	}
#tree_rate_info[order(tree_rate_info$index),]

dd <- max(hist(vector_tree[vector_tree>0],breaks=(min(vector_tree[vector_tree>0])-1):max(vector_tree),plot=F)$counts);
nnodes <- max(tree_rate_info$index)-notu;
matrix_tree_indexed <- transform_vector_tree_to_matrix_tree(vector_tree);
rownames(matrix_tree_indexed) <- rep(0,nnodes);
for (i in 1:nnodes)	{
	if (matrix_tree_indexed[i,1]<=notu)	{
		nn <- match(matrix_tree_indexed[i,1],tree_rate_info$index)
		} else	{
		zeros <- floor(log10(notu))-floor(log10(matrix_tree_indexed[i,1]-notu));
		jj <- paste("node_",rep(0,zeros),matrix_tree_indexed[i,1]-notu,sep="");
		nn <- match(jj,tree_rate_info$taxon)
		}
	ii <- tree_rate_info$index[(nn+1):nrow(tree_rate_info)][tree_rate_info$index[(nn+1):nrow(tree_rate_info)]>notu][1]
	rownames(matrix_tree_indexed)[i] <- ii;
	}

ancestral_node <- vector(length=htu);
matrix_tree_index_string <- matrix_tree_indexed;
for (i in 1:nnodes)	{
	for (j in 1:ncol(matrix_tree_indexed))	{
		if (matrix_tree_indexed[i,j]>0)
			matrix_tree_index_string[i,j] <- tree_rate_info$taxon[match(matrix_tree_index_string[i,j],tree_rate_info$index)];
		}
	}

output <- list(tree_rate_info,matrix_tree_indexed,matrix_tree_index_string)
names(output) <- c("branch_rate_information","matrix_tree","indexed_matrix_tree");
return(output);
}

read_newick_string_mcmc_indexed <- function(newick_string_full,otu_names)	{
notu <- length(otu_names);
atomized_newick <- strsplit(newick_string_full,"")[[1]];
nnodes <- sum(atomized_newick==")");
htu <- nnodes+notu;
get_basics <- read_newick_string_mcmc(newick_string_full,otu_names);
vector_tree <- read_newick_string(get_basics$newick_modified);
matrix_tree <- transform_vector_tree_to_matrix_tree(vector_tree);

indexed_taxa <- c();
for (i in 1:notu)	{
	altered_newick <- strsplit(newick_string_full,otu_names[i])[[1]];
	this_index <- strsplit(altered_newick[2],"\\]")[[1]][1];
	this_index <- gsub("\\[&","",this_index);
	this_index_no <- as.numeric(strsplit(this_index,"=")[[1]][2]);
	zeros <- floor(log10(htu))-floor(log10(this_index_no));
	this_index <- paste("index_",rep(0,zeros),this_index_no,sep="");
	node_no <- as.numeric(which(matrix_tree==i,arr.ind = TRUE)[,1]+notu);
	indexed_taxa <- rbind(indexed_taxa,c(otu_names[i],this_index,node_no));
	}
indexed_nodes <- c();
for (i in nnodes:1)	{
	j <- nnodes+1-i;
	htu <- i+notu;
	altered_newick <- strsplit(newick_string_full,")")[[1]][j+1];
	this_index <- strsplit(altered_newick,"\\]")[[1]][1];
	this_index <- gsub("\\[&","",this_index);
	this_index_no <- as.numeric(strsplit(this_index,"=")[[1]][2]);
	zeros <- floor(log10(htu))-floor(log10(this_index_no));
	this_index <- paste("index_",rep(0,zeros),this_index_no,sep="");
	zeros2 <- floor(log10(htu))-floor(log10(i));
	node_name <- paste("node_",rep(0,zeros),i,sep="");
	node_no <- as.numeric(which(matrix_tree==htu,arr.ind = TRUE)[,1]+notu);
	if (length(node_no)==0 || is.null(node_no) || is.na(node_no))	node_no <- -1;
	indexed_nodes <- rbind(c(node_name,this_index,node_no),indexed_nodes);
	}

indexed_taxa <- data.frame(rbind(indexed_taxa,indexed_nodes));
colnames(indexed_taxa) <- c("taxonomic_unit","index","ancestral_node")
indexed_taxa$index_no <- as.numeric(gsub("index_","",indexed_taxa$index));
return(indexed_taxa);
}
 
#beast_tree <- mctree1;
fix_beast_trees <- function(beast_tree)	{
beast_tree <- gsub("\t","",beast_tree);
begin_taxa <- match("begin taxa;",tolower(beast_tree));
notu <- as.numeric(strsplit(gsub(";","",beast_tree[begin_taxa+1]),split = "=")[[1]][2]);
taxon_st <- match("taxlabels",tolower(beast_tree))+1;
taxon_en <- match("taxlabels",tolower(beast_tree))+notu;
otu_names_beast <- beast_tree[taxon_st:taxon_en];
otu_names_revby <- sort(otu_names_beast);
if (beast_tree[length(beast_tree)]=="End;")
	rb_tree <- beast_tree[length(beast_tree)-1];	
rb_tree <- gsub(" \\(","\\(",rb_tree);

for (sp in notu:1)	{
	molecularized <- strsplit(rb_tree,paste(sp,"\\[",sep=""))[[1]];
	rb_tree <- paste(molecularized[1],paste(otu_names_beast[sp],"[",sep=""),molecularized[2],sep="");
	}
rb_output <- "#NEXUS";
rb_output <- c(rb_output,"");
rb_output <- c(rb_output,"\tBegin taxa;");
rb_output <- c(rb_output,paste("\tDimensions ntax=",notu,";",sep=""));
rb_output <- c(rb_output,"\tTaxlabels");

for (sp in 1:notu)
	rb_output <- c(rb_output,paste("\t\t",otu_names_revby[sp],sep=""));
rb_output <- c(rb_output,"\t\t;");
rb_output <- c(rb_output,"End;");
rb_output <- c(rb_output,"");
rb_output <- c(rb_output,"Begin trees;");
rb_output <- c(rb_output,rb_tree);
rb_output <- c(rb_output,"End;");
return(rb_output);
}

ladderize_newick_string <- function(newick_string)	{
venn_tree <- transform_newick_string_to_venn_tree(newick_string)
mat_tree <- transform_venn_tree_to_matrix_tree(venn_tree);
atomized_newick_string <- strsplit(newick_string,"")[[1]];
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_string);
a_n_s <- length(atomized_newick_string);
rparens <- (1:a_n_s)[atomized_newick_string %in% "("];
lparens <- (1:a_n_s)[atomized_newick_string %in% ")"];
commas <- (1:a_n_s)[atomized_newick_string %in% ","];
numbers <- (1:a_n_s)[atomized_newick_string %in% 0:9];
comma_number <- commas[(commas+1) %in% numbers];
comma_number <- comma_number[(comma_number-1) %in% lparens];

#((((1,((2,(((((3,(((((11,(10,(9,(8,(4,(7,(5,6))))))),((12,13),(((16,(14,15)),(17,18)),((((21,(19,20)),(22,(23,(24,(27,(25,26)))))),28),29)))),30),31),32)),33),34),35),36)),37)),38),((((39,(40,41)),(42,43)),(((44,45),(46,47)),48)),(49,50))),51)
cn <- 1;
while (length(comma_number)>=1)	{
	clade_end_dists <- clade_boundaries$rp-comma_number[cn];
	this_clade <- clade_boundaries[match(min(clade_end_dists[clade_end_dists>0]),clade_end_dists),]
	otu_no <- atomized_newick_string[numbers[numbers>comma_number[cn]][1]:(lparens[lparens>comma_number[cn]][1]-1)];
	atomized_newick_string[comma_number[cn]:(comma_number[cn]+length(otu_no))]
	transpose2 <- atomized_newick_string[(this_clade$lp+1):(comma_number[cn]-1)]
	transpose1 <- c(otu_no,",");
	rewrite <- (this_clade$lp+1):(this_clade$rp-1);
	atomized_newick_string[rewrite] <- c(transpose1,transpose2);
	
#	local_lparents <- atomized_newick_string[(1+max(numbers[numbers<comma_number[cn]])):(comma_number[cn]-1)];
#	nsclades <- length(local_lparents); # number of sister clades;
#	prec_rparens <- length(rparens[rparens<comma_number[cn]]); # right boundaries of clades so far;
#	sis_clade_boundary <- rparens[(prec_rparens+1-nsclades):prec_rparens];
#	transpose_1 <- atomized_newick_string[sis_clade_boundary[1]:(comma_number[cn]-1)];
#	transpose_2 <- c(otu_no,",");
#	rewrite <- sis_clade_boundary[1]:(comma_number[cn]+length(otu_no));
#	atomized_newick_string[rewrite] <- c(transpose_2,transpose_1);
	rparens <- (1:a_n_s)[atomized_newick_string %in% "("];
	lparens <- (1:a_n_s)[atomized_newick_string %in% ")"];
	commas <- (1:a_n_s)[atomized_newick_string %in% ","];
	numbers <- (1:a_n_s)[atomized_newick_string %in% 0:9];
	comma_number <- commas[(commas+1) %in% numbers];
	comma_number <- comma_number[(comma_number-1) %in% lparens];
	ladderized_newick_string <- paste(atomized_newick_string,collapse="");
	}
return(ladderized_newick_string);
}

#### MODIFY TREES STORED IN MEMORY ####
transform_matrix_tree_to_vector_tree <- function (matrix_tree)	{
Nnode <- dim(matrix_tree)[1]
notus <- max(matrix_tree)-Nnode
ttus <- max(matrix_tree)
vector_tree <- vector(length=ttus)
vector_tree[notus+1] <- -1
for (n in Nnode:1)
	vector_tree[matrix_tree[n,]] <- n+notus
return(vector_tree)
}

transform_matrix_tree_to_venn_tree <- function(matrix_tree) {
nNodes <- nrow(matrix_tree);
notu <- max(matrix_tree)-nNodes;
venn_tree <- array(0,dim=c(nNodes,notu));
venn_tree[,1:ncol(matrix_tree)] <- matrix_tree;
#for (nd in nNodes:16)  {
nd <- nNodes;
while (nd > 0)  {
  if (max(venn_tree[nd,])>notu) {
    daughter_clades <- venn_tree[nd,venn_tree[nd,]>notu];
    daughter_nodes <- daughter_clades-notu;
    dclades <- length(daughter_clades);
    for (dc in dclades:1) {
      f1c <- match(daughter_clades[dc],venn_tree[nd,]);
      daughter_div <- sum(venn_tree[daughter_nodes[dc],]>0);
      replace_cells <- f1c:(f1c+daughter_div-1);
      grandkids <- venn_tree[daughter_nodes[dc],venn_tree[daughter_nodes[dc],]>0];
      displaced_cuz <- venn_tree[nd,(1+f1c):notu][venn_tree[nd,(1+f1c):notu]>0]
      if (length(displaced_cuz)>0)  {
        displaced_cells <- (f1c+daughter_div-1)+1:length(displaced_cuz);
        venn_tree[nd,displaced_cells] <- displaced_cuz;
        }
      venn_tree[nd,replace_cells] <- grandkids;
      #venn_tree[,1:6]
      #venn_tree[nd,]f1c:notu] <- c(venn_tree[daughter_nodes[dc],1:daughter_div],venn_tree[nd,(f1c+1):(notu+1-daughter_div)]);
      }
    }
  nd <- nd-1;
  }
return(venn_tree);
}

# routine to extract vector tree from matrix giving total progeny of a node
transform_venn_tree_to_vector_tree <- function (venn_tree)	{
n_Nodes <- nrow(venn_tree);
notus <- ncol(venn_tree);
max_otus <- max(venn_tree);
base <- max_otus+1
new_vector_tree <- vector(length=(max_otus+n_Nodes));
otus <- sort(as.integer(venn_tree[1,]));
for (s in 1:notus)	{
	spc <- otus[s];
	new_vector_tree[spc] <- max_otus+sort(which(venn_tree==spc,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
#	new_vector_tree[spc] <- max_otus+sort(which(venn_tree==spc,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
new_vector_tree[base] <- -1
new_vector_tree[base+1] <- base
for (n in 3:n_Nodes)	{
	htu <- max_otus+n
	lead <- venn_tree[n,1]
	new_vector_tree[htu] <- max_otus+sort(which(venn_tree[1:(n-1),]==lead,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
return(new_vector_tree)
}

transform_venn_tree_to_matrix_tree <- function(venn_tree)  {
return(transform_vector_tree_to_matrix_tree(transform_venn_tree_to_vector_tree(venn_tree)));
}

transform_vector_tree_to_matrix_tree <- function(vector_tree)	{
node_rosetta <- sort(unique(vector_tree[vector_tree>0]))
Nnodes <- length(node_rosetta)
maxtomy <- max((hist(vector_tree[vector_tree>1],breaks=((min(vector_tree[vector_tree>1])-1):max(vector_tree[vector_tree>1])),plot=FALSE)$counts))
#order(vector_tree)[2:length(vector_tree)]
node_rich <- vector(length=Nnodes)
matrix_tree <- matrix(0,Nnodes,maxtomy)
for (i in 1:length(vector_tree))	{
	node <- match(vector_tree[i],node_rosetta)
	if(!is.na(node))	{
		node_rich[node] <- node_rich[node]+1
		matrix_tree[node,node_rich[node]] <- i
		}
#	if (vector_tree[i]>=node_rosetta[1])	{
#		node <- match(vector_tree[i],node_rosetta)
#		node_rich[node] <- node_rich[node]+1
#		matrix_tree[node,node_rich[node]] <- i
#		}
	}
return(matrix_tree)
}

transform_vector_tree_to_venn_tree <- function(vector_tree)	{
ohtu <- length(vector_tree);
base <- match(-1,vector_tree);
otu <- base-1
htu <- ohtu-otu
venn_tree <- matrix(0,ohtu,otu)
for (i in 1:otu)	venn_tree[base,i] <- i

node_rich <- vector(length=ohtu)
for (sp in otu:1)	if (vector_tree[sp]!=0)			node_rich[vector_tree[sp]] <- node_rich[vector_tree[sp]]+1
for (nd in ohtu:(base+1))	if (vector_tree[nd]>0)	node_rich[vector_tree[nd]] <- node_rich[vector_tree[nd]]+node_rich[nd]
node_div <- vector(length=ohtu)
for (sp in 1:otu)	{
	node_div[vector_tree[sp]] <- node_div[vector_tree[sp]]+1
	venn_tree[vector_tree[sp],node_div[vector_tree[sp]]] <- sp
	}

for (nd in ohtu:(base+1))	{
	anc <- vector_tree[nd]
	for (i in 1:node_div[nd])	{
		node_div[anc] <- node_div[anc]+1
		venn_tree[anc,node_div[anc]] <- venn_tree[nd,i]
		}
	}
#venn_tree[base:ohtu,1:15]

return(venn_tree[base:ohtu,])
}

create_phylo_class_from_nexus_tree <- function(vector_tree,tip.label)	{
htu1 <- min(vector_tree[vector_tree>0])
Nnode <- 1+(max(vector_tree)-htu1)
otus <- htu1-1
edges <- matrix(0,otus+Nnode-1,2)
j <- 0
for (i in 1:length(vector_tree))	{
	if (vector_tree[i]!=-1)	{
		j <- j+1
		edges[j,1] <- vector_tree[i]
		edges[j,2] <- i
		}
	}
output_tree <- list(edges,tip.label,Nnode)
names(output_tree) <- c("edge","tip.label","Nnode")
class(output_tree) <- "phylo"
#str(output_tree)
return(output_tree)
}

create_phylo_class_from_nexus_tree_file <- function(nexustreefile,tip.label)	{
vector_tree <- read_newick_tree_from_file(nexustreefile)
htu1 <- min(vector_tree[vector_tree>0])
Nnode <- 1+(max(vector_tree)-htu1)
otus <- htu1-1
edges <- matrix(0,otus+Nnode-1,2)
j <- 0
for (i in 1:length(vector_tree))	{
	if (vector_tree[i]!=-1)	{
		j <- j+1
		edges[j,1] <- vector_tree[i]
		edges[j,2] <- i
		}
	}
output_tree <- list(edges,tip.label,Nnode)
names(output_tree) <- c("edge","tip.label","Nnode")
class(output_tree) <- "phylo"
#str(output_tree)
return(output_tree)
}

transmogrify_additive_dependents_to_multistate_old <- function(ind_char,dep_chars,chmatrix,secondary_dependencies,INAP=-22,UNKNOWN=-11,theoretical=TRUE)	{
# theortical: if TRUE, then Q matrix & recoding allows for all theoretically possible combinations
# unknown_inap_sep: code unknowns separately; this is a bad idea that I regret and that we'll cut
notu <- nrow(chmatrix);
combos <- chmatrix[,c(ind_char,dep_chars)];
rstates <- count_states(combos,UNKNOWN,INAP)
combos <- combos[!(rowMaxs(combos)==UNKNOWN & rowMins(combos)==UNKNOWN),];
combos <- combos[combos[,2]!=INAP,];
# get rid of missing examples
combos <- combos[!(1:nrow(combos)) %in% unique(which(combos==UNKNOWN,arr.ind = TRUE)[,1]),];
# get secondary dependents
secondaries <- unique(which(combos==INAP,arr.ind = TRUE)[,2]);
semi_indies <- secondary_dependencies[secondary_dependencies!=ind_char];

nchars <- ncol(combos);
dchars <- length(dep_chars);
missing_combos <- unique(combos)[sort(unique(which(unique(combos)==UNKNOWN,arr.ind = TRUE)[,1])),];
nstates <- count_states(combos);
nstates[secondaries] <- nstates[secondaries]+1;
#combos <- combos[order(combos[,1],combos[,2],combos[,3]),];
for (cn in ncol(combos):1)	combos <- combos[order(combos[,cn]),];
all_obs_complex_combos <- unique(combos);
all_obs_complex_combos <- all_obs_complex_combos[all_obs_complex_combos[,1]>=0,];
ind_states <- sort(unique(chmatrix[,ind_char][chmatrix[,ind_char]>=0]))
key_states <- unique(all_obs_complex_combos[,1]);
wrong_states <- ind_states[!ind_states %in% key_states];

if (!theoretical)	{
	all_poss_combos <- unique(combos);
	for (cu in ncol(all_poss_combos):1)	all_poss_combos <- all_poss_combos[order(abs(all_poss_combos[,cu])),];
	all_poss_combos <- all_poss_combos[(1:nrow(all_poss_combos))[!(1:nrow(all_poss_combos)) %in% which(all_poss_combos==UNKNOWN,arr.ind=TRUE)[,1]],];
	} else	{
#	unq_combos <- array(0,dim=c(prod(nstates),length(nstates)));
#	unq_combos[,1] <- unique(unique(combos)[,1]);
	# make sure that secondary dependents get inapplicable in the combos
	# should the Q-matrix for the "wrong" independent character reflect this?
	#which(unique(combos)==INAP,arr.ind=TRUE);
#	for (ds in 1:dchars)	{
	dc <- nchars+1;
	while (dc >1)	{
		dc <- dc-1;
		char_states <- sort(unique(combos[,dc]));
		char_states <- char_states[!char_states %in% c(UNKNOWN,INAP)];
		cs <- 1;
		while (char_states[cs]<0)	{
			char_states <- c(char_states,unravel_polymorph_badass(char_states[cs]))
			char_states <- unique(sort(char_states[2:length(char_states)]));
			}
#		if (dc %in% secondaries)	unq_combos[unq_combos[,dc]==(nstates[dc]-1),dc] <- INAP;
		if (dc %in% secondaries)	char_states <- c(char_states,INAP);
		nstates_ch <- length(char_states);
		if (dc==nchars)	{
			all_poss_combos <- array(char_states,dim=c(nstates_ch,1));
			} else	{
			cs <- 1;
			existing_combos <- nrow(all_poss_combos)
			added_char_states <- rep(char_states[cs],existing_combos);
			unq_combos_orig <- all_poss_combos;
			while (cs<nstates_ch)	{
				cs <- cs+1;
				all_poss_combos <- rbind(all_poss_combos,unq_combos_orig);
				added_char_states <- c(added_char_states,rep(char_states[cs],existing_combos));
				}
			all_poss_combos <- cbind(added_char_states,all_poss_combos);
			}
		}
#	while (dc > 1)	{
#		dc <- dc-1;
#		if (dc<nchars)	{
#			this_char_states <- c();
#			for (dss in 1:nstates[dc])	#{
#				this_char_states <- c(this_char_states,rep(dss-1,prod(nstates[(dc+1):nchars])));
#				if (dss==nstates[dc] && (dc %in% secondaries))	{
#					this_char_states <- c(this_char_states,rep(INAP,prod(nstates[(dc+1):dchars])));
#					this_char_states <- c(this_char_states,rep(INAP,nstates[dc-1]));
#					} else	{
#					this_char_states <- c(this_char_states,rep(dss-1,nstates[dc-1]));
#					}
#				}
#			} else	{
#			this_char_states <- (1:nstates[dc])-1;
#			}
#		ttl_code <- round(nrow(unq_combos)/length(this_char_states),0);
#		unq_combos[,dc] <- rep(this_char_states,ttl_code);
#		}
	}
all_poss_combos <- all_poss_combos[all_poss_combos[,1] %in% key_states,];
# alternative for secondary independents: replace all pairs with -22 & then reduce via unique(all_poss_combos);

sc <- 0;
while (sc < length(secondaries))	{
	sc <- sc+1;
	independents <- 2:ncol(all_poss_combos);
	semi_indy <- find_independent_character(dchar=secondaries[sc],independents=2:ncol(all_poss_combos),chmatrix=combos,UNKNOWN,INAP);
	this_pair <- unique(combos[,c(semi_indy,secondaries[sc])]);
	semi_key_state <- unique(this_pair[this_pair[,2]!=INAP,][,1]);
	semi_wrn_state <- unique(this_pair[,1][!this_pair[,1] %in% semi_key_state]);
	
	# get rid of keystate:inapplicable pairs
	zz <- (1:nrow(all_poss_combos));
	xx <- zz[all_poss_combos[,semi_indy] %in% semi_key_state]
	yy <- zz[all_poss_combos[,secondaries[sc]]==INAP];
	xx <- xx[xx %in% yy];
	all_poss_combos <- all_poss_combos[!zz %in% xx,];
	
	# get rid of nonkeystate:state-pairs
	zz <- (1:nrow(all_poss_combos));
	xx <- zz[all_poss_combos[,semi_indy] %in% semi_wrn_state]
	yy <- zz[all_poss_combos[,secondaries[sc]]!=INAP];
	xx <- xx[xx %in% yy];
	all_poss_combos <- all_poss_combos[!zz %in% xx,];
	}

wrstates <- length(wrong_states);
null_combos <- cbind(wrong_states,array(INAP,dim=c(wrstates,ncol(combos)-1)));
colnames(null_combos) <- colnames(all_poss_combos)
all_poss_combos <- rbind(null_combos,all_poss_combos);

ucombos <- rcombos <- ttl_states <- nrow(all_poss_combos);
#if (unknown_inap_sep)	{
#	unq_combos <- rbind(unq_combos,missing_combos);
#	ttl_states <- nrow(unq_combos);
#	rcombos <- ttl_states-(wrstates+length(missing_combos));
#	ucombos <- ttl_states-length(missing_combos);
#	} else	{
rcombos <- ttl_states-wrstates;
#	}

state_combos <- rep("",ttl_states);
for (uc in 1:ttl_states)	{
	this_combo <- paste(all_poss_combos[uc,],collapse="");
	this_combo <- gsub(as.character(INAP),"-",this_combo);
	this_combo <- gsub(as.character(UNKNOWN),"?",this_combo);
	state_combos[uc] <- this_combo;
	}
state_combos <- gsub(as.character(INAP),"",state_combos);
rownames(all_poss_combos) <- state_combos;
colnames(all_poss_combos) <- c(ind_char,dep_chars);
# make sure that secondaries dependencies are weeded out.
#	if we see only •2-, then make sure that •20 & •21 are eliminated
# prepare Q-Matrix
if (nrow(all_poss_combos)>100)	{
	print("Getting basic distances among state combinations")
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=TRUE);
	} else	{
	Q <- pairwise_differences_discrete(all_poss_combos,UNKNOWN=UNKNOWN,INAP=INAP,progress_bar=F);
	}
# if there are secondaries, the figure out how to weight them here!!!
colnames(Q) <- rownames(Q) <- state_combos;


Q[wrstates,(wrstates+1):ttl_states] <- 1/rcombos;
si <- 0;
while (si < length(semi_indies))	{
	si <- si+1;
	sc <- match(semi_indies[si],c(ind_char,dep_chars));
	cs <- (1:length(c(ind_char,dep_chars)))[secondary_dependencies==semi_indies[si]];
	relv_combos <- unique(all_poss_combos[,c(sc,cs)]);
	relv_combos <- relv_combos[(1:nrow(relv_combos))[!(1:nrow(relv_combos)) %in% unique(which(relv_combos==INAP,arr.ind=TRUE)[,1])],];
	relv_combos_all <- all_poss_combos[,c(sc,cs)];
	if (length(cs)>1)	{
		key_combos <- which(all_poss_combos[,cs]==INAP,arr.ind=TRUE)[,1]
		} else	{
		key_combos <- as.numeric(which(all_poss_combos[,cs]==INAP,arr.ind=TRUE))
		}
	# fix revl_combos to get all of the right matches & not just those for one possibility
	key_combos <- key_combos[!key_combos %in% key_combos[all_poss_combos[,sc]==INAP]];
	relv_combos_2 <- all_poss_combos[,c(sc,cs)];
	relv_combos_2 <- relv_combos_2[!(1:ttl_states) %in% which(relv_combos_2==INAP,arr.ind=TRUE)[,1],];
	combos_key <- match(rownames(relv_combos_2),rownames(Q));
#	Q[key_combos,combos_key]==1;
	for (i in 1:length(key_combos))
		Q[key_combos[i],combos_key][Q[key_combos[i],combos_key]==1] <- 1/nrow(relv_combos);
#	1/nrow(relv_combos)
	}
Q[Q>1] <- 0;
for (qq in 1:ttl_states)	{
	Q[qq,] <- Q[qq,]/sum(Q[qq,1:ucombos]);
	Q[qq,qq] <- -1;
#	Q[qq,qq] <- -sum(Q[qq,1:ucombos])
	}
#write.csv(Q,"Q.csv",row.names = TRUE)
new_multistate <- vector(length=notu);
for (nn in 1:notu)
	new_multistate[nn] <- match_vector_to_matrix_row(test_vector=chmatrix[nn,c(ind_char,dep_chars)],test_matrix=all_poss_combos)-1;
all_states <- c(0:9,letter_states,more_letter_states);
if (nrow(all_poss_combos)>10)
	if (is.numeric(new_multistate))
		new_multistate[!is.na(new_multistate)] <- all_states[1+new_multistate[!is.na(new_multistate)]]
#cbind(new_multistate,chmatrix[,c(ind_char,dep_chars)])
#if (unknown_inap_sep)	{
#	new_multistate[is.na(new_multistate)] <- UNKNOWN;
#	} else	{
prob_child <- (1:notu)[is.na(new_multistate)];
combos <- chmatrix[,c(ind_char,dep_chars)];
pc <- 0;
polymorphs <- unique(combos[combos<0])
polymorphs <- polymorphs[!polymorphs %in% c(UNKNOWN,INAP)];
while (pc < length(prob_child))	{
	pc <- pc+1;
	pcc <- prob_child[pc];
	if (combos[pcc,1]==UNKNOWN)	{
		new_multistate[pcc] <- UNKNOWN;
		} else	{
			# make it polymorphic for all possible states
		this_combo <- combos[pcc,(0:dchars)+1];
#		doofi <- (1:length(this_combo))[this_combo==UNKNOWN];
		# this if/else probably is unneeded now!
		if (sum(this_combo %in% polymorphs)==0)	{
			set_chars <- (1:length(this_combo))[this_combo!=UNKNOWN];
#			set_chars <- set_chars[set_chars>1];
			set_states <- this_combo[set_chars];
			poss_combos <- all_poss_combos[all_poss_combos[,1] %in% key_states,];
			ss <- 0;
			while (ss < length(set_chars))	{
				ss <- ss+1;
				poss_combos <- subset(poss_combos,poss_combos[,set_chars[ss]]==set_states[ss]);
				}
			polymorph <- 0;
			if (nrow(all_poss_combos)>10)	polymorph <- c();
			for (cp in 1:nrow(poss_combos))	{
				this_state <- (row.match(poss_combos[cp,],as.data.frame(all_poss_combos))-1);
				if (nrow(all_poss_combos)>10)	{
					polymorph <- c(polymorph,all_states[this_state]);
#					base64encode(this_state);
#					for (i in 1:5)	print(base64encode(10^(i-1)));
					} else	{
					polymorph <- polymorph-((10^(cp-1))*this_state);
					}
				}
			} else	{
			set_chars <- (1:nchars)[!this_combo %in% c(UNKNOWN,polymorphs)];				unset_chars <- (1:nchars)[this_combo %in% c(UNKNOWN,polymorphs)];
			un <- 0;
			missing <- unset_chars[this_combo[unset_chars]==UNKNOWN];
			polys <- unset_chars[!this_combo[unset_chars] %in% UNKNOWN];
			while (un < length(missing))	{
				un <- un+1;
				mchar <- missing[un];
				rstates <- (1:nstates[mchar])-1;
				this_combo[mchar] <- 0;
				for (rs in 1:length(rstates))
					this_combo[mchar] <- this_combo[rs]-rstates[rs]*(10^(rs-1));
				}
#				set_chars <- set_chars[set_chars>1];
			set_states <- this_combo[set_chars];
			poss_combos <- all_poss_combos;
#			print(nrow(poss_combos));
			# reduce the possible combinations to those consistent
			sc <- length(set_chars);
#			for (sc in 1:length(set_chars))	{
			while (sc > 1)	{
				poss_combos <- subset(poss_combos,poss_combos[,set_chars[sc]]==set_states[sc]);
#				print(nrow(poss_combos));
				sc <- sc-1;
				}
			uc <- 0;
			while (uc < length(polys))	{
				uc <- uc+1;
				poss_combos <- poss_combos[poss_combos[,polys[uc]] %in% unravel_polymorph_badass(this_combo[polys[uc]]),];
				}
			pstates <- match(rownames(poss_combos),rownames(all_poss_combos));
			polymorph <- 0;
			if (nrow(all_poss_combos)>10)	polymorph <- c();
			for (cp in 1:nrow(poss_combos))	{
				this_state <- (row.match(poss_combos[cp,],as.data.frame(all_poss_combos))-1);
				if (nrow(all_poss_combos)>10)	{
					polymorph <- c(polymorph,all_states[this_state]);
#					base64encode(this_state);
#					for (i in 1:5)	print(base64encode(10^(i-1)));
					} else	{
					polymorph <- polymorph-((10^(cp-1))*this_state);
					}
				}
			}
#		if (nrow(unq_combos)>10)	{
#			if (is.numeric(new_multistate))
#				new_multistate[!is.na(new_multistate)] <- all_states[1+new_multistate[!is.na(new_multistate)]]
			new_multistate[pcc] <- paste("(",paste(polymorph,collapse=""),")",sep="");
#			} else	{
#			new_multistate[pcc] <- polymorph;
#			}
		}
	}
#	}
#new_multistate[new_multistate==UNKNOWN] <- "?";

#cbind(new_multistate,chmatrix[,c(ind_char,dep_chars)])
output <- list(all_poss_combos,Q,new_multistate);
names(output) <- c("unique_combinations","Q","new_character");
return(output);
}
