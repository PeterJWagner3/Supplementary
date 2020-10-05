newick_verbotten <- c(".","?","\"","\'");
letter_states <- LETTERS[!LETTERS %in% c("I","O")];
zzzz <- 0.25;
MAXNO <- 1.797693e+308

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Routines for handling character matrices ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
mundify_nexus_text <- function(nexus_line)	{
nexus_line <- gsub("\xd4","",nexus_line);
nexus_line <- gsub("\xd5","",nexus_line);
nexus_line <- gsub("\x87","a",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x93\xa0\xbc","ae",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x83\xa0\xbc","c",nexus_line);
nexus_line <- gsub("\x8e","e",nexus_line);
nexus_line <- gsub("\x8f","e",nexus_line);
nexus_line <- gsub("\x92","i",nexus_line);
nexus_line <- gsub("\xbf","o",nexus_line);
nexus_line <- gsub("\x9a","o",nexus_line);
nexus_line <- gsub("\x97","o",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8e\x93\xa4\xbc","o",nexus_line);
nexus_line <- gsub("\x9f","ue",nexus_line);
nexus_line <- gsub("\xd0","-",nexus_line);
nexus_line <- gsub("\xd2","\"",nexus_line);
nexus_line <- gsub("\xd3","\"",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x90\xbc","\'",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x8c\xbc","ƒ",nexus_line);
nexus_line <- gsub("\xa7","ß",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8c\xa6\x88\xbc","≤",nexus_line);
nexus_line <- gsub("\xb3","≥",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x96\x8c\xbc","≈",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\xa6\x98\xbc","˚",nexus_line);
nexus_line <- gsub("\xb6","∂",nexus_line);
nexus_line <- gsub("\xc6","∆",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\xb6\x88\xbc","∑",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x86\x88\xbc","Ω",nexus_line);
nexus_line <- gsub("\xa5"," ",nexus_line);
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
nexus_line <- gsub("\x88","a",nexus_line);
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

accersi_data_from_chosen_nexus_file <- function(polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partition="", trend_partition="")	{
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
ml <- 0;
#i <- 1

for (i in 1:length(nexus))  {
	nexus[i] <- mundify_nexus_text(nexus_line = nexus[i])
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	if (length(j)>ml) ml <- length(j)
	}
ml <- ml+1;	# LENGTH OF LONGEST LINE
	
# file is now a vector of characters.  Turn it into a matrix with one char per cell
nexusfile <- matrix("\n",length(nexus),ml)
for (i in 1:length(nexus))  {
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	for (k in 1:length(j))      nexusfile[i,k] <- j[k]
	if ((length(j)+2)<ml)
		for (k in (length(j)+2):ml) nexusfile[i,k] <- ""
	}

top <- 0;
ln <- 1;		# this is the row with the word "Matrix": character data starts next.
while (top==0)	{
	em_nexus <- gsub("\t","",nexus[ln]);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	if (!is.na(match("matrix",tolower(nexus_words))))	{
		top <- ln;
		}
	else	ln <- ln+1;
	}
top <- top+1;	# this will give the first row of data
# skip the comment text denoting character numbers (if present)
while(nexusfile[top,1]=="[" || nexusfile[top,1]==" ") top <- top+1

missing <- "?";
gap <- "-";
notu <- nchars <- strat <- range <- geog <- 0
for (i in 2:top)  {
	while ((nexusfile[i,1]=="[" || nexusfile[i,1]=="\n") && i<top)	i <- i+1;
	em_nexus <- gsub("\t","",nexus[i]);
	em_nexus <- gsub("="," = ",em_nexus);
	em_nexus <- gsub(";"," ; ",em_nexus);
	em_nexus <- gsub(",","",em_nexus);
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
		j <- 1+match("gap",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		gap <- nexus_words[j];
		}
	if (!is.na(match("missing",tolower(nexus_words))))	{
		j <- 1+match("missing",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		missing <- nexus_words[j];
		}
	if (!is.na(match("fa",tolower(nexus_words))))	{
		strat <- as.numeric(nexus_words[match("fa",tolower(nexus_words))-1]);
		}
	if (!is.na(match("la",tolower(nexus_words))))	{
		range <- as.numeric(nexus_words[match("la",tolower(nexus_words))-1]);
		}
	if (!is.na(match("geog",tolower(nexus_words))) || !is.na(match("geography",tolower(nexus_words))))	{
		geog <- c(nexus_words[match("geog",tolower(nexus_words))-1],nexus_words[match("geography",tolower(nexus_words))-1]);
		geog <- as.numeric(geog[!is.na(geog)]);
		}
	}

extra <- 0;
if (strat>0)	{
	if (range>0)	{
		nchars <- nchars-2
		extra <- 2
		} else {
		nchars <- nchars-1
		extra <- 1
		} 
	strat_ranges <- data.frame(FA=as.numeric(rep(0,notu)),LA=as.numeric(rep(0,notu)));
	}
if (geog>0)	{
	nchars <- nchars-1
	geography <- vector(length=notu)
	extra <- extra+1
	}
	
taxa <- vector(length=notu);
nstates <- array(0,dim=nchars);
chmatrix <- matrix(0,notu,nchars);
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
		if (length(xxx)==0 || is.na(xxx))
			xxx <- "";
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

# look for rate variation partitions
if (rate_partition!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(rate_partition,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			nexus[ln] <- gsub(":"," : ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  :"," :",nexus[ln]);
			nexus[ln] <- gsub(":  ",": ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(rate_partition,breakup_this_line):length(breakup_this_line)];
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
			got_splits<- T;
			}
		}
	} else	character_rate_partitions <- rep("imagine",nchars);

if (trend_partition!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(trend_partition,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(trend_partition,breakup_this_line):length(breakup_this_line)];
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
			got_splits<- T;
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
		if (!is.na(match("ord:",tolower(xxx))) && !is.na(match("revbayes",tolower(xxx))))	{
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
	if (nexusfile[s,1]=="'" || nexusfile[s,2]=="'")	{
		jj <- ((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"]);
		i <- max((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"])
		taxa[tx] <- pracma::strcat(nexusfile[s,(jj[1]+1):(jj[2]-1)])
		i <- i+1
		while (nexusfile[s,i]==" " && i<ncol(nexusfile))	i <- i+1
		}	else	{
		i <- 1
		if (nexusfile[s,1]!="\"")  {
			while (nexusfile[s,i]=="\t")	i <- i+1
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
			taxa[tx] <- nexusfile[s,2]
			i <- 3
			while (nexusfile[s,i]!=" " && nexusfile[s,i+1]!=" " && i<ncol(nexusfile))	{
				taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
				i <- i+1
				}
			}
		# now, get to characters
		while ((nexusfile[s,i]==" " || nexusfile[s,i]=="\t") && i<ncol(nexusfile))
			i <- i+1
		}
	k <- i;
	endline <- match("\n",nexusfile[s,])
	if (is.na(endline))	endline <- length(nexusfile[s,])
	if ((endline-k)==(nchars+extra))	{
		# true if there are no polymorphic characters for the taxon
		dummy <- nexusfile[s,k:(endline-1)]
		dummy[dummy==missing] <- UNKNOWN
		dummy[dummy==gap] <- INAP
		letterstate <- dummy
		dummy <- sapply(letterstate,switch_letter_state_to_numeric)
		chmatrix[tx,] <- as.numeric(dummy[1:nchars]);
		if (strat>0)	{
			strat_ranges$FA[tx] <- strat_ranges$LA[tx] <- as.numeric(dummy[strat])
			if (range>0)	strat_ranges$LA[tx] <- as.numeric(dummy[range])
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
			if (c<=nchars)	{
				if (nexusfile[s,i]=="(" || nexusfile[s,i]=="{")	{
					if (polymorphs==TRUE || polymorphs==1)	{
						i <- i+1
						w <- as.numeric(nexusfile[s,i])
						chmatrix[tx,c] <- -1*as.numeric(nexusfile[s,i])
						if ((1+w)>nstates[c])  nstates[c] <- 1+w
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
					chmatrix[tx,c] <- switch_letter_state_to_numeric(nexusfile[s,i]);
					}	else if (nexusfile[s,i]>="0" && nexusfile[s,i]<="9") {
					chmatrix[tx,c] <- as.numeric(nexusfile[s,i]);
					}
				if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
				i <- i+1
				}  else {
				if (c==strat)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,1]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,1]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					if (range==0)	strat_ranges[tx,2] <- strat_ranges[tx,1]
					i <- i+1
					} else if (c==range)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,2]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,2]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					i <- i+1
					} else if (c==geog)	{
						if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
							geography[tx]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
							geography[tx]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					}
				} # end non-morphological data
#			print(nexusfile[s,k:83]);
#			print(chmatrix[tx,])
			if (nexusfile[s,i+1]=="\n" || i==(mxln-1)) c <- nchars+extra;
			}
		}
#	chmatrix[tx,];
#	tx <- tx+1;
#	s <- s+1
	}

if (trend_partition!="")	print("Choose a file giving the order in which taxa appear (with low numbers = early): ");
chmatrix <- mundify_character_matrix(chmatrix,minst=0,UNKNOWN,INAP);	# clean up coding
if (trend_partition!="")	{
	appearance_order_file <- file.choose();
	appearance_order <- read.csv(appearance_order_file,header=T);
	otu_fas <- match(appearance_order$appearance,sort(unique(appearance_order$appearance)));
	for (nch in 1:ncol(chmatrix))	{
		otu_states <- chmatrix[,nch];
		chmatrix[,nch] <- rescore_states_by_first_appearances(otu_states,otu_fas,UNKNOWN,INAP);
		}
	}
nstates <- count_states(chmatrix,UNKNOWN,INAP);

tree_found <- 0;
while (s<length(nexus) && tree_found==0)	{
	while (nexus[s]!= "BEGIN TREES; " && s<length(nexus))
		s <- s+1;
	if (s<length(nexus))	{
		while (tree_found==0 && s<length(nexus))	{
			s <- s+1
			jj <- strsplit(nexus[s],split=c("\t"," "),fixed=TRUE)[[1]];
			jj <- paste(jj,collapse="")
			jj <- strsplit(jj,split=" ",fixed=TRUE)[[1]];
			if (sum(jj=="TREE")>0 || sum(jj=="tree")>0)	tree_found <- 1;
			}
		newick_string <- jj[length(jj)];
		tree <- read_newick_string(newick_string);
		tree_found <- 1
		s <- length(nexus);
		}
	}

row.names(chmatrix) <- taxa

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
	character_trend_partitions <- character_trend_partitions[keepers];
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

mundify_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix)	# replaces spc to standardize coding.
nchars <- ncol(chmatrix)
for (ch in 1:nchars)	{
#ch <- 0
#while (ch<=nchars)	{
#	ch <- ch+1;
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP])
	if (length(rem)>0)	{
		test <- chmatrix[-rem,ch]
		}	else test <- chmatrix[,ch]
	# check polymorphics for anything that needs to be changed
	if (length(rem) < notu)	{
		polys <- sum(test<0)
		coded <- sort(unique(test[test>=0]))
		if (polys>0)	{
			examps <- test[test<0]
			polycoded <- sort(unique(test[test<0]))
			for (i in 1:length(examps))	{
				polystates <- unravel_polymorph(examps[i])
				coded <- sort(unique(c(coded,polystates)))
#				if (min(polystates)<minstch)	minstch <- min(polystates)
				}
			} else	{
			polycoded <- c();
			}
		minstch <- min(coded)
	# eliminate gaps in states
		if ((max(coded)-(minstch-1)) > length(coded))	{
			gaps_in_codes <- c(1,coded[2:length(coded)]-coded[1:(length(coded)-1)])-1
			new_codes <- coded
			for (st in 2:length(coded))	{
				new_codes[st:length(coded)] <- coded[st:length(coded)]-gaps_in_codes[st]
				}
		
			for (st in 2:length(coded))	{
				# start at 2 because the first state cannot have a gap before it!
#				exp_st <- (minstch+(st-1))
#				if (coded[st]!=exp_st)	{
				if (coded[st]!=new_codes[st])	{
					rec <- (1:notu)[chmatrix[,ch]==coded[st]]
					chmatrix[rec,ch] <- new_codes[st]
					redo_poly <- 1
					while (redo_poly <= length(polycoded))	{
						polystates <- unravel_polymorph(polycoded[redo_poly])
						polystates[polystates==coded[st]] <- new_codes[st]
						newpolystates <- ravel_polymorph(polystates)
						if (newpolystates != polycoded[redo_poly])	{
							chmatrix[(1:notu)[chmatrix[,ch]==polycoded[redo_poly]],ch] <- newpolystates
							}
#						polycoded[redo_poly] <- chmatrix[,ch][chmatrix[,ch] %in% polycoded[redo_poly]] <- newpolystates
						redo_poly <- redo_poly+1
						}
					coded[st] <- new_codes[st]
					}
				}
			}
		# standardize minimum state
		if (minstch!=minst)	{
			adj <- minst-minstch
			test2 <- (1:notu)[chmatrix[,ch]>=0]
			chmatrix[test2,ch] <- chmatrix[test2,ch]+adj
			if (polys>0)	{
				examps2 <- examps
				for (i in 1:polys)	{
					doh <- unravel_polymorph(examps[i])
					doh <- doh+adj
					examps2[i] <- ravel_polymorph(doh)
					}
				chmatrix[(chmatrix[,ch] %in% examps),ch] <- examps2
				}
			}
		}
	}
return(chmatrix)
}

switch_letter_state_to_numeric <- function(state)  {
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

count_states <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
nstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- unique(c(char_states,unravel_polymorph(char_states[1])))[2:length(char_states)];
			}
		}
	nstates <- c(nstates,length(char_states));
	}
return(nstates);
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Routines to read tree information from Newick files & format ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
read_newick_tree_from_file <- function(newicktree_file) {
newick_tree <- scan(file=newicktree_file,what=character(),sep="\n")
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
		sp <- as.numeric(nexus_string[i])+(sp*(10^p))
		p <- p+1
		if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Character reconstruction routines (likelihood, parsimony, etc.) ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
accersi_minimum_steps_history <- function(mtree,cdata,types,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# mtree: matrix tree, where each row gives the branches stemming from a node
# cdata: character data
# types: 1 for unordered, 0 for ordered
# UNKNOWN: numeric representation for "?"
# INAP: numeric representation for inapplicable or gap
# outgroup: the number of the outgroup taxon
otus <- nrow(cdata)
Nnode <- nrow(mtree)
ncharss <- ncol(cdata)
steps <- vector(length=ncharss)
for (c in 1:ncharss)	{
	cvector <- cdata[,c];
	type <- types[c];
	char_evolution <- Sankoff_character(mtree,cvector,type,UNKNOWN,INAP,outgroup);
	steps[c] <- char_evolution$Steps
	if (c==1)	{
		full_matrix <- char_evolution$States
		changes_matrix <- char_evolution$Derivation
		}	else	{
		full_matrix <- cbind(full_matrix,char_evolution$States)
		changes_matrix <- cbind(changes_matrix,char_evolution$Derivation)
		}
	}
char_labels <- vector(length=ncharss)
for (c in 1:ncharss)	{
	if (c<10)	{
		if (ncharss<10)	{
			char_labels[c] <- paste("ch_",c,sep="")
			}	else if (ncharss>9 && ncharss<100)	{
			char_labels[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_labels[c] <- paste("ch_00",c,sep="")
			}
		}	else if (c<100)	{
			if (ncharss<100)	{
			char_labels[c] <- paste("ch_",c,sep="")
			}	else	{
			char_labels[c] <- paste("ch_0",c,sep="")
			}
		}
	}

ttl_br <- otus+Nnode
nonzero <- 1
rbr <- 1
for (br in 1:ttl_br)	{
	nonzero <- sum(changes_matrix[br,]>0)
	if (nonzero>0)	{
		dch <- (1:ncharss)[changes_matrix[br,]>0]
		dst <- changes_matrix[br,changes_matrix[br,]>0]
		dbr <- rep(br,length(dch))
		if (rbr==1)	{
			branch_changes <- cbind(dbr,dch,dst)
			} else	{
			branch_changes <- rbind(branch_changes,cbind(dbr,dch,dst))
			}
		rbr <- rbr+1
		}
	}
rownames(branch_changes) <- rep("",dim(branch_changes)[1])
colnames(full_matrix) <- colnames(changes_matrix) <- char_labels
output <- list(steps,full_matrix[((otus+1):(otus+Nnode)),],branch_changes)
names(output) <- c("Steps","Ancestral_Reconstructions","Changes")
return(output)
}

Sankoff_character <- function(mtree,cvector,type=1,UNKNOWN=-11,INAP=-22,outgroup=1)	{
# method for reconstructing ancestral conditions.
obs_states <- sort(unique(cvector[cvector>=0]))
ttl_states <- length(obs_states)
Nnode <- nrow(mtree)
otus <- length(cvector)
scored <- (1:otus)
sankoff_matrix <- matrix(1,otus+Nnode,(length=ttl_states))
for (s in 1:otus)	{
	if (cvector[s] >= 0)	{
		st <- match(cvector[s],obs_states)
		sankoff_matrix[s,st] <- 0
		}	else	sankoff_matrix[s,] <- 0
	}
cvector <- c(cvector,rep(UNKNOWN,Nnode))
node_rich <- vector(length=Nnode);
for (n in Nnode:1)	{
	missing <- gap <- 0
	node_rich[n] <- f1 <- length(mtree[n,mtree[n,]>0])
	# if all taxa are scored
	# mtree[n,]
#	if (sum(mtree[n,(1:f1)] %in% scored)==f1)	{
	ht <- n+otus
	sankoff_matrix[ht,] <- 0
	for (s in 1:f1)	{
		#### add something to deal with inapplicables here.
		# list states that demand more than minimal change
		sp <- mtree[n,s]
		if (cvector[sp]!=INAP && cvector[sp]!=UNKNOWN)	{
			#	we will add a step to each of these because if sankoff is:
			#		0 1 1 for states 0, 1 & 2
			#	then we need 1 step from either 1 or 2
			n_s <- obs_states[sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
			a_s <- obs_states[!sankoff_matrix[sp,] %in% min(sankoff_matrix[sp,])]
			no_step <- match(n_s,obs_states)
			add_step <- match(a_s,obs_states)
			sankoff_matrix[ht,no_step] <- sankoff_matrix[ht,no_step]+min(sankoff_matrix[sp,])
			if (length(no_step) < length(obs_states))	{
				sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+min(sankoff_matrix[sp,])+1
	#			sankoff_matrix[ht,add_step] <- sankoff_matrix[ht,add_step]+sankoff_matrix[sp,add_step]
	#			sankoff_matrix[ht,]
				}
			}	else if (cvector[sp]==INAP)	{
			gap <- gap+1
			}	else if (cvector[sp]==UNKNOWN)	{
			missing <- missing+1
			}
		}
	if (missing==f1)	{
		cvector[ht] <- UNKNOWN
		scored <- c(scored,ht)
		}	else if (gap==f1)	{
		cvector[ht] <- INAP
		scored <- c(scored,ht)
		}	else if (sum(sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,]))==1)	{
		# see if we can fix a state
		ast <- match(min(sankoff_matrix[ht,]),sankoff_matrix[ht,])
		cvector[ht] <- obs_states[ast]
		scored <- c(scored,ht)
		}	else	{
		cvector[ht] <- ravel_polymorph(obs_states[sankoff_matrix[ht,] %in% min(sankoff_matrix[ht,])])
		}
	}

base <- otus+1
ch_steps <- min(sankoff_matrix[base,])

if (cvector[base]<0 && (cvector[base]!=UNKNOWN && cvector[base]!=INAP))	{
	#poss_starts <- (1:ttl_states)[sankoff_matrix[base,] %in% min(sankoff_matrix[base,])]
	poss_starts <- unravel_polymorph(cvector[base])
	# IF the designated outgroup is attached to the first node AND if it 
	#	has one of the possible nodal states, then assign that to basal node
	if (!is.na(match(outgroup,mtree[1,])) && !is.na(match(cvector[outgroup],poss_starts)))	{
		cvector[base] <- cvector[outgroup]
		}	else	{
		# otherwise, just grab one of them at random: it really doesn't matter
		grab <- ceiling(runif(1)/(1/length(poss_starts)))
		cvector[base] <- poss_starts[grab]
		}
	}

### start here: work up the tree, using cvector[htu] to set the state
changes_above <- vector(length=Nnode)
for (n in 1:Nnode)	{
	ht <- n+otus
	if (cvector[ht]!=INAP && cvector[ht]!=UNKNOWN)	{
		f1 <- node_rich[n]
		# if all taxa are scored in mtree[n,]
		if (sum(mtree[n,(1:f1)] %in% scored)<f1)	{
			anc_st <- cvector[ht]
			unscored <- mtree[n,!mtree[n,(1:f1)] %in% scored]
			unscored <- unscored[unscored>0]
			for (u in 1:length(unscored))	{
				if (length(unscored)>1)	{
					f2 <- unscored[u]
					}	else	{
					f2 <- unscored
					}
				# if parent value matches one of the possible values for the daughter node
				#	then go with that value
				if(anc_st>=0 && (!is.na(match(anc_st,obs_states[sankoff_matrix[f2,] %in% min(sankoff_matrix[f2,])]))))	{
					cvector[f2] <- anc_st
					}	else if (anc_st<0 && (anc_st!=UNKNOWN && anc_st!=INAP))	{
					anc_poss <- unravel_polymorph(anc_st)
					f2_poss <- unravel_polymorph(cvector[f2])
					reduced_poss <- f2_poss[f2_poss %in% anc_poss]
					if (length(reduced_poss)==1)	{
						cvector[ht] <- cvector[f2] <- reduced_poss
						}	else	{
	#					cvector[ht] <- cvector[f2] <- ravel_polymorph(reduced_poss)
						# if it is still up in the air at this point, then you might as
						#	well roll the dice!
						cvector[f2] <- reduced_poss[ceiling(runif(1)/(1/length(reduced_poss)))]
						}
					}	else	{
				# it really does not matter at this point what state we pick: it will all be the same.
					f2_poss <- unravel_polymorph(cvector[f2])
					cvector[f2] <- f2_poss[ceiling(runif(1)/(1/length(f2_poss)))]
					}
					# end case where we can assign ancestral condition to daughter node
				}	# end search of unscored nodes
			}	# I could add a routine here to compare sets of equally parsimonious states
		changes_above[n] <- min(sankoff_matrix[n+otus,])
		}		# if node has (0,1) and daughter node has (1,2), then go with 1
	}
cchanges <- rep(0,length(cvector))
for (n in 1:Nnode)	{
	ht <- otus+n
	if (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP)	{
		f1 <- node_rich[n]
		for (f in 1:f1)	{
			f2 <- mtree[n,f]
			if (cvector[f2] != cvector[ht])
				if ((cvector[f2]!=UNKNOWN && cvector[f2]!=INAP) && (cvector[ht]!=UNKNOWN && cvector[ht]!=INAP))	{
					# if both taxa are resolved, then this is easy
					if (cvector[f2]<0)	{
						ply <- unravel_polymorph(cvector[f2])
						if (is.na(match(cvector[ht],ply)))	cchanges[f2] <- match(min(ply),obs_states)
						}	else if (length(cvector[f2])==1 && length(cvector[ht])==1)	{
						cchanges[f2] <- match(cvector[f2],obs_states)
						}
					}
			}
		}
	}
output <- list(ch_steps,cvector,cchanges)
names(output) <- c("Steps","States","Derivation")
return(output)
}

#optimo_rate_and_tree_divergences_posterior_scaled_to_base
optimo_rate_and_tree_divergences_posterior_scaled_to_base <- function(base,m_alpha,divergence_times,vtree,notu,marginals,states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vtree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,notu=notu,marginals=marginals,states=states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vtree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times)
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times")
return(results)
}

#base=base_rep[bb];m_alpha_1=m_alpha_bang_rep;m_alpha_2=m_alpha_post_rep;divergence_times=init_divergence_times_rep;branch_rates=branch_rates_rep;vtree=vtree_ord;notu=notu_ord;marginals=marginals;states=states;BDS=BDS;rate_quants=rate_quants
optimo_rate_and_tree_divergences_posterior_scaled_to_base_two_rates_new <- function(base,m_alpha_1,m_alpha_2,divergence_times,branch_rates,vtree,notu,marginals,states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vtree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,notu=notu,marginals=marginals,states=states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
#ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
m_alphas <- c(m_alpha_1,m_alpha_2);
ddd <- optim(m_alphas,fn=likelihood_of_early_and_late_alphas_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,branch_rates=branch_rates,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vtree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times);
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times");
return(results)
}

optimo_rate_and_tree_divergences_posterior_scaled_to_base_rate_shift <- function(base,m_alpha_1,m_alpha_2,divergence_times,rate_shift,vtree,notu,marginals,states,BDS,rate_quants=1)	{
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
# node_ages <- -470+runif(nNodes)*30
#print(base);
first_sample <- min(divergence_times[1:notu]);
ur_b <- notu+1	# htu number of basal node
ttu <- max(vtree);
nNodes <- ttu-notu;
predates <- (1:ttu)[divergence_times<first_sample];
orig_base <- divergence_times[ur_b];
divergence_times[predates] <- base+abs(base-first_sample)*(1-(divergence_times[predates]-first_sample)/(orig_base-first_sample));

#sapply(m_alpha,accersi_tree_log_likelihood_simple,)
#ddd <- optim(m_alpha,fn=accersi_tree_log_likelihood_simple,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,notu=notu,marginals=marginals,states=states,lower=m_alpha/10,upper=1.0,control=cl)
cl <- list(fnscale=-1);
#ddd <- optim(m_alpha,fn=likelihood_of_alpha_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha/10,upper=2*m_alpha,control=cl);
m_alphas <- c(m_alpha_1,m_alpha_2);
ddd <- optim(m_alphas,fn=likelihood_of_rate_shift_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,rate_shift=rate_shift,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
#ddd <- optim(m_alphas,fn=likelihood_of_early_and_late_alphas_given_divergences,method="L-BFGS-B",divergence_times=divergence_times,branch_rates=branch_rates,vtree=vtree,marginals=marginals,states=states,rate_quants=rate_quants,lower=m_alpha_2/10,upper=2*m_alpha_1,control=cl);
ml_alpha <- ddd$par;
mlgln <- ddd$value;
mlgprior <- prior_probability_of_evolutionary_history_given_basal_divergence(divergence_times,vtree,BDS);
mlgposter <- mlgln+mlgprior;
results <- list(ml_alpha,mlgln,mlgprior,mlgposter,divergence_times);
names(results) <- c("ml_Rate","lnL","lnPrior","lnPosterior","divergence_times");
return(results)
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
likelihood_of_alpha_given_divergences <- function(m_alpha,divergence_times,vtree,marginals,states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
# states: number fo state per character
#
if (m_alpha < 0)
	m_alpha <- abs(m_alpha);
#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha <- m_alpha*rate_quants;
log_like <- sapply(alpha,calculate_rate_likelihood_given_tree_divergences_and_characters,vtree,divergence_times,marginals,states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# good for optimizing a single rate or a mean rate given a distribution (summarized in rate_quants)
# created 2020-07-11
likelihood_of_rate_shift_given_divergences <- function(m_alphas,divergence_times,rate_shift,vtree,marginals,states,rate_quants=1)	{
# m_alpha: typical rate for character change
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
# states: number fo state per character
#
m_alphas <- abs(m_alphas);
m_alpha_1 <- m_alphas[1];
m_alpha_2 <- m_alphas[2];

#alpha_sts <- alpha*1/((1:mxstates)-1);
alpha_1 <- m_alpha_1*rate_quants;
alpha_2 <- m_alpha_2*rate_quants;
alphas <- cbind(m_alpha_1*rate_quants,m_alpha_2*rate_quants);
#log_like_e <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_1,vtree=vtree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,states=states);
#log_like_u <- calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alpha_2,vtree=vtree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,states=states);

log_like <- apply(alphas,MARGIN=1,FUN=calculate_rate_shift_likelihood_given_tree_divergences_and_characters,vtree=vtree,divergence_times=divergence_times,rate_shift=rate_shift,marginals=marginals,states=states);
#lll <- vector(length=length(rate_quants));
#for (ii in 1:length(rate_quants))	{
#	lll[ii] <- calculate_rate_shift_likelihood_given_tree_divergences_and_characters(alphas=alphas[ii,],vtree,divergence_times,rate_shift,marginals,states);
#	}

#log_like <- c();
#for (rq in 1:nrow(alphas))
#	log_like <- c(log_like,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters(alphas=alphas[rq,],vtree=vtree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,states=states));
#	log_like <- sapply(alphas,calculate_early_and_late_rate_likelihood_given_tree_divergences_and_characters,vtree=vtree,divergence_times=divergence_times,branch_rates=branch_rates,marginals=marginals,states=states);
return(BayesFactor::logMeanExpLogs(log_like))
}

# created 2020-07-25
calculate_rate_shift_likelihood_given_tree_divergences_and_characters <- function(alphas,vtree,divergence_times,rate_shift,marginals,states,rate_quants=1)	{
# alphas: the first and second rate
# NOTE: Do this separately for all 4 quartiles if we have rate variation!!!
notu <- match(-1,vtree)-1;
tree_base <- notu+1;
ttu <- length(vtree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- c();
for (aa in 1:length(alphas))
	alpha_sts <- rbind(alpha_sts,alphas[aa]*1/((1:mxstates)-1));

for (nn in nNodes:1)	{
	htu <- nn+notu;
#	b_rate <- branch_rates[htu];
	daughters <- (1:ttu)[vtree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	node_onset <- divergence_times[htu];
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[states==k];
		mk <- array(0,dim=c(2,length(daughters)))
		for (dd in 1:length(daughters))	{
			f1 <- daughters[dd];
			f1_onset <- divergence_times[f1];
			if (node_onset < rate_shift && f1_onset > rate_shift)	{
				r1 <- abs(rate_shift-node_onset)/(f1_onset-node_onset);
				r2 <- 1-r1;
				} else if (f1_onset <= rate_shift)	{
				r1 <- 1;
				r2 <- 0;
				} else	{
				r1 <- 0;
				r2 <- 1;
				}
			br_alphas <- (r1*alpha_sts[1,k]) + (r2*alpha_sts[2,k]);
#			mk[,dd] <- 0
			mk[,dd] <- Mk_continuous(bd=bd[dd],k,br_alphas);
#				mk[,dd] <- (mk[,dd] + Mk_continuous(bd=bd[dd],k,br_alphas[bra])/length(br_alphas));
			}
		rc <- 0;
		while (rc < length(rchars))	{
			rc <- rc + 1;
			ch <- rchars[rc];
			marginals[htu,ch,1:k] <- 1;
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k]);
#			print(logmarginals[nn,ch,1:k])
			}
		}
#	print(sum(logmarginals[nn,,]));
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like);
}

# routine to get likelihood of one rate given tree & divergence times.  Good for optimizing one rate or calculating rates over lognormal/gamma distributions
calculate_rate_likelihood_given_tree_divergences_and_characters <- function(alpha,vtree,divergence_times,marginals,states)	{
notu <- match(-1,vtree)-1;
tree_base <- notu+1;
ttu <- length(vtree);
nNodes <- ttu-notu;
tnchars <- dim(marginals)[2];
mxstates <- dim(marginals)[3];
logmarginals <- array(0,dim=c(nNodes,tnchars,mxstates));
gap_lgprobs <- 0;
alpha_sts <- alpha*1/((1:mxstates)-1);
for (nn in nNodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vtree==htu];
	tomy <- length(daughters);
	bd <- abs(divergence_times[daughters]-divergence_times[htu]);
	for (k in 2:mxstates)	{
		rchars <- (1:tnchars)[states==k];
		mk <- sapply(bd,Mk_continuous,k,alpha=alpha_sts[k]);	# mk model for each node, with first value being P[net stasis]
		rc <- 0;
		while (rc < length(rchars))	{
#		for (rc in 1:length(rchars))	{
			rc <- rc+1;
			ch <- rchars[rc]
			marginals[htu,ch,1:k] <- 1
			for (d in 1:tomy)	{
#			d <- 0
#			while (d<=tomy)	{
#				d <- d+1
				dd <- daughters[d]	# dth taxon attached to node
				for (st in 1:k)	{
					# prob net stasis
					mg_st <- mk[1,d]*marginals[dd,ch,st]
					# prob net change
					mg_ch <- mk[2,d]*marginals[dd,ch,(1:k)[!(1:k) %in% st]]
					marginals[htu,ch,st] <- marginals[htu,ch,st]*(mg_st+sum(mg_ch))
					}	
				}
			logmarginals[nn,ch,1:k] <- log(marginals[htu,ch,1:k])
			}
		}
	}
basal_char_likelihoods <- rowSums(marginals[tree_base,,]);
log_like <- sum(log(basal_char_likelihoods));
return(log_like)
}

prior_probability_of_evolutionary_history_given_basal_divergence <- function(divergence_times,vtree,BDS)	{
# modified 2019-08-26
# divergence_times: vector giving 
# vtree: vector giving htu number of node from which each htu & otu evolved
#
notu <- match(-1,vtree)-1;
nNodes <- length(vtree)-notu;
gap_lgprobs <- 0;
for (nn in nNodes:1)	{
	htu <- nn+notu;
	daughters <- (1:ttu)[vtree==htu];
	for (d in 1:length(daughters))	{
		stgap <- (1:nrow(BDS))[BDS$ma > divergence_times[htu]];		# all intervals after divergence
		stgap <- stgap[BDS$ma[stgap] < divergence_times[daughters[d]]];	# intervals between divergence & appearance
		gap_lgprobs <- gap_lgprobs+sum(log(BDS$pgap[stgap]));
		}
	}
return(gap_lgprobs);
}

allocate_node_marginals <- function(chmatrix,nttu,states)	{
# chmatrix: original character matrix;
# nttu: total taxa + nodes;
# states: states per character
mxstates <- max(states);
tnchars <- ncol(chmatrix);
notu <- nrow(chmatrix);
marginals <- array(0,dim=c(nttu,tnchars,mxstates));
states[states<2] <- 2;
for (n in 1:notu)	{
	for (c in 1:tnchars)	{
		if (chmatrix[n,c]>=0)	{
			marginals[n,c,chmatrix[n,c]] <- 1
			} else if (chmatrix[n,c]==UNKNOWN || chmatrix[n,c]==INAP)	{
			marginals[n,c,1:states[c]]	<- 1/states[c];
			} else {
			pc <- -1*chmatrix[n,c];
			polys <- vector(length=ceiling(log10(pc+0.5)));
			ex <- length(polys);
			for (i in 1:length(polys))	{
				polys[i] <- floor(pc/(10^(ex-1)));
				pc <- pc-(polys[i]*10^(ex-1));
				ex <- ex-1;
				}
#			polys <- polys+1;
			marginals[n,c,polys] <- 1/length(polys)
			}
		}
	}
return(marginals)
}

Mk_continuous <- function(bd,k,alpha)	{
# bd: time (branch duration)
# k: states
# alpha: instantaneous rate
# from Lewis 2001
pstasis <- (1/k)+((k-1)*exp(-k*alpha*bd)/k)
pchange <- (1/k)-(exp(-k*alpha*bd)/k)
return(c(pstasis,pchange))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### basic "stratocladistic" routines ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# get minimum implied gap along each branch
accersi_gaps_from_vector_tree <- function(vtree,strat_ranges,apos=NULL)	{
if (is.matrix(strat_ranges))	strat_ranges <- data.frame(strat_ranges,stringsAsFactors = F);
if (ncol(strat_ranges)==2)	{
	colnames(strat_ranges) <- c("FA","LA");
	} else if (ncol(strat_ranges)==3)	{
	colnames(strat_ranges) <- c("taxon","FA","LA");
	}
ohtu <- length(vtree);
notu <- nrow(strat_ranges);
if (is.null(apos))
	apos <- rep(1,ohtu);
base <- notu+1;
dates <- date_clades_from_vector_tree(vtree,FA=strat_ranges$FA);
gaps <- vector(length=ohtu);
node_poss_anc <- accersi_poss_ancestors_for_nodes(vtree,FA=strat_ranges$FA,apos);

for (tx in 1:ohtu)	{
	if (tx==base)	tx <- tx+1
	if (node_poss_anc[vtree[tx]]==0)	{
		gaps[tx] <- abs(dates[tx]-dates[vtree[tx]])
		} else	{
		gaps[tx] <- max(0,dates[tx]-strat_ranges$LA[node_poss_anc[vtree[tx]]])
		}
	}
return(gaps);
}

# collect node ages from vector tree
date_clades_from_vector_tree <- function(vtree,FAs)	{
ohtu <- length(vtree)
notu <- length(FAs)
base <- notu+1
dates <- vector(length=ohtu)
end <- max(FAs)
for (nd in base:ohtu)	dates[nd] <- end
for (sp in 1:notu)	{
	if (vtree[sp]>0)	{
		dates[sp] <- FAs[sp]
		anc <- vtree[sp]
		if (dates[anc]>dates[sp])	dates[anc] <- dates[sp]
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
for (nd in ohtu:(base+1))	{
	if (vtree[nd]>0)	{
		anc <- vtree[nd]
		if (dates[anc]>dates[nd])	dates[anc] <- dates[nd]
		}	# this is to make sure that species excluded from the tree do not mess up analysis
	}
return(dates)
}

# get possible ancestors given apomorphies
accersi_poss_ancestors_for_nodes <- function(vector_tree,FA,apos)	{
ohtu <- length(vector_tree);
notu <- length(FA);
node_poss_anc <- rep(0,ohtu);
if (!is.na(apos[1]))
	for (n in 1:notu)
		if (apos[n]==0)	
			if (node_poss_anc[vector_tree[n]]==0)	{
				node_poss_anc[vector_tree[n]] <- n;
				} else if (FA[node_poss_anc[vector_tree[n]]] > FA[n])
				FA[node_poss_anc[vector_tree[n]]] <- FA[n];
return(node_poss_anc)
}

accersi_starter_divergences <- function(vtree,simple_clock,strat_ranges,prec=0.1)	{
# vtree: vector giving htu from which each otu & htu evovled
# simple_clock: vector giving expected branch durations given overall rates
# strat_ranges: first & last appearances of species
if (is.matrix(strat_ranges))
	strat_ranges <- data.frame(strat_ranges,stringsAsFactors = F);
if (sum(colnames(strat_ranges) %in% c("fa","FA","fad","FAD","Ma_LB","onset"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("fa","FA","fad","FAD","Ma_LB","onset")];
	col_titles[fix_col] <- "ma_lb";
	colnames(strat_ranges) <- col_titles;
	}
if (sum(colnames(strat_ranges) %in% c("la","LA","Lad","LAD","Ma_UB","end"))>0)	{
	col_titles <- colnames(strat_ranges);
	fix_col <- (1:ncol(strat_ranges))[colnames(strat_ranges) %in% c("la","LA","Lad","LAD","Ma_UB","end")];
	col_titles[fix_col] <- "ma_ub";
	colnames(strat_ranges) <- col_titles;
	}

notu <- min(vtree[vtree>0])-1;
nttu <- max(vtree);
nNodes <- nttu-notu;
node_ages_starters <- c(-abs(strat_ranges$ma_lb[1:notu]),rep(0,nNodes));
nn <- nNodes;
#while (nn>0)	{
for (nn in nNodes:1)	{
	htu <- notu+nn;							# get vector number for the node's age
	daughters <- (1:nttu)[vtree==htu];		# get descendants of this node
	if (max(simple_clock[daughters])==0)
		simple_clock[daughters[length(daughters)]] <- prec;
	local_clock <- node_ages_starters[daughters]-simple_clock[daughters];	# get divergence times implied by simple clock
	# make sure that local clock is at least as old minimum divergence times
	local_clock[local_clock>min(node_ages_starters[daughters])] <- min(node_ages_starters[daughters]);
	node_ages_starters[htu] <- prec*round(mean(local_clock)/prec,0);
#	nn <- nn-1
#	daughters
#	node_ages_starters[daughters]
#	node_ages_starters
	}
return(node_ages_starters);
}

# get node ages using minimum necessary ages
date_taxa_on_tree_simple <- function(vtree,FAs)	{
ttus <- max(vtree)	# total htus+otus
dates <- vector(length=ttus)
notu <- length(FAs)	# number of otus
if (length(dim(as.array(vtree)))==1)	{
	# routine if tree given as a single vector with ancestral
	base <- notu+1
	end <- max(FAs)
	for (nd in base:ttus)	dates[nd] <- end
	for (sp in 1:notu)	{
		if (vtree[sp]>0)	{
			dates[sp] <- FAs[sp]
			anc <- vtree[sp]
			if (dates[anc]>dates[sp])	dates[anc] <- dates[sp]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	for (nd in ttus:(base+1))	{
		if (vtree[nd]>0)	{
			anc <- vtree[nd]
			if (dates[anc]>dates[nd])	dates[anc] <- dates[nd]
			}	# this is to make sure that species excluded from the tree do not mess up analysis
		}
	} else	{
	Nnode <- dim(vtree)[1]
	dates[1:notu] <- FAs
	for (n in Nnode:1)	{
		f1 <- sum(vtree[n,]>0)
		ht <- n+notu
		dates[ht] <- min(dates[vtree[n,(1:f1)]])
		}
	}
return(dates)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Distribution Analysis ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# find expectations of this gamma at this sample size
# exp will give the expected number of taxa found 1…ncoll times
expected_abundances <- function(rel_ab_dist, nspec, S)	{
expected <- vector(length=nspec)
for (t in 1:S)	{
	for (i in 1:nspec)	{
		expected[i] <- expected[i]+dbinom(i,nspec,rel_ab_dist[t])
		}
	}
return(expected)
}

optimize_uniform_abundance <- function(counts)	{
# written 2017-01-28
hS <- minS <- length(counts)
nspec <- sum(counts)
oS <- sum(counts>0);
#observed <- fisher_plot(counts);
observed <- hist(counts[counts>0],breaks=0:max(counts),plot=F)$counts

mxlnl <- lnl <- -1*MAXNO
while (lnl == mxlnl)	{
	rel_ab_dist <- rep(1/hS,hS)
	raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	if (mxlnl < lnl)	mxlnl <- lnl
	hS <- hS+1
	}
bH <- c(hS-1,round(mxlnl,3),round(modified_AIC(mxlnl,1,nspec),3))
names(bH) <- c("Uniform_S","Uniform_log-likelihood","Uniform_AICc")
return(bH)
}

# generate geometric (= exponential) distribution
geometric_distribution <- function(decay)	{
if (decay>1)	decay <- (1/decay)
S <- round(1+((log(10^-9)-log(1-decay))/log(decay)))
ranks <- (1:S)
prop <- (1-decay)*(decay^(ranks-1))
return(prop)
}

# get log-likelihood of geometric given decay
loglikelihood_geometric_rad <- function(decay,nspec,oS,observed)	{
# p0[1]: r	# p0[2]: decay	# p0[3]: S
rel_ab_dist <- geometric_distribution(decay)		# basic exponential distribution
hS <- length(rel_ab_dist)
#print(c(round(decay,5),hS))		# for debugging
if (hS>=oS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,hS)
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	}	else	{
	lnl <- oS*log(MINNO)
	}
#print(c(round(decay,5),round(lnl,10)))		# for debugging
return(lnl)
}

# get best-fit geometric
optimize_geometric_abundance <- function(counts)	{
oS <- length(counts)				# observed taxa
#observed <- fisher_plot(counts);
observed <- hist(counts[counts>0],breaks=0:max(counts),plot=F)$counts
decay <- (min(counts)/max(counts))^(1/(length(counts)-1))
max_decay <- exp(log(decay)/2)
min_decay <- exp(log(decay)*2)

cl <- list(fnscale=-1)
nspec <- sum(counts)
w <- optim(decay,fn=loglikelihood_geometric_rad,method="L-BFGS-B",nspec=nspec,oS=oS,observed=observed,lower=min_decay,upper=max_decay,control=cl)
bH <- c(round(w$par,6),round(w$value,3),round(modified_AIC(w$value,1,nspec),3))
names(bH) <- c("Geometric_decay","Geometric_log-likelihood","Geometric_AICc")
return(bH)
}

# get log-likelihood for lognormal
loglikelihood_lognormal_rad_for_optim <- function(oS,nspec,rand_no,hS,observed,min_mag,max_mag)	{
mag <- min_mag + rand_no*(max_mag - min_mag)
#print(c(round(mag,9)))
rel_ab_dist <- lognormal_distribution(mag=mag,S=hS)		# basic lognormal distribution
raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
return(round(lnl,2))
}

# get minimum and maximum magnitude of increase for a lognormal given hS
accersi_min_and_max_lognormal_mag_given_hS <- function(observed,counts,hS)	{
oS <- sum(observed)				# observed taxa
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
rank_shifts <- sort(qnorm((1:hS)/(hS+1)),decreasing=TRUE)[1:oS]
local_m <- mag_shifts <- vector(length=(oS-1))
for (i in 1:(oS-1))		{
	mag_shifts[i] <- counts[i]/counts[i+1]
	local_m[i] <- exp(log(mag_shifts[i])/(rank_shifts[i]-rank_shifts[i+1]))
	}
if (max(observed)==1)	{
	min_mag <- min(local_m)
	}	else {
	min_mag <- 1.01
	}
max_mag <- prod(local_m)^(1/(length(local_m)))
if (max_mag < min_mag)	{
	m <- min_mag
	min_mag <- max_mag
	max_mag <- m
	}
return(c(min_mag,max_mag))
}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_lognormal_abundance_given_hS <- function(observed,counts,hS)	{
#print(hS)		# for debugging
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
mxfinds <- length(observed)
mnfinds <- min((1:mxfinds)[!observed %in% 0])
mm <- accersi_min_and_max_lognormal_mag_given_hS(observed,counts,hS)
mag <- prod(mm)^0.5
#inev <- mag <- exp(log(mxfinds/mnfinds)/(qnorm((hS-1)/(hS+1))-qnorm((hS-oS)/(hS+1))))
cl <- list(fnscale=-1)
rand_no <- (mag - mm[1])/(mm[2] - mm[1])
#w <- optim(rand_no,loglikelihood_lognormal_rad_for_optim,,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,min_mag=mm[1],max_mag=mm[2],lower=0,upper=1,control=cl)
w <- optim(rand_no,loglikelihood_lognormal_rad_for_optim,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,min_mag=mm[1],max_mag=mm[2],lower=0,upper=1,control=cl);
w$par <- mm[1] + (w$par*(mm[2] - mm[1]))
#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=mm[1],upper=mm[2],control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("Lognormal_magnitude", "Lognormal_S","Lognormal_log-likelihood")
return(bH)
}

# get most likely lognormal distribution for a population of some sort
# modified 2020-07-01 to allow unobserved but present 0 finds
optimize_lognormal_abundance <- function(counts)	{
minS <- stS <- length(counts)				# observed taxa
oS <- sum(counts>0);
observed <- hist(counts[counts>0],breaks=0:max(counts),plot=F)$counts
nspec <- sum(counts)
span <- 5
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
while (incr>0)	{
	hS <- seq(stS,enS,by=incr)
	results <- sapply(hS,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts[counts>0])
	mlnl <- max(results[3,])
	mlSc <- match(mlnl,results[3,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[3,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				stS <- max(minS,hS[spn]);
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]
					incr <- 1
					}
				stS <- max(minS,(hS[spn]-(span*incr)));
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- max(minS,hS[mlSc]);
			enS <- hS[mlSc+1]-1;
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- max(minS,hS[mlSc-1]+1);
				enS <- hS[mlSc]
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- max(minS,hS[mlSc]);
				enS <- hS[mlSc+1]-1
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	}
bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
names(bH)[4] <- "Lognormal_AICc"
return(bH)
}

# get gamma distribution
gamma_distribution <- function(a, b, S)	{
p0 <- vector(length=S)
p0[1] <- S/(S+1)
for(i in 2:S) p0[i] <- p0[i-1]-(1/(S+1))
prop <- qgamma(p0,a,b)/sum(qgamma(p0,a,b))
return(prop)
}

# get minimum alpha
accersi_minimum_alpha_for_gamma_one <- function(hS)	{
min_alpha <- 1
iS <- hS
while (iS==hS)	{
	rel_ab_dist <- gamma_distribution(a=min_alpha,b=min_alpha,S=hS)		# basic lognormal distribution
	iS <- sum(rel_ab_dist>MINNO)
	if (iS == hS)	min_alpha <- min_alpha/2
	}
return(2*min_alpha)
}

# get log-likelihood for gamma
loglikelihood_gamma_one_rad <- function(oS,nspec,alpha,hS,observed)	{
#print(c(log10(alpha),hS))
rel_ab_dist <- gamma_distribution(a=alpha,b=alpha,S=hS)		# basic lognormal distribution
#print(rel_ab_dist)
if (sum(rel_ab_dist>MINNO) == hS)	{
	raw_expected <- expected_abundances(rel_ab_dist,nspec,S=hS)
	#lnl <- abundance_distribution_loglikelihood_suf(observed,expected=raw_expected,oS=oS,hS=hS)
	lnl <- distribution_loglikelihood_mul(observed,expected=raw_expected,oS,hS)
	}	else	lnl <- -1*MAXNO
return(round(lnl,2))
}

# get most likely magnitude of increase for a lognormal distribution of hS entities
optimize_gamma_one_abundance_given_hS <- function(observed,counts,hS)	{
#print(hS)		# for debugging
oS <- sum(observed)				# observed taxa
nspec <- sum((1:length(observed))*observed)
cl <- list(fnscale=-1)
alpha <- 2
min_alpha <- accersi_minimum_alpha_for_gamma_one(hS)
w <- optim(alpha,loglikelihood_gamma_one_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=min_alpha,upper=10000,control=cl)
#w <- optim(mag,fn=loglikelihood_lognormal_rad,method="L-BFGS-B",oS=oS,nspec=nspec,hS=hS,observed=observed,lower=mm[1],upper=mm[2],control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("Gamma_magnitude", "Gamma_S","Gamma_log-likelihood")
return(bH)
}

# get most likely gamma distribution for a population of some sort
optimize_gamma_one_abundance <- function(counts)	{
#stS <- oS <- length(counts)				# observed taxa
#observed <- fisher_plot(counts)
minS <- stS <- length(counts)				# observed taxa
oS <- sum(counts>0);
observed <- hist(counts[counts>0],breaks=0:max(counts),plot=F)$counts
nspec <- sum(counts)
span <- 5
enS <- stS+((span-1)*oS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
while (incr>0)	{
	hS <- seq(stS,enS,by=incr)
	results <- sapply(hS,optimize_gamma_one_abundance_given_hS,observed=observed,counts=counts)
	mlnl <- max(results[3,])
	mlSc <- match(mlnl,results[3,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[3,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				stS <- max(minS,hS[spn]);
				enS <- stS+((spn-1)*oS)
				}	else	{
				enS <- hS[spn]
				if (span<incr)	{
					incr <- floor((enS-hS[spn-1])/span)
					}	else {
					span <- enS-hS[spn-1]
					incr <- 1
					}
				stS <- max(minS,(hS[spn]-(span*incr)));
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- max(minS,hS[mlSc]);
			enS <- hS[mlSc+1]-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- floor((enS-stS)/span)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_abundance_given_hS,observed=observed,counts=counts)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- max(minS,(hS[mlSc-1]+1));
				enS <- hS[mlSc]
				if (span<incr)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- max(minS,hS[mlSc]);
				enS <- hS[mlSc+1]-1
				if (incr>span)	{
					incr <- floor((enS-stS)/span)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	}
bH <- c(round(bH[1],6),bH[2],bH[3],modified_AIC(mlnl,2,nspec))
names(bH)[4] <- "Gamma_AICc"
return(bH)
}



