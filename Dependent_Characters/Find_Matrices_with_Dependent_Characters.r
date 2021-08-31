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
matrices_w_inaps <- c();
for (bm in 1:length(brachiopod_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=brachiopod_matrices[bm],character_database);
	inap_cases <- unique(which(character_info$Matrix==INAP,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(brachiopod_matrices[bm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
brachiopod_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (bm in 1:length(brachiopod_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=brachiopod_matrices_inaps[bm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==brachiopod_matrices_inaps[bm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		for (cd in (dependent_chars[dc]-1):1)	{
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}

if (is.null(inapplicable_examples$independent_character))	{
	inapplicable_examples <- tibble::add_column(inapplicable_examples,independent_character=as.numeric(independents),.after = match("dependent_character",colnames(inapplicable_examples)));
	} else	{
	inapplicable_examples$independent_character=as.numeric(independents);
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Brachiopods_w_Inapplicables.csv",sep=""),row.names = F);

trilobite_matrices <- names(character_database$Panarthropods$Trilobites);
matrices_w_inaps <- c();
for (tm in 1:length(trilobite_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=trilobite_matrices[tm],character_database);
	inap_cases <- unique(which(character_info$Matrix==-22,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(trilobite_matrices[tm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
trilobite_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (tm in 1:length(trilobite_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=trilobite_matrices_inaps[tm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==trilobite_matrices_inaps[tm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		for (cd in (dependent_chars[dc]-1):1)	{
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}

if (is.null(inapplicable_examples$independent_character))	{
	inapplicable_examples <- tibble::add_column(inapplicable_examples,independent_character=as.numeric(independents),.after = match("dependent_character",colnames(inapplicable_examples)));
	} else	{
	inapplicable_examples$independent_character=as.numeric(independents);
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Trilobites_w_Inapplicables.csv",sep=""),row.names = F);

cephalopod_matrices <- names(character_database$Molluscs$Cephalopods);
matrices_w_inaps <- c();
for (cm in 1:length(cephalopod_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=cephalopod_matrices[cm],character_database);
	inap_cases <- unique(which(character_info$Matrix==-22,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(cephalopod_matrices[cm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
cephalopod_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (cm in 1:length(cephalopod_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=cephalopod_matrices_inaps[cm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==cephalopod_matrices_inaps[cm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		for (cd in (dependent_chars[dc]-1):1)	{
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}

if (is.null(inapplicable_examples$independent_character))	{
	inapplicable_examples <- tibble::add_column(inapplicable_examples,independent_character=as.numeric(independents),.after = match("dependent_character",colnames(inapplicable_examples)));
	} else	{
	inapplicable_examples$independent_character=as.numeric(independents);
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Cephalopods_w_Inapplicables.csv",sep=""),row.names = F);

bivalve_matrices <- names(character_database$Molluscs$Bivalves);
matrices_w_inaps <- c();
for (cm in 1:length(bivalve_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=bivalve_matrices[cm],character_database);
	inap_cases <- unique(which(character_info$Matrix==-22,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(bivalve_matrices[cm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
bivalve_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (cm in 1:length(bivalve_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=bivalve_matrices_inaps[cm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==bivalve_matrices_inaps[cm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		for (cd in (dependent_chars[dc]-1):1)	{
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}

if (is.null(inapplicable_examples$independent_character))	{
	inapplicable_examples <- tibble::add_column(inapplicable_examples,independent_character=as.numeric(independents),.after = match("dependent_character",colnames(inapplicable_examples)));
	} else	{
	inapplicable_examples$independent_character=as.numeric(independents);
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Bivalves_w_Inapplicables.csv",sep=""),row.names = F);

gastropod_matrices <- names(character_database$Molluscs$Gastropods);
matrices_w_inaps <- c();
for (gm in 1:length(gastropod_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=gastropod_matrices[gm],character_database);
	inap_cases <- unique(which(character_info$Matrix==-22,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(gastropod_matrices[gm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
gastropod_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (gm in 1:length(gastropod_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=gastropod_matrices_inaps[gm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==gastropod_matrices_inaps[gm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		for (cd in (dependent_chars[dc]-1):1)	{
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}

if (is.null(inapplicable_examples$independent_character))	{
	inapplicable_examples <- tibble::add_column(inapplicable_examples,independent_character=as.numeric(independents),.after = match("dependent_character",colnames(inapplicable_examples)));
	} else	{
	inapplicable_examples$independent_character=as.numeric(independents);
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Gastropods_w_Inapplicables.csv",sep=""),row.names = F);

mammal_matrices <- names(character_database$Vertebrates$Mammals);
matrices_w_inaps <- c();
for (bm in 1:length(mammal_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=mammal_matrices[bm],character_database);
	inap_cases <- unique(which(character_info$Matrix==INAP,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(mammal_matrices[bm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
mammal_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (bm in 1:length(mammal_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=mammal_matrices_inaps[bm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==mammal_matrices_inaps[bm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		cd <- dependent_chars[dc];
		while (cd > 1)	{
#		for (cd in (dependent_chars[dc]-1):1)	{
			cd <- cd-1;
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Mammals_w_Inapplicables.csv",sep=""),row.names = F);

dinosaur_matrices <- names(character_database$Vertebrates$Dinosaurs);
matrices_w_inaps <- c();
for (bm in 1:length(dinosaur_matrices))	{
	character_info <- accersi_data_from_RData(matrix_name=dinosaur_matrices[bm],character_database);
	inap_cases <- unique(which(character_info$Matrix==INAP,arr.ind=T)[,2]);
	if (length(inap_cases)>0)
		matrices_w_inaps <- rbind(matrices_w_inaps,cbind(rep(dinosaur_matrices[bm],length(inap_cases)),inap_cases));
	}
inapplicable_examples <- data.frame(matrix=as.character(matrices_w_inaps[,1]),dependent_character=as.numeric(matrices_w_inaps[,2]))
dinosaur_matrices_inaps <- unique(inapplicable_examples$matrix);
independents <- c();
for (bm in 1:length(dinosaur_matrices_inaps))	{
	character_info <- accersi_data_from_RData(matrix_name=dinosaur_matrices_inaps[bm],character_database);
	chmatrix <- character_info$Matrix;
	dependent_chars <- inapplicable_examples$dependent_character[inapplicable_examples$matrix==dinosaur_matrices_inaps[bm]];
	for (dc in 1:length(dependent_chars))	{
		dcc <- dependent_chars[dc];
		independents <- c(independents,0);
		cd <- dependent_chars[dc];
		while (cd > 1)	{
#		for (cd in (dependent_chars[dc]-1):1)	{
			cd <- cd-1;
			unique_combos <- unique(chmatrix[,c(cd,dcc)]);
			unique_combos <- unique_combos[!unique_combos[,1] %in% c(UNKNOWN,INAP),];
			if (is.matrix(unique_combos))
				unique_combos <- unique_combos[!unique_combos[,2] %in% UNKNOWN,];
			if (is.matrix(unique_combos))	{
				g_c <- unique_combos[unique_combos[,2]==INAP,];
				gapped_combos <- array(g_c,dim=c((length(g_c)/2),2));
				ungapped_combos <- array(unique_combos[!unique_combos[,2] %in% INAP,],dim=c((length(unique_combos[!unique_combos[,2] %in% INAP,])/2),2));
				if (length(unique(gapped_combos[,1]))==1 && sum(ungapped_combos[,1]==gapped_combos[1,1])==0)	{
					independents[length(independents)] <- cd;
					cd <- 1;
#					} else if (is.matrix(gapped_combos) && length(unique(gapped_combos[,1])))	{
					}
				}
			}
		}
	}
write.csv(inapplicable_examples,paste(dependent_directory,"Dinosaurs_Inapplicables.csv",sep=""),row.names = F,fileEncoding = "UTF-8");
write.csv(unique(inapplicable_examples$matrix),paste(dependent_directory,"Dinosaurs_w_Inapplicables.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

{}