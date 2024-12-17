# brute force kluges to make R what I effin' well want it to do..... 
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal

#### send NAs to hell where they belong..... ####
clear_na_from_matrix <- function(data, replacement)  {
#size <- dim(data)
#for (i in 1:size[1])	{
#	for (j in 1:size[2]) if (is.na(data[i,j]))	data[i,j] <- replacement
#	}
for (i in 1:ncol(data))	data[is.na(data[,i]),i] <- replacement;
#	if(sum(is.na(data[,i]))>0)	{
#		duds <- (1:nrow(data))[is.na(data[,i])];
#		data[duds,i] <- replacement;
#		}
#	}
return(data)
}

clear_matrix_na_with_another_cell_value <- function(data,j, k)	{
size <- dim(data)
for (i in 1:size[1])	{
	if (is.na(data[i,j]))	data[i,j] <- data[i,k]
	}
return(data)
}

clear_na_from_vector <- function(data, replacement)	{
size <- length(data)
for (i in 1:size[1])	if (is.na(data[i]))	data[i] <- replacement
return(data)
}

#### Non-Asshole conversion therapy ####
matrix_to_list_by_columns <- function(data_mat)	{
output <- tapply(data_mat,rep(1:ncol(data_mat),each=nrow(data_mat)),function(i)i);
if (length(colnames(data_mat))>0)	
	names(output) <- colnames(data_mat);
return(output);
}

matrix_to_list_by_rows <- function(data_mat)	{
output <- tapply(data_mat,rep(1:nrow(data_mat),each=ncol(data_mat)),function(i)i);
if (length(colnames(data_mat))>0)	
	names(output) <- rownames(data_mat);
return(output);
}

#### Text Cleaning funcitons ####
mundus_web_text <- function(web_text)	{
web_text <- gsub("–","-",web_text);
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
web_text <- gsub("é","e",web_text);
web_text <- gsub("√©","e",web_text);
#web_text <- gsub("ö","&ouml;",web_text);
#web_text <- gsub("ü","&uuml;",web_text);
web_text <- gsub("ü","u",web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("ý","y",web_text);
web_text <- gsub("√Ω","y",web_text);
return(web_text);
}

frak_it <- function(x)	{
return(as.numeric(as.character(x)))
}

frak_it_character <- function(x)	{
if (is.factor(x)) {
	return(as.character(x));
	}	else {
	return(x);
	}
}

transmogrify_to_title_case <- function(name) {
name_part <- strsplit(name," ")[[1]];
return (paste(toupper(substring(name_part, 1,1)), substring(name_part, 2),sep="", collapse=" "));
}

transmogrify_diacritics <- function(funky_text)	{
j <- strsplit(as.character(funky_text),split="",fixed=TRUE)[[1]];
j <- j[!j %in% c("́")];
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
eek <- c("ć","ĉ","ċ","č","ƈ","ç","ḉ","ȼ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "c";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "C";
eek <- c("ḋ","ɗ","ḍ","ḏ","ḑ","ḓ","ď");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "d";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "D";
eek <- c("è","é","ê","ẽ","ē","ĕ","ė","ë","ẻ","ě","ȅ","ȇ","ẹ","ȩ","ę","ḙ","ḛ","ề","ế","ễ","ể","ḕ","ḗ","ệ","ḝ");
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
eek <- c("ĺ","ľ","ŀ","ł","ƚ");
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
eek <- c("ò","ó","ô","õ","ō","ŏ","ȯ","ö","ỏ","ő","ǒ","ȍ","ȏ","ơ","ǫ","ọ","ɵ","ø");
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
eek <- c("ù","ú","û","ũ","ū","ŭ","ü","ủ","ű","ǔ","ȕ","ȗ","ụ","ṳ","ų","ṷ","ṵ");
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

expello_na_from_matrix <- function(data, replacement="")  {
#i <- 0;
for (i in 1:ncol(data))	{
#	i <- i+1;
	if(sum(is.na(data[,i]))>0)	{
		duds <- (1:nrow(data))[is.na(data[,i])]
		data[duds,i] <- replacement
		}
#	print(paste(i,sum(is.na(data[,i]))))
	}
return(data)
}

expello_na_from_vector <- function(data, replacement="")	{
if(sum(is.na(data))>0)	{
	duds <- (1:length(data))[is.na(data)]
	data[duds] <- replacement
	}
return(data)
}

# written 2021-04-04 to help make faster rock search function
find_row_with_value <- function(go_seek,hide)	{
return(unique(which(hide==go_seek,arr.ind = T)[,1]));
}

match_value_to_matrix_cell <- function(go_seek,hide)	{
row_matches <- c()
for (cc in 1:ncol(hide))	{
	rc <- match(go_seek,hide[,cc])
	if (!is.na(rc))
		row_matches <- c(row_matches,rc)
	}
if (is.null(row_matches))
	row_matches <- "";	
row_matches <- unique(row_matches);
return(row_matches);
}

insert_cell_into_vector_x <- function(x,new_value,cell_no)	{
# routine to put take vector x and add a new value at some point (cell_no)
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (cell_no==1)	{
	return(c(new_value,x));
	} else if (cell_no==(length(x)+1))	{
	return(c(x,new_value));
	} else	{
	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

insert_row_into_matrix_x <- function(x,new_row,row_no)	{
# routine to put take matrix x and add a new row at that bumps other rows down one
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (row_no==1)	{
	xx <- rbind(new_row,x);
	colnames(xx) <- colnames(x);
	return(xx);
	} else if (row_no==(nrow(x)+1))	{
	return(rbind(x,new_row));
	} else	{
	return(rbind(x[1:(row_no-1),],new_row,x[(row_no:nrow(x)),]));
#	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

#find_me <- paleodb_clean_member_basic_no_rock[need_info_member[3]]
matrix_row_match <- function(find_me,examined_matrix)	{
xxx <- which(examined_matrix==find_me,arr.ind=TRUE);
#examined_matrix[unique(xxx[,1]),]
if (nrow(xxx)==0)	{
	return(-1);
	} else if (length(unique(xxx[,1]))==1)	{
	return(xxx[1,1])
	} else {
	vvv <- unique(xxx[,1]);
	mx <- 0;
	mj <- 0;
	for (j in 1:length(vvv))	{
		if (sum(vvv[j]==xxx[,1])>mx)	{
			mj <- j;
			mx <- sum(vvv[j]==xxx[,1])
			} else if (sum(vvv[j]==xxx[,1])==mx)	{
			mj <- c(mj,j)
			}
		}
	if (length(mj)==1)	{
		return(vvv[mj]);
		} else {
		return(-2);
		}
	}
}

#match_matrix <- genus_taxonomy; match_vector <- reef_builder_clades;
match_vector_to_matrix <- function(match_matrix,match_vector)	{
all_matches <- c();
for (mv in 1:length(match_vector))  {
  add_these <- which(match_matrix == match_vector[mv],arr.ind=T);
  if (length(add_these)==0) {
    add_these <- c(mv,0,0);
    add_these <- array(c(mv,0,0),dim=c(1,3));
    rownames(add_these) <- match_vector[mv];
    } else if (length(add_these)==1)  {
    add_these <- array(c(mv,add_these),dim=c(1,3));
    rownames(add_these) <- match_vector[mv];
    } else {
    add_these <- cbind(rep(mv,nrow(add_these)),add_these);
    rownames(add_these) <- rep(match_vector[mv],nrow(add_these));
    }
  all_matches <- rbind(all_matches,add_these);
  }
colnames(all_matches) <- c("vector_cell","matrix_row","matrix_column");
return(data.frame(all_matches));
}

#character_array <- strsplit(rock_unit_name,split = "");
simplify2vector <- function(character_array)	{
return(simplify2array(character_array)[,1])
}

revelare_file_type <- function(filename)	{
file_info <- stringr::str_split(filename,pattern="")[[1]];
f_i <- length(file_info);
f_i_s <- 1+(1:f_i)[file_info %in% "."];
file_type <- c();
for (i in f_i_s:f_i)	file_type <- paste(file_type,file_info[i],sep="");
return(file_type);
}

# transform dataframe into a list based on some criterion.
#	list_criterion gives the column name with the variables that will be used to divide the data.frame into a list
listifor <- function(dataframe_to_transform,list_criterion)	{
criterion_col <- match(list_criterion,colnames(dataframe_to_transform));
list_categories <- sort(unique(dataframe_to_transform[,criterion_col]));
list_categories <- list_categories[list_categories!=""];

listed_dataframe <- list();
for (lc in 1:length(list_categories))	{
	this_subset <- subset(dataframe_to_transform,dataframe_to_transform[,criterion_col]==list_categories[lc]);
	listed_dataframe <- rlist::list.append(listed_dataframe,this_subset);
	}
names(listed_dataframe) <- list_categories;
return(listed_dataframe);
}

choose_with_logs <- function(N,n)	{
if (n==0 || n==N)	{
	return(1);
	} else	{
	return(exp(sum(log(1:N))-(sum(log(1:n))+sum(log(1:(N-n))))));
	}
}

delete_row_from_matrix <- function(matrx,rw)	{
if (rw==1)	{
	return(matrx[2:nrow(matrx),])
	} else if (rw==nrow(matrx))	{
	return(matrx[1:(nrow(matrx)-1),])
	} else	{
	z <- rw+1;
	new_matrx <- matrx[1:(rw-1),];
	new_matrx <- rbind(new_matrx,matrx[(rw+1):nrow(matrx),]);
	return(new_matrx);
	}
}

#dataframe <- wagner_rocks;
tolower_dataframe <- function(dataframe)	{
for (nc in 1:ncol(dataframe))	{
	if (is.character(dataframe[,nc]))	dataframe[,nc] <- tolower(dataframe[,nc]);
	}
return(dataframe)
}

#### Matrix Summaries ####
find_gaps_in_vector <- function(vv)	{
gaps <- c();
for (v in 2:length(vv))	gaps <- c(gaps,(vv[v]-vv[v-1])-1)
return(gaps);
}

# matching vectors to rows in a matrix
#row1 <- 202; matrix2 <- cooccr_mat; matrix1 <- reduced_cooccr_mat;
accersi_matching_row <- function(row2,matrix2,matrix1)	{
test_vector <- matrix2[row2,];
return(match_vector_to_matrix_row(test_vector,matrix1))
}

#test_matrix <- reduced_cooccr_mat; test_vector <- cooccr_mat[19,]
match_vector_to_matrix_row <- function(test_vector,test_matrix)	{
rn <- 1:nrow(test_matrix);
xx <- sapply(rn,vector_matches_in_matrix,test_vector,test_matrix);
return(match(ncol(test_matrix),xx));
}

# vector2 <- test_vector; vector1 <- test_matrix[1,]
vector_matches_in_matrix <- function(rn,test_vector,test_matrix)	{
return(sum(test_vector==test_matrix[rn,]))
}

