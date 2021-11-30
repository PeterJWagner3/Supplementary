common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
source('~/Documents/R_Projects/Common_R_Source_Files/General_Plot_Templates.r');	#
source(paste(common_source_folder,"Nexus_File_Routines.r",sep=""));
source(paste(common_source_folder,"paleophylogeny_routines.r",sep=""));

ancestral_a <- c(0,-1);
ancestral_b <- c(1,0);
ancestral_c <- c(1,1);

Q2 <- construct_Q_matrix_unordered(2);
Qinap <- array(1,dim=c(3,3));
Qinap[1,2:3] <- 1;
nn <- rowSums(Qinap)-1;
Qinap[1,] <- 0.5;
for (i in 1:3)	Qinap[i,i] <- -(sum(Qinap[i,])-Qinap[i,i]);
#for (i in 1:3)	Qinap[i,] <- Qinap[i,]/abs(Qinap[i,i]);

alpha <- 0.05;
tmat_2 <- expm::expm(alpha*Q2);
tmat_inap <- expm::expm(alpha*Qinap);

ttl_reps <- 1000000;
changes1 <- rpois(ttl_reps,alpha);

ttl_change_counts <- hist(changes1,breaks=-1:max(changes1))$counts;
max_changes <- max(changes1);
names(ttl_change_counts) <- 0:max_changes;
combo_cases <- vector(length=3);
even_changes <- seq(1,max_changes+1,by=2);
odd_changes <- seq(2,max_changes+1,by=2);
names(combo_cases) <- c("0-","10","11");
combo_cases[1] <- sum(ttl_change_counts[even_changes]);
poss_states <- c(2,3);
nposs_states <- length(poss_states);
# case where all changes from tailless to tail generate a new tail color. Thus, a 3rd change can generate blue or red regardless of what the first change generated.
for (j in 1:length(odd_changes))	{
	steps <- odd_changes[j]-1;
	for (i in 1:ttl_change_counts[odd_changes[j]])	{
		dependent_state <- min(poss_states)+floor(nposs_states*runif(1));	# 50:50 that we get "red" or "blue"; 33:33:33 we get "red","blue" or "green"
		change_time <- runif(steps);	# when does it change? That effects how long it has to then change color from the original;
		dependent_changes <- rpois(1,alpha*min(change_time));
		if ((dependent_changes%%2)==1)	dependent_state <- 	poss_states[!poss_states %in% dependent_state];
		combo_cases[dependent_state] <- combo_cases[dependent_state]+1;
		}
	}
combo_cases <- combo_cases/ttl_reps;

#mxx=nrow(Qinap); mnx=0; mxy=1; mny=0; main="P Final Condition | 0- Ancestor";subtitle="New Color Each Derivation";abcissa=""; ordinate="Runs"; xsize=3; ysize=3; cexaxis=1; cexlab=1; cexmain=1; cexsub=1
tail_colors <- c("gray80","red","dodgerblue");
specify_basic_plot(mxx=nrow(Qinap), mnx=0, mxy=1, mny=0, main="P Final Condition | 0- Ancestor",subtitle="New Color Each Derivation",abcissa="", ordinate="Runs", xsize=3, ysize=3, cexaxis=1, cexlab=1, cexmain=1, cexsub=1)
specified_axis_w_labels(axe=1,max_val=nrow(Qinap),min_val=0,maj_break=1,med_break=0,min_break=0,axis_labels=names(combo_cases),linewd=4/3,label_pos="mid",font_size=1.25,orient=1,print_label=T)
specified_axis(axe=2,max_val=1,min_val=0,maj_break=0.10,med_break=0.05,min_break=0.025,orient=2);
for (i in 1:length(combo_cases))	rect(i-1,0,i,combo_cases[i],col=tail_colors[i]);
segments(0,tmat_inap[1,1],nrow(Qinap),tmat_inap[1,1],lty=2,lwd=2)
segments(0,tmat_inap[1,2],nrow(Qinap),tmat_inap[1,2],lty=3,lwd=1.5)

# case where all changes from tailless to tail generate a new tail color. Thus, a 3rd change can generate blue or red regardless of what the first change generated.
combo_cases_2 <- vector(length=3);
combo_cases_2[1] <- ttl_reps*combo_cases[1];
for (j in 1:length(odd_changes))	{
	steps <- odd_changes[j]-1;
	for (i in 1:ttl_change_counts[odd_changes[j]])	{
		dependent_state <- min(poss_states)+floor(nposs_states*runif(1));	# 50:50 that we get "red" or "blue"; 33:33:33 we get "red","blue" or "green"
		change_time <- runif(steps);	# when does it change? That effects how long it has to then change color from the original;
		dependent_changes <- rpois(1,alpha*max(change_time));
		if ((dependent_changes%%2)==1)	dependent_state <- 	poss_states[!poss_states %in% dependent_state];
		combo_cases_2[dependent_state] <- combo_cases_2[dependent_state]+1;
		}
	}
combo_cases_2 <- combo_cases_2/ttl_reps;

specify_basic_plot(mxx=nrow(Qinap), mnx=0, mxy=1, mny=0, main="P Final Condition | 0- Ancestor",subtitle="2ndary Tail Retains Ancestral Condition",abcissa="", ordinate="Runs", xsize=3, ysize=3, cexaxis=1, cexlab=1, cexmain=1, cexsub=1)
specified_axis_w_labels(axe=1,max_val=nrow(Qinap),min_val=0,maj_break=1,med_break=0,min_break=0,axis_labels=names(combo_cases),linewd=4/3,label_pos="mid",font_size=1.25,orient=1,print_label=T)
specified_axis(axe=2,max_val=1,min_val=0,maj_break=0.10,med_break=0.05,min_break=0.025,orient=2);
for (i in 1:length(combo_cases))	rect(i-1,0,i,combo_cases_2[i],col=tail_colors[i]);
segments(0,tmat_inap[1,1],nrow(Qinap),tmat_inap[1,1],lty=2,lwd=2)
segments(0,tmat_inap[1,2],nrow(Qinap),tmat_inap[1,2],lty=3,lwd=1.5)

# start with red tail: inherit ancestral tail;
combo_cases_3 <- vector(length=3);
combo_cases_3[1] <- sum(ttl_change_counts[odd_changes]);
even_changes_2 <- seq(2,max_changes,by=2);
changes2 <- rpois(ttl_change_counts[1],alpha);
ttl_changes_2ndary <- hist(changes2,breaks=-1:max(changes2),plot=F)$counts;
combo_cases_3[2] <- sum(ttl_changes_2ndary[even_changes]);
combo_cases_3[3] <- sum(ttl_changes_2ndary[odd_changes]);
combo_cases_3 <- combo_cases_3/ttl_reps;

specify_basic_plot(mxx=nrow(Qinap), mnx=0, mxy=1, mny=0, main="P Final Condition | 10 Ancestor",subtitle="2ndary Tail Retains Ancestral Condition",abcissa="", ordinate="Runs", xsize=3, ysize=3, cexaxis=1, cexlab=1, cexmain=1, cexsub=1)
specified_axis_w_labels(axe=1,max_val=nrow(Qinap),min_val=0,maj_break=1,med_break=0,min_break=0,axis_labels=names(combo_cases),linewd=4/3,label_pos="mid",font_size=1.25,orient=1,print_label=T)
specified_axis(axe=2,max_val=1,min_val=0,maj_break=0.10,med_break=0.05,min_break=0.025,orient=2);
for (i in 1:length(combo_cases_3))	rect(i-1,0,i,combo_cases_3[i],col=tail_colors[i]);
segments(0,tmat_inap[1,1],nrow(Qinap),tmat_inap[1,1],lty=2,lwd=2)
segments(0,tmat_inap[1,2],nrow(Qinap),tmat_inap[1,2],lty=3,lwd=1.5)


for (j in 1:length(odd_changes))	{
	steps <- odd_changes[j]-1;
	for (i in 1:ttl_change_counts[odd_changes[j]])	{
		dependent_state <- min(poss_states)+floor(nposs_states*runif(1));	# 50:50 that we get "red" or "blue"; 33:33:33 we get "red","blue" or "green"
		change_time <- runif(steps);	# when does it change? That effects how long it has to then change color from the original;
		dependent_changes <- rpois(1,alpha*max(change_time));
		if ((dependent_changes%%2)==1)	dependent_state <- 	poss_states[!poss_states %in% dependent_state];
		combo_cases_2[dependent_state] <- combo_cases_2[dependent_state]+1;
		}
	}


right_changes_ok <- changes[changes[,1]>0,];
right_changes_ok2 <- right_changes_ok[right_changes_ok[,2]!=1,];

{}