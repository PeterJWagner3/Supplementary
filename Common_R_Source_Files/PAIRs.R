#### Standard error ####

st_er <- function(x){ sd(x) / sqrt(length(x))}


#### Co-occurence functions (Toth et al. 2019) ####

simpairs <- function(x){ # FETmP function
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  
  occs = array()
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #FETmP Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      if(z[i,j]>=1) {z[i,j] <- 0.99999999999999994}
      z[i,j] = qnorm(z[i,j])
      z[j,i] = z[i,j]
    }
  }
  return(as.dist(z, diag = F, upper = F))
} 


simpairs.score <- function(x){ # FETmP function to return original co-occurrence score (prior to Z-transformation)
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  
  occs = array()
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #FETmP Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      if(z[i,j]>=1) {z[i,j] <- 0.99999999999999994}
      z[j,i] = z[i,j]
    }
  }
  return(as.dist(z, diag = F, upper = F))
} 

simpairspd <- function(A, B, nsites){
  cscore <- rep(0, min(A, B)+1)
  for(a in max(0, A+B-nsites):min(A,B)){
    for (k in max(0, A+B-nsites):a){
      cscore[a+1] = 
        choose(B , k) * choose(nsites - B , A - k) / choose(nsites , A)}
  }
  return(cscore)
}

resamp <- function(PA, sites = 20, reps = 100, replace = F){
  samp <- list()
  for(j in 1:reps){
    samp[[j]] <- PA[, sample(1:ncol(PA), sites, replace = replace)]
    samp[[j]] <- samp[[j]][which(rowSums(samp[[j]]) > 0),]
  }
  return(samp)}

resamp.MC <- function(PA,sites,reps,replace = F){
  samp <- list()
  for(j in 1:reps){
    samp[[j]] <- PA[, sample(1:ncol(PA), sites, replace = replace)]
    samp[[j]] <- samp[[j]][which(rowSums(samp[[j]]) > 0),]
  }
  return(samp)}

#### Redefine environmental categories ####

# Convert PBDB env. categories to shallow deep, corresponding to coast->shelf vs. slope/basinal
Redefine.Environ = function(x) { x %>% mutate(environ = case_when(
  environ == "deltaic indet." ~ "Shallow",
  environ == "intertidal" ~ "Shallow",
  environ == "deltaic indet" ~ "Shallow",
  environ == "peritidal" ~ "Shallow",
  environ == "reef, buildup or bioherm" ~ "Shallow",
  environ == "perireef or subreef" ~ "Shallow",
  environ == "shoreface" ~ "Shallow",
  environ == "delta front" ~ "Shallow",
  environ == "sand shoal" ~ "Shallow",
  environ == "foreshore" ~ "Shallow",
  environ == "coastal indet." ~ "Shallow",
  environ == "Shallow" ~ "Shallow",
  environ == "Delta" ~ "Shallow",
  environ == "lagoonal/restricted shallow subtidal" ~ "Shallow",
  environ == "shallow subtidal indet." ~ "Shallow",
  environ == "shallow subtidal indet" ~ "Shallow",
  environ == "open shallow subtidal"  ~ "Shallow",
  environ == "intrashelf/intraplatform reef" ~ "Shallow",
  environ == "transition zone/lower shoreface" ~ "Shallow",
  environ == "deep subtidal shelf" ~ "Shallow",
  environ == "deep subtidal indet." ~ "Shallow",
  environ == "deep subtidal indet" ~ "Shallow",
  environ == "deep subtidal ramp" ~ "Shallow",
  environ == "marginal marine indet." ~ "Shallow",
  environ == "platform/shelf-margin reef" ~ "Shallow",
  environ == "offshore shelf" ~ "Shallow",
  environ == "offshore" ~ "Shallow",
  environ == "offshore indet."  ~ "Shallow",
  
  environ == "prodelta" ~ "Deep",
  environ == "submarine fan" ~ "Deep",
  environ == "basinal (carbonate)" ~ "Deep",
  environ == "basinal (siliciclastic)" ~ "Deep",
  environ == "basinal (siliceous)" ~ "Deep",
  environ == "deep-water indet." ~ "Deep",
  environ == "slope" ~ "Deep",
  environ == "Deep" ~ "Deep",
  environ == "carbonate indet." ~  "",
  environ == "marine indet." ~ "",
  environ == "glacial" ~ "",
  environ == "terrestrial indet." ~ "Terrestrial"
))}

Redefine.Environ2 = function(x) {
  for(l in 1:nrow(x)){
  if(x[l,"environ"]== "deltaic indet."){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "intertidal"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "deltaic indet"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "peritidal"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "reef, buildup or bioherm"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "perireef or subreef"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "shoreface"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "delta front"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "sand shoal"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "foreshore"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "coastal indet."){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "coastal indet"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "Shallow"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "Delta"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "delta"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "lagoonal/restricted shallow subtidal"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "shallow subtidal indet."){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "shallow subtidal indet"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "open shallow subtidal" ){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "intrashelf/intraplatform reef"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "transition zone/lower shoreface"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "deep subtidal shelf"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "deep subtidal indet."){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "deep subtidal indet"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "deep subtidal ramp"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "marginal marine indet."){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "marginal marine indet"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "platform/shelf-margin reef"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "offshore shelf"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "offshore"){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "offshore indet." ){x[l,"environ"]<-"Shallow"}
  if(x[l,"environ"]== "offshore indet" ){x[l,"environ"]<-"Shallow"}
  
  if(x[l,"environ"]== "prodelta"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "submarine fan"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "basinal (carbonate)"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "basinal (siliciclastic)"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "basinal (siliceous)"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "deep-water indet."){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "slope"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "Deep"){x[l,"environ"]<-"Deep"}
  if(x[l,"environ"]== "carbonate indet."){x[l,"environ"]<- "Missing"}
  if(x[l,"environ"]== "carbonate indet"){x[l,"environ"]<- "Missing"}
  if(x[l,"environ"]== "marine indet."){x[l,"environ"]<-"Missing"}
  if(x[l,"environ"]== "marine indet"){x[l,"environ"]<-"Missing"}
  if(x[l,"environ"]== "glacial"){x[l,"environ"]<-"Missing"}
  if(x[l,"environ"]== "terrestrial indet."){x[l,"environ"]<-"Terrestrial"}
  if(x[l,"environ"]== "terrestrial indet"){x[l,"environ"]<-"Terrestrial"}
  }}




#### Custom cooccurence summary functions ####

# Function for the proportion of aggregated pairs over significant pairs
Prop.Agg = function(x) {x$positive/(x$positive+x$negative)}
# Function for the proportion of significant pairs over random pairs
Prop.Sig = function(x) {x$percent_sig/100}

# Generate summary table for resampled cooccurence analysis
Cooccur.Subsample=function(Data,Interval.Name){
  Data.Sub=resamp(Data)
  Cooccur.summary=data.frame()
  for(i in 1:100){
    cooc.run=cooccur(mat=Data.Sub[[i]],type="spp_site",spp_names=TRUE,prob="comb",thresh=FALSE)
    Cooccur.summary[1,i]=cooc.run$positive
    Cooccur.summary[2,i]=cooc.run$negative
    Cooccur.summary[3,i]=cooc.run$random
    Cooccur.summary[4,i]=cooc.run$unclassifiable
    Cooccur.summary[5,i]=sum(cooc.run$random,cooc.run$positive,cooc.run$negative,cooc.run$unclassifiable)
  }
  rownames(Cooccur.summary)=c("Positive","Negative","Random","Unclassifiable","Total Pairs")
  write.csv(Cooccur.summary, file=paste("Results/",Interval.Name, " - Subsample Summary.csv", sep=""), row.names=T)
}
#### Modified cooccurence matrix plot function ####

# plot custom plot
plot.cooccur.custom <-
  function(x, ...){
    
    ##
    allargs <- match.call(expand.dots = TRUE)
    plotrand <- allargs$plotrand
    plotrand <- ifelse(test = is.null(plotrand),yes = FALSE,no = plotrand)
    randsummary<- allargs$randsummary
    randsummary <- ifelse(test = is.null(randsummary),yes = FALSE,no = randsummary)
    
    ##
    
    dim <- x$species
    comat_pos <- comat_neg <- matrix(nrow=dim,ncol=dim)
    
    co_tab <- x$result
    for (i in 1:nrow(co_tab)){
      comat_pos[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_gt"]
      comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_gt"]
      
      row.names(comat_pos[co_tab[i,"sp2"],co_tab[i,"sp1"]])
      
    }
    for (i in 1:nrow(co_tab)){
      comat_neg[co_tab[i,"sp1"],co_tab[i,"sp2"]] <- co_tab[i,"p_lt"]
      comat_neg[co_tab[i,"sp2"],co_tab[i,"sp1"]] <- co_tab[i,"p_lt"]
    }
    comat <- ifelse(comat_pos>=0.05,0,1) + ifelse(comat_neg>=0.05,0,-1)
    colnames(comat) <- 1:dim
    row.names(comat) <- 1:dim
    
    if ("spp_key" %in% names(x)){
      
      sp1_name <- merge(x=data.frame(order=1:length(colnames(comat)),sp1=colnames(comat)),y=x$spp_key,by.x="sp1",by.y="num",all.x=T)
      sp2_name <- merge(x=data.frame(order=1:length(row.names(comat)),sp2=row.names(comat)),y=x$spp_key,by.x="sp2",by.y="num",all.x=T)
      
      colnames(comat) <- sp1_name[with(sp1_name,order(order)),"spp"]  
      row.names(comat) <- sp2_name[with(sp2_name,order(order)),"spp"]
      
    }  
    
    #ind <- apply(comat, 1, function(x) all(is.na(x)))
    #comat <- comat[!ind,]
    #ind <- apply(comat, 2, function(x) all(is.na(x)))
    #comat <- comat[,!ind]
    
    comat[is.na(comat)] <- 0
    
    origN <- nrow(comat)
    
    # SECTION TO REMOVE SPECIES INTERACTION WITH NO OTHERS
    
    #rmrandomspp <- function(orimat,plotrand = FALSE,randsummary = FALSE){
    if(plotrand == FALSE){
      ind <- apply(comat, 1, function(x) all(x==0))
      comat <- comat[!ind,]    
      ind <- apply(comat, 2, function(x) all(x==0))
      comat <- comat[,!ind]
      #ind <- apply(orimat, 1, function(x) all(x==0))
      #orimat <- orimat[!ind,]    
      #ind <- apply(orimat, 2, function(x) all(x==0))
      #orimat <- orimat[,!ind]
    }
    #return(orimat)
    #}
    
    #comat <- rmrandomspp(orimat = comat, dots)
    
    postN <- nrow(comat)
    
    
    comat <- comat[order(rowSums(comat)),]
    comat <- comat[,order(colSums(comat))]
    
    #comat <- rmrandomspp(orimat = comat, ...)
    
    #ind <- apply(comat, 1, function(x) all(x==0))
    #comat <- comat[!ind,]
    #ind <- apply(comat, 2, function(x) all(x==0))
    #comat <- comat[,!ind]
    
    ind <- apply(comat, 1, function(x) all(x==0))
    comat <- comat[names(sort(ind)),]
    ind <- apply(comat, 2, function(x) all(x==0))
    comat <- comat[,names(sort(ind))]
    
    #comat
    data.m = melt(comat)
    colnames(data.m) <- c("X1","X2","value")
    data.m$X1 <- as.character(data.m$X1)
    data.m$X2 <- as.character(data.m$X2)
    
    meas <- as.character(unique(data.m$X2))
    
    dfids <- subset(data.m, X1 == X2)
    
    X1 <- data.m$X1
    X2 <- data.m$X2
    
    df.lower = subset(data.m[lower.tri(comat),],X1 != X2)
    
    ##### testing the rand summary
    if(randsummary == FALSE){  
    }else{
      dim <- nrow(comat)
      ext.dim <- round(dim*0.2,digits = 0)
      if(ext.dim<0){ext.dim<-1}
      placehold <- paste("ext_", rep(c(1:ext.dim),each = dim), sep="")
      
      randcol.df <- data.frame(
        X1 = placehold,
        X2 = rep(meas,times = ext.dim),
        value = rep(x = c(-2), times = dim*ext.dim))
      
      df.lower <- rbind(df.lower,randcol.df)
      meas <- c(meas,unique(placehold))
    }
    
    X1 <- df.lower$X1
    X2 <- df.lower$X2
    value <- df.lower$value

    if(randsummary == FALSE){  
      p <- ggplot(df.lower, aes(X1, X2)) + 
        geom_raster(aes(fill = factor(value,levels=c(-1,0,1)))) +
        scale_fill_manual(values = c("red","dark gray","blue"), name = "", labels = c("negative","random","positive"),drop=FALSE) + 
        theme(axis.text.y.left=element_blank(),axis.text.x=element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=15, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.2, 0.8),legend.text=element_text(size=12)) + 
        xlab("<< Higher Abundance Genus") + ylab("<< Higher Abundance Genus") + 
        scale_x_discrete(limits=meas, expand = c(0.05, 0),drop=FALSE) + 
        scale_y_discrete(limits=meas, expand = c(0.05, 0),drop=FALSE) 
    }else{
  
      p <- ggplot(df.lower, mapping=aes(X1, X2)) + geom_tile(aes(fill = factor(value,levels=c(-1,0,1,-2))), colour ="white") 
      p <- p + scale_fill_manual(values = c("#FFCC66","dark gray","light blue","light gray"), name = "", labels = c("negative","random","positive","random"),drop=FALSE) + 
        theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),plot.title = element_text(vjust=-4,size=20, face="bold"),panel.background = element_rect(fill='white', colour='white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9, 0.5),legend.text=element_text(size=18)) + 
        ggtitle("Species Co-occurrence Matrix") + 
        xlab("") + ylab("") + 
        scale_x_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) + 
        scale_y_discrete(limits=meas, expand = c(0.3, 0),drop=FALSE) 
      p <- p + geom_text(data=dfids,aes(label=X1),hjust=1,vjust=0,angle = -22.5)#, color="dark gray")
      
      dim <- nrow(comat)
      ext_x <- dim + 0.5 #(ext.dim/2)
      ext_y <- dim + 1
      nrem <- origN - postN
      randtext <- paste(nrem, " completely\nrandom species")
      ext_dat <- data.frame(ext_x=ext_x,ext_y=ext_y,randtext=randtext)
      
      p <- p + geom_text(data=ext_dat,aes(x = ext_x,y = ext_y,label=randtext),hjust=0,vjust=0, color="dark gray")
    }
    ####
    
    p
    
  }

#### Blois 2014 ####

# classify site type (01,10,00,11)
classifySites<- function(sppDat, spp1, spp2){
  
  # identify the row number of each species in the original data matrix  
  row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
  row2<- which(rownames(sppDat)==spp2)
  
  # Determine the four classes of sites (10, 01, 11, 00) for each pair
  siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
  
  # which sites have spp1 but not spp2? (10)
  sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
  siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
  
  # which sites have spp2 but not spp1? (01)
  sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
  siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
  
  # which sites have both spp1 and spp2? (11)
  sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
  siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
  
  # which sites have neither spp1 or spp2? (00)
  sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
  siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
  
  return(siteType)  
}

# manova
performManova<- function(manovaDat){
  
  Y<- cbind(manovaDat[,2], manovaDat[,3])
  X<- as.factor(manovaDat[,1])
  maovResults<- manova(Y ~ X)
  
  return(maovResults)
}


#### Analysis of abiotic factors (Blois 2014) ####
Blois.Factors<-function(time,pairsDat,sppDat,siteDat,randomPairs,segPairs,aggPairs){
  
  # path to the output files
  outputFile<- paste("Results/",time,"/",sep="")
  
  # prepare data files
  pairsDat <- pairsDat %>% select(sp1_name, sp2_name,sp1,sp2,prob_cooccur,exp_cooccur) 
  
  # ARRANGE EACH FILE
  
  #check site names and order
  all(rownames(siteDat)==colnames(sppDat))
  
  # SITES
  # retain sites from appropriate time 
  siteDatOrig<- siteDat
  siteDat<- siteDat[match(colnames(sppDat), rownames(siteDat)),]
  
  #check site names and order
  all(rownames(siteDat)==colnames(sppDat))
  
  # Set the number of variables we are testing for differences among sites
  numVar <- ncol(siteDat) #how many variables are we testing for?
  
  # Check for no aggregated or segregated pairs
  if (sum(aggPairs$prob_cooccur)==0 & sum(segPairs$prob_cooccur)>0){
    #### Blois: Segregated - Random loop ####
    # set the three different types of analyses
    analyses<- c("random", "segregated")
    
    # create summary table of analyses
    finalSummary<- matrix(nrow=2, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    # Analyse each pair type separately
    for (k in 1:2){  
      #set the analysis type
      if (k==1){dat<- randomPairs}
      if (k==2){dat<- segPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        num.run=as.numeric(i)
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[-NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar} #included for completeness.
        if (k ==2){siteChar<- siteChar[which((siteChar$siteType=="10") | (siteChar$siteType=="01")),]} # seg
        
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique(manovaDat$siteType))<2) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            manovaSignificance[i,]<- round(summary(maovResults)$stats[1,6],4)
          } #DS
       
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                  
                  cat("### CHI SQUARE results", fill=T)
                  chiResults<- chisq.test(subDat[,1], subDat[,2])
                  
                  print(chiResults)    
                  cat(" ", fill=T)      
                  
                  varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1, analysis2)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(randomPairs)/nrow(pairsDat),
                  nrow(aggPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)
    
    }else{
  if (sum(aggPairs$prob_cooccur)==0 & sum(segPairs$prob_cooccur)==0){
    #### Blois: Random loop ####
    # set the three different types of analyses
    analyses<- c("random")

    # create summary table of analyses
    finalSummary<- matrix(nrow=1, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    # Analyse each pair type separately
    for (k in 1){  
      
      #set the analysis type
      if (k==1){dat<- randomPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        num.run=as.numeric(i)
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[-NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar} #included for completeness.
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique(manovaDat$siteType))<2) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            manovaSignificance[i,]<- round(summary(maovResults)$stats[1,6],4)
          } #DS
          
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                  
                  cat("### CHI SQUARE results", fill=T)
                  chiResults<- chisq.test(subDat[,1], subDat[,2])
                  
                  print(chiResults)    
                  cat(" ", fill=T)      
                  
                  varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(randomPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)
    
  }else{if (sum(segPairs$prob_cooccur)==0 & sum(aggPairs$prob_cooccur)>0){
    #### Blois: Aggregated - Random loop ####
    # set the three different types of analyses
    analyses<- c("random", "aggregated")

    # create summary table of analyses
    finalSummary<- matrix(nrow=2, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    # Analyse each pair type separately
    for (k in 1:2){  

      #set the analysis type
      if (k==1){dat<- randomPairs}
      if (k==2){dat<- aggPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        num.run=as.numeric(i)
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[-NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar} #included for completeness.
        if (k ==2){siteChar<- siteChar[which((siteChar$siteType=="00") | (siteChar$siteType=="11")),]} # agg
        
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique(manovaDat$siteType))<2) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            manovaSignificance[i,]<- round(summary(maovResults)$stats[1,6],4)
          } #DS
          
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA
                  }else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                  
                  cat("### CHI SQUARE results", fill=T)
                  chiResults<- chisq.test(subDat[,1], subDat[,2])
                  
                  print(chiResults)    
                  cat(" ", fill=T)      
                  varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1, analysis2)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(randomPairs)/nrow(pairsDat),
                  nrow(aggPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)
  }else{
    #### Blois: Aggregated - Segregated - Random loop ####
    # set the three different types of analyses
    analyses<- c("random", "segregated", "aggregated")

    # create summary table of analyses
    finalSummary<- matrix(nrow=3, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    # Analyse each pair type separately
    for (k in 1:3){  
      
      #set the analysis type
      if (k==1){dat<- randomPairs}
      if (k==2){dat<- segPairs}
      if (k==3){dat<- aggPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[nas,]
        # eliminatedSites<- siteChar[-nas,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar} #included for completeness. Random
        if (k ==2){siteChar<- siteChar[which((siteChar$siteType=="10") | (siteChar$siteType=="01")),]} # seg
        if (k ==3){siteChar<- siteChar[which((siteChar$siteType=="00") | (siteChar$siteType=="11")),]} # ag
        
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique(manovaDat$siteType))<2) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            manovaSignificance[i,]<- round(summary(maovResults)$stats[1,6],4)
          } #DS
          
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                  cat("### CHI SQUARE results", fill=T)
                  chiResults<- chisq.test(subDat[,1], subDat[,2])
                  
                  print(chiResults)    
                  cat(" ", fill=T)      
                  varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1, analysis2, analysis3)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(randomPairs)/nrow(pairsDat),
                  nrow(segPairs)/nrow(pairsDat),
                  nrow(aggPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)
  }}} # else 2: random only and segregated only
} # closes function


#### Blois 2014 modified for FETmP results ####
Blois.FETmP <- function(time,sppDat,siteDat,Pairs.Results){
  # Process data
  
  Cooccur.TS<-cooccur(mat=sppDat,type="spp_site",spp_names=TRUE,prob="comb",thresh=FALSE)
  pairsDat <- as.data.frame(Cooccur.TS[["results"]])
  aggPairs <- Pairs.Results %>% filter(Zone==time & Pair.Type=="Agg") %>% select(sp1_name,sp2_name,Mean.FETmP);colnames(aggPairs)<-c("sp1_name","sp2_name","prob_cooccur")
  segPairs <- Pairs.Results %>% filter(Zone==time & Pair.Type=="Seg") %>% select(sp1_name,sp2_name,Mean.FETmP);colnames(segPairs)<-c("sp1_name","sp2_name","prob_cooccur")
  randPairs <- Pairs.Results %>% filter(Zone==time & Pair.Type=="Random") %>% select(sp1_name,sp2_name,Mean.FETmP);colnames(randPairs )<-c("sp1","sp2","prob_cooccur")
  
  # path to the output files
  outputFile<- paste("Results/",time,"/",sep="")
  
  # prepare data files
  pairsDat <- pairsDat %>% select(sp1_name, sp2_name,sp1,sp2,prob_cooccur,exp_cooccur) 
  
  # ARRANGE EACH FILE
  
  #check site names and order
  all(rownames(siteDat)==colnames(sppDat))
  
  # SITES
  # retain sites from appropriate time 
  siteDatOrig<- siteDat
  siteDat<- siteDat[match(colnames(sppDat), rownames(siteDat)),]
  
  #check site names and order
  all(rownames(siteDat)==colnames(sppDat))
  
  # Set the number of variables we are testing for differences among sites
  numVar <- ncol(siteDat) #how many variables are we testing for?
  
  #### Analysis variant: Agg pairs ####
  if (sum(segPairs$prob_cooccur)==0 & sum(aggPairs$prob_cooccur)>0){
    
    analyses<- c("aggregated")
    
    # create summary table of analyses
    finalSummary<- matrix(nrow=1, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    for (k in 1){
      #set the analysis type
      if (k==1){dat<- aggPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar[which((siteChar$siteType=="00") | (siteChar$siteType=="11")),]} # ag
        
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          
          if(length(unique((manovaDat$siteType)))<2 | nrow(manovaDat)<3 ) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            X1<-Y[,1] # MC
            X2<-Y[,2] # MC
            Y1<-as.numeric(X) # MC
            corr.P<-lm(Y1 ~ X1*X2)
            if(maovResults$df.residual<=3|summary(corr.P)$adj.r.squared==1){manovaSignificance[i,]<-NA # Skip MANOVA with insufficient DF or assumption violation
            #round(summary(corr.P)$coefficients[,"Pr(>|t|)"][2],4)# Perform a MLR when df is insufficient for MANOVA
            }else{manovaSignificance[i,]<-round(summary(maovResults)$stats[1,6],4)}
            if(is.na(manovaSignificance[i,])==TRUE){manovaSignificance[i,]<-1} # MC: set p value to an NA ID and interpret as non-significant for skipped MANOVAs
          } #DS
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    cat("### CHI SQUARE results", fill=T)
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    print(chiResults)    
                    cat(" ", fill=T)      
                    varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                  }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(aggPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)}# end agg analysis
  
  #### Analysis variant: Agg+Seg pairs ####
  if(sum(segPairs$prob_cooccur)<0 & sum(aggPairs$prob_cooccur)>0){
    analyses<- c("segregated","aggregated")
    
    # create summary table of analyses
    finalSummary<- matrix(nrow=2, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    for (k in 1:2){ 
      #set the analysis type
      if (k==1){dat<- segPairs}
      if (k==2){dat<- aggPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar[which((siteChar$siteType=="10") | (siteChar$siteType=="01")),]} # seg
        if (k ==2){siteChar<- siteChar[which((siteChar$siteType=="00") | (siteChar$siteType=="11")),]} # ag
        
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique((manovaDat$siteType)))<2 | nrow(manovaDat)<3 ) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            X1<-Y[,1] # MC
            X2<-Y[,2] # MC
            Y1<-as.numeric(X) # MC
            corr.P<-lm(Y1 ~ X1*X2)
            if(maovResults$df.residual<=3|summary(corr.P)$adj.r.squared==1){manovaSignificance[i,]<-NA # Skip MANOVA with insufficient DF or assumption violation
            #round(summary(corr.P)$coefficients[,"Pr(>|t|)"][2],4)# Perform a MLR when df is insufficient for MANOVA
            }else{manovaSignificance[i,]<-round(summary(maovResults)$stats[1,6],4)}
            if(is.na(manovaSignificance[i,])==TRUE){manovaSignificance[i,]<-1} # MC: set p value to an NA ID and interpret as non-significant for skipped MANOVAs
          } #DS
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    cat("### CHI SQUARE results", fill=T)
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    print(chiResults)    
                    cat(" ", fill=T)      
                    varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                  }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1, analysis2)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(segPairs)/nrow(pairsDat),
                  nrow(aggPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)}# end analysis variant
  
  #### Analysis variant: Seg pairs ####
  if(sum(segPairs$prob_cooccur)<0 & sum(aggPairs$prob_cooccur)==0){
    analyses<- c("segregated")
    
    # create summary table of analyses
    finalSummary<- matrix(nrow=1, ncol=(numVar+1)*2)
    rownames(finalSummary)<- analyses
    
    for (k in 1){ 
      #set the analysis type
      if (k==1){dat<- segPairs}
      
      # write the outputs to an ANOVA file for each type of analysis
      sink(file=paste(outputFile, "ANOVAResults-", analyses[k], ".txt", sep="")) 
      
      cat("######### Analysis Type = ", analyses[k], " #########", fill=T) #print the analysis to file
      
      # Set result matrices
      sumDat<- matrix(ncol=4, nrow=nrow(dat)) #matrix of 4 columns: number of sites with 10, 01, 11, 00 
      varSignificance<- matrix(nrow=nrow(dat), ncol=numVar) #ANOVA or Chi-squared p-values for each site variable
      manovaSignificance<- matrix(nrow=nrow(dat), ncol=1) #MANOVA on latitude and longitude
      
      # SCROLL THROUGH EACH PAIR AND RUN THE SITE-CHARACTERISTICS ANALYSIS
      for (i in 1:nrow(dat)){
        
        # identify the species names in the pairs file
        spp1<- as.character(dat$sp1_name[i])  
        spp2<- as.character(dat$sp2_name[i])
        cat("######### Species Pair = ", spp1, " and ", spp2, " #########", fill=T)  #write the species pair to the output file
        
        # identify the row number of each species in the original data matrix  
        row1<- which(rownames(sppDat)==spp1)  #pairs output truncates the spp names to 18 char
        row2<- which(rownames(sppDat)==spp2)
        
        # Determine the four classes of sites (10, 01, 11, 00) for each pair
        ###  rename this siteType (to separate from pairType)
        siteType<- vector(length=ncol(sppDat)) # create a new vector to store the site type
        
        # which sites have spp1 but not spp2? (10)
        sites10<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,]==0))]<- "10"
        
        # which sites have spp2 but not spp1? (01)
        sites01<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,] > 0))]<- "01"
        
        # which sites have both spp1 and spp2? (11)
        sites11<- colnames(sppDat)[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]
        siteType[ which((sppDat[row1,] > 0) & (sppDat[row2,] > 0))]<- "11"
        
        # which sites have neither spp1 or spp2? (00)
        sites00<- colnames(sppDat)[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]
        siteType[ which((sppDat[row1,]==0) & (sppDat[row2,]==0))]<- "00"
        
        # Create a new site characteristics file based on siteDat, by removing site name column and adding siteType
        ## siteChar is a matrix with each site in a row, with columns indicating the site type (11, 10, 01, 00) and the characteristics of the site (precip, temp, etc)
        siteChar<- cbind(as.character(siteType), siteDat) 
        colnames(siteChar)[1]<- "siteType"
        
        # Do any of the sites have missing values?  If so, remove them.  #Note 6/17.  Missing values should be ok
        # siteChar<- siteChar[NA,]
        # eliminatedSites<- siteChar[NA,] 
        
        # fill in the sumDat results file
        # sumDat shows the number of sites in each category for each pair
        sumDat[i,]<- c(length(which(siteChar$siteType=="10")), 
                       length(which(siteChar$siteType=="01")), 
                       length(which(siteChar$siteType=="11")), 
                       length(which(siteChar$siteType=="00")))
        
        # trim siteChar depending on the type of analysis:
        # 1.  For "random pairs" analysis, compare all 4 site types
        # 2.  For "segregated pairs" analysis, compare sites (01) & (10)
        # 3.  For "aggregated pairs" analysis, compare sites (00) & (11)
        
        if (k ==1){siteChar<- siteChar[which((siteChar$siteType=="10") | (siteChar$siteType=="01")),]} # seg
        
        # scroll through each site characteristic and determine if there are significant differences in the characteristics among the 4 types of sites
        if (length(unique(siteChar$siteType))==1){  #skip this step if there is only one type of site, after trimming for missing data and analysis type
          varSignificance[i,]<- NA
          cat("Only one type of site remaining", fill=T)
        }else{
          
          # dispersal test section
          manovaDat<- siteChar[, c(1, match(c("LONG", "LAT"), colnames(siteChar)))]
          
          ##DS correction: check to see if there are multiple siteTypes left, if not, skip
          manovaDat=na.omit(manovaDat) #DS
          if(length(unique((manovaDat$siteType)))<2 | nrow(manovaDat)<3 ) { #DS
            manovaSignificance[i,] = varSignificance[i,]= NA #DS
          }else{ #DS
            Y<- cbind(manovaDat[,2], manovaDat[,3])
            X<- factor(manovaDat[,1])
            maovResults<- manova(Y ~ X)
            X1<-Y[,1] # MC
            X2<-Y[,2] # MC
            Y1<-as.numeric(X) # MC
            corr.P<-lm(Y1 ~ X1*X2)
            if(maovResults$df.residual<=3|summary(corr.P)$adj.r.squared==1){manovaSignificance[i,]<-NA # Skip MANOVA with insufficient DF or assumption violation
            #round(summary(corr.P)$coefficients[,"Pr(>|t|)"][2],4)# Perform a MLR when df is insufficient for MANOVA
            }else{manovaSignificance[i,]<-round(summary(maovResults)$stats[1,6],4)}
            if(is.na(manovaSignificance[i,])==TRUE){manovaSignificance[i,]<-1} # MC: set p value to an NA ID and interpret as non-significant for skipped MANOVAs
          } #DS
          # anova or chi-square single variable analyses
          for (j in 1:numVar){   
            col<- j+1
            char<- colnames(siteChar)[col]  #name of the site characteristic
            subDat<- siteChar[,c(1, col)]
            
            ##DS remove NAs, if only one siteType left, skip
            subDat=na.omit(subDat)
            if(length(unique(subDat$siteType))<2){
              manovaSignificance[i,] = varSignificance[i,]= NA
            }else{
              #print summary stats  
              cat("##### Variable =", char, " #####", fill=T)  
              if ((k>1) && (char=="DRIFTLESS")) {
                varSignificance[i,(j)]<- NA
              }else{       
                if (is.numeric(siteChar[,col])){  #if the variable is numeric, do an anova
                  cat("## Mean for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) mean(na.omit(x)))) #mean value for each site type
                  
                  cat("## Standard deviation for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), function(x) sd(na.omit(x)))) #sd value for each site type
                  
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type
                  
                  #boxplot(get(char)~siteType, data=subDat, main=char)
                  
                  cat("### ANOVA results", fill=T)
                  ss<- tapply(subDat[,2], list(subDat$siteType), length)
                  ss<- as.data.frame(ss)
                  ss<- na.omit(ss)
                  
                  if ((dim(ss)[1]==2) && (min(ss)==1)){ #if only two categories and one of them only has one site, then skip
                    varSignificance[i,(j)]<- NA
                    
                  }else{
                    aovResults<- aov(get(char)~siteType, data=subDat)
                    
                    print(summary(aovResults))    
                    cat(" ", fill=T)      
                    
                    varSignificance[i,(j)]<- round(summary(aovResults)[[1]][5][1,],4)
                  }
                }else{  #if the variable is categorical, do a chi-square test
                  cat("## Count of each site and habitat type", fill=T)
                  print(table(subDat)) #count for each site type, broken up into habitat types
                  cat("## Sample size for each site type", fill=T)
                  print(tapply(subDat[,2], list(subDat$siteType), length)) #sample size of each site type       
                  
                  if(length(unique(subDat[,2]))<2){varSignificance[i,(j)]<-NA} else{
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    cat("### CHI SQUARE results", fill=T)
                    chiResults<- chisq.test(subDat[,1], subDat[,2])
                    
                    print(chiResults)    
                    cat(" ", fill=T)      
                    varSignificance[i,(j)]<- round(chiResults$ p.value,4)
                  }}
              }
            }
          } #close j    
        } #close ifelse
      } #close i
      # clean up the varSignificance and sumDat matrices
      colnames(varSignificance)<- colnames(siteChar)[-1]
      colnames(manovaSignificance)<- "distance"
      varSignificance<- cbind(manovaSignificance, varSignificance)
      colnames(sumDat)<- c("nSites10", "nSites01", "nSites11", "nSites00") 
      finalSummary[k,1:(ncol(varSignificance))]<- apply(varSignificance, 2, function(x) length(which(x<=0.05)))
      finalSummary[k,(ncol(varSignificance)+1):ncol(finalSummary)]<-apply(varSignificance, 2, function(x) length(which(x>0.05)))
      # bind all the relevant information together
      pairsType<- rep(analyses[k], nrow(varSignificance))
      datMod<- cbind(dat, pairsType, sumDat, varSignificance)
      
      # write csv file for posterity and store in different variable name for the future
      write.csv(datMod, paste(outputFile, "R_output-", analyses[k], "_pairs.csv", sep=""), row.names=F)
      assign(paste("analysis", k, sep=""), datMod)
      sink() #close the sink file for the analysis  
    } #close k
    
    #  bind the results of the three analyses together and save as csv
    finalAnalysis<- rbind(analysis1)
    write.csv(finalAnalysis, file=paste(outputFile, "R_output-allPairs.csv", sep=""), row.names=F)
    
    # tidy up finalSummary table and save as csv
    colnames(finalSummary)<- c(paste(colnames(varSignificance), "-sig", sep=""), paste(colnames(varSignificance), "-ns", sep=""))
    finalSummary<- finalSummary[,order(colnames(finalSummary))]
    
    # add frequence of each type of pair to the summary file
    Frequency<- c(nrow(segPairs)/nrow(pairsDat))
    finalSummary<- cbind(Frequency, finalSummary)
    
    write.csv(finalSummary, file=paste(outputFile, "R_output-finalSummary.csv", sep=""), row.names=T)}# end analysis variant
}

#### Mega Cooccurence Analysis ####

Cooccur.Mega=function(Genus.Dat,Env.Dat,Interval,sites,reps){
  Data.Sub=resamp.MC(Genus.Dat,sites,reps) # Subsample sites (20)
  Env.TS = Env.Dat
  Cooccur.Summary=data.frame() # Creates summary table for cooccurence
  Factor.Summary=data.frame() # Creates summary table for Blois anlaysis
  Agg.Total<-data.frame()
  Seg.Total<-data.frame()
  for(i in 1:reps){
    print(i)
    Data.TS = Data.Sub[[i]]
    f = c(colnames(Data.TS))
    Env.Sub = as.data.frame(t(subset(t(Env.TS),select=f)))
    Cooccur.TS = cooccur(mat=Data.TS,type="spp_site",spp_names=TRUE,prob="comb",thresh=FALSE)
    Genus.TS = as.data.frame(Cooccur.TS[["results"]])
    Genus.TS$OE <- abs(Genus.TS$obs_cooccur-Genus.TS$exp_cooccur)
    Genus.TS$p_gt[Genus.TS$OE<1] <-0.995
    Genus.TS$p_lt[Genus.TS$OE<1] <-0.995
    Agg.Genus.TS = Genus.TS %>% filter(p_gt<0.05) %>% dplyr::select(sp1_name,sp2_name,prob_cooccur)
    Seg.Genus.TS = Genus.TS %>% filter(p_lt<0.05) %>% dplyr::select(sp1_name,sp2_name,prob_cooccur)
    Un.Genus.TS = Genus.TS %>% filter(p_lt==0.995 & p_gt==0.995) %>% dplyr::select(sp1_name,sp2_name,prob_cooccur)
    Rand.Genus.TS = Genus.TS %>% filter(p_lt>0.05 & p_gt>0.05) %>% dplyr::select(sp1_name,sp2_name,prob_cooccur)
    Genus.TS.Pairs=pair.attributes(Cooccur.TS)
    
    # Create Summary Tables for co-ooccurence
    Cooccur.Summary[1,i]=nrow(Agg.Genus.TS)
    Cooccur.Summary[2,i]=nrow(Seg.Genus.TS)
    Cooccur.Summary[3,i]=(nrow(Rand.Genus.TS)-nrow(Un.Genus.TS))
    Cooccur.Summary[4,i]=nrow(Un.Genus.TS)
    Cooccur.Summary[5,i]=Cooccur.TS$pairs
    
    Agg.Total<-rbind(Agg.Total,Agg.Genus.TS)
    Seg.Total<-rbind(Seg.Total,Seg.Genus.TS)
    
    #Factor Analysis
    Blois.Factors(Interval,Genus.TS,Data.TS,Env.Sub,Rand.Genus.TS,Seg.Genus.TS,Agg.Genus.TS)
    
    # Summary Table
    finalAnalysis=as.data.frame(read.csv(paste("Results/",Interval,"/R_output-allPairs.csv", sep="")))
    finalAnalysis[is.na(finalAnalysis)]=1  # clear NA's from tabulation
    
    #Both
    F.Both.Agg=sum((finalAnalysis$pairsType=="aggregated")==TRUE &
                     (finalAnalysis$environ<0.05)==TRUE & 
                     (finalAnalysis$distance<0.05)==TRUE)
    F.Both.Seg=sum((finalAnalysis$pairsType=="segregated")==TRUE &
                     (finalAnalysis$environ<0.05)==TRUE & 
                     (finalAnalysis$distance<0.05)==TRUE)
    # Neither
    F.Neither.Agg=sum((finalAnalysis$pairsType=="aggregated")==TRUE &
                        (finalAnalysis$environ>0.05)==TRUE & 
                        (finalAnalysis$distance>0.05)==TRUE)
    F.Neithter.Seg=sum((finalAnalysis$pairsType=="segregated")==TRUE &
                         (finalAnalysis$environ>0.05)==TRUE & 
                         (finalAnalysis$distance>0.05)==TRUE)
    # Distance
    F.Dist.Agg=sum((finalAnalysis$pairsType=="aggregated")==TRUE &
                     (finalAnalysis$environ>0.05)==TRUE & 
                     (finalAnalysis$distance<0.05)==TRUE)
    F.Dist.Seg=sum((finalAnalysis$pairsType=="segregated")==TRUE &
                     (finalAnalysis$environ>0.05)==TRUE & 
                     (finalAnalysis$distance<0.05)==TRUE)
    # Environment
    F.Env.Agg=sum((finalAnalysis$pairsType=="aggregated")==TRUE &
                    (finalAnalysis$environ<0.05)==TRUE & 
                    (finalAnalysis$distance>0.05)==TRUE)
    F.Env.Seg=sum((finalAnalysis$pairsType=="segregated")==TRUE &
                    (finalAnalysis$environ<0.05)==TRUE & 
                    (finalAnalysis$distance>0.05)==TRUE)
    
    Factor.Summary[1,i]=as.numeric(F.Both.Agg)
    Factor.Summary[2,i]=as.numeric(F.Neither.Agg)
    Factor.Summary[3,i]=as.numeric(F.Dist.Agg)
    Factor.Summary[4,i]=as.numeric(F.Env.Agg)
    Factor.Summary[5,i]=as.numeric(F.Both.Seg)
    Factor.Summary[6,i]=as.numeric(F.Neithter.Seg)
    Factor.Summary[7,i]=as.numeric(F.Dist.Seg)
    Factor.Summary[8,i]=as.numeric(F.Env.Seg)
    Factor.Summary[8,i]=as.numeric(F.Env.Seg)
  }
  rownames(Cooccur.Summary)=c("Positive","Negative","Random","Unclassifiable","Total Pairs")
  write.csv(Cooccur.Summary, file=paste("Results/",Interval,"/100 Cooccurence Summary.csv", sep=""), row.names=T)
  
  rownames(Factor.Summary)=c("Both.Agg","Neither.Agg","Dist.Agg","Env.Agg",
                             "Both.Seg","Neither.Seg","Dist.Seg","Env.Seg")
  write.csv(Factor.Summary, file=paste("Results/",Interval,"/100 Factor Summary.csv", sep=""), row.names=T)
  
  if(nrow(Agg.Total)>0){
    Agg.Total.C <- Agg.Total %>% select(sp1_name,sp2_name)
    Agg.Total.C$Var1<-paste(Agg.Total.C$sp1_name,"-",Agg.Total.C$sp2_name,sep="")
    Agg.Total.Table<-as.data.frame(table(Agg.Total.C$Var1))
    Agg.Total.Table$Pair.Type<-"Aggregated"
    Agg.Total.Table=left_join(as.data.frame(Agg.Total.Table),as.data.frame(Agg.Total.C),by="Var1")
    Agg.Total.Table=unique(Agg.Total.Table[,1:ncol(Agg.Total.Table)])
  }else{Agg.Total.Table<-data.frame()}
  
  if(nrow(Seg.Total)>0){
    Seg.Total.C <- Seg.Total %>% select(sp1_name,sp2_name)
    Seg.Total.C$Var1<-paste(Seg.Total.C$sp1_name,"-",Seg.Total.C$sp2_name,sep="")
    Seg.Total.Table<-as.data.frame(table(Seg.Total.C$Var1))
    Seg.Total.Table$Pair.Type<-"Segregated"
    Seg.Total.Table=left_join(as.data.frame(Seg.Total.Table),Seg.Total.C,by="Var1")
    Seg.Total.Table=unique(Seg.Total.Table[,1:ncol(Seg.Total.Table)])
  }else{Seg.Total.Table<-data.frame()}

  Species.Pairs<-rbind(Agg.Total.Table,Seg.Total.Table)
  write.csv(Species.Pairs,file=paste("Results/",Interval,"/Species Summary ",Interval,".csv",sep=""))
}


#### Veech co-occurrence with specified subsampling ####

Cooc.Combined <- function(Genus.Dat,Interval,sites,reps){
  Data.Sub<-resamp.MC(Genus.Dat,sites,reps) # Subsample sites over N reps
  # FETmP-based Cooccurrence
  St.TS<-as.matrix(simpairs.score(Genus.Dat)) # Base FETmP
  colnames(St.TS)<-rownames(Genus.Dat)
  rownames(St.TS)<-rownames(Genus.Dat)
  St.Table.1 <- PairsTable(St.TS) # Create core of summary table
  St.Table.1$sp_ID <- paste(St.Table.1$Sp1,"-",St.Table.1$Sp2,sep="")
  for(i in 1:100){
    Data.TS=Data.Sub[[i]]
    St.TS=as.matrix(simpairs.score(Data.TS))
    colnames(St.TS)<-rownames(Data.TS)
    rownames(St.TS)<-rownames(Data.TS)
    St.Table.2 <- PairsTable(St.TS)
    St.Table.2$id <- paste(St.Table.2$Sp1,"-",St.Table.2$Sp2,sep="")
    colnames(St.Table.2) <- c("sp1_name","sp2_name","FETmP","sp_ID")
#    St.Table.2 <- St.Table.2[,c("sp_ID","FETmP")]
    St.Table.2 <- St.Table.2 %>% select(sp_ID,FETmP)
#    St.Table.2[1,]
    St.Table.1 <- left_join(as.data.frame(St.Table.1),as.data.frame(St.Table.2),by="sp_ID") # Join collection info fields to species matrix
  }
  colnames(St.Table.1) <- c("sp1_name","sp2_name","Base.FETmP","sp_ID",seq(1:100))
  St.Table.1=unique(St.Table.1[,1:ncol(St.Table.1)])
  St.Table.1$Mean.FETmP <- rowMeans(St.Table.1[,5:(ncol(St.Table.1))],na.rm=TRUE)
  TS.Summary <- St.Table.1 %>% select(sp1_name,sp2_name,sp_ID,Base.FETmP,Mean.FETmP)
  
  for(i in 1:nrow(TS.Summary)){
    if(is.na(TS.Summary[i,5])==TRUE){TS.Summary[i,5]<-0.5}
    if(TS.Summary[i,5]==0.5){TS.Summary[i,6] <-"Random"}
    if(TS.Summary[i,5]>0.5){TS.Summary[i,6]<-"Agg"}
    if(TS.Summary[i,5]<0.5){TS.Summary[i,6] <-"Seg"}
  }
  
  colnames(TS.Summary) <- c("sp1_name","sp2_name","sp_ID","Base.FETmP","Mean.FETmP","Pair.Type")
  TS.Summary$Zone <- Interval

  Cooccur.Summary=as.data.frame(matrix(nrow=4,ncol=reps)) # Creates summary table for cooccurence
  Agg.Total<-data.frame()
  Seg.Total<-data.frame()
  for(i in 1:reps){
    print(i)
    Data.TS = Data.Sub[[i]]
    f = c(colnames(Data.TS))
    Cooccur.TS = cooccur(mat=Data.TS,type="spp_site",spp_names=TRUE,prob="comb",thresh=FALSE)
    Genus.TS = as.data.frame(Cooccur.TS[["results"]])
    Agg.Genus.TS = Genus.TS %>% filter(p_gt<0.05) %>% dplyr::select(sp1_name,sp2_name,p_gt)
    Seg.Genus.TS = Genus.TS %>% filter(p_lt<0.05) %>% dplyr::select(sp1_name,sp2_name,p_lt)
    Rand.Genus.TS = Genus.TS %>% filter(p_lt>0.05 & p_gt>0.05) %>% dplyr::select(sp1_name,sp2_name,prob_cooccur)
    Genus.TS.Pairs=pair.attributes(Cooccur.TS)
    
    # Create Summary Tables for co-ooccurence
    Cooccur.Summary[1,i]=nrow(Agg.Genus.TS)
    Cooccur.Summary[2,i]=nrow(Seg.Genus.TS)
    Cooccur.Summary[3,i]=nrow(Rand.Genus.TS)
    Cooccur.Summary[4,i]=Cooccur.TS$pairs
    
    Agg.Total<-rbind(Agg.Total,Agg.Genus.TS)
    Seg.Total<-rbind(Seg.Total,Seg.Genus.TS)
  }
  colnames(Agg.Total)<-c("sp1_name","sp2_name","p_value")
  colnames(Seg.Total)<-c("sp1_name","sp2_name","p_value")
  rownames(Cooccur.Summary)=c("Positive","Negative","Random","Total Pairs")
  write.csv(Cooccur.Summary, file=paste("Results/",Interval,"/100 Cooccurence Summary.csv", sep=""), row.names=T)
  
  if(nrow(Agg.Total)>0){
    Agg.Total.C <- Agg.Total %>% select(sp1_name,sp2_name)
    Agg.Total.C$Var1<-paste(Agg.Total.C$sp1_name,"-",Agg.Total.C$sp2_name,sep="")
    Agg.Total.Table<-as.data.frame(table(Agg.Total.C$Var1))
    Agg.Total.Table$Pair.Veech<-"Aggregated"
    Agg.Total.Table=left_join(as.data.frame(Agg.Total.Table),as.data.frame(Agg.Total.C),by="Var1")
    Agg.Total.Table=unique(Agg.Total.Table[,1:ncol(Agg.Total.Table)])
  }else{Agg.Total.Table<-data.frame()}
  
  if(nrow(Seg.Total)>0){
    Seg.Total.C <- Seg.Total %>% select(sp1_name,sp2_name)
    Seg.Total.C$Var1<-paste(Seg.Total.C$sp1_name,"-",Seg.Total.C$sp2_name,sep="")
    Seg.Total.Table<-as.data.frame(table(Seg.Total.C$Var1))
    Seg.Total.Table$Pair.Veech<-"Segregated"
    Seg.Total.Table=left_join(as.data.frame(Seg.Total.Table),Seg.Total.C,by="Var1")
    Seg.Total.Table=unique(Seg.Total.Table[,1:ncol(Seg.Total.Table)])
  }else{Seg.Total.Table<-data.frame()}
  
  
  rownames(Cooccur.Summary)=c("Positive","Negative","Random","Total Pairs")
  
  if(nrow(Agg.Total.Table)>0 | nrow(Seg.Total.Table)>0){ # if all pairs were random, create blank table
    Pairs.Table.V<-rbind(Agg.Total.Table,Seg.Total.Table)}else{
      Pairs.Table.V<-as.data.frame(matrix(ncol=5,nrow=1))
      Pairs.Table.V[1,]<-"No Pairs"}
  
  colnames(Pairs.Table.V)<-c("sp_ID","Freq","Pair.Veech","sp1_name","sp2_name")
  write.csv(Pairs.Table.V,file=paste("Results/",Interval,"/Species Summary ",Interval,".csv",sep=""))
  
  # Combine results from FETmP with Veech
  Pairs.Table.V <- Pairs.Table.V %>% select("sp_ID","Pair.Veech","Freq") # Select data from Veech analysis
  if(nrow(Agg.Total.Table)>0 | nrow(Seg.Total.Table)>0){ # if there were pairs
    TS.Summary  <- left_join(as.data.frame(TS.Summary),as.data.frame(Pairs.Table.V),by="sp_ID") # Join Veech to FETmP summary
    TS.Summary  <- unique(TS.Summary [,1:ncol(TS.Summary)])}else{# Remove duplicate row from Join
      TS.Summary$Pair.Veech<-NA
      TS.Summary$Freq<-NA}
  TS.Summary[is.na(TS.Summary)]<-"Null" # Designate non-significant pairs
  write.csv(TS.Summary,file=paste("Results/",Interval,"/Combined Pairs Summary.csv",sep=""))
}


Cooc.Bodge <- function(Genus.Dat,Interval,sites,reps){ # code to get around memory limitations of my PC and the Veech 2013 analysis code
  Data.Sub<-resamp.MC(Genus.Dat,sites,reps) # Subsample sites over N reps
  
  # FETmP-based Cooccurrence
  St.TS<-as.matrix(simpairs.score(Genus.Dat)) # Base FETmP
  colnames(St.TS)<-rownames(Genus.Dat)
  rownames(St.TS)<-rownames(Genus.Dat)
  St.Table.1 <- PairsTable(St.TS) # Create core of summary table
  St.Table.1$sp_ID <- paste(St.Table.1$Sp1,"-",St.Table.1$Sp2,sep="")
  for(i in 1:100){
    Data.TS=Data.Sub[[i]]
    St.TS=as.matrix(simpairs.score(Data.TS))
    colnames(St.TS)<-rownames(Data.TS)
    rownames(St.TS)<-rownames(Data.TS)
    St.Table.2 <- PairsTable(St.TS)
    St.Table.2$id <- paste(St.Table.2$Sp1,"-",St.Table.2$Sp2,sep="")
    colnames(St.Table.2) <- c("sp1_name","sp2_name","FETmP","sp_ID")
    St.Table.2 <- St.Table.2 %>% select(sp_ID,FETmP)
    St.Table.1 <- left_join(as.data.frame(St.Table.1),as.data.frame(St.Table.2),by="sp_ID") # Join collection info fields to species matrix
  }
  colnames(St.Table.1) <- c("sp1_name","sp2_name","Base.FETmP","sp_ID",seq(1:100))
  St.Table.1=unique(St.Table.1[,1:ncol(St.Table.1)])
  St.Table.1$Mean.FETmP <- rowMeans(St.Table.1[,5:(ncol(St.Table.1))],na.rm=TRUE)
  TS.Summary <- St.Table.1 %>% select(sp1_name,sp2_name,sp_ID,Base.FETmP,Mean.FETmP)
  
  for(i in 1:nrow(TS.Summary)){
    if(is.na(TS.Summary[i,5])==TRUE){TS.Summary[i,5]<-0.5}
    if(TS.Summary[i,5]==0.5){TS.Summary[i,6] <-"Random"}
    if(TS.Summary[i,5]>0.5){TS.Summary[i,6]<-"Agg"}
    if(TS.Summary[i,5]<0.5){TS.Summary[i,6] <-"Seg"}
  }
  
  colnames(TS.Summary) <- c("sp1_name","sp2_name","sp_ID","Base.FETmP","Mean.FETmP","Pair.Type")
  TS.Summary$Zone <- Interval
  
  Cooccur.Summary=as.data.frame(matrix(nrow=4,ncol=reps)) # Creates summary table for cooccurence
  Agg.Total<-data.frame()
  Seg.Total<-data.frame()

  if(nrow(Agg.Total)>0){
    Agg.Total.C <- Agg.Total %>% select(sp1_name,sp2_name)
    Agg.Total.C$Var1<-paste(Agg.Total.C$sp1_name,"-",Agg.Total.C$sp2_name,sep="")
    Agg.Total.Table<-as.data.frame(table(Agg.Total.C$Var1))
    Agg.Total.Table$Pair.Veech<-"Aggregated"
    Agg.Total.Table=left_join(as.data.frame(Agg.Total.Table),as.data.frame(Agg.Total.C),by="Var1")
    Agg.Total.Table=unique(Agg.Total.Table[,1:ncol(Agg.Total.Table)])
  }else{Agg.Total.Table<-data.frame()}
  
  if(nrow(Seg.Total)>0){
    Seg.Total.C <- Seg.Total %>% select(sp1_name,sp2_name)
    Seg.Total.C$Var1<-paste(Seg.Total.C$sp1_name,"-",Seg.Total.C$sp2_name,sep="")
    Seg.Total.Table<-as.data.frame(table(Seg.Total.C$Var1))
    Seg.Total.Table$Pair.Veech<-"Segregated"
    Seg.Total.Table=left_join(as.data.frame(Seg.Total.Table),Seg.Total.C,by="Var1")
    Seg.Total.Table=unique(Seg.Total.Table[,1:ncol(Seg.Total.Table)])
  }else{Seg.Total.Table<-data.frame()}
  
  
  rownames(Cooccur.Summary)=c("Positive","Negative","Random","Total Pairs")
  
  if(nrow(Agg.Total.Table)==0 & nrow(Seg.Total.Table)==0){ # if all pairs were random, create blank table
    Pairs.Table.V<-as.data.frame(matrix(ncol=5,nrow=1))
    Pairs.Table.V[1,]<-"No Pairs"
  }else{Pairs.Table.V<-rbind(Agg.Total.Table,Seg.Total.Table)}
  
  colnames(Pairs.Table.V)<-c("sp_ID","Freq","Pair.Veech","sp1_name","sp2_name")
  write.csv(Pairs.Table.V,file=paste("Results/",Interval,"/Species Summary ",Interval,".csv",sep=""))
  
  # Combine results from FETmP with Veech
  Pairs.Table.V <- Pairs.Table.V %>% select("sp_ID","Pair.Veech","Freq") # Select data from Veech analysis
  if(nrow(Agg.Total.Table)>0 | nrow(Seg.Total.Table)>0){ # if there were pairs
    TS.Summary  <- left_join(as.data.frame(TS.Summary),as.data.frame(Pairs.Table.V),by="sp_ID") # Join Veech to FETmP summary
    TS.Summary  <- unique(TS.Summary [,1:ncol(TS.Summary)])}else{# Remove duplicate row from Join
      TS.Summary$Pair.Veech<-NA
      TS.Summary$Freq<-NA}
  TS.Summary[is.na(TS.Summary)]<-"Null" # Designate non-significant pairs
  write.csv(TS.Summary,file=paste("Results/",Interval,"/Combined Pairs Summary.csv",sep=""))
}

#### Co-Occurence Strength ####

# Onoco Co-occurence analysis
simpairs <- function(x){ # FETmP function
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  
  occs = array()
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #FETmP Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      if(z[i,j]>=1) {z[i,j] <- 0.99999999999999994}
      z[i,j] = qnorm(z[i,j])
      z[j,i] = z[i,j]
    }
  }
  return(as.dist(z, diag = F, upper = F))
} 

# Summarize co-occurence strengths with resampling
Pairs.Strength.Summary<-function(Data.Sub,Interval,reps){
  Onoco.Summary<-data.frame()
  for(i in 1:reps){
    Data.Test<-Data.Sub[[i]]
    Onoco.Pairs<-as.matrix(simpairs(Data.Test))
    colnames(Onoco.Pairs)<-rownames(Data.Test)
    rownames(Onoco.Pairs)<-rownames(Data.Test)
    
    PT.AV=PairsTable(Onoco.Pairs)
    Onoco.Summary[i,1] <- mean(PT.AV$FETmP[PT.AV$FETmP>0])
    Onoco.Summary[i,2] <- sd(PT.AV$FETmP[PT.AV$FETmP>0])
    
    if(is.na(mean(PT.AV$FETmP[PT.AV$FETmP>0]))){Onoco.Summary[i,3:4]<-0}else{
      Onoco.Summary[i,3:4] <- range(PT.AV$FETmP[PT.AV$FETmP>0])}
    
    Onoco.Summary[i,5] <- skewness(PT.AV$FETmP[PT.AV$FETmP>0])
    Onoco.Summary[i,6] <- mean(PT.AV$FETmP[PT.AV$FETmP<(0)])
    Onoco.Summary[i,7] <- sd(PT.AV$FETmP[PT.AV$FETmP<(0)])
    
    if(is.na(mean(PT.AV$FETmP[PT.AV$FETmP<(-0)]))){Onoco.Summary[i,8:9]<-0}else{
      Onoco.Summary[i,8:9] <- range(PT.AV$FETmP[PT.AV$FETmP<(0)])}
    
    Onoco.Summary[i,10] <- skewness(PT.AV$FETmP[PT.AV$FETmP<(0)])
    colnames(Onoco.Summary) <- c("Mean.Agg","SD.Agg","LB.Agg","UB.Agg","Skew.Agg",
                                 "Mean.Seg","SD.Seg","LB.Seg","UB.Seg","Skew.Seg",
                                 "Num.Agg","Num.Seg")
  }
  print(colMeans(Onoco.Summary))
  write.csv(Onoco.Summary,file=paste("Results/",Interval,"/Pairs Strength Summary.csv",sep=""))
}

# Convert Species X species matrix to table
PairsTable <- function(x) {data.frame(Sp1=rownames(x)[row(x)[upper.tri(x)]], 
                                      Sp2=colnames(x)[col(x)[upper.tri(x)]], 
                                      FETmP=x[upper.tri(x)])}
# Convert distance matrix to table
DistTable <- function(x) {data.frame(Site1 = rep(colnames(x), each = nrow(x)), 
                                     Site2 = rep(rownames(x), ncol(x)), 
                                     Distance = as.vector(x))}

DivTable <- function(x) {data.frame(Site1 =rownames(x)[row(x)[upper.tri(x)]], 
                                     Site2 = colnames(x)[col(x)[upper.tri(x)]], 
                                     Sor.Dis = x[upper.tri(x)])}


#### Beta Diversity ####
# Calculate diversity with distance-based subsampling
Diversity.Dist <- function(Data.TS,Dist.T,Band.D,Interval,D.List,i.timer) {

  Dist.TS<-as.data.frame(Data.TS)
  Site.List <- c(levels(Dist.T$Site1))
  set.seed(1)
  samp.i <- sample(Site.List,100,replace=TRUE)
  Div.Table.Sor <- data.frame()
  Div.List <- as.vector("A")
  for(i in 1:100){
    if(i.timer==1){print(i)}
    Dist.I<- Dist.T %>% filter(Dist.T$Site1==samp.i[i] & Dist.T$Distance<Band.D)
    Dist.I$Site2 <- factor(Dist.I$Site2)
    Sub.Sites <- c(levels(Dist.I$Site2))
    Sub.Dat <- subset(Data.TS,select=Sub.Sites)
    
    if(ncol(as.data.frame(Sub.Dat))==1){
      Div.Table.Sor[1,i]<-NA
      Div.Table.Sor[2,i]<-NA
      Div.Table.Sor[3,i]<-sum(Sub.Dat[,1])
      Div.Table.Sor[4,i]<-samp.i[i]} else{
      Div.Table.Sor[1,i]<-mean(ecol.dist(Sub.Dat,method=sorenson,type="dis"))
      Div.Table.Sor[2,i]<-fossil::bootstrap(Sub.Dat,taxa.row=TRUE,abund=FALSE)
      Div.Table.Sor[3,i]<-length(Sub.Sites)
      Div.Table.Sor[4,i]<-samp.i[i]
      if(D.List==1){
      Div.List <- c(Div.List,ecol.dist(Sub.Dat,method=sorenson,type="dis"))}
      }
    
    
    rownames(Div.Table.Sor) <- c("Sorenson.Dis","Boot.Rich","Sites","Center")
    
    if(D.List==1){write.csv(Div.List,file=paste("Results/",Interval,"/Diversity List - ",Interval,".csv",sep=""))}
    write.csv(Div.Table.Sor,file=paste("Results/",Interval,"/Diversity Table - ",Interval,".csv",sep=""))
    }}
