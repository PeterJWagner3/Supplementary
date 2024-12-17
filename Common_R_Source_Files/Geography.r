# use this for simulated expectations (Supplementary figures)
library(gtools)
earth_radius <- 6378.16
equatorial_circumferance <- 40075
polar_circumferance <- 40008

earth_circumference_at_latitude <- function(latitude)	{
return(2*pi*(earth_radius*cos(abs(latitude/180)*pi)));
}

#lat <- plat1; lng <- plng1; km_north <- -72.41374; km_east <- -267.1877;
new_coordinates_given_km_from_lat_and_long <- function(lat,lng,km_north,km_east)	{
lat_circumf <- earth_circumference_at_latitude(lat);
new_lat <- lat+(km_east/(lat_circumf/360));
new_lng <- lng+(km_north/(polar_circumferance/360));
return(data.frame(lat=as.numeric(new_lat),lng=as.numeric(new_lng),stringsAsFactors = F));
}

#dodgy_site_info <- data.frame(lat=as.numeric(dodgy_site_info$lat),lng=as.numeric(mlngØ))
reestimate_paleocoordinates <- function(suspect_site,relv_pbdb_sites)	{
relv_pbdb_info <- data.frame(lat=as.numeric(relv_pbdb_sites$lat),
							 lng=as.numeric(relv_pbdb_sites$lng),
							 paleolat=as.numeric(relv_pbdb_sites$paleolat),
							 paleolng=as.numeric(relv_pbdb_sites$paleolng));
relv_pbdb_info <- unique(relv_pbdb_info);

poss_indentical <- relv_pbdb_info[relv_pbdb_info$lng==suspect_site$lng & relv_pbdb_info$lat==suspect_site$lat,]
if (nrow(poss_indentical)==0)	{
	estimated_coords <- data.frame(lat=as.numeric(),lng=as.numeric(),stringsAsFactors = F);
	rsites <- nrow(relv_pbdb_info);
	for (rs1 in 1:(rsites-1))	{
		mlat1 <- relv_pbdb_info$lat[rs1];
		mlng1 <- relv_pbdb_info$lng[rs1];
		plat1 <- relv_pbdb_info$paleolat[rs1];
		plng1 <- relv_pbdb_info$paleolng[rs1];
		dm1Ø <- distance_in_km_between_coordinates(suspect_site$lat,suspect_site$lng,mlat1,mlng1);
		if (suspect_site$lng!=mlng1)	{
			dm1Øew <- ((suspect_site$lng-mlng1)/abs(suspect_site$lng-mlng1))*distance_in_km_between_coordinates(suspect_site$lat,suspect_site$lng,suspect_site$lat,mlng1);
			} else	{
			dm1Øew <- 0;
			}
		if (suspect_site$lat!=mlat1)	{
			dm1Øns <- ((suspect_site$lat-mlat1)/abs(suspect_site$lat-mlat1))*distance_in_km_between_coordinates(suspect_site$lat,suspect_site$lng,mlat1,suspect_site$lng);
			} else	{
			dm1Øns <- 0;
			}
		polar_coord1Ø <- cartesian_to_polar_coordinates(dm1Øew,dm1Øns);
		angle1Øm <- polar_coord1Ø$theta;
		for (rs2 in (rs1+1):rsites)	{
			mlat2 <- relv_pbdb_info$lat[rs2];
			mlng2 <- relv_pbdb_info$lng[rs2];
			if (mlat2!=mlat1 || mlng1!=mlng2)	{
				plat2 <- relv_pbdb_info$paleolat[rs2];
				plng2 <- relv_pbdb_info$paleolng[rs2];
				dØ2 <- distance_in_km_between_coordinates(suspect_site$lat,suspect_site$lng,mlat2,mlng2);
				d12 <- distance_in_km_between_coordinates(mlat1,mlng1,mlat2,mlng2);
				if (mlng2!=mlng1)	{
					dm12ew <- ((mlng2-mlng1)/abs(mlng2-mlng1))*distance_in_km_between_coordinates(mlat1,mlng1,mlat1,mlng2);
					} else	{
					dm12ew <- 0;
					}
				if (mlat2!=mlat1)	{
					dm12ns <- ((mlat2-mlat1)/abs(mlat2-mlat1))*distance_in_km_between_coordinates(mlat1,mlng1,mlat2,mlng1);
					} else	{
					dm12ns <- 0;
					}
				if (plng2!=plng1)	{
					dp12ew <- ((plng2-plng1)/abs(plng2-plng1))*distance_in_km_between_coordinates(plat1,plng1,plat1,plng2);
					} else	{
					dp12ew <- 0;
					}
				if (plat2!=plat1)	{
					dp12ns <- ((plat2-plat1)/abs(plat2-plat1))*distance_in_km_between_coordinates(plat1,plng1,plat2,plng1); 
					} else	{
					dp12ns <- 0;
					}
#				print(c(rs1,rs2,dm12ew,dm12ns,dp12ew,dp12ns));
				polar_coord12m <- cartesian_to_polar_coordinates(x=dm12ew,y=dm12ns);
				angle12m <- polar_coord12m$theta;
				polar_coord12p <- cartesian_to_polar_coordinates(x=dp12ew,y=dp12ns);
				angle12p <- polar_coord12p$theta;
				angle1Øp <- angle12p+(angle12m-angle1Øm);
				km_north <- cos(angle1Øp)*dm1Ø;
				km_east <- sin(angle1Øp)*dm1Ø;
				estimated_coords <- rbind(estimated_coords,new_coordinates_given_km_from_lat_and_long(plat1,plng1,km_north,km_east));
				}
			}
		}
	} else	{
	estimated_coords <- data.frame(lat=as.numeric(poss_indentical$paleolat),lng=as.numeric(poss_indentical$paleolng),stringsAsFactors = F);	
	}
return(data.frame(lat=as.numeric(median(estimated_coords$lat)),lng=as.numeric(median(estimated_coords$lng)),stringsAsFactors = F));
}
#lat1 <- 45.783
#lat2 <- -30.652
#lng1 <- 30.1
#lng2 <- 50.0
distance_in_km_between_coordinates <- function(lat1, lng1, lat2, lng2)	{
#double	X1, X2, Y1, Y2, Z1, Z2;
#double	dotproduct, theta, distancekm=0.0f;
PI <- 3.14159265358979
RAD <- 6370.9
if (lat1!=lat2 || lng1!=lng2)	{
	X1 <- (cos(lat1*(PI/180)))*(cos(lng1*(PI/180)));
	Y1 <- (cos(lat1*(PI/180)))*(sin(lng1*(PI/180)));
	Z1 <- sin(lat1*(PI/180));

	X2 <- (cos(lat2*(PI/180)))*(cos(lng2*(PI/180)));
	Y2 <- (cos(lat2*(PI/180)))*(sin(lng2*(PI/180)));
	Z2 <- sin(lat2*(PI/180));

	dotproduct <- (X1*X2)+(Y1*Y2)+(Z1*Z2);

	theta=acos(dotproduct)*180/PI;

	if (theta>180)	theta <- 360-theta;	# if you go over halfway around the world, then reverse!

	distancekm <- (theta/360)*2*PI*RAD;
	}	else	distancekm <- 0.0

return (distancekm)
}

distances_in_km_among_many_coordinates <- function(lats,longs)	{
sites <- length(lats);
distances <- array(0,dim=c(sites,sites));
for (s1 in 1:(sites-1))	{
	for (s2 in (s1+1):sites)	{
		distances[s2,s1] <- distances[s1,s2] <- distance_in_km_between_coordinates(lat1=lats[s1],lng1=longs[s1],lat2=lats[s2],lng2=longs[s2]);
		}
	}
return(distances);
}

mean_distance_in_km_between_coordinates <- function(lats,longs)	{
#double	X1, X2, Y1, Y2, Z1, Z2;
#double	dotproduct, theta, distancekm=0.0f;
avekm <- 0
ttlcoll <- length(lats)
for (p1 in 1:(ttlcoll-1))	{
	for (p2 in (p1+1):ttlcoll)	{
		avekm <- avekm+distance_in_km_between_coordinates(lats[p1],longs[p1],lats[p2],longs[p2])
		}
	}
avekm <- avekm/(((ttlcoll^2)-ttlcoll)/2)
return(avekm)	
}

max_distance_in_km_between_coordinates <- function(lats,longs)	{
#double	X1, X2, Y1, Y2, Z1, Z2;
#double	dotproduct, theta, distancekm=0.0f;
maxkm <- 0
ttlcoll <- length(lats)
for (p1 in 1:(ttlcoll-1))	{
	for (p2 in (p1+1):ttlcoll)	{
		km <- distance_in_km_between_coordinates(lats[p1],longs[p1],lats[p2],longs[p2])
		if (km>maxkm)	maxkm <- km
		}
	}
return(maxkm)	
}

#latlong1 <- coordinates_i; latlong2 <- coordinates_j;
distances_in_km_between_two_areas <- function(latlong1,latlong2)	{
# latlong1: coordinates for site in one area / rock unit / whatever
# latlong2: coordinates for site in other area / rock unit / whatever
area_dists <- c();
for (i in 1:nrow(latlong1))
	for (j in 1:nrow(latlong2))
		area_dists <- c(area_dists,distance_in_km_between_coordinates(lat1=latlong1[i,1],lng1 = latlong1[i,2],lat2=latlong2[j,1],lng2=latlong2[j,2]));
return(area_dists);
}

convert_degrees_minutes_seconds_to_degrees_fraction <- function(degrees,minutes,seconds)	{
return(degrees+(minutes/60)+(seconds/3600))
}

convert_degrees_fraction_to_degrees_minutes_seconds <- function(degrees)	{
sign <- degrees/abs(degrees)
d <- sign*floor(abs(degrees))
m <- floor(60*(abs(degrees)-abs(d)))
s <- round(60*(60*(abs(degrees)-abs(d))-m),0)
return(c(d,m,s))
}

contemporaneous_collections <- function(lb1,ub1,lb2,ub2)	{
range1 <- (10*lb1):(10*ub1)
range2 <- (10*lb2):(10*ub2)
overlap <- range1[range1 %in% range2];
if (length(overlap)>0)	{
	return(TRUE);
	} else	{
	return(FALSE);
	}
}

contemporaneity_among_collections <- function(lbs,ubs)	{
sites <- length(lbs);
contemporaneity <- array(TRUE,dim=c(sites,sites));
for (s1 in 1:(sites-1))	{
	lb1 <- lbs[s1];
	ub1 <- ubs[s1];
	for (s2 in (s1+1):sites)	{
		lb2 <- lbs[s2];
		ub2 <- ubs[s2];
		contemporaneity[s1,s2] <- contemporaneity[s2,s1] <- contemporaneous_collections(lb1,ub1,lb2,ub2);
		}
	}
return(contemporaneity);
}

# routine to lump collections & within a certain distance ("lump")
binge_sampling_lumping <- function(coll_geog,LOC,FORM,LAT,LNG,LB,UB,lump)	{
# Provides two columns: column 1 is original locality #, column 2 is the post-lumping number.  (These usually are the same.)
# this will be easiest if localities are sorted
# coll_geog: matrix giving information about localities
# LOC: column number giving collection number
# FORM: column number giving rock-unit number
#	NOTE: if all values for FORM are the same, then all localities are lumped within a span regardlesss of rock-unit
# LAT: column number giving collection latitude
# LNG: column number giving collection longitude
# LB: oldest possible age in myr
# UB: youngest possible age in myr
# lump: radius within which to lump collections

coll_geog <- coll_geog[order(coll_geog[,LOC],decreasing=FALSE),]
lumped_coll <- coll_geog[,LOC]	# this will be the "final" number for each collection
ttl_coll <- nrow(coll_geog);		# total collections considered here

rock_nos <- sort(unique(coll_geog[,FORM]));
#rock_nos <- rock_nos[rock_nos>0];
coll_nos <- data.frame(collection_no=as.numeric(coll_geog$collection_no),
					   collection_no_sr=as.numeric(coll_geog$collection_no));

for (rr in 1: length(rock_nos))	{
	rock_collections <- subset(coll_geog,coll_geog[,FORM]==rock_nos[rr]);
	
	if (nrow(rock_collections)>1)	{
		lats <- rock_collections[,LAT];
		longs <- rock_collections[,LNG];
		distances <- distances_in_km_among_many_coordinates(lats,longs);
		lbs <- rock_collections[,LB];
		ubs <- rock_collections[,UB];
		contemporaneity <- contemporaneity_among_collections(lbs,ubs);

		rock_coll_nos <-rock_collections$collection_no;
		to_lump <- distances<=lump;
		colnames(to_lump) <- rownames(to_lump) <- rock_collections$collection_no;
		to_lump <- to_lump*contemporaneity;		
		for (cc in 1:nrow(distances))	{
			rcoll <- match(rock_coll_nos[cc],coll_nos$collection_no);
			if (coll_nos$collection_no[rcoll]==coll_nos$collection_no_sr[rcoll])	{
				binged <- rock_coll_nos[to_lump[cc,]==TRUE];
				coll_nos$collection_no_sr[match(binged,coll_nos$collection_no)] <- binged[1];
				}
			}
		}	
	}

return(coll_nos)
}

# routine to lump collections & within a certain distance ("lump")
binge_sampling_correction <- function(coll_geog,LOC,FORM,LAT,LNG,lump)	{
# Provides two columns: column 1 is original locality #, column 2 is the post-lumping number.  (These usually are the same.)
# this will be easiest if localities are sorted
# coll_geog: matrix giving information about localities
# LOC: column number giving collection number
# FORM: column number giving rock-unit number
#	NOTE: if all values for FORM are the same, then all localities are lumped regardlesss of rock-unit
# LAT: column number giving collection latitude
# LNG: column number giving collection longitude
# lump: radius within which to lump collections

coll_geog <- coll_geog[order(coll_geog[,LOC],decreasing=FALSE),]
lumped_coll <- coll_geog[,LOC]	# this will be the "final" number for each collection
ttl_coll <- nrow(coll_geog);		# total collections considered here

for (c1 in 1:(ttl_coll-1))	{
	# skip already lumped collections
	while (lumped_coll[c1]!=coll_geog[c1,LOC] && c1<(ttl_coll-1))	c1 <- c1+1;
	if (c1<=(ttl_coll-1))	{
		rem_coll <- coll_geog[(c1+1):ttl_coll,];
		if (c1<(ttl_coll-1))	{
			same_form <- subset(rem_coll,rem_coll[,FORM]==coll_geog[c1,FORM])
			rem_form <- nrow(same_form);
			}	else	{
			if (rem_coll[FORM]==coll_geog[c1,FORM])	{
				same_form <- matrix(0,1,4)
				same_form[1,] <- coll_geog[c1+1,]
				rem_form <- 1
				}	else {
				rem_form <- 0		# MAKE SURE ABOUT THIS!!!
				}
			}
		if (rem_form>0)	{
			poss_matches <- match(same_form[,LOC],coll_geog[,LOC])	# get vector of collections that might match c1
			for (cc in 1:rem_form)	{
				c2 <- poss_matches[cc];
				# only bother if this collection has not already been lumped
				if (lumped_coll[c2]==coll_geog[c2,LOC])	{
					if (coll_geog[c1,LAT]==coll_geog[c2,LAT] && coll_geog[c1,LNG]==coll_geog[c2,LNG])	{
						km <- 0
						}	else	{
						km <- distance_in_km_between_coordinates(coll_geog[c1,LAT],coll_geog[c1,LNG],coll_geog[c2,LAT],coll_geog[c2,LNG])	
						}
					if (km<=lump)		lumped_coll[c2] <- lumped_coll[c1]
					}	# end test to see if the following collection (c2) is close tenough to be lumped with c1	
				}	# end going through other localities matching the same criterion
			}	# end case where there are other localities from the same "group" (formation, etc.)
		}	# end case where there are still localities that might need to be lumped
	}	# end search of collections

x <- coll_geog[,1]	# total kluge.....
return(cbind(x,lumped_coll))
}

# unnormalized sinc function
sinc <- function(x)	{
return(sin(x)/x)
}

#lat <- 0;
#lng <- 180;
winkel_projection <- function(lat,lng,dataframe=F)	{
lat_rad <- pi*lat/180;
lng_rad <- pi*lng/180;
winkel_phi <- acos(2/pi);
winkel_alpha <- acos(cos(lat_rad)*cos(lng_rad/2));
winkel_lng_rad <- ((lng_rad*cos(winkel_phi)) + 2*cos(lat_rad)*sin(lng_rad/2))/2;
winkel_lat_rad <- (lat_rad + (sin(lat_rad)/sinc(winkel_alpha)))/2;
if (dataframe)	{
	winkel_coords <- data.frame(lat=as.numeric(winkel_lat_rad),lng=as.numeric(winkel_lng_rad));
	return(winkel_coords)
	} else	{
	return(c(winkel_lat_rad,winkel_lng_rad));
	}
}

winkel_cartography <- function(main_title="",azimi=5,globe_lwd=1,map_lwd=0.25)	{
wink_lng_max <- winkel_projection(0,180)[2];
wink_lng_min <- -wink_lng_max;
wink_lat_max <- winkel_projection(90,0)[1];
wink_lat_min <- -wink_lat_max;
plot(NA,type='n',axes=FALSE,main=main_title,xlab="",ylab="",xlim=c(wink_lng_min,wink_lng_max),ylim=c(wink_lat_min,wink_lat_max));
longitude_markers <- seq(-180,180,by=azimi)
#for (lm in 1:length(longitude_markers))	{
lm <- 0;
while (lm < length(longitude_markers))	{
	lm <- lm+1;
	lng <- longitude_markers[lm];
	long_line <- c();
	lat_markers <- (-180:180)/2;
	for (ml in 1:length(lat_markers))	{
		lat <- lat_markers[ml];
		long_line <- rbind(long_line,winkel_projection(lat,lng));
		}
	if (abs(round(lng,0))==180)	{
		lines(long_line[,2],long_line[,1],lwd=globe_lwd);
		} else if ((abs(round(lng,0))) %% 30 == 0)	{
		lines(long_line[,2],long_line[,1],lwd=(globe_lwd+map_lwd)/2);
		} else	{
		lines(long_line[,2],long_line[,1],lwd=map_lwd);
		}
	}
latitude_markers <- seq(-90,90,by=azimi)
ml <- 0;
while (ml < length(latitude_markers))	{
	ml <- ml+1;
	lat <- latitude_markers[ml];
	lng_markers <- (-360:360)/2;
	lat_line <- c();
	for (lm in 1:length(lng_markers))	{
		lng <- lng_markers[lm];
		lat_line <- rbind(lat_line,winkel_projection(lat,lng));
		}
	if (round(abs(lat),0)==90)	{
		lines(lat_line[,2],lat_line[,1],lwd=globe_lwd);
		} else if ((round(lat,0)) %% 30 ==0)	{
		lines(lat_line[,2],lat_line[,1],lwd=(globe_lwd+map_lwd)/2);
		} else	{
		lines(lat_line[,2],lat_line[,1],lwd=map_lwd);
		}
	}
}

aux_angle_estimator <- function(lat,theta_inc=0.005)	{
# 2θ + 2sin(θ) = π*sin(lat)
try0 <- 0;
goal <- pi*sin(lat);
poss_theta <- theta_inc;
try1 <- (2*poss_theta) + sin(2*poss_theta);
while (try1 < goal)	{
	poss_theta <- poss_theta+theta_inc;
	try0 <- try1;
	try1 <- (2*poss_theta) + sin(2*poss_theta);
#	print(c(try1,try0));
	}
if (abs(goal-try1)>abs(goal-try0))
	poss_theta <- poss_theta-theta_inc;
if (poss_theta > lat)	poss_theta <- lat;
return(poss_theta);
}

aux_angle_estimator_1 <- function(lat,long,prev_aux_angle=0)	{
if (prev_aux_angle==0)	prev_aux_angle <- lat;
if (lat==abs(pi/2))	{
	return (lat);
	} else	{
	return(prev_aux_angle - ((2*prev_aux_angle) + sin(2*prev_aux_angle) - pi*sin(lat)) / (2 + 2*cos(2*prev_aux_angle)));
	}
}

# lat <- 90*pi/180; lng <- 0*pi/180
mollweide_projection <- function(lat,lng,radius=1,meridian=0,dataframe=F)	{
lat_rad <- pi*lat/180;
lng_rad <- pi*lng/180;
aux_angle <- aux_angle_estimator(lat_rad);
mollweide_lng <- (radius*2*sqrt(2)/pi)*(lng_rad-meridian)*cos(aux_angle);
mollweide_lat <- (radius*sqrt(2)/pi)*sin(aux_angle);
if (dataframe)	{
	mollweide <- data.frame(lat=as.numeric(mollweide_lat),lng=as.numeric(mollweide_lng));
	return(mollweide)
	} else	{
	return(c(mollweide_lat,mollweide_lng));
	}
}

mollweide_cartography <- function(main_title="",azimi=5,radius=1,xsize=6,ysize=4.5,globe_lwd=1,map_lwd=0.25,font="")	{
moll_lng_max <- mollweide_projection(0,180,radius=radius)[2];
moll_lng_min <- -moll_lng_max;
moll_lat_max <- mollweide_projection(90,0,radius=radius)[1];
moll_lat_min <- -moll_lat_max;
par(pin=c(xsize,ysize));
if (font=="")	{
	plot(NA,type='n',axes=FALSE,main=main_title,xlab="",ylab="",xlim=c(moll_lng_min,moll_lng_max),ylim=c(moll_lat_min,moll_lat_max));
	} else	{
	plot(NA,type='n',axes=FALSE,main=main_title,xlab="",ylab="",xlim=c(moll_lng_min,moll_lng_max),ylim=c(moll_lat_min,moll_lat_max),family=font);
	}
longitude_markers <- seq(-180,180,by=azimi)
#for (lm in 1:length(longitude_markers))	{
lm <- 0;
while (lm < length(longitude_markers))	{
	lm <- lm+1;
	lng <- longitude_markers[lm];
	long_line <- c();
	lat_markers <- (-180:180)/2;
	for (ml in 1:length(lat_markers))	{
		lat <- lat_markers[ml];
		long_line <- rbind(long_line,mollweide_projection(lat,lng,radius=radius));
		}
	if (abs(round(lng,0))==180)	{
		lines(long_line[,2],long_line[,1],lwd=globe_lwd);
		} else if ((abs(round(lng,0))) %% 30 == 0)	{
		lines(long_line[,2],long_line[,1],lwd=(globe_lwd+map_lwd)/2);
		} else	{
		lines(long_line[,2],long_line[,1],lwd=map_lwd);
		}
	}
latitude_markers <- seq(-90,90,by=azimi)
ml <- 0;
while (ml < length(latitude_markers))	{
	ml <- ml+1;
	lat <- latitude_markers[ml];
	lng_markers <- (-360:360)/2;
	lat_line <- c();
	for (lm in 1:length(lng_markers))	{
		lng <- lng_markers[lm];
		lat_line <- rbind(lat_line,mollweide_projection(lat,lng,radius=radius));
		}
	if (round(abs(lat),0)==90)	{
		lines(lat_line[,2],lat_line[,1],lwd=globe_lwd);
		} else if ((round(lat,0)) %% 30 ==0)	{
		lines(lat_line[,2],lat_line[,1],lwd=(globe_lwd+map_lwd)/2);
		} else	{
		lines(lat_line[,2],lat_line[,1],lwd=map_lwd);
		}
	}
}

#x <- dm12ew; y <- dm12ns;
#x <- -1; y <- -1;
cartesian_to_polar_coordinates <- function(x,y)	{
# counter clockwise with 3 o'clock set to zero
radius <- sqrt((x^2)+(y^2));
angle1 <- 180*acos(x/radius)/pi;
angle2 <- 180*asin(y/radius)/pi;
#print(c(x,y)); 2023-02-21: UGH! I had x & y backwards for some reason...
if (x>0 && y>0)	{
	theta <- atan(y/x);
	} else if (x<0 && y>0)	{
	theta <- pi+atan(y/x);
	} else if (x>0 && y<0)	{
	theta <- (2*pi)+atan(y/x);
	} else if (x<0 && y<0)	{
	theta <- pi+atan(y/x);
	} else if (x==0 && y>0)	{
	theta <- pi/2;
	} else if (x==0 && y<0)	{
	theta <- 3*pi/4;
	} else if (x>0 && y==0)	{
	theta <- 0;
	} else if (x<0 && y==0)	{
	theta <- pi;
	} else if (x==0 && y==0)	{
	theta <- 0;
	}

#theta <- pi*(360+180*atan(y/x)/pi)/180;
return(data.frame(radius=as.numeric(radius),theta=as.numeric(theta)));
}

accersi_gplates_coastline_reconstructions <- function(ma,model="MERDITH2021",print_update=F)	{
options(warn=-1);
geo_url <- paste("https://gws.gplates.org/reconstruct/coastlines/?&time=",ma,"&model=",model,sep="");
info_raw <- read.table(geo_url);
if (length(info_raw)>6)	{
	if (print_update)	print(paste("Getting coastal outlines for",ma,"Ma"));
	info_raw <- gsub("\\[","",info_raw);
	info_raw <- gsub("\\]","",info_raw);
	info_raw <- gsub(",","",info_raw);
	info_raw <- gsub("\\{","",info_raw);
	info_raw <- gsub("\\}","",info_raw);
	info_raw <- gsub("\"","",info_raw);
	i_r <- length(info_raw);
	new_plates_types <- (1:i_r)[info_raw %in% "type:"];
	new_plates_coords <- (1:i_r)[info_raw %in% "coordinates"];
	new_plates_shapes <- (1:i_r)[info_raw %in% "Polygon"];
	new_plates_types <- c(max(new_plates_types[new_plates_types<min(new_plates_coords)]),new_plates_types[new_plates_types>min(new_plates_coords)]);
	n_p_c <- length(new_plates_coords);
#	specify_basic_plot(180,-180,90,-90)
	coastline_coords <- list();
	for (np in 1:(n_p_c-1))	{
		coords <- as.numeric(info_raw[new_plates_coords[np]:new_plates_coords[np+1]]);
		coords <- coords[!is.na(coords)];
		coast <- cbind(coords[odd(1:length(coords))],coords[even(1:length(coords))]);
		colnames(coast) <- c("paleolng","paleolat");
		coastline_coords <- rlist::list.append(coastline_coords,coast);
	#	polygon(coast[,1],coast[,2])
		}
	np <- np+1;
	coords <- as.numeric(info_raw[new_plates_coords[np]:i_r]);
	coords <- coords[is.numeric(coords)];
	coords <- coords[!is.na(coords)];
	coast <- cbind(coords[odd(1:length(coords))],coords[even(1:length(coords))]);
	colnames(coast) <- c("paleolng","paleolat");
	coastline_coords <- rlist::list.append(coastline_coords,coast);
	options(warn=1);
	return(coastline_coords);
	} else	{
	return(NULL);
	}
}

accersi_gplates_static_plate_reconstructions <- function(ma,model="MULLER2022",print_update=F)	{
options(warn=-1);
geo_url <- paste("https://gws.gplates.org/reconstruct/static_polygons/?&time=",ma,"&model=",model,sep="");
info_raw <- read.table(geo_url);
if (length(info_raw)>6)	{
	if (print_update)	print(paste("Getting static plate polygons for",ma,"Ma"));
	info_raw <- gsub("\\[","",info_raw);
	info_raw <- gsub("\\]","",info_raw);
	info_raw <- gsub(",","",info_raw);
	info_raw <- gsub("\\{","",info_raw);
	info_raw <- gsub("\\}","",info_raw);
	info_raw <- gsub("\"","",info_raw);
	i_r <- length(info_raw);
	new_plates_types <- (1:i_r)[info_raw %in% "type:"];
	new_plates_coords <- (1:i_r)[info_raw %in% "coordinates"];
	new_plates_shapes <- (1:i_r)[info_raw %in% "Polygon"];
	new_plates_types <- c(max(new_plates_types[new_plates_types<min(new_plates_coords)]),new_plates_types[new_plates_types>min(new_plates_coords)]);
	n_p_c <- length(new_plates_coords);
#	specify_basic_plot(180,-180,90,-90)
	polygon_coordinates <- list();
	for (np in 1:(n_p_c-1))	{
		coords <- as.numeric(info_raw[new_plates_coords[np]:new_plates_coords[np+1]]);
		coords <- coords[!is.na(coords)];
		coast <- cbind(coords[odd(1:length(coords))],coords[even(1:length(coords))]);
		colnames(coast) <- c("paleolng","paleolat");
		polygon_coordinates <- rlist::list.append(polygon_coordinates,coast);
	#	polygon(coast[,1],coast[,2])
		}
	np <- np+1;
	coords <- as.numeric(info_raw[new_plates_coords[np]:i_r]);
	coords <- coords[is.numeric(coords)];
	coords <- coords[!is.na(coords)];
	coast <- cbind(coords[odd(1:length(coords))],coords[even(1:length(coords))]);
	colnames(coast) <- c("paleolng","paleolat");
	polygon_coordinates <- rlist::list.append(polygon_coordinates,coast);
	options(warn=1);
	return(polygon_coordinates);
	} else	{
	return(NULL);
	}
}

#specify_basic_plot(180,-180,90,-90,xsize=6,ysize=4.5)
#i <- 1100;
#while (i < length(coastline_coords))	{
#i <- i+1
#while (i %% 100 !=99)	{
#	polygon(coastline_coords[[i]][,1],coastline_coords[[i]][,2]);
#	i <- i+1;
#	}
