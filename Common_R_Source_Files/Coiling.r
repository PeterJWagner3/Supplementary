accersi_overlie <- function(theta_loc,incr,up=TRUE)	{
# get "click" number of part of shell directly behind this part of the whorl
#		up=TRUE gives it for the subsequent half whorl
#			if this click is 160˚, then we want 200˚ (180+[180-160])
#		up=FALSE gives it for the prior whorl
#			if this click is 520˚, the we want 160˚
aps <- 360/incr
### get increments away from whorl border
#rot <- (((theta_loc-1)*incr)%%360)/incr
rot <- ((theta_loc*incr)%%360)/incr	# no. of clicks past onset of whorl
if (up)	{
	return (theta_loc + 2* (aps - rot))
	} else	{
	return(theta_loc + 2* (aps - rot) - aps)	
	}
#(180-(theta_loc*incr)%%360)
### get increments away from whorl border
#rot_start <- aps*(floor((theta_loc-1)/aps))
#rot_start <- aps*(floor(theta_loc/aps))
#next_rot_start <- rot_start+aps
#overlie <- next_rot_start - rot
}

cartesian_to_polar <- function(x, y)	{
if (x>=0)	{
	if (y>=0)	{
		azimuth <- atan(y/x)
		}	else {
		azimuth <- (2*pi)+atan(y/x)
		}
	}	else	{
	if (y>=0)	{
		azimuth <- pi+atan(y/x)
		}	else	{
		azimuth <- pi+atan(y/x)	
		}
	}
polar_coords <- c(sqrt(x^2 + y^2),azimuth)
names(polar_coords) <- c("radius","azimuth")
return(polar_coords)
}

polar_to_cartesian <- function(r, theta)	{
x <- r*cos(theta);
y <- r*sin(theta);
cartesian_coords <- c(x,y);
names(cartesian_coords) <- c("x","y");
return(cartesian_coords);
}

ellipse_coordinates <- function(a,b,radius,circum,center)	{
ellipse <- matrix(0,length(circum),2)
for (d in 1:length(circum)) {
	ellipse[d,1] <- center[1]+a*cos(circum[d])*radius
	ellipse[d,2] <- center[2]+b*sin(circum[d])*radius 
	}
return(ellipse)
}

relative_radius_at_any_azimuth_around_ellipse <- function(azimuth,a,b,mean_radius)	{
# mean_radius: radius of a circle with the same area
x <- a*cos(azimuth)
y <- b*sin(azimuth)
radial_dist <- mean_radius*sqrt((x^2)+(y^2))
names(radial_dist) <- ""
return(radial_dist)
}

logarithmic_spiral  <- function(a,b,theta)	{
return(a*exp(b*theta));
}

logarithmic_spiral_from_W  <- function(W,theta,a=1)	{
b <- log(W)/(2*pi);
return(a*exp(b*theta));
#a*exp(theta*log(W)/(2*pi))
}

arc_length <- function(b,theta=2*pi,a=1)	{
W <- exp(b*2*pi).
return(((a*sqrt(1+(b^2))*exp(b*theta))/b)-(((a/2)*sqrt(1+(b^2))*exp(b*theta))/b));
}

arc_length_given_W <- function(W,theta=2*pi,a=1)	{
#W <- exp(b*2*pi);
b <- log(W)/(2*pi);
return(((a*sqrt(1+(b^2))*exp(b*theta))/b)-(((a/2)*sqrt(1+(b^2))*exp(b*theta))/b));
}

curvature <- function(b,a=1,theta=2*pi)	{
return(exp(-b*theta)/(a*sqrt(1+(b^2))));
}

curvature_given_W <- function(W,a=1,theta=2*pi)	{
return(exp(-(log(W)/(2*pi))*theta)/(a*sqrt(1+((log(W)/(2*pi))^2))));
}

#length_of_spiral <- function(W,theta=2*pi,r=1,T=0,rotations=1)	{
#if (rotations!=1)	theta <- 2*pi*rotations;
#((2*pi*W^(theta/(2*pi)))/log(W));
#dx <- 2*pi/360;
#x <- 0;
#recti
#while (x <= theta)	{
#	2*pi*r*W^x;
#	x <- x+dx;
#	}
#}
