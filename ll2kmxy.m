function [kmx, kmy]=ll2kmxy(slat,slon,olat,olon)
%Function was stolen from Whoi.edu and converted by PETER NOLAN @ Virginia
%Tech, 2016. Converts latitude and longitude into kilometer distances.

%latitude/longitude
%slat=1;
%slon=1;
      
%Origin latitude/longitude
%olat=0;
%olon=0;

%origin=[slat,slon,coord_system,olat,olon,xoffset_mtrs,yoffset_mtrs,rotation_angle_degs,rms_error];
 origin=[slat,slon,1,           olat,olon,0,           0,           0,                  0];

sxy = translate_coordinates('LL_TO_XY',origin); 

kmx=METERS_TO_KLMS(sxy(1));
kmy=METERS_TO_KLMS(sxy(2));
end

function [sxy]=translate_coordinates(trans_option,porg)


    angle = DEG_TO_RADIANS(porg(8));
    if(trans_option == 'XY_TO_LL')
      %{
         /* X,Y to Lat/Lon Coordinate Translation  */
         pxpos_mtrs = porg.x;  
         pypos_mtrs = porg.y;
	 xx = pxpos_mtrs - porg.xoffset_mtrs;
	 yy = pypos_mtrs - porg.yoffset_mtrs;
	 r = sqrt(xx*xx + yy*yy);

	 if(r)
	 {
            ct = xx/r;
	    st = yy/r;
	    xx = r * ( (ct * cos(angle))+ (st * sin(angle)) );
	    yy = r * ( (st * cos(angle))- (ct * sin(angle)) );
	 }

	 var plon = porg.olon + xx/METERS_DEGLON(olat);
	 var plat = porg.olat + yy/METERS_DEGLAT(olat);

         var sll={};
         sll={slat:plat, slon:plon};
         return(sll);
      }
      %}
    elseif(trans_option == 'LL_TO_XY')
        xx = (porg(2) - porg(5))*METERS_DEGLON(porg(5));
        yy = (porg(1) - porg(4))*METERS_DEGLAT(porg(4));

        r = sqrt(xx*xx + yy*yy);


        if(r)
            ct = xx/r;
            st = yy/r;
            xx = r * ( (ct * cos(angle)) + (st * sin(angle)) );
            yy = r * ( (st * cos(angle)) - (ct * sin(angle)) );
        end
            
            
        pxpos_mtrs = xx + porg(6);
        pypos_mtrs = yy + porg(7);

        sxy=[pxpos_mtrs,pypos_mtrs];

    end
end

function meters=METERS_DEGLON(x)
    d2r=DEG_TO_RADIANS(x);
    meters=((111415.13 * cos(d2r))- (94.55 * cos(3.0*d2r)) + (0.12 * cos(5.0*d2r)));
end

function meters=METERS_DEGLAT(x)
    d2r=DEG_TO_RADIANS(x);
    meters=(111132.09 - (566.05 * cos(2.0*d2r))+ (1.20 * cos(4.0*d2r)) - (0.002 * cos(6.0*d2r)));
end

function rad=DEG_TO_RADIANS(x)
    RADIANS = 57.2957795;
    rad=(x/RADIANS); 
end

function kms=METERS_TO_KLMS(x)
    kms=x*0.001;
end