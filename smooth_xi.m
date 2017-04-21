function [xi1_1, xi1_2]=smooth_xi1_2(x,y,xi,yi,v1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Code generated by Mattia Serra                     %
%                  Mattia Serra (ETH Zurich)                        %                               
%                http://www.zfm.ethz.ch/~serra/                     % 
%                    email: serram@ethz.ch                          % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help
% locally smoothes a vector field v1 defined on a grid (xi,yi)
% Then interpolates the vector field to find the value at (x,y)
% xi is a matrix with the x-coordinares in the meshgrid format
% yi is a matrix with the y-coordinares in the meshgrid format
% v1 is the direction field to be integrated in the form of a 3D matrix:
%  v1(i1,i2,i3): i1 and i2 are the idices related to the location in the
%  grid, and i3 (which can be 1 or 2) denotes the entries of the v1 vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin=min(min(xi)); ymin=min(min(yi));
xmax=max(max(xi)); ymax=max(max(yi));
[n1 n2]=size(xi);
dx=(xmax-xmin)/(n2-1); dy=(ymax-ymin)/(n1-1);

n=length(x);
xi1_1=zeros(n,1); xi1_2=zeros(n,1);
for j=1:n
    i2=floor((x(j)-xmin)/dx)+1;
    i1=floor((y(j)-ymin)/dy)+1;
    if (i1==n1 || i2==n2 || i1==0 || i2==0)
        xi1_1(j,1)=0; xi1_2(j,1)=0; % identifies where you are on the grid
    else
        % local smoothing
        if (v1(i1,i2,1)*v1(i1,i2+1,1)+v1(i1,i2,2)*v1(i1,i2+1,2))<0
           v1(i1,i2+1,:)=-v1(i1,i2+1,:); 
        end
        if (v1(i1,i2,1)*v1(i1+1,i2,1)+v1(i1,i2,2)*v1(i1+1,i2,2))<0
           v1(i1+1,i2,:)=-v1(i1+1,i2,:); 
        end
        if (v1(i1,i2,1)*v1(i1+1,i2+1,1)+v1(i1,i2,2)*v1(i1+1,i2+1,2))<0
           v1(i1+1,i2+1,:)=-v1(i1+1,i2+1,:); 
        end
        % Bilinear interpolation for v1(:,:,1)
        tmp1 = (xi(i1,i2+1)-x(j))*v1(i1,i2,1)/dx + (x(j)-xi(i1,i2))*v1(i1,i2+1,1)/dx;
        tmp2 = (xi(i1,i2+1)-x(j))*v1(i1+1,i2,1)/dx + (x(j)-xi(i1,i2))*v1(i1+1,i2+1,1)/dx;
        xi1_1(j)= (yi(i1+1,i2)-y(j))*tmp1/dy + (y(j)-yi(i1,i2))*tmp2/dy;
        % Bilinear interpolation for v1(:,:,2)        
        tmp1 = (xi(i1,i2+1)-x(j))*v1(i1,i2,2)/dx + (x(j)-xi(i1,i2))*v1(i1,i2+1,2)/dx;
        tmp2 = (xi(i1,i2+1)-x(j))*v1(i1+1,i2,2)/dx + (x(j)-xi(i1,i2))*v1(i1+1,i2+1,2)/dx;
        xi1_2(j)= (yi(i1+1,i2)-y(j))*tmp1/dy + (y(j)-yi(i1,i2))*tmp2/dy;
        % normalize
        nrm_xi1 = sqrt(xi1_1(j,1)^2+xi1_2(j,1)^2);
        if (nrm_xi1 ~= 0)
            xi1_1(j,1)=xi1_1(j,1)/nrm_xi1;
            xi1_2(j,1)=xi1_2(j,1)/nrm_xi1;
        end    
    end
end