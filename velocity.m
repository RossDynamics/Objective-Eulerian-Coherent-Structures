function [x,y,u,v] = velocity(xmin,xmax,ymin,ymax,resolution)
%DOUBLEGYREVELOCITY Summary of this function goes here
%   Detailed explanation goes here
%   xdot = -pi*sin(pi*x)*cos(pi*y)
%   ydot = pi*cos(pi*x)*sin(pi*y)
MAG=1e0;
[x,y] = meshgrid(xmin:1/resolution:xmax,ymin:1/resolution:ymax);
%{
%Autonomous Double Gyre
u=-pi.*sin(pi.*x).*cos(pi.*y);
v=pi.*cos(pi.*x).*sin(pi.*y);
%}

%Modified Autonomous Double Gyre
u=-pi.*sin(pi.*x).*cos(pi.*y)+MAG.*randn(size(x));%+1.*x;
v=pi.*cos(pi.*x).*sin(pi.*y)+MAG.*randn(size(x));%-1.*y;

%u=x.*y;
%v=x.*y;
%u = -pi*sin(pi*f(x,1)).*cos(pi*y);
%v = pi*cos(pi*f(x,1)).*sin(pi*y).*2.*sin(1).*x+(1-2.*sin(1));
%u=x+x^2-y;
%v=-y-y^3+x;
%figure
%surf(x,y,u)
%figure
%surf(x,y,v)
%figure
%quiver(x,y,u,v,'k')
%axis equal
end

function value=f(x,t)
    value=sin(t).*x.^2+(1-2*sin(t)).*x
end
