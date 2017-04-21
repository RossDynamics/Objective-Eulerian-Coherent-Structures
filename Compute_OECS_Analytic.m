function Compute_OECS
clear all

%Global Variables
global tol x y u v dxcheck dycheck Eigenvalue1 Eigenvalue2 Eigenvector1 Eigenvector2 index Eigencheck xmin xmax ymin ymax dxev1 dxev2 dyev1 dyev2

%Input Parameters
%{
xmin=-3;
xmax=3;
ymin=-3;
ymax=3;
%}
%
xmin=-1;
xmax=1;
ymin=-1;
ymax=1;
% ideally tol = 0
tol = 10
%}
res=30;
dist=1/res;
ForwardLength=[0 200];
BackwardLength=[0 -200];
maxits=500;
%opts=odeset('MaxStep',0.1,'RelTol',1e-9,'Events', @eventfun);
%opts=odeset('MaxStep',0.1,'Events', @eventfun);
opts=odeset('Events', @eventfun);

%Calculate Velocity Field & Coordinate Mesh
[x,y,u,v]=velocity(xmin,xmax,ymin,ymax,res);

%Calculate Gradients of Velocity Field
[dxu,dyu]=gradient(u,dist);
[dxv,dyv]=gradient(v,dist);

figure
surface(x,y,dxu)
title('d/dx u')
view(30,30)

figure
surface(x,y,dyu)
title('d/dy u')
view(30,30)

figure
surface(x,y,dxv)
title('d/dx v')
view(30,30)

figure
surface(x,y,dyv)
title('d/dy v')
view(30,30)

%Compute Eulerian Rate of Strain Tensor for each point in the Coordinate Mesh
S = ComputeEulerianStrainTensor(dxu,dyu,dxv,dyv);

%Calculate Eigenvalue and Eigenvector fields for Eulerian Strain Tensor
[Eigenvalue1 Eigenvalue2 Eigenvector1 Eigenvector2]=CalculateEigenFields(S);

[dxev1,dyev1]=gradient(abs(Eigenvalue1),dist);
%{
figure
surface(x,y,dxev1)

figure
surface(x,y,dyev1)
%}
[dxev2,dyev2]=gradient(abs(Eigenvalue2),dist);
%{
figure
surface(x,y,dxev2)

figure
surface(x,y,dyev2)
%}

%Compute Hyperbolic OECS
%Set up grid vectors
x1=x(1,:);
y1=y(:,1);
[s1xmin s1ymin]=Find2DPeak(abs(Eigenvalue1(:,:)),x1,y1,'maxima');

[s2xmax s2ymax]=Find2DPeak(abs(Eigenvalue2(:,:)),x1,y1,'maxima');


figure
hold on
quiver(x,y,u,v,'k')
for i=[1:1:length(s2xmax)]
    Ax0=s2xmax(i);
    Ay0=s2ymax(i);
    A0=[Ax0,Ay0,0];
    clear Eigencheck;
    index=0;
    [T A1]=ode113(@Attracting,ForwardLength,A0,opts);
    clear Eigencheck;
    index=0;
    [T A2]=ode113(@Attracting,BackwardLength,A0,opts);
    %[T A2]=ode113(@Attracting2,timeA,A0,opts);
    plot(A1(:,1),A1(:,2),'b','LineWidth',2)
    plot(A2(:,1),A2(:,2),'b','LineWidth',2)
end    
%
for i=[1:1:length(s1xmin)]
    Rx0=s1xmin(i);
    Ry0=s1ymin(i);
    R0=[Rx0,Ry0,0];
    clear Eigencheck;
    index=0;
    [T R1]=ode113(@Repelling,ForwardLength,R0,opts);
    clear Eigencheck;
    index=0;
    [T R2]=ode113(@Repelling,BackwardLength,R0,opts);
    plot(R1(:,1),R1(:,2),'r','LineWidth',2)
    plot(R2(:,1),R2(:,2),'r','LineWidth',2)
end
%}
scatter(s2xmax, s2ymax,'m')%,'filled')
scatter(s1xmin, s1ymin,'m')%,'filled')
axis equal

%{
function dzdt=Attracting(t,z)
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        index=index+1;
        Eigencheck(index)=abs(interp2(x,y,Eigenvalue1,z(1),z(2)));
        [Ix2 Iy2]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
        dzdt(1) = Ix2; 
        dzdt(2) = Iy2;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
                dzdt(3)=0;
            else
                Eigencheck(index)=Eigencheck(index-1);
                dzdt(3)=1;
            end
        end
    else
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
end

function dzdt=Repelling(t,z)
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        index=index+1;
        Eigencheck(index)=abs(interp2(x,y,Eigenvalue2,z(1),z(2)));
        [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
        dzdt(1) = Ix1; 
        dzdt(2) = Iy1;
        if index > 1
            if Eigencheck(index-1)>=Eigencheck(index)
                dzdt(3)=0;
            else
                Eigencheck(index)=Eigencheck(index-1);
                dzdt(3)=1;
            end
        end
    else
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
end
%}

function dzdt=Attracting(t,z)
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        index=index+1;
        dxcheck=interp2(x,y,dxev1,z(1),z(2));
        dycheck=interp2(x,y,dyev1,z(1),z(2));
        if t<0
            dxcheck=-dxcheck;
            dycheck=-dycheck;
        end
        [Ix2 Iy2]=smooth_xi(z(1),z(2),x,y,Eigenvector2);
        dzdt(1) = Ix2; 
        dzdt(2) = Iy2;
        if index > 1
            if and(dxcheck<=tol,dycheck<=tol)
                dzdt(3)=0;
            else
                dzdt(3)=1;
            end
        end
    else
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
    if index>maxits
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
    %}
end

function dzdt=Repelling(t,z)
    dzdt = zeros(3,1);
    if and(and(z(1)<=xmax,z(1)>=xmin),and(z(2)<=ymax,z(2)>=ymin))
        index=index+1;
        dxcheck=interp2(x,y,dxev2,z(1),z(2));
        dycheck=interp2(x,y,dyev2,z(1),z(2));
        if t<0
            dxcheck=-dxcheck;
            dycheck=-dycheck;
        end
        [Ix1 Iy1]=smooth_xi(z(1),z(2),x,y,Eigenvector1);
        dzdt(1) = Ix1; 
        dzdt(2) = Iy1;
        if index > 1
            if and(dxcheck<=tol,dycheck<=tol)
                dzdt(3)=0;
            else
                dzdt(3)=1;
            end
        end
    else
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
    if index>maxits
        dzdt(1) = 0;
        dzdt(2) = 0;
        dzdt(3) = 1;
    end
    %}
end

function [value,isterminal,direction]=eventfun(t,z)
    value=[z(3)-eps,z(3)+eps];
    isterminal=[1,1];
    direction=[0,0];
end

end